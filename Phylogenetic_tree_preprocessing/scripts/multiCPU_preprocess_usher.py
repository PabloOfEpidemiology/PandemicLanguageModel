# Copyright Contributors to the Pyro-Cov project.
# SPDX-License-Identifier: Apache-2.0

import argparse
import datetime
import logging
import math
import pickle
import re
from collections import Counter, defaultdict

import pandas as pd
import torch
import tqdm

from pyrocov.geo import get_canonical_location_generator, gisaid_normalize
from pyrocov.mutrans import START_DATE
from pyrocov.sarscov2 import nuc_mutations_to_aa_mutations
from pyrocov.usher import (
    FineToMeso,
    load_mutation_tree,
    load_proto,
    prune_mutation_tree,
    refine_mutation_tree,
)
from pyrocov.util import gzip_open_tqdm
import multiprocessing
import os

logger = logging.getLogger(__name__)
logging.basicConfig(format="%(relativeCreated) 9d %(message)s", level=logging.INFO)

DATE_FORMATS = {7: "%Y-%m", 10: "%Y-%m-%d"}


def try_parse_date(string):
    fmt = DATE_FORMATS.get(len(string))
    if fmt is not None:
        try:
            return datetime.datetime.strptime(string, fmt)
        except ValueError:
            return


def try_parse_genbank(strain):
    match = re.search(r"([A-Z]+[0-9]+)\.[0-9]", strain)
    if match:
        return match.group(1)


def try_parse_gisaid(string, public_to_gisaid):
    match = re.search(r"\bEPI_ISL_[0-9]+\b", string)
    if match:
        return match.group(0)
    for part in string.split("|"):
        result = public_to_gisaid.get(part)
        if result is not None:
            return result


def load_nextstrain_metadata(args):
    """
    Returns a dict of dictionaries from genbank_accession to metadata.
    """
    logger.info("Loading nextstrain metadata")
    get_canonical_location = get_canonical_location_generator(
        args.recover_missing_usa_state
    )
    df = pd.read_csv("results/nextstrain/metadata.tsv", sep="\t", dtype=str)
    result = defaultdict(dict)
    for row in tqdm.tqdm(df.itertuples(), total=len(df)):
        # Key on genbank accession.
        key = row.genbank_accession
        if not isinstance(key, str) or not key:
            continue

        # Extract date.
        date = row.date
        if isinstance(date, str) and date and date != "?":
            date = try_parse_date(date)
            if date is not None:
                if date < args.start_date:
                    date = args.start_date  # Clip rows before start date.
                result["day"][key] = (date - args.start_date).days

        # Extract a standard location.
        location = get_canonical_location(
            row.strain, row.region, row.country, row.division, row.location
        )
        if location is not None:
            location = gisaid_normalize(location)
            result["location"][key] = location

        # Extract pango lineage.
        lineage = row.pango_lineage
        if isinstance(lineage, str) and lineage and lineage != "?":
            result["lineage"][key] = lineage

    logger.info("Found metadata:\n{}".format({k: len(v) for k, v in result.items()}))
    return result


# These are merely fallback locations in case nextrain is missing genbank ids.
USHER_LOCATIONS = {
    "England": "Europe / United Kingdom / England",
    "Wales": "Europe / United Kingdom / Wales",
    "Scotland": "Europe / United Kingdom / Scotland",
    "Northern Ireland": "Europe / United Kingdom / Northern Ireland",
    "China": "Asia / China",
    "Pakistan": "Asia / Pakistan",
}


def load_usher_metadata(args):
    """
    Returns a dict of dictionaries from usher strain id to metadata.
    """
    logger.info("Loading usher metadata")
    df = pd.read_csv("results/usher/metadata.tsv", sep="\t", dtype=str)
    result = defaultdict(dict)
    for row in tqdm.tqdm(df.itertuples(), total=len(df)):
        # Key on usher strain.
        key = row.strain
        assert isinstance(key, str) and key

        # Extract date.
        date = row.date
        if isinstance(date, str) and date and date != "?":
            date = try_parse_date(date)
            if date is not None:
                if date < args.start_date:
                    date = args.start_date  # Clip rows before start date.
                result["day"][key] = (date - args.start_date).days

        # Extract a standard location.
        if isinstance(row.strain, str) and row.strain:
            prefix = re.split(r"[/|_]", row.strain)[0]
            location = USHER_LOCATIONS.get(prefix)
            if location is not None:
                location = gisaid_normalize(location)
                result["location"][key] = location

        # Extract pango lineage.
        lineage = row.pangolin_lineage
        if isinstance(lineage, str) and lineage and lineage != "?":
            result["lineage"][key] = lineage

    logger.info("Found metadata:\n{}".format({k: len(v) for k, v in result.items()}))
    return result


def load_gisaid_metadata(args):
    """
    Returns a dict of dictionaries from gisaid accession to metadata.
    """
    filename = args.gisaid_metadata_file_in
    logger.info(f"Loading gisaid metadata from {filename}")
    assert filename.endswith(".tsv.gz") or filename.endswith(".tsv")

    # Open the file depending on its format:
    if filename.endswith(".gz"):
        f = gzip_open_tqdm(filename, "rt")
    else:
        # For a plain .tsv file, use a normal open and wrap with tqdm manually if needed
        f = open(filename, "rt")

    result = defaultdict(dict)
    header = ()
    for i, line in enumerate(f):
        line = line.strip().split("\t")
        if i == 0:
            header = tuple(line)
            continue
        row = dict(zip(header, line))

        # Key on gisaid accession id.
        key = row.get("Accession ID")
        if not key:
            continue

        # Extract date.
        date = row.get("Collection date")
        if date and date != "?":
            date = try_parse_date(date)
            if date is not None:
                if date < args.start_date:
                    date = args.start_date  # Clip rows before start date.
                result["day"][key] = (date - args.start_date).days

        # Extract location.
        location = row.get("Location")
        if location:
            location = gisaid_normalize(location)
            result["location"][key] = location

        # Extract pango lineage.
        lineage = row.get("Pango lineage")
        if lineage:
            result["lineage"][key] = lineage

    f.close()  # Close the file handle

    logger.info("Found metadata:\n{}".format({k: len(v) for k, v in result.items()}))
    return result


def load_metadata(args):
    # Load metadata.
    public_to_gisaid = {}
    if args.gisaid_metadata_file_in:
        metadata = load_gisaid_metadata(args)
        with open("results/epiToPublicAndDate.latest", "rt") as f:
            for line in f:
                row = line.strip().split()
                if row:
                    public_to_gisaid[row[1]] = row[0]
    else:
        # Use nextstrain metadata when available; otherwise fallback to usher.
        usher_metadata = load_usher_metadata(args)  # keyed on strain
        nextstrain_metadata = load_nextstrain_metadata(args)  # keyed on genbank
        metadata = usher_metadata
        for field, usher_col in metadata.items():
            nextstrain_col = nextstrain_metadata[field]
            for strain in usher_col:
                genbank = try_parse_genbank(strain)
                if genbank:
                    value = nextstrain_col.get(genbank)
                    if value:
                        usher_col[strain] = value

    # Load usher tree.
    # Collect all names appearing in the usher tree.
    # if hasattr(args, 'refined_proto') and hasattr(args, 'refined_tree'):
    #     proto, tree = args.refined_proto, args.refined_tree
    # else:
    #     proto, tree = load_proto(args.tree_file_out)
    proto, tree = load_proto(args.tree_file_out)
    nuc_mutations_by_clade, proto, tree = load_mutation_tree(proto, tree)
    assert nuc_mutations_by_clade

    # Collect background mutation statistics.
    stats = defaultdict(Counter)
    aa_substitutions = stats["aaSubstitutions"]
    for mutations in nuc_mutations_by_clade.values():
        aa_substitutions.update(nuc_mutations_to_aa_mutations(mutations))

    # Collect condensed samples.
    condensed_nodes = {}
    for node in proto.condensed_nodes:
        condensed_nodes[node.node_name] = list(node.condensed_leaves)

    # Propagate fine clade names downward.
    node_to_clade = {}
    for clade, meta in zip(tree.find_clades(), proto.metadata):
        if meta.clade:
            node_to_clade[clade] = meta.clade
        # Propagate down to descendent clones.
        for child in clade.clades:
            node_to_clade[child] = node_to_clade[clade]

    # Collect info from each node in the tree.
    fields = "day", "location", "lineage"
    columns = defaultdict(list)
    sample_keys = set()
    skipped = stats["skipped"]
    skipped_by_day = Counter()
    nodename_to_count = Counter()
    for node, meta in zip(tree.find_clades(), proto.metadata):
        keys = condensed_nodes.get(node.name, [node.name])
        for key in keys:
            if key is None:
                continue
            sample_keys.add(key)

            if args.gisaid_metadata_file_in:
                key = try_parse_gisaid(key, public_to_gisaid)
                if key is None:
                    skipped["no gisaid id"] += 1
                    continue

            row = {k: metadata[k].get(key) for k in fields}
            if row["day"] is None:
                skipped["no date"] += 1
                continue
            if row["location"] is None:
                skipped["no location"] += 1
                skipped_by_day["no location", row["day"]] += 1
                continue

            columns["clade"].append(node_to_clade[node])
            columns["index"].append(key)
            for k, v in row.items():
                columns[k].append(v)
            nodename_to_count[node.name] += 1
    logger.info(f"Found {len(sample_keys)} samples in the usher tree")
    logger.info(f"Skipped {sum(skipped.values())} nodes because:\n{skipped}")
    columns = dict(columns)
    assert columns
    assert len(set(map(len, columns.values()))) == 1, "columns have unequal length"
    assert sum(skipped.values()) < args.max_skippage, f"suspicious skippage:\n{skipped}"
    logger.info(f"Kept {len(columns['day'])} rows")
    stats["skipped_by_day"] = skipped_by_day

    # Save columns.pkl in the new directory
    columns_path = os.path.join(args.output_dir, "columns.pkl")
    with open(columns_path, "wb") as f:
        pickle.dump(columns, f)
    logger.info(f"Saved {columns_path}")

    # Stats file is already using args.stats_file_out which was modified
    with open(args.stats_file_out, "wb") as f:
        pickle.dump(stats, f)
    logger.info(f"Saved {args.stats_file_out}")

    return columns, nodename_to_count


def prune_tree(args, max_num_clades, coarse_to_fine, nodename_to_count):
    proto, tree = load_proto(args.tree_file_out)

    # Add weights for leaves.
    cum_weights = {}
    mrca_weights = {}
    for clade in tree.find_clades():
        count = nodename_to_count[clade.name]
        cum_weights[clade] = count
        mrca_weights[clade] = count**2

    # Add weights of MRCA pairs.
    reverse_clades = list(tree.find_clades())
    reverse_clades.reverse()  # from leaves to root
    for parent in reverse_clades:
        for child in parent.clades:
            cum_weights[parent] += cum_weights[child]
        for child in parent.clades:
            mrca_weights[parent] += cum_weights[child] * (
                cum_weights[parent] - cum_weights[child]
            )
    num_samples = sum(nodename_to_count.values())
    assert cum_weights[tree.root] == num_samples
    assert sum(mrca_weights.values()) == num_samples**2

    # Aggregate among clones to basal representative.
    weights = defaultdict(float)
    for meta, parent in zip(reversed(proto.metadata), reverse_clades):
        for child in parent.clades:
            mrca_weights[parent] += mrca_weights.get(child, 0)
        assert isinstance(meta.clade, str)
        if meta.clade:
            weights[meta.clade] = mrca_weights.pop(parent)
    assert sum(weights.values()) == num_samples**2

    # To ensure pango lineages remain distinct, set their weights to infinity.
    for fine in coarse_to_fine.values():
        weights[fine] = args.current_pango_loss
    assert "" not in weights

     # Prune the tree, minimizing the number of incorrect mutations.
    pruned_tree_filename = os.path.join(args.output_dir, 
                                       f"lineageTree_{args.current_pango_loss}_{max_num_clades}.pb")
    
    meso_set, pruned_proto, pruned_tree = prune_mutation_tree(
        proto,
        tree,
        args.tree_file_out,
        max_num_clades,
        weights
    )
    
     # Save the pruned tree
    with open(pruned_tree_filename, "wb") as f:
        f.write(pruned_proto.SerializeToString())
    logger.info(f"Saved pruned tree to {pruned_tree_filename}")
    

    assert len(meso_set) == max_num_clades
    return FineToMeso(meso_set), (pruned_proto, pruned_tree)


def extract_features(
    args,
    max_num_clades,
    fine_to_coarse,
    coarse_to_fine,
    nodename_to_count,
    columns,
):
    logger.info(f"Extracting features with {max_num_clades} clades and pango loss {args.current_pango_loss}")

    # Prune tree, updating data structures to use meso-scale clades.
    fine_to_meso, (pruned_proto, pruned_tree)  = prune_tree(
        args, max_num_clades, coarse_to_fine, nodename_to_count
    )
    fine_to_coarse = {fine_to_meso(f): c for f, c in fine_to_coarse.items()}
    coarse_to_fine = {c: fine_to_meso(f) for c, f in coarse_to_fine.items()}

    # Save finer columns.
    columns = dict(columns)
    columns["clade"] = [fine_to_meso(f) for f in columns["clade"]]
    clade_set = set(columns["clade"])
    assert len(clade_set) <= max_num_clades

    columns_file_out = os.path.join(args.output_dir, 
                                   f"{args.current_pango_loss}_columns_{max_num_clades}.pkl")
    with open(columns_file_out, "wb") as f:
        pickle.dump(columns, f)
    logger.info(f"Saved {columns_file_out}")

    if not args.skip_features:
        # Convert from nucleotide mutations to amino acid mutations.
        nuc_mutations_by_clade = load_mutation_tree(pruned_proto, pruned_tree)[0]

        assert nuc_mutations_by_clade
        aa_mutations_by_clade = {
            clade: nuc_mutations_to_aa_mutations(mutations)
            for clade, mutations in nuc_mutations_by_clade.items()
        }

        # Create dense aa features.
        clades = sorted(nuc_mutations_by_clade)
        clade_ids = {k: i for i, k in enumerate(clades)}
        aa_mutations = sorted(set().union(*aa_mutations_by_clade.values()))
        logger.info(f"Found {len(aa_mutations)} amino acid mutations")
        mutation_ids = {k: i for i, k in enumerate(aa_mutations)}
        aa_features = torch.zeros(len(clade_ids), len(mutation_ids), dtype=torch.bool)
        for clade, ms in aa_mutations_by_clade.items():
            i = clade_ids[clade]
            for m in ms:
                j = mutation_ids.get(m)
                aa_features[i, j] = True

        # Save features.
        features = {
            "clades": clades,
            "clade_to_lineage": fine_to_coarse,
            "lineage_to_clade": coarse_to_fine,
            "aa_mutations": aa_mutations,
            "aa_features": aa_features,
        }
        features_file_out = os.path.join(args.output_dir, 
                                        f"{args.current_pango_loss}_features_{max_num_clades}.pt")
        torch.save(features, features_file_out)
        logger.info(f"Saved {features_file_out}")



def process_combination(args, max_num_clades, pango_loss, preprocessed_data):
    args.current_max_num_clades = max_num_clades
    args.current_pango_loss = pango_loss
    extract_features(
        args,
        max_num_clades,
        preprocessed_data['fine_to_coarse'],
        preprocessed_data['coarse_to_fine'],
        preprocessed_data['nodename_to_count'],
        preprocessed_data['columns'],
    )


def load_or_process_data(args):
    """Load preprocessed data from cache if available, otherwise process from scratch."""
    cache_file = 'preprocessed_data.pkl'
    
    # Try to load from cache first
    if os.path.exists(cache_file):
        logger.info(f"Loading preprocessed data from {cache_file}")
        try:
            with open(cache_file, 'rb') as f:
                preprocessed_data = pickle.load(f)
            
            # Validate that the cache contains all required keys
            required_keys = {'fine_to_coarse', 'coarse_to_fine', 'columns', 'nodename_to_count'}
            if all(key in preprocessed_data for key in required_keys):
                logger.info("Successfully loaded preprocessed data from cache")
                return preprocessed_data
            else:
                logger.warning("Cache file is missing required data, reprocessing...")
        except Exception as e:
            logger.warning(f"Failed to load cache: {e}, reprocessing...")
    
    # If cache doesn't exist or is invalid, process the data
    logger.info("Processing data from scratch...")
    
    # Extract mappings between coarse lineages and fine clades
    coarse_proto = args.tree_file_in
    fine_proto = args.tree_file_out
    fine_to_coarse, proto, tree = refine_mutation_tree(coarse_proto, fine_proto)

    # Create columns
    columns, nodename_to_count = load_metadata(args)
    columns["lineage"] = [fine_to_coarse[f] for f in columns["clade"]]

    # Choose the basal representative
    coarse_to_fines = defaultdict(list)
    for fine, coarse in fine_to_coarse.items():
        coarse_to_fines[coarse].append(fine)
    coarse_to_fine = {c: min(fs) for c, fs in coarse_to_fines.items()}

    # Create preprocessed data dictionary
    preprocessed_data = {
        'fine_to_coarse': fine_to_coarse,
        'coarse_to_fine': coarse_to_fine,
        'columns': columns,
        'nodename_to_count': nodename_to_count
    }

    # Save to cache
    logger.info(f"Saving preprocessed data to {cache_file}")
    with open(cache_file, 'wb') as f:
        pickle.dump(preprocessed_data, f)
    
    return preprocessed_data


def create_run_directory(args):
    """Create a directory name based on all parameters"""
    mnc_params = args.max_num_clades.replace(",", "_")
    pl_params = args.pango_loss.replace(",", "_")
    dir_name = f"MNC_{mnc_params}_PL_{pl_params}"
    full_path = os.path.join("results", dir_name)
    os.makedirs(full_path, exist_ok=True)
    logger.info(f"Created output directory: {full_path}")
    return full_path

def modify_save_paths(args):
    """Modify the default save paths to use the new directory"""
    output_dir = create_run_directory(args)
    
    # Modify the default paths in args
    args.stats_file_out = os.path.join(output_dir, "stats.pkl")
    args.output_dir = output_dir  # Add this to args for other functions to use
    
    return args




def main(args):
    # Modify save paths before any processing
    args = modify_save_paths(args)

    # Load or process data using cache
    preprocessed_data = load_or_process_data(args)

    # Create a list of all combinations
    combinations = [
        (int(max_num_clades), float(pango_loss))
        for max_num_clades in args.max_num_clades.split(",")
        for pango_loss in args.pango_loss.split(",")
    ]

    # Use multiprocessing to process combinations in parallel
    with multiprocessing.Pool(processes=args.num_cpus) as pool:
        pool.starmap(
            process_combination,
            [(args, max_num_clades, pango_loss, preprocessed_data) for max_num_clades, pango_loss in combinations]
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess pangolin mutations")
    parser.add_argument(
        "--usher-metadata-file-in", default="results/usher/metadata.tsv"
    )
    parser.add_argument(
        "--nextstrain-metadata-file-in", default="results/nextstrain/metadata.tsv"
    )
    parser.add_argument("--gisaid-metadata-file-in", default="")
    parser.add_argument("--tree-file-in", default="results/usher/all.masked.pb")
    parser.add_argument("--tree-file-out", default="results/lineageTree.fine.pb")
    parser.add_argument("--stats-file-out", default="results/stats.pkl")
    parser.add_argument("--recover-missing-usa-state", action="store_true")
    parser.add_argument("-s", "--max-skippage", type=float, default=1e7)
    parser.add_argument("-c", "--max-num-clades", default="2000,3000,5000,10000")
    parser.add_argument("--start-date", default=START_DATE)
    # The new arguments, added by Trojanskis et al.
    parser.add_argument('--pango-loss', type=str, default="0", help='Comma-separated list of weights for Pango lineage loss')
    parser.add_argument('--skip-features', action='store_true', help='Skip feature extraction')
    parser.add_argument('--num-cpus', type=int, default=1, help='Number of CPUs to use for parallel processing')


    args = parser.parse_args()
    args.start_date = try_parse_date(args.start_date)
    main(args)
