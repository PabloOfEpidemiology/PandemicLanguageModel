# Set-up ----
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  } else {
    library(package, character.only = TRUE)
  }
}

# List of packages to install and load
packages <- c("tidyverse", "data.table", "rbin", "cmdstanr", "patchwork", 
              "RColorBrewer", "nnet", "broom", "seqinr", "Biostrings", 
              "tictoc", "MASS", "foreach", "parallel", "doParallel", "dplyr", "ggplot2", "lubridate")

# Install and load each package
invisible(sapply(packages, install_and_load))

# Set global bitmap type to cairo
options(bitmapType = "cairo")

print("All libraries have been set up")

# Load necessary packages for parallelization
library(doParallel)
library(foreach)


# Initialize GLOBAL to FALSE by default (equivalent to --countries)
GLOBAL <- FALSE  # Set to TRUE if --global flag is provided

# Initialize num.cores to 1 by default
num.cores <- 1


# Parse command-line arguments ----
args <- commandArgs(trailingOnly=TRUE)
count.mut_thresholds <- c(1500, 3000, 7000, 10000, 20000, 40000, 70000, 100000)  # default values

# Determine the number of tasks based on count.mut_thresholds
num_tasks <- length(count.mut_thresholds)

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--count.mut" && i + 1 <= length(args)) {
      count_mut_values <- args[i + 1]
      count.mut_thresholds <- as.numeric(unlist(strsplit(count_mut_values, ",")))
    }
    # Check for the --global flag
    if (args[i] == "--global") {
      GLOBAL <- TRUE
    }
    # Optionally, handle --countries explicitly if needed
    if (args[i] == "--countries") {
      GLOBAL <- FALSE
    }
    # Add parsing for --num.cores
    if (args[i] == "--num.cores" && i + 1 <= length(args)) {
      num.cores <- as.integer(args[i + 1])
    }
  }
}


print(paste("Using count.mut_threshold values:", paste(count.mut_thresholds, collapse=", ")))
print(paste("GLOBAL set to:", GLOBAL))
print(paste("Number of cores set to:", num.cores))

CLUSTER <- TRUE  # Set to FALSE to use test inputs
count.country_threshold <- 100000  # Set to 100,000 by default
bin.size <- 1
date.max <- "2024-10-28"  # Changed to 28th October 2024

# Define the output prefix based on the GLOBAL switch
output_prefix_base <- if (GLOBAL) {
  "Spike_Global"
} else {
  "Spike_Countries"
}

# Set the working directory to where the data files are located ----
setwd("/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing")
print(paste("Working directory set to:", getwd()))

# Define input directory ----
input_dir_spike_only <- "/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing/Input_data/Spike_only/"

input_dir_raw_gisaid <- "/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing/Input_data/raw_GISAID_data/"

# Define output directory ----
output_dir <- "Simulations/Spike"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
print(paste("Output directory created/confirmed:", output_dir))

print("Skipping cache loading and processing data from scratch.")

# Load data files before the loop ----
print("Loading data files before the loop.")


# Define UK countries vector
UK.v <- c("England", "Scotland", "Wales")

# Load and preprocess the metadata
metadata.original.name <- file.path(input_dir_raw_gisaid, "metadata.tsv")
metadata.original <- fread(metadata.original.name, header=TRUE, sep="\t", quote="", check.names=TRUE)
print("Original metadata loaded.")
print(paste("Loaded original metadata with", nrow(metadata.original), "rows"))

# Filter for human hosts
metadata.original <- metadata.original %>% filter(Host == "Human")
print(paste("Filtered original metadata for human hosts:", nrow(metadata.original), "rows remaining"))

# Separate virus names
metadata.original <- metadata.original %>% mutate(Virus.name2 = Virus.name) %>% separate(Virus.name2, c("virus", "country", "name", "year"), sep = "/")



# Load additional metadata
metadata.name <- file.path(input_dir_spike_only, "Spike_full.nextclade.tsv")
metadata <- fread(metadata.name, header=TRUE, sep="\t", quote="", check.names=TRUE)
print("Additional metadata loaded.")
print(paste("Loaded nextclade metadata with", nrow(metadata), "rows"))

# Load and preprocess N count data
data.count_N.name <- file.path(input_dir_spike_only, "gene_S.translation.count_table.txt")
data.count_N <- fread(data.count_N.name, header=TRUE, sep="\t", quote="", check.names=TRUE)
print(paste("Loaded gene count data with", nrow(data.count_N), "rows"))


print("Data files loaded.")



# Filter translation data
data.count_N.filtered <- data.count_N %>% 
  dplyr::rename(seqName = name) %>% 
  mutate(prop.N = count_N / length) %>% 
  filter(prop.N < 0.01, count_stop == 1)
print(paste("Filtered gene count data to", nrow(data.count_N.filtered), "rows"))

# Merge with metadata
metadata.merged <- metadata %>% 
  inner_join(data.count_N.filtered, by = "seqName")
metadata.merged <- metadata.merged %>% distinct(seqName, .keep_all = TRUE)
print(paste("Merged metadata and gene count data to", nrow(metadata.merged), "rows"))

# Additional processing of metadata
metadata.merged <- metadata.merged %>% 
  mutate(seqName2 = seqName) %>% 
  separate(seqName2, c("virus", "country", "name", "time"), sep = "/")

metadata.merged <- metadata.merged %>% 
  mutate(time2 = time) %>% 
  separate(time2, c("year", "Collection.date", "Submission.date"), sep = "\\|")

metadata.merged <- metadata.merged %>% 
  filter(!is.na(Collection.date), 
         str_length(Collection.date) == 10, 
         virus == "hCoV-19")
metadata.merged <- metadata.merged %>% 
  mutate(Collection.date = as.Date(Collection.date)) %>% 
  filter(Collection.date <= as.Date(date.max))

metadata.merged <- metadata.merged %>% 
  inner_join(metadata.original %>% dplyr::select(name), by = "name")


print(paste("Further processed metadata to", nrow(metadata.merged), "rows"))

# Group by country and summarize counts
metadata.merged <- metadata.merged %>% 
  mutate(country = ifelse(country %in% UK.v, "UK", as.character(country)))

# Set 'country' to 'global' if GLOBAL is TRUE
if (GLOBAL) {
  metadata.merged <- metadata.merged %>% mutate(country = "global")
}

# Proceed with grouping and summarizing
country.interest.df <- metadata.merged %>% 
  group_by(country) %>%
  summarize(count.country = n()) %>%
  arrange(desc(count.country)) %>%
  filter(count.country > count.country_threshold)
country.interest.v <- country.interest.df$country
print("Countries with significant data:")
print(country.interest.df)

# Filtering metadata based on country of interest.
print("Filtering metadata based on country of interest.")
metadata.merged <- metadata.merged %>% filter(country %in% country.interest.df$country)
print(paste("Number of entries after filtering metadata:", nrow(metadata.merged)))


# Load mutation information and display a summary
Spike_mut.info.name <- file.path(input_dir_spike_only, "Spike_full.nextclade.mut_long.tsv")
Spike_mut.info_pre_loop <- fread(Spike_mut.info.name, header = TRUE, sep = "\t", quote = "", check.names = TRUE)
print("Loaded mutation information. First few rows:")
print(head(Spike_mut.info_pre_loop))

# Filter mutations with "S:"
Spike_mut.info_pre_loop <- Spike_mut.info_pre_loop %>% filter(grepl("S:", mut))
print("Filtered for 'S:' mutations. First few rows:")
print(head(Spike_mut.info_pre_loop))


print("Data pre-processing completed outside of the loop.")

# Set up the parallel backend using the specified number of cores
registerDoParallel(cores = num.cores)

# Start parallel processing for selected count.mut_threshold values ----
print("Starting parallel processing for selected count.mut_threshold values")

# Use foreach to run simulations in parallel
summary_stats_list <- foreach(count.mut_threshold = count.mut_thresholds, .packages = c("tidyverse", "data.table", "lubridate")) %dopar% {
  
  print(paste("Starting iteration with count.mut_threshold =", count.mut_threshold))
  
  # Define the output prefix based on the GLOBAL switch and count.mut_thresholds
  output_prefix <- paste0(output_prefix_base, "_", paste(count.mut_thresholds, collapse="_"))
  
  # Pre-processing data ----
  print("Starting data pre-processing inside the loop...")
  
  # Note: data.count_N.filtered, metadata.merged, Spike_mut.info, etc., are already loaded and preprocessed
  
  # Join with metadata.merged
  Spike_mut.info <- Spike_mut.info_pre_loop %>% right_join(metadata.merged %>% dplyr::select(seqName), by = "seqName")
  print("Joined mutation info with filtered metadata. Rows in merged data:")
  print(nrow(Spike_mut.info))
  
  # Count mutations and filter based on count.mut_threshold
  Spike_count.df.mut <- Spike_mut.info %>% 
    group_by(mut) %>% 
    summarize(count.mut = n())
  Spike_count.df.mut <- Spike_count.df.mut %>% 
    arrange(desc(count.mut)) %>% 
    filter(count.mut > count.mut_threshold)
  print(paste("Mutations with counts >", count.mut_threshold, ":"))
  print(Spike_count.df.mut)
  
  # Further refine the mutation information
  Spike_mut.info.iter <- Spike_mut.info %>% 
    inner_join(Spike_count.df.mut %>% dplyr::select(mut), by = "mut")
  Spike_mut.info.iter <- Spike_mut.info.iter %>% unique()
  Spike_mut.info.iter <- Spike_mut.info.iter %>% mutate(mut = as.factor(mut))
  print(paste("Refined mutation information. Unique mutations count:", length(unique(Spike_mut.info.iter$mut))))
  
  # Continue with the rest of the processing inside the loop
  # Summarize haplotype information
  print("Summarizing haplotype information.")
  
  # Convert mut.info to a data.table
  Spike_mut.info.dt <- as.data.table(Spike_mut.info.iter)
  # Summarize haplotype information with the correct delimiter
  Spike_mut.info.sum.dt <- Spike_mut.info.dt[, .(haplotype = paste(mut, collapse = ", ")), by = seqName]
  # Convert back to tibble if needed
  Spike_mut.info.sum <- as_tibble(Spike_mut.info.sum.dt)
  print("Haplotype information summarized using data.table with correct delimiter.")
  
  # Create haplotype IDs and join with metadata
  Spike_mut.info.sum <- Spike_mut.info.sum %>% 
    mutate(hap_Id = paste("hap_", as.character(as.numeric(as.factor(haplotype))), sep=""))
  Spike_mut.info.sum <- Spike_mut.info.sum %>% 
    inner_join(metadata.merged %>% dplyr::select(seqName, country), by="seqName")
  print("Haplotype information summarized and joined with metadata.")
  print(head(Spike_mut.info.sum))
  
  # Unique haplotype information
  Spike_hap.info.df <- Spike_mut.info.sum %>% 
    ungroup() %>% 
    dplyr::select(hap_Id, haplotype) %>% 
    unique()
  print("Unique haplotypes information:")
  print(Spike_hap.info.df)
  
  # Group by haplotype and country, filter by count >= 20
  Spike_count.df.hap <- Spike_mut.info.sum %>% group_by(hap_Id, country) %>% summarize(count.hap = n()) %>% arrange(desc(count.hap))
  Spike_count.df.hap <- Spike_count.df.hap %>% filter(count.hap >= 20)
  print("Haplotypes with counts >= 20:")
  print(Spike_count.df.hap)
  
  # Filter for haplotypes in count.df.hap
  Spike_mut.info.sum <- Spike_mut.info.sum %>% filter(hap_Id %in% as.character(Spike_count.df.hap$hap_Id))
  
  # Merge haplotypes back with metadata
  metadata.merged.iter <- metadata.merged %>% 
    inner_join(Spike_mut.info.sum %>% dplyr::select(seqName, hap_Id), by = "seqName")
  metadata.merged.iter <- metadata.merged.iter %>% distinct(seqName, .keep_all = TRUE)
  metadata.merged.iter <- metadata.merged.iter %>% filter(!is.na(Collection.date))
  print(paste("Metadata merged with haplotype information. Rows remaining:", nrow(metadata.merged.iter)))
  
  # Further filter metadata
  Spike_metadata.filtered.interest <- metadata.merged.iter %>%
    mutate(Collection.date = as.Date(Collection.date),
           date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date)) + 1,
           date.bin = cut(date.num, seq(0, max(date.num), bin.size)),
           date.bin.num = as.numeric(date.bin))
  Spike_metadata.filtered.interest <- Spike_metadata.filtered.interest %>% filter(!is.na(date.bin.num))
  print(paste("Filtered metadata with date bins. Rows remaining:", nrow(Spike_metadata.filtered.interest)))
  
  print("Data pre-processing completed.")
  
  # Create Spike_data.coef.df with dummy Re values, including filtering step
  print("Creating Spike_data.coef.df with dummy Re values, including filtering hap_Id-country combinations with counts >= 20.")
  
  # Initialize an empty list to store data frames for each country
  list_of_combinations <- list()
  
  # Loop over each country to filter hap_Id-country combinations with counts >= 20
  for (country.interest in country.interest.df$country) {
    print(paste("Processing country:", country.interest))
    
    # Filter data for the current country
    Spike_metadata_country <- Spike_metadata.filtered.interest %>%
      filter(country == country.interest)
    
    if (nrow(Spike_metadata_country) == 0) {
      print(paste("No data available for country:", country.interest))
      next
    }
    
    # Group by hap_Id and count occurrences
    hap_counts <- Spike_metadata_country %>%
      group_by(hap_Id) %>%
      summarize(total = n()) %>%
      filter(total >= 20)
    
    if (nrow(hap_counts) == 0) {
      print(paste("No haplotypes with counts >= 20 for country:", country.interest))
      next
    }
    
    # Keep only hap_Id with counts >= 20
    Spike_metadata_country_filtered <- Spike_metadata_country %>%
      filter(hap_Id %in% hap_counts$hap_Id)
    
    # Get unique combinations of hap_Id and country
    unique_combinations_country <- Spike_metadata_country_filtered %>%
      dplyr::select(hap_Id, country) %>%
      unique()
    
    # Append to the list
    list_of_combinations[[country.interest]] <- unique_combinations_country
  }
  
  # Combine all unique combinations
  unique_combinations <- bind_rows(list_of_combinations)
  
  # Create hap_Id_country combinations in Spike_metadata.filtered.interest
  Spike_metadata.filtered.interest <- Spike_metadata.filtered.interest %>%
    mutate(hap_Id_country = paste(hap_Id, country, sep="_"))
  
  # Assign dummy Re values between 0 and 2
  set.seed(123)  # For reproducibility
  unique_combinations$relative_Re <- runif(nrow(unique_combinations), min=0, max=2)
  
  Spike_data.coef.df <- unique_combinations
  
  print("Dummy Re values assigned to hap_Id-country combinations.")
  print("Filtered hap_Id-country combinations:")
  print(Spike_data.coef.df)
  
  # Post-processing ----
  print("Starting post-processing.")
  
  Spike_mut.info <- Spike_mut.info_pre_loop %>% 
    right_join(Spike_metadata.filtered.interest %>% dplyr::select(seqName, hap_Id), by = "seqName")
  print("Joined mutation info with metadata.")
  
  Spike_mut.info <- Spike_mut.info %>% unique()
  Spike_mut.info <- Spike_mut.info %>% mutate(mut = as.factor(mut))
  print("Duplicates removed and mutations converted to factor.")
  
  Spike_mut.info.dt <- as.data.table(Spike_mut.info)
  Spike_mut.info.sum.all_S <- Spike_mut.info.dt[, .(haplotype.all = toString(mut)), by = seqName]
  print("Haplotype information summarized.")
  
  Spike_mut.info.sum.all_S <- Spike_mut.info.sum.all_S[, hap_Id.all_mut := paste("hap_all_", as.character(as.numeric(as.factor(haplotype.all))), sep="")]
  print("Haplotype ID created.")
  
  Spike_mut.info.sum.all_S <- as_tibble(Spike_mut.info.sum.all_S)
  print("Converting Spike_mut.info.sum.all_S back to tibble.")
  
  Spike_metadata.filtered.interest <- Spike_metadata.filtered.interest %>% 
    inner_join(Spike_mut.info.sum.all_S %>% dplyr::select(seqName, hap_Id.all_mut), by="seqName")
  print("Joined summarized haplotype info with metadata.")
  
  date.quantile.df <- Spike_metadata.filtered.interest %>% 
    mutate(Collection.date = as.Date(Collection.date)) %>% 
    group_by(hap_Id) %>% 
    summarize(date.first = quantile(Collection.date, 0.01, type=1))
  date.quantile.df <- date.quantile.df %>% arrange(date.first)
  print("First collection date quantile calculated.")
  
  # Determine major haplotype
  all_S.major.hap.df <- Spike_metadata.filtered.interest %>%
    dplyr::group_by(hap_Id, hap_Id.all_mut) %>%
    dplyr::summarize(count = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(hap_Id) %>%
    dplyr::arrange(desc(count)) %>%
    dplyr::slice(1)  # Select the first (most frequent) haplotype
  
  print("Major haplotype determined.")
  print(all_S.major.hap.df)
  
  Spike_metadata.filtered.interest.representative <- Spike_metadata.filtered.interest %>% 
    inner_join(all_S.major.hap.df %>% dplyr::select(hap_Id, hap_Id.all_mut), by=c("hap_Id", "hap_Id.all_mut")) %>% 
    group_by(hap_Id) %>% 
    sample_n(1)
  
  Spike_metadata.filtered.interest.representative <- Spike_metadata.filtered.interest.representative %>% 
    dplyr::select(seqName, clade, Nextclade_pango, hap_Id)
  
  Spike_metadata.filtered.interest.representative <- Spike_metadata.filtered.interest.representative %>% 
    inner_join(Spike_hap.info.df, by="hap_Id")
  print("Representative sequences selected.")
  
  Spike_data.coef.df.merged <- Spike_data.coef.df %>% 
    inner_join(Spike_metadata.filtered.interest.representative, by="hap_Id")
  Spike_data.coef.df.merged <- Spike_data.coef.df.merged %>% 
    arrange(desc(relative_Re))
  print("Coefficient data merged with representative metadata.")
  
  Spike_data.coef.df.merged <- Spike_data.coef.df.merged %>% 
    left_join(date.quantile.df, by="hap_Id")
  print("Joined with date quantile data.")
  
  # Add hap_Id_country column
  Spike_data.coef.df.merged <- Spike_data.coef.df.merged %>%
    mutate(hap_Id_country = paste(hap_Id, country, sep="_"))
  
  print("Post-processing completed.")
  
  # Generate summary statistics ----
  
  print("Generating summary statistics.")
  
  # 1) Number of unique haplotype names
  num_unique_haplotypes <- length(unique(Spike_data.coef.df.merged$hap_Id))
  print(paste("Number of unique haplotype names:", num_unique_haplotypes))
  
  # 1.2) Number of mutations that pass the count.mut_threshold
  num_mutations_pass_threshold <- nrow(Spike_count.df.mut)
  print(paste("Number of mutations that pass the count.mut_threshold:", num_mutations_pass_threshold))
  
  # 2) Average haplotype size and standard deviation
  # Modified to include only haplotype-country combinations with counts >= 20
  haplotype_sizes <- Spike_metadata.filtered.interest %>%
    mutate(hap_Id_country = paste(hap_Id, country, sep="_")) %>%
    filter(hap_Id_country %in% Spike_data.coef.df.merged$hap_Id_country) %>%
    group_by(hap_Id) %>%
    summarize(count = n())
  
  average_haplotype_size <- mean(haplotype_sizes$count)
  std_haplotype_size <- sd(haplotype_sizes$count)
  
  print(paste("Average haplotype size:", average_haplotype_size))
  print(paste("Standard deviation of haplotype size:", std_haplotype_size))
  
  # 4) Number of unique hap_Id-country combinations
  num_unique_hap_country_combinations <- length(unique(Spike_data.coef.df.merged$hap_Id_country))
  print(paste("Number of unique hap_Id-country combinations:", num_unique_hap_country_combinations))
  
  # 5) Average size of unique hap_Id-country combinations and standard deviation
  hap_country_sizes <- Spike_metadata.filtered.interest %>%
    mutate(hap_Id_country = paste(hap_Id, country, sep="_")) %>%
    filter(hap_Id_country %in% Spike_data.coef.df.merged$hap_Id_country) %>%
    group_by(hap_Id, country) %>%
    summarize(count = n())
  
  average_hap_country_size <- mean(hap_country_sizes$count)
  std_hap_country_size <- sd(hap_country_sizes$count)
  
  print(paste("Average size of unique hap_Id-country combinations:", average_hap_country_size))
  print(paste("Standard deviation of unique hap_Id-country combinations:", std_hap_country_size))
  
  # 6) Calculate final number of sequences in the dataset
  total_sequences_in_final_analysis <- num_unique_haplotypes * average_haplotype_size
  
  # Save the summary statistics for this run
  summary_stats <- data.frame(
    count.mut_threshold = count.mut_threshold,
    num_unique_haplotypes = num_unique_haplotypes,
    num_mutations_pass_threshold = num_mutations_pass_threshold,
    average_haplotype_size = average_haplotype_size,
    std_haplotype_size = std_haplotype_size,
    num_unique_hap_country_combinations = num_unique_hap_country_combinations,
    average_hap_country_size = average_hap_country_size,
    std_hap_country_size = std_hap_country_size,
    total_sequences_in_final_analysis = total_sequences_in_final_analysis
  )
  
  print("Summary statistics for this run:")
  print(summary_stats)
  
  # Return the summary statistics
  return(summary_stats)
}

# After the loop, combine all summary statistics into one dataframe
summary_stats_df <- do.call(rbind, summary_stats_list)

print("Combined summary statistics:")
print(summary_stats_df)

# Write the combined summary statistics to simulations.csv in the output directory
output_file <- file.path(output_dir, paste0(output_prefix_base, "_", paste(count.mut_thresholds, collapse="_"), "_simulations.csv"))
write.csv(summary_stats_df, file = output_file, row.names = FALSE)

print(paste("Summary statistics saved to", output_file))

print("Script execution completed.")