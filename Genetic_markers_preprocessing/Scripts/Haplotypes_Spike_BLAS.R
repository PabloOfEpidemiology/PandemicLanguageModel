#!/usr/bin/env R

# Set-up ----
install_and_load <- function(package) {
  if (!require(package, character.only = TRUE)) {
    install.packages(package, repos = "https://cran.r-project.org")
    library(package, character.only = TRUE)
  } else {
    library(package, character.only = TRUE)
  }
}

# List of packages to install and load
packages <- c("tidyverse", "data.table", "rbin", "cmdstanr", "patchwork", 
              "RColorBrewer", "nnet", "broom", "seqinr", "Biostrings", 
              "tictoc", "MASS", "foreach", "parallel", "doParallel", "dplyr")

# Install and load each package
invisible(sapply(packages, install_and_load))

# Check the version of 'dplyr' and update if necessary
if (packageVersion("dplyr") < "0.8.0") {
  install.packages("dplyr", repos = "https://cran.r-project.org")
  library(dplyr)
}




# Verify BLAS configuration
cat("R is using the following BLAS:\n")
print(sessionInfo())

# Alternatively, check the BLAS library directly
cat("Linking to BLAS:\n")
print(R.version$BLAS)


# Command line arguments parsing ----
args <- commandArgs(trailingOnly = TRUE)

# Set default values
CLUSTER <- TRUE  # Set to FALSE to use test inputs
GLOBAL <- FALSE  # Default is FALSE, can be overridden by CLI
count.mut_threshold <- 7000
count.country_threshold <- 100000
bin.size <- 1

# Parse command line arguments
if ("--Global" %in% args) {
  GLOBAL <- TRUE
}

if ("--Countries" %in% args) {
  GLOBAL <- FALSE
}

threshold_idx <- which(args == "--count.mut_threshold")
if (length(threshold_idx) > 0 && length(args) >= threshold_idx + 1) {
  count.mut_threshold <- as.numeric(args[threshold_idx +1])
}

# Define the output prefix based on the GLOBAL switch and count.mut_threshold
output_prefix <- if (GLOBAL) {
  paste0("Global_", count.mut_threshold)
} else {
  paste0("Countries_", count.mut_threshold)
}

# Number of cores to use ----
num_cores_available <- detectCores()
num_cores_to_use <- min(num_cores_available - 6, 1)
num_cores_to_use <- max(1, num_cores_to_use)
print(paste("Number of available cores:", num_cores_available))
print(paste("Number of cores that will be used for parallel processing:", num_cores_to_use))

# Transmissibility parameters ----
generation_time <- 2.1
date.max <- "2024-10-28"
ref_hap_clade <- "22B"

# Set the working directory ----
if (CLUSTER) {
  setwd("/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing/")
} else {
  setwd("/Users/user/OneDrive/CovidGPT/Jumpei_multinomial_model")
}
print(paste("Working directory set to:", getwd()))

# Define input directories ----
input_dir_spike_only <- if (CLUSTER) {
  "/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing/Input_data/Spike_only/"
} else {
  "test_Inputs/"
}

input_dir_raw_gisaid <- if (CLUSTER) {
  "/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing/Input_data/raw_GISAID_data/"
} else {
  "test_Inputs/"
}

# Define output directory ----
output_base_dir <- if (CLUSTER) {
  "/data/zool-paleovirology/ball6168/PandemicLanguageModel/Genetic_markers_preprocessing/Haplotypes"
} else {
  "test_Outputs"
}

output_subdir <- paste0("Spike_", output_prefix)
output_dir <- file.path(output_base_dir, output_subdir)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
print(paste("Output directory created/confirmed:", output_dir))

# Caching code ----

# Define cache directory
cache_base_dir <- if (CLUSTER) {
  file.path(getwd(), "cache")
} else {
  file.path(getwd(), "test_cache")
}

cache_subdir <- "Spike"
cache_dir <- file.path(cache_base_dir, cache_subdir)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

# Function to check if cached data exists and is recent enough
is_cache_valid <- function(cache_file, max_age_days = 40) {
  if (file.exists(cache_file)) {
    file_age <- difftime(Sys.time(), file.mtime(cache_file), units = "days")
    return(file_age <= max_age_days)
  }
  return(FALSE)
}

# Function to get cache file path
get_cache_file <- function(threshold, GLOBAL) {
  global_str <- ifelse(GLOBAL, "_global", "_countries")
  file.path(cache_dir, paste0("Spike_preprocessed_data_cache_", threshold, global_str, ".rds"))
}

cache_file <- get_cache_file(count.mut_threshold, GLOBAL)

# Check if valid cache exists for the current count.mut_threshold value
process_from_scratch <- TRUE

if (is_cache_valid(cache_file)) {
  print(paste("Loading cached preprocessed data for count.mut_threshold =", count.mut_threshold))
  load(cache_file)
  
  # Verify loaded objects
  required_vars <- c("Spike_metadata.filtered.interest", "Spike_count.hap.mat", "Spike_count.df.hap2", 
                     "Spike_count.df.hap2.min", "Spike_hap_Id.ref", "Spike_count.df.mut", 
                     "country.interest.df", "bin.size", "generation_time",
                     "country.interest.v", "Spike_mut.info.name", "Spike_hap.info.df","Spike_count.df.hap")
  
  missing_vars <- required_vars[!sapply(required_vars, exists)]
  if (length(missing_vars) > 0) {
    warning(paste("Missing variables in cache:", paste(missing_vars, collapse = ", ")))
    print("Processing data from scratch due to missing variables in cache.")
    process_from_scratch <- TRUE
  } else {
    print("All required variables loaded from cache.")
    process_from_scratch <- FALSE
  }
} else {
  print(paste("Processing data from scratch for count.mut_threshold =", count.mut_threshold))
}

if (process_from_scratch) {
  # Pre-processing data ----
  
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
  
  # Filter translation data
  data.count_N.filtered <- data.count_N %>% dplyr::rename(seqName = name) %>% mutate(prop.N = count_N / length) %>% filter(prop.N < 0.01, count_stop == 1)
  print(paste("Filtered gene count data to", nrow(data.count_N.filtered), "rows"))
  
  # Merge with metadata
  metadata.merged <- metadata %>% inner_join(data.count_N.filtered, by = "seqName")
  metadata.merged <- metadata.merged %>% distinct(seqName, .keep_all = TRUE)
  print(paste("Merged metadata and gene count data to", nrow(metadata.merged), "rows"))
  
  # Additional processing of metadata
  metadata.merged <- metadata.merged %>% mutate(seqName2 = seqName) %>% separate(seqName2, c("virus", "country", "name", "time"), sep = "/")
  metadata.merged <- metadata.merged %>% mutate(time2 = time) %>% separate(time2, c("year", "Collection.date", "Submission.date"), sep = "\\|")
  metadata.merged <- metadata.merged %>% filter(!is.na(Collection.date), str_length(Collection.date) == 10, virus == "hCoV-19")
  metadata.merged <- metadata.merged %>% mutate(Collection.date = as.Date(Collection.date)) %>% filter(Collection.date <= as.Date(date.max))
  metadata.merged <- metadata.merged %>% 
    inner_join(metadata.original %>% dplyr::select(name, Accession.ID), by = "name")
  print(paste("Further processed metadata to", nrow(metadata.merged), "rows"))
  
  
  print(colnames(metadata.original))
  print(colnames(metadata.merged))
  
  
  # Group by country and summarize counts
  metadata.merged <- metadata.merged %>% mutate(country = ifelse(country %in% UK.v, "UK", as.character(country)))
  
  # Set 'country' to 'global' if GLOBAL is TRUE
  if (GLOBAL) {
    metadata.merged <- metadata.merged %>% mutate(country = "global")
  }
  
  # Proceed with grouping and summarizing
  country.interest.df <- metadata.merged %>% group_by(country) %>%
    summarize(count.country = n()) %>%
    arrange(desc(count.country)) %>%
    filter(count.country > count.country_threshold)
  country.interest.v <- country.interest.df$country
  print("Countries with significant data:")
  print(country.interest.df)
  
  # Filtering metadata based on country of interest
  print("Filtering metadata based on country of interest.")
  metadata.merged <- metadata.merged %>% filter(country %in% country.interest.df$country)
  print(paste("Number of entries after filtering metadata:", nrow(metadata.merged)))
  
  # Load mutation information and display a summary
  Spike_mut.info.name <- file.path(input_dir_spike_only, "Spike_full.nextclade.mut_long.tsv")
  Spike_mut.info <- fread(Spike_mut.info.name, header = TRUE, sep = "\t", quote = "", check.names = TRUE)
  print("Loaded mutation information. First few rows:")
  print(head(Spike_mut.info))
  
  # Filter mutations with "S:"
  Spike_mut.info <- Spike_mut.info %>% filter(grepl("S:", mut))
  print("Filtered for 'S:' mutations. First few rows:")
  print(head(Spike_mut.info))
  
  # Join with metadata.merged
  Spike_mut.info <- Spike_mut.info %>% right_join(metadata.merged %>% dplyr::select(seqName), by = "seqName")
  print("Joined mutation info with filtered metadata. Rows in merged data:")
  print(nrow(Spike_mut.info))
  
  # Count mutations and filter based on count
  Spike_count.df.mut <- Spike_mut.info %>% group_by(mut) %>% summarize(count.mut = n())
  Spike_count.df.mut <- Spike_count.df.mut %>% arrange(desc(count.mut)) %>% filter(count.mut > count.mut_threshold)
  print(paste("Mutations with counts >", count.mut_threshold, ":"))
  print(Spike_count.df.mut)
  
  # Further refine the mutation information
  Spike_mut.info <- Spike_mut.info %>% inner_join(Spike_count.df.mut %>% dplyr::select(mut), by = "mut")
  Spike_mut.info <- Spike_mut.info %>% unique()
  Spike_mut.info <- Spike_mut.info %>% mutate(mut = as.factor(mut))
  print("Refined mutation information. Unique mutations count:")
  print(length(unique(Spike_mut.info$mut)))
  
  # Summarize haplotype information
  
  # Convert mut.info to a data.table
  Spike_mut.info.dt <- as.data.table(Spike_mut.info)
  # Summarize haplotype information with the correct delimiter
  Spike_mut.info.sum.dt <- Spike_mut.info.dt[, .(haplotype = paste(mut, collapse = ", ")), by = seqName]
  # Convert back to tibble if needed
  Spike_mut.info.sum <- as_tibble(Spike_mut.info.sum.dt)
  print("Haplotype information summarized using data.table with correct delimiter.")
  
  Spike_mut.info.sum <- Spike_mut.info.sum %>% mutate(hap_Id = paste("hap_", as.character(as.numeric(as.factor(haplotype))), sep=""))
  Spike_mut.info.sum <- Spike_mut.info.sum %>% inner_join(metadata.merged %>% dplyr::select(seqName, country), by="seqName")
  print("Haplotype information summarized and joined with metadata.")
  print(head(Spike_mut.info.sum))
  
  # Unique haplotype information
  Spike_hap.info.df <- Spike_mut.info.sum %>% ungroup() %>% dplyr::select(hap_Id, haplotype) %>% unique()
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
  metadata.merged <- metadata.merged %>% inner_join(Spike_mut.info.sum %>% dplyr::select(seqName, hap_Id), by = "seqName")
  metadata.merged <- metadata.merged %>% distinct(seqName, .keep_all = TRUE)
  metadata.merged <- metadata.merged %>% filter(!is.na(Collection.date))
  print("Metadata merged with haplotype information. Rows remaining:")
  print(nrow(metadata.merged))
  
  # Further filter metadata
  Spike_metadata.filtered.interest <- metadata.merged %>%
    mutate(Collection.date = as.Date(Collection.date),
           date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date)) + 1,
           date.bin = cut(date.num, seq(0, max(date.num), bin.size)),
           date.bin.num = as.numeric(date.bin))
  Spike_metadata.filtered.interest <- Spike_metadata.filtered.interest %>% filter(!is.na(date.bin.num))
  print("Filtered metadata with date bins. Rows remaining:")
  print(nrow(Spike_metadata.filtered.interest))
  
  # Create count matrix for haplotypes
  Spike_count.hap.mat <- Spike_metadata.filtered.interest %>%
    filter(clade == ref_hap_clade) %>%
    group_by(country, hap_Id) %>%
    summarize(count.hap = n()) %>%
    spread(key = country, value = count.hap) %>% na.omit()
  print("Count matrix for haplotypes (after filtering by clade):")
  print(Spike_count.hap.mat)
  
  # Gather data for haplotype count comparisons
  Spike_count.df.hap2 <- Spike_count.hap.mat %>% gather(key = country, value = count, -hap_Id)
  Spike_count.df.hap2.min <- Spike_count.df.hap2 %>% group_by(hap_Id) %>% summarize(count.min = min(count)) %>% arrange(desc(count.min))
  Spike_hap_Id.ref <- Spike_count.df.hap2.min$hap_Id[1]
  print(paste("Reference haplotype ID:", Spike_hap_Id.ref))
  print("Count matrix for reference haplotype:")
  print(Spike_count.hap.mat %>% filter(hap_Id == Spike_hap_Id.ref))
  
  # Save preprocessed data to cache
  print(paste("Saving preprocessed data to cache for count.mut_threshold =", count.mut_threshold))
  save(list = c("Spike_metadata.filtered.interest", "Spike_count.hap.mat", "Spike_count.df.hap2", 
                "Spike_count.df.hap2.min", "Spike_hap_Id.ref", "Spike_count.df.mut", 
                "country.interest.df", "bin.size", "generation_time",
                "country.interest.v", "Spike_mut.info.name", "Spike_hap.info.df","Spike_count.df.hap"), 
       file = cache_file, envir = .GlobalEnv)
}

print("Preprocessed data ready for parallelization.")

# Cleanup
print("Cleaning up memory.")
rm(metadata.original); gc(); gc()
rm(metadata); gc(); gc()
rm(Spike_mut.info); gc(); gc()

# Register the parallel backend to use many processors
print("Starting parallel processing...")

registerDoParallel(cores=num_cores_to_use)

Spike_data.coef.df <- data.frame()
# Convert to data.table for faster processing
setDT(Spike_metadata.filtered.interest)

# Extract the desired columns: country, date.bin.num, hap_Id
Spike_metadata_to_save <- Spike_metadata.filtered.interest[, .(country, date.bin.num, hap_Id)]

# MULTINOMIAL MODEL FITTING ----
print(paste("Number of countries to process:", length(country.interest.df$country)))
print("Countries to process:")
print(country.interest.df$country)

registerDoParallel(cores = num_cores_to_use)
print(paste("Registered parallel backend workers:", getDoParWorkers()))

results <- foreach(country.interest = country.interest.df$country,
                   .combine = rbind,
                   .packages = c("data.table", "nnet"),
                   .export = c("country.interest.df", "Spike_metadata.filtered.interest", "Spike_hap_Id.ref", "bin.size", "generation_time")) %dopar% {
                     tryCatch({
                       start_time <- Sys.time()
                       print(paste("[", format(start_time, "%Y-%m-%d %H:%M:%S"), "]", "Starting processing for country:", country.interest))
                       
                       # Filter data for the current country
                       Spike_metadata.filtered.interest.country <- Spike_metadata.filtered.interest[country == country.interest]
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Filtered data size for", country.interest, ":", nrow(Spike_metadata.filtered.interest.country)))
                       
                       if(nrow(Spike_metadata.filtered.interest.country) == 0){
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "No data available after filtering for country:", country.interest))
                         return(NULL)
                       }
                       
                       # Group by date.bin.num and hap_Id, count occurrences
                       Spike_metadata.filtered.interest.bin <- Spike_metadata.filtered.interest.country[, .(count = .N), by = .(date.bin.num, hap_Id)]
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Grouped data by date.bin.num and hap_Id"))
                       
                       # Calculate total counts for each hap_Id, filter those with count >= 20, order by descending total
                       hap.count.df <- Spike_metadata.filtered.interest.bin[, .(total = sum(count)), by = hap_Id][total >= 20][order(-total)]
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Number of haplotypes with count >= 20:", nrow(hap.count.df)))
                       
                       if(nrow(hap.count.df) == 0) {
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "No haplotypes with count >= 20 for", country.interest))
                         return(NULL)
                       }
                       
                       # Spread the data, filling NA with 0
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Spreading data into wide format"))
                       Spike_metadata.filtered.interest.bin.spread <- dcast(Spike_metadata.filtered.interest.bin[hap_Id %in% hap.count.df$hap_Id],
                                                                            date.bin.num ~ hap_Id,
                                                                            value.var = "count",
                                                                            fill = 0)
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Data spread completed with dimensions:", paste(dim(Spike_metadata.filtered.interest.bin.spread), collapse = "x")))
                       
                       # Create X matrix for the model
                       X <- cbind(1, Spike_metadata.filtered.interest.bin.spread$date.bin.num)
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "X matrix created with dimensions:", paste(dim(X), collapse = "x")))
                       
                       # Create Y matrix for the model
                       Y <- as.matrix(Spike_metadata.filtered.interest.bin.spread[, -1, with = FALSE])
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Y matrix created with dimensions:", paste(dim(Y), collapse = "x")))
                       
                       # Check for NA or infinite values in Y
                       if(any(is.na(Y)) || any(is.infinite(Y))){
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "Y matrix contains NA or infinite values. Skipping country:", country.interest))
                         return(NULL)
                       }
                       
                       # Order hap_Id
                       hap_Id.order.v <- c(Spike_hap_Id.ref, setdiff(hap.count.df$hap_Id, Spike_hap_Id.ref))
                       Y_cols <- colnames(Y)
                       hap_Id.order.v <- intersect(hap_Id.order.v, Y_cols)
                       
                       if(length(hap_Id.order.v) > 0) {
                         Y <- Y[, hap_Id.order.v, drop = FALSE]
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "Y matrix reordered. Current haplotypes:", paste(hap_Id.order.v, collapse = ", ")))
                       } else {
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "No matching haplotypes found for this country"))
                         return(NULL)
                       }
                       
                       # Calculate max date and date vector
                       max.date <- max(Spike_metadata.filtered.interest.country$date.bin.num)
                       date.v <- X[,2] / max.date
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Max date:", max.date))
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Date vector (first 5 values):", paste(head(date.v, 5), collapse = ", ")))
                       
                       # Check for variance in date.v
                       if(length(unique(date.v)) == 1) {
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "No variation in date.v for", country.interest, "- skipping"))
                         return(NULL)
                       }
                       
                       # Estimate the number of parameters
                       num_params <- (ncol(Y) - 1) * ncol(X)
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Number of parameters to estimate:", num_params))
                       
                       # Adjust MaxNWts if necessary
                       max_nwts <- max(1000, num_params * 10)
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Setting MaxNWts to:", max_nwts))
                       
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Fitting multinomial model..."))
                       # Fit multinomial model
                       fit.multinom <- multinom(Y ~ date.v, maxit = 10000, MaxNWts = max_nwts, trace = FALSE)
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Multinomial model fitted"))
                       
                       # Check for convergence
                       if(fit.multinom$convergence != 0){
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "Model did not converge for country:", country.interest))
                       } else {
                         print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                     "Model converged successfully for country:", country.interest))
                       }
                       
                       # Extract coefficients
                       fit.coef <- coef(fit.multinom)
                       colnames(fit.coef) <- c("intercept", "slope")
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Coefficients extracted"))
                       
                       # Convert coefficients to data.table
                       fit.coef <- as.data.table(fit.coef, keep.rownames = "hap_Id")
                       
                       # Display coefficients for debugging
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Model coefficients:\n"))
                       print(fit.coef)
                       
                       # Calculate relative_Re
                       fit.coef[, relative_Re := exp(((slope / max.date) / bin.size) * generation_time)]
                       fit.coef <- fit.coef[, .(hap_Id, relative_Re)]
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Relative_Re calculated"))
                       
                       # Add reference hap_Id
                       fit.coef <- rbindlist(list(data.table(hap_Id = Spike_hap_Id.ref, relative_Re = 1), fit.coef))
                       
                       # Add country column
                       fit.coef[, country := country.interest]
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Results for", country.interest, "completed"))
                       
                       # Clear large objects to free memory
                       rm(Spike_metadata.filtered.interest.country,
                          Spike_metadata.filtered.interest.bin,
                          Spike_metadata.filtered.interest.bin.spread,
                          X, Y)
                       gc()
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Memory cleared for country:", country.interest))
                       
                       end_time <- Sys.time()
                       elapsed_time <- difftime(end_time, start_time, units = "secs")
                       print(paste("[", format(end_time, "%Y-%m-%d %H:%M:%S"), "]",
                                   "Finished processing for country:", country.interest,
                                   "| Time elapsed:", round(elapsed_time, 2), "seconds"))
                       
                       fit.coef
                     }, error = function(e) {
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Error occurred while processing country:", country.interest))
                       print(paste("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]",
                                   "Error message:", e$message))
                       return(NULL)
                     })
                   }

print("All countries processed")

print("Finish parallel processing")

# CREATING INCIDENCE FILES----
print("Creating incidence files...")

# Filter for countries of interest
Spike_metadata.filtered.interest.country <- Spike_metadata.filtered.interest[country %in% country.interest.df$country]
print("Filtered data for countries of interest")

# Create bin counts
Spike_metadata.filtered.interest.bin <- Spike_metadata.filtered.interest.country[, .(count = .N), by = .(country, date.bin.num, hap_Id)]
print("Created binned counts")

# Create hap_country column
Spike_metadata.filtered.interest.bin[, hap_country := paste(hap_Id, country, sep="_")]
print("Created hap_country column")

Spike_metadata.filtered.interest.bin <- Spike_metadata.filtered.interest.bin[, .SD[sum(count) >= 20], by = hap_country]
print("Filtered out hap_country combinations with fewer than 20 counts")

# Get min and max dates
min_date <- min(Spike_metadata.filtered.interest.bin$date.bin.num)
max_date <- max(Spike_metadata.filtered.interest.bin$date.bin.num)

# Create a complete sequence of dates
all_dates <- data.table(Collection.date = min_date:max_date)
print(paste("Created sequence of all dates from", min_date, "to", max_date))

# Get unique hap_country combinations
unique_hap_country <- unique(Spike_metadata.filtered.interest.bin$hap_country)
print(paste("Number of unique hap_country combinations:", length(unique_hap_country)))

# Reshape the data from long to wide format
hap_counts_wide <- dcast(Spike_metadata.filtered.interest.bin,
                         date.bin.num ~ hap_country,
                         value.var = "count",
                         fill = 0)

# Merge with all_dates
hap_counts_daily <- merge(all_dates, hap_counts_wide,
                          by.x = "Collection.date", by.y = "date.bin.num",
                          all.x = TRUE)

# Replace NA with 0
setnafill(hap_counts_daily, fill = 0)

print("Filled in counts for each hap_country combination")

# Replace NA with 0
for(col in names(hap_counts_daily)[-1]) {
  set(hap_counts_daily, i = which(is.na(hap_counts_daily[[col]])), j = col, value = 0)
}
print("Replaced NA values with 0")

# Create weekly and 3-daily aggregations
hap_counts_daily[, Week_N := as.integer((Collection.date - min_date) / 7) + 1]
hap_counts_daily[, ThreeDay_N := as.integer((Collection.date - min_date) / 3) + 1]

hap_counts_weekly <- hap_counts_daily[, lapply(.SD, sum), by = Week_N, .SDcols = -c("Collection.date", "ThreeDay_N")]
hap_counts_3daily <- hap_counts_daily[, lapply(.SD, sum), by = ThreeDay_N, .SDcols = -c("Collection.date", "Week_N")]

setorder(hap_counts_weekly, Week_N)
setorder(hap_counts_3daily, ThreeDay_N)
print("Created weekly and 3-daily aggregations")

# Print summary information
print(paste("Dimensions of hap_counts_daily:", paste(dim(hap_counts_daily), collapse = "x")))
print(paste("Dimensions of hap_counts_weekly:", paste(dim(hap_counts_weekly), collapse = "x")))
print(paste("Dimensions of hap_counts_3daily:", paste(dim(hap_counts_3daily), collapse = "x")))

# Save count files to output directory with prefixed filenames
fwrite(hap_counts_daily, file.path(output_dir, paste0(output_prefix, "_hap_counts_daily.csv")))
fwrite(hap_counts_weekly, file.path(output_dir, paste0(output_prefix, "_hap_counts_weekly.csv")))
fwrite(hap_counts_3daily, file.path(output_dir, paste0(output_prefix, "_hap_counts_3daily.csv")))
print("Saved count files to output directory with prefixed filenames")

print("Incidence files creation completed")

# Post-processing ----

# Combine all results into a single dataframe
print("Combining all results into a single dataframe.")
Spike_data.coef.df <- as.data.frame(results)
print("Combined results.")

# Stop the parallel cluster after the job is done
print("Stopping the parallel cluster.")
stopImplicitCluster()
print("Parallel cluster stopped.")

# Plotting results
print("Plotting results.")
g <- ggplot(results, aes(x=country, y=relative_Re)) +
  geom_violin()
print("Plot created. Saving to PDF.")
pdf(file.path(output_dir, paste0(output_prefix, "_violin_plot.pdf")))
print(g)
dev.off()
print(paste("Plot saved to", file.path(output_dir, paste0(output_prefix, "_violin_plot.pdf"))))

print("Loading mutation information.")
Spike_mut.info <- fread(Spike_mut.info.name, header=TRUE, sep="\t", quote="", check.names=TRUE)

Spike_mut.info <- Spike_mut.info %>% filter(grepl("S:",mut))
print("Filtered for S gene mutations.")

Spike_mut.info <- Spike_mut.info %>% right_join(Spike_metadata.filtered.interest %>% dplyr::select(seqName,hap_Id), by = "seqName")
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

Spike_metadata.filtered.interest <- Spike_metadata.filtered.interest %>% inner_join(Spike_mut.info.sum.all_S %>% dplyr::select(seqName,hap_Id.all_mut),by="seqName")
print("Joined summarized haplotype info with metadata.")

date.quantile.df <- Spike_metadata.filtered.interest %>% mutate(Collection.date = as.Date(Collection.date)) %>% group_by(hap_Id) %>% summarize(date.first = quantile(Collection.date,0.01,type=1))
date.quantile.df <- date.quantile.df %>% arrange(date.first)
print("First collection date quantile calculated.")

# Determining major haplotypes
all_S.major.hap.df <- Spike_metadata.filtered.interest %>%
  dplyr::group_by(hap_Id, hap_Id.all_mut) %>%
  dplyr::summarize(count = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(hap_Id) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::top_n(1, count) %>%
  dplyr::filter(count == dplyr::first(count))

print("Major haplotype determined.")
print(all_S.major.hap.df)

Spike_metadata.filtered.interest.representative <- Spike_metadata.filtered.interest %>% inner_join(all_S.major.hap.df %>% dplyr::select(hap_Id,hap_Id.all_mut),by=c("hap_Id","hap_Id.all_mut")) %>% group_by(hap_Id) %>% sample_n(1)
Spike_metadata.filtered.interest.representative <- Spike_metadata.filtered.interest.representative %>% dplyr::select(seqName,clade,Nextclade_pango,hap_Id,Accession.ID)
Spike_metadata.filtered.interest.representative <- Spike_metadata.filtered.interest.representative %>% inner_join(Spike_hap.info.df,by="hap_Id")
print("Representative sequences selected.")

print("Colnames of Spike_metadata.filtered.interest")
print(colnames(Spike_metadata.filtered.interest))

print("Colnames of Spike_metadata.filtered.interest.representative")
print(colnames(Spike_metadata.filtered.interest.representative))



Spike_data.coef.df.merged <- Spike_data.coef.df %>% inner_join(Spike_metadata.filtered.interest.representative,by="hap_Id")
Spike_data.coef.df.merged <- Spike_data.coef.df.merged %>% arrange(desc(relative_Re))
print("Coefficient data merged with representative metadata.")

Spike_data.coef.df.merged <- Spike_data.coef.df.merged %>% left_join(date.quantile.df,by="hap_Id")
print("Joined with date quantile data.")

print("Writing the output to a file.")
out.name <- file.path(output_dir, paste0(output_prefix, "_metadata_representative_all_countries_with_date.txt"))
write.table(Spike_data.coef.df.merged, out.name, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
print(paste("Output written to", out.name))

print("Saving Spike_metadata")
# Filter Spike_metadata_to_save based on the combination of country and hap_Id
Spike_metadata_to_save_filtered <- Spike_metadata_to_save %>%
  dplyr::semi_join(Spike_data.coef.df.merged, by = c("country", "hap_Id"))

# Save the extracted data to a CSV file
fwrite(
  Spike_metadata_to_save_filtered,
  file.path(output_dir, paste0(output_prefix, "_metadata_filtered_interest.csv"))
)

print("Spike_metadata saved")





# Plotting results from original Jumpei Ito et al. paper ----


Spike_count.df.hap.all_regions <- Spike_metadata.filtered.interest %>% group_by(hap_Id) %>% summarize(count.total = n())
print("Count dataframe created.")

Spike_data.coef.df.merged <- Spike_data.coef.df.merged %>% inner_join(Spike_count.df.hap,by=c("hap_Id","country"))
print("Coefficient data merged.")

print("Creating and saving time_vs_fitness_each_country_log.pdf.")
g <- ggplot(Spike_data.coef.df.merged, aes(x=date.first,y=log(relative_Re), size = count.hap, fill = clade))
g <- g + geom_point(shape = 21, alpha = 0.7)
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + geom_hline(yintercept=c(-1,-0.5,0,0.25,0.5,1),color="gray70")
g <- g + facet_wrap(~country)
g <- g + xlab("Date") + ylab("Relative Re (log)")
#g <- g + scale_y_continuous(lim=c(0,2))

pdf.name <- file.path(output_dir, paste0(output_prefix, "Spike_time_vs_fitness_each_country_log.pdf"))
pdf(pdf.name, width = 20, height = 10)
print(g)
dev.off()
print(paste("First plot saved to", pdf.name))

g <- ggplot(Spike_data.coef.df.merged, aes(x=date.first,y=log(relative_Re+0.1), size = count.hap, fill = clade))
g <- g + geom_point(shape = 21, alpha = 0.7)
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + geom_hline(yintercept=c(-1,-0.5,0,0.25,0.5,1),color="gray70")
g <- g + facet_wrap(~country)
g <- g + xlab("Date") + ylab("Relative Re (log)")
#g <- g + scale_y_continuous(lim=c(0,2))

pdf.name <- file.path(output_dir, paste0(output_prefix, "Spike_time_vs_fitness_each_country_log_with_pseudo.pdf"))
pdf(pdf.name, width = 20, height = 10)
print(g)
dev.off()
print(paste("Second plot saved to", pdf.name))

g <- ggplot(Spike_data.coef.df.merged, aes(x=date.first,y=relative_Re, size = count.hap, fill = clade))
g <- g + geom_point(shape = 21, alpha = 0.7)
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + geom_hline(yintercept=c(1,2,3,4),color="gray70")
g <- g + facet_wrap(~country)
g <- g + xlab("Date") + ylab("Relative Re")
#g <- g + scale_y_continuous(lim=c(0,2.5))

pdf.name <- file.path(output_dir, paste0(output_prefix, "Spike_time_vs_fitness_each_country.pdf"))
pdf(pdf.name, width = 20, height = 10)
print(g)
dev.off()
print(paste("Third plot saved to", pdf.name))

# Separate haplotype into separate rows
haplotype.mut_profile <- Spike_metadata.filtered.interest.representative %>%
  dplyr::select(hap_Id, haplotype) %>%
  tidyr::separate_rows(haplotype, sep = ",\\s*")
print("Step 1: Haplotype rows separated.")
print(head(haplotype.mut_profile))

# Extract numerical position data from the haplotype string
haplotype.mut_profile <- haplotype.mut_profile %>%
  mutate(pos = gsub("[^0-9]", "", haplotype), pos = as.numeric(pos))
print("Step 2: Numeric positions extracted.")
print(head(haplotype.mut_profile))

# Join additional data from the same dataset
haplotype.mut_profile <- haplotype.mut_profile %>%
  dplyr::inner_join(Spike_metadata.filtered.interest.representative %>%
                      dplyr::select(hap_Id, clade), by = "hap_Id")
print("Step 3: Data joined with clade information.")
print(head(haplotype.mut_profile))

# Final data manipulation
haplotype.mut_profile <- haplotype.mut_profile %>%
  dplyr::ungroup() %>%
  mutate(hap_Id = factor(hap_Id, levels = date.quantile.df$hap_Id)) %>%
  arrange(hap_Id)
print("Step 4: Final manipulation and arrangement.")
print(head(haplotype.mut_profile))

print("Haplotype mutation profile analyzed.")

# Existing plotting code
count.hap_with_mut_each_pos <- haplotype.mut_profile %>% group_by(pos) %>% summarize(count.hap = n()) %>% arrange(pos)
print("Count of haplotypes with mutations at each position calculated.")

g2 <- ggplot(count.hap_with_mut_each_pos,aes(x=pos,y=log(count.hap,10)))
g2 <- g2 + geom_point(size=0.01)
g2 <- g2 + xlab("") + ylab("log10(count)")

g1 <- ggplot(haplotype.mut_profile,aes(x=pos,y=hap_Id, color = clade))
g1 <- g1 + geom_point(size=0.01)

g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g1 <- g1 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g1 <- g1 + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 strip.text = element_text(size=8)
)
g1 <- g1 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
g1 <- g1 + xlab("Amino acid position") + ylab(paste("haplotype (",as.character(nrow(date.quantile.df)),")",sep=""))

pdf.name <- file.path(output_dir, "Spike_mut_heatmap.pdf")
pdf(pdf.name, width = 10, height = 8)
print(g2 / g1 + plot_layout(ncol = 1, heights = c(1, 3)))
dev.off()
print(paste("Mutation heatmap saved to", pdf.name))




# Boxplot
ggplot(Spike_data.coef.df.merged, aes(x = clade, y = relative_Re, fill = clade)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Relative Re by Clade",
       x = "Clade",
       y = "Relative Re") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Save the plot
ggsave(file.path(output_dir, paste0(output_prefix, "RelativeRe_by_Clade_Boxplot.pdf")), width = 12, height = 6)






# Calculate the number of mutations per haplotype
mutation_counts <- haplotype.mut_profile %>%
  group_by(hap_Id) %>%
  summarize(mutation_count = n())

# Merge with fitness data
fitness_vs_mutations <- Spike_data.coef.df.merged %>%
  inner_join(mutation_counts, by = "hap_Id")

# Plot
ggplot(fitness_vs_mutations, aes(x = mutation_count, y = relative_Re)) +
  geom_point(alpha = 0.6, aes(color = clade)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal() +
  labs(title = "Relative Re vs. Mutation Count",
       x = "Number of Mutations in Haplotype",
       y = "Relative Re") +
  theme(legend.position = "right")

# Save the plot
ggsave(file.path(output_dir, paste0(output_prefix, "RelativeRe_vs_MutationCount.pdf")), width = 12, height = 6)



print("Script execution completed.")

