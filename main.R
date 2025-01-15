# Load required libraries
library(dplyr)       # For data manipulation
library(tidyr)       # For data tidying
library(R.matlab)    # To read MATLAB (.mat) files

# Setting the current working directory
wd <- getwd()       
print(wd)           

# Load custom R scripts for Flux Balance Analysis (FBA) functions and utilities
source(paste0(wd, "/epimod_FBAfunctions/R/FBAgreatmodeClass.R"))    # Core FBA class functions
source(paste0(wd, "/epimod_FBAfunctions/R/class_generation.R"))     # Functions for generating FBA models
source(paste0(wd, "/epimod_FBAfunctions/R/readMat.R"))              # MATLAB file reading utilities

# Define the models to be processed
models <- list(
  list(name = "Escherichia_coli_str_K_12_substr_MG1655", type = "ecoli"),
  list(name = "Anaerostipes_caccae_DSM_14662", type = "Anaerostipes_caccae"),
  list(name = "Bacteroides_thetaiotaomicron_VPI_5482", type = "Bacteroides_thetaiotaomicron"),
  list(name = "Bifidobacterium_longum_NCC2705", type = "Bifidobacterium_longum"),
  list(name = "Blautia_producta_DSM_2950", type = "Blautia_producta"),
  list(name = "Clostridium_butyricum_DSM_10702", type = "Clostridium_butyricum"),
  list(name = "Clostridium_ramosum_VPI_0427_DSM_1402", type = "Clostridium_ramosum"),
  list(name = "Lactobacillus_plantarum_subsp_plantarum_ATCC_14917", type = "Lactobacillus_plantarum")
)

# Set the working directory to the current directory
setwd(wd)

# Initialize total time for tracking execution duration
total_time <- 0  

# Create a results directory if it does not already exist
dir.create("./results", showWarnings = FALSE)

# Define the output directory for results
dest_dir <- paste0(wd, "/results/")
print(dest_dir)  # Print the destination directory for verification

# Loop through each model for FBA processing
for (model_info in models) {
  model_name <- model_info$name    # Extract the model name
  model_type <- model_info$type    # Extract the model type
  
  # Construct file paths for FBA input
  fba_fname <- paste0(model_name, ".txt")
  matfile <- paste0(wd, '/input/', model_type, "/", model_name, ".mat")
  print(matfile)  # Print the .mat file path for debugging
  
  # Check if the .mat file exists
  if (!file.exists(matfile)) {
    stop(paste("File not found:", matfile))  # Stop execution if file is missing
  }
  
  # Generate the model for the specified .mat file
  cat("Generating model for:", model_name, "\n")
  model <- FBA4Greatmod.generation(fba_mat = matfile)
  
  # Set biomass parameters for the FBA model
  model <- setBiomassParameters(model, bioMax = 1.172, bioMean = 0.489, bioMin = 0.083)
  
  # Optionally, set the pFBA flag:
  #   geneOption = 0: Minimizes all reactions
  #   geneOption = 1: Minimizes only gene-associated reactions
  #   geneOption = 2: Minimizes only non-gene-associated reactions
  #   Default = -1: pFBA is not applied
  # model <- setPFbaGeneOption(model, geneOption = 1)
  
  # Prepare the output path for .rds files
  output_rds_path <- paste0(wd, "/input/models/", model_type)
  files_to_move <- list.files(pattern = "\\.rds$", full.names = TRUE)
  
  # Move .rds files to the appropriate directory
  if (length(files_to_move) > 0) {
    sapply(files_to_move, function(f) {
      file.rename(from = f, to = file.path(output_rds_path, basename(f)))
    })
  } else {
    cat("No .rds files found for", model_name, "\n")
  }
  
  # Measure the time taken to write the FBA file
  start_time <- Sys.time()
  writeFBAfile(model, model_name, dest_dir)
  end_time <- Sys.time()
  
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  total_time <- total_time + elapsed_time
  cat(paste("Time to write", model_name, ".txt:", elapsed_time, "seconds\n"))
}

# Print the total time elapsed for processing all models
cat("Total time elapsed:", total_time, "seconds\n")

### Setting Upper Bounds for FBA reactions and Not FBA reactions ###

# Load the script for setting bounds for exchange reactions
source(paste0(wd, "/epimod_FBAfunctions/R/ex_bounds_module.R"))

# Define the FBA models and associated counts
bacteria_models <- c(
  "./results/Clostridium_butyricum_DSM_10702.txt",
  "./results/Escherichia_coli_str_K_12_substr_MG1655.txt",
  "./results/Bacteroides_thetaiotaomicron_VPI_5482.txt",
  "./results/Clostridium_ramosum_VPI_0427_DSM_1402.txt",
  "./results/Lactobacillus_plantarum_subsp_plantarum_ATCC_14917.txt",
  "./results/Blautia_producta_DSM_2950.txt",
  "./results/Bifidobacterium_longum_NCC2705.txt",
  "./results/Anaerostipes_caccae_DSM_14662.txt"
)

bacteria_counts <- c(23796, 377148, 2325473, 4231665, 2435, 3439626, 2860957, 969418)

# This function extracts all the specified EX_ reactions from the provided models (bacteria_files)
# and calculates their upper bound as `non_fba_base_bound / bacteria_counts` for the corresponding position.
# Additionally, `fba_upper_bound` is applied to the specific reactions listed in `fba_reactions`.
run_full_ex_bounds(
  extraction_output  = "extracted_ex_reactions.txt",  # Output file for extracted reactions
  bacteria_files     = bacteria_models,              # List of FBA models to process
  output_dir         = "results_ex_reactions",       # Directory to store results
  fba_reactions      = c("EX_but_e", "EX_ac_e", "EX_ppa_e", "EX_lac_L_e"), # Specific FBA reactions
  bacteria_counts    = bacteria_counts,              # Count of bacteria used for bound calculation
  non_fba_base_bound = 1000,                         # Upper bound for non-FBA reactions
  fba_upper_bound    = 0.015                         # Upper bound for specified FBA reactions
)

