# Load necessary libraries
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(ggplot2)
# library(ggpubr)
# library(grid)
# library(data.table) # fread

# devtools::document()

# From the package root directory
# devtools::build()
# devtools::install()
devtools::load_all() # Automatically points to the current directory as the package root

# varsome ----
# LE = less than equal to, GE = greater than equal to
# varsome <- read.csv(file = "../../data/singlecase/varsome_calibrated_insilico_thresholds.tsv", sep="\t")

# user variables ----
af_threshold <- 0.1  # Allele frequency threshold
gnomad_max <- 1e-6
# metadataclass_file_path <- "../output/custom_reference_metadata.csv"
samples_file_path <- "../data/samples.tsv" # phenotype data
output_path <- "../output/"

# input single file
input_path <- "../data/phrtmma_v1_chr21_40411318_41411317.csv"

# input list of files
input_path <- c(
"../data/phrtmma_v1_chr21_40411318_41411317.csv",
"../data/phrtmma_v1_chr21_41411318_42411317.csv",
"../data/phrtmma_v1_chr21_42411318_43411317.csv"
)

# input all files
input_path <- "../data/"

# start analysis ----

processed_data_list <- process_genetic_data(input_path, samples_file_path, af_threshold)

# Check that all imported data chucks detected the correct column classes
# compare_column_classes_and_output_csv(processed_data_list, metadataclass_file_path)
compare_column_classes_and_output_csv(processed_data_list)

# Apply the corrected column classes to each dataset in processed_data_list
# processed_data_list <- apply_column_classes_to_processed_data(processed_data_list, metadataclass_file_path)

processed_data_list <- apply_column_classes_to_processed_data(processed_data_list)

# Merge all dataframes in the list into a single dataframe
all_data <- bind_rows(processed_data_list)

df <- all_data |> head()

# processed_data_list[[1]] |> nrow()
# processed_data_list[[2]] |> nrow()
# processed_data_list[[1]] |> ncol()
# processed_data_list[[2]] |> ncol()

# df1 <- processed_data_list[[1]]
# df2 <- processed_data_list[[2]]

# Ensure all_data is indeed a dataframe
if (!is.data.frame(all_data)) {
  stop("all_data is not a dataframe")
}

# Now, you have a single dataframe 'all_data' with data from all processed files

# Apply ACMG labels preparation on the aggregated data
all_data <- prepare_acmg_labels(all_data)

# Define a generic suffix or use a specific identifier for your output files
file_suffix <- "aggregate"

# Generate plots based on the aggregated data
plot_criteria_count_each_gene(all_data, file_suffix)
plot_criteria_gene_total(all_data, file_suffix)
plot_variants_per_criteria(all_data, file_suffix)

# Optionally, you can save the aggregated dataframe for further analysis
# write.csv(all_data, paste0("../output/aggregated_data_", Sys.Date(), ".csv"), row.names = FALSE)


df_test <- apply_acmg_pp3_with_external_scores(df, varsome_data)

# select ----
# use dt to make selection using manual list or indexes
# manually select columns to report

columns_to_select <- c(
  "sample",
  "comp_het_flag",
  "SYMBOL",
  "ACMG_PVS1",
  "ACMG_PS1",
  "ACMG_PS5",
  "ACMG_PM2",
  "ACMG_PM3",
  "ACMG_highest",
  "ACMG_count",
  "AF.x",
  "CLIN_SIG",
  "IMPACT",
  "genotype",
  # "Inheritance",
  "gnomAD_AF",
  # "Engine",
  "Strong_pathogenic_GE",
  "Moderate_pathogenic_GE",
  "Supporting_pathogenic_GE",
  "BayesDel_addAF_score",
  "BayesDel_noAF_score",
  # "CADD_PHRED",
  "DANN_score",
  "Eigen.raw_coding",
  "Eigen.PC.phred_coding",
  "FATHMM_score",
  "fathmm.MKL_coding_score",
  "fathmm.XF_coding_score",
  "LRT_score",
  "M.CAP_score",
  "MetaLR_score",
  "MetaSVM_score",
  "MetaRNN_score",
  "MutPred_score",
  "MutationAssessor_score",
  "MutationTaster_score",
  "phastCons100way_vertebrate",
  "Polyphen2_HDIV_score",
  "Polyphen2_HVAR_score",
  "PROVEAN_score",
  "REVEL_score",
  "SIFT_score"
)

# df_selected <- df %>% dplyr::select(dplyr::all_of(columns_to_select))
library(data.table)
dt <- as.data.table(all_data)
first_50_indices <- 1:50
first_50_names <- names(dt)[first_50_indices]
all_columns <- unique(c(first_50_names, columns_to_select)) # all column names
dt_selected <- dt[, ..all_columns] # Select columns in data.table

# Report ----

library(dplyr)
dt_selected <- dt_selected |> dplyr::filter(ACMG_count >= 2)
dt_test <- dt_selected |> head(1)

dt_test <- dt_test |> select(sample, HGVSc, HGVSp, ACMG_count)
# write.csv(dt_test, "../output/dt_test.csv", row.names = FALSE)
write.csv(dt_test, "../output/dt_test.csv", row.names = FALSE, quote = FALSE)

names(dt_test)

# sample ID
dt_test$sample
# has coding variant:
dt_test$HGVSc
# which encodes the protein variant:
dt_test$HGVSp
# and is given a ACMG score:
dt_test$ACMG_count

# Assuming dt_test is already prepared and contains at least one row
dt_test <- dt_selected |> head(1) |> dplyr::select(sample, HGVSc, HGVSp, ACMG_count)


# Function to escape LaTeX special characters in text
escape_latex <- function(text) {
  text <- gsub("_", "\\\\_", text)
return(text)
}

# Assuming dt_test is prepared
dt_test <- dt_selected |> head(1) |> dplyr::select(sample, HGVSc, HGVSp, ACMG_count)

# Apply the escape function
dt_test <- data.frame(lapply(dt_test, escape_latex), stringsAsFactors = FALSE)

# Prepare LaTeX commands
latex_commands <- sprintf("
\\newcommand{\\SampleID}{%s}
\\newcommand{\\CodingVariant}{%s}
\\newcommand{\\ProteinVariant}{%s}
\\newcommand{\\ACMGScore}{%s}
",
dt_test$sample,
dt_test$HGVSc,
dt_test$HGVSp,
dt_test$ACMG_count)

# Write LaTeX commands to a file
writeLines(latex_commands, "../output/variables.tex")
