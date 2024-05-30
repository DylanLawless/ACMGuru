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
input_path <- "../data/study_v1_chr21_40411318_41411317.csv"

# input list of files
input_path <- c(
"../data/study_v1_chr21_40411318_41411317.csv",
"../data/study_v1_chr21_41411318_42411317.csv",
"../data/study_v1_chr21_42411318_43411317.csv"
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


# df_test <- apply_acmg_pp3_with_external_scores(df, varsome_data)

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


# acmg tally  ----
# SEE:  ./archive/AMCGuru_singlecase_vcurrent.R
# SEE:  ./archive/AMCGuru_singlecase_vcurrent.R
# # SEE:  ./archive/AMCGuru_singlecase_vcurrent.R
# # SEE:  ./archive/AMCGuru_singlecase_vcurrent.R
# List of all ACMG labels
acmg_labels <- c("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5",
                 "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6",
                 "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4")

# Check if each ACMG column exists, if not create it and fill with NA
for (acmg_label in acmg_labels) {
  if (!acmg_label %in% names(df)) {
    df[[acmg_label]] <- NA
  }
}

# Then use coalesce to find the first non-NA ACMG label
df$ACMG_highest <- dplyr::coalesce(!!!df[acmg_labels])
df <- df %>% dplyr::select(ACMG_highest, everything())

# Count the number of non-NA values across the columns
df$ACMG_count <- rowSums(!is.na(df[, acmg_labels ]))
df <- df %>% dplyr::select(ACMG_count, everything())
# df$ACMG_count[df$ACMG_count == 0] <- NA

p.criteria_count_each_gene <- df |>
  filter(ACMG_count > 1) |>
  ggplot(aes(y = ACMG_count, x = SYMBOL)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x  = element_text(angle=45, hjust=1, vjust=1)) +
  xlab("Gene symbol") +
  ylab("ACMG criteria count (>1)")
p.criteria_count_each_gene
# ggsave(paste("../../data/singlecase/", file_suffix, "criteria_count_each_gene.pdf", sep = "") ,plot = p.criteria_count_each_gene )

# as table
df |>
  filter(ACMG_count > 1) |>
  dplyr::select(sample, SYMBOL, ACMG_count) |>
  arrange(desc(ACMG_count))

p.criteria_gene_total <- df %>%
  group_by(SYMBOL) |>
  summarise(acmg_count_per_symbol = sum(ACMG_count)) |>
  na.omit() |>
  ggplot(aes(x = acmg_count_per_symbol, fill=..x..) ) +
  geom_histogram(stat="count", binwidth = 1, color="black"
  ) +
  theme_minimal() +
  xlab("No. ACMG criteria (P) variants per gene") +
  ylab("Number of genes") +
  geom_text(stat='count', aes(label=..count.., y=..count..+20), color = "black") +
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.criteria_gene_total
# ggsave(paste("../../data/singlecase/", file_suffix, "criteria_gene_total.pdf", sep = "") ,plot = p.criteria_gene_total )

# as table
df |>
  group_by(SYMBOL) |>
  summarise(acmg_count_per_symbol = sum(ACMG_count)) |>
  na.omit() |>
  arrange(desc(acmg_count_per_symbol))

p.variants_per_criteria <- df |>
  ggplot(aes(x = ACMG_count, fill=..x..)) +
  geom_histogram(binwidth = 1, color="black") +
  xlab("No. ACMG criteria\nassigned (P)") +
  ylab("No. variants") +
  theme_minimal() +
  geom_text(stat='count', aes(label=..count.., y=..count..+300), color = "black") +
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.variants_per_criteria
# ggsave(paste("../../data/singlecase/", file_suffix, "variants_per_criteria.pdf", sep = "") ,plot = p.variants_per_criteria , width = 9, height = 5)

# Check we only have approx. 1 "casual" variant per sample
p.criteria_per_sample <- df %>%
  group_by(sample) %>%
  summarise(ACMG_count = max(ACMG_count, na.rm = TRUE))  %>%
  ggplot(aes(x = ACMG_count, fill=..x..)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(x = "No. ACMG criteria\nassigned (P)", y = "No. samples") +
  theme_minimal() +
  geom_text(stat='count', aes(label=..count.., y=..count..+20), color = "black") +
  guides(fill=FALSE) +
  scale_fill_scico(palette = 'acton', direction = 1) # batlowK, acton, lajolla, lapaz, turku
p.criteria_per_sample
# ggsave(paste("../../data/singlecase/", file_suffix, "criteria_per_sample.pdf", sep = "") ,plot = p.criteria_per_sample )

# as table
df |>
  group_by(sample, ACMG_count) |>
  tally(n = "count_per_sample") |>
  ungroup() |>
  dplyr::select(-sample) |>
  group_by(ACMG_count) |>
  tally(n = "count_per_sample")

# ACMG Verdict----
# Rules are combined using the point system described in PMID:32720330
# Each rule triggered is assigned a number of points based on the strength of the evidence provided:
#
# Supporting: 1 point
# Moderate: 2 points
# Strong: 4 points
# Very Strong: 8 points
# A total score is computed as the sum of the points from the pathogenic rules, minus the sum of the points from benign rules.
#
# The total score is then compared to thresholds to assign the final verdict:
#
# Pathogenic if greater than or equal to 10,
# Likely Pathogenic if between 6 and 9 inclusive,
# Uncertain Significance if between 0 and 5,
# Likely Benign if between -6 and -1,
# Benign if less than or equal to -7.

# Benign if less than or equal to -7.

df <-  df |> dplyr::select("ACMG_PVS1", "ACMG_PS1", "ACMG_PS2", "ACMG_PS3", "ACMG_PS4", "ACMG_PS5",
                           "ACMG_PM1", "ACMG_PM2", "ACMG_PM3", "ACMG_PM4", "ACMG_PM5", "ACMG_PM6",
                           "ACMG_PM7", "ACMG_PP1", "ACMG_PP2", "ACMG_PP3", "ACMG_PP4",
                           everything())


# Define scores for each ACMG label
acmg_scores <- c("PVS1" = 8,
					  "PS1" = 4, "PS2" = 4, "PS3" = 4, "PS4" = 4, "PS5" = 4,
					  "PM1" = 2, "PM2" = 2, "PM3" = 2, "PM4" = 2, "PM5" = 2,
					  "PP3" = 1)

# Create ACMG_score column by looking up ACMG_highest in acmg_scores
df$ACMG_score <- acmg_scores[df$ACMG_highest]

# If there are any ACMG labels that don't have a corresponding score, these will be NA. You may want to set these to 0.
df$ACMG_score[is.na(df$ACMG_score)] <- 0
df <- df |> dplyr::select(ACMG_score, everything())













# Total ACMG score ----
# # List of all ACMG columns
# acmg_columns <- grep("ACMG_P", colnames(df), value = TRUE)
# acmg_columns
#
# # Mutate all ACMG columns
# df <- df %>%
#   mutate_at(acmg_columns, function(x) acmg_scores[x])
#
# # Replace NAs with 0 in ACMG columns only
# df[acmg_columns] <- lapply(df[acmg_columns], function(x) ifelse(is.na(x), 0, x))
#
# # Calculate total ACMG score
# df$ACMG_total_score <- rowSums(df[acmg_columns])
#
# df <- df |> dplyr::select(ACMG_total_score, everything())
#
# p.acmg_score <- df |>
# 	ggplot(aes(x = as.character(ACMG_total_score), fill= as.numeric(ACMG_total_score) )) +
# 	geom_histogram(stat='count', bins = length(acmg_scores), color="black") +
# 	theme_minimal() +
# 	xlab("ACMG score") +
# 	ylab("No. variants") +
# 	geom_text(stat='count', aes(label=..count.., y=..count..+300), color = "black") +
# 	guides(fill=FALSE) +
# 	scale_fill_scico(palette = 'bamako', direction = 1) # batlowK, acton, lajolla, lapaz, turku
# p.acmg_score
# # ggsave(paste("../../data/singlecase/", file_suffix, "acmg_score.pdf", sep = "") ,plot = p.acmg_score )
#
