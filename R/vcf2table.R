# This script converts the complex annotation data structure from VCF to a
# simple dataframe. The main package dependency is: VariantAnnotation.

# -----------
# set up ----
# -----------

args <- commandArgs(trailingOnly = TRUE)
version <- args[1]
output_dir <- args[2]
file_suffix <- args[3]
sample_count <- as.numeric(args[4])
vcf_out_path <- args[5]
file_idx <- as.numeric(args[6])  # This captures the SLURM_ARRAY_TASK_ID

sink(paste("../../data/log/filter_maf04/annotated/annotated_launch_", "chr", file_idx, "_", version, "_", file_suffix,  ".log", sep = ""), append = TRUE)
print(paste("Arguements:", args))

# check sample count
if(sample_count == 0) {
        warning("The number of samples (sample_count) is 0.")
        stop("The number of samples in the dataset must be greater than 0.")
} else {
        cat("Number of samples:", sample_count, "\n")
}

library(dplyr)
library(tidyr)
library(stringr)
library(VariantAnnotation)

# --------------------
# function set up ----
# --------------------

read_and_convert_vcf <- function(vcfFile) {
        cat("\n------------------")
        cat("\nread and convert start")
        cat("\n------------------")

        cat("\ngetwd()\n")
        print(getwd())
        command <- paste("du -sh", vcfFile)
        result <- system(command, intern = TRUE)
        print(result)

        cat("\ntabix\n")
        tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

        cat("\nvcf\n")
        vcf <- VariantAnnotation::readVcf(file = tabixVcf)

        cat("\nevcf\n")
        evcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)
        rm(vcfFile, tabixVcf, vcf)
        gc()

        cat("\ncsq\n")
        csq <- ensemblVEP::parseCSQToGRanges(x = evcf)
        df_csq <- as.data.frame(csq, row.names = NULL)
        rm(csq)
        gc()

        cat("\nrownames\n")
        df_csq$rownames <- ensemblVEP::parseCSQToGRanges(x = evcf) |> names()
        df_geno <- as.data.frame(geno(evcf)[["GT"]])
        df_info <- as.data.frame(info(evcf))
        rm(evcf)
        gc()

        df_geno$rownames <- df_geno |> rownames()
        df_info$rownames <- df_geno$rownames # note the reuse of rownames since package version updated

        cat("\nmerge geno info")
        df_merge <- merge(df_geno, df_info, by = "rownames")
        df_merge <- df_merge |> dplyr::select(-CSQ)
        rm(df_geno, df_info)
        gc()

        cat("\nmerge csq")
        df_main <- merge(df_merge, df_csq, by = "rownames")
        rm(df_merge, df_csq)
        gc()

        return(df_main)
}

# filter_and_gather <- function(df_main, sample_count, output_dir, file_suffix, version, f) {
filter_and_gather <- function(df_main, sample_count, output_dir, file_suffix, version, f, split_file) {
        cat("\n------------------")
        cat("\nfilter and gather start")
        cat("\n------------------")

        # Temp debug - output before gather ///////////////////////////////////////
        # # Extract the ${chromosome}_${start}_${end} from split_file
        # pattern <- paste("chr", chr_names[f], "_(\\d+)_(\\d+)", sep="")
        # matches <- regmatches(split_file, regexpr(pattern, split_file))
        # region_info <- matches[1]  # This extracts the matched pattern from the filename
        # # Clone df_main to df_temp
        # df_temp <- df_main
        # # Identify list columns
        # list_columns <- names(df_temp)[sapply(df_temp, class) == "list"]
        # # Replace list columns in df_temp with a placeholder
        # df_temp[list_columns] <- "PLACEHOLDER_VALUE"
        # # Construct the file path for df_temp using region_info
        # before_gather_file_path_temp <- paste(output_dir, file_suffix, "_", version, "_", region_info, "_before_gather_TEMP.csv", sep = "")
        # # Write df_temp to the constructed file path
        # utils::write.table(df_temp, before_gather_file_path_temp, row.names = FALSE, sep = ",")
        # print(paste("df_temp before gather for chr", chr_names[f], "saved to:", before_gather_file_path_temp))
        # Temp debug - output before gather ///////////////////////////////////////

        df_main <- df_main |> dplyr::filter(IMPACT %in% c("HIGH", "MODERATE"))

        df_main$AC <- as.numeric(df_main$AC)
        df_main$AF.x <- as.numeric(df_main$AF.x)
        df_main$MLEAC <- as.numeric(df_main$MLEAC)
        df_main$MLEAF <- as.numeric(df_main$MLEAF)

        n_sample_col_start <- 1+1
        n_sample_col_end <- 1+sample_count

        df <- tidyr::gather(df_main, all_of(n_sample_col_start:n_sample_col_end),
                                                key = "sample", value = "genotype_call") |>
                        dplyr::select(sample, genotype_call, SYMBOL, HGVSp, HGVSc, Consequence, IMPACT, everything())

        rm(df_main)
        gc()

        df$SNP <- df$"rownames"
        rm(n_sample_col_start, n_sample_col_end)

        cat("Class of df : ", class(df), "\n")
        return(df)
}

clean_genotype <- function(df) {
        cat("\n------------------")
        cat("\ngenotype_clean start")
        cat("\n------------------")

        # Create new column "genotype"
        df$genotype <- df$genotype_call
        df$genotype[df$genotype_call == "0/0"] <- "0"
        df$genotype[df$genotype_call == "./0"] <- "0"
        df$genotype[df$genotype_call == "0/."] <- "0"
        df$genotype[df$genotype_call == "0|0"] <- "0"
        df$genotype[df$genotype_call == "0/1"] <- "1"
        df$genotype[df$genotype_call == "1/0"] <- "1"
        df$genotype[df$genotype_call == "0|1"] <- "1"
        df$genotype[df$genotype_call == "./1"] <- "1"
        df$genotype[df$genotype_call == "1/."] <- "1"
        df$genotype[df$genotype_call == ".|1"] <- "1"
        df$genotype[df$genotype_call == "1/1"] <- "2"
        df$genotype[df$genotype_call == "1|1"] <- "2"
        df$genotype[df$genotype_call == "./."] <- "0"
        df$genotype[df$genotype_call == "."] <- "0"

        cat("\nunique genotypes\n")
        genotype_unique <- unique(df$genotype)

        if (all(genotype_unique %in% c("0", "1", "2"))) {
                cat("\nGenotypes found: ", paste(genotype_unique), "\n")
        } else {
                cat("\nGenotypes found: ", paste(genotype_unique))
                cat("\nError: Invalid genotype found.\nGenotypes must be format 0,1,2. See clean_genotype function for details.\n")
                stop("\nStopping analysis.\n")
        }

        # cat("Removing all genotype: 0.\n")
        df$genotype <- as.numeric(df$genotype)
        df <- df |> dplyr::filter(genotype > 0)
        if (!is.data.frame(df)) {
                stop("Unexpected object type after genotype cleaning")
        }

        cat("\nClass of df : ", class(df), "\n")
        return(df)
}

process_and_save <- function(split_file, df, output_dir, file_suffix, version, f) {
        cat("\n------------------")
        cat("\nprocess and save start")
        cat("\n------------------")

        df <- df |> distinct()
        cat("\nClass of df : ", class(df), "\n")

        # Extract the ${chromosome}_${start}_${end} from split_file
        pattern <- paste("chr", chr_names[f], "_(\\d+)_(\\d+)", sep="")
        matches <- regmatches(split_file, regexpr(pattern, split_file))
        region_info <- matches[1]  # This extracts the matched pattern from the filename

        # Update the file_path to include the region_info
        file_path <- paste(output_dir, file_suffix, "_", version, "_", region_info, ".csv", sep = "")

        cat("\nprint(output_dir, file_suffix, version, f):")
        print(output_dir)
        print(file_suffix)
        print(version)
        print(f)
        cat("\nprint(file_path):")
        print(file_path)

        cat("\ncheck for list columns:")
        list_cols <- sapply(df, is.list)
        if (any(list_cols)) {
          cat("Columns that are lists:", names(df)[list_cols], "\n")
        } else {
          cat("No columns are lists.\n")
        }

        write.csv(df, file_path)
        print(paste("chr", f, "filtered and saved to:", file_path))
}

progress_bar <- function(idx, file_list) {
        # percent_complete <- round(idx / length(file_list) * 100, 2)
        percent_complete <- round((idx / length(file_list)) * 100, 2)
        num_symbols <- round(percent_complete/2)
        cat("\nFinished processing sample: ", idx, "of ", length(file_list), "\n",
                "\r[", paste(rep("=", num_symbols), collapse = ""), ">",
                paste(rep(" ", 50 - num_symbols), collapse = ""), "]",
                percent_complete, "%\n")
}

# --------------------------------
# run analysis with functions ----
# --------------------------------

print(args)
cat("\n------------------")
cat("\nvcf_to_table start")
cat("\n------------------")
cat("\n")

# Since we launch an array of one job per chr file,
# the following analysis is run once for each.

# Get the chromosome-level file names
chr_names <- c(as.character(1:22), "X", "Y", "M")
f <- file_idx # variable f used in functions to capture the chromosome ID

# To get chromosome-level index for progress bar:
file_list <- paste0(vcf_out_path, "/bcftools_gatk_decom_maf04_chr", chr_names[f], ".recode_vep.vcf.gz")

# We have indivudal chunks per chromosome since files were split previously
split_files_pattern <- paste0("^bcftools_gatk_decom_maf04_chr", chr_names[f], "\\.recode_vep.*\\.vcf\\.gz$")
split_files_list <- list.files(path = vcf_out_path, pattern = split_files_pattern, full.names = TRUE)

print(paste0("Looking for files in:", vcf_out_path))
print(paste0("Using pattern:", split_files_pattern))
print(split_files_list)

cat("\nStarting loop on all chr files\n")
for(idx in 1:length(split_files_list)) {
        split_file <- split_files_list[idx]
        print(split_file)
        df_main <- read_and_convert_vcf(split_file)
        # df <- filter_and_gather(df_main, sample_count)
        df <- filter_and_gather(df_main, sample_count, output_dir, file_suffix, version, f, split_file)
        df <- clean_genotype(df)
        process_and_save(split_file, df, output_dir, file_suffix, version, f)
        progress_bar(idx, split_files_list)
}

sink() # stop console output to log file
~
~
~
~
~
~
