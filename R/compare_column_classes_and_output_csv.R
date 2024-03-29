#' Compare Column Classes and Output CSV
#'
#' Compares column classes across all dataframes in a list and outputs a CSV file with the most common class for each column where inconsistencies are found.
#'
#' @param processed_data_list A list of dataframes to be analyzed.
#' @param metadataclass_file_path Optional. Path where the output CSV file will be saved.
#' If not provided, uses the default file included in the package.
#' @export
compare_column_classes_and_output_csv <- function(processed_data_list, metadataclass_file_path = NULL) {
  # Use the default file if no path is provided
  # if (is.null(metadataclass_file_path)) {
  #   metadataclass_file_path <- system.file("extdata", "suggested_reference_metadata.csv", package = "ACMGuru")
  #   if (metadataclass_file_path == "") {
  #     stop("Default metadata file not found in package 'ACMGuru'. Please provide a valid metadataclass_file_path.")
  #   }
  # }

  # Use the default file if no path is provided
  if (is.null(metadataclass_file_path)) {
    metadataclass_file_path <- system.file("extdata", "suggested_reference_metadata.csv", package = "ACMGuru")
    if (metadataclass_file_path == "") {
      stop("Default metadata file not found in package 'ACMGuru'. Please provide a valid metadataclass_file_path. (Function: compare_column_classes_and_output_csv)")
    } else {
      cat("Using default metadata file located at:", metadataclass_file_path, "\n (Function: compare_column_classes_and_output_csv).")
    }
  }


  column_classes_list <- lapply(processed_data_list, function(df) sapply(df, class))
  all_column_names <- unique(unlist(lapply(column_classes_list, names)))

  suggested_reference_metadata <- data.frame(column_name = character(), most_common_class = character(), stringsAsFactors = FALSE)

  for (column_name in all_column_names) {
    column_classes <- sapply(column_classes_list, `[`, column_name)
    unique_classes <- unique(na.omit(column_classes))

    if (length(unique_classes) > 1) {
      class_counts <- table(unlist(column_classes))
      most_common_class <- names(which.max(class_counts))

      suggested_reference_metadata <- rbind(suggested_reference_metadata, data.frame(column_name = column_name, most_common_class = most_common_class, stringsAsFactors = FALSE))
    }
  }

  if(nrow(suggested_reference_metadata) > 0) {
    write.csv(suggested_reference_metadata, metadataclass_file_path, row.names = FALSE)
    cat("Suggested reference metadata classes have been output to:", metadataclass_file_path, "\n")
  } else {
    cat("All columns have consistent classes across dataframes. No output file created.\n")
  }
}
