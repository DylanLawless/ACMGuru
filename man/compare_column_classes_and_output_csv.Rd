% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_column_classes_and_output_csv.R
\name{compare_column_classes_and_output_csv}
\alias{compare_column_classes_and_output_csv}
\title{Compare Column Classes and Output CSV}
\usage{
compare_column_classes_and_output_csv(
  processed_data_list,
  metadataclass_file_path = NULL
)
}
\arguments{
\item{processed_data_list}{A list of dataframes to be analyzed.}

\item{metadataclass_file_path}{Optional. Path where the output CSV file will be saved.
If not provided, uses the default file included in the package.}
}
\description{
Compares column classes across all dataframes in a list and outputs a CSV file with the most common class for each column where inconsistencies are found.
}
