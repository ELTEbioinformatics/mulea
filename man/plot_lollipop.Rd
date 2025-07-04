% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plotting.R
\name{plot_lollipop}
\alias{plot_lollipop}
\title{Plot Lollipop}
\usage{
plot_lollipop(
  reshaped_results,
  ontology_id_colname = "ontology_id",
  plot_names_colname = NULL,
  selected_rows_to_plot = NULL,
  p_value_type_colname = "eFDR",
  p_value_max_threshold = 0.05
)
}
\arguments{
\item{reshaped_results}{data.table in relaxed form, obtained as the output
of the \code{reshape_results} function. The data source for generating the
barplot.}

\item{ontology_id_colname}{Character, specifies the column name that contains
ontology IDs in the input data.}

\item{plot_names_colname}{Column name used for representing labels of
ontologies.}

\item{selected_rows_to_plot}{A numeric vector specifying which rows of the
reshaped results data frame should be included in the plot. Default is NULL.
frame should be included in the plot?}

\item{p_value_type_colname}{Character, specifies the column name for p-values
in the input data. Default is 'eFDR'.}

\item{p_value_max_threshold}{Numeric, representing the maximum p-value
threshold for filtering data. Default is 0.05.}
}
\value{
Returns a lollipop plot
}
\description{
Plots lollipop plot of p-values.
}
\details{
Create a customized  lollipop plot of p-values, facilitating visual
exploration and analysis of statistical significance within ontology
categories.
}
\examples{
library(mulea)

# loading and filtering the example ontology from a GMT file
tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
    "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
    max_nr_of_elements = 400)

# loading the example data
sign_genes <- readLines(system.file(package = "mulea", "extdata", 
    "target_set.txt"))
background_genes <- readLines(system.file(
    package="mulea", "extdata", "background_set.txt"))

# creating the ORA model
ora_model <- ora(gmt = tf_gmt_filtered,
    # the test set variable
    element_names = sign_genes, 
    # the background set variable
    background_element_names = background_genes, 
    # the p-value adjustment method
    p_value_adjustment_method = "eFDR", 
    # the number of permutations
    number_of_permutations = 10000,
    # the number of processor threads to use
    nthreads = 2)
# running the ORA
ora_results <- run_test(ora_model)

# reshaping results for visualisation
ora_reshaped_results <- reshape_results(
    model = ora_model,
    model_results = ora_results,
    # choosing which column to use for the indication of significance
    p_value_type_colname = "eFDR")
        
# Plot lollipop
plot_lollipop(reshaped_results = ora_reshaped_results,
    # the column containing the names we wish to plot
    ontology_id_colname = "ontology_id",
    # upper threshold for the value indicating the significance
    p_value_max_threshold = 0.05,
    # column that indicates the significance values
    p_value_type_colname = "eFDR")

}
\seealso{
\code{\link{reshape_results}}
}
