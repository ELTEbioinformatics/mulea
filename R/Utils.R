#' Read GMT File
#' @description Reads gene set or ontology data from a Gene Matrix Transposed
#' (GMT) file and parse into a `data.frame`.
#'
#' @param file Character, a path to a file.
#' @return Returns a `data.frame` with three columns:
#'
#' * 'ontology_id': Ontology identifier that uniquely identifies the element 
#'     within the referenced ontology.
#' * 'ontology_name': Ontology name or description that provides a
#'     user-friendly label or textual description for the 'ontology_id'.
#' * 'list_of_values': Associated genes or proteins that is a vector of
#'     identifiers of genes or proteins belonging to the 'ontology_id'.
#'
#' @export
#'
#' @examples
#' # import example gene set
#' library(mulea)
#' tf_gmt <- read_gmt(file = system.file(
#'     package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))

read_gmt <- function(file) { fileConnection <- file(file)
    tryCatchRes <- tryCatch(
        lines <- readLines(fileConnection),
        warning = function(w) { })
    close(fileConnection)
    lines <- lines[!grepl('^#+', lines, fixed = FALSE)]
    lines <- lines["" != lines]
    gmtAsDF <- plyr::adply(.data = lines, .margins = 1, .fun = function(line) {
        fields <- strsplit(line, split = "\t")[[1]]
        category <- fields[1]
        description <- fields[2]
        list_of_values <- fields[3:length(fields)]
        data.frame("ontology_id" = category, "ontology_name" = description,
            "list_of_values" = I(list(list_of_values)), 
            stringsAsFactors = FALSE)})
    gmtAsDF[c("ontology_id", "ontology_name", "list_of_values")]
}


#' Write GMT file
#' @description Writes gene set or ontology `data.frame` with specific
#' formatting (columns representing ontology identifiers, descriptions, and
#' associated lists of values) and writes it to a file in a standardized Gene 
#' Matrix Transposed (GMT) file format.
#'
#' @param gmt A `data.frame` containing the data to be written, 
#'   imported from a GMT file with the `read_gmt` function.
#' @param file Character, a path to a file.
#' @return Returns the input as a GMT file at a specific location.
#' @export
#' @examples
#' library(mulea)
#'
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(
#'     package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#'
#' # writing the filtered ontology to a GMT file
#' \dontrun{
#' write_gmt(
#'     gmt = tf_gmt, 
#'     file = "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt")
#' }

write_gmt <- function(gmt, file) {
    vectorOfModel <- plyr::daply(.data = gmt, .variables = c("ontology_id"),
        .fun = function(dataFrameRow) {
        collapsedlist_of_values <- paste( dataFrameRow[, 3][[1]], 
            collapse = "\t")
        paste(dataFrameRow[1], dataFrameRow[2], collapsedlist_of_values, 
            sep = "\t")})
    fileConnection <- file(file)
    writeLines(vectorOfModel, con = fileConnection, sep = "\n", 
        useBytes = FALSE)
    close(fileConnection)
}

#' Filter Ontology
#' @description Filtering ontology to contain entries having number of elements
#'   (genes or proteins) between a given range. The reason for this is
#'   enrichment analysis results can sometimes be skewed by overly specific or
#'   broad entries. Filtering ontologies allows you to customize the size of
#'   ontology entries, ensuring your analysis aligns with your desired scope.
#' @param gmt A `data.frame` which contains the entries 
#'   (gene or protein sets), imported from a GMT file with the 
#'   `read_gmt` function.
#' @param min_nr_of_elements Minimum number of elements. Ontology entries
#'   containing as many or fewer elements (genes or proteins) will be excluded.
#' @param max_nr_of_elements Maximum number of elements. Ontology entries
#'   containing as many or more elements (genes or proteins) will be excluded.
#' @return Return a `data.frame`which contains the entries (gene or protein
#'   sets) in a similar format that produced by the `read_gmt` function.
#' @export
#' @examples
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(
#'     package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt,
#'         min_nr_of_elements = 3,
#'         max_nr_of_elements = 400)

filter_ontology <- function(gmt, min_nr_of_elements = NULL,
    max_nr_of_elements = NULL) {
    if (is.null(min_nr_of_elements)) {
        terms_sizes <- plyr::laply(.data = gmt$list_of_values,
            .fun = function(term) { length(term) })
            min_nr_of_elements <- 3}
    if (is.null(max_nr_of_elements)) {
        terms_sizes <- plyr::laply(.data = gmt$list_of_values,
            .fun = function(term) { length(term) })
            max_nr_of_elements <- 400}
    filtered_input_gmt <- plyr::ddply(.data = gmt, 
        .variables = c("ontology_id"),
        .fun = function(df_row) {
            if (length(df_row$list_of_values[[1]]) > min_nr_of_elements) {
                df_row } else { df_row[-1, ] }})
    filtered_input_gmt <- plyr::ddply(.data = filtered_input_gmt,
        .variables = c("ontology_id"),
        .fun = function(df_row) {
            if (length(df_row$list_of_values[[1]]) < max_nr_of_elements) {
                df_row } else { df_row[-1, ]}})
    filtered_input_gmt
}


#' Convert a list to ontology (GMT) `data.frame`.
#' @description Converts a list of ontology elements (gene sets) to an ontology
#' (GMT) `data.frame` object.
#' @param gmt_list A list with named character vectors. The name will become the
#'   'ontology_id', and the elements in the vector will become the
#'   'list_of_values' in the ontology (GMT) `data.frame`.
#' @examples
#' library(mulea)
#' 
#' # creating a list of gene sets
#' ontology_list <- list(gene_set1 = c("gene1", "gene2", "gene3"),
#'         gene_set2 = c("gene4", "gene5", "gene6"))
#'
#' # converting the list to a ontology (GMT) object
#' new_ontology_object <- list_to_gmt(ontology_list)
#' @return Returns ontology (GMT) `data.frame` where the 'ontology_name'
#'   contains random unique strings.
#' @export
list_to_gmt <- function(gmt_list) {
    listAsGmtDataFrame <- plyr::ldply(.data = gmt_list, .id = c('ontology_id'),
        .fun = function(element) {
        ontology_name <- stringi::stri_rand_strings(length = 5, n = 1)
        data.frame("ontology_name" = ontology_name, 
            "list_of_values" = I(list(element)), stringsAsFactors = FALSE)})
    return(listAsGmtDataFrame)
}


