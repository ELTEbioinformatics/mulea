#' Gene Set Enrichment Analysis (GSEA)
#'
#' An S4 class to represent the gsea tests in mulea.
#'
#' @slot gmt  A `data.frame` representing the ontology GMT.
#' @slot element_names A vector of elements names 
#'   (gene or protein names or identifiers) to include in the analysis.
#' @slot element_scores A vector of numeric values representing a 
#'   score (*e.g.* *p*-value, *z*-score, log fold change) for each 
#'  'element_name', in the same number and order as element_name.
#' @slot gsea_power A power of weight. Default value is 1.
#' @slot element_score_type Defines the GSEA score type.
#' * 'pos': Only positive element_scores
#' * 'neg': Only negative element_scores
#' * 'std': standard, containing both positive and negative scores
#'   Default value is 'std'.
#' @slot number_of_permutations The number of permutations used in 
#'   `gsea` test. Default value is 1000.
#' @slot test character
#' @return GSEA object. This object represents the result of the 
#'   `gsea` tests.
#' @export
#' @examples
#' library(mulea)
#'
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(
#'     package="mulea", "extdata",
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'     max_nr_of_elements = 400)
#'
#' # loading the example `data.frame`
#' scored_gene_tab <- read.delim(file = system.file(package = "mulea", "extdata", 
#'         "ordered_set.tsv"))
#'
#' # creating the GSEA model
#' gsea_model <- gsea(gmt = tf_gmt_filtered,
#'         # the names of elements to test
#'         element_names = scored_gene_tab$Gene.symbol,
#'         # the logFC-s of elements to test
#'         element_scores = scored_gene_tab$logFC,
#'         # consider elements having positive logFC values only
#'         element_score_type = "pos",
#'         # the number of permutations
#'         number_of_permutations = 10000)

gsea <- setClass(
    "gsea",
    slots = list(
        gmt = "data.frame",
        element_names = "character",
        element_scores = "numeric",
        gsea_power = "numeric",
        element_score_type = "character",
        number_of_permutations = "numeric",
        test = "function"))

setMethod("initialize", "gsea",
    function(.Object,
        gmt = data.frame(),
        element_names = character(),
        element_scores = numeric(),
        gsea_power = 1,
        element_score_type = "std",
        number_of_permutations = 1000,
        test = NULL,
        ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@element_scores <- element_scores
            .Object@gsea_power <- gsea_power
            .Object@element_score_type <- element_score_type
            .Object@number_of_permutations <- number_of_permutations
            .Object@test <- function(rankedBasemodel) {
                rankedTestRes <- NULL
                subramanianTest <- SubramanianTest(
                    gmt = rankedBasemodel@gmt,
                    element_names = rankedBasemodel@element_names,
                    element_scores = rankedBasemodel@element_scores,
                    gsea_power = rankedBasemodel@gsea_power,
                    element_score_type = rankedBasemodel@element_score_type)
                rankedTestRes <- run_test(subramanianTest)
            rankedTestRes
            }
            .Object
})

#' @describeIn gsea runs test calculations.
#' @param model Object of S4 class representing the mulea test.
#' @return run_test method for GSEA object. Returns results of
#' the enrichment analysis.
#' @examples
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#'         "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'         max_nr_of_elements = 400)
#' 
#' # loading the example `data.frame`
#' scored_gene_tab <- read.delim(file = system.file(package = "mulea", "extdata", 
#'         "ordered_set.tsv"))
#'
#' # creating the GSEA model
#' gsea_model <- gsea(gmt = tf_gmt_filtered,
#'         # the names of elements to test
#'         element_names = scored_gene_tab$Gene.symbol,
#'         # the logFC-s of elements to test
#'         element_scores = scored_gene_tab$logFC,
#'         # consider elements having positive logFC values only
#'         element_score_type = "pos",
#'         # the number of permutations
#'         number_of_permutations = 10000)
#'
#' # running the test
#' gsea_results <- run_test(gsea_model)

setMethod("run_test",
    signature(model = "gsea"),
    function(model) {
        model@test(model)})
