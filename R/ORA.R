#'An S4 class to represent a set based tests in mulea.
#'
#'@slot method The overrepresentation (ora) method. Possible values:
#'  "Hypergeometric", "SetBasedEnrichment".
#'@slot gmt A `data.frame` representing the ontology GMT.
#'@slot element_names A vector of elements names 
#'  (gene or protein names or identifiers) representing the target set to 
#'  analyse. For example
#'  differentially expressed genes.
#'@slot background_element_names A vector of elements names 
#'  (gene or protein names or identifiers) representing all the 
#'  elements involved in the previous analyses For example all 
#'  genes that were measured in differential expression analysis.
#'@slot p_value_adjustment_method A character string representing the 
#'  type of the *p*-value adjustment method. Possible values:
#' * 'eFDR': empirical false discovery rate correction method
#' * all `method` options from `stats::p.adjust` documentation.
#'@slot number_of_permutations A numeric value representing the number of
#'  permutations used to calculate the eFDR values. Default value is 10000.
#'@slot nthreads Number of processor's threads to use in calculations.
#'@return ora object. This object represents the result of the
#'  overrepresentation test in mulea.
#'@export ora
#' @examples
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(
#'     package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'         max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(package = "mulea", "extdata", 
#'         "sign_genes.csv"))
#' background_genes <- readLines(system.file(package="mulea", "extdata", 
#'         "background_genes.csv"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered, 
#'         # the test set variable
#'         element_names = sign_genes, 
#'         # the background set variable
#'         background_element_names = background_genes, 
#'         # the p-value adjustment method
#'         p_value_adjustment_method = "eFDR", 
#'         # the number of permutations
#'         number_of_permutations = 10000,
#'         # the number of processor threads to use
#'         nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)

ora <- setClass("ora", slots = list(
    gmt = "data.frame",
    element_names = "character",
    background_element_names = "character",
    p_value_adjustment_method = "character",
    number_of_permutations = "numeric",
    nthreads = "numeric",
    test = "function"))

setMethod("initialize", "ora", function(.Object, gmt = data.frame(), 
    element_names = character(),
    background_element_names = character(),
    p_value_adjustment_method = "eFDR",
    number_of_permutations = 10000,
    test = NULL,
    nthreads = 4,
    ...) {
        adjustMethod <- NULL
        .Object@gmt <- gmt
        .Object@element_names <- element_names
        .Object@background_element_names <- background_element_names
        .Object@p_value_adjustment_method <- p_value_adjustment_method
        .Object@number_of_permutations <- number_of_permutations
        .Object@nthreads <- nthreads
        .Object@test <- function(setBasemodel) {
            setBasedTestRes <- computeSetBasedTestRes(setBasemodel)
            setBasedTestRes}
        .Object})

performSetBasedEnrichmentTest <- function(setBasemodel) {
    muleaSetBaseEnrichmentTest <- SetBasedEnrichmentTest(gmt = setBasemodel@gmt,
        element_names = setBasemodel@element_names,
        pool = setBasemodel@background_element_names,
        number_of_permutations = setBasemodel@number_of_permutations,
        nthreads = setBasemodel@nthreads)
    muleaSetBaseEnrichmentTest <- run_test(muleaSetBaseEnrichmentTest)
    muleaSetBaseEnrichmentTest <- merge(setBasemodel@gmt[c("ontology_id",
        "ontology_name")], muleaSetBaseEnrichmentTest, by.x = "ontology_id",
        by.y = "DB_names", all = TRUE)
    return(muleaSetBaseEnrichmentTest)}

performSetBasedEnrichmentTestElse <- function(setBasemodel) {
    MuleaHypergeometricTest <- MuleaHypergeometricTest(gmt = setBasemodel@gmt,
        element_names = setBasemodel@element_names,
        pool = setBasemodel@background_element_names,
        nthreads = setBasemodel@nthreads)
    setBasedTestRes <- run_test(MuleaHypergeometricTest)
    muleaSetBaseEnrichmentTest <- merge(setBasemodel@gmt[c("ontology_id",
        "ontology_name")], setBasedTestRes, by.x = "ontology_id",
        by.y = "ontology_name", all = TRUE)
    res <- list("setBasedTestRes"=setBasedTestRes,
                "muleaSetBaseEnrichmentTest"=muleaSetBaseEnrichmentTest)
    return(res)
}

computeSetBasedTestRes <- function(setBasemodel) {setBasedTestRes <- NULL
    if (!identical(setBasemodel@p_value_adjustment_method, character(0)) &&
        setBasemodel@p_value_adjustment_method == "eFDR") {
        muleaSetBaseEnrichmentTest<-performSetBasedEnrichmentTest(setBasemodel)
        for (i in seq_along(muleaSetBaseEnrichmentTest$FDR)) {
            if (!is.nan(muleaSetBaseEnrichmentTest$FDR[i])
                && muleaSetBaseEnrichmentTest$FDR[i] > 1.0) {
                muleaSetBaseEnrichmentTest$FDR[i] <- 1.0e+00}}
            names(muleaSetBaseEnrichmentTest) <- c("ontology_id", 
                "ontology_name", "nr_common_with_tested_elements", 
                "nr_common_with_background_elements", "Genes_in_DB", "p_value", 
                "P_adj_Bonf", "adjustedPValue", "R_obs", "R_exp", "eFDR")
            setBasedTestRes <- muleaSetBaseEnrichmentTest[, 
                !names(muleaSetBaseEnrichmentTest) %in% c("Genes_in_DB",
                    "P_adj_Bonf", "R_obs", "R_exp", "adjustedPValue")]
            } else { res <- performSetBasedEnrichmentTestElse(setBasemodel)
                muleaSetBaseEnrichmentTest <- res$muleaSetBaseEnrichmentTest
                setBasedTestRes <- res$setBasedTestRes
                names(muleaSetBaseEnrichmentTest) <- c("ontology_id",
                    "ontology_name", "list_of_values", "p_value")
                if (!identical(setBasemodel@p_value_adjustment_method, 
                    character(0)) && 
                        setBasemodel@p_value_adjustment_method != "eFDR") {
                        muleaSetBaseEnrichmentTest <- data.frame(
                            muleaSetBaseEnrichmentTest, 
                            "adjusted_p_value" = stats::p.adjust(
                            muleaSetBaseEnrichmentTest$p_value,
                            method = setBasemodel@p_value_adjustment_method))
                        setBasedTestRes <- muleaSetBaseEnrichmentTest[, 
                            !names(muleaSetBaseEnrichmentTest) %in% c(
                                "list_of_values")]}}
    return(setBasedTestRes)
}

#' @describeIn run_test ora test.
#' @param model Object of S4 class representing the mulea test.
#' @return run_test method for ora object. Returns the results of the 
#' overrepresentation analysis.
#' @examples
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(
#'         package="mulea", "extdata", 
#'         "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'         max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(package = "mulea", "extdata", 
#'         "sign_genes.csv"))
#' background_genes <- readLines(system.file(package="mulea", "extdata", 
#'         "background_genes.csv"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered, 
#'         # the test set variable
#'         element_names = sign_genes, 
#'         # the background set variable
#'         background_element_names = background_genes, 
#'         # the p-value adjustment method
#'         p_value_adjustment_method = "eFDR", 
#'         # the number of permutations
#'         number_of_permutations = 10000,
#'         # the number of processor threads to use
#'         nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)

setMethod("run_test", signature(model = "ora"), function(model) {
    model@test(model)})
