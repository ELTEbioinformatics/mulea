################################################
# Function for FDR corrected hypergeometric enrichment test
################################################

#' PRIVATE class : An S4 class to represent a Hypergeometric tests in 
#'   mulea.
#'
#' @slot gmt A data.frame representing GMT's representation of model.
#' @slot element_names A data from experiment to analyse across model.
#' @slot pool A background data to count test.
#' @slot nthreads Number of processor's threads used in calculations.
#' @return SetBasedEnrichmentTest object. Used as private function.
SetBasedEnrichmentTest <- setClass(
    "SetBasedEnrichmentTest",
    slots = list(
        gmt = "data.frame",
        element_names = "character",
        pool = "character",
        number_of_permutations = "numeric",
        only_hyper_geometric_test = "logical",
        nthreads = "numeric",
        test = "function"))

setMethod("initialize", "SetBasedEnrichmentTest",
    function(.Object,
        gmt = data.frame(),
        element_names = character(),
        pool = character(),
        number_of_permutations = 10000,
        test = NULL,
        only_hyper_geometric_test = FALSE,
        nthreads = 2,
        ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@pool <- pool
            .Object@number_of_permutations <- number_of_permutations
            .Object@only_hyper_geometric_test <- only_hyper_geometric_test
            .Object@nthreads <- nthreads
            .Object@test <- function(model) {
                pool <- NULL
                if (0 == length(model@pool)) {
                    pool <- unique(unlist(.Object@gmt[, 'list_of_values']))
                } else {
                    pool <- unique(model@pool)
                }
                element_names <- .Object@element_names
                if (!all(element_names %in% pool)) {
                    element_names <- element_names[element_names %in% pool]
                    warning("Not all elements of element_names (sample) 
                        are from pool.", "TestData vector is automatically 
                        cut off to pool vector.")
                }
                DB <- .Object@gmt[, 'list_of_values']
                names(DB) <- .Object@gmt$ontology_id
                testResults <- set.based.enrichment.test.wrapper(
                    steps = .Object@number_of_permutations,
                    pool = pool, select = element_names, DB = DB,
                    only_hyper_geometric_test = 
                      model@only_hyper_geometric_test,
                    nthreads = model@nthreads)
                testResults
            }
            .Object
})

#' @describeIn SetBasedEnrichmentTest runs test calculations.
#' @param model Object of s4 class represents mulea Test.
#' @return run_test method for SetBasedEnrichmentTest object. Used as private 
#' function.
setMethod("run_test",
    signature(model = "SetBasedEnrichmentTest"),
    function(model) {
        model@test(model)
    })


set.based.enrichment.test.wrapper <- function(
    steps, pool, select, DB, nthreads = 2, 
    only_hyper_geometric_test = FALSE) {

    setEnrTestRes <- set.based.enrichment.test(
        steps = steps,
        pool = pool,
        select = select,
        DB = DB,
        nthread = nthreads)
    return(setEnrTestRes)
}
