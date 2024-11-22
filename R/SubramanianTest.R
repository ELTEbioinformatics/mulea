#' PRIVATE class : An S4 class to represent a ranked based tests in mulea.
#'
#' @slot gmt A `data.frame` representing the ontology GMT.
#' @slot element_names A vector of elements names (gene or protein names or
#'   identifiers) to include in the analysis.
#' @slot element_scores A vector of numeric values representing a 
#'   score (*e.g.* *p*-value, *z*-score, log fold change) for each 
#'   'element_name', in the same number and order as element_name.
#' @slot p A power of weight.
#' @slot element_score_type Defines the GSEA score type.
#' * "pos": Only positive element_scores
#' * "neg": only negative element_scores - "neg" and mixed
#' * "std": standard -- containing both positive and negative scores
#'   Default value is "std".
#' @return `data.frame` with presented columns 'ontology_id', 'ontology_name',
#' 'nr_common_with_tested_elements',
#' 'p_value', 'adjusted_p_value'
SubramanianTest <- setClass(
    "SubramanianTest",
    slots = list(
        gmt = "data.frame",
        element_names = "character",
        element_scores = "numeric",
        gsea_power = "numeric",
        element_score_type = "character",
        test = "function"))

#' @importFrom fgsea fgsea
setMethod("initialize", "SubramanianTest",
    function(.Object,
        gmt = data.frame(),
        element_names = character(),
        element_scores = numeric(),
        gsea_power = 1,
        element_score_type = "std",
        test = NULL,
        ...) {
            .Object@gmt <- gmt
            .Object@element_names <- element_names
            .Object@element_scores <- element_scores
            .Object@gsea_power <- gsea_power
            .Object@element_score_type <- element_score_type
            .Object@test <- function(model) {
                listmodelDfFromFile <- lapply(model@gmt$list_of_values, 
                                              strsplit, split=" ", fixed=T)
                names(listmodelDfFromFile) <- model@gmt$ontology_id
                samplesToAnalisys <- model@element_scores
                names(samplesToAnalisys) <- model@element_names
                fgseaRes <- fgsea::fgsea(pathways = listmodelDfFromFile,
                    stats = samplesToAnalisys,
                    gseaParam = model@gsea_power,
                    scoreType = model@element_score_type)
                resultDf <- merge(model@gmt, fgseaRes, by.x = "ontology_id",
                    by.y = "pathway", all = TRUE)
                resultDf <- plyr::ddply(.data = resultDf,
                    .variables = c('ontology_id'),
                    .fun = function(df_row) {
                        nr_common_genes_ontology_background <-
                            length(df_row[, 'list_of_values'][[1]])
                        cbind(df_row, nr_common_genes_ontology_background)
                    })[c("ontology_id", "ontology_name", "size", "pval", 
                        "padj")]
                colnames(resultDf) <- c("ontology_id", "ontology_name",
                    "nr_common_with_tested_elements", "p_value",
                    "adjusted_p_value")
                resultDf
                }
            .Object
})

#' @describeIn SubramanianTest runs test calculations.
#' @param model Object of s4 class represents mulea Test.
#' @return run_test method for SubramanianTest object. Used as private 
#'   function.
setMethod("run_test",
    signature(model = "SubramanianTest"),
    function(model) {
        model@test(model)})
