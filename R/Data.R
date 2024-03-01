#' Ontology Example
#'
#' Data frame parsed from a Gene Matrix Transposed (GMT) 
#' file containing 7 gene
#' sets. Each row represents a gene set. 
#'
#' @format A data frame with 7 rows and 3 variables:
#' \describe{
#'   \item{ontology_id}{Ontology ID of the gene set}
#'   \item{ontology_name}{Description of the gene set}
#'   \item{list_of_values}{Genes in the gene set}
#' }
#' @source \url{http://geneontology.org/}
"geneSet"

#' Genes and scores for ora and gsea analyses
#'
#' Data frame containing 13 genes in the "select" column 
#' and their scores in the
#' "score" column.
#'
#' @format A data frame with 13 rows and 3 variables:
#' \describe{
#'   \item{X}{Row names}
#'   \item{select}{Flybase gene IDs of 13 genes from Drosophyla melanogaster}
#'   \item{score}{Log fold change scores for each gene}
#' }
#' @source \url{https://flybase.org/}
"selectDf"

#' Background genes for ora analyses
#'
#' Data frame containing 13 genes in the "background_element_names" column.
#'
#' @format A data frame with 33 rows and 2 variables:
#' \describe{
#'   \item{X}{Row names}
#'   \item{background_element_names}{Flybase gene IDs of 33 genes from 
#'   Drosophyla melanogaster}
#' }
#' @source \url{https://flybase.org/}
"poolDf"

