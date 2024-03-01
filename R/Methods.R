checkIfPoolIncludeSample <-
    function(model, sampleVector, poolVector = NULL) {
        # Chcking if experiment data all in model data.
        if (0 != length(poolVector)) {
            if (0 != sum(!(sampleVector %in% poolVector))) {
                to_msg <- paste(setdiff(sampleVector, unique(poolVector)), 
                                collapse = ", ")
                warning("element_names are outside of pool.", " ", to_msg)
        return(setdiff(sampleVector, 
                       setdiff(sampleVector, unique(poolVector))))
        }
    } else {
        if (0 != sum(!(sampleVector %in% unique(unlist(
            model$list_of_values
            ))))) {
          to_msg <- paste(setdiff(sampleVector, unique(
            unlist(model$list_of_values))), collapse = ", ")
        warning("element_names are outside of gmt.", " ", to_msg)
        return(setdiff(sampleVector, setdiff(sampleVector, unique(
            unlist(model$list_of_values)
        ))))
            }
        }
    return(sampleVector)
    }

cutGmtToPool <- function(gmt, pool) {
    cutDF <- plyr::ddply(
        .data = gmt,
        .variables = c("ontology_id"),
        .fun = function(dfRow) {
            dfRow$list_of_values[[1]] <- intersect(
                dfRow$list_of_values[[1]], pool)
        dfRow
        }
    )
    cutDF
}
