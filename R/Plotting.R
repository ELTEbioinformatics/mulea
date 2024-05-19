validate_column_names_and_function_args <- function(data, ...) {
    arguments <- list(...)
    if (!all(unlist(arguments) %in% names(data))) {
        stop('Wrongly set data column names.')
    }
}

filterRelaxedResultsForPlotting <- function(reshaped_results,
    p_value_type_colname = 'ontologyStatValue',
    p_value_max_threshold = 0.05) {
        include <- !is.na(reshaped_results[[p_value_type_colname]])
    reshaped_results_filtered_na <- reshaped_results[include,]
    include <- reshaped_results_filtered_na[[
        p_value_type_colname]] <= p_value_max_threshold
    reshaped_results_filtered_cutoff <- 
        reshaped_results_filtered_na[include, ]
    reshaped_results_filtered_cutoff
}

#' Reshape Results
#' @description This function takes a model and model_results data, reshapes
#'   them into a suitable format for plotting, and returns the resulting data
#'   frame, which can be used for further analysis or visualization.
#'
#' @param model a mulea model, created by the 
#'   `ora` or the `gsea` functions.
#' @param model_results Result `data.frame` 
#'   returned by the `run_test` function.
#' @param model_ontology_col_name Character, specifies the column name in the
#'   model that contains ontology IDs. It defines which column in the model
#'   should be used for matching ontology IDs. Possible values are 'ontology_id'
#'   and 'ontology_name'. The default value is 'ontology_id'.
#' @param ontology_id_colname Character, specifies the column name for ontology
#'   IDs in the model results. It indicates which column in the model results
#'   contains ontology IDs for merging. Possible values are 'ontology_id' and
#'   'ontology_name'. The default value is 'ontology_id'.
#' @param p_value_type_colname Character, specifies the column name 
#'   for the type or raw or adjusted *p*-value in the result 
#'   `data.frame` returned by the `run_test` function. 
#'   The default value is 'eFDR'.
#' @param p_value_max_threshold Logical, indicating whether to apply a 
#'   *p*-value threshold when filtering the resulting data. If TRUE, 
#'   the function filters the data based on a *p*-value threshold.
#' @seealso \code{\link{plot_graph}}, \code{\link{plot_barplot}},
#'   \code{\link{plot_heatmap}}
#' @importFrom data.table :=
#' @import tidyverse
#' @import magrittr
#' @export
#'
#' @return Return detailed and relaxed `data.table` where model and results are
#'   merged for plotting purposes.
#'
#' @examples
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'     max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(
#'     package = "mulea", "extdata", "target_set.txt"))
#' background_genes <- readLines(
#'     system.file(package="mulea", "extdata", "background_set.txt"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered,
#'     # the test set variable
#'     element_names = sign_genes, 
#'     # the background set variable
#'     background_element_names = background_genes, 
#'     # the p-value adjustment method
#'     p_value_adjustment_method = "eFDR", 
#'     # the number of permutations
#'     number_of_permutations = 10000,
#'     # the number of processor threads to use
#'     nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)
#' 
#' # reshaping results for visualisation
#' ora_reshaped_results <- reshape_results(model = ora_model,
#'     model_results = ora_results,
#'     # choosing which column to use for the indication of significance
#'     p_value_type_colname = "eFDR")

reshape_results <- function(model = NULL, model_results = NULL,
    model_ontology_col_name = 'ontology_id', 
    ontology_id_colname = 'ontology_id',
    p_value_type_colname = 'eFDR', p_value_max_threshold = TRUE) {
    genIdInOntology <- NULL
    model_with_res <- merge(x = model@gmt, y = model_results,
        by.x = model_ontology_col_name, by.y = ontology_id_colname, all = TRUE)
    model_with_res_dt <- data.table::setDT(model_with_res)
    model_with_res_dt_size <- 0
    for (i in seq_len(nrow(model_with_res_dt))) { model_with_res_dt_size <-
        model_with_res_dt_size + length(model_with_res_dt[[i, 
                'list_of_values']])}
    model_with_res_dt_relaxed <- data.table::data.table(
        ontology_id = rep('a', length.out = model_with_res_dt_size),
        genIdInOntology = rep('a',length.out = model_with_res_dt_size),
        ontologyStatValue = rep(1.0, length.out = model_with_res_dt_size))
    p_value_type_colname <- match.arg(p_value_type_colname,
        choices = names(model_with_res_dt))
    model_with_res_dt_relaxed_counter <- 1
    for (i in seq_len(nrow(model_with_res_dt))) {
        category_name <- model_with_res_dt[[i, 'ontology_id']]
        category_p_stat <- model_with_res_dt[[i, p_value_type_colname]]
        for (item_name in model_with_res_dt[[i, 'list_of_values']]) {
            model_with_res_dt_relaxed[model_with_res_dt_relaxed_counter,
                c("ontology_id", "genIdInOntology", 
                    "ontologyStatValue") := list(category_name, item_name, 
                        category_p_stat)]
            model_with_res_dt_relaxed_counter  <- 
                model_with_res_dt_relaxed_counter + 1}}
    if (p_value_max_threshold) {
        model_with_res_dt_relaxed <- model_with_res_dt_relaxed[
            genIdInOntology %in% model@element_names]}
    names(model_with_res_dt_relaxed) <- c("ontology_id", 
        "element_id_in_ontology", p_value_type_colname)
    model_with_res_dt_relaxed
}


#' Plot Graph (Network)
#'
#' @description Plots graph representation of enrichment results.
#'
#' @details This function generates a graph (network) visualization of the
#'   enriched ontology entries. On the plot each node represents an ontology
#'   entry below a given *p*-value threshold, and is coloured based on its
#'   significance level. A connection (edge) is drawn between two nodes if they
#'   share at least one common element (gene) belonging to the target set -- in
#'   the case of ORA results -- or all analysed elements -- in the case of GSEA
#'   results.
#' @param reshaped_results Character, the input `data.table` containing the
#'   reshaped results.
#' @param shared_elements_min_threshold Numeric, threshold specifying the
#'   minimum number of shared elements required between two ontologies to
#'   consider them connected by an edge on the graph. Default value is 0.
#' @param p_value_type_colname Character, the name of the column in the reshaped
#'   results that contains the type of *p*-values associated with the ontology
#'   elements. Default value is 'eFDR'.
#' @param ontology_id_colname Character, the name of the column in the reshaped
#'   results that contains ontology identifiers or names. Default value is
#'   'ontology_id'.
#' @param ontology_element_colname Character, the name of the column in the
#'   reshaped results that contains element identifiers within the ontology.
#'   Default value is 'element_id_in_ontology'.
#' @param p_value_max_threshold Numeric, a threshold value for filtering 
#'   rows in the reshaped results based on the 
#'   *p*-values. Rows with *p*-values greater
#'   than this threshold will be filtered out. 
#'   Default value is 0.05.
#' @return Returns a graph plot.
#' @importFrom data.table :=
#' @importFrom rlang .data
#' @seealso \code{\link{reshape_results}}
#' @export
#'
#' @examples
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'     max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(package = "mulea", "extdata", 
#'     "target_set.txt"))
#' background_genes <- readLines(system.file(
#'     package="mulea", "extdata", "background_set.txt"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered,
#'     # the test set variable
#'     element_names = sign_genes, 
#'     # the background set variable
#'     background_element_names = background_genes, 
#'     # the p-value adjustment method
#'     p_value_adjustment_method = "eFDR", 
#'     # the number of permutations
#'     number_of_permutations = 10000,
#'     # the number of processor threads to use
#'     nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)
#' 
#' # reshaping results for visualisation
#' ora_reshaped_results <- reshape_results(model = ora_model,
#'     model_results = ora_results,
#'     # choosing which column to use for the indication of significance
#'     p_value_type_colname = "eFDR")
#'
#' # Plot graph
#' plot_graph(reshaped_results = ora_reshaped_results,
#'     # the column containing the names we wish to plot
#'     ontology_id_colname = "ontology_id",
#'     # upper threshold for the value indicating the significance
#'     p_value_max_threshold = 0.05,
#'     # column that indicates the significance values
#'     p_value_type_colname = "eFDR")

plot_graph <- function(reshaped_results, ontology_id_colname='ontology_id',
    ontology_element_colname='element_id_in_ontology',
    shared_elements_min_threshold = 0,
    p_value_type_colname = 'eFDR', p_value_max_threshold = 0.05) {
    ontology_id <- ontology_id_colname
    edges <- NULL
    validate_column_names_and_function_args(
        data = reshaped_results, p_value_type_colname, 
        ontology_id_colname, ontology_element_colname)
    reshaped_results <- data.table::setDT(reshaped_results)
    model_with_res_dt_relaxed <- filterRelaxedResultsForPlotting(
        reshaped_results = reshaped_results,
        p_value_type_colname = p_value_type_colname,
        p_value_max_threshold = p_value_max_threshold)
    ontologies <- unique(
        model_with_res_dt_relaxed[, ontology_id_colname, with = FALSE])
    ontologies_graph_edges <- getOntologiesGraphEdgesTemplate(ontologies)
    ontologies_graph_edges <- fillOntologiesGrapgEdge(
        ontologies = ontologies, ontology_id = ontology_id, 
        model_with_res_dt_relaxed = model_with_res_dt_relaxed, 
        ontology_element_colname = ontology_element_colname, 
        ontologies_graph_edges = ontologies_graph_edges,
        shared_elements_min_threshold = shared_elements_min_threshold)
    nodes_ids <- model_with_res_dt_relaxed[[ontology_id]]
    nodes_p_stat <- model_with_res_dt_relaxed[[p_value_type_colname]]
    ontologies_graph_nodes <- data.table::data.table(id = nodes_ids,
        label = nodes_ids, p_stat = nodes_p_stat)
    ontologies_graph_nodes <- unique(ontologies_graph_nodes)
    routes_tidy <- tidygraph::tbl_graph(nodes = ontologies_graph_nodes,
        edges = ontologies_graph_edges, directed = TRUE)
    graph_plot <- ggraph::ggraph(routes_tidy, layout = "linear", 
        circular = TRUE)
    if (0 != nrow(tibble::as_tibble(tidygraph::activate(routes_tidy, 
        edges)))) {
        graph_plot <- graph_plot + 
        ggraph::geom_edge_arc(aes(width = .data$weight), alpha = 0.5)}
    graph_plot <- addExtraParamsToPlt(graph_plot=graph_plot, 
        p_value_max_threshold=p_value_max_threshold,
        p_value_type_colname=p_value_type_colname)
    graph_plot
}

getOntologiesGraphEdgesTemplate <- function(ontologies){
    ontologies_graph_edges_num <- sum(seq_len(nrow(ontologies) - 1))
    ontologies_graph_edges <- data.table::data.table(
        from = rep('a', length.out = ontologies_graph_edges_num),
        to = rep('a', length.out = ontologies_graph_edges_num),
        weight = rep(0, length.out = ontologies_graph_edges_num))
    if (0 == ontologies_graph_edges_num) {
        stop("No edges at all. Wrong data.table or manipulate 
            p_value_max_threshold please.")}
    return(ontologies_graph_edges)
}

fillOntologiesGrapgEdge <- function(ontologies, ontology_id, 
    model_with_res_dt_relaxed, ontology_element_colname, 
    ontologies_graph_edges, shared_elements_min_threshold) {
    ontologies_graph_edges_counter <- 1
    for (i in seq_len(nrow(ontologies) - 1)) {
        ontology_name_i <- ontologies[[i, ontology_id]]
        filter_model_row_select <- model_with_res_dt_relaxed[[
            ontology_id]] == ontology_name_i
        genes_in_ontology_i <- model_with_res_dt_relaxed[
            filter_model_row_select][[ontology_element_colname]]
        for (j in (i + 1):nrow(ontologies)) {
            ontology_name_j <- ontologies[[j, ontology_id]]
            filter_model_row_select_j <- model_with_res_dt_relaxed[[
                ontology_id]] == ontology_name_j
            genes_in_ontology_j <- model_with_res_dt_relaxed[
                filter_model_row_select_j][[ontology_element_colname]]
            genes_in_ontology_i_j_intersection_num <- length(intersect(
                genes_in_ontology_i, genes_in_ontology_j))
            if (shared_elements_min_threshold < 
                genes_in_ontology_i_j_intersection_num) {
                ontologies_graph_edges[ontologies_graph_edges_counter,
                    c("from", "to", "weight") := list(
                        ontology_name_i, ontology_name_j,
                        genes_in_ontology_i_j_intersection_num)]
                ontologies_graph_edges_counter <- 
                    ontologies_graph_edges_counter + 1}}}
    ontologies_graph_edges <- ontologies_graph_edges[
        0:(ontologies_graph_edges_counter - 1), ]
    return(ontologies_graph_edges)
}

addExtraParamsToPlt <- function(graph_plot, p_value_max_threshold,
    p_value_type_colname) {graph_plot <- graph_plot + 
        ggraph::scale_edge_width(breaks = scales::pretty_breaks(n = 5), 
            range = c(0, 3), name="Nr. of shared elements") +
        ggraph::geom_node_point(aes(color = .data$p_stat)) +
        ggraph::geom_node_point(aes(color = .data$p_stat, 
            size = (1 - .data$p_stat)), show.legend = FALSE) +
        scale_size_area(max_size = 10) +
        scale_color_gradient2(mid =  '#ff6361', high = 'grey90',
            limits = c(0.0,p_value_max_threshold),
            name = p_value_type_colname) + 
        ggraph::geom_node_text(aes(label = .data$label), repel = TRUE) +
        ggraph::theme_graph()+ 
        theme(text=element_text(family="Helvetica"))
    return(graph_plot)
}

#' Plot Barplot
#' 
#' @description 
#' Plots barplot of p-values.
#'
#' @details 
#' Create a customized barplot of *p*-values, facilitating visual exploration 
#' and analysis of statistical significance within ontology categories.
#' 
#' @param reshaped_results  'data.table' in relaxed form, obtained as the 
#' output of the `reshape_results` function. 
#' The data source for generating the barplot.
#' @param selected_rows_to_plot A numeric vector specifying which rows of the 
#' reshaped results 'data.frame' should be included in the plot. 
#' Default is 'NULL'.
#' @param ontology_id_colname Character, specifies the column name that contains
#' ontology IDs in the input data.
#' @param p_value_type_colname Character, specifies the column name for 
#' *p*-values in the input data. 
#' Default is 'eFDR'.
#' @param p_value_max_threshold Numeric, representing the maximum *p*-value
#' threshold for filtering data. 
#' Default is 0.05.
#' @import ggplot2
#' @importFrom dplyr arrange
#' @seealso \code{\link{reshape_results}}
#' @export
#'
#' @return Returns a barplot.
#' 
#' @examples 
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'     max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(package = "mulea", "extdata", 
#'     "target_set.txt"))
#' background_genes <- readLines(system.file(
#'     package="mulea", "extdata", "background_set.txt"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered,
#'     # the test set variable
#'     element_names = sign_genes, 
#'     # the background set variable
#'     background_element_names = background_genes, 
#'     # the p-value adjustment method
#'     p_value_adjustment_method = "eFDR", 
#'     # the number of permutations
#'     number_of_permutations = 10000,
#'     # the number of processor threads to use
#'     nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)
#' 
#' # reshaping results for visualisation
#' ora_reshaped_results <- reshape_results(model = ora_model,
#'     model_results = ora_results,
#'     # choosing which column to use for the indication of significance
#'     p_value_type_colname = "eFDR")
#' 
#' # Plot barplot
#' plot_barplot(reshaped_results = ora_reshaped_results,
#'     # the column containing the names we wish to plot
#'     ontology_id_colname = "ontology_id",
#'     # upper threshold for the value indicating the significance
#'     p_value_max_threshold = 0.05,
#'     # column that indicates the significance values
#'     p_value_type_colname = "eFDR")

plot_barplot <- function(reshaped_results, ontology_id_colname = 'ontology_id',
    selected_rows_to_plot = NULL, p_value_type_colname = 'eFDR',
    p_value_max_threshold = 0.05) {
    validate_column_names_and_function_args(data = reshaped_results,
        p_value_type_colname, ontology_id_colname)
    reshaped_results <- filterRelaxedResultsForPlotting(
        reshaped_results = reshaped_results,
        p_value_type_colname = p_value_type_colname,
        p_value_max_threshold = p_value_max_threshold)
    if (is.null(selected_rows_to_plot)) {
        selected_rows_to_plot <- seq_len(nrow(reshaped_results))}
    unique_reshaped_results <- unique(reshaped_results[selected_rows_to_plot, 
        c(ontology_id_colname, p_value_type_colname), with = FALSE])
    unique_reshaped_results <- unique_reshaped_results %>%
        dplyr::arrange(dplyr::desc((!!as.name(p_value_type_colname))))
    unique_reshaped_results_df <- as.data.frame(unique_reshaped_results)
    unique_reshaped_results_df[, 1] <- factor(unique_reshaped_results_df[[1]],
        levels = unique_reshaped_results_df[[1]])
    mulea_gg_plot <- ggplot(unique_reshaped_results_df, 
        aes_string(x = ontology_id_colname, y = p_value_type_colname,
            fill = p_value_type_colname)) +
    geom_bar(stat = "identity") +
    scale_fill_gradient2(mid =  '#ff6361', high = 'grey90',
        limits = c(0.0,p_value_max_threshold),
        name = p_value_type_colname) +
    coord_flip() +
    theme_light() + 
    theme(text = element_text(family = "Helvetica"))
    mulea_gg_plot
}

#' Plot Lollipop
#' 
#' @description 
#' Plots lollipop plot of p-values.
#'
#' @details 
#' Create a customized  lollipop plot of p-values, facilitating visual 
#' exploration and analysis of statistical significance within ontology 
#' categories.
#' 
#' @param reshaped_results  data.table in relaxed form, obtained as the output 
#' of the `reshape_results` function. The data source for generating the 
#' barplot.
#' @param selected_rows_to_plot A numeric vector specifying which rows of the 
#' reshaped results data frame should be included in the plot. Default is NULL.
#' frame should be included in the plot?
#' @param ontology_id_colname Character, specifies the column name that contains
#' ontology IDs in the input data.
#' @param p_value_type_colname Character, specifies the column name for p-values
#' in the input data. Default is 'eFDR'.
#' @param p_value_max_threshold Numeric, representing the maximum p-value 
#' threshold for filtering data. Default is 0.05.
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom readr read_tsv
#' @import ggplot2
#' @seealso \code{\link{reshape_results}}
#' @export
#'
#' @return Returns a lollipop plot
#' 
#' @examples 
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'     max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(package = "mulea", "extdata", 
#'     "target_set.txt"))
#' background_genes <- readLines(system.file(
#'     package="mulea", "extdata", "background_set.txt"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered,
#'     # the test set variable
#'     element_names = sign_genes, 
#'     # the background set variable
#'     background_element_names = background_genes, 
#'     # the p-value adjustment method
#'     p_value_adjustment_method = "eFDR", 
#'     # the number of permutations
#'     number_of_permutations = 10000,
#'     # the number of processor threads to use
#'     nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)
#' 
#' # reshaping results for visualisation
#' ora_reshaped_results <- reshape_results(
#'     model = ora_model,
#'     model_results = ora_results,
#'     # choosing which column to use for the indication of significance
#'     p_value_type_colname = "eFDR")
#'         
#' # Plot lollipop
#' plot_lollipop(reshaped_results = ora_reshaped_results,
#'     # the column containing the names we wish to plot
#'     ontology_id_colname = "ontology_id",
#'     # upper threshold for the value indicating the significance
#'     p_value_max_threshold = 0.05,
#'     # column that indicates the significance values
#'     p_value_type_colname = "eFDR")
#' 
plot_lollipop <- function(reshaped_results, ontology_id_colname = 'ontology_id',
    selected_rows_to_plot = NULL, p_value_type_colname = 'eFDR',
    p_value_max_threshold = 0.05) {
    validate_column_names_and_function_args(data = reshaped_results,
        p_value_type_colname, ontology_id_colname)
    reshaped_results <- filterRelaxedResultsForPlotting(
        reshaped_results = reshaped_results,
        p_value_type_colname = p_value_type_colname,
        p_value_max_threshold = p_value_max_threshold)
    if (is.null(selected_rows_to_plot)) {
        selected_rows_to_plot <- seq_len(nrow(reshaped_results))}
    unique_reshaped_results <- unique(
        reshaped_results[selected_rows_to_plot, 
            c(ontology_id_colname, p_value_type_colname), with = FALSE])
    unique_reshaped_results <- unique_reshaped_results %>%
        dplyr::arrange(dplyr::desc((!!as.name(p_value_type_colname))))
    unique_reshaped_results_df <- as.data.frame(unique_reshaped_results)
    unique_reshaped_results_df[, 1] <- factor(
        unique_reshaped_results_df[[1]],
        levels = unique_reshaped_results_df[[1]])
    mulea_gg_plot <- ggplot(unique_reshaped_results_df,
        aes_string(x = ontology_id_colname, y = p_value_type_colname)) +
    geom_segment(aes(x = get(ontology_id_colname), 
        xend = get(ontology_id_colname), y = 0, 
        yend = as.numeric(get(p_value_type_colname))), color = 'black') +
    geom_point(aes(size = 5, color = get(p_value_type_colname))) +
    guides(size = 'none')+
    scale_color_gradient2(mid = '#ff6361', high = 'grey90',
        limits = c(0.0, p_value_max_threshold),
        name = p_value_type_colname) +
    coord_flip() +
    theme_light() + 
    theme(text = element_text(family="Helvetica"))
    mulea_gg_plot
}


#' Plot Heatmap
#' 
#' @description 
#' Plots heatmap of enriched terms and obtained p-values.
#' 
#' @details 
#'  The `plot_heatmap` function provides a convenient way to create a ggplot2 
#'  heatmap illustrating the significance of enriched terms within ontology 
#'  categories based on their associated p-values.
#' @param reshaped_results  data.table in relaxed form, obtained as the output 
#' of the `reshape_results` function. The data source for generating the 
#' barplot.
#' @param ontology_id_colname Character, specifies the column name that contains
#' ontology IDs in the input data.
#' @param p_value_type_colname Character, specifies the column name for p-values
#' in the input data. Default is 'eFDR'.
#' @param ontology_element_colname Character, specifying the column name that
#' contains ontology elements or terms in the input data. Default: 
#' 'element_id_in_ontology'.
#' @param p_value_max_threshold Numeric, representing the maximum p-value 
#' threshold for filtering data. Default is 0.05.
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom readr read_tsv
#' @import ggplot2
#' @seealso \code{\link{reshape_results}}
#' @export
#'
#' @return Returns a ggplot2 heatmap.
#' 
#' @examples 
#' library(mulea)
#' 
#' # loading and filtering the example ontology from a GMT file
#' tf_gmt <- read_gmt(file = system.file(package="mulea", "extdata", 
#'     "Transcription_factor_RegulonDB_Escherichia_coli_GeneSymbol.gmt"))
#' tf_gmt_filtered <- filter_ontology(gmt = tf_gmt, min_nr_of_elements = 3, 
#'     max_nr_of_elements = 400)
#' 
#' # loading the example data
#' sign_genes <- readLines(system.file(package = "mulea", "extdata", 
#'     "target_set.txt"))
#' background_genes <- readLines(system.file(
#'     package="mulea", "extdata", "background_set.txt"))
#'
#' # creating the ORA model
#' ora_model <- ora(gmt = tf_gmt_filtered,
#'     # the test set variable
#'     element_names = sign_genes, 
#'     # the background set variable
#'     background_element_names = background_genes, 
#'     # the p-value adjustment method
#'     p_value_adjustment_method = "eFDR", 
#'     # the number of permutations
#'     number_of_permutations = 10000,
#'     # the number of processor threads to use
#'     nthreads = 2)
#' # running the ORA
#' ora_results <- run_test(ora_model)
#' 
#' # reshaping results for visualisation
#' ora_reshaped_results <- reshape_results(
#'     model = ora_model,
#'     model_results = ora_results,
#'     # choosing which column to use for the indication of significance
#'     p_value_type_colname = "eFDR")
#' 
#' # Plot heatmap                                        
#' plot_heatmap(reshaped_results = ora_reshaped_results,
#'     # the column containing the names we wish to plot
#'     ontology_id_colname = "ontology_id",
#'     # column that indicates the significance values
#'     p_value_type_colname = "eFDR")

plot_heatmap <- function(reshaped_results,
    ontology_id_colname = 'ontology_id',
    ontology_element_colname = 'element_id_in_ontology',
    p_value_type_colname = 'eFDR',
    p_value_max_threshold = 0.05) {
    validate_column_names_and_function_args(data = reshaped_results,
        p_value_type_colname, ontology_element_colname, ontology_id_colname)
    model_with_res_dt_relaxed <- filterRelaxedResultsForPlotting(
        reshaped_results = reshaped_results,
        p_value_type_colname = p_value_type_colname,
        p_value_max_threshold = p_value_max_threshold)
    model_with_res_dt_relaxed_sort_pval <- model_with_res_dt_relaxed %>%
        dplyr::arrange(dplyr::desc((!!rlang::sym(p_value_type_colname))), 
            .by_group = FALSE)
    model_with_res_dt_relaxed_sort_pval[, 1] <- factor(
        model_with_res_dt_relaxed_sort_pval[[1]],
        levels = unique(model_with_res_dt_relaxed_sort_pval[[1]]))
    model_with_res_dt_relaxed_sort_pval[, 2] <- factor(
        model_with_res_dt_relaxed_sort_pval[[2]],
        levels = unique(rev(model_with_res_dt_relaxed_sort_pval[[2]])))
    ggplot(model_with_res_dt_relaxed_sort_pval, 
        aes(!!rlang::sym(ontology_element_colname),
            !!rlang::sym(ontology_id_colname),
        fill = !!rlang::sym(p_value_type_colname))) +
    geom_tile() +
    scale_fill_gradient2(mid = '#ff6361', high = 'grey90',
        limits = c(0.0,p_value_max_threshold)) +
    coord_fixed() +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90), 
        text = element_text(family = "Helvetica"))
}
