##########################################################
# Function for FDR corrected hypergeometric enrichment test
##########################################################
#' @import parallel
#' @import mulea

initialize_result_df <- function(pool, select, DB) {
    DB_names <- names(DB)
    number_of_DB_categories <- length(DB)  # Number of DB categories
    size_pool <- length(pool)
    size_select <- length(select)
    # init empty data frame with the row-size of DB 
    # (number of DB categories) to store the result
    result_df1 <- data.frame(DB_names = DB_names, DB_in_select = as.integer(NA),
        DB_in_pool = as.integer(NA), Genes_in_DB = as.integer(NA),
        P_val = as.double(NA), R_obs = as.integer(NA))
    # for every DB entity in the DB list
    for (i in seq_len(number_of_DB_categories)) {
        # create a vector of genes connected to the i-th DB category
        current_DB_category <- DB[[i]]
        # hypergometric test
        # q: number of common genes between a DBterm and select
        result_df1$DB_in_select[i] <- length(intersect(select,
            current_DB_category))
        # m: number of common genes between DBterm and BackGround
        result_df1$DB_in_pool[i] <- length(intersect(pool,
            current_DB_category))
        result_df1$Genes_in_DB[i] <- length(current_DB_category)
    }
    # compute the p values (raw p values without any corrections)
    # n: number of non-pool genes among DB
    # k: number of genes in select
    # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
    result_df1$P_val <- 1-stats::phyper(result_df1$DB_in_select-1, 
        result_df1$DB_in_pool,size_pool-result_df1$DB_in_pool,size_select) 
    # can change the digits, this is important for 
    # the precision of '0' is R
    P_val_round <- round(result_df1$P_val, digits = 15)
    result_df1$R_obs <- rank(P_val_round, na.last = TRUE, ties.method ="max")
    result_df1$P_adj_Bonf <- stats::p.adjust(result_df1$P_val, 
        method = "bonferroni")
    result_df1$P_adj_BH <- stats::p.adjust(result_df1$P_val, method = "BH")
    return(result_df1)
}

do_the_simulation <- function(list_of_all_genes, pool, select,DB, steps,
    nthread) {if(nthread==1) { # simple call op C++ part of the code
        simulation_result_tbl <- tryCatch(enrichment_test_simulation(DB, 
        list_of_all_genes, pool, length(select), steps, 0), error = warning)  
    } else if (nthread>1 && round(nthread) == nthread) {# parallel call of C++ 
        vv <- floor(steps / nthread)
        cc <- steps %% nthread
        steps_per_thread <- c(rep(vv, nthread-cc), rep(vv+1 , cc))
        seeds_per_thread <- sample(2^15, nthread)
        stopifnot(sum(steps_per_thread)==steps)
        rm(cc, vv)
        requireNamespace("parallel")
        if (interactive()) {cl_outfile <- paste(tempdir(), 'paralell.log', 
            sep = "\\")
            cl <- makeCluster(spec=nthread, type="PSOCK", outfile=cl_outfile)
        } else { cl <- makeCluster(spec = nthread, type = "PSOCK")}
        clusterExport(cl, c("DB","list_of_all_genes", "pool", "select",
            "seeds_per_thread", "steps_per_thread"), envir = environment())
        clusterEvalQ(cl, devtools::load_all())
        result_of_paralel <- clusterApplyLB(cl = cl, seq_len(nthread), 
            function(idx) {simulation_result_tbl <- tryCatch(
                enrichment_test_simulation(DB, list_of_all_genes, pool, 
                    length(select), steps_per_thread[[idx]], 
                    seeds_per_thread[[idx]]), error = warning)  # RCPP
                return(simulation_result_tbl)})
        stopCluster(cl)
        simulation_result_tbl <- do.call(rbind, result_of_paralel)
    } else {  stop("nthred parameter must be a positive integer number")  }
    names(simulation_result_tbl) <- c("DB_in_pool", "intersect.size",
        "multiplicity") # set the column names
    simulation_result_tbl$p <- 1-stats::phyper(
        simulation_result_tbl$intersect.size-1,
        simulation_result_tbl$DB_in_pool,
        length(pool)-simulation_result_tbl$DB_in_pool, length(select))
    # test consistency
    stopifnot(steps*length(DB) == sum(simulation_result_tbl$multiplicity))
    if(! all(is.finite(simulation_result_tbl$p)) ) {
        stopifnot(all(is.finite(simulation_result_tbl$p)), "ERROR_002")  }
    # group by the equal p-values and summarise the multiplicities 
    simulation_result_tbl <- simulation_result_tbl[,c("multiplicity", "p")]
    simulation_result_tbl$p <- round(simulation_result_tbl$p, digits = 15)
    simulation_result_tbl <- simulation_result_tbl[order(
        simulation_result_tbl$p),]
    list11 <- split(simulation_result_tbl, simulation_result_tbl$p)
    for(i in seq(list11)) { if(nrow(list11[[i]])>1) { list11[[i]] <- data.frame(
        multiplicity = sum(list11[[i]]$multiplicity), p = list11[[i]]$p[1])}}
    simulation_result_tbl <- do.call(rbind, list11)
    rownames(simulation_result_tbl) <- NULL
    return(simulation_result_tbl)  }


set.based.enrichment.test <- function(steps, pool, select, DB, nthread=1, 
    debug=FALSE) {
    ## convert the database and the select and pollt to sorted integer lists
    list_of_all_genes <- unique(c(unlist(DB),pool))
    select <- intersect(select, pool)
    number_of_DB_categories <- length(DB)  # Number of DB categories
    result_df1 <- initialize_result_df(pool,select, DB)
    P_val_round <- round(result_df1$P_val, digits = 15)
    ############# simualtion ############ 
    if(length(select)==0 ) {result_df1$R_exp <- NaN
        result_df1$FDR <- NaN
        return(result_df1)}
    names(DB) <- NULL
    simulation_result_tbl <- do_the_simulation(list_of_all_genes, pool, 
        select,DB, steps, nthread)
    ##############################################
    # some preparation to help the binary search
    if(simulation_result_tbl$p[1]!=0) {
        simulation_result_tbl <- rbind(data.frame(multiplicity = 0, p = 0), 
            simulation_result_tbl)}
    if(simulation_result_tbl$p[nrow(simulation_result_tbl)]!=1) {
        simulation_result_tbl <- rbind(simulation_result_tbl, 
            data.frame(multiplicity = 0, p = 1))}
    simulation_result_tbl$cum_sum_multiplicity <- cumsum(
        simulation_result_tbl$multiplicity)
    # init an empty vector. I  will store here the result
    R_exp <- integer(number_of_DB_categories)
    # this id the greatest/last in cumsum.multiplicity
    NN <- simulation_result_tbl$cum_sum_multiplicity[nrow(
        simulation_result_tbl)] 
    for (i in seq_len(number_of_DB_categories)) {if(P_val_round[i]>=1) {
        R_exp[i] <- NN  # it makes things much faster 
        } else {target <- P_val_round[i]
            # binary search
            a1 <- 1
            a2 <- nrow(simulation_result_tbl) #-cnt.of.ones
            while(a1+1 < a2) {a3 <- floor((a1+a2)/2)
                current <- simulation_result_tbl$p[a3]
                if( current <= target) {a1 <- a3} else { a2<-a3 }}
            # at the end of the binary search the 
            # simulation_result_tbl$cum_sum_multiplicity[[a1]] 
            # is the count of  p where p<=target  and a2=a1+1
            R_exp[i] <- simulation_result_tbl$cum_sum_multiplicity[[a1]]}}
    result_df1$R_exp <- R_exp/steps
    result_df1$FDR <- result_df1$R_exp/result_df1$R_obs
    return(result_df1)
}
