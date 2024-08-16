create_random_db <- function() {
  DB <- list()
  for (cat_i in seq_len(10)) {
    DB_cat_values <- c()
    for (el_i in seq_len(floor(runif(1, min = 2, max = 10)))) {
      DB_cat_values <- append(DB_cat_values, paste0("el_", el_i))
    }
    cat_id <- paste0("cat_", cat_i)
    DB[[cat_id]]  <- DB_cat_values
  }
  return(DB)
}

test_that("set.based.enrichment.test.", {
  steps <- 100
  DB <- create_random_db()
  pool <- unique(unlist(DB))
  select <- c("el_4", "el_5", "el_6")
  nthread <- 2
  res <- set.based.enrichment.test(
    steps = steps,
    pool = pool,
    select = select,
    DB = DB,
    nthread = nthread
  )
  
  testthat::expect_gt(res[["FDR"]][1], 0)
})


test_that("do_the_simulation", {
  steps <- 100
  DB <- create_random_db()
  pool <- unique(unlist(DB))
  select <- c("el_4", "el_5", "el_6")
  nthread <- 2
  list_of_all_genes <- unique(c(unlist(DB), pool))
  res_1 <- mulea:::do_the_simulation(list_of_all_genes = list_of_all_genes, 
                                   pool = pool, 
                                   select = select,
                                   DB = DB, 
                                   steps = steps, 
                                   nthread = nthread, 
                                   random_seed = 123)
  res_2 <- mulea:::do_the_simulation(list_of_all_genes = list_of_all_genes, 
                                   pool = pool, 
                                   select = select,
                                   DB = DB, 
                                   steps = steps, 
                                   nthread = nthread, 
                                   random_seed = 123)
  testthat::expect_equal(res_1, res_2)
})
