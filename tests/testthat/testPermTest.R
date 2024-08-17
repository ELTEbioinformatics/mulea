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

test_that("do_the_simulation_same_seed", {
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

test_that("do_the_simulation_diff_seed", {
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
                                     random_seed = 321)
  testthat::expect_failure(testthat::expect_equal(res_1, res_2))
})

test_that("ora_same_seed", {
  gmtMock1 <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontology_id = "GO:0000002",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("e", "f", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("b", "d", "e", "f")
  poolMock <- c("a","b","c","d","e","f","g","h","i","j",
                "k","l","m","n","o","p","q","r","s","t","u","w","x","y"
  )
  
  mulea_ora_model_1 <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2,
    random_seed = 123)
  
  mulea_test_res_1 <- run_test(mulea_ora_model_1)
  
  mulea_ora_model_2 <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2,
    random_seed = 123)
  
  mulea_test_res_2 <- run_test(mulea_ora_model_2)
  
  testthat::expect_equal(mulea_test_res_1, mulea_test_res_2)
})

test_that("ora_diff_seed", {
  gmtMock1 <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontology_id = "GO:0000002",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("e", "f", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("b", "d", "e", "f")
  poolMock <- c("a","b","c","d","e","f","g","h","i","j",
      "k","l","m","n","o","p","q","r","s","t","u","w","x","y"
  )
  
  mulea_ora_model_1 <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2,
    random_seed = 123)
  
  mulea_test_res_1 <- run_test(mulea_ora_model_1)
  
  mulea_ora_model_2 <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2,
    random_seed = 321)
  
  mulea_test_res_2 <- run_test(mulea_ora_model_2)
  
  testthat::expect_failure(
    testthat::expect_equal(mulea_test_res_1, mulea_test_res_2))
})
