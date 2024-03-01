test_that("ora : object creation test.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "c", "d")
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2)
  
  testthat::expect_equal(mulea_ora_model@gmt, gmtMock)
  testthat::expect_equal(mulea_ora_model@element_names, c("a", "b", "c"))
  testthat::expect_equal(mulea_ora_model@background_element_names, 
                         c("a", "c", "d"))
  testthat::expect_equal(mulea_ora_model@p_value_adjustment_method, "eFDR")
})

test_that("ora : object creation test : adjustment type.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "c", "d")
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "BH", 
    nthreads = 2)
  
  testthat::expect_equal(mulea_ora_model@p_value_adjustment_method, "BH")
})

test_that("ora : element_names out of DB model.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "d")
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    p_value_adjustment_method = "eFDR",
    nthreads = 2
  )
  
  testthat::expect_warning(muleaTestRes <- run_test(mulea_ora_model))
                             
  testthat::expect_equal(muleaTestRes$p_value, 1)
})

test_that("ora : element_names out of pool.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  poolMock <- c("a", "b", "d")
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2)
  
  testthat::expect_warning(muleaTestRes <- run_test(mulea_ora_model))
  testthat::expect_equal(muleaTestRes$p_value, 1 / 3)
})

test_that("ora : matrix 2,2,2,2.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h")
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2)
  
  muleaTestRes <- run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$p_value, 53 / 70)
})

test_that("ora : pool >> var + DBi, matrix 2,2,2,18.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <-
    c(
      "a",
      "b",
      "c",
      "d",
      "e",
      "f",
      "g",
      "h",
      "i",
      "j",
      "k",
      "l",
      "m",
      "n",
      "o",
      "p",
      "q",
      "r",
      "s",
      "t",
      "u",
      "w",
      "x",
      "y"
    )
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2)
  
  muleaTestRes <- run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$p_value, 37 / 322)
})

test_that("ora : DBi not include pool, matrix 2,0,2,2.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "e", "f", "g", "h")
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2)
  
  muleaTestRes <- run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$p_value, 0.4)
})

test_that("ora : DB1 + DB2 => pool, matrix 1,3,2,2 and 2,2,1,3.", {
  gmtMock1 <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontology_id = "GO:0000002",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("e", "f", "g", "h"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2)
  testDataMock <- c("d", "e", "f")

  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    p_value_adjustment_method = "eFDR",
    nthreads = 2
  )
  
  muleaTestRes <- run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$p_value, c(13 / 14, 0.5))
})

test_that("ora : DB1 + DB2 => pool, matrix 2,2,2,0 and 2,2,1,3.", {
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
  poolMock <-
    c(
      "a",
      "b",
      "c",
      "d",
      "e",
      "f",
      "g",
      "h",
      "i",
      "j",
      "k",
      "l",
      "m",
      "n",
      "o",
      "p",
      "q",
      "r",
      "s",
      "t",
      "u",
      "w",
      "x",
      "y"
    )
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names = poolMock,
    p_value_adjustment_method = "eFDR", 
    nthreads = 2)
  
  muleaTestRes <- run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$p_value, c(37 / 322, 27 / 3542))
})


test_that("ora : DB1 + DB2 => pool, matrix 2,2,2,0 and 2,2,1,3.", {
  gmtMock <- read_gmt(file = system.file(
    package="mulea", "extdata", "model.gmt"))
  testDataMock <- c("FBgn0004407", "FBgn0010438", "FBgn0037044", 
                    "FBgn0002887", "FBgn0028434", "FBgn0030170", 
                    "FBgn0263831")
  poolMock <- unique(c(c("FBgn0033690", "FBgn0261618", "FBgn0004407", 
                         "FBgn0010438", "FBgn0032154", "FBgn0039930", 
                         "FBgn0040268", "FBgn0013674", "FBgn0037008", 
                         "FBgn0003116", "FBgn0037743", "FBgn0035401", 
                         "FBgn0037044", "FBgn0051005", "FBgn0026737", 
                         "FBgn0026751", "FBgn0038704", "FBgn0002887", 
                         "FBgn0028434", "FBgn0030170", "FBgn0263831", 
                         "FBgn0000579"),
                       c("FBgn0066666", "FBgn0000000", "FBgn0099999", 
                         "FBgn0011111", "FBgn0022222", "FBgn0777777", 
                         "FBgn0333333")))
  
  mulea_ora_model <- ora(
    gmt = gmtMock,
    element_names = testDataMock,
    background_element_names=poolMock,
    nthreads = 2
  )
  
  muleaTestRes <- run_test(mulea_ora_model)
  testthat::expect_equal(muleaTestRes$p_value[c(1,2)], c(1, 19205/26423))
})

test_that("ora : private : matrix 2,2,2,2.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "c", "d", "e", "f", "g", "h")
  
  mulea_hyper_test <- mulea:::MuleaHypergeometricTest(
    gmt = gmtMock,
    element_names = testDataMock,
    pool = poolMock, 
    nthreads = 2)
  
  muleaTestRes <- run_test(mulea_hyper_test)
  testthat::expect_equal(muleaTestRes$p.value, 53 / 70)
})


test_that("ora : private : DBi not include pool, matrix 2,0,2,2.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  poolMock <- c("a", "b", "e", "f", "g", "h")
  
  mulea_hyper_test <- MuleaHypergeometricTest(
    gmt = gmtMock,
    element_names = testDataMock,
    pool = poolMock, 
    nthreads = 2)
  
  muleaTestRes <- run_test(mulea_hyper_test)
  testthat::expect_equal(muleaTestRes$p.value, 0.4)
})


