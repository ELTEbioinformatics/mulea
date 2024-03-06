test_that("GSEA : object creation test.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "c")
  scoreDataMock <- c(0.1, 0.5, 1)
  
  mulea_ranked_based_test_model <- gsea(
    gmt = gmtMock,
    element_names = testDataMock,
    element_scores = scoreDataMock)

  testthat::expect_equal(mulea_ranked_based_test_model@gmt, gmtMock)
  testthat::expect_equal(
    mulea_ranked_based_test_model@element_names, c("a", "b", "c"))
  testthat::expect_equal(
    mulea_ranked_based_test_model@element_scores, c(0.1, 0.5, 1))
})

test_that("GSEA : no element_scores vector.", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "d")
  
  mulea_ranked_based_test_model <- gsea(
    gmt = gmtMock, element_names = testDataMock)
  testthat::expect_error(muleaTestRes <-
                           run_test(mulea_ranked_based_test_model))
})
