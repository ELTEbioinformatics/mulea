test_that("Methods : checkIfPoolIncludeSample false", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "e", "f")
  
  testthat::expect_warning(checkIfPoolIncludeSample(
    model = gmtMock, sampleVector = testDataMock))
})

test_that("Methods : checkIfPoolIncludeSample false", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  testDataMock <- c("a", "b", "d")
  
  testthat::expect_equal(
    checkIfPoolIncludeSample(model = gmtMock, sampleVector = testDataMock),
    c("a", "b", "d")
  )
})


test_that("Methods : cutGmtToPool", {
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  poolMock <- c("a", "b", "e", "f", "g", "h")
  
  testthat::expect_equal(cutGmtToPool(
    gmt = gmtMock, pool = poolMock)[['list_of_values']][[1]],
                         list(c("a", "b"))[[1]])
})


test_that("Utils : read_gmt", {
  filePath <- tempfile(fileext=".gmt")
  fileConn <- file(filePath)
  writeLines(
    c('ID:0000001	"mitochondrion inheritance"	FBgn0033690	FBgn0261618'), 
    fileConn)
  close(fileConn)
  modelDfFromFile <- read_gmt(file = filePath)
  
  testthat::expect_equal(modelDfFromFile$ontology_id, "ID:0000001")
  testthat::expect_equal(modelDfFromFile$ontology_name, 
                         "\"mitochondrion inheritance\"")
  testthat::expect_equal(modelDfFromFile$list_of_values[[1]], 
                         c("FBgn0033690", "FBgn0261618"))
})


test_that("Utils : write_gmt", {
  filePath <- tempfile(fileext=".gmt")
  gmtMock <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d"))),
    stringsAsFactors = FALSE
  )
  
  write_gmt(gmt = gmtMock, file = filePath)
  modelDfFromFile <- read_gmt(file = filePath)
  
  testthat::expect_equal(modelDfFromFile$ontology_id, "GO:0000001")
  testthat::expect_equal(modelDfFromFile$ontology_name, 
                         "Imagin gen ontology to tests.")
  testthat::expect_equal(modelDfFromFile$list_of_values[[1]], 
                         c("a", "b", "c", "d"))
})


test_that("Utils : filter_ontology default", {
  gmtMock1 <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontology_id = "GO:0000002",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("e", "f", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock3 <- data.frame(
    ontology_id = "GO:0000003",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d", "e", "f", "g", "h", "i"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2, gmtMock3)
  
  gmt_filtered_model <- filter_ontology(
    gmt = gmtMock,
    min_nr_of_elements = 5,
    max_nr_of_elements = 350
  )
  
  testthat::expect_equal(nrow(gmt_filtered_model), 1)
  testthat::expect_equal(gmt_filtered_model$ontology_id, "GO:0000003")
  testthat::expect_equal(gmt_filtered_model$list_of_values[[1]], 
                         c("a", "b", "c", "d", "e", "f", "g", "h", "i"))
})

test_that("Utils : filter_ontology between min and max", {
  gmtMock1 <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontology_id = "GO:0000002",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("e", "f", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock3 <- data.frame(
    ontology_id = "GO:0000003",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d", "e", "f", "g", "h", "i"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2, gmtMock3)
  
  gmt_filtered_model <- filter_ontology(gmt = gmtMock, 
                                        min_nr_of_elements = 3, 
                                        max_nr_of_elements = 5)
  
  testthat::expect_equal(nrow(gmt_filtered_model), 1)
  testthat::expect_equal(gmt_filtered_model$ontology_id, "GO:0000002")
  testthat::expect_equal(gmt_filtered_model$list_of_values[[1]], 
                         c("e", "f", "c", "d"))
})

test_that("Utils : filter_ontology above min", {
  gmtMock1 <- data.frame(
    ontology_id = "GO:0000001",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a"))),
    stringsAsFactors = FALSE
  )
  gmtMock2 <- data.frame(
    ontology_id = "GO:0000002",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("e", "f", "c", "d"))),
    stringsAsFactors = FALSE
  )
  gmtMock3 <- data.frame(
    ontology_id = "GO:0000003",
    ontology_name = "Imagin gen ontology to tests.",
    list_of_values = I(list(c("a", "b", "c", "d", "e", "f", "g", "h", "i"))),
    stringsAsFactors = FALSE
  )
  gmtMock <- rbind(gmtMock1, gmtMock2, gmtMock3)
  
  gmt_filtered_model <- filter_ontology(
    gmt = gmtMock,
    min_nr_of_elements = 3,
    max_nr_of_elements = 350
  )
  
  testthat::expect_equal(nrow(gmt_filtered_model), 2)
  testthat::expect_equal(
    gmt_filtered_model$ontology_id, c("GO:0000002", "GO:0000003"))
})


test_that("Utils : helper_set_up_namespace_params dll loader helper", {
  sink(tempfile())
  callRes <- helper_set_up_namespace_params()
  sink()
  testthat::expect_equal(callRes, "Not implemented helper func.")
})
