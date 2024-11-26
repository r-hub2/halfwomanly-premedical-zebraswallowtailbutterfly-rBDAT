if (interactive()) pkgload::load_all()
context("getVolume")

# Fails on CRAN running stale versions under windows and r-devel on new M1-mac
# hence allow some tolerance on comparison
test_equal <- function(...) testthat::expect_equal(..., tolerance =  1e-06)

test_that("test getVolume incl. bark", {
  result <- getVolume(
    tree = list(spp = 1, D1 = 30, H = 25),
    AB = list(A = 1, B = 2, sl = 0.1),
    iAB = "H",
    bark = TRUE,
    mapping = NULL
  )
  if (FALSE) {
    dump("result")
  }
  expectation <- 0.0686297789216042
  test_equal(result, expectation)
})

test_that("test getVolume excl. bark", {
  result <- getVolume(
    tree = list(spp = 1, D1 = 30, H = 25),
    AB = list(A = 1, B = 2, sl = 0.1),
    iAB = "H",
    bark = FALSE,
    mapping = NULL
  )
  if (FALSE) {
    dump("result")
  }
  expectation <- 0.0615232065320015
  test_equal(result, expectation)
})

test_that("test getAssortment", {
  result <- getAssortment(tree = list(spp = 1, D1 = 30, H = 25))
  if (FALSE) {
    dump("result")
  }
  expectation <- structure(list(
    tree = c(1L, 1L, 1L, 1L, 1L),
    No = 1:5,
    Sort = c("X", "Sth", "Ab", "Ind", "nvDh"),
    height = c(0, 0.25, 0, 0, 20.0459995269775),
    length = c(0, 19.6000003814697, 0, 0, 2.20400047302246),
    midD = c(0, 20.8709373474121, 0, 0, 8.20349884033203),
    topD = c(0, 10.5018358230591, 0, 0, 6.01526641845703),
    Vol = c(0, 0.670548021793365, 0, 0, 0.0116492994129658)
  ),
  row.names = c(NA, 5L), class = "data.frame"
  )
  test_equal(result, expectation)
})
