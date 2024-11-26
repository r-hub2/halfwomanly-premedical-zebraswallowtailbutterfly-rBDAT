testthat::context("Testing rBDAT:::throw()")
testthat::test_that("throw the rBDAT exception", {
  error_message <- "hello, testthat"
  string <- "hello, testthat"
  testthat::expect_error(
    rBDAT:::throw(string),
    error_message
  )
})
