context("getHeight")

test_that("output is a vector", {
  expect_is(getHeight(
    tree = list(spp = 1, D1 = 30, H = 25),
    Dx = list(Dx = 7),
    bark = TRUE,
    mapping = NULL
  ), "numeric")
})
