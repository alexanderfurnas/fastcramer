test_that("statistic matches cramer::cramer.test for every kernel", {
  skip_if_not_installed("cramer")
  set.seed(1)
  x <- matrix(rnorm(120 * 6), 120); y <- matrix(rnorm(80 * 6, mean = 0.3), 80)
  for (k in c("cramer", "bahr", "log", "fraca", "fracb")) {
    ref <- cramer::cramer.test(x, y, just.statistic = TRUE,
                               kernel = switch(k, cramer = "phiCramer", bahr = "phiBahr", log = "phiLog",
                                               fraca = "phiFracA", fracb = "phiFracB"))$statistic
    expect_equal(cramer_statistic(x, y, kernel = k), ref, tolerance = 1e-6, label = k)
  }
})

test_that("chunking does not change the result", {
  set.seed(2)
  x <- matrix(rnorm(300 * 4), 300); y <- matrix(rnorm(250 * 4), 250)
  expect_equal(cramer_statistic(x, y, chunk = 7L), cramer_statistic(x, y, chunk = 1000L))
})

test_that("vectors are treated as one-dimensional samples and inputs are validated", {
  v1 <- rnorm(50); v2 <- rnorm(40)
  expect_equal(cramer_statistic(v1, v2), cramer_statistic(matrix(v1), matrix(v2)))
  expect_error(cramer_statistic(matrix(1:6, 3), matrix(1:6, 2)), "same number of columns")
  expect_error(cramer_statistic(matrix(c(1, NA), 2), matrix(1:2, 2)), "missing")
})
