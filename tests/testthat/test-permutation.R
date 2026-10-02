test_that("permutation test detects a shift and not a null", {
  set.seed(3)
  x <- matrix(rnorm(150 * 5), 150); y <- matrix(rnorm(150 * 5, mean = 0.6), 150)
  shifted <- cramer_test(x, y, replicates = 199, seed = 1)
  expect_lt(shifted$p.value, 0.05)
  same <- cramer_test(x, matrix(rnorm(150 * 5), 150), replicates = 199, seed = 1)
  expect_gt(same$p.value, 0.05)
  expect_s3_class(shifted, "cramer_test")
  expect_length(shifted$null, 199)
})

test_that("blocked shuffling preserves each stratum's label composition", {
  n <- 400; block <- sample(letters[1:5], n, replace = TRUE)
  g <- c(rep(TRUE, 150), rep(FALSE, 250))
  shuffle <- fastcramer:::make_shuffler(block, n)
  for (i in 1:5) {
    s <- shuffle(g)
    expect_equal(tapply(s, block, sum), tapply(g, block, sum))
    expect_equal(sum(s), sum(g))
  }
})

test_that("bootstrap design returns a sensible summary", {
  set.seed(4)
  x <- matrix(rnorm(500 * 6), 500); y <- matrix(rnorm(300 * 6, mean = 0.4), 300)
  b <- cramer_bootstrap(x, y, replicates = 30, draw = 150, seed = 1)
  expect_s3_class(b, "cramer_bootstrap")
  expect_equal(nrow(b$replicates), 60)
  expect_gt(b$summary$exceedance, 0.9)
  expect_equal(b$draw, 150)
  expect_equal(cramer_bootstrap(x, y, replicates = 2)$draw, 300)   # min(5000, m, n)
})
