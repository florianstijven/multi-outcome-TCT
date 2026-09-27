testthat::test_that("build_summing_contrast_multi_outcome() works with strata", {
  contrast_matrix <- build_summing_contrast_multi_outcome(times = list(1:3, 1:3, 1:5, 1:5),
                                                          strata = c(1, 1, 2, 2))
  testthat::expect_equal(dim(contrast_matrix), c(6, 32))
  testthat::expect_equal(sum(contrast_matrix), 0)
  testthat::expect_equal(contrast_matrix[1, 1:8], c(0, -1, 0, 0, 1, 0, 0, -1))
})
