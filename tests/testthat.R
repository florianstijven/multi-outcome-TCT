library(testthat)

# Run all tests from the project root.
testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
