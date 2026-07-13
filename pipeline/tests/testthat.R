#!/usr/bin/env Rscript

if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_dir("tests/testthat", reporter = "summary")
} else {
  stop("Install testthat: install.packages('testthat')")
}
