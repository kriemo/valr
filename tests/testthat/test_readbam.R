context("read_bam")

test_that("reading bam works", {
  res <- bam_to_df()
  expect_equal(nrow(res), 9976)
})
