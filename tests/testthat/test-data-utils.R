testthat::test_that("demo data follows the fine-mapping contract", {
  demo <- make_demo_finemap_data(n = 100, seed = 1)
  testthat::expect_silent(validate_finemap_data(demo))
  testthat::expect_equal(names(normalize_finemap_data(demo)), required_columns)
  testthat::expect_true(all(demo$PIP >= 0 & demo$PIP <= 1))
})

testthat::test_that("invalid fine-mapping data fails with useful errors", {
  invalid <- data.frame(CHR = 22, POS = 10, SNP = "rs1", PIP = 1.2)
  testthat::expect_error(validate_finemap_data(invalid), "between 0 and 1")

  wrong_chromosome <- data.frame(CHR = 21, POS = 10, SNP = "rs1", PIP = 0.2)
  testthat::expect_error(
    validate_finemap_data(wrong_chromosome),
    "chromosome 22"
  )
})

testthat::test_that("filtering and ranking are deterministic", {
  data <- data.frame(
    CHR = 22,
    POS = c(10, 20, 30),
    SNP = c("rs1", "rs2", "rs3"),
    PIP = c(0.2, 0.9, 0.5)
  )
  filtered <- filter_finemap_data(data, c(15, 30))
  testthat::expect_equal(filtered$SNP, c("rs2", "rs3"))
  testthat::expect_equal(top_finemap_variants(filtered, 1)$SNP, "rs2")
})
