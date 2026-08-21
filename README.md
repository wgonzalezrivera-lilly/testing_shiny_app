# testing_shiny_app

A structured Shiny fine-mapping practice app for chromosome 22. It plots
posterior inclusion probabilities (PIP), rather than GWAS p-values, and can
load a real fine-mapping result file or deterministic demo data.

## Project structure

- `app.R`: UI and server composition.
- `R/data_utils.R`: input contract, validation, normalization, filtering, and ranking.
- `R/demo_data.R`: reproducible simulated chromosome 22 PIP data.
- `R/plot_utils.R`: reusable PIP and credible-set plot.
- `tests/testthat/`: unit tests for the pure data functions.

## Run the app

Install Shiny once if needed:

```r
install.packages("shiny")
```

Then run this from the project directory:

```r
shiny::runApp("app.R")
```

`run_finemap_app.R` remains as a compatibility launcher for the organized app.

## Input format

Uploaded CSV files must contain exactly these fields (additional fields are
ignored):

```text
CHR,POS,SNP,PIP
22,123456,rs123,0.83
```

`CHR` must be `22` or `chr22`, `POS` must be a positive base-pair position,
and `PIP` must be between 0 and 1. The app limits uploads to 100 MB and shows
validation errors in the app when a file is invalid.

## Tests and deployment

Install `testthat` once to run the unit tests:

```r
install.packages(c("shiny", "testthat"))
testthat::test_dir("tests/testthat")
```

For Posit Connect or shinyapps.io, deploy the project directory containing
`app.R` and the `R/` directory. The app has no database or external service
dependency, so it is also suitable for a Docker or on-premises Shiny Server
deployment.
