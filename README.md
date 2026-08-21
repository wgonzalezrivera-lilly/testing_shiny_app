# testing_shiny_app

A basic Shiny practice app for exploring chromosome 22 GWAS results.

## Run the app

Install Shiny once if needed:

```r
install.packages("shiny")
```

Then run this from the project directory:

```r
shiny::runApp("hello_world_shiny.R")
```

The app displays a Manhattan-style plot, highlights genome-wide significant
variants, filters by chromosome position, and shows the most significant
variants in a table. It uses simulated chromosome 22 data until a CSV is
uploaded. Uploaded files must contain `CHR`, `POS`, `SNP`, and `P` columns.
