# Compatibility launcher for the organized app.
source("R/data_utils.R", local = TRUE)
source("R/demo_data.R", local = TRUE)
source("R/plot_utils.R", local = TRUE)
source("app.R", local = TRUE)
shinyApp(ui = ui, server = server)
