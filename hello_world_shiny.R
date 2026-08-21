# My first Shiny app in R

library(shiny)

ui <- fluidPage(
	titlePanel("Chromosome 22 GWAS Explorer"),
	sidebarLayout(
		sidebarPanel(
			fileInput("gwas_file", "Upload GWAS results (CSV):", accept = ".csv"),
			helpText("Required columns: CHR, POS, SNP, P. Leave empty to use simulated data."),
			sliderInput(
				"position_range",
				"Chromosome 22 position:",
				min = 1,
				max = 51e6,
				value = c(1, 51e6),
				step = 1e5,
				sep = ""
			),
			numericInput(
				"top_n",
				"Variants to show:",
				value = 10,
				min = 1,
				max = 100,
				step = 1
			),
				actionButton("reset_view", "Reset position range")
		),
		mainPanel(
			plotOutput("manhattan_plot", height = "550px"),
			h4("Most significant variants"),
			tableOutput("top_variants"),
			verbatimTextOutput("data_status")
		)
	)
)

server <- function(input, output, session) {
	set.seed(22)
	demo_data <- data.frame(
		CHR = 22,
		POS = sort(sample(seq(1, 51e6), 5000)),
		SNP = sprintf("rs%07d", seq_len(5000)),
		P = runif(5000, min = 1e-8, max = 1)
	)

	demo_data$P[sample(seq_len(nrow(demo_data)), 12)] <- 10^runif(12, -12, -6)

	observeEvent(input$reset_view, {
		updateSliderInput(session, "position_range", value = c(1, 51e6))
	})

	gwas_data <- reactive({
		if (is.null(input$gwas_file)) {
			return(demo_data)
		}

		data <- read.csv(input$gwas_file$datapath, stringsAsFactors = FALSE)
		required_columns <- c("CHR", "POS", "SNP", "P")
		validate(need(all(required_columns %in% names(data)), "CSV must contain CHR, POS, SNP, and P columns."))
		data <- data[, required_columns]
		data$CHR <- as.numeric(data$CHR)
		data$POS <- as.numeric(data$POS)
		data$P <- as.numeric(data$P)
		validate(need(all(is.finite(data$POS)) && all(is.finite(data$P)), "POS and P must contain only numeric values."))
		validate(need(all(data$P > 0 & data$P <= 1), "P values must be greater than 0 and no greater than 1."))
		data[data$CHR == 22, ]
	})

	filtered_data <- reactive({
		data <- gwas_data()
		data[data$POS >= input$position_range[1] & data$POS <= input$position_range[2], ]
	})

	output$manhattan_plot <- renderPlot({
		data <- filtered_data()
		validate(need(nrow(data) > 0, "No chromosome 22 variants are in the selected range."))
		plot(
			data$POS,
			-log10(data$P),
			pch = 20,
			col = ifelse(data$P < 5e-8, "firebrick", "steelblue"),
			xlab = "Position on chromosome 22",
			ylab = expression(-log[10](P)),
			main = "Chromosome 22 Manhattan plot"
		)
		abline(h = -log10(5e-8), col = "firebrick", lty = 2)
		legend("topright", legend = "Genome-wide significance (P = 5e-8)", col = "firebrick", lty = 2, bty = "n")
	})

	output$top_variants <- renderTable({
		data <- filtered_data()
		validate(need(nrow(data) > 0, "No variants to display."))
		data <- data[order(data$P), , drop = FALSE]
		data <- head(data, input$top_n)
		data$P <- format(data$P, scientific = TRUE, digits = 3)
		data
	}, striped = TRUE, bordered = TRUE, spacing = "s")

	output$data_status <- renderPrint({
		data <- filtered_data()
		cat("Variants in view:", nrow(data), "\n")
		if (is.null(input$gwas_file)) cat("Source: simulated chromosome 22 data\n") else cat("Source: uploaded CSV\n")
	})
}

shinyApp(ui = ui, server = server)
