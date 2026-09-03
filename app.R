library(shiny)

MAX_POSITION <- 51e6
DEFAULT_POSITION_RANGE <- c(1, MAX_POSITION)

ui <- fluidPage(
	tags$head(
		tags$style(HTML("
			.sidebar-toggle-btn { margin-right: 8px; text-decoration: none; }
			#sidebar_wrapper > .row > div:first-child { transition: width 0.15s ease; }
			#sidebar_wrapper.collapsed > .row > div:first-child { display: none; }
			#sidebar_wrapper.collapsed > .row > div:nth-child(2) { width: 100%; flex: 0 0 100%; max-width: 100%; }
		")),
		tags$script(HTML("
			$(document).on('click', '#toggle_sidebar', function(e) {
				e.preventDefault();
				$('#sidebar_wrapper').toggleClass('collapsed');
			});
		"))
	),
	titlePanel(
		tagList(
			actionLink("toggle_sidebar", label = icon("bars"), class = "sidebar-toggle-btn"),
			"Chromosome 22 Fine-mapping Explorer"
		)
	),
	div(
		id = "sidebar_wrapper",
		class = "collapsed",
		navlistPanel(
			id = "main_nav",
		well = TRUE,
		widths = c(2, 10),
		tabPanel(
			"Explorer",
			fluidRow(
				column(
					width = 3,
					wellPanel(
						fileInput(
							"finemap_file",
							"Upload fine-mapping results (CSV)",
							accept = c(".csv", "text/csv")
						),
						helpText("Required columns: CHR, POS, SNP, PIP."),
						sliderInput(
							"position_range",
							"Position range (bp)",
							min = DEFAULT_POSITION_RANGE[1],
							max = DEFAULT_POSITION_RANGE[2],
							value = DEFAULT_POSITION_RANGE,
							step = 1e5,
							sep = ""
						),
						sliderInput(
							"pip_threshold",
							"PIP highlight threshold",
							min = 0,
							max = 1,
							value = 0.1,
							step = 0.01
						),
						sliderInput(
							"credible_set_threshold",
							"Credible set coverage",
							min = 0.5,
							max = 0.99,
							value = 0.95,
							step = 0.01
						),
						numericInput("top_n", "Rows in variant table", value = 10, min = 1, max = 100, step = 1),
						actionButton("reset_view", "Reset view", class = "btn-primary"),
						downloadButton("download_variants", "Download visible variants")
					)
				),
				column(
					width = 9,
					plotOutput("pip_plot", height = "560px"),
					tags$hr(),
					h4("Highest-PIP variants"),
					tableOutput("variant_table"),
					verbatimTextOutput("data_status")
				)
			)
		),
		tabPanel(
			"About",
			h4("About this app"),
			p("This app explores statistical fine-mapping results for chromosome 22, highlighting variants with high posterior inclusion probability (PIP) and their credible sets.")
		),
		tabPanel(
			"Settings",
			h4("Settings"),
			p("App-wide settings will live here.")
		),
		"More to come \U0001F642"
	)
	)
)

server <- function(input, output, session) {
	demo_data <- make_demo_finemap_data()

	data_state <- reactive({
		if (is.null(input$finemap_file)) {
			return(list(data = demo_data, error = NULL, source = "simulated demo data"))
		}
		if (input$finemap_file$size > MAX_UPLOAD_BYTES) {
			return(list(data = NULL, error = "The uploaded file is larger than 100 MB.", source = "uploaded CSV"))
		}
		tryCatch(
			list(data = read_finemap_csv(input$finemap_file$datapath), error = NULL, source = input$finemap_file$name),
			error = function(error) list(data = NULL, error = conditionMessage(error), source = input$finemap_file$name)
		)
	})

	current_data <- reactive({
		state <- data_state()
		validate(need(is.null(state$error), state$error))
		state$data
	})

	observeEvent(data_state(), {
		state <- data_state()
		if (!is.null(state$error)) return()
		data <- state$data
		data_max <- max(data$POS)
		updateSliderInput(
			session,
			"position_range",
			min = 1,
			max = max(data_max, 1),
			value = c(1, data_max),
			step = max(1, round(data_max / 500))
		)
	}, ignoreInit = TRUE)

	observeEvent(input$reset_view, {
		data <- current_data()
		data_max <- max(data$POS)
		updateSliderInput(session, "position_range", min = 1, max = data_max, value = c(1, data_max))
	})

	visible_data <- reactive({
		filter_finemap_data(current_data(), input$position_range)
	})

	output$pip_plot <- renderPlot({
		plot_finemap_pip(current_data(), input$position_range, input$pip_threshold, input$credible_set_threshold)
	})

	output$variant_table <- renderTable({
		data <- visible_data()
		validate(need(nrow(data) > 0, "No variants are in the selected position range."))
		result <- top_finemap_variants(data, input$top_n)
		result$PIP <- formatC(result$PIP, format = "f", digits = 4)
		result
	}, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

	output$data_status <- renderPrint({
		state <- data_state()
		if (!is.null(state$error)) {
			cat("Data error:", state$error, "\n")
			return(invisible(NULL))
		}
		data <- visible_data()
		cat("Source:", state$source, "\n")
		cat("Variants in view:", nrow(data), "of", nrow(state$data), "\n")
		cat("Total PIP in view:", formatC(sum(data$PIP), format = "f", digits = 4), "\n")
	})

	output$download_variants <- downloadHandler(
		filename = function() paste0("chr22_visible_finemapping_", Sys.Date(), ".csv"),
		content = function(file) write.csv(visible_data(), file, row.names = FALSE)
	)
}

shinyApp(ui = ui, server = server)
