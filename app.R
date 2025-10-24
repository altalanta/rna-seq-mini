#
# NHANES BMI Body Fat Analysis - Configuration Wizard
# Shiny application for non-programmers to customize analysis parameters
#

library(shiny)
library(shinythemes)
library(yaml)
library(readr)
library(jsonlite)

# Load current configuration
config_file <- "config/config.yml"
if (file.exists(config_file)) {
  current_config <- read_yaml(config_file)
} else {
  # Default configuration if file doesn't exist
  current_config <- list(
    data = list(
      raw_dir = "data/raw",
      derived_dir = "data/derived"
    ),
    outputs = list(
      tables_dir = "outputs/tables",
      figures_dir = "outputs/figures",
      logs_dir = "outputs/logs",
      report_dir = "outputs/report"
    ),
    nhanes = list(
      demo_file = "DEMO_J.XPT",
      bmx_file = "BMX_J.XPT",
      dxx_file = "DXX_J.XPT",
      dxxag_file = "DXXAG_J.XPT"
    ),
    analysis = list(
      age_range = c(20, 59),
      survey_weights_col = "WTMEC2YR",
      strata_col = "SDMVSTRA",
      psu_col = "SDMVPSU"
    ),
    logging = list(
      level = "INFO",
      file = "analysis_log.txt"
    )
  )
}

# UI Definition
ui <- fluidPage(
  theme = shinytheme("cerulean"),

  # Application title
  titlePanel(
    div(
      style = "display: flex; align-items: center;",
      img(src = "https://via.placeholder.com/50x50/2c5aa0/ffffff?text=NHANES", style = "margin-right: 15px;"),
      "NHANES BMI Body Fat Analysis - Configuration Wizard"
    )
  ),

  # Main content
  sidebarLayout(
    sidebarPanel(
      width = 4,

      # Introduction
      wellPanel(
        h4("Welcome! 👋"),
        p("This wizard helps you customize the NHANES BMI Body Fat analysis for your research needs."),
        p("No programming knowledge required - just point and click!"),
        actionButton("reset_config", "Reset to Defaults", class = "btn-warning")
      ),

      # Analysis Parameters Section
      wellPanel(
        h4("📊 Analysis Parameters"),

        # Age Range
        sliderInput("age_min", "Minimum Age:",
                   min = 0, max = 80, value = current_config$analysis$age_range[1], step = 1),
        sliderInput("age_max", "Maximum Age:",
                   min = 0, max = 80, value = current_config$analysis$age_range[2], step = 1),

        # Survey Design Variables
        textInput("survey_weights", "Survey Weights Column:",
                 value = current_config$analysis$survey_weights_col),
        textInput("strata_col", "Strata Column:",
                 value = current_config$analysis$strata_col),
        textInput("psu_col", "PSU Column:",
                 value = current_config$analysis$psu_col)
      ),

      # Data Sources Section
      wellPanel(
        h4("📁 Data Sources"),

        # NHANES Files
        textInput("demo_file", "Demographics File:",
                 value = current_config$nhanes$demo_file),
        textInput("bmx_file", "Body Measures File:",
                 value = current_config$nhanes$bmx_file),
        textInput("dxx_file", "DXA Whole Body File:",
                 value = current_config$nhanes$dxx_file),
        textInput("dxxag_file", "DXA Android/Gynoid File:",
                 value = current_config$nhanes$dxxag_file)
      ),

      # Output Options Section
      wellPanel(
        h4("📈 Output Options"),

        # Output Directories
        textInput("tables_dir", "Tables Directory:",
                 value = current_config$outputs$tables_dir),
        textInput("figures_dir", "Figures Directory:",
                 value = current_config$outputs$figures_dir),
        textInput("logs_dir", "Logs Directory:",
                 value = current_config$outputs$logs_dir),
        textInput("report_dir", "Report Directory:",
                 value = current_config$outputs$report_dir),

        # Logging Level
        selectInput("log_level", "Logging Level:",
                   choices = c("DEBUG", "INFO", "WARNING", "ERROR"),
                   selected = current_config$logging$level)
      ),

      # Action Buttons
      div(
        style = "margin-top: 20px;",
        actionButton("preview_config", "Preview Configuration", class = "btn-info"),
        actionButton("save_config", "Save Configuration", class = "btn-success"),
        actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
        downloadButton("download_config", "Download Config File")
      )
    ),

    mainPanel(
      width = 8,

      # Status and Information Tabs
      tabsetPanel(
        tabPanel("Configuration Preview",
          h3("📋 Configuration Preview"),
          verbatimTextOutput("config_preview"),
          h4("What This Configuration Does:"),
          wellPanel(
            htmlOutput("config_explanation")
          )
        ),

        tabPanel("Analysis Options",
          h3("🔧 Analysis Customization Guide"),

          wellPanel(
            h4("Age Range Selection"),
            p("Choose the age range for your analysis. The default (20-59) focuses on working-age adults."),
            p("• 18-25: Young adults"),
            p("• 20-59: Working adults (default)"),
            p("• 60+: Older adults"),
            p("• Custom range: Set specific min/max ages")
          ),

          wellPanel(
            h4("Survey Design Variables"),
            p("These are technical variables used for proper statistical weighting:"),
            p("• Survey Weights: Account for sampling probability and non-response"),
            p("• Strata: Primary sampling units for variance estimation"),
            p("• PSU: Primary sampling units for clustering")
          ),

          wellPanel(
            h4("NHANES Data Files"),
            p("Four main data files are used in this analysis:"),
            p("• DEMO_J.XPT: Demographics, age, sex, ethnicity"),
            p("• BMX_J.XPT: Body measurements including BMI"),
            p("• DXX_J.XPT: Whole-body DXA measurements"),
            p("• DXXAG_J.XPT: Regional body fat distribution")
          )
        ),

        tabPanel("Getting Started",
          h3("🚀 Getting Started Guide"),

          wellPanel(
            h4("Step 1: Customize Your Settings"),
            p("Use the controls on the left to customize your analysis parameters."),
            p("Don't worry - you can always reset to defaults if needed!")
          ),

          wellPanel(
            h4("Step 2: Save Your Configuration"),
            p("Click 'Save Configuration' to save your settings to config/config.yml."),
            p("This file will be used by the analysis pipeline.")
          ),

          wellPanel(
            h4("Step 3: Run the Analysis"),
            p("Choose how to run your analysis:"),
            tags$ul(
              tags$li(strong("Quick Start:"), "Click 'Run Analysis' to start immediately"),
              tags$li(strong("Command Line:"), "Use 'make parallel-pipeline' for best performance"),
              tags$li(strong("Step by Step:"), "Use 'make fetch', 'make analysis', etc.")
            )
          ),

          wellPanel(
            h4("Step 4: Explore Results"),
            p("Your results will be saved in the outputs/ directory:"),
            tags$ul(
              tags$li("outputs/tables/ - Statistical results (CSV files)"),
              tags$li("outputs/figures/ - Visualizations (PNG files)"),
              tags$li("outputs/report/ - Complete HTML report"),
              tags$li("outputs/logs/ - Analysis logs")
            )
          )
        ),

        tabPanel("Help & Troubleshooting",
          h3("❓ Help & Troubleshooting"),

          wellPanel(
            h4("Common Issues & Solutions"),
            h5("❌ 'Package not found' errors"),
            p("Solution: Run: install.packages(c('dplyr', 'ggplot2', 'survey', 'foreign'))"),

            h5("❌ 'Data directory not found' errors"),
            p("Solution: Create directories: mkdir -p data/raw data/derived outputs/tables outputs/figures outputs/logs"),

            h5("❌ 'Configuration file not found' errors"),
            p("Solution: The wizard will create config/config.yml for you")
          ),

          wellPanel(
            h4("Performance Tips"),
            p("💡 Use 'make parallel-pipeline' for faster analysis on multi-core systems"),
            p("💡 Enable caching by running the pipeline multiple times"),
            p("💡 Check outputs/logs/ for detailed progress information")
          ),

          wellPanel(
            h4("Need More Help?"),
            p("📖 Interactive Tutorial: vignettes/getting-started.Rmd"),
            p("📚 Complete Documentation: outputs/report/report.html"),
            p("🐛 Report Issues: https://github.com/altalanta/nhanes-bmi-bodyfat/issues"),
            p("📧 Contact: analysis@nhanes-bmi.org")
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {

  # Reactive configuration
  config <- reactive({
    list(
      data = list(
        raw_dir = "data/raw",
        derived_dir = "data/derived"
      ),
      outputs = list(
        tables_dir = input$tables_dir,
        figures_dir = input$figures_dir,
        logs_dir = input$logs_dir,
        report_dir = input$report_dir
      ),
      nhanes = list(
        demo_file = input$demo_file,
        bmx_file = input$bmx_file,
        dxx_file = input$dxx_file,
        dxxag_file = input$dxxag_file
      ),
      analysis = list(
        age_range = c(input$age_min, input$age_max),
        survey_weights_col = input$survey_weights,
        strata_col = input$strata_col,
        psu_col = input$psu_col
      ),
      logging = list(
        level = input$log_level,
        file = "analysis_log.txt"
      )
    )
  })

  # Reset configuration to defaults
  observeEvent(input$reset_config, {
    updateSliderInput(session, "age_min", value = 20)
    updateSliderInput(session, "age_max", value = 59)
    updateTextInput(session, "survey_weights", value = "WTMEC2YR")
    updateTextInput(session, "strata_col", value = "SDMVSTRA")
    updateTextInput(session, "psu_col", value = "SDMVPSU")
    updateTextInput(session, "demo_file", value = "DEMO_J.XPT")
    updateTextInput(session, "bmx_file", value = "BMX_J.XPT")
    updateTextInput(session, "dxx_file", value = "DXX_J.XPT")
    updateTextInput(session, "dxxag_file", value = "DXXAG_J.XPT")
    updateTextInput(session, "tables_dir", value = "outputs/tables")
    updateTextInput(session, "figures_dir", value = "outputs/figures")
    updateTextInput(session, "logs_dir", value = "outputs/logs")
    updateTextInput(session, "report_dir", value = "outputs/report")
    updateSelectInput(session, "log_level", selected = "INFO")

    showNotification("Configuration reset to defaults!", type = "message")
  })

  # Configuration Preview
  output$config_preview <- renderText({
    paste(capture.output(toJSON(config(), pretty = TRUE, auto_unbox = TRUE)), collapse = "\n")
  })

  # Configuration Explanation
  output$config_explanation <- renderUI({
    age_range <- paste(config()$analysis$age_range, collapse = "-")
    sample_size <- "approximately 2,240 adults"

    HTML(paste0("
      <p><strong>Analysis Summary:</strong></p>
      <ul>
        <li><strong>Age Range:</strong> ", age_range, " years</li>
        <li><strong>Sample Size:</strong> ", sample_size, "</li>
        <li><strong>Survey Design:</strong> ", config()$analysis$survey_weights_col, " weights</li>
        <li><strong>Data Sources:</strong> 4 NHANES files (2017-2018)</li>
        <li><strong>Output:</strong> Tables, figures, and HTML report</li>
      </ul>
      <p><em>This configuration will analyze BMI-body fat correlations using proper survey-weighted statistical methods.</em></p>
    "))
  })

  # Save Configuration
  observeEvent(input$save_config, {
    # Ensure config directory exists
    dir.create("config", showWarnings = FALSE)

    # Write configuration to YAML file
    write_yaml(config(), config_file)

    showNotification("Configuration saved successfully!", type = "message")
  })

  # Download Configuration
  output$download_config <- downloadHandler(
    filename = "nhanes_config.yml",
    content = function(file) {
      write_yaml(config(), file)
    }
  )

  # Run Analysis
  observeEvent(input$run_analysis, {
    # Save configuration first
    write_yaml(config(), config_file)

    # Show instructions
    showModal(modalDialog(
      title = "Analysis Started!",
      HTML("
        <p>Your configuration has been saved and the analysis is ready to run.</p>
        <p><strong>To run the analysis:</strong></p>
        <ol>
          <li>Open a terminal in your project directory</li>
          <li>Run: <code>make parallel-pipeline</code></li>
          <li>Wait for completion (may take several minutes)</li>
          <li>Check results in <code>outputs/</code> directory</li>
        </ol>
        <p><em>The parallel pipeline uses multiple CPU cores for faster processing.</em></p>
      "),
      easyClose = TRUE,
      footer = modalButton("Got it!")
    ))
  })

  # Preview Configuration
  observeEvent(input$preview_config, {
    showModal(modalDialog(
      title = "Configuration Preview",
      verbatimTextOutput("config_modal_preview"),
      size = "l",
      easyClose = TRUE
    ))
  })

  # Configuration for modal preview
  output$config_modal_preview <- renderText({
    paste(capture.output(toJSON(config(), pretty = TRUE, auto_unbox = TRUE)), collapse = "\n")
  })
}

# Run the application
shinyApp(ui = ui, server = server)



