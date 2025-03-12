library(shiny)
library(ggplot2)
library(dplyr)
library(ggiraph)
library(readr)

# Sample Data (Replace with actual data)
df <- req_fpt[1:10000,] %>% 
  filter(hits_ref < 100)

# UI
ui <- fluidPage(
  titlePanel("Filtered Polar Plot"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("tm_range", "Tm Range:",
                  min = min(df$tm), max = max(df$tm),
                  value = range(df$tm), step = 0.1),
      sliderInput("hits_ref_range", "Hits Ref Range:",
                  min = min(df$hits_ref), max = max(df$hits_ref),
                  value = range(df$hits_ref), step = 1),
      sliderInput("access_range", "Class Range:",
                  min = min(df$access), max = max(df$access),
                  value = range(df$access), step = 1),
      checkboxGroupInput("target_filter", "Select Targets:",
                         choices = unique(df$target),
                         selected = unique(df$target)),
      downloadButton("downloadData", "Download Selected Rows as TSV")
    ),
    mainPanel(
      girafeOutput("polar_plot"),
      girafeOutput("temp_plot")
    )
  )
)

# Server
server <- function(input, output, session) {
  filtered_data <- reactive({
    df %>%
      filter(
        tm >= input$tm_range[1], tm <= input$tm_range[2],
        hits_ref >= input$hits_ref_range[1], hits_ref <= input$hits_ref_range[2],
        access >= input$access_range[1], access <= input$access_range[2],
        target %in% input$target_filter
      )
  })
  
  # Reactive value to store selected rows (as a list of data_id)
  selected_rows <- reactiveVal(data.frame())
  
  # First plot (Polar plot)
  output$polar_plot <- renderGirafe({
    filtered_df <- filtered_data()
    
    # If filtered data is empty, return a blank plot
    if (nrow(filtered_df) == 0) {
      return(girafe(ggobj = ggplot() + theme_void()))
    }
    
    filtered_df$target_num <- as.numeric(factor(filtered_df$target))
    
    p <- ggplot(filtered_df, aes(x = factor(target), y = hits_ref, color = target)) +
      geom_jitter_interactive(aes(tooltip = as.character(probe_id), 
                                  data_id = as.character(probe_id))) +
      coord_polar() +
      theme_minimal() +
      scale_y_continuous(limits = c(-1, max(filtered_df$hits_ref))) +
      labs(title = "Filtered Polar Plot", x = "Target", y = "Hits Ref")
    
    girafe(ggobj = p)
  })
  
  # Second plot (Temperature plot)
  output$temp_plot <- renderGirafe({
    filtered_df <- filtered_data()
    
    # If filtered data is empty, return a blank plot
    if (nrow(filtered_df) == 0) {
      return(girafe(ggobj = ggplot() + theme_void()))
    }
    
    filtered_df$target_num <- as.numeric(factor(filtered_df$target))
    
    p <- ggplot(filtered_df, aes(x = factor(target), y = tm, color = target)) +
      geom_jitter_interactive(aes(tooltip = as.character(probe_id), 
                                  data_id = as.character(probe_id))) +
      theme_minimal() +
      labs(title = "Filtered Temperature Plot", x = "Target", y = "Tm")
    
    girafe(ggobj = p)
  })
  
  # Observe click events on both plots
  observe({
    click_data_polar <- input$polar_plot_click
    click_data_temp <- input$temp_plot_click
    
    selected_points <- NULL
    if (!is.null(click_data_polar)) {
      selected_points <- click_data_polar$id
    }
    if (!is.null(click_data_temp)) {
      selected_points <- click_data_temp$id
    }
    
    # Store the selected rows based on clicked data_id
    if (!is.null(selected_points)) {
      selected_rows(filtered_data() %>% filter(probe_id %in% selected_points))
    }
  })
  
  # Download handler to export selected rows as .tsv
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("selected_rows_", Sys.Date(), ".tsv", sep = "")
    },
    content = function(file) {
      # If rows are selected, write them to a .tsv file
      selected_data <- selected_rows()
      if (nrow(selected_data) > 0) {
        write_tsv(selected_data, file)
      }
    }
  )
}

# Run the app
shinyApp(ui, server)
