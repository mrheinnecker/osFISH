library(shiny)
library(ggplot2)
library(dplyr)
library(ggiraph)

# Sample Data (Replace with actual data)
df <- req_fpt[1:10000,] %>% 
  filter(hits_ref<100)
  
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
      sliderInput("access_range", "class Range:",
                  min = min(df$access), max = max(df$access),
                  value = range(df$access), step = 1),
      checkboxGroupInput("target_filter", "Select Targets:",
                         choices = unique(df$target),
                         selected = unique(df$target))
    ),
    mainPanel(
      girafeOutput("polar_plot"),
      girafeOutput("temp_plot"),
    )
  )
)

# Server
server <- function(input, output) {
  filtered_data <- reactive({
    df %>%
      filter(
        tm >= input$tm_range[1], tm <= input$tm_range[2],
        hits_ref >= input$hits_ref_range[1], hits_ref <= input$hits_ref_range[2],
        access >= input$access_range[1], access <= input$access_range[2],
        target %in% input$target_filter
      )
  })
  
  output$polar_plot <- renderGirafe({
    filtered_df <- filtered_data()
    
    # If filtered data is empty, return a blank plot
    if (nrow(filtered_df) == 0) {
      return(girafe(ggobj = ggplot() + theme_void()))
    }
    
    filtered_df$target_num <- as.numeric(factor(filtered_df$target))
    
    p <- ggplot(filtered_df, aes(x = factor(target), y = hits_ref, color=target)) +
      geom_jitter_interactive(aes(tooltip = as.character(probe_id))) +
      coord_polar() +
      theme_minimal() +
      #scale_color_gradientn(colors = c("blue", "red"))+
      scale_y_continuous(limits = c(-1, max(filtered_df$hits_ref)))+
      labs(title = "Filtered Polar Plot", x = "Target", y = "Hits Ref")
    
    girafe(ggobj = p)
  })
  
  output$temp_plot <- renderGirafe({
    filtered_df <- filtered_data()
    
    # If filtered data is empty, return a blank plot
    if (nrow(filtered_df) == 0) {
      return(girafe(ggobj = ggplot() + theme_void()))
    }
    
    filtered_df$target_num <- as.numeric(factor(filtered_df$target))
    
    p <- ggplot(filtered_df, aes(x = factor(target), y = tm, color=target)) +
      geom_jitter_interactive(aes(tooltip = as.character(probe_id))) +
     # coord_polar() +
      theme_minimal() +
      #scale_color_gradientn(colors = c("blue", "red"))+
      #scale_y_continuous(limits = c(-1, max(filtered_df$hits_ref)))+
      labs(title = "Filtered Polar Plot", x = "Target", y = "Hits Ref")
    
    girafe(ggobj = p)
  })
  
  
}

# Run the app
shinyApp(ui, server)
