# 22JAN24 - Creating an advanced interactive plot 

# Check and set the correct working directory
getwd()
# setwd("/Users/evaedwards/Final-Year-Project/Datasets")

# Read in data that previously made
data <- read.csv("interactiveplot_updated.csv")
head(data)

# Define UI
ui <- fluidPage(
  titlePanel("Interactive Mutation Frequency Plot"),
  sidebarLayout(
    sidebarPanel(
      # selectInput("mutation", "Select Mutation ID:", 
      #             choices = unique(data$mutation_ID), 
      #             selected = unique(data$mutation_ID)[1],
      #             multiple = TRUE),
      selectInput("type", "Filter by Type:", 
                  choices = c("All", unique(data$type)), 
                  selected = "All",
                  multiple = TRUE),
      selectInput("strain", "Filter by Strain:",
                  choices = c("All", unique(data$Strain)),
                  selected = "All",
                  multiple = TRUE),
      selectInput("Control", "Filter by Presence of control:",
                  choices = c("All", unique(data$Control)),
                  selected = "All",
                  multiple = TRUE),
      selectInput("Count", "Filter by amount of schools this gene is mutated in:",
                  choices = c("All", unique(data$Count)),
                  selected = "All",
                  multiple = TRUE),
      selectInput("mut_type", "Filter by Mutation Type:", 
                  choices = c("All", unique(data$mut_type)), 
                  selected = "All",
                  multiple = TRUE),
      selectInput("gene", "Filter by Gene:", 
                  choices = c("All", unique(data$gene)), 
                  selected = "All",
                  multiple = TRUE),
      selectInput("school", "Filter by School:", 
                  choices = c("All", unique(data$School)), 
                  selected = "All",
                  multiple = TRUE),
      selectInput("background", "Filter by Background:", 
                  choices = c("All", unique(data$Background)), 
                  selected = "All",
                  multiple = TRUE),
    ),
    mainPanel(
      plotOutput("freqPlot")
    )
  )
)

server <- function(input, output) {
  
  filtered_data <- reactive({
    filtered <- data
    
    # Apply all filters using AND logic
    if (!"All" %in% input$gene) {
      filtered <- filtered[filtered$gene %in% input$gene, ]
    }
    if (!"All" %in% input$type) {
      filtered <- filtered[filtered$type %in% input$type, ]
    }
    if (!"All" %in% input$Control) {
      filtered <- filtered[filtered$Control %in% input$Control, ]
    }
    if (!"All" %in% input$mut_type) {
      filtered <- filtered[filtered$mut_type %in% input$mut_type, ]
    }
    if (!"All" %in% input$strain) {
      filtered <- filtered[filtered$Strain %in% input$strain, ]
    }
    if (!"All" %in% input$school) {
      filtered <- filtered[filtered$School %in% input$school, ]
    }
    if (!"All" %in% input$background) {
      filtered <- filtered[filtered$Background %in% input$background, ]
    }
    if (!"All" %in% input$Count) {
      filtered <- filtered[filtered$Count %in% input$Count, ]
    }
    if (length(input$mutation) > 0) {
      filtered <- filtered[filtered$mutation_ID %in% input$mutation, ]
    }
    
    filtered  # Return the filtered data
  })
  
  
  # Render Plot
  output$freqPlot <- renderPlot({
    plot_data <- filtered_data()
    
    # Generate unique colors for each gene
    gene_colors <- scales::hue_pal()(length(unique(plot_data$gene)))
    names(gene_colors) <- unique(plot_data$gene)  # Assign colors to genes
    
    # Add a color column to the plot data for gene coloring
    plot_data$color <- gene_colors[plot_data$gene]
    
    # Plot with gene colors
    par(mar = c(20, 6, 4, 2)) # Making the margins bigger
    ylim_max <- max(plot_data$frequency) * 1.2  # Add 20% space above bars
    
    barplot(
      plot_data$frequency, 
      names.arg = plot_data$mutation_ID, 
      las = 2, # Rotate x-axis labels
      cex.names = 0.5,
      col = plot_data$color,  # Color bars by gene
      main = "Frequency of Mutations",
      xlab = "Mutation ID",
      ylab = "Frequency",
      ylim = c(0, ylim_max)
    )
  }, height = 600, width = 900)
}

# Run the application
shinyApp(ui = ui, server = server)



