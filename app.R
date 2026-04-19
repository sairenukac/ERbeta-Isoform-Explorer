library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(DT)
library(ggpubr) # For p-values

# 1. DATA LOADING & PREP
# Load main data
data <- read_tsv("GSE104296_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
data_df <- as.data.frame(data)
rownames(data_df) <- data_df$GeneID
data_clean <- data_df[, -1]

# Load mapping file
gene_map <- read_csv("mapped_genes.csv")

# Metadata mapping
meta <- data.frame(
  sample = colnames(data_clean),
  condition = factor(c("Control", "Control", "ERb-KO", "ERb-KO", 
                       "ERb1", "ERb1", "ERb5", "ERb5"),
                     levels = c("Control", "ERb-KO", "ERb1", "ERb5"))
)

# 2. UI
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"),
  titlePanel("Glioblastoma: ERβ Isoform Explorer"),
  
  sidebarLayout(
    sidebarPanel(
      # Dropdown menu with search capability
      selectizeInput("gene_choice", "Search Gene Symbol or ID:", 
                     choices = setNames(gene_map$GeneID, 
                                        paste0(gene_map$Symbol, " (", gene_map$GeneID, ")")),
                     selected = 653635),
      
      checkboxInput("log", "Log2 Scale (TPM)", TRUE),
      checkboxInput("show_p", "Show P-values (vs Control)", TRUE),
      hr(),
      helpText("Select a gene to see how ERβ isoforms modulate its expression.")
    ),
    
    mainPanel(
      plotOutput("expressionPlot", height = "500px"),
      hr(),
      h4("Detailed Sample Data"),
      DTOutput("rawTable") # Interactive table
    )
  )
)

# 3. SERVER
server <- function(input, output) {
  
  # Reactive logic to fetch and join data
  plot_data <- reactive({
    req(input$gene_choice)
    
    # Extract row and format
    df <- data.frame(
      expression = as.numeric(data_clean[as.character(input$gene_choice), ]),
      sample = colnames(data_clean)
    ) %>%
      left_join(meta, by = "sample")
    
    if (input$log) {
      df$expression <- log2(df$expression + 1)
    }
    df
  })
  
  # Generate the Plot
  output$expressionPlot <- renderPlot({
    df <- plot_data()
    gene_info <- gene_map %>% filter(GeneID == input$gene_choice)
    
    p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "black") +
      geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.7) +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = paste("Expression of", gene_info$Symbol),
        subtitle = gene_info$`Gene Name / Type`,
        y = if(input$log) "log2(TPM + 1)" else "TPM",
        x = "Condition"
      ) +
      theme_minimal(base_size = 16) +
      theme(legend.position = "none")
    
    # Add P-values comparing everything to the 'Control' group
    if (input$show_p) {
      p <- p + stat_compare_means(method = "t.test", 
                                  ref.group = "Control", 
                                  label = "p.signif", 
                                  label.y.npc = 0.9)
    }
    
    p
  })
  
  # Generate Interactive Table
  output$rawTable <- renderDT({
    datatable(plot_data(), 
              options = list(pageLength = 8, dom = 't'), 
              rownames = FALSE,
              colnames = c("Expression Value", "GSM Sample ID", "Group")) %>%
      formatRound("expression", 3)
  })
}

shinyApp(ui, server)