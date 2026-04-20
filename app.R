library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(DT)
library(ggpubr) # For p-values

# Loading data
expression_data <- read_csv("expression_data_clean.csv") %>% 
  column_to_rownames("GeneID")

gene_map <- read_csv("annotated_with_stats.csv")

# Metadata mapping
meta <- data.frame(
  sample = colnames(data_clean),
  condition = factor(c("Control", "Control", "ERb-KO", "ERb-KO", 
                       "ERb1", "ERb1", "ERb5", "ERb5"),
                     levels = c("Control", "ERb-KO", "ERb1", "ERb5"))
)

# 2. UI
ui <- navbarPage(
  title = "Glioblastoma: ERβ Isoform Explorer",
  theme = bslib::bs_theme(bootswatch = "flatly"),
  
  # Tab 1: The Gene Explorer
  tabPanel("Gene Explorer",
           sidebarLayout(
             sidebarPanel(
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
               DT::DTOutput("rawTable")
             )
           )
  ),
  
  # Tab 2: Volcano Plot
  tabPanel("Volcano Plot",
           sidebarLayout(
             sidebarPanel(
               numericInput("logFC_cut", "Log2FC Threshold", value = 1, min = 0, step = 0.5),
               numericInput("p_cut", "P-value Threshold", value = 0.05, min = 0, max = 1),
               helpText("Significant genes are highlighted based on the thresholds above.")
             ),
             
             mainPanel(
               plotOutput("volcanoPlot", height = "500px", click = "plot_click"),
               hr(),
               h4("Clicked Gene Details"),
               tableOutput("clickedPoints")
             )
           )
  )
)

# 3. SERVER
server <- function(input, output) {
  
  # Reactive logic to fetch and join data
  plot_data <- reactive({
    req(input$gene_choice)
    
    # Extracting row and format
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
  
  # Generating the Boxplot
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
    
    if (input$show_p) {
      p <- p + stat_compare_means(method = "t.test", 
                                  ref.group = "Control", 
                                  label = "p.signif", 
                                  label.y.npc = 0.9)
    }
    
    p
  })
  
  # Generating the Clicked Table
  output$clickedPoints <- renderTable({
    req(input$plot_click)
    
    df_click <- gene_map %>%
      mutate(log10_p = -log10(pvalue))
    
    res <- nearPoints(df_click, input$plot_click, 
                      xvar = "log2FC", 
                      yvar = "log10_p",
                      threshold = 10, 
                      maxpoints = 1)
    
    res %>% 
      dplyr::select(Symbol, GeneID, log2FC, pvalue)
  })
  
  # Generating Interactive Table
  output$rawTable <- renderDT({
    datatable(plot_data(), 
              options = list(pageLength = 8, dom = 't'), 
              rownames = FALSE,
              colnames = c("Expression Value", "GSM Sample ID", "Group")) %>%
      formatRound("expression", 3)
  })
  
  # Generating Comprehensive Volcano Plot
  output$volcanoPlot <- renderPlot({
    req(gene_map)
    
    df <- gene_map %>%
      filter(!is.na(log2FC), !is.na(pvalue)) %>%
      mutate(status = case_when(
        log2FC > input$logFC_cut & pvalue < input$p_cut ~ "Up",
        log2FC < -input$logFC_cut & pvalue < input$p_cut ~ "Down",
        TRUE ~ "Not Significant"
      ))
    
    ggplot(df, aes(x = log2FC, y = -log10(pvalue), color = status)) +
      geom_point(alpha = 0.4, size = 1.5) +
      scale_color_manual(values = c("Up" = "firebrick", "Down" = "dodgerblue", "Not Significant" = "grey80")) +
      # Setting range between -10 to 10
      coord_cartesian(xlim = c(-10, 10)) + 
      theme_minimal(base_size = 14) + # Sets font size for everything at once
      geom_vline(xintercept = c(-input$logFC_cut, input$logFC_cut), linetype = "dashed") +
      geom_hline(yintercept = -log10(input$p_cut), linetype = "dashed") +
      labs(
        title = "Differential Expression: ERβ-KO vs Control",
        subtitle = "Dashed lines indicate user-defined thresholds",
        x = "log2(Fold Change)", 
        y = "-log10(P-value)"
      )
  })
}

shinyApp(ui, server)