library(shiny)
library(plotly)
library(ggplot2)
library(ggpubr)
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
library(shinyWidgets)
library(ComplexHeatmap)

options(scipen = 999)

#UI Function
ui = navbarPage(
  "Expression plots",
  tabPanel("Volcano plot",
           sidebarLayout(
             sidebarPanel(
               div(
                 style = "font-size: 12px;",
                 helpText(
                   "Text must be in CSV format with 3 columns. First column Gene names, Second column Log2 expression values and third column with p-values"
                 )
               ),
               fileInput(
                 "df",
                 "Choose text File",
                 multiple = FALSE,
                 accept = c(
                   "text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv",
                   ".txt"
                 )
               ),
               checkboxInput("pvalue_log", "P Value (-Log10)", value = T),
               textInput("plot_title", "Plot Title"),
               textInput("x_axis_title", "X Axis Title"),
               textInput("y_axis_title", "Y Axis Title"),
               width = 3
             ),
             mainPanel(tabsetPanel(
               tabPanel("plot", plotlyOutput("vp")),
               tabPanel("Input Data", dataTableOutput("table"))
             ))
           )),
  tabPanel("Heatmap plot",
           sidebarLayout(
             sidebarPanel(
               div(
                 style = "font-size: 12px;",
                 helpText(
                   "Text must be in CSV format. First row must have Sample names and first column must have gene names/symbols"
                 )
               ),
               fileInput(
                 "hm_df",
                 "Choose text File",
                 multiple = FALSE,
                 accept = c(
                   "text/csv",
                   "text/comma-separated-values,text/plain",
                   ".csv",
                   ".txt"
                 )
               ),
               textInput("plot_title", "Plot Title"),
               textInput("row_title", "Row Title"),
               textInput("column_title", "Column Title"),
               pickerInput(
                 "other_params",
                 "Other Parameters (click to select)",
                 choices = c(
                   "Row Dendrogram" = "rdend",
                   "Column Dendrogram" = "cdend",
                   "Scale Rows" = "Rnorm",
                   "Scale columns" = "Cnorm"
                 ),
                 multiple = TRUE
               ),
               width = 3
             ),
             mainPanel(tabsetPanel(
               tabPanel("Heat map", plotOutput("hm_plot")),
               tabPanel("Input Data", dataTableOutput("hm_table"))
             ))
           ))
)
##################################################
##############server function#####################
##################################################
server = function(input, output) {
  #####################################
  ### For volcano plot ################
  ### #################################
  ## Create data frane
  df <- reactive({
    df = read.csv(
      input$df$datapath,
      strip.white = T,
      stringsAsFactors = F,
      header = T
    ) %>%
      rename("genes" = 1,
             "fold_change" = 2 ,
             "p_value" = 3) %>%
      mutate(status = ifelse(
        fold_change > 0.6,
        "Up regulated",
        ifelse (fold_change < -0.6, "Down regulated", "Did not change")
      )) %>%
      drop_na() %>%
      arrange(p_value)


    if (input$pvalue_log)
    {
      df$p_value = -log10(df$p_value)
    }

    return (df)
  })
  ## Variable for p-value cutoff, plot, x and y axis titles
  plot_title = reactive({
    if (input$plot_title == "")
    {
      plot_title = "Expression values"
    }
    if (input$plot_title != "")
    {
      plot_title = input$plot_title
    }
    return (plot_title)
  })
  title_y = reactive({
    if (input$pvalue_log && input$y_axis_title == "")
    {
      title_y = paste0("P Value (-Log", tags$sub("10"), ")")
    }
    else if (input$pvalue_log == F && input$y_axis_title == "")
    {
      title_y = "P Value"
    }
    else{
      title_y = input$y_axis_title
    }
    return (title_y)
  })
  title_x = reactive({
    if (input$x_axis_title == "")
    {
      title_x = paste0("Fold change (-Log", tags$sub("2"), ")")
    }
    else {
      title_x = input$x_axis_title
    }
    return (title_x)
  })
  p_cutoff = reactive ({
    if (input$pvalue_log) {
      p_cutoff = -log10(0.01)
    }
    else{
      p_cutoff = 0.01
    }
    return (p_cutoff)
  })
  ## Create ggplot object
  p = reactive({
    ggplot(df(),
           aes(
             x = fold_change,
             y = p_value,
             col = status,
             text = paste(
               df()[, 1],
               "<br>Fold change: ",
               round(fold_change, digits = 2),
               "<br>P value: ",
               round(p_value, digits = 2),
               "<br>status: ",
               status
             )
           )) +
      geom_point() +
      geom_vline(
        xintercept = c(-0.6, 0.6),
        color = "grey",
        linetype = "dashed"
      ) +
      geom_hline(yintercept = p_cutoff(),
                 color = "grey",
                 linetype = "dashed") +
      xlab(title_x()) +
      ylab(title_y()) +
      ggtitle(plot_title()) +
      scale_colour_manual(values = c("grey", "darkred", "darkgreen")) +
      theme_pubr() +
      theme(legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5))
  })
  ## output ggplot in plot tab
  output$vp = renderPlotly({
    req(input$df)
    validate(
      need(ncol(input$df)==3, 'Error!! Input file should have three columns only')
    )
    ggplotly(p(), tooltip = c("text"))
  })
  ## output table in data tab
  output$table = renderDataTable({
    req(input$df)
    validate(
      need(ncol(input$df)==3, 'Error!! Input file should have three columns only')
      )

    df() %>% mutate(p_value = round(p_value, 3),
                    fold_change = round(fold_change, 2))
  })
  #####################################
  ##For heatmap tab####################
  #####################################
  ## Create variables################
  ### Create data frame to load
  hm_df <- reactive({
    read.csv(
      input$hm_df$datapath,
      strip.white = T,
      row.names = 1,
      stringsAsFactors = F,
      header = T
    )
  })
  ### Create heatmap
  ### Output table in data tab
  output$hm_table = renderDataTable({
    req(input$hm_df)
    hm_df()
  })
  ### output heatmap in heatmap tab
  output$hm_plot = renderPlot({
    req(input$hm_df)
    Heatmap(as.matrix(hm_df())[1:25, ])
  })

}

shinyApp(ui = ui, server = server)