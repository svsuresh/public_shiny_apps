library(shiny)
library(plotly)
library(ggplot2)
library(ggpubr)
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))

options(scipen = 999)

#UI.R
ui = navbarPage("Expression plots",
                tabPanel("Volcano plot",
                         sidebarLayout(
                           sidebarPanel(
                             div(
                               style = "font-size: 11px;",
                               wellPanel(
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
                             checkboxInput("pvalue_log", "Log (-10) P-values", value = T),
                             textInput("plot_title", "Plot Title"),
                             textInput("x_axis_title", "X Axis Title"),
                             textInput("y_axis_title", "Y Axis Title"),
                             width = 3
                           ),
                           mainPanel(tabsetPanel(
                             tabPanel("plot", plotlyOutput("vp")),
                             tabPanel("Input Data", dataTableOutput("table"))
                           ))
                         )))
server = function(input, output) {
  df <- reactive({
    read.csv(
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
  })
  title_y = reactive({
    if (input$y_axis_title == "")
    {
      title_y = paste0('P-value (-log', tags$sub("10"), ')')
    }
    if (input$x_axis_title != "")
    {
      title_y = as.character(input$y_axis_title)
    }
    return (title_y)
  })

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

  title_x = reactive({
    if (input$x_axis_title == "")
    {
      title_x = paste0('Fold change (log', tags$sub("2"), ')')
    }
    if (input$x_axis_title != "")
    {
      title_x = as.character(input$x_axis_title)
    }
    return (title_x)
  })

  p = reactive({
    ggplot(df(),
           aes(
             x = fold_change,
             y = -log10(p_value),
             col = status,
             text = paste(
               df()[, 1],
               "<br>Fold change: ",
               round(fold_change, digits = 2),
               "<br>P value: ",
               round(-log10(p_value), digits = 2),
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
      geom_hline(yintercept = 2,
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

  output$vp = renderPlotly({
    req(input$df)
    ggplotly(p(), tooltip = c("text"))
  })

  output$table = renderDataTable({
    req(input$df)
    df()
  })
}
shinyApp(ui = ui, server = server)