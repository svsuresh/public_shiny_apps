library(shiny)
library(plotly)
library(ggplot2)
library(ggpubr)
library(DT)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))



options(scipen = 999)

#UI.R
ui = navbarPage("Expression plots",
                tabPanel("Volcano plot",
                         sidebarLayout(
                           sidebarPanel(
                             fileInput(
                               "df",
                               "Choose text File",
                               multiple = FALSE,
                               accept = c(
                                 "text/tsv",
                                 "text/tab-separated-values,text/plain",
                                 "text/comma-separated-values,text/plain",
                                 ".tsv",
                                 ".csv",
                                 ".txt"
                               )
                             ),
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

  output$vp = renderPlotly({
    req(input$df)

    p=ggplot(df(), aes(
      x = fold_change,
      y = -log10(p_value),
      col = status,
      text = paste(df()[,1], "<br>Fold change: ", round(fold_change, digits = 2), "<br>P value: ", round(-log10(p_value), digits = 2), "<br>status: ", status))) +
      geom_point() +
      geom_vline(xintercept = c(-0.6, 0.6),
                 color = "grey",
                 linetype = "dashed") +
      geom_hline(yintercept = 2,
                 color = "grey",
                 linetype = "dashed") +
      xlab("Fold change (log2)") +
      ylab("Adjusted p-value (-log10)") +
      scale_colour_manual(values = c("grey", "darkred", "darkgreen"))+
      theme_pubr() +
      theme(legend.title = element_blank() )

    ggplotly(p, tooltip=c("text"))

  })

  output$table = renderDataTable({
    req(input$df)
    df()
  })
}
shinyApp(ui = ui, server = server)
