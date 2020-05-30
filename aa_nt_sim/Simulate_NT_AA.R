options(repos = BiocManager::repositories())
library(Biostrings)
library(shiny)
library(shinythemes)


##UI.R
ui = fluidPage(
    theme = shinytheme("flatly"),
    titlePanel("DNA/Protein Sequence simulation"),
    sidebarLayout(
        sidebarPanel(
            selectInput(
                "bio_molecule",
                "Select DNA/Protein",
                choices = c("DNA", "Protein")
            ),
            numericInput("seq_num", "Number of sequences you want", 2, min = 1),
            numericInput("seq_len", "Length of each sequence you want", 10, min = 1),
            width = 3
        ),
        mainPanel(
            verbatimTextOutput("sequence"),
            tags$head(
                tags$style("#sequence{overflow-y:scroll; max-height: 350px;}")
            )
            ,
            downloadLink('download_seq', 'Download')
        )
    )
)

##Server.R
server = function(input, output, session) {
    stringset = reactive({
        if (input$bio_molecule == "DNA") {
            stringset = DNAStringSet
        }
        if (input$bio_molecule == "Protein") {
            stringset = AAStringSet
        }
        return (stringset)
    })

    baseset = reactive({
        if (input$bio_molecule == "DNA") {
            baseset = DNA_BASES
        }
        if (input$bio_molecule == "Protein") {
            baseset = AA_STANDARD
        }
        return (baseset)
    })


    name = reactive({
        if (input$bio_molecule == "DNA") {
            name = "DNA_Sequence"
        }
        if (input$bio_molecule == "Protein") {
            name = "AA_Sequence"
        }
        return (name)
    })

    DNAseq = reactive({
        function(seqnum, seqlen) {
            test = stringset()(sapply (1:seqnum, function (x)
                paste(
                    sample(baseset(), seqlen, replace = T), collapse = ""
                )))
            names(test) = paste(name(), letters[1:seqnum], sep = "_")
            return(test)
        }
    })


    output$sequence = renderText({
        session$userData$text1 = DNAseq()(input$seq_num, input$seq_len)
        paste(
            paste(">", names(session$userData$text1)),
            session$userData$text1,
            sep = "\n",
            collapse = "\n"
        )
    })

    output$download_seq <- downloadHandler(
        filename = function() {
            paste0(name(), ".fa")
        },
        content = function(file) {
            writeXStringSet(session$userData$text1, file, format = "fasta")
        }
    )

}


shinyApp(ui, server)
#rm(list=ls())
