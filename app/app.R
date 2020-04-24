library(shiny)
library(shinythemes)
library(shinydashboard)
library(knitr)
library(stringr)
library(pdftools)
library(magick)
library(digest)
library(DT)
library(IRanges)
library(GenomicRanges)
library(GenomeInfoDb)
library(sangerseqR)
source("./tool.R")
toolDir <- "."
dbPrefix <- "./human_g1k_v37_decoy"

# default files ----
outputDir <- "/project/response"
log <- "/project/access.log"
if (! dir.exists(outputDir)) dir.create(outputDir)

# Define UI for data upload app ----
ui <- fluidPage( 
  
    theme = shinytheme("united"),
    # App title ----
    headerPanel(title = "Sanger Analysis"),
    # Sidebar panel for inputs ----
    sidebarPanel(
    
        tags$h4('1. Upload Data:'),
        textInput("position", label = "Version: GRCh37", value = "", placeholder = "1:1000 or 1:100-1000", width = "70%"),
        helpText(tags$span(style="color: red", tags$small("Note: To match best alignment."))),
        fileInput("abif", "Choose chromatogram file (.ab1)", multiple = FALSE, accept = c(".ab1")),
        tags$hr(),
        tags$h4('2. Set Chromatogram Options:'),
        sliderInput("x", strong("Approx. number of bases per row"), min = 10, max = 200, value = 50, step = 10, ticks = FALSE),
        helpText(tags$small("Higher numbers fit more peaks on a single line.")),
        tags$div(style="display: inline-block;vertical-align:top; width: 150px;", numericInput("trim5", strong("5' Trim"), 50)),
        tags$div(style="display: inline-block;vertical-align:top; width: 150px;", numericInput("trim3", strong("3' Trim"), 50)),
        tags$br(style="clear:both"),
        checkboxInput('showtrim', 'Show Trimmed Region', TRUE),
        helpText(tags$small("Removes low quality bases from the ends.")),
        numericInput("ratio", strong('Signal Ratio'), '0.33', min = 0, max = 1, step = .01, width = "30%"),
        helpText(tags$small("Signal Ratio is calculated as peak signal/max signal for each position. Signals above this ratio are called as alternate bases.")),
        tags$hr(),
    
        tags$h4('3. Screenshot:'),
        numericInput("winSize", strong('Window size'), '10', min = 1, max = 100, step = 1, width = "30%"),
        downloadButton("downloadData", "Download"),
        tags$hr(),
    
        tags$small('Created by ', tags$a(href="https://github.com/sanadamakomi", "sanadamakomi", target="_blank")),
        tags$br(),
        tags$small('Hosted by', tags$a(href="http://www.rstudio.com", target="_blank", "Rstudio")),
        tags$br(),
        tags$small('Using sangerseqR version ', textOutput('ssRversion'))
    ),
  
  
    # Main panel for displaying outputs ----
    mainPanel(
      
        # tableset ----
        tabsetPanel(
        
            # chromatogram
            tabPanel("Chromatogram", 
            
                conditionalPanel(
                    condition = "!output$fileUploaded", 
                    tags$p("A Chromatogram will be shown here when a valid sequencing file has been uploaded.")
                ), 
    
                conditionalPanel(
                    condition = "output$fileUploaded", 
                    verbatimTextOutput('checkposMess'),
                    plotOutput('chromatogram')
                )
                
            ),
    
            # variants
            tabPanel("Variants",
            
                conditionalPanel(
                    condition = "!output$fileUploaded", 
                    tags$p("A variant calling result will be shown here when a valid sequencing file has been uploaded.")
                ), 
      
                conditionalPanel(
                    condition = "output$fileUploaded", 
                    DT::dataTableOutput("variantMatrix")
                )
                
            ),
    
            # screenshot
            tabPanel("Screenshot",
      
                conditionalPanel(
                    condition = "!output$filedownload", 
                    tags$p("A screenshot of validated position will be shown here when is validated")
                ), 
             
                conditionalPanel(
                    condition = "output$filedownload", 
                    plotOutput('screenshot')
                )        
            )
    
        )
    )
)


# Define server logic to read selected file ----
server <- function(input, output, session) {

    inputdata <- reactive({
        if (!is.null(input$abif)) {
            return(makeBaseCalls(readsangerseq(input$abif$datapath), input$ratio))
        } else return(NULL)
    })
  
    inputpos <- reactive({ return(formatPos(input$position)) })
  
    h <- reactive({
        if (!is.null(input$abif)) {
            figheight(inputdata(), input$trim5, input$trim3, width=input$x, showtrim=input$showtrim)
        } else return("auto")
    })
  
    calldata <- reactive({
        if (!is.null(input$abif)) {
            return(batchFunc(inputdata(), pos = inputpos(), toolDir = toolDir, dbPrefix = dbPrefix))
        } else return(list(variants=NULL, positions=NULL))
    })
  
    hetcallUpdate <- reactive({
        if (! is.null(input$abif) & ! is.null(calldata()$positions)) {
            return(updateHetcall(inputdata(), winSize = input$winSize, position.df = calldata()$positions))
        } else return(NULL)
    })
  
    ## chromatogram plot
    chroplot <- reactive({
        if( !is.null(input$abif)) {
            chromatogram(inputdata(), showcalls = "both", trim3 = input$trim3, trim5 = input$trim5, width = input$x, showtrim = input$showtrim, showhets = TRUE, cex.base = 2, cex.mtext = 1)
        } else return(NULL)
    })
  
    ## download plot
    screenshotplot <- reactive({
        if (! is.null(input$abif) & ! is.null(calldata()$positions) ) {
            obj <- hetcallUpdate()
            seq.len <- length(obj@primarySeq)
            position.df <- calldata()$positions
            realwin <- position.df$query.loc.end - position.df$query.loc.start + 1 + input$winSize * 2
            plotwin <- 50
            if (input$winSize > 50) {
                plotwin <- realwin + 1 
            }
            format.trim5 <- position.df$trim5 - input$winSize
            format.trim3 <- position.df$trim3 - input$winSize
            if (format.trim5 < 0) format.trim5 <- 0
            if (format.trim3 < 0) format.trim3 <- 0
            chromatogram(obj, width = plotwin, height = 2, trim5 = format.trim5, trim3 = format.trim3, showcalls = "both", cex.base = 2, cex.mtext = 1, showhets = TRUE)
        } else return(NULL)
    })
  
    screenshotdownload <- reactive({
        if (! is.null(input$abif) & ! is.null(calldata()$positions)) {
            obj <- hetcallUpdate()
            seq.len <- length(obj@primarySeq)
            position.df <- calldata()$positions
            realwin <- position.df$query.loc.end - position.df$query.loc.start + 1 + input$winSize * 2
            plotwin <- 50
            if (input$winSize > 50) {
                plotwin <- realwin + 1 
            }
            fileName <- paste(get_time_human(), digest(hetcallUpdate(), algo = "md5"), input$winSize, sep = "_")
            pdfName <- paste0(fileName, ".pdf")
            pngName <- paste0(fileName, ".png")
            chromatogram(obj, width = plotwin, height = 2, trim5 = position.df$trim5 - input$winSize, trim3 = position.df$trim3 - input$winSize, showcalls = "both", cex.base = 2, cex.mtext = 1, showhets = TRUE, filename = file.path(outputDir, pdfName))
            addRect(file.path(outputDir,pdfName), file.path(outputDir, pngName), total.width = seq.len, bin.width = input$winSize)
            return(file.path(outputDir, pngName))
        } else return(NULL)
    })
  
    ## ouput
    ## chromatogram 
  
    output$fileUploaded <- reactive({return(!is.null(input$abif))})
    output$chromatogram <- renderPlot(chroplot(), height = reactive(h()), res = 120)
    output$checkposMess <- renderPrint({
        if(is.null(calldata()$positions)) return("Please input position.")
        else {return(paste0("Position locate in ", calldata()$positions$query.loc.start, " to ", calldata()$positions$query.loc.end))}
    })
    output$h <- reactive(h())
  
    ## screenshot
    output$filedownload <- reactive({return(!is.null(screenshotplot()))})
    output$screenshot <- renderPlot(screenshotplot(), height = 250, res = 120)
    output$downloadData <- downloadHandler(
    
        filename <- function() {
            paste0("screenshot", "_", input$winSize, "_", get_time_human(), ".png")
        },
        
        content <- function(file) {
            file.copy(screenshotdownload(), file)
        },
        
        contentType = "png"
    )
  
    ## data table
    output$variantMatrix <- DT::renderDataTable({
        if (is.null(calldata()$variants)) return(NULL)
        else return(DT::datatable(calldata()$variants[,1:12,drop=FALSE],
            extensions = 'Buttons',
            rownames = FALSE,
            options = list(
                orderClasses = TRUE,
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
            )
        ))
    })
  
    ## observe
    observe({
        c(input$ratio, input$trim5, input$trim3)
        if(! is.null(input$abif))
        updateTabsetPanel(session, "maintabset", selected = "Chromatogram")
    })
  
    observe({
        variants <- calldata()$variants
        if(! is.null(variants)) {
            updateTabsetPanel(session, "maintabset", selected = "Variants")
        }
    })
  
    # observe({
        # sch <- screenshotplot()
        # if(! is.null(sch)) {
            # updateTabsetPanel(session, "maintabset", selected = "Screenshot")
        # }
    # })
  
    #logging
    observe({
        if(!is.null(inputdata())) {
            isolate({
                alog <- file(log, "a")
                cat(format(Sys.time(), "%Y-%b-%d_%H-%M-%S"), file=alog)
                cat("\t", file=alog)
                cat(input$seq$name, file=alog)
                cat("\n", file=alog)
                close(alog)
            })
        }
    })
    
    output$ssRversion <- renderText(as.character(packageVersion("sangerseqR")))
}

# Create Shiny app ----
shinyApp(ui, server)