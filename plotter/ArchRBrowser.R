library(ArchR)
library(ggplot2)
library(shiny)
library(argparse)

isColor <- function(x = NULL) {
                unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), 
                  error = function(e) FALSE)))
}

MyArchRBrowser <- function(
  ArchRProj = NULL,
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  minCells = 25,
  baseSize = 10,
  borderWidth = 0.5,
  tickWidth = 0.5,
  facetbaseSize = 12,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  browserTheme = "cosmo",
  threads = getArchRThreads(),
  verbose = TRUE,
  logFile = createLogFile("ArchRBrowser"),
  host="10.100.44.33",
  port=5013
  ){

  geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)

  ArchR:::.startLogging(logFile=logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "ArchRBrowser Input-Parameters", logFile = logFile)

  ArchR:::.requirePackage("shiny", installInfo = 'install.packages("shiny")')
  ArchR:::.requirePackage("rhandsontable", installInfo = 'install.packages("rhandsontable")')

  #Determine Grouping Methods
  ccd <- getCellColData(ArchRProj)
  discreteCols <- lapply(seq_len(ncol(ccd)), function(x){
   ArchR:::.isDiscrete(ccd[, x])
  }) %>% unlist %>% {colnames(ccd)[.]}
  if("Clusters" %in% discreteCols){
    selectCols <- "Clusters"
  }else{
    selectCols <- "Sample"
  }

  #Extend where upstream can be negative for browser
  extendGR2 <-  function(gr = NULL, upstream = NULL, downstream = NULL){
    ArchR:::.validInput(input = gr, name = "gr", valid = c("GRanges"))
    ArchR:::.validInput(input = upstream, name = "upstream", valid = c("integer"))
    ArchR:::.validInput(input = downstream, name = "downstream", valid = c("integer"))
    #Get Info From gr
    st <- start(gr)
    ed <- end(gr)
    #https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
    isMinus <- BiocGenerics::which(strand(gr) == "-")
    isOther <- BiocGenerics::which(strand(gr) != "-")
    #Forward
    st[isOther] <- st[isOther] - upstream
    ed[isOther] <- ed[isOther] + downstream
    #Reverse
    ed[isMinus] <- ed[isMinus] + upstream
    st[isMinus] <- st[isMinus] - downstream
    #If Any extensions now need to be flipped.
    end(gr) <- pmax(st, ed)
    start(gr) <- pmin(st, ed)
    return(gr)
  }


  #####################
  #Shiny App UI
  #####################
  if(!requireNamespace("shinythemes", quietly = TRUE)){
    ArchR:::.logMessage("shinythemes not found! To see a nice theme use :\n\tinstall.packages('shinythemes')\nContinuing wihtout shinythemes!", verbose = verbose, logFile = logFile)
    theme <- NULL
  }else{
    theme <- shinythemes::shinytheme(browserTheme)
  }

  ui <- fluidPage(
    theme = theme,
    titlePanel(
        h1(div(HTML(paste0("<b>ArchR Browser v1 : nCells = ", formatC(nCells(ArchRProj), format="f", big.mark = ",", digits=0), "</b>"))), align = "left")
    ),
    sidebarLayout(
      sidebarPanel(
        tags$head(
          tags$style(HTML("hr {border-top: 1px solid #000000;}"))
        ),
        tags$head(
          tags$style(HTML('#exitButton{background-color:#D60000}'))
        ), 
        tags$head(
          tags$style(HTML('#restartButton{background-color:#02A302}'))
        ),
        tags$head(
          tags$style(HTML('#plot_height{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#plot_width{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#ymax{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#tile_size{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#range_min{height: 35px}'))
        ),
        tags$head(
          tags$style(HTML('#range_max{height: 35px}'))
        ),
        actionButton(inputId = "exitButton", label = "Exit Browser", icon = icon("times-circle")),
        br(),
        br(),
        actionButton(inputId = "restartButton", label = "Plot Track!", icon = icon("play-circle")),
        #div(style="display:inline-block;width:50%;text-align: center;",actionButton("exitButton", label = "Exit Browser", icon = icon("paper-plane"))),
        #br(),
        #div(style="display:inline-block;width:50%;text-align: center;",actionButton("restartButton", label = "Plot Track!", icon = icon("paper-plane"))),
        br(),
        br(),
        selectizeInput("name",
                       label = "Gene Symbol",
                       choices = as.vector(geneAnnotation$genes$symbol)[!is.na(as.vector(geneAnnotation$genes$symbol))],
                       multiple = FALSE,
                       options = list(placeholder = 'Select a Center'),
                       selected = "CD4"
                      ),
        selectizeInput("grouping",
           label = "groupBy",
           choices = discreteCols,
           multiple = FALSE,
           options = list(placeholder = 'Select Grouping'),
           selected = selectCols
        ),
        sliderInput("range", "Distance From Center (kb):", min = -250, max = 250, value = c(-50,50)),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("range_min", "Distance (-kb):", min = -250, max = 250, value = -50),
          numericInput("range_max", "Distance (+kb):", min = -250, max = 250, value = 50)
        ),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("tile_size", "TileSize:", min = 10, max = 5000, value = 250),
          numericInput("ymax", "Y-Max (0,1):", min = 0, max = 1, value = 0.99)
        ),
        hr(),
        downloadButton(outputId = "down", label = "Download the Track!"),
        br(),
        br(),
        splitLayout(cellWidths = c("50%","50%"),
          numericInput("plot_width", "Width", min = 0, max = 250, value = 8),
          numericInput("plot_height", "Height", min = 0, max = 250, value = 12)
        ),
        width = 2, height = "750px", position = "left"
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Plot",
            plotOutput(outputId = "ATAC", width= "800px", height = "725px")
          ),
          tabPanel("Additional Params",
            br(),
            br(),
            selectizeInput("normATAC",
               label = "normMethod",
               choices = c("ReadsInTSS", "ReadsInPromoter", "nFrags", "None"),
               multiple = FALSE,
               options = list(placeholder = 'Select NormMethod'),
               selected = "ReadsInTSS"
            ),
            br(),
            div(HTML("<b>Group Metadata</b>")),
            rHandsontableOutput("Metadata", width= "1085px",height = "800px")
          )
        )
      )
    )
  )

  #####################
  #Shiny App Server
  #####################
  server <- function(
        input = input, 
        output = output, 
        session = session
    ){

    output$Metadata <- renderRHandsontable({
        groups <- gtools::mixedsort(unique(ccd[,input$grouping]))
        mdata <- data.frame(
          groupBy = input$grouping,
          include = rep(TRUE,length(groups)), 
          group = groups, 
          color = paletteDiscrete(values = groups)[groups], 
          nCells = as.vector(table(ccd[,input$grouping])[groups]),
          medianTSS = getGroupSummary(ArchRProj = ArchRProj, select = "TSSEnrichment", summary = "median", groupBy = input$grouping)[groups],
          medianFragments = getGroupSummary(ArchRProj = ArchRProj, select = "nFrags", summary = "median", groupBy = input$grouping)[groups],
          stringsAsFactors = FALSE
        )
        rownames(mdata) <- NULL
        rhandsontable(mdata)
      })

    #Update Sliders
    observeEvent(input$range_min, {
      updateSliderInput(session, "range",
                        value = c(input$range_min,max(input$range)))
    })
    
    observeEvent(input$range_max, {
      updateSliderInput(session, "range",
                        value = c(input$range_min,input$range_max))
    })

    observeEvent(input$range , {

      updateNumericInput(session, "range_min", value = min(input$range))
      updateNumericInput(session, "range_max", value = max(input$range))

    }, priority = 200)

    output$checkbox <- renderUI({
      choice <- gtools::mixedsort(unique(ccd[,input$grouping,drop=TRUE]))
      checkboxGroupInput("checkbox","Select Groups", choices = choice, selected = choice)
    })    

    #################################
    # Inputs that cause re-plotting
    #################################
    #toListen <- reactive({
    #  list(inputrestartButton,inputname)
    #})

    restartFN <- observeEvent(input$restartButton, {

      if (input$name == ""){
          
          output$ATAC <- renderPlot({
            p <- ggplot() +
                xlim(c(-5,5)) + ylim(c(-5,5)) +
              geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
            print(p)
          })

      }else{
        output$ATAC <- renderPlot({

          withProgress(message = 'Plotting', style = "notification", value = 0, {

            #Get Region if Gene Symbol
            region <- geneAnnotation$genes

            if(tolower(input$name) %ni% tolower(mcols(region)$symbol)){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Please Supply\nA Valid Gene!")) + theme_void()
              return(print(p))
            }

            region <- region[which(tolower(mcols(region)$symbol) %in% tolower(input$name))]
            region <- region[order(match(tolower(mcols(region)$symbol), tolower(input$name)))]
            region1 <- resize(region, 1, "start")
            strand(region1) <- "*"

            #Extend Region
            #region <- extendGR(region, upstream = -min(inputrange)∗1000,downstream=max(inputrange) * 1000)
            #Pre-Load full window for even faster plotting
            region <- extendGR2(region1, upstream = 250000, downstream = 250000)
            tmpArchRRegion <<- extendGR2(region1, 
              upstream = -min(isolate(input$range)) * 1000, 
              downstream = max(isolate(input$range)) * 1000
            )
            region <- tmpArchRRegion

            setProgress(0.1)

            #User Inputs
            groupBy <- isolate(input$grouping)

            groupDF <- tryCatch({
              isolate(hot_to_r(input$Metadata))
            },error=function(x){
              groups <- gtools::mixedsort(unique(ccd[,isolate(input$grouping)]))
              mdata <- data.frame(
                groupBy = input$grouping,
                include = rep(TRUE,length(groups)), 
                group = groups, 
                color = paletteDiscrete(values = groups)[groups], 
                nCells = as.vector(table(ccd[,input$grouping])[groups]),
                medianTSS = getGroupSummary(ArchRProj = ArchRProj, select = "TSSEnrichment", summary = "median", groupBy = input$grouping)[groups],
                medianFragments = getGroupSummary(ArchRProj = ArchRProj, select = "nFrags", summary = "median", groupBy = input$grouping)[groups],
                stringsAsFactors = FALSE
              )
              rownames(mdata) <- NULL
              mdata
            })

            if(groupDF$groupBy[1] != groupBy){
              groups <- gtools::mixedsort(unique(ccd[,isolate(input$grouping)]))
              groupDF <- data.frame(
                groupBy = groupBy,
                include = rep(TRUE,length(groups)), 
                group = groups, 
                color = paletteDiscrete(values = groups)[groups], 
                stringsAsFactors = FALSE
              )
              rownames(groupDF) <- NULL
            }

            useGroups <- groupDF[groupDF[,"include"],"group"]


            if(!all(isColor(groupDF[groupDF[,"include"], "color"]))){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Colors from Metadata\n not real R colors!")) + theme_void()
              return(print(p))               
            }

            if(any(useGroups %ni% ccd[, groupBy])){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Groups from Metadata\n not present in groupBy!")) + theme_void()
              return(print(p)) 
            }

            pal <- groupDF[groupDF[,"include"], "color"]
            names(pal) <- groupDF[groupDF[,"include"], "group"]

            ylim <- c(0, isolate(input$ymax))
            normMethod <- isolate(input$normATAC)
            tileSize <- isolate(input$tile_size)

            p <- ArchR:::.bulkTracks(
                ArchRProj = ArchRProj, 
                region = region, 
                tileSize = tileSize, 
                useGroups = useGroups,
                groupBy = groupBy,
                threads = threads, 
                minCells = minCells,
                ylim = ylim,
                baseSize = baseSize,
                borderWidth = borderWidth,
                tickWidth = tickWidth,
                facetbaseSize = facetbaseSize,
                normMethod = normMethod,
                geneAnnotation = geneAnnotation,
                title = "",
                pal = pal, 
                tstart = NULL,
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))

            #p <- p + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

            tmpArchRP <<- p

            setProgress(0.5)

            if(!is.null(features)){

              f <- ArchR:::.featureTracks(
                  features = features, 
                  region = tmpArchRRegion,
                  facetbaseSize = facetbaseSize, 
                  hideX = TRUE, 
                  title = "Peaks",
                  logFile = logFile
                ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

              #f <- f + .suppressAll(scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

            }

            if(!is.null(loops)){

              l <- ArchR:::.loopTracks(
                loops = loops, 
                region = tmpArchRRegion, 
                facetbaseSize = facetbaseSize,
                hideX = TRUE, 
                hideY = TRUE,
                title = "Loops",
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
            }

            setProgress(0.6)

            g <- ArchR:::.geneTracks(
              geneAnnotation = geneAnnotation, 
              region = tmpArchRRegion, 
              facetbaseSize = facetbaseSize,
              labelSize = 3,
              title = "Genes",
              logFile = logFile
            ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

            #g <- .suppressAll(g + scale_x_continuous(limits = c(start(tmpArchRRegion), end(tmpArchRRegion)), expand = c(0,0)))

            setProgress(0.8)

            if(!is.null(loops)){
              if(!is.null(features)){
                suppressWarnings(print(ggAlignPlots(p, f, l, g, sizes = c(10, 1.5, 3, 4),type = "v", draw = TRUE)))
              }else{
                suppressWarnings(print(ggAlignPlots(p, l, g, sizes = c(10, 3, 4),type = "v", draw = TRUE)))
              }
            }else{
              if(!is.null(features)){
                suppressWarnings(print(ggAlignPlots(p, f, g, sizes = c(10, 2, 4),type = "v", draw = TRUE)))
              }else{
                suppressWarnings(print(ggAlignPlots(p, g, sizes = c(10, 4),type = "v", draw = TRUE)))
              }
            }

            setProgress(1)

          })

        })
      }    

    })

    #######################################
    # When Download Is Initiated
    #######################################

    # downloadHandler contains 2 arguments as functions, namely filename, content
    output$down <- downloadHandler(

      filename <- function(){
        paste0("ArchRBrowser-",input$name,"-",seqnames(tmpArchRRegion)[1],":",start(tmpArchRRegion)[1],"-",end(tmpArchRRegion)[1],".pdf")
      },

      # content is a function with argument file. content writes the plot to the device
      content = function(file) {
        
        withProgress(message = 'Creating PDF', style = "notification", value = 0, {
         
          if(!exists("tmpArchRP")){

            #User Inputs
            groupBy <- isolate(input$grouping)
            groupDF <- isolate(hot_to_r(input$Metadata))
            useGroups <- groupDF[groupDF[,"include"],"group"]

            isColor <- function(x = NULL) {
                unlist(lapply(x, function(y) tryCatch(is.matrix(col2rgb(y)), 
                  error = function(e) FALSE)))
            }

            if(!all(isColor(groupDF[groupDF[,"include"], "color"]))){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Colors from Metadata\n not real R colors!")) + theme_void()
              return(print(p))               
            }

            if(any(useGroups %ni% ccd[, groupBy])){
              p <- ggplot() +
                  xlim(c(-5,5)) + ylim(c(-5,5)) +
                geom_text(size=20, aes(x = 0, y = 0, label = "Error Groups from Metadata\n not present in groupBy!")) + theme_void()
              return(print(p)) 
            }

            pal <- groupDF[groupDF[,"include"], "color"]
            names(pal) <- groupDF[groupDF[,"include"], "group"]

            ylim <- c(0, isolate(input$ymax))
            normMethod <- isolate(input$normATAC)
            tileSize <- isolate(input$tile_size)

            p <- ArchR:::.bulkTracks(
                ArchRProj = ArchRProj, 
                region = tmpArchRRegion, 
                tileSize = tileSize, 
                useGroups = useGroups,
                groupBy = groupBy,
                threads = threads, 
                minCells = minCells,
                ylim = ylim,
                baseSize = baseSize,
                borderWidth = borderWidth,
                tickWidth = tickWidth,
                facetbaseSize = facetbaseSize,
                normMethod = normMethod,
                geneAnnotation = geneAnnotation,
                title = "",
                pal = pal, 
                tstart = NULL,
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))

          }else{

            print("Using previous ggplot")

            p <- tmpArchRP

          }

          setProgress(0.5)

          if(!is.null(features)){

            f <-ArchR:::.featureTracks(
                features = features, 
                region = tmpArchRRegion,
                facetbaseSize = facetbaseSize, 
                hideX = TRUE, 
                title = "Peaks",
                logFile = logFile
              ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

          }

          setProgress(0.6)

          if(!is.null(loops)){

            l <- ArchR:::.loopTracks(
              loops = loops, 
              region = tmpArchRRegion, 
              facetbaseSize = facetbaseSize,
              hideX = TRUE, 
              hideY = TRUE,
              title = "Loops",
              logFile = logFile
            ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
          }

          setProgress(0.7)

          g <- ArchR:::.geneTracks(
            geneAnnotation = geneAnnotation, 
            region = tmpArchRRegion, 
            facetbaseSize = facetbaseSize,
            labelSize = 3,
            title = "Genes",
            logFile = logFile
          ) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))

          setProgress(0.8)

          pdf(file = file, width = input$plot_width, height = input$plot_height)
          
          if(!is.null(loops)){
            if(!is.null(features)){
              a <- suppressWarnings(ggAlignPlots(p, f, l, g, sizes = c(10, 1.5, 3, 4),type = "v", draw = FALSE))
            }else{
              a <- suppressWarnings(ggAlignPlots(p, l, g, sizes = c(10, 3, 4),type = "v", draw = FALSE))
            }
          }else{
            if(!is.null(features)){
              a <- suppressWarnings(ggAlignPlots(p, f, g, sizes = c(10, 2, 4),type = "v", draw = FALSE))
            }else{
              a <- suppressWarnings(ggAlignPlots(p, g, sizes = c(10, 4),type = "v", draw = FALSE))
            }
          }
          
          suppressWarnings(grid::grid.draw(a))
          dev.off()

          setProgress(1)

        })

    })


    exitFN <- observeEvent(input$exitButton, {
      if(exists("tmpArchRRegion")){
	      ArchR:::.suppressAll(rm(tmpArchRRegion))
      }
      if(exists("tmpArchRP")){
	      ArchR:::.suppressAll(rm(tmpArchRP))
      }
      shiny::stopApp()
    })

  }

  shiny::runApp(list(ui = ui, server = server),host=host,port=port, launch.browser = FALSE)

}

###########################################################
#parser <- ArgumentParser(description='To show browser track in web')
#parser$add_argument("--project",
#                    type="character",
#                    default=NULL,
#                    help="the project path  of ArchR")

#parser$add_argument("--host",
#                    type="character",
#                    default="10.100.44.33",
#                    help="host ip address")


#parser$add_argument("--port",
#                    type="integer",
#                    default="5013",
#                    help="address port")

#parser$add_argument("--jobs",
#                    type="integer",
#                    default=16,
#                    help="number threads")

#args <- parser$parse_args()
#message("INFO : Loading dataset ...")
#addArchRThreads(args$jobs)
#projHeme=loadArchRProject(args$project)
#message("INFO : Lauch browser ...")
#MyArchRBrowser(projHeme,host=args$host,port=args$port)
#message("INFo : Done !")


