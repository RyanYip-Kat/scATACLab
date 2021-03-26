library(ArchR)
library(ggplot2)
library(ggrepel)

getBrowserTrack <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  facet=FALSE,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  pal = NULL,
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = "BrowserTrack"
  ){
  
  plotSummary <- tolower(plotSummary)
  names(sizes) <- plotSummary
  sizes <- sizes[order(plotSummary)]
  plotSummary <- plotSummary[order(plotSummary)]

  geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)
  tstart <- Sys.time()
  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      print(region)
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- .validGRanges(region)
  names(region)=paste(seqnames(region),start(region),end(region),sep="-")
  if(is.null(geneSymbol)){
    useMatrix <- NULL
  }

  if(!is.null(useMatrix)){
    featureMat <- ArchR:::.getMatrixValues(
      ArchRProj = ArchRProj,
      matrixName = useMatrix,
      name = mcols(region)$symbol
    )
    if(log2Norm){
      featureMat <- log2(featureMat + 1) 
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[rownames(featureMat), 1]
  }
  dfList <- list()
  ggList <- list()
  for(x in seq_along(region)){
    plotList<-list()
    DFList<-list()
    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){
       cat(sprintf("INFO : %d of %d --- bulktrack\n",x,length(region)))
       bulktrack <- .bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = minCells,
        pal = pal,
	facet=facet,
        ylim = ylim,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        facetbaseSize = facetbaseSize,
        normMethod = normMethod,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) #+ theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
       plotList[["bulktrack"]]=bulktrack[[1]]+theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
       DFList[["bulktrack"]]=bulktrack[[2]]
    }
    
    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("sctrack" %in% tolower(plotSummary)){
      cat(sprintf("INFO : %d of %d --- sctrack\n",x,length(region)))
      sctrack <- .scTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = 5,
        maxCells = scCellsMax,
        pal = pal,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        scTileSize = scTileSize,
        facetbaseSize = facetbaseSize,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) #+ theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
      plotList[["sctrack"]]=sctrack[[1]]+theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
      DFList[["sctrack"]]=sctrack[[2]]
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){
      if(!is.null(features)){
        featuretrack <- .featureTracks(
            features = features, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "Peaks",
            logFile = logFile) #+ theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
        plotList[["featuretrack"]]=featuretrack[[1]]+theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
	DFList[["featuretrack"]]=featuretrack[[2]]
      }
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){
      if(!is.null(loops)){
	cat(sprintf("INFO : %d of %d --- looptrack\n",x,length(region)))
        looptrack <- .loopTracks(
            loops = loops, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            hideY = TRUE,
            title = "Loops",
            logFile = logFile) #+ theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
        plotList[["looptrack"]]=looptrack[[1]]+theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
	DFList[["looptrack"]]=looptrack[[2]]
      }
    }

    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){
      cat(sprintf("INFO : %d of %d --- genetrack\n",x,length(region)))
      genetrack <- .geneTracks(
        geneAnnotation = geneAnnotation, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes",
        logFile = logFile)# + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      plotList[["genetrack"]]=genetrack[[1]]+theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      DFList[["genetrack"]]=genetrack[[2]]
    }
    ##########################################################
    # Time to plot
    ##########################################################
    #plotSummary <- tolower(plotSummary)
    #names(sizes) <- plotSummary
    #sizes <- sizes[order(plotSummary)]
    #plotSummary <- plotSummary[order(plotSummary)]

    # nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    # if(any(nullSummary)){
    #   sizes <- sizes[-which(nullSummary)]
    # }
    sizes <- sizes[tolower(names(plotList))]
    p=ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE)
    dfList[[names(region)[x]]]=DFList
    ggList[[names(region)[x]]]=p 
  }
  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
    names(dfList)=names(ggList)
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }

  list("ggList"=ggList,"DFList"=dfList)

}
    
#######################################################
# Bulk Aggregated ATAC Track Methods
#######################################################
.bulkTracks <- function(
  ArchRProj = NULL, 
  region = NULL, 
  tileSize = 100, 
  minCells = 25,
  groupBy = "Clusters",
  useGroups = NULL,
  normMethod = "ReadsInTSS",
  threads = 1, 
  ylim = NULL,
  baseSize = 7,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  facet=TRUE,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){

  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  df <- .groupRegionSumArrows(
    ArchRProj = ArchRProj, 
    groupBy = groupBy, 
    normMethod = normMethod,
    useGroups = useGroups,
    minCells = minCells,
    region = region, 
    tileSize = tileSize, 
    threads = threads,
    verbose = verbose,
    logFile = logFile
  )

  ######################################################
  # Plot Track
  ######################################################
  if(!is.null(ylim)){
    ylim <- quantile(df$y, ylim)
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }else{
    ylim <- c(0,quantile(df$y, probs=c(0.999)))
    df$y[df$y < ylim[1]] <- ylim[1]
    df$y[df$y > ylim[2]] <- ylim[2]
  }
  uniqueGroups <- gtools::mixedsort(unique(paste0(df$group)))
  if(!is.null(useGroups)){
    uniqueGroups <- unique(useGroups)
  }
  df$group <- factor(df$group, levels = uniqueGroups)
  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)

  allGroups <- gtools::mixedsort(unique(getCellColData(ArchRProj = ArchRProj, select = groupBy, drop = TRUE)))

  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = allGroups))
  }
  
  #Plot Track
  if(facet){
    p <- ggplot(df, aes_string("x","y", color = "group", fill = "group")) + 
      geom_area(stat = "identity") + 
      facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
      ylab(sprintf("Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_y_continuous(limits = ylim, expand = c(0,0)) +
      theme_ArchR(baseSize = baseSize,
                baseRectSize = borderWidth,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
      theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize, 
            color = "black", 
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color="black")) +
       guides(fill = FALSE, colour = FALSE) + ggtitle(title)
  }else{
    p <- ggplot(df, aes_string("x","y", group="group",color = "group", fill = "group")) +
      geom_area(stat = "identity") +
      #facet_wrap(facets = ~group, strip.position = 'right', ncol = 1) +
      ylab(sprintf("Coverage\n(Norm. ATAC Signal Range (%s-%s) by %s)", round(min(ylim),2), round(max(ylim),2), normMethod)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_y_continuous(limits = ylim, expand = c(0,0)) +
      theme_ArchR(baseSize = baseSize,
                baseRectSize = borderWidth,
                baseLineSize = tickWidth,
                legendPosition = "right",
                axisTickCm = 0.1) +
      theme(panel.spacing= unit(0, "lines"),
          axis.title.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text = element_text(
            size = facetbaseSize,
            color = "black",
            margin = margin(0,0.35,0,0.35, "cm")),
            strip.text.y = element_text(angle = 0),
          strip.background = element_rect(color="black")) + ggtitle(title)
       #guides(fill = FALSE, colour = FALSE) + ggtitle(title)
  }

  return(list("bulkTracksPlot"=p,"bulkTracksDF"=df))

}

##############################################################################
# Create Average Tracks from Arrows
##############################################################################
.groupRegionSumArrows <- function(
  ArchRProj = NULL,
  useGroups = NULL,
  groupBy = NULL,
  region = NULL,
  tileSize = NULL,
  normMethod = NULL,
  verbose = FALSE,
  minCells = 25,
  maxCells = 500,
  threads = NULL,
  logFile = NULL
  ){

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    cat(sprintf("Getting Region From Arrow Files %s of %s\n", i, length(ArrowFiles)), logFile = logFile)
    tryCatch({
      .regionSumArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    })
  }, threads = threads) %>% Reduce("+" , .)

  #Plot DF
  df <- data.frame(which(groupMat > 0, arr.ind=TRUE))
  df$y <- groupMat[cbind(df[,1], df[,2])]

  #Minus 1 Tile Size
  dfm1 <- df
  dfm1$row <- dfm1$row - 1
  dfm1$y <- 0

  #Plus 1 Size
  dfp1 <- df
  dfp1$row <- dfp1$row + 1
  dfp1$y <- 0

  #Create plot DF
  df <- rbind(df, dfm1, dfp1)
  df <- df[!duplicated(df[,1:2]),]
  df <- df[df$row > 0,]
  df$x <- regionTiles[df$row]
  df$group <- uniqueGroups[df$col]

  #Add In Ends
  dfs <- data.frame(
    col = seq_along(uniqueGroups), 
    row = 1, 
    y = 0,
    x = start(region),
    group = uniqueGroups
  )

  dfe <- data.frame(
    col = seq_along(uniqueGroups),
    row = length(regionTiles),
    y = 0,
    x = end(region),
    group = uniqueGroups
  )
  
  #Final output
  plotDF <- rbind(df,dfs,dfe)
  plotDF <- df[order(df$group,df$x),]
  plotDF <- df[,c("x", "y", "group")]
  
  #Normalization 
  g <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(tolower(normMethod) == "readsintss"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "readsinpromoter"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "nfrags"){
      v <- getCellColData(ArchRProj, normMethod, drop = TRUE)
      groupNormFactors <- unlist(lapply(split(v, g), sum))
  }else if(tolower(normMethod) == "ncells"){
      groupNormFactors <- table(g)
  }else if(tolower(normMethod) == "none"){
      groupNormFactors <- rep(10^4, length(g))
      names(groupNormFactors) <- g
  }else{
    stop("Norm Method Not Recognized : ", normMethod)
  }

  #Scale with Norm Factors
  scaleFactors <- 10^4 / groupNormFactors
  matchGroup <- match(paste0(plotDF$group), names(scaleFactors))
  plotDF$y <- plotDF$y * as.vector(scaleFactors[matchGroup])

  return(plotDF)

}

.regionSumArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){
  
  cellFragsRegion <- ArchR:::.getFragsFromArrow(
      ArrowFile = ArrowFile, 
      chr = paste0(seqnames(region)), 
      cellNames = cellNames, 
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)
  
  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)
  
  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)
  
  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)
  
  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = c(matchID[ids], matchID[ide]),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames
  
  mat@x[mat@x > 1] <- 1

  #Create Group Matrix
  groupMat <- matrix(0, nrow = length(regionTiles), ncol = length(uniqueGroups))
  colnames(groupMat) <- uniqueGroups
  uniqueGroups <- uniqueGroups[uniqueGroups %in% unique(cellGroups)]
  for(i in seq_along(uniqueGroups)){
    groupMat[,uniqueGroups[i]] <- Matrix::rowSums(mat[,which(cellGroups == uniqueGroups[i]),drop=FALSE])
  }

  return(groupMat)

}

#######################################################
# Gene Tracks
#######################################################
.geneTracks <- function(
  geneAnnotation = NULL, 
  region = NULL, 
  baseSize = 9, 
  borderWidth = 0.4, 
  title = "Genes",
  geneWidth = 2, 
  exonWidth = 4, 
  labelSize = 2,
  facetbaseSize,
  colorMinus = "dodgerblue2",
  colorPlus = "red",
  logFile = NULL
  ){


  #only take first region
  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  genes <- sort(sortSeqlevels(geneAnnotation$genes), ignore.strand = TRUE)
  exons <- sort(sortSeqlevels(geneAnnotation$exons), ignore.strand = TRUE)
  genesO <- data.frame(subsetByOverlaps(genes, region, ignore.strand = TRUE))

  if(nrow(genesO) > 0){

    #Identify Info for Exons and Genes
    exonsO <- data.frame(subsetByOverlaps(exons, region, ignore.strand = TRUE))
    exonsO <- exonsO[which(exonsO$symbol %in% genesO$symbol),]
    genesO$facet = title
    genesO$start <- matrixStats::rowMaxs(cbind(genesO$start, start(region)))
    genesO$end <- matrixStats::rowMins(cbind(genesO$end, end(region)))

    #Collapse Iteratively
    #backwards iteration so that the last value chosen is the lowest cluster possible to fit in.
    genesO$cluster <- 0
    for(i in seq_len(nrow(genesO))){
      if(i==1){
        genesO$cluster[i] <- 1
      }else{
        for(j in seq_len(max(genesO$cluster))){
          jEnd <- rev(genesO$end)[match(rev(seq_len(max(genesO$cluster)))[j], rev(genesO$cluster))]
          if(genesO$start[i] > jEnd + median(genesO$width)){
            genesO$cluster[i] <- rev(genesO$cluster)[match(rev(seq_len(max(genesO$cluster)))[j],rev(genesO$cluster))]
          }
        }
        if(genesO$cluster[i]==0){
          genesO$cluster[i] <- genesO$cluster[i-1] + 1
        }
      }
    }
    exonsO$cluster <- genesO$cluster[match(exonsO$symbol, genesO$symbol)]
    pal <- c("-"=colorMinus,"+"=colorPlus,"*"=colorPlus)
    
    p <- ggplot(data = genesO, aes(color = strand, fill = strand)) +
      facet_grid(facet~.) +
      #################################################
      #Limits
      #################################################
      ylim(c(0.5, max(genesO$cluster) + 0.5)) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) + 
      #################################################
      #Segment for Not Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)!="-"),], 
        aes(x = start, xend = end, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segment for Minus Stranded
      #################################################
      geom_segment(data = genesO[which(as.character(genesO$strand)=="-"),], 
        aes(x = end, xend = start, y = cluster, yend = cluster, color = strand),size=geneWidth) +
      #################################################
      #Segement for Exons
      #################################################
      geom_segment(data = exonsO, aes(x = start, xend = end, y = cluster, 
        yend = cluster, color = strand),size=exonWidth) +
      #################################################
      #Colors
      #################################################
      scale_color_manual(values = pal, guide = FALSE) + 
      scale_fill_manual(values = pal) +
      #################################################
      #Theme
      #################################################
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
      theme(legend.text = element_text(size = baseSize), strip.text.y = element_text(size = facetbaseSize, angle = 0)) +
      guides(fill = guide_legend(override.aes = list(colour = NA, shape = "c", size=3)), color = FALSE) + 
      theme(legend.position="bottom") +
      theme(legend.title=element_text(size=5), legend.text=element_text(size=7),
        legend.key.size = unit(0.75,"line"), legend.background = element_rect(color =NA), strip.background = element_blank())

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand!="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand!="-"),], 
        aes(x = start, y = cluster, label = symbol, color = strand, fill = NA), 
          segment.color = "grey", nudge_x = -0.01*(end(region) - start(region)), nudge_y = -0.25, 
          size = labelSize, direction = "x")
    }

    #Add Labels if There are Genes with this orientation!
    if(length(which(genesO$strand=="-")) > 0){
      p <- p + ggrepel::geom_label_repel(data=genesO[which(genesO$strand=="-"),], 
        aes(x = end, y = cluster, label = symbol, color = strand, fill = NA), 
          segment.color = "grey", nudge_x = +0.01*(end(region) - start(region)), nudge_y = 0.25, 
          size = labelSize, direction = "x")
    }

    p <- p + theme(legend.justification = c(0, 1), 
      legend.background = element_rect(colour = NA, fill = NA), legend.position="none")

  }else{

    #create empty plot
    df <- data.frame(facet = "GeneTrack", start = 0, end = 0, strand = "*", symbol = "none")
    pal <- c("*"=colorPlus)
    p <- ggplot(data = df, aes(start, end, fill = strand)) + geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_color_manual(values = pal) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(!is.ggplot(p)){
    message("INFO : geneTrack is not a ggplot!")
  }

  return(list("geneTracksPlot"=p,"geneTracksDF"=df))

}

#######################################################
# Feature Tracks
#######################################################
.featureTracks <- function(
  features = NULL, 
  region = NULL, 
  title = "FeatureTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = NULL,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  region <- .validGRanges(region)
  region <- .subsetSeqnamesGR(region[1], as.character(seqnames(region[1])))

  if(!is.null(features)){

    if(!ArchR:::.isGRList(features)){
      features <- .validGRanges(features)
      featureList <- SimpleList(FeatureTrack = features)
      hideY <- TRUE
    }else{
      featureList <- features
      hideY <- FALSE
    }
    featureList <- featureList[rev(seq_along(featureList))]

    featureO <- lapply(seq_along(featureList), function(x){
      featurex <- featureList[[x]]
      namex <- names(featureList)[x]
      mcols(featurex) <- NULL
      sub <- subsetByOverlaps(featurex, region, ignore.strand = TRUE)
      if(length(sub) > 0){
        data.frame(sub, name = namex)
      }else{
        empty <- GRanges(as.character(seqnames(region[1])), ranges = IRanges(0,0))
        data.frame(empty, name = namex)
      }

    })

    featureO <- Reduce("rbind", featureO)
    featureO$facet <- title

    if(is.null(pal)){
      pal <- paletteDiscrete(set = "stallion", values = rev(unique(paste0(featureO$name))))
    }
    
    featureO$name <- factor(paste0(featureO$name), levels=names(featureList))

    p <- ggplot(data = featureO, aes(color = name)) +
      facet_grid(facet~.) +
      geom_segment(data = featureO, aes(x = start, xend = end, y = name, yend = name, color = name), size=featureWidth) +
      ylab("") + xlab("") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      scale_color_manual(values = pal) +
      theme(legend.text = element_text(size = baseSize)) + 
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      guides(color = FALSE, fill = FALSE) + theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank())

  }else{

    #create empty plot
    df <- data.frame(facet = "FeatureTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    message("featureTrack is not a ggplot!")
  }

  return(list("featureTracksPlot"=p,"featureTracksDF"=df))

}

#######################################################
# Loop Tracks
#######################################################
.loopTracks <- function(
  loops = NULL, 
  region = NULL, 
  title = "LoopTrack", 
  pal = NULL,
  baseSize = 9, 
  facetbaseSize = 9,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){

  getArchDF <- function(lp, r = 100){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * (rx/max(rx))
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$value)){
      mcols(lp)$value <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), value = mcols(lp)$value[z])
    }) %>% Reduce("rbind",.)
    return(df)
  }

  if(!is.null(loops)){

    if(is(loops, "GRanges")){
      loops <- SimpleList(Loops = loops)
    }else if(ArchR:::.isGRList(loops)){
    }else{
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }

    valueMin <- min(unlist(lapply(loops, function(x) min(x$value))))
    valueMax <- max(unlist(lapply(loops, function(x) max(x$value))))

    loopO <- lapply(seq_along(loops), function(x){
       subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within") 
       if(length(subLoops)>0){
         dfx <- getArchDF(subLoops)
         dfx$name <- Rle(paste0(names(loops)[x]))
         dfx
       }else{
         NULL
       }
    }) %>% Reduce("rbind",.)

    testDim <- tryCatch({
      if(is.null(loopO)){
        FALSE
      }
      if(nrow(loopO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })

    if(testDim){

      loopO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
      }

      p <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = value)) + 
        geom_line() +
        facet_grid(name ~ .) +
        ylab("") + 
        coord_cartesian(ylim = c(-100,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(valueMin, valueMax)) +
        theme(legend.text = element_text(size = baseSize)) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
        guides(color= guide_colorbar(barwidth = 0.75, barheight = 3))

    }else{

      #create empty plot
      df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
      p <- ggplot(data = df, aes(start, end)) + 
        geom_point() +
        facet_grid(facet~.) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

    }

  }else{

    #create empty plot
    df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
	  ArchR:::.logError("loopTracks is not a ggplot!", fn = ".loopTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(list("loopTracksPlot"=p,"loopTracksDF"=df))

}

.subsetSeqnamesGR <- function(gr = NULL, names = NULL){
  gr <- gr[which(as.character(seqnames(gr)) %in% names),]
  seqlevels(gr) <- as.character(unique(seqnames(gr)))
  return(gr)
}

#######################################################
# scATAC Track Methods
#######################################################

.scTracks <- function(
  ArchRProj = NULL,
  region = NULL,
  tileSize = 100,
  minCells = 5,
  maxCells = 100,
  groupBy = "Clusters",
  useGroups = NULL,
  threads = 1,
  baseSize = 7,
  scTileSize = 0.5,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  pal = NULL,
  tstart = NULL,
  verbose = FALSE,
  logFile = NULL
  ){


  if(is.null(tstart)){
    tstart <- Sys.time()
  }

  #Group Info
  cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  if(!is.null(minCells)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% names(table(cellGroups)[table(cellGroups) >= minCells]),,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  if(!is.null(useGroups)){
    ArchRProj@cellColData <- ArchRProj@cellColData[cellGroups %bcin% useGroups,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
  }
  tabGroups <- table(cellGroups)
  
  if(any(tabGroups > maxCells)){
    cellGroups2 <- getCellColData(ArchRProj, groupBy, drop = FALSE)
    splitGroups <- split(rownames(cellGroups2), cellGroups2[,1])
    useCells <- lapply(seq_along(splitGroups), function(x){
      if(length(splitGroups[[x]]) > maxCells){
        sample(splitGroups[[x]], maxCells)
      }else{
        splitGroups[[x]]
      }
    }) %>% unlist
    ArchRProj@cellColData <- ArchRProj@cellColData[useCells,,drop=FALSE]
    cellGroups <- getCellColData(ArchRProj, groupBy, drop = TRUE)
    tabGroups <- table(cellGroups)
  }

  cellsBySample <- split(rownames(getCellColData(ArchRProj)), getCellColData(ArchRProj, "Sample", drop = TRUE))
  groupsBySample <- split(cellGroups, getCellColData(ArchRProj, "Sample", drop = TRUE))
  uniqueGroups <- gtools::mixedsort(unique(cellGroups))
  
  #Tile Region
  regionTiles <- seq(trunc(start(region) / tileSize), trunc(end(region) / tileSize) + 1) * tileSize
  ArrowFiles <- getArrowFiles(ArchRProj)
  ArrowFiles <- ArrowFiles[names(cellsBySample)]

  groupMat <- .safelapply(seq_along(ArrowFiles), function(i){
    tryCatch({
      .regionSCArrows(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    }, error = function(e){
      errorList <- list(
        ArrowFile = ArrowFiles[i], 
        region = region, 
        regionTiles = regionTiles,
        tileSize = tileSize,
        cellNames = cellsBySample[[names(ArrowFiles)[i]]],
        cellGroups = groupsBySample[[names(ArrowFiles)[i]]],
        uniqueGroups = uniqueGroups
      )
    })
  }, threads = threads) %>% Reduce("cbind" , .)

  groupDF <- data.frame(Matrix::summary(groupMat))
  groupDF$group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[colnames(groupMat)[groupDF$j], 1]
  groupDF <- lapply(split(groupDF, groupDF$group), function(z){
    nz <- tabGroups[z$group[1]]
    nc <- length(unique(z$j))
    idx <- sort(sample(seq_len(nz), nc))
    idx[1] <- 1
    idx[length(idx)] <- nz
    z$y <- idx[match(z$j, unique(z$j))]
    z
  }) %>% Reduce("rbind", .)
  groupDF$bp <- regionTiles[groupDF$i]
  
  if(is.null(pal)){
    pal <- suppressWarnings(paletteDiscrete(values = names(tabGroups)))
  }

  nn <- paste0(names(tabGroups), ":", tabGroups)
  names(nn) <- names(tabGroups)
  groupDF$group2 <- nn[groupDF$group]
  names(pal) <- nn[names(pal)]

  title <- paste0(as.character(seqnames(region)),":", start(region)-1, "-", end(region), " ", title)
  
  #Re-Order
  groupDF$group2 <- factor(
    paste0(groupDF$group2), 
    levels = gtools::mixedsort(unique(paste0(groupDF$group2)))
  )

  p <- ggplot(groupDF, aes(x=bp, y=y, width = tileSize, fill = group2, color = group2)) + 
      geom_tile(size = scTileSize) + 
      facet_grid(group2 ~ ., scales="free_y") + 
      theme_ArchR() + 
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      ylab("Binarized SC Coverage") + 
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme_ArchR(baseSize = baseSize,
                  baseRectSize = borderWidth,
                  baseLineSize = tickWidth,
                  legendPosition = "right",
                  axisTickCm = 0.1) +
      theme(panel.spacing= unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.text = element_text(
              size = facetbaseSize, 
              color = "black", 
              margin = margin(0,0.35,0,0.35, "cm")),
              strip.text.y = element_text(angle = 0),
            strip.background = element_rect(color="black")) +
      guides(fill = FALSE, colour = FALSE) + ggtitle(title)

   return(list("scTracksPlot"=p,"scTracksDF"=df))

}

.safelapply<-function (..., threads = 1, preschedule = FALSE)
{
    if (tolower(.Platform$OS.type) == "windows") {
        threads <- 1
    }
    if (threads > 1) {
        o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)
        errorMsg <- list()
        for (i in seq_along(o)) {
            if (inherits(o[[i]], "try-error")) {
                capOut <- utils::capture.output(o[[i]])
                capOut <- capOut[!grepl("attr\\(\\,|try-error",
                  capOut)]
                capOut <- head(capOut, 10)
                capOut <- unlist(lapply(capOut, function(x) substr(x,
                  1, 250)))
                capOut <- paste0("\t", capOut)
                errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ",
                  i, " : "), capOut), "\n")
            }
        }
        if (length(errorMsg) != 0) {
            errorMsg <- unlist(errorMsg)
            errorMsg <- head(errorMsg, 50)
            errorMsg[1] <- paste0("\n", errorMsg[1])
            stop(errorMsg)
        }
    }
    else {
        o <- lapply(...)
    }
    o
}

.validGRanges<-function (gr = NULL)
{
    stopifnot(!is.null(gr))
    if (inherits(gr, "GRanges")) {
        return(gr)
    }
    else {
        stop("Error cannot validate genomic range!")
    }
}

.regionSCArrows <- function(
  ArrowFile = NULL,
  region = NULL,
  regionTiles = NULL,
  tileSize = NULL,
  cellNames = NULL,
  cellGroups = NULL,
  uniqueGroups = NULL,
  logFile = NULL
  ){

  cellFragsRegion <- ArchR:::.getFragsFromArrow(
      ArrowFile = ArrowFile,
      chr = paste0(seqnames(region)),
      cellNames = cellNames,
      out = "GRanges"
    ) %>% subsetByOverlaps(., region, ignore.strand = TRUE)

  #Starts
  ts <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ids <- which(ts > 0)

  #Ends
  te <- match(trunc(start(cellFragsRegion)/tileSize) * tileSize, regionTiles, nomatch = 0)
  ide <- which(te > 0)

  #Match
  matchID <- S4Vectors::match(mcols(cellFragsRegion)$RG, cellNames)

  #Sparse Matrix
  mat <- Matrix::sparseMatrix(
    i = c(ts[ids], te[ide]),
    j = as.vector(c(matchID[ids], matchID[ide])),
    x = rep(1,  length(ids) + length(ide)),
    dims = c(length(regionTiles), length(cellNames))
  )
  colnames(mat) <- cellNames

  mat@x[mat@x > 1] <- 1

  return(mat)

}

.combinedFeaturePlot <- function(
  plotList = NULL,
  useMatrix = NULL,
  featureMat = NULL,
  log2Norm = FALSE,
  feature = NULL,
  pal = NULL,
  sizes = NULL,
  baseSize = NULL,
  facetbaseSize = NULL,
  borderWidth = NULL,
  tickWidth = NULL
  ){


  if(is.null(pal)){
    pal <- paletteDiscrete(values=featureMat$Group, set = "stallion")
  }

  if(log2Norm){
    title <- paste0("Log2 ", useMatrix, " : ", feature)
  }else{
    title <- paste0("Raw ", useMatrix, " : ", feature)
  }

  featurePlot <- ggGroup(
      x = featureMat$Group,
      y = featureMat[,feature],
      groupOrder = gtools::mixedsort(paste0(unique(featureMat$Group))),
      pal = pal
    ) +
    facet_wrap(x~., ncol=1,scales="free_y",strip.position="right") +
    guides(fill = FALSE, colour = FALSE) +
    theme_ArchR(baseSize = baseSize,
              baseRectSize = borderWidth,
              baseLineSize = tickWidth,
              legendPosition = "right",
              axisTickCm = 0.1) +
    theme(panel.spacing= unit(0, "lines"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text = element_text(
          size = facetbaseSize,
          color = "black",
          margin = margin(0,0.35,0,0.35, "cm")),
          strip.text.y = element_text(angle = 0),
        strip.background = element_rect(color="black")) +
    theme(plot.margin = unit(c(0.35, 0.15, 0.35, 0.15), "cm")) +
    ggtitle(title)

  if(any(tolower(names(plotList)) %in% "bulktrack")){

    idx <- which(tolower(names(plotList)) == "bulktrack")

    p <- plotList[[idx]] + featurePlot + plot_spacer()

    plotList[idx] <- NULL

    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] + plot_spacer() + plot_spacer()
    }

    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2),
      heights = sizes
    )

  }else{


    idx <- which(tolower(names(plotList)) == "sctrack")

    p <- plotList[[idx]] + featurePlot + plot_spacer()

    plotList[idx] <- NULL

    for(i in seq_along(plotList)){
      p <- p + plotList[[i]] + plot_spacer() + plot_spacer()
    }

    p <- p + plot_layout(
      ncol = 3,
      widths = c(3, 1, 0.2),
      heights = sizes
    )

  }

  p

}
