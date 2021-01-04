library(data.table)
library(Matrix)
library(GenomicRanges)
library(magrittr)
library(SummarizedExperiment)
library(yaml)
library(Rcpp)
#Helper function for summing sparse matrix groups
binarizeMat <- function(mat){
    mat@x[mat@x > 0] <- 1
    mat
}

keepFilteredChromosomes <- function(x, remove = c("chrM"), underscore = TRUE, standard = TRUE, pruning.mode="coarse"){
        #first we remove all non standard chromosomes
        if(standard){
                x <- GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = pruning.mode)
        }
        #Then check for underscores or specified remove
        seqNames <- seqlevels(x)
        chrRemove <- c()
        #first we remove all chr with an underscore
        if(underscore){
                chrRemove <- c(chrRemove, which(grepl("_", seqNames)))
        }
        #next we remove all chr specified in remove
        chrRemove <- c(chrRemove, which(seqNames %in% remove))
        if(length(chrRemove) > 0){
                chrKeep <- seqNames[-chrRemove]
        }else{
                chrKeep <- seqNames
        }
        #this function restores seqlevels
        seqlevels(x, pruning.mode=pruning.mode) <- chrKeep
        return(x)
}


grToFeature <- function(gr){
    peakinfo <- data.frame(
        row.names = paste(seqnames(gr),start(gr),end(gr),sep="_"),
        site_name = paste(seqnames(gr),start(gr),end(gr),sep="_"),
        chr = gsub("chr","",as.character(seqnames(gr))),
        bp1 = start(gr),
        bp2 = end(gr)
    )
    return(peakinfo)
}

#############################
featureToGR <- function(feature){
    featureSplit <- stringr::str_split(paste0(feature), pattern = "_", n = 3, simplify = TRUE)
    gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
    return(gr)
}

getGeneGTF <- function(file){
    #Import
    message("Reading in GTF...")
    importGTF <- rtracklayer::import(file)
    #Exon Info
    message("Computing Effective Exon Lengths...")
    exonGTF <- importGTF[importGTF$type=="exon",]
    exonList <- reduce(split(exonGTF, mcols(exonGTF)$gene_id))
    exonReduced <- unlist(exonList, use.names=TRUE)
    mcols(exonReduced)$gene_id <- names(exonReduced)
    mcols(exonReduced)$widths <- width(exonReduced)
    exonSplit <- split(exonReduced$widths, mcols(exonReduced)$gene_id)
    exonLengths <- lapply(seq_along(exonSplit), function(x) sum(exonSplit[[x]])) %>%
        unlist %>% data.frame(row.names=names(exonSplit), effLength=.)
    #Gene Info
    message("Constructing gene GTF...")
    geneGTF1 <- importGTF[importGTF$type=="gene",]
    geneGTF2 <- GRanges(
            seqnames=paste0("chr",seqnames(geneGTF1)),
            ranges=ranges(geneGTF1),
            strand=strand(geneGTF1),
            gene_name=geneGTF1$gene_name,
            gene_id=geneGTF1$gene_id
        ) %>% keepFilteredChromosomes %>% sortSeqlevels %>% sort(.,ignore.strand=TRUE)
    mcols(geneGTF2)$exonLength <- exonLengths[geneGTF2$gene_id,]
    return(geneGTF2)
}
# gtfFile <- "data/genes.gtf"
#genes <- getGeneGTF(gtfFile) %>% resize(1,"start") %>% resize(tssWindow * 2 + 1, "center")
