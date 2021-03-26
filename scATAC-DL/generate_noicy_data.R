library(ArchR)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(chromVARmotifs)
library(SummarizedExperiment)
library(argparse)
parser <- ArgumentParser(description='ChromVAR handle motif from ArchR')
parser$add_argument("--project",
                    type="character",
                    default=NULL,
                    help="ArchR Project Path")


parser$add_argument("--column",
                    type="character",
                    default="Clusters",
                    help="the column in cellcol to be exported")

parser$add_argument("--subset",
                    nargs="+",
                    type="character",
                    default=NULL,
                    help="the subset of column to be extracted")

parser$add_argument("--ratio",
		    type="double",
		    default=0.1,
		    help="the ratio use for subsample")

parser$add_argument("--groupBy",
		    type="character",
		    default=NULL,
		    help="goupby use for exporting BW file in function getGroupBW")

parser$add_argument("--outdir",
                    type="character",
                    default="output")
args <- parser$parse_args()

############################### funciton
makedir<-function(path){
        if(!dir.exists(path)){
                dir.create(path,recursive=TRUE)
        }
}

outDir=args$outdir
column=args$column
colSub=args$subset
ratio=args$ratio
groupBy=ifelse(is.null(args$groupBy),column,args$groupBy)
makedir(outDir)

message("INFO : Loading dataset")
proj=loadArchRProject(args$project)

#################  wether subset
if(!is.null(column) & !is.null(colSub)){
	metadata=as.data.frame(getCellColData(proj))
	cells=getCellNames(proj)
	message("INFO : Subset ")
	stopifnot(column%in%colnames(metadata))
	target=metadata[[column]]
	idxPass= which(target%in%colSub)
	cellPass=cells[idxPass]
        proj=subsetCells(proj,cellNames=cellPass)

}

#################  subsample
metadata=as.data.frame(getCellColData(proj))
cells=getCellNames(proj)
Groups=unique(metadata[[groupBy]])

message("INFO : Subsample")
cellSample=c()
message("---------")
for(g in Groups){
	cat(sprintf("subsample group--- [ %s ]\n",g))
	idx=which(metadata[[groupBy]]==g)
	cell=rownames(metadata[idx,])
	cell=sample(cell,ceiling(length(cell)*ratio)) # subsample
	cat(sprintf("subsample size --- [ %f ]\n",length(cell)))
	message("--------")
	cellSample=c(cellSample,cell)
}

proj=subsetCells(proj,cellNames=cellSample)

#################
message("INFO : Call Peak")
proj <- addGroupCoverages(ArchRProj = proj, groupBy = groupBy,force = TRUE)
proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    additionalParams = "--bdg --nomodel --nolambda",
    extendSummits=250,     # control peaks width
    groupBy = groupBy,
    pathToMacs2 = findMacs2(),
    force=TRUE
)

message("INFO : add PeaksMatrix")
proj <- addPeakMatrix(proj)

message("INFO : Save!")
saveArchRProject(ArchRProj = proj,
		 outputDirectory = file.path(outDir,"Noicy-PeakSetMacs"),
		 load = FALSE)
