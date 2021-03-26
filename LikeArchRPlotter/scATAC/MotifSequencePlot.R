library(argparse)
library(stringr)
library(ggplot2)
library(Signac)

parser <- ArgumentParser(description='Motif Plot in sequence weight')
parser$add_argument("--object",
                    type="character",
                    default=NULL,
                    help="Object from ArchRCreateMotif.R ")



parser$add_argument("--motifs",
		    nargs="+",
                    type="character",
                    default=NULL,
                    help="motifs to be plot in MotifPlot")

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

################################
outDir=file.path(args$outdir,"Motifs")
makedir(outDir)

################################
message("INFO : Loading dataset ...")
motif.matrix=readRDS(args$object)
motifs=str_to_lower(args$motifs)
#motif.names=colnames(motif.matrix)
motif.names=GetMotifData(object=motif.matrix,slot = "motif.names")
all.motifs=as.character(unlist(motif.names))
input_motif_names=names(motif.names)[str_to_lower(all.motifs)%in%motifs]

message("INFO : Plot ...")
for(i in seq_along(input_motif_names)){
	mo=input_motif_names[i]
	name=as.character(motif.names[[mo]])
	cat(sprintf("INFO : %d [ %s ] of %d\n",i,name,length(input_motif_names)))
	p=MotifPlot(motif.matrix,mo)+theme(plot.title=element_text(size=16))
	ggsave(file.path(outDir,paste0(name,"-motif.pdf")),plot=p,width=10,height=8)
}
message("INFO : Done!")



