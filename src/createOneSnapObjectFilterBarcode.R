library(SnapATAC)
library(GenomicRanges)
library(viridisLite)
library(ggplot2)
library(argparse)
library(stringr)
############################
mkdir=function(path){
	if(!dir.exists(path)){
		dir.create(path,recursive=TRUE)
	}
}

createSnapObject<-function(snap_file,
			   sample_name="pbmc",
			   barcode_file=NULL,
			   umi=c(3.2,5),
			   promoter_ratio=c(0.15,0.5),
			   outDir="./Results"){

	umi_l=umi[1]
	umi_u=umi[2]

	promoter_l=promoter_ratio[1]
	promoter_u=promoter_ratio[2]

        mkdir(outDir)	
	message("INFO : Create snapObject ...")
	x.sp = createSnap(file=snap_file, sample=sample_name,num.cores=16)
	if(!is.null(barcode_file)){
		barcodes = read.csv(barcode_file,head=TRUE,stringsAsFactors=FALSE)
		barcodes = barcodes[2:nrow(barcodes),];
                promoter_ratio = (barcodes$promoter_region_fragments+1) / (barcodes$passed_filters + 1);
                UMI = log(barcodes$passed_filters+1, 10);
                data = data.frame(UMI=UMI, promoter_ratio=promoter_ratio);
                barcodes$promoter_ratio = promoter_ratio

		message("INFO : Plot metric ...")
		metricDir=file.path(outDir,"metric")
                mkdir(metricDir)

		p1 = ggplot(data,
			    aes(x= UMI, y= promoter_ratio)) +
                            geom_point(size=0.1, col="grey") +
                            theme_classic() +
                            ggtitle(sample_name) +
                            ylim(0, 1) + xlim(0, 6) +
                            labs(x = "log10(UMI)", y="promoter ratio")
                ggsave(file.path(metricDir,paste0(sample_name,".pdf")),plot=p1,width=12,height=12)

                barcodes.sel = barcodes[which(UMI >= umi_l & UMI <= umi_u & promoter_ratio >= promoter_l & promoter_ratio <= promoter_u),];
                barcodes.sel$barcode=paste(str_to_upper(sample_name),as.character(barcodes.sel$barcode),sep="_")
                rownames(barcodes.sel) = barcodes.sel$barcode;
		x.sp = x.sp[which(x.sp@barcode%in%barcodes.sel$barcode),];
                x.sp@metaData = barcodes.sel[x.sp@barcode,]
	}
	message("INFO : addBmat ...")
	x.sp = addBmatToSnap(x.sp, bin.size=5000)
	x.sp@metaData$sample=x.sp@sample
	message("INFO : Save ...")
	saveRDS(x.sp,file.path(outDir,"snapObject.rds"))
	message("INFO : Done!")
}


############################
parser <- ArgumentParser(description='Create snap Object from mutiple sample ...')
parser$add_argument("--snap",
                    type="character",
                    default=NULL,
                    help="path of *.snap file")

parser$add_argument("--name",
		    type="character",
		    default=NULL,
		    help="sample name for snap file")

parser$add_argument("--barcode",
		    type="character",
                    default=NULL,
		    help="10x barcode,*singlecell.csv file")

parser$add_argument("--umi",
		    type="double",
		    nargs="+",
		    default=NULL,
		    help="umi interval,length must be == 2")

parser$add_argument("--promoter",
                    type="double",
                    nargs="+",
                    default=NULL,
                    help="promoter_ratio interval,length must be == 2")


parser$add_argument("--outdir",
                    type="character",
                    default="Snap_result")

args <- parser$parse_args()

#############################
outDir = args$outdir
metricDir=file.path(outDir,"metric")
mkdir(outDir)
mkdir(metricDir)


#############################
message("INFO : loading dataset ...")
sample_name=as.character(args$name)
snap_file=as.character(args$snap)
barcode_file=as.character(args$barcode)
umi_interval=args$umi
promoter_interval=args$promoter


stopifnot(length(umi_interval)==2) # c(umi_lower,umi_upper)
stopifnot(length(promoter_interval)==2)  # c(promoter_lower,promoter_upper)

message("INFO : create snp ...")
createSnapObject(snap_file=snap_file,
		 sample_name=sample_name,
		 barcode_file=barcode_file,
		 umi=umi_interval,
		 promoter=promoter_interval,
		 outDir=outDir)

message("INFO : Done!")
