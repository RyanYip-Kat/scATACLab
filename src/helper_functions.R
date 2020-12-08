library(data.table)
library(stringr)
bgzip="/home/ye/anaconda3/envs/BulkBio/bin/bgzip"
tabix="/home/ye/anaconda3/envs/BulkBio/bin/tabix"
ZipAndBuildTbi=function(gzipFile){
        fileName=file.path(dirname(gzipFile),paste0("sorted-",basename(gzipFile)))
        cmd=paste0("sort -k1,1d -k2,2n -k3,3n"," ",gzipFile," > ",fileName)
        system(cmd)
        cmd=paste0(bgzip," ",fileName)
        message("INFO : Zip File ...")
        system(cmd)

        cmd=paste0(tabix," ","--preset=bed"," ",fileName,".gz")
        message("INFO : Build tabix Index  ...")
        system(cmd)
        system(paste0("rm ",gzipFile))
        return(paste0(fileName,".gz"))
}

ReadFragments<- function(frag_gz_file, prefix,sep="_"){
  name="fragments.tsv"
  filename=file.path(dirname(frag_gz_file),paste0("reFormat-",name))
  if(!file.exists(filename)){
	  cat(sprintf("INFO : Loading fragment.tsv.gz from [ %s ] \n",frag_gz_file))
          fragments<- data.table::fread(cmd = paste0("zcat < ", frag_gz_file))
          colnames(fragments)=c("chrom","chromStart","chromEnd","barcode","duplicateCount")
          barcode=paste(prefix,fragments[["barcode"]],sep=sep)
          fragments[["barcode"]]=barcode

	  cat(sprintf("INFO : Writing reFromat fragment.tsv.gz from [ %s ] \n",filename))
	  write.table(fragments,filename,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  }

  message("INFO : Zip and Build Index ...")
  z=ZipAndBuildTbi(filename)
  return(z)
}

