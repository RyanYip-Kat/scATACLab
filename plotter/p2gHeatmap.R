library(ArchR)
p2gheatmap=function(p2gMatrices,
		    limitsATAC = c(-2, 2),
		    limitsRNA = c(-2, 2),
		    palGroup = NULL,
		    palATAC = paletteContinuous("solarExtra"),
                    palRNA = paletteContinuous("blueYellow")){

       ATAC=p2gMatrices[["ATAC"]]
       RNA=p2gMatrices[["RNA"]]
       p2g=p2gMatrices[["Peak2GeneLinks"]]

       mATAC=ATAC[["matrix"]]
       colATAC=ATAC[["colData"]]
       mRNA=RNA[["matrix"]]
       colRNA=RNA[["colData"]]

       KDF_3=RNA[["kmeansId"]]

       htATAC <- ArchR:::.ArchRHeatmap(
				mat = mATAC,
                                scale = FALSE,
                                limits = limitsATAC,
                                color = palATAC,
                                colData = colATAC,
                                clusterCols = FALSE,
                                clusterRows = FALSE,
                                split = KDF_3,
                                labelRows = FALSE,
                                labelCols = FALSE,
                                draw = FALSE,
                                name = paste0("ATAC Z-Scores\n", nrow(mATAC), " P2GLinks"))
	
	htRNA <- ArchR:::.ArchRHeatmap(
			       mat = mRNA,
                               scale = FALSE,
                               limits = limitsRNA,
                               color = palRNA,
			       colData = colRNA,
                               clusterCols = FALSE,
                               clusterRows = FALSE,
                               split = KDF_3,
                               labelRows = FALSE,
                               labelCols = FALSE,
                               draw = FALSE,
			       name = paste0("RNA Z-Scores\n", nrow(mRNA), " P2GLinks")
			       )
        return(htATAC + htRNA)
}
