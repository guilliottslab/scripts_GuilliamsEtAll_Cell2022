##### Function get colors of ggplot
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

if(exists("listLabels")){
  listColors<-ggplotColours(length(listLabels))
}else{
  listColors<-ggplotColours(6)
}

##### Function drawVlnPlot
drawVlnPlot<-function(toPlot, fileName, colsToColor){
  toPlot<-toPlot[order(toPlot[,colsToColor[1]]),]
  p_nGene <- ggplot(toPlot, aes(staticNr, nGene)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()
  
  toPlot<-toPlot[order(toPlot[,colsToColor[2]]),]
  p_nUMI <- ggplot(toPlot, aes(staticNr, nUMI)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()
  
  toPlot<-toPlot[order(toPlot[,colsToColor[3]]),]
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mito)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic()
  
  grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
  if(fileName != ""){
    ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file=fileName, dpi=200)
  }
}

##### Function drawVlnPlot_split
drawVlnPlot_split<-function(toPlot, fileName, colsToColor){
  p_nGene <- ggplot(toPlot, aes(orig.ident, nGene, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_nUMI <- ggplot(toPlot, aes(orig.ident, nUMI, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  p_mito <- ggplot(toPlot, aes(orig.ident, percent.mito, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
  if(fileName != ""){
    ggsave(grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3),file=fileName, dpi=200, width=25, height=12)
  }
}

##### Function drawVlnPlotSeurat_split
drawVlnPlotSeurat_split<-function(metaDataTable, fileName){
  p1<-ggplot(metaDataTable, aes(orig.ident, nGene, fill=orig.ident)) +
    geom_jitter(height = 0, width = 0.3, alpha=0.4) + 
    geom_violin(alpha=0.6) + 
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('nGene')
  p2<-ggplot(metaDataTable, aes(orig.ident, nUMI, fill=orig.ident)) +
    geom_jitter(height = 0, width = 0.3, alpha=0.4) + 
    geom_violin(alpha=0.6) + 
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('nUMI')
  p3<-ggplot(metaDataTable, aes(orig.ident, percent.mito, fill=orig.ident)) +
    geom_jitter(height = 0, width = 0.3, alpha=0.4) + 
    geom_violin(alpha=0.6) + 
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('percent.mito')
  
  grid.arrange(p1,p2,p3, ncol=3)
  if(fileName != ""){
    ggsave(grid.arrange(p1,p2,p3, ncol=3),file=fileName, dpi=200)
  }
}

##### Function drawPCAplotSplit
drawPCAplotSplit<-function(pcaTable, fileName){
  pdf(file=fileName, width = 10)
  p <- ggplot()+
    geom_point(aes(x=PC1,y=PC2, colour=orig.ident), data=pcaTable, size=2, shape=20) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  print(p)
  
  p <- ggplot()+
    geom_point(aes(x=PC1,y=PC3, colour=orig.ident), data=pcaTable, size=2, shape=20) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  print(p)
  dev.off()
}

##### Function drawFeaturePlot
drawFeaturePlot<-function(group, geneSymbol, reductionType){
  p <- ggplot()+
    geom_point(aes(x=tSNE_1,y=tSNE_2, colour=logTable[geneSymbol,]), data=tsneTable, size=1, shape=20)
  
  if(reductionType=="umap"){
    p <- ggplot()+
      geom_point(aes(x=UMAP_1,y=UMAP_2, colour=logTable[geneSymbol,]), data=umapTable, size=1, shape=20)
  }
  
  p <- p +
    scale_colour_gradient(low = "gray", high="blue") +
    ggtitle(paste0(group,": ",geneSymbol)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face="bold"), legend.position = "none")
  return(p)
}

##### Function drawMultipleFeaturePlot
drawMultipleFeaturePlot<-function(groupName, markers, ncol, reductionType){
  myList=list()
  for(i in 1:length(markers)){
    p<-drawFeaturePlot(groupName,markers[i], reductionType)
    myList[[i]]<-p
  }
  return(grid.arrange(grobs = myList, ncol = ncol, nrow = 2))
}

##### Function drawUMI_mitoPlot
drawUMI_mitoPlot<-function(coordsTable, reductionType, clusterMatrix, columnName, titleInfo){
  
  columnNr<-which(colnames(clusterMatrix)==columnName)

  p <- ggplot()+
    geom_point(aes(x=tSNE_1,y=tSNE_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  
  if(reductionType=="umap"){
    p <- ggplot()+
      geom_point(aes(x=UMAP_1,y=UMAP_2, colour=clusterMatrix[,columnNr]), data=coordsTable, size=2, shape=20) 
  }
  
  p<-p +
    scale_colour_gradientn(colours = c("darkblue","cyan","green","yellow","orange","darkred")) +
    ggtitle(paste0(titleInfo," (",reductionType,")")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  
  return(p)
}


##### Function drawTSNEplot_colorSample
drawTSNEplot_colorSample<-function(highlight, clusterMatrix, columnName, reductionType){
  
  columnNr<-which(colnames(clusterMatrix)==columnName)
  
  p<-ggplot()+
    geom_point(aes(x=tSNE_1,y=tSNE_2, colour=clusterMatrix[,columnNr], alpha=clusterMatrix[,columnNr]), data=tsneTable[rownames(clusterMatrix),], 
               size=1)
  
  if(reductionType=="umap"){
    p<-ggplot()+
      geom_point(aes(x=UMAP_1,y=UMAP_2, colour=clusterMatrix[,columnNr], alpha=clusterMatrix[,columnNr]), data=umapTable[rownames(clusterMatrix),], 
                 size=1)
  }

  p<-p +
    scale_color_manual(values=highlight$color) +
    scale_alpha_manual(values=highlight$alpha) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.text=element_text(size=12))
  
  return(p)

}


##### Function getTheColGaps
getTheColGaps<-function(clusterMatrix, resColumn, clusters){
  columnID<-which(colnames(clusterMatrix)==resColumn)
  
  theGaps<-nrow(clusterMatrix[clusterMatrix[,columnID] == clusters[1],])
  for(i in clusters[-1]){
    lengthCluster<-nrow(clusterMatrix[clusterMatrix[,columnID] == i,])
    toAdd<-theGaps[length(theGaps)]+lengthCluster
    theGaps<-c(theGaps, toAdd)
  }
  return(theGaps)
}


##### Function colorSomeCells
colorSomeCells<-function(clusterMatrix, coordsTable, cellsToColor, title=NULL){
  clusterMatrix$toColor=FALSE
  clusterMatrix[cellsToColor,which(colnames(clusterMatrix)=="toColor")]<-TRUE
  clusterMatrix<-clusterMatrix[order(clusterMatrix$toColor),]
  coordsTable<-coordsTable[rownames(clusterMatrix),]
  
  p1 <- ggplot()+
    geom_point(aes_string(x=colnames(coordsTable)[1],y=colnames(coordsTable)[2], colour=clusterMatrix$toColor), data=coordsTable, size=1) +
    scale_color_manual(values=c('TRUE'="orangered",'FALSE'="lightgray")) +
    theme_classic() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.text=element_text(size=12))
  
  return(p1)
}

##### Function getRandomCells
getRandomCells<-function(maxNrClusters, maxNrCells){
  listSelectedCells<-list()
  for(i in 0:maxNrClusters){
    cells<-WhichCells(seuratObj, idents = i)
    IDs<-sample(1:length(cells), maxNrCells, replace=F)
    selectedCells<-cells[IDs]
    
    j<-i+1
    listSelectedCells[[j]]<-selectedCells
  }
  names(listSelectedCells)<-paste0('cl',0:maxNrClusters)
  return(listSelectedCells)
}


##### Function printDiagnostics
printDiagnostics<-function(sampleName, fileName, diagnostics){
  cat("------------------------------------------------------------",file=fileName,sep="\n")
  cat(paste0("\t",sampleName),file=fileName,sep="\n", append=TRUE)
  cat("------------------------------------------------------------",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("The raw dataset contains ",diagnostics[['nrGenes']]," genes and ",diagnostics[['nrCells']]," cells."),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("Minimum number of genes per cell detected: ",diagnostics[['minGenesPerCell']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Maximum number of genes per cell detected: ",diagnostics[['maxGenesPerCell']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Mean number of genes per cell detected: ",diagnostics[['meanGenesPerCell']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Median number of genes per cell detected: ",diagnostics[['medianGenesPerCell']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Cells with < 200 genes: ",diagnostics[['cellsLess200genes']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("Number of genes expressed in none of the cells: ",diagnostics[['genesNotExpressed']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Number of genes expressed in < 3 cells: ",diagnostics[['genesLess3cells']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("% zero inflated: ",diagnostics[['zeroInflation']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("Cells per sample: ",diagnostics[['splitSamples']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)

  cat("--------------- Filtering ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("nmad for nGene low: ",diagnostics[['nmad.low.feature']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("nmad for nGene high: ",diagnostics[['nmad.high.feature']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for nGene low: ",diagnostics[['feature.drop.low']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for nGene high: ",diagnostics[['feature.drop.high']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for nGene: ",diagnostics[['feature.drop']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="", append=TRUE)

  cat(paste0("nmad for nUMI low: ",diagnostics[['nmad.low.libsize']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("nmad for nUMI high: ",diagnostics[['nmad.high.libsize']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for nUMI low: ",diagnostics[['libsize.drop.low']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for nUMI high: ",diagnostics[['libsize.drop.high']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for nUMI: ",diagnostics[['libsize.drop']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="", append=TRUE)

  cat(paste0("nmad for %mito low: ",diagnostics[['nmad.low.mito']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("nmad for %mito high: ",diagnostics[['nmad.high.mito']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for %mito low: ",diagnostics[['mito.drop.low']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for %mito high: ",diagnostics[['mito.drop.high']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Outlier cells for %mito: ",diagnostics[['mito.drop']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("===> Removed outlier cells: ",diagnostics[['firstRemove']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("Outlier cells based on PCA: ",diagnostics[['pcaRemove']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("===> Total removed cells: ",diagnostics[['totalRemove']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat(paste0("Cells per sample after filtering: ",diagnostics[['splitSamplesAfterFiltering']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat("--------------- Dimensions ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim raw data:\t\t\t",diagnostics[['dimRawData']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim sce object:\t\t\t",diagnostics[['dimSce']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim after 1st filtering:\t",diagnostics[['dimFirstRemove']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim after PCA:\t\t\t",diagnostics[['dimAfterPCA']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim before seuratObj:\t\t",diagnostics[['dimBeforeSeuratObj']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim after 2nd filtering:\t",diagnostics[['dimAfterSeuratObj']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat("--------------- Analysis ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("Highly variable genes: ",diagnostics[['varGenes']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("HVG: x low cutoff: ",diagnostics[['varGenes_xLow']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("HVG: x high cutoff: ",diagnostics[['varGenes_xHigh']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("HVG: y cutoff: ",diagnostics[['varGenes_y']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)

  cat(paste0("PCs: ",diagnostics[['dimsPC']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Resolution: ",diagnostics[['res']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat("--------------- DE genes ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("Markers per cluster: ",diagnostics[['markersPerCluster']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
}


##### Function printDiagnosticsHarmony
printDiagnosticsHarmony<-function(sampleName, fileName, diagnostics){
  cat("------------------------------------------------------------",file=fileName,sep="\n")
  cat(paste0("\tHarmony: ",sampleName),file=fileName,sep="\n", append=TRUE)
  cat("------------------------------------------------------------",file=fileName,sep="\n", append=TRUE)
  
  cat("--------------- Dimensions ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim before seuratObj:\t\t",diagnostics[['dimBeforeSeuratObj']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Dim after seuratObj:\t",diagnostics[['dimAfterSeuratObj']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat("--------------- Filtering ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("Cells per sample after filtering: ",diagnostics[['splitSamplesAfterFiltering']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
  
  cat("--------------- Analysis ---------------",file=fileName,sep="\n", append=TRUE)
  cat(paste0("Highly variable genes: ",diagnostics[['varGenes']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("PCs: ",diagnostics[['dimsPC']]),file=fileName,sep="\n", append=TRUE)
  cat(paste0("Resolution: ",diagnostics[['res']]),file=fileName,sep="\n", append=TRUE)
  cat("\n",file=fileName,sep="\n", append=TRUE)
}


