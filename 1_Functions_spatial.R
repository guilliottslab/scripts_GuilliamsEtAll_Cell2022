##### Function colorSomeCells
colorSomeCells<-function(clusterMatrix, coordsTable, cellsToColor, title=NULL){
  clusterMatrix$toColor=FALSE
  clusterMatrix[cellsToColor,which(colnames(clusterMatrix)=="toColor")]<-TRUE
  clusterMatrix<-clusterMatrix[order(clusterMatrix$toColor),]
  coordsTable<-coordsTable[rownames(clusterMatrix),]
  
  p1 <- ggplot()+
    geom_point(aes_string(x=colnames(coordsTable)[1],y=colnames(coordsTable)[2], colour=clusterMatrix$toColor), data=coordsTable, size=0.2) +
    scale_color_manual(values=c('TRUE'="orangered",'FALSE'="lightgray")) +
    theme_classic() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right", legend.text=element_text(size=12))
  
  return(p1)
}

##### Function drawVlnPlot
drawVlnPlot<-function(toPlot,  colsToColor){
  toPlot$staticNr=1
  
  toPlot<-toPlot[order(toPlot[,colsToColor[1]]),]
  p_nGene <- ggplot(toPlot, aes(staticNr, nFeature_Spatial)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5, size=0.8) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank())
  
  toPlot<-toPlot[order(toPlot[,colsToColor[2]]),]
  p_nUMI <- ggplot(toPlot, aes(staticNr, nCount_Spatial)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5, size=0.8) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank())
  
  toPlot<-toPlot[order(toPlot[,colsToColor[3]]),]
  p_mito <- ggplot(toPlot, aes(staticNr, percent.mt)) + 
    geom_violin(fill="gray80") + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5, size=0.8) +
    scale_color_manual(values=c("#00bfc4", "#F8766D")) +
    theme_classic() +
    theme(legend.position = "none", axis.title.x = element_blank())
  
  return(grid.arrange(p_nGene,p_nUMI,p_mito,ncol=3))
}


##### Function drawVlnPlot_split
drawVlnPlot_split<-function(toPlot, colsToColor){
  toPlot$staticNr=1
  
  columnNr<-which(colnames(toPlot)==colsToColor[1])
  toPlot<-toPlot[order(toPlot[,columnNr]),]
  p_nGene <- ggplot(toPlot, aes(orig.ident, nFeature_Spatial, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[1]), alpha=0.5, size=0.8) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank())
  
  columnNr<-which(colnames(toPlot)==colsToColor[2])
  toPlot<-toPlot[order(toPlot[,columnNr]),]
  p_nUMI <- ggplot(toPlot, aes(orig.ident, nCount_Spatial, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[2]), alpha=0.5, size=0.8) +
    geom_violin() + 
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size=8), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank())
  
  columnNr<-which(colnames(toPlot)==colsToColor[3])
  toPlot<-toPlot[order(toPlot[,columnNr]),]
  p_mito <- ggplot(toPlot, aes(orig.ident, percent.mt, fill=orig.ident)) + 
    geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor[3]), alpha=0.5, size=0.8) +
    geom_violin() +
    scale_color_manual(values=c("#00BFC4", "#F8766D")) +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
          axis.text.y = element_text(size = 8), axis.title.x = element_blank())
  
  thePlot<-grid.arrange(p_nGene, p_nUMI,p_mito, ncol=3)
  return(thePlot)
}

