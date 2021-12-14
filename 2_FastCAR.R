
################################################################################
########## LOAD PACKAGES
################################################################################

library('dplyr')
library('gridExtra')
library('scater')
library('qlcMatrix')
library('FastCAR')
library('ggplot2')
library('Matrix')

################################################################################
########## FUNCTIONS
################################################################################
source('/PATH/TO/1_Functions.R')


##### Function to do calculations for ambient plot
calculateForPlot<-function(fullMatrix, stopLimit, byStep){
  print('Getting ambient profile')
  ambProfile = describe.ambient.RNA.sequence(fullCellMatrix = fullMatrix, 
                                             start = 0, 
                                             stop = stopLimit, 
                                             by = byStep, 
                                             contaminationChanceCutoff = 0.05)
  return(ambProfile)
}


##### Function to make ambient plot
makePlot<-function(ambProfile, theCutOff, sampleName){
  ### Make plot
  print('Making plot')
  png(file=paste0(sampleFolder,'results_fastCAR/plotAmbientProfile_',sampleName,'.png'), width = 900, height = 800)
  par(mfrow = c(3, 1))
  plot(as.numeric(rownames(ambProfile)), ambProfile[,1], 
       main = "Total number of empty droplets at cutoffs", 
       xlab = "empty droplet UMI cutoff", ylab = "Number of empty droplets")
  abline(v=theCutOff, col="red")
  
  plot(as.numeric(rownames(ambProfile)), ambProfile[,2], 
       main = "Number of genes in ambient RNA", xlab = "empty droplet UMI cutoff", 
       ylab = "Genes in empty droplets")
  abline(v=theCutOff, col="red")
  
  plot(as.numeric(rownames(ambProfile)), ambProfile[,3], 
       main = "number of genes to correct", xlab = "empty droplet UMI cutoff", 
       ylab = "Genes identified as contamination")
  abline(v=theCutOff, col="red")
  dev.off()
}

##### Function to get ambient profile
getAmbientProfile<-function(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff){
  ### Determine ambient RNA
  ambientProfile = determine.background.to.remove(fullMatrix, cellMatrix, emptyDropletCutoff, contaminationChanceCutoff)
  topGenes<-tail(sort(ambientProfile), 50)
  
  print(paste0(sum(ambientProfile>0),' genes using for correction'))
  
  ##Create barplot of top genes
  toPlot<-as.data.frame(topGenes)
  toPlot$gene<-rownames(toPlot)
  p<-ggplot(data=toPlot, aes(x=reorder(gene, -topGenes), y=topGenes)) +
    geom_bar(stat="identity",color='gray20',fill='steelblue') +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y =  element_text(size=6))
  ggsave(p, file=paste0(sampleFolder,"results_fastCAR/barplot_topAmbientGenes_",sampleName,".png"))
  
  ##Write all genes used for correction
  corrGenes<-ambientProfile[ambientProfile>0]
  corrGenes<-corrGenes[order(corrGenes, decreasing = T)]
  write.table(as.matrix(corrGenes), 
              file=paste0(sampleFolder,'results_fastCAR/correctionGenes_',sampleName,'.txt'),
              sep='\t', col.names = F)
  
  ###Return ambientProfile
  return(ambientProfile)
}


################################################################################
########## CALCULATE AMBIENT RNA
################################################################################
sampleName <- "SampleName"

########## Read data ##########
cellMatrix = read.cell.matrix("/PATH/TO/SampleName/filtered_feature_bc_matrix")
## ADD next line when having RNA and ADT data to extract only RNA data
## cellMatrix = cellMatrix$`Gene Expression`
dim(cellMatrix)


fullMatrix = read.full.matrix("/PATH/TO/SampleName/raw_feature_bc_matrix")
## ADD next line when having RNA and ADT data to extract only RNA data
## fullMatrix = fullMatrix$`Gene Expression`
dim(fullMatrix)

### Filter raw folder
min(colSums(cellMatrix))

the_colSums<-colSums(fullMatrix)
neededCells<-names(the_colSums[the_colSums > 10])
fullMatrix<-fullMatrix[,neededCells]
dim(fullMatrix)


########## Make plot ##########
ambProfile<-calculateForPlot(fullMatrix, 1200, 10)

### Determine cut off from the Cellranger web summary plot (for every sample separate)
theCutOff<-80
makePlot(ambProfile, theCutOff, sampleName)

### How many empty droplets
tmp<-names(the_colSums[the_colSums < theCutOff])
length(intersect(tmp, neededCells))


########## Get ambient genes ##########
ambientProfile = getAmbientProfile(fullMatrix, cellMatrix, theCutOff, 0.05)


########## Remove ambient genes ##########
cellMatrix_new = remove.background(cellMatrix, ambientProfile)
dim(cellMatrix_new)

########## Save new count matrix ##########
saveRDS(cellMatrix_new, file="/PATH/TO/cellMatrixNew_",sampleName,".rds")




