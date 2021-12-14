
################################################################################
########## LOAD PACKAGES
################################################################################

library('Seurat')
library('dplyr')
library('ggplot2')
library('Matrix')
library('R.utils')
library('gridExtra')
library('grid')
library("reticulate")

################################################################################
########## FUNCTIONS
################################################################################
source('/PATH/TO/1_Functions.R')


################################################################################
########## LOAD TOTALVI/ScVI RESULTS
################################################################################
py <- import_builtins()
pickle <- import("pickle")

##### Load post_adata #####
post_adata <- pickle$load(py$open("/PATH/TO/post_adata.pkl"), "rb")

##### Load HVG dataset #####
all_dataset <- pickle$load(py$open("/PATH/TO/HVG_dataset.pkl"), "rb")
length(all_dataset$barcodes)
length(all_dataset$gene_names)

##### Load Full dataset #####
all_dataset_full <- pickle$load(py$open("/PATH/TO/all_dataset.pkl"), "rb")
length(all_dataset_full$barcodes)
length(all_dataset_full$gene_names)

##### Get umap coords #####
umapTable <- post_adata$obsm['X_umap']
dim(umapTable)
rownames(umapTable)<-all_dataset$barcodes
colnames(umapTable)<-c('UMAP_1','UMAP_2')
umapTable<-as.data.frame(umapTable)

##### Get HVG #####
hvgGenes<-rownames(post_adata$var)
length(hvgGenes)
# 4000
head(hvgGenes)

##### Get cluster IDs #####
obs<-post_adata$obs
dim(obs)

clusterIDs <- obs$louvain
names(clusterIDs)<-all_dataset$barcodes

batchIDs<-obs$sample
names(batchIDs)<-all_dataset$barcodes

typeIDs<-obs$type
names(typeIDs)<-all_dataset$barcodes

##### Create metadata variable #####
clusterMatrixVI<-umapTable
clusterMatrixVI$cluster<-as.factor(clusterIDs[rownames(clusterMatrixVI)])
clusterMatrixVI$sample<-as.factor(batchIDs[rownames(clusterMatrixVI)])
clusterMatrixVI$type<-as.factor(typeIDs[rownames(clusterMatrixVI)])

##### Get denoised proteins (Only for TotalVI) #####
badColIDs<-c(which(colnames(obs) %in% c('louvain','sample','type')),
             grep('fore_prob', colnames(obs)),
             grep('observed', colnames(obs)))
badColIDs<-badColIDs[order(badColIDs)]

denoisedProteins<-t(obs[,-badColIDs])
dim(denoisedProteins)

colnames(denoisedProteins)<-all_dataset$barcodes

##### Get denoised genes (For TotalVI and ScVI) #####
denoisedGenes<-t(post_adata$layers['norm_genes'])
dim(denoisedGenes)

colnames(denoisedGenes)<-all_dataset$barcodes
rownames(denoisedGenes)<-hvgGenes


##### Get normalised values full dataset #####
normData<-readMM(file="PATH/TO/normData.mtx")
normData<-t(normData)
dim(normData)

rownames(normData)<-all_dataset_full$gene_names
colnames(normData)<-all_dataset_full$barcodes
normData[1:5,1:5]


########################################
##### Create umap plot
########################################

### Make umap
centers <- clusterMatrixVI %>% dplyr::group_by(cluster) %>% summarize(x = median(x = UMAP_1), y = median(x = UMAP_2))
ggplot(aes(x=UMAP_1,y=UMAP_2, colour=cluster), data=clusterMatrixVI) +
  geom_point(size=0.2, shape=20) +
  geom_text(data = centers, aes(x = x, y = y, label = cluster), color='black', fontface='bold', size=6) +
  theme_classic() +
  theme(legend.position = "none")


########################################
##### Check AB expression - raw counts
########################################

myABs<-rownames(rawDataADT)
length(myABs)

for(i in 1:length(myABs)){
  ### Color the cells where an AB is 'nicely expressed'.
  ### 'nicely expressed' is defined as: the cells where the AB counts are higher that the 98% quantile cutoff
  ### (not taking the cells with 0 count into account)
  ABname<-myABs[i]
  print(paste0(i,'. Working on ',ABname))
  countsGene<-rawDataADT[ABname,]
  cutoff<-quantile(countsGene[countsGene>0],0.98)
  if(! is.na(cutoff)){
    wantedCells<-names(countsGene[countsGene>=cutoff])
    p1<-colorSomeCells(clusterMatrixVI, umapTable, wantedCells,ABname)
    fileName<-paste0('plot_',ABname,'_countQuantile.png')
    ggsave(p1, file=paste0("/PATH/TO/featurePlot_AB/",fileName))
  }
}



########################################
##### Check for AB pattern
########################################

### Check if there are cells that express each AB

meanPerCell<-apply(denoisedProteins,2,mean)
meanPerCell[meanPerCell>=quantile(meanPerCell,0.99)]<-quantile(meanPerCell,0.99)
meanPerCell<-meanPerCell[order(meanPerCell, decreasing = T)]

toPlot<-as.data.frame(meanPerCell)
toPlot$staticNr<-1
myCutoff<-5  ## this cutoff changes per sample.
cutOff<-median(toPlot$meanPerCell)+(myCutoff*mad(toPlot$meanPerCell))
toPlot$tooHigh<-FALSE
toPlot[toPlot$meanPerCell >= cutOff,'tooHigh']<-TRUE
colsToColor<-'tooHigh'
ggplot(toPlot, aes(staticNr, meanPerCell)) + 
  geom_jitter(height = 0, width = 0.3, aes_string(col=colsToColor), alpha=0.5) +
  geom_violin(fill="gray80") + 
  scale_color_manual(values=c("#00bfc4", "#F8766D")) +
  theme_classic()


badCells<-rownames(toPlot[toPlot$tooHigh==TRUE,])
length(badCells)
colorSomeCells(clusterMatrixVI, umapTable, badCells,'Bad cells')



