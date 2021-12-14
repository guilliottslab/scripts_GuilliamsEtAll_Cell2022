
################################################################################
########## LOAD PACKAGES
################################################################################

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('DisneyTools')

################################################################################
########## FUNCTIONS
################################################################################
source('/PATH/TO/1_Functions.R')

################################################################################
########## LOAD DATA
################################################################################

##### Sample1 #####
rawData1 <- readRDS("/PATH/TO/cellMatrixNew_Sample1.rds")
dim(rawData1)

## Change colnames
colnames(rawData1) <- paste0(colnames(rawData1), "-1")


##### Sample2 #####
rawData2 <- readRDS("/PATH/TO/cellMatrixNew_Sample2.rds")
dim(rawData2)

## Change colnames
colnames(rawData2) <- paste0(colnames(rawData2), "-2")


##### Combine samples #####
rawData <- cbind(rawData1,rawData2)
dim(rawData)


##### Save rawData #####
rawData <- saveRDS(rawData, file="PATH/TO/rawData.rds")


################################################################################
########## QC: CELLS
################################################################################

########################################
########## Calculate metrics
########################################

##### Create object #####
library("scater")
sce<-SingleCellExperiment(list(counts=rawData))
dim(sce)

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)

rownames(sce)[is.mito]

##### Calculate QC metrics #####
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
dim(colData(sce))


##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(rawData),"orig.ident"=listLabels[[1]],
                     "nGene"=sce$total_features_by_counts,"nUMI"=sce$total_counts,"percent.mito"=sce$pct_counts_Mt, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1


##### Add sample names to metaData #####
for(i in 2:length(listLabels)){
  toSearch<-paste0("-",i)
  metaData[grep(toSearch,rownames(metaData)), which(colnames(metaData)=="orig.ident")]<-listLabels[[i]]
}


########################################
########## Get outliers
########################################

##### Set number of MADs (This varies between the different datasets) #####
nmad_low_feature<-3
nmad_high_feature<-3

nmad_low_UMI<-3
nmad_high_UMI<-3

nmad_high_mito<-3

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop.low <- isOutlier(sce$total_features_by_counts, nmads=nmad_low_feature, type="lower", log=TRUE)
sum(feature.drop.low)

feature.drop.high <- isOutlier(sce$total_features_by_counts, nmads=nmad_high_feature, type="higher", log=TRUE)
sum(feature.drop.high)

feature.drop<-as.logical(feature.drop.low + feature.drop.high)
sum(feature.drop)

##same as UMI in Seurat pipeline
libsize.drop.low <- isOutlier(sce$total_counts, nmads=nmad_low_UMI, type="lower", log=TRUE)
sum(libsize.drop.low)

libsize.drop.high <- isOutlier(sce$total_counts, nmads=nmad_high_UMI, type="higher", log=TRUE)
sum(libsize.drop.high)

libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
sum(libsize.drop)

##% mitochondrial genes
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=nmad_high_mito, type="higher")
sum(mito.drop)

##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop

########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))

toPlot<-metaData
savePlots<-TRUE

##nGene
if(savePlots==TRUE){png(file="/PATH/TO/1a_nGene.png", width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

##nUMI
if(savePlots==TRUE){png(file="/PATH/TO/1b_nUMI.png", width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

##percent.mito
if(savePlots==TRUE){png(file="/PATH/TO/1c_percMito.png", width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

########################################
########## Create violinPlots
########################################

### Before filtering
toPlot<-metaData
drawVlnPlot(toPlot, fileName = "/PATH/TO/2a_beforeFiltering.png", colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, fileName = "/PATH/TO/2a_beforeFiltering_splitted.png", colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))


drawVlnPlot(toPlot, fileName = "/PATH/TO/2a_beforeFiltering_nGene.png", 
            colsToColor = c('nGene.drop','nGene.drop','nGene.drop'))
drawVlnPlot(toPlot, fileName = "/PATH/TO/2a_beforeFiltering_nUMI.png", 
            colsToColor = c('nUMI.drop','nUMI.drop','nUMI.drop'))
drawVlnPlot(toPlot, fileName = "/PATH/TO/2a_beforeFiltering_mito.png", 
            colsToColor = c('mito.drop','mito.drop','mito.drop'))


### After filtering
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = "/PATH/TO/2b_afterFiltering.png", colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, fileName = "/PATH/TO/2b_afterFiltering_splitted.png", colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))


########################################
########## Remove outliers
########################################

sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
dim(sce)

### Number of cells removed
nrow(metaData)-ncol(sce)

################################################################################
########## FINALIZE QC
################################################################################

##### Filter rawData #####
rawDataFiltered<-rawData[rownames(sce),colnames(sce)]
dim(rawDataFiltered)

##### Save rawDataFiltered #####
saveRDS(rawDataFiltered, file="/PATH/TO/rawDataFiltered.rds")

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataFiltered, project = "seuratObj", min.cells = 3, min.features = 200)


################################################################################
########## ADD extra information to the Seurat metaData
################################################################################

##### Add % mitochondrial genes #####
seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

##### Create violin plot for 3 metrics (ngene,nUMI,%mito) #####
png(file="/PATH/TO/4_vlnPlotSeurat.png", width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"))
dev.off()

##### Add sample names to orig.ident #####
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
for(i in 1:length(listLabels)){
  toSearch<-paste0('-',i)
  metaDataTable[grep(toSearch,rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-listLabels[[i]]
}
seuratObj@meta.data<-metaDataTable


################################################################################
########## NORMALIZE SCT
################################################################################
seuratObj <- SCTransform(seuratObj,verbose = FALSE)


##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, file = "/PATH/TO/5_afterNorm_splitted.png")


################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

########################################
########## PCA PLOT
########################################
pdf(file="/PATH/TO/8a_PCA.pdf", width = 10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3))
dev.off()


########################################
########## HEATMAP OF PCs
########################################

### Create heatmap of PC 1-40
pdf(file="/PATH/TO/9a_selectPC.pdf")
PCHeatmap(seuratObj, dims = 1:12, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 13:24, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 25:36, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 37:40, cells = 500, balanced = TRUE)
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

##### Create PCElbowplot #####
png(file="/PATH/TO/9b_selectPC.png", width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 40)
dev.off()

## This varies between the different datasets

################################################################################
########## CLUSTER THE CELLS
################################################################################

## Number of PCs and resolution used
dimsToTry<-c(20)
resToUse<-0.8

for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse)
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0("/PATH/TO/10a_tSNE_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 21, height = 9)
  
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n_neighbors = 30)
  umapPlot<-DimPlot(seuratObj, reduction.use = "umap", label = T, label.size = 8, pt.size = 2)
  umapPlotSplit<-DimPlot(seuratObj, reduction.use = "umap", label = F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0("/PATH/TO/10b_UMAP_",min(dimsToUse),"-",max(dimsToUse),".png"), width = 21, height = 9)
  
}

##### Save object
saveRDS(seuratObj, file="/PATH/TO/seuratObj.rds")




################################################################################
########## GET DE GENES
################################################################################

########################################
##### all clusters vs all clusters
########################################

### Find markers for every cluster compared to all remaining cells, report only the positive ones
allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10, min.diff.pct=0.20, logfc.threshold = 0.20, return.thresh = 0.01, only.pos = F)
table(allMarkers$cluster)


### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1

  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC

  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel
library('openxlsx')
write.xlsx(markersList, file ="/PATH/TO/allClusters.xlsx")

### Test
FeaturePlot(object = seuratObj, features = c("ALB"), cols = c("grey", "blue"),
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98')


