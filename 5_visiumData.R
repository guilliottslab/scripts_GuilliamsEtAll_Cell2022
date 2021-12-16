
################################################################################
########## LOAD PACKAGES
################################################################################

library('Seurat')
library('ggplot2')
library('patchwork')
library('dplyr')
library('gridExtra')

################################################################################
########## FUNCTIONS
################################################################################
source('/PATH/TO/1_Functions_spatial.R')

################################################################################
########## LOAD DATA
################################################################################

##### Sample1 #####
seuratObjSample1<-Load10X_Spatial(data.dir = '/PATH/TO/sample1/')
dim(seuratObjSample1)

##### Sample2 #####
seuratObjSample2<-Load10X_Spatial(data.dir = '/PATH/TO/sample2/')
dim(seuratObjSample2)

##### Sample3 #####
seuratObjSample3<-Load10X_Spatial(data.dir = '/PATH/TO/sample3/')
dim(seuratObjSample3)

########################################
##### Check for outlier spots
########################################

##### Sample1: no outlier spots #####
coordsSample1<-seuratObjSample1@images$slice1@coordinates

p1<-SpatialDimPlot(seuratObjSample1)
p2 <- ggplot(aes(x=imagecol, y=imagerow),data=coordsSample1)+
  geom_point(size=0.5) +
  theme_classic()
grid.arrange(p1,p2,ncol=2)
ggsave(grid.arrange(p1,p2,ncol=2), file='/PATH/TO/results/image_sample1.png', width=10, height = 5)


##### Sample2: remove outlier spots #####
coordsSample2<-seuratObjSample2@images$slice1@coordinates

p1<-SpatialDimPlot(seuratObjSample2)
p2 <- ggplot(aes(x=imagecol, y=imagerow),data=coordsSample2)+
  geom_point(size=0.5) +
  theme_classic()
grid.arrange(p1,p2,ncol=2)

# Select outlier spots
badSpots<-coordsSample2 %>% mutate('spot'=rownames(.)) %>% 
  dplyr::filter(imagecol < 5000 | imagerow < 4500) %>%
  pull(spot)

# Visualize bad spots
coordsSample2$badSpot<-FALSE
coordsSample2[badSpots,'badSpot']<-TRUE
p2 <- ggplot(aes(x=imagecol, y=imagerow,color=badSpot),data=coordsSample2)+
  geom_point(size=0.5) +
  theme_classic()
grid.arrange(p1,p2,ncol=2)

# Remove bad spots
seuratObjSample2<-seuratObjSample2[,setdiff(colnames(seuratObjSample2),badSpots)]
dim(seuratObjSample2)
# 31053  1274
p3<-SpatialDimPlot(seuratObjSample2)
ggsave(grid.arrange(p1,p2,p3,ncol=3), file='/PATH/TO/results/image_sample2.png', 
       width=15, height = 5)


##### Sample3: no outlier spots #####
coordsSample3<-seuratObjSample3@images$slice1@coordinates

p1<-SpatialDimPlot(seuratObjSample3)
p2 <- ggplot(aes(x=imagecol, y=imagerow),data=coordsSample3)+
  geom_point(size=0.5) +
  theme_classic()
grid.arrange(p1,p2,ncol=2)
ggsave(grid.arrange(p1,p2,ncol=2), file='/PATH/TO/results/image_sample3.png', width=10, height = 5)



########################################
##### Normalize
########################################
seuratObjSample1 <- SCTransform(seuratObjSample1, assay = "Spatial", verbose = FALSE)
seuratObjSample2 <- SCTransform(seuratObjSample2, assay = "Spatial", verbose = FALSE)
seuratObjSample3 <- SCTransform(seuratObjSample3, assay = "Spatial", verbose = FALSE)


########################################
##### Merge
########################################
seuratObjSample1@meta.data$orig.ident<-'sample1'
seuratObjSample2@meta.data$orig.ident<-'sample2'
seuratObjSample3@meta.data$orig.ident<-'sample3'

seuratObjMerge <- merge(seuratObjSample1, list(seuratObjSample2,seuratObjSample3))
dim(seuratObjMerge)


head(seuratObjMerge@meta.data)
tail(seuratObjMerge@meta.data)


table(seuratObjMerge@meta.data$orig.ident)




################################################################################
########## FILTERING
################################################################################
seuratObjMerge[["percent.mt"]] <- PercentageFeatureSet(seuratObjMerge, pattern = "^mt-")

### Create vlnPlot
plot1 <- VlnPlot(seuratObjMerge, 
                 features = c("nFeature_Spatial","nCount_Spatial","percent.mt"), 
                 pt.size = 0.1, group.by = 'orig.ident') + NoLegend()
print(plot1)
ggsave(plot1,file='/PATH/TO/results/vlnPlot_1a.png')

### Plot on image
plot2 <- SpatialFeaturePlot(seuratObjMerge, features = "percent.mt") + 
  theme(legend.position = "right")
print(plot2)
ggsave(plot2,file='/PATH/TO/results/colorMito_onImage.png', width = 10, height = 6)


### Color spots to be removed
toPlot<-seuratObjMerge@meta.data
toPlot$badCell1<-FALSE
toPlot[toPlot$nFeature_Spatial<1000,'badCell1']<-TRUE
toPlot$badCell2<-FALSE
toPlot[toPlot$nCount_Spatial<5000,'badCell2']<-TRUE
toPlot[toPlot$nCount_Spatial>60000,'badCell2']<-TRUE
toPlot$badCell3<-FALSE
toPlot[toPlot$percent.mt>15,'badCell3']<-TRUE


p1<-drawVlnPlot(toPlot,colsToColor = c('badCell1','badCell2','badCell3'))
print(p1)
ggsave(p1,file='/PATH/TO/results/vlnPlot_1b.png')

p2<-drawVlnPlot_split(toPlot,colsToColor = c('badCell1','badCell2','badCell3'))
print(p2)
ggsave(p2,file='/PATH/TO/results/vlnPlot_1c.png', width = 8)



### Do filtering
seuratObjSlice<-subset(seuratObjMerge, 
                       nFeature_Spatial > 1000 & nCount_Spatial > 5000 & 
                         nCount_Spatial < 60000 & percent.mt < 15)
dim(seuratObjSlice)


seuratObjMerge<-seuratObjSlice


################################################################################
########## UMAP
################################################################################
DefaultAssay(seuratObjMerge) <- "SCT"
VariableFeatures(seuratObjMerge) <- c(VariableFeatures(seuratObjSample1), 
                                      VariableFeatures(seuratObjSample2),
                                      VariableFeatures(seuratObjSample3))

seuratObjMerge <- RunPCA(seuratObjMerge, verbose = FALSE)

seuratObjMerge <- FindNeighbors(seuratObjMerge, dims = 1:30)
seuratObjMerge <- FindClusters(seuratObjMerge, resolution = 0.6)

seuratObjMerge <- RunUMAP(seuratObjMerge, dims = 1:30)


seuratObjMerge[["percent.mt"]] <- PercentageFeatureSet(seuratObjMerge, pattern = "^mt-")
seuratObjMerge[["percent.rbc"]] <- PercentageFeatureSet(seuratObjMerge, pattern = "Hb[ab]-")
plotMito<-FeaturePlot(seuratObjMerge,c("percent.mt"))
plotRBC<-FeaturePlot(seuratObjMerge,c("percent.rbc"))

########################################
##### Plot on umap
########################################

### Plot clusters
umapPlot<-DimPlot(seuratObjMerge, reduction = "umap", label=T)
umapPlotSplit<-DimPlot(seuratObjMerge, reduction = "umap", group.by = 'orig.ident')
grid.arrange(umapPlot, umapPlotSplit, ncol=2)

### Plot mito
plotMito<-FeaturePlot(seuratObjMerge,'percent.mt', min.cutoff = 'q2', max.cutoff = 'q98')
grid.arrange(umapPlot, umapPlotSplit, plotMito, ncol=3)

### Plot RBC
plotRBC<-FeaturePlot(seuratObjMerge,'percent.rbc', min.cutoff = 'q2', max.cutoff = 'q98')
grid.arrange(umapPlot, umapPlotSplit, plotMito, plotRBC, ncol=4)
ggsave(grid.arrange(umapPlot, umapPlotSplit, plotMito, plotRBC, ncol=4),
       file='/PATH/TO/results/umap.png', width = 24, height = 6)


### Plot some genes
FeaturePlot(seuratObjMerge,c("Glul", "Cyp2f2"), min.cutoff = 'q2', max.cutoff = 'q98')





########################################
##### Plot on image
########################################

### Plot clusters
SpatialDimPlot(seuratObjMerge)
SpatialDimPlot(seuratObjMerge,label = TRUE, label.size = 3)


#Plot certain clusters
clustersOI<-c(6,9,5,8,3,11)

SpatialDimPlot(seuratObjMerge, 
               cells.highlight = CellsByIdentities(object = seuratObjMerge, idents = clustersOI), 
               facet.highlight = TRUE, images = 'slice1',ncol=6)
SpatialDimPlot(seuratObjMerge, 
               cells.highlight = CellsByIdentities(object = seuratObjMerge, idents = clustersOI), 
               facet.highlight = TRUE, images = 'slice1.1',ncol=6)
SpatialDimPlot(seuratObjMerge, 
               cells.highlight = CellsByIdentities(object = seuratObjMerge, idents = clustersOI), 
               facet.highlight = TRUE, images = 'slice1.2',ncol=6)

### Plot mito
SpatialFeaturePlot(seuratObjMerge, features = "percent.mt") + 
  theme(legend.position = "right")


### Plot genes
SpatialFeaturePlot(seuratObjMerge, features = c("Glul", "Cyp2f2"))


################################################################################
########## GET DE GENES
################################################################################

allMarkers <- FindAllMarkers(seuratObjMerge)


### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  markersList[[i]]<-tmp[order(tmp$avg_logFC, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel
library('openxlsx')
write.xlsx(markersList, file ="/PATH/TO/results/DEgenes.xlsx")


################################################################################
########## SAVE
################################################################################

### Save seuratObj
saveRDS(seuratObjMerge, file='/PATH/TO/Robjects/seuratObj_sample1-3.rds')

### Save metaData table
toPlot<-seuratObjMerge@meta.data[,c('orig.ident','seurat_clusters')]
toPlot$spot<-rownames(toPlot)
umapTable<-as.data.frame(seuratObjMerge@reductions$umap@cell.embeddings)
toPlot<-cbind(toPlot,umapTable)
dim(toPlot)





