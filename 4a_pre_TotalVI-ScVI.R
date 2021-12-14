
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

################################################################################
########## FUNCTIONS
################################################################################
source('/PATH/TO/1_Functions.R')

################################################################################
########## Create list with sample names
################################################################################

listLabels <- c('Sample1','Sample2')

################################################################################
########## LOAD SEURATOBJ
################################################################################

seuratObj <- readRDS(file="/PATH/TO/seuratObj.rds")

## Get Counts matrix
rawData <- seuratObj@assays$RNA@counts
dim(rawData)

################################################################################
########## WRITE NEEDED FILES FOR TOTALVI/SCVI: RNA
################################################################################

########## Create RNA file ##########
### Function to write files
writeToMtx<-function(theData, sampleName){
  print(paste0('Writing mtx files for ',sampleName))
  dir.create(paste0("/PATH/TO/",sampleName))
  
  writeMM(obj = theData, file=paste0("/PATH/TO/",sampleName,"/matrix.mtx"))
  write(x = rownames(theData), file = paste0("/PATH/TO/",sampleName,"/features.tsv"))
  write(x = colnames(theData), file = paste0("/PATH/TO/",sampleName,"/barcodes.tsv"))
  gzip(paste0("/PATH/TO/",sampleName,"/features.tsv"))
  gzip(paste0("/PATH/TO/",sampleName,"/barcodes.tsv"))
  gzip(paste0("/PATH/TO/",sampleName,"/matrix.mtx"))
  print('done')
}



### Sample1
neededCells<-colnames(rawData)[grep('-1',colnames(rawData))]
rawData1<-rawData1[neededGenes,neededCells]
writeToMtx(rawData1, 'Sample1')

### Sample2
neededCells<-colnames(rawData)[grep('-2',colnames(rawData))]
rawData2<-rawData2[neededGenes,neededCells]
writeToMtx(rawData2, 'Sample2')


################################################################################
########## WRITE NEEDED FILES FOR TOTALVI: ADT
################################################################################

## Some of the samples were run to an older ADT whitelist than others.
## We changed to old AB names to the new AB names.
## How to change the old AB names into the correct AB name, can be found in the 4e_TranslateABs.txt file 

########## Sample1 ##########
rawDataSparse <- Read10X("PATH/TO/Sample1/filtered_feature_bc_matrix/")
rawDataADT1<-rawDataSparse$`Antibody Capture`

## Add ID number to the cellnames
theID<-which(unlist(listLabels)=='Sample1')
colnames(rawDataADT1)<-paste0(colnames(rawDataADT1),'-',theID)

## Filter rawdatADT
rawDataADT1<-rawDataADT1[,intersect(colnames(rawDataADT1),colnames(rawData))]
dim(rawDataADT1)


## Calculate in how many cells each AB is expressed
exprsInNrCells<-apply(rawDataADT1, 1, function(x){sum(x>0)})
exprsInNrCells<-exprsInNrCells[order(exprsInNrCells, decreasing = T)]

## Show in how many cells each AB expressed
barplot(exprsInNrCells)
cutoff<-200  ## Clear drop is visible at a certain cutoff, this should match with the actually added ABs
expressedABs1<-names(exprsInNrCells[exprsInNrCells>=cutoff])
length(expressedABs1)

### Only keep the ABs that were actually added
rawDataADT1<-rawDataADT1[expressedABs1,]
dim(rawDataADT1)


### Change underscores
rownames(rawDataADT1)<-gsub('_','-',rownames(rawDataADT1))
rownames(rawDataADT1)<-paste0('adt-',rownames(rawDataADT1))
write.csv(rawDataADT1, file="PATH/TO/rawDataADT1.csv")


########## Sample2 ##########
rawDataSparse <- Read10X("PATH/TO/Sample2/filtered_feature_bc_matrix/")
rawDataADT2<-rawDataSparse$`Antibody Capture`

## Add ID number to the cellnames
theID<-which(unlist(listLabels)=='Sample2')
colnames(rawDataADT2)<-paste0(colnames(rawDataADT2),'-',theID)

## Filter rawdatADT
rawDataADT2<-rawDataADT2[,intersect(colnames(rawDataADT2),colnames(rawData))]
dim(rawDataADT2)


## Calculate in how many cells each AB is expressed
exprsInNrCells<-apply(rawDataADT2, 1, function(x){sum(x>0)})
exprsInNrCells<-exprsInNrCells[order(exprsInNrCells, decreasing = T)]

## Show in how many cells each AB expressed
barplot(exprsInNrCells)
cutoff<-210  ## Clear drop is visible at a certain cutoff, this should match with the actually added ABs
expressedABs2<-names(exprsInNrCells[exprsInNrCells>=cutoff])
length(expressedABs2)

### Only keep the ABs that were actually added
rawDataADT2<-rawDataADT2[expressedABs2,]
dim(rawDataADT2)


### Change underscores and add -adt to the ABs names
rownames(rawDataADT2)<-gsub('_','-',rownames(rawDataADT2))
rownames(rawDataADT2)<-paste0('adt-',rownames(rawDataADT2))
write.csv(rawDataADT2, file="/PATH/TO/rawDataADT2.csv")

################################################################################
########## In case of non-CiteSeq + CiteSeq samples
################################################################################

### EXAMPLE
### Sample1 = non-CiteSeq
### Sample2 = CiteSeq 100 ABs added (WhiteList = 180 ABs)
### Sample3 = CiteSeq 120 ABs added (WhiteList = 180 ABs)


### 1) Preparing RNA data stays the same as shown above for the 3 samples
### 2) load the rawDataADT and remove per sample the antibodies (ABs) that weren't added (See above) and save. 
###    For Sample2 100 ABs were added and for Sample3 120 ABs were added

write.csv(rawDataADT2, file="PATH/TO/rawDataADT2.csv")
write.csv(rawDataADT3, file="PATH/TO/rawDataADT3.csv")

### 3) Create txt file containing the ABs overlap between the two CiteSeq samples (needed for step '' from the 4c_TotalVI_CiteSeq+non-CiteSeq.ipynb script) 

# All ABs names in the count table / WhiteList
allABs <- rownames(rawDataADT2)
length(allABs)
# 180

# ADT names expressed in all CiteSeq samples
overlapABs<-Reduce(intersect, list(expressedABs_Sample2,expressedABs_Sample3))
length(overlapABs)
# 100

write.table(overlapABs, file="/PATH/TO/overlapABs.txt", row.names = F,
            col.names = F, quote = F)


### 4) Create txt file containing the ABs difference between the two CiteSeq samples

ABs_Sample3 <- setdiff(expressedABs_Sample3,overlapABs)
length(ABs_Sample3)
# 20

write.table(ABs_Sample3, file="/PATH/TO/ABs_Sample3.txt", row.names = F,
            col.names = F, quote = F)


### 5) Create 'Dummy' ABs count table for non-CiteSeq Sample1

rawDataADT1<-rawData1[1:2,] ## Take First 2 rows from RNA count table for all the cells (All counts will be set to 0 in TotalVI)
rownames(rawDataADT1)<-c("adt-CD4","adt-CD8a") ## Assign random ADT names to the 2 rows. These counts will be changed into 0 in the 4c_TotalVI_CiteSeq+non-CiteSeq.ipynb script.

write.csv(rawDataADT1, file="PATH/TO/rawDataADT1.csv")







