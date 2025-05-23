
###======================================================================================
###Part I: Analyze the scRNA-seq data from HCC cohorts independently and show tSNE plots.
###======================================================================================

library(Seurat)
library(dplyr)
library(magrittr)
library(limma)

##---------------------------------------------------------------------GSE149614
#Normalization, dimensionality reduction and clustering
GSE149614.counts <- read.table("C:/Users/10104/Desktop/GSE149614_HCC.scRNAseq.S71915.count.txt")

GSE149614.metadata <- read.table("C:/Users/10104/Desktop/GSE149614_HCC.metadata.updated.txt")

GSE149614_data <- CreateSeuratObject(counts = GSE149614.counts, 
                                project = "GSE149614", 
                                meta.data = GSE149614.metadata,
                                min.cells = 3, min.features = 200)

GSE149614_data <- NormalizeData(GSE149614_data)

GSE149614_data <- FindVariableFeatures(object = GSE149614_data,selection.method = 'vst', nfeatures = 2000)

GSE149614_data <- ScaleData(GSE149614_data)    

GSE149614_data<- RunPCA(GSE149614_data, features = VariableFeatures(object = GSE149614_data))

GSE149614_data <- FindNeighbors(GSE149614_data,reduction = "pca", dims = 1:30)

GSE149614_data <- FindClusters(GSE149614_data, verbose = F, resolution = 2.0) 

GSE149614_data <- RunTSNE(GSE149614_data, reduction = "pca", dims = 1:30)

#Cell type annotation
GSE149614_data$cell.type.identification.step1 <- "Hepatocyte"
for(i in rownames(GSE149614_data@meta.data)){
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==0){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==1){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==3){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==11){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==16){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==21){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==24){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==28){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==29){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==54){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==56){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==58){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==61){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==17){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==31){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==35){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==55){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==62){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==64){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==68){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==69){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==10){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==18){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==41){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==49){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==60){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==25){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==32){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==65){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==67){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==4){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==6){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==13){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==14){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==22){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==23){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==26){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==30){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==37){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==38){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==40){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==46){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==47){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==48){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==50){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==51){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==52){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==57){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE149614_data@meta.data[i,"RNA_snn_res.2"]==66){GSE149614_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"}
}

#Show tSNE plot
DimPlot(GSE149614_data, group.by = "cell.type.identification.step1", pt.size = 0.01, label = F)

#Show THY1 distribution
FeaturePlot(GSE149614_data, features = "THY1", pt.size = 0.01)


##---------------------------------------------------------------------GSE156337
#Normalization, dimensionality reduction and clustering
GSE156337_data <- NormalizeData(GSE156337_data)

GSE156337_data <- FindVariableFeatures(object = GSE156337_data,selection.method = 'vst', nfeatures = 2000)

GSE156337_data <- ScaleData(GSE156337_data)    

GSE156337_data<- RunPCA(GSE156337_data, features = VariableFeatures(object = GSE156337_data))

GSE156337_data<- FindNeighbors(GSE156337_data, dims = 1:30)

GSE156337_data<- FindClusters(GSE156337_data, resolution = 2.0)

GSE156337_data<- RunTSNE(GSE156337_data, dims = 1:30)

#Cell type annotation
GSE156337_data$cell.type.identification.step1 <- "Hepatocyte"
for(i in rownames(GSE156337_data@meta.data)){
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==0){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==1){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==2){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==3){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==4){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==5){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==6){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==7){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==8){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==9){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==10){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==11){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==12){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==13){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==14){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==15){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==16){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==17){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==18){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==19){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==20){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==21){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==22){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==23){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==24){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==25){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==26){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==27){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==28){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==29){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==30){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==31){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==32){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==33){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==34){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==35){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==36){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==37){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==38){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==39){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==40){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==41){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==42){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==43){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==44){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==45){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==46){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Epithelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==47){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==48){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==49){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==50){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE156337_data@meta.data[i,"RNA_snn_res.2"]==51){GSE156337_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"}
}

#Show tSNE plot
DimPlot(GSE156337_data, group.by = "cell.type.identification.step1", pt.size = 0.01, label = F)

#Show THY1 distribution
FeaturePlot(GSE156337_data, features = "THY1", pt.size = 0.01)


##---------------------------------------------------------------------GSE202642
GSE202642.counts <- Read10X("/GSE202642/")

GSE202642_data <- CreateSeuratObject(counts = GSE202642.counts, 
                                project = "GSE202642", 
                                min.cells = 3, min.features = 200)

#Quality control
VlnPlot(GSE202642_data, features = c("nFeature_RNA"))

VlnPlot(GSE202642_data, features = c("nCount_RNA"))

GSE202642_data[["percent.mt"]] <- PercentageFeatureSet(GSE202642_data, pattern = "^MT-")

VlnPlot(GSE202642_data, features = c("percent.mt"))

GSE202642_data <- subset(GSE202642_data, subset = nFeature_RNA > 200 & nCount_RNA < 100000 & percent.mt < 10)

#Normalization, dimensionality reduction and clustering
GSE202642_data <- NormalizeData(GSE202642_data, normalization.method = "LogNormalize", scale.factor = 10000)

GSE202642_data <- FindVariableFeatures(GSE202642_data, selection.method = "vst", nfeatures = 2000)

GSE202642_data <- ScaleData(GSE202642_data)

GSE202642_data <- RunPCA(GSE202642_data, features = VariableFeatures(object = GSE202642_data))

GSE202642_data <- FindNeighbors(GSE202642_data,reduction = "pca", dims = 1:20)

GSE202642_data <- FindClusters(GSE202642_data, verbose = F, resolution = 1.0)   

GSE202642_data <- RunTSNE(GSE202642_data, reduction = "pca", dims = 1:20)

#Cell type annotation
GSE202642_data$cell.type.identification.step1<-0
for(i in rownames(GSE202642_data@meta.data)){
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==0){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==1){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==2){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==3){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==4){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==5){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==6){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==7){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==8){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==9){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==10){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==11){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==12){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==13){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==14){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==15){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Fibroblast"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==16){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==17){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==18){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==19){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==20){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==21){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==22){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==23){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==24){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==25){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==26){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==27){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==28){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==29){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"T/NK.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==30){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"B.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==31){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Endothelial.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==32){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==33){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==34){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==35){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==36){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Hepatocyte"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==37){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"};
  if(GSE202642_data@meta.data[i,"RNA_snn_res.1"]==38){GSE202642_data@meta.data[i,"cell.type.identification.step1"]<-"Myeloid.cell"}
}

#Show tSNE plot
DimPlot(GSE202642_data, group.by = "cell.type.identification.step1", pt.size = 0.01, label = F)

#Show THY1 distribution
FeaturePlot(GSE202642_data, features = "THY1", pt.size = 0.01)



###=================================================================================
###Part II: GO analysis of the DEGs of THY1.high cells compared with THY1.low cells.
###=================================================================================

library(org.Hs.eg.db)
library(clusterProfiler)
library(patchwork)
library(tidyverse)

##---------------------------------------------------------------------GSE149614
#Extract fibroblast and endothelial populations
GSE149614_data_fibroblast <- subset(GSE149614_data, cell.type.identification.step1 %in% c("Fibroblast"))

GSE149614_data_endothelial <- subset(GSE149614_data, cell.type.identification.step1 %in% c("Endothelial.cell"))

#Re-analyze the fibroblast population
GSE149614_data_fibroblast <- NormalizeData(GSE149614_data_fibroblast)

GSE149614_data_fibroblast <- FindVariableFeatures(object = GSE149614_data_fibroblast,selection.method = 'vst', nfeatures = 2000)

GSE149614_data_fibroblast <- ScaleData(GSE149614_data_fibroblast)    

GSE149614_data_fibroblast<- RunPCA(GSE149614_data_fibroblast, features = VariableFeatures(object = GSE149614_data_fibroblast))

GSE149614_data_fibroblast<- FindNeighbors(GSE149614_data_fibroblast, dims = 1:30)

GSE149614_data_fibroblast<- FindClusters(GSE149614_data_fibroblast, resolution = 2)

GSE149614_data_fibroblast<- RunTSNE(GSE149614_data_fibroblast, dims = 1:30)

DimPlot(GSE149614_data_fibroblast, reduction = "tsne", pt.size=0.1, label = F, label.size = 5)

#THY1 expression status annotation
GSE149614_data_fibroblast$THY1.status <- "THY1.high"
for(i in rownames(GSE149614_data_fibroblast@meta.data)){
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==3){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==4){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==9){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==11){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==13){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==15){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==16){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"};
  if(GSE149614_data_fibroblast@meta.data[i,"RNA_snn_res.2"]==22){GSE149614_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.low"}
}

#DEGs screening
GSE149614.fibroblast_THY1.HIGH_markers <- FindMarkers(GSE149614_data_fibroblast, 
                                                      ident.1 = "THY1.high",
                                                      group.by = "THY1.status",
                                                      only.pos = T)

filter_GSE149614.fibroblast_THY1.HIGH_markers <-subset(GSE149614.fibroblast_THY1.HIGH_markers, p_val_adj<0.01)

#Perform pathway enrichment analysis
ego_BP <- enrichGO(gene          = row.names(filter_GSE149614.fibroblast_THY1.HIGH_markers),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 

ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)

barplot(ego_BP,showCategory = 30)

#Re-analyze the endothelial population
GSE149614_data_endothelial <- NormalizeData(GSE149614_data_endothelial)

GSE149614_data_endothelial <- FindVariableFeatures(object = GSE149614_data_endothelial,selection.method = 'vst', nfeatures = 2000)

GSE149614_data_endothelial <- ScaleData(GSE149614_data_endothelial)    

GSE149614_data_endothelial<- RunPCA(GSE149614_data_endothelial, features = VariableFeatures(object = GSE149614_data_endothelial))

GSE149614_data_endothelial<- FindNeighbors(GSE149614_data_endothelial, dims = 1:30)

GSE149614_data_endothelial<- FindClusters(GSE149614_data_endothelial, resolution = 2)

GSE149614_data_endothelial<- RunTSNE(GSE149614_data_endothelial, dims = 1:30)

DimPlot(GSE149614_data_endothelial, reduction = "tsne", pt.size=0.1, label = F, label.size = 5)

#THY1 expression status annotation
GSE149614_data_endothelial$THY1.status <- "THY1.low"
for(i in rownames(GSE149614_data_endothelial@meta.data)){
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==0){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==1){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==17){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==20){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==22){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==23){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE149614_data_endothelial@meta.data[i,"RNA_snn_res.2"]==24){GSE149614_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"}
}


#DEGs screening
GSE149614.endothelial_THY1.HIGH_markers <- FindMarkers(GSE149614_data_endothelial, 
                                                       ident.1 = "THY1.high",
                                                       group.by = "THY1.status",
                                                       only.pos = T)

filter_GSE149614.endothelial_THY1.HIGH_markers <-subset(GSE149614.endothelial_THY1.HIGH_markers, p_val_adj<0.01)

#Perform pathway enrichment analysis
genelist<-bitr(rownames(filter_GSE149614.endothelial_THY1.HIGH_markers),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = 'org.Hs.eg.db')

genelist <- pull(genelist,ENTREZID)  

ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')

barplot(ekegg, showCategory=20)



##---------------------------------------------------------------------GSE156337
#Extract fibroblast and endothelial populations
GSE156337_data_fibroblast <- subset(GSE149614_data, cell.type.identification.step1 %in% c("Fibroblast"))

GSE156337_data_endothelial <- subset(GSE149614_data, cell.type.identification.step1 %in% c("Endothelial.cell"))

#Re-analyze the fibroblast population
GSE156337_data_fibroblast <- NormalizeData(GSE156337_data_fibroblast)

GSE156337_data_fibroblast <- FindVariableFeatures(object = GSE156337_data_fibroblast,selection.method = 'vst', nfeatures = 2000)

GSE156337_data_fibroblast <- ScaleData(GSE156337_data_fibroblast)    

GSE156337_data_fibroblast<- RunPCA(GSE156337_data_fibroblast, features = VariableFeatures(object = GSE156337_data_fibroblast))

GSE156337_data_fibroblast<- FindNeighbors(GSE156337_data_fibroblast, dims = 1:30)

GSE156337_data_fibroblast<- FindClusters(GSE156337_data_fibroblast, resolution = 2)

GSE156337_data_fibroblast<- RunTSNE(GSE156337_data_fibroblast, dims = 1:30)

DimPlot(GSE156337_data_fibroblast, reduction = "tsne", pt.size=0.1, label = F, label.size = 5)

#THY1 expression status annotation
GSE156337_data_fibroblast$THY1.status <- "THY1.low"
for(i in rownames(GSE156337_data_fibroblast@meta.data)){
  if(GSE156337_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==4){GSE156337_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE156337_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==7){GSE156337_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE156337_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==10){GSE156337_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"}
}

#DEGs screening
GSE156337.fibroblast_THY1.HIGH_markers <- FindMarkers(GSE156337_data_fibroblast, 
                                                      ident.1 = "THY1.high",
                                                      group.by = "THY1.status",
                                                      only.pos = T)

filter_GSE156337.fibroblast_THY1.HIGH_markers <-subset(GSE156337.fibroblast_THY1.HIGH_markers, p_val_adj<0.01)

#Perform pathway enrichment analysis
ego_BP <- enrichGO(gene          = row.names(filter_GSE156337.fibroblast_THY1.HIGH_markers),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 

ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)

barplot(ego_BP,showCategory = 30)

#Re-analyze the endothelial population
GSE156337_data_endothelial <- NormalizeData(GSE156337_data_endothelial)

GSE156337_data_endothelial <- FindVariableFeatures(object = GSE156337_data_endothelial,selection.method = 'vst', nfeatures = 2000)

GSE156337_data_endothelial <- ScaleData(GSE156337_data_endothelial)    

GSE156337_data_endothelial<- RunPCA(GSE156337_data_endothelial, features = VariableFeatures(object = GSE156337_data_endothelial))

GSE156337_data_endothelial<- FindNeighbors(GSE156337_data_endothelial, dims = 1:30)

GSE156337_data_endothelial<- FindClusters(GSE156337_data_endothelial, resolution = 2)

GSE156337_data_endothelial<- RunTSNE(GSE156337_data_endothelial, dims = 1:30)

DimPlot(GSE156337_data_endothelial, reduction = "tsne", pt.size=0.1, label = F, label.size = 5)

#THY1 expression status annotation
GSE156337_data_endothelial$THY1.status <- "THY1.low"
for(i in rownames(GSE156337_data_endothelial@meta.data)){
  if(GSE156337_data_endothelial@meta.data[i,"RNA_snn_res.2"]==9){GSE156337_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE156337_data_endothelial@meta.data[i,"RNA_snn_res.2"]==10){GSE156337_data_endothelial@meta.data[i,"THY1.status"]<-"THY1.high"}
}

#DEGs screening
GSE156337.endothelial_THY1.HIGH_markers <- FindMarkers(GSE156337_data_endothelial, 
                                                       ident.1 = "THY1.high",
                                                       group.by = "THY1.status",
                                                       only.pos = T)

filter_GSE156337.endothelial_THY1.HIGH_markers <-subset(GSE156337.endothelial_THY1.HIGH_markers, p_val_adj<0.01)

#Perform pathway enrichment analysis
genelist<-bitr(rownames(filter_GSE156337.endothelial_THY1.HIGH_markers),
               fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = 'org.Hs.eg.db')

genelist <- pull(genelist,ENTREZID)  

ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')

barplot(ekegg, showCategory=20)


##---------------------------------------------------------------------GSE202642
#Extract fibroblast and endothelial populations
GSE202642_data_fibroblast <- subset(GSE149614_data, cell.type.identification.step1 %in% c("Fibroblast"))

GSE202642_data_endothelial <- subset(GSE149614_data, cell.type.identification.step1 %in% c("Endothelial.cell"))

#Re-analyze the fibroblast population
GSE202642_data_fibroblast <- NormalizeData(GSE202642_data_fibroblast)

GSE202642_data_fibroblast <- FindVariableFeatures(object = GSE202642_data_fibroblast,selection.method = 'vst', nfeatures = 2000)

GSE202642_data_fibroblast <- ScaleData(GSE202642_data_fibroblast)    

GSE202642_data_fibroblast<- RunPCA(GSE202642_data_fibroblast, features = VariableFeatures(object = GSE202642_data_fibroblast))

GSE202642_data_fibroblast<- FindNeighbors(GSE202642_data_fibroblast, dims = 1:30)

GSE202642_data_fibroblast<- FindClusters(GSE202642_data_fibroblast, resolution = 1)

GSE202642_data_fibroblast<- RunTSNE(GSE202642_data_fibroblast, dims = 1:30)

DimPlot(GSE202642_data_fibroblast, reduction = "tsne", pt.size=0.1, label = F, label.size = 5)

#THY1 expression status annotation
GSE202642_data_fibroblast$THY1.status <- "THY1.low"
for(i in rownames(GSE202642_data_fibroblast@meta.data)){
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==1){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==3){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==5){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==8){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==9){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==10){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"};
  if(GSE202642_data_fibroblast@meta.data[i,"RNA_snn_res.1"]==11){GSE202642_data_fibroblast@meta.data[i,"THY1.status"]<-"THY1.high"}
}

#DEGs screening
GSE202642.fibroblast_THY1.HIGH_markers <- FindMarkers(GSE202642_data_fibroblast, 
                                                      ident.1 = "THY1.high",
                                                      group.by = "THY1.status",
                                                      only.pos = T)

filter_GSE202642.fibroblast_THY1.HIGH_markers <-subset(GSE202642.fibroblast_THY1.HIGH_markers, p_val_adj<0.01)

#Perform pathway enrichment analysis
ego_BP <- enrichGO(gene          = row.names(filter_GSE202642.fibroblast_THY1.HIGH_markers),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 

ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)

barplot(ego_BP,showCategory = 30)



###=======================================================
###Part III: Cellular interaction analysis using CellChat.
###=======================================================

library(CellChat)
library(ggplot2)
library(ggalluvial)
library(Seurat)
library(dplyr)
library(patchwork)
library(EnhancedVolcano) 

##---------------------------------------------------------------------GSE149614
#Data re-construction to reduce calculation 

#Cell type extraction
GSE149614_data_B <- subset(GSE149614_data, cell.type.identification.step1 %in% "B.cell")
GSE149614_data_Endothelial <- subset(GSE149614_data, cell.type.identification.step1 %in% "Endothelial.cell")
GSE149614_data_Hepatocyte <- subset(GSE149614_data, cell.type.identification.step1 %in% "Hepatocyte")
GSE149614_data_Myeloid <- subset(GSE149614_data, cell.type.identification.step1 %in% "Myeloid.cell")
GSE149614_data_T.NK <- subset(GSE149614_data, cell.type.identification.step1 %in% "T/NK.cell")

#Randomly select 1000 cells out of each cell type
GSE149614_data_B_ids <- colnames(GSE149614_data_B)
GSE149614_data_B_source <- factor(GSE149614_data_B$cell.type.identification.step1)
GSE149614_data_B_list <- lapply(split(GSE149614_data_B_ids, GSE149614_data_B_source),function(x){sample(x,1000)})            
GSE149614_data_B_id <- unlist(GSE149614_data_B_list)
GSE149614_data_B_sml <- GSE149614_data_B[ ,GSE149614_data_B_id]

GSE149614_data_Endothelial_ids <- colnames(GSE149614_data_Endothelial)
GSE149614_data_Endothelial_source <- factor(GSE149614_data_Endothelial$cell.type.identification.step1)
GSE149614_data_Endothelial_list <- lapply(split(GSE149614_data_Endothelial_ids, GSE149614_data_Endothelial_source),function(x){sample(x,1000)})            
GSE149614_data_Endothelial_id <- unlist(GSE149614_data_Endothelial_list)
GSE149614_data_Endothelial_sml <- GSE149614_data_Endothelial[ ,GSE149614_data_Endothelial_id]

GSE149614_data_Hepatocyte_ids <- colnames(GSE149614_data_Hepatocyte)
GSE149614_data_Hepatocyte_source <- factor(GSE149614_data_Hepatocyte$cell.type.identification.step1)
GSE149614_data_Hepatocyte_list <- lapply(split(GSE149614_data_Hepatocyte_ids, GSE149614_data_Hepatocyte_source), function(x){sample(x,1000)})            
GSE149614_data_Hepatocyte_id <- unlist(GSE149614_data_Hepatocyte_list)
GSE149614_data_Hepatocyte_sml <- GSE149614_data_Hepatocyte[ ,GSE149614_data_Hepatocyte_id]

GSE149614_data_Myeloid_ids <- colnames(GSE149614_data_Myeloid)
GSE149614_data_Myeloid_source <- factor(GSE149614_data_Myeloid$cell.type.identification.step1)
GSE149614_data_Myeloid_list <- lapply(split(GSE149614_data_Myeloid_ids, GSE149614_data_Myeloid_source),function(x){sample(x,1000)})            
GSE149614_data_Myeloid_id <- unlist(GSE149614_data_Myeloid_list)
GSE149614_data_Myeloid_sml <- GSE149614_data_Myeloid[ ,GSE149614_data_Myeloid_id]

GSE149614_data_T.NK_ids <- colnames(GSE149614_data_T.NK)
GSE149614_data_T.NK_source <- factor(GSE149614_data_T.NK$cell.type.identification.step1)
GSE149614_data_T.NK_list <- lapply(split(GSE149614_data_T.NK_ids, GSE149614_data_T.NK_source),function(x){sample(x,1000)})            
GSE149614_data_T.NK_id <- unlist(GSE149614_data_T.NK_list)
GSE149614_data_T.NK_sml <- GSE149614_data_T.NK[ ,GSE149614_data_T.NK_id]

#Edit metadata for each cell type
GSE149614_data_B_sml$THY1.status <- GSE149614_data_B_sml$cell.type.identification.step1
GSE149614_data_Endothelial_sml$THY1.status <- GSE149614_data_Endothelial_sml$cell.type.identification.step1
GSE149614_data_Hepatocyte_sml$THY1.status <- GSE149614_data_Hepatocyte_sml$cell.type.identification.step1
GSE149614_data_Myeloid_sml$THY1.status <- GSE149614_data_Myeloid_sml$cell.type.identification.step1
GSE149614_data_T.NK_sml$THY1.status <- GSE149614_data_T.NK_sml$cell.type.identification.step1

#Dataset re-construction
GSE149614_data_remerge.for.cellchat <- merge(x = GSE149614_data_fibroblast,
                                             y = c(GSE149614_data_B_sml,
                                                   GSE149614_data_Endothelial_sml,
                                                   GSE149614_data_Hepatocyte_sml,
                                                   GSE149614_data_Myeloid_sml,
                                                   GSE149614_data_T.NK_sml))

#Cellular interaction analysis
cellchat <- createCellChat(object = GSE149614_data_remerge.for.cellchat, group.by = "THY1.status")
cellchat
summary(cellchat)
meta <- data.frame(cellType = GSE149614_data_remerge.for.cellchat$THY1.status, row.names =  Cells(GSE149614_data_remerge.for.cellchat))
cellchat <- addMeta(cellchat, meta = meta,  meta.name = "THY1.status")
cellchat <- setIdent(cellchat, ident.use = "THY1.status")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human          
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
cellchat <- updateCellChat(cellchat)


##---------------------------------------------------------------------GSE156337
#Data re-construction to reduce calculation 

#Cell type extraction
GSE156337_data_B <- subset(GSE156337_data, cell.type.identification.step1 %in% "B.cell")
GSE156337_data_Endothelial <- subset(GSE156337_data, cell.type.identification.step1 %in% "Endothelial.cell")
GSE156337_data_Hepatocyte <- subset(GSE156337_data, cell.type.identification.step1 %in% "Hepatocyte")
GSE156337_data_Myeloid <- subset(GSE156337_data, cell.type.identification.step1 %in% "Myeloid.cell")
GSE156337_data_T.NK <- subset(GSE156337_data, cell.type.identification.step1 %in% "T/NK.cell")

#Randomly select 1000 cells out of each cell type
GSE156337_data_B_ids <- colnames(GSE156337_data_B)
GSE156337_data_B_source <- factor(GSE156337_data_B$cell.type.identification.step1)
GSE156337_data_B_list <- lapply(split(GSE156337_data_B_ids, GSE156337_data_B_source),function(x){sample(x,1000)})            
GSE156337_data_B_id <- unlist(GSE156337_data_B_list)
GSE156337_data_B_sml <- GSE156337_data_B[ ,GSE156337_data_B_id]

GSE156337_data_Endothelial_ids <- colnames(GSE156337_data_Endothelial)
GSE156337_data_Endothelial_source <- factor(GSE156337_data_Endothelial$cell.type.identification.step1)
GSE156337_data_Endothelial_list <- lapply(split(GSE156337_data_Endothelial_ids, GSE156337_data_Endothelial_source),function(x){sample(x,1000)})            
GSE156337_data_Endothelial_id <- unlist(GSE156337_data_Endothelial_list)
GSE156337_data_Endothelial_sml <- GSE156337_data_Endothelial[ ,GSE156337_data_Endothelial_id]

GSE156337_data_Hepatocyte_ids <- colnames(GSE156337_data_Hepatocyte)
GSE156337_data_Hepatocyte_source <- factor(GSE156337_data_Hepatocyte$cell.type.identification.step1)
GSE156337_data_Hepatocyte_list <- lapply(split(GSE156337_data_Hepatocyte_ids, GSE156337_data_Hepatocyte_source), function(x){sample(x,1000)})            
GSE156337_data_Hepatocyte_id <- unlist(GSE156337_data_Hepatocyte_list)
GSE156337_data_Hepatocyte_sml <- GSE156337_data_Hepatocyte[ ,GSE156337_data_Hepatocyte_id]

GSE156337_data_Myeloid_ids <- colnames(GSE156337_data_Myeloid)
GSE156337_data_Myeloid_source <- factor(GSE156337_data_Myeloid$cell.type.identification.step1)
GSE156337_data_Myeloid_list <- lapply(split(GSE156337_data_Myeloid_ids, GSE156337_data_Myeloid_source),function(x){sample(x,1000)})            
GSE156337_data_Myeloid_id <- unlist(GSE156337_data_Myeloid_list)
GSE156337_data_Myeloid_sml <- GSE156337_data_Myeloid[ ,GSE156337_data_Myeloid_id]

GSE156337_data_T.NK_ids <- colnames(GSE156337_data_T.NK)
GSE156337_data_T.NK_source <- factor(GSE156337_data_T.NK$cell.type.identification.step1)
GSE156337_data_T.NK_list <- lapply(split(GSE156337_data_T.NK_ids, GSE156337_data_T.NK_source),function(x){sample(x,1000)})            
GSE156337_data_T.NK_id <- unlist(GSE156337_data_T.NK_list)
GSE156337_data_T.NK_sml <- GSE156337_data_T.NK[ ,GSE156337_data_T.NK_id]

#Edit metadata for each cell type
GSE156337_data_B_sml$THY1.status <- GSE156337_data_B_sml$cell.type.identification.step1
GSE156337_data_Endothelial_sml$THY1.status <- GSE156337_data_Endothelial_sml$cell.type.identification.step1
GSE156337_data_Hepatocyte_sml$THY1.status <- GSE156337_data_Hepatocyte_sml$cell.type.identification.step1
GSE156337_data_Myeloid_sml$THY1.status <- GSE156337_data_Myeloid_sml$cell.type.identification.step1
GSE156337_data_T.NK_sml$THY1.status <- GSE156337_data_T.NK_sml$cell.type.identification.step1

#Dataset re-construction
GSE156337_data_remerge.for.cellchat <- merge(x = GSE156337_data_fibroblast,
                                             y = c(GSE156337_data_B_sml,
                                                   GSE156337_data_Endothelial_sml,
                                                   GSE156337_data_Hepatocyte_sml,
                                                   GSE156337_data_Myeloid_sml,
                                                   GSE156337_data_T.NK_sml))

#Cellular interaction analysis
cellchat <- createCellChat(object = GSE156337_data_remerge.for.cellchat, group.by = "THY1.status")
cellchat
summary(cellchat)
meta <- data.frame(cellType = GSE156337_data_remerge.for.cellchat$THY1.status, row.names =  Cells(GSE156337_data_remerge.for.cellchat))
cellchat <- addMeta(cellchat, meta = meta,  meta.name = "THY1.status")
cellchat <- setIdent(cellchat, ident.use = "THY1.status")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human          
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
cellchat <- updateCellChat(cellchat)


##---------------------------------------------------------------------GSE202642
#Data re-construction to reduce calculation 

#Cell type extraction
GSE202642_data_B <- subset(GSE202642_data, cell.type.identification.step1 %in% "B.cell")
GSE202642_data_Endothelial <- subset(GSE202642_data, cell.type.identification.step1 %in% "Endothelial.cell")
GSE202642_data_Hepatocyte <- subset(GSE202642_data, cell.type.identification.step1 %in% "Hepatocyte")
GSE202642_data_Myeloid <- subset(GSE202642_data, cell.type.identification.step1 %in% "Myeloid.cell")
GSE202642_data_T.NK <- subset(GSE202642_data, cell.type.identification.step1 %in% "T/NK.cell")

#Randomly select 1000 cells out of each cell type
GSE202642_data_B_ids <- colnames(GSE202642_data_B)
GSE202642_data_B_source <- factor(GSE202642_data_B$cell.type.identification.step1)
GSE202642_data_B_list <- lapply(split(GSE202642_data_B_ids, GSE202642_data_B_source),function(x){sample(x,1000)})            
GSE202642_data_B_id <- unlist(GSE202642_data_B_list)
GSE202642_data_B_sml <- GSE202642_data_B[ ,GSE202642_data_B_id]

GSE202642_data_Endothelial_ids <- colnames(GSE202642_data_Endothelial)
GSE202642_data_Endothelial_source <- factor(GSE202642_data_Endothelial$cell.type.identification.step1)
GSE202642_data_Endothelial_list <- lapply(split(GSE202642_data_Endothelial_ids, GSE202642_data_Endothelial_source),function(x){sample(x,1000)})            
GSE202642_data_Endothelial_id <- unlist(GSE202642_data_Endothelial_list)
GSE202642_data_Endothelial_sml <- GSE202642_data_Endothelial[ ,GSE202642_data_Endothelial_id]

GSE202642_data_Hepatocyte_ids <- colnames(GSE202642_data_Hepatocyte)
GSE202642_data_Hepatocyte_source <- factor(GSE202642_data_Hepatocyte$cell.type.identification.step1)
GSE202642_data_Hepatocyte_list <- lapply(split(GSE202642_data_Hepatocyte_ids, GSE202642_data_Hepatocyte_source), function(x){sample(x,1000)})            
GSE202642_data_Hepatocyte_id <- unlist(GSE202642_data_Hepatocyte_list)
GSE202642_data_Hepatocyte_sml <- GSE202642_data_Hepatocyte[ ,GSE202642_data_Hepatocyte_id]

GSE202642_data_Myeloid_ids <- colnames(GSE202642_data_Myeloid)
GSE202642_data_Myeloid_source <- factor(GSE202642_data_Myeloid$cell.type.identification.step1)
GSE202642_data_Myeloid_list <- lapply(split(GSE202642_data_Myeloid_ids, GSE202642_data_Myeloid_source),function(x){sample(x,1000)})            
GSE202642_data_Myeloid_id <- unlist(GSE202642_data_Myeloid_list)
GSE202642_data_Myeloid_sml <- GSE202642_data_Myeloid[ ,GSE202642_data_Myeloid_id]

GSE202642_data_T.NK_ids <- colnames(GSE202642_data_T.NK)
GSE202642_data_T.NK_source <- factor(GSE202642_data_T.NK$cell.type.identification.step1)
GSE202642_data_T.NK_list <- lapply(split(GSE202642_data_T.NK_ids, GSE202642_data_T.NK_source),function(x){sample(x,1000)})            
GSE202642_data_T.NK_id <- unlist(GSE202642_data_T.NK_list)
GSE202642_data_T.NK_sml <- GSE202642_data_T.NK[ ,GSE202642_data_T.NK_id]

#Edit metadata for each cell type
GSE202642_data_B_sml$THY1.status <- GSE202642_data_B_sml$cell.type.identification.step1
GSE202642_data_Endothelial_sml$THY1.status <- GSE202642_data_Endothelial_sml$cell.type.identification.step1
GSE202642_data_Hepatocyte_sml$THY1.status <- GSE202642_data_Hepatocyte_sml$cell.type.identification.step1
GSE202642_data_Myeloid_sml$THY1.status <- GSE202642_data_Myeloid_sml$cell.type.identification.step1
GSE202642_data_T.NK_sml$THY1.status <- GSE202642_data_T.NK_sml$cell.type.identification.step1

#Dataset re-construction
GSE202642_data_remerge.for.cellchat <- merge(x = GSE202642_data_fibroblast,
                                             y = c(GSE202642_data_B_sml,
                                                   GSE202642_data_Endothelial_sml,
                                                   GSE202642_data_Hepatocyte_sml,
                                                   GSE202642_data_Myeloid_sml,
                                                   GSE202642_data_T.NK_sml))

#Cellular interaction analysis
cellchat <- createCellChat(object = GSE202642_data_remerge.for.cellchat, group.by = "THY1.status")
cellchat
summary(cellchat)
meta <- data.frame(cellType = GSE202642_data_remerge.for.cellchat$THY1.status, row.names =  Cells(GSE202642_data_remerge.for.cellchat))
cellchat <- addMeta(cellchat, meta = meta,  meta.name = "THY1.status")
cellchat <- setIdent(cellchat, ident.use = "THY1.status")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human          
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
cellchat <- updateCellChat(cellchat)



###==============================================
###Part IV: Subcluster analysis of myeloid cells.
###==============================================

library(harmony)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(magrittr)
library(dplyr)

#Integration of datasets
Integrated.HCC.dataset <- merge(x = GSE149614_data_remerge.for.cellchat, y = c(GSE156337_data_remerge.for.cellchat, GSE202642_data_remerge.for.cellchat))

Integrated.HCC.dataset <- NormalizeData(Integrated.HCC.dataset) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

system.time({Integrated.HCC.dataset <- RunHarmony(Integrated.HCC.dataset, group.by.vars = "orig.ident", assay.use = "RNA")})

Integrated.HCC.dataset <- FindNeighbors(Integrated.HCC.dataset,reduction = "pca", dims = 1:30)

Integrated.HCC.dataset <- FindClusters(Integrated.HCC.dataset, verbose = F, resolution = 1.0)

Integrated.HCC.dataset <- RunTSNE(Integrated.HCC.dataset, reduction = "pca", dims = 1:30)

DimPlot(Integrated.HCC.dataset, reduction = "tsne", pt.size = 0.8, label = T, label.size = 6, raster=FALSE)

#Extract myeloid cell population
Integrated.HCC.dataset_myeloid <- subset(Integrated.HCC.dataset, cell.type.identification.step1 %in% c("Myeloid.cell"))

#Re-analysis of myeloid cells
Integrated.HCC.dataset_myeloid <- NormalizeData(Integrated.HCC.dataset_myeloid) %>% FindVariableFeatures() %>% ScaleData()

#NMF dimensionality reduction
library(NMF)
vm <- Integrated.HCC.dataset_myeloid@assays$RNA@scale.data
vm[vm<0] <- 0
res <- nmf(vm, 15, method = "snmf/r", seed = 'nndsvd')
nmf.embeddings <- t(coef(res))
colnames(nmf.embeddings) <- paste0("NMF_", 1:15)  #1:15指的是有15列，把这15列都按照规定改名字
nmf.embeddings <- as.matrix(nmf.embeddings)
Integrated.HCC.dataset_myeloid[["nmf"]] <- CreateDimReducObject(embeddings = nmf.embeddings, key = "NMF_")

#Clustering of myeloid cell population
Integrated.HCC.dataset_myeloid <- FindNeighbors(Integrated.HCC.dataset_myeloid,reduction = "nmf", dims = 1:15)

Integrated.HCC.dataset_myeloid <- FindClusters(Integrated.HCC.dataset_myeloid, verbose = F, resolution = 1)

Integrated.HCC.dataset_myeloid <- RunTSNE(Integrated.HCC.dataset_myeloid, reduction = "nmf", dims = 1:15)

DimPlot(Integrated.HCC.dataset_myeloid, reduction = "tsne", pt.size = 0.1, label = T, label.size = 6)

#Renaming of myeloid cell subclusters
Integrated.HCC.dataset_myeloid$myeloid.subclusters <- "0"
for(i in rownames(Integrated.HCC.dataset_myeloid@meta.data)){
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==0){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M1"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==1){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M2"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==2){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M3"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==3){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M4"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==4){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M5"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==5){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M6"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==6){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M7"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==7){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M8"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==8){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M9"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==9){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M10"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==10){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M11"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==11){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M12"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==12){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M13"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==13){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M14"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==14){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M15"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==15){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M16"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==16){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M17"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==17){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M18"};
  if(Integrated.HCC.dataset_myeloid@meta.data[i,"seurat_clusters"]==18){Integrated.HCC.dataset_myeloid@meta.data[i,"myeloid.subclusters"]<-"M19"}
}


#Marker gene identification
myeloid.markers <- FindAllMarkers(Integrated.HCC.dataset_myeloid, only.pos = TRUE)

myeloid.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

DoHeatmap(Integrated.HCC.dataset_myeloid, features = top5$gene)

#Dot plot showing the marker genes of each myeloid cell subcluster
markers.to.plot <- c('SELENOP', 'TENT2', 
                     'FCN1', 'EREG', 
                     'SPP1', 'CSTB', 
                     'APOC1','CTSC', 
                     'RACK1', 'SNHG29', 
                     'AREG', 'CLEC10A',
                     'GZMA', 'CD3E', 
                     'DONSON',
                     'FCGR3B', 'S100A9', 
                     'TIMD4', 'NDST3', 
                     'S100P', 'ADGRG3', 
                     'IDO1', 'CLEC9A', 
                     'IGFBP7', 'SPARC',
                     'APOM', 'AGXT',
                     'APOC3', 'ALB',
                     'MARCO', 'CYP2E1',
                     'MKI67', 'UBE2C',
                     'AKR1C2', 'HULC',
                     'ETS1', "IKZF3")

DotPlot(Integrated.HCC.dataset_myeloid, 
        features = markers.to.plot, 
        group.by = "myeloid.subclusters",
        dot.scale = 7,
        cols = c("grey95", "mediumblue")) + RotatedAxis()



###========================================================================
###Part V: Cellular interaction between fibroblast and myeloid subclusters.
###========================================================================
#Extract fibroblast component out of the intergrated dataset
Integrated.HCC.dataset_fibroblast <- subset(Integrated.HCC.dataset, cell.type.identification.step1 %in% c("Fibroblast"))

Integrated.HCC.dataset_fibroblast$myeloid.subclusters <- Integrated.HCC.dataset_fibroblast$THY1.status

#Construct a dataset that composed of fibroblast and myeloid sub-clusters
Merged_fibroblast_myeloid <- merge(x = Integrated.HCC.dataset_fibroblast,
                                   y = Integrated.HCC.dataset_myeloid)

#Cellular interaction between fibroblast and myeloid sub-clusters
cellchat <- createCellChat(object = Merged_fibroblast_myeloid, group.by = "myeloid.subclusters")
meta <- data.frame(cellType = Merged_fibroblast_myeloid$myeloid.subclusters, row.names =  Cells(Merged_fibroblast_myeloid))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "myeloid.subclusters")
cellchat <- setIdent(cellchat, ident.use = "myeloid.subclusters")
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human          
colnames(CellChatDB$interaction)
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

#Show the cellular interaction network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 10,height = 18)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 10,height = 18)
ht1 + ht2

#Show the ADGRE signaling pathway
netVisual_bubble(cellchat, sources.use = 20, targets.use = c(1:19), signaling = c("ADGRE"), remove.isolate = FALSE)



###============================================================
###Part VI: De-convolution of the spatial transcriptomics data.
###============================================================

library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)

##Marker gene identification from scRNA-seq data
sce <- SingleCellExperiment(assays = list(counts = Integrated.HCC.dataset@assays$RNA@counts))

sce = scater::logNormCounts(sce)

colData(sce) <- DataFrame(Integrated.HCC.dataset@meta.data)

sce <- runPCA(sce)

sce <- runUMAP(sce)

genes <- !grepl(pattern = "^Rp[l|s]|Mt", x = rownames(sce))

dec <- modelGeneVar(sce, subset.row = genes)

hvg <- getTopHVGs(dec, n = 3000)

colLabels(sce) <- colData(sce)$THY1.status

mgs <- scoreMarkers(sce, subset.row = genes)

mgs_fil <- lapply(names(mgs), function(i) {
  x <- mgs[[i]]
  x <- x[x$mean.AUC > 0.8, ]
  x <- x[order(x$mean.AUC, decreasing = TRUE), ]
  x$gene <- rownames(x)
  x$cluster <- i
  data.frame(x)
})

mgs_df <- do.call(rbind, mgs_fil)

idx <- split(seq(ncol(sce)), sce$THY1.status)

n_cells <- 500 

cs_keep <- lapply(idx, function(i) {
  n <- length(i)
  if (n < n_cells)
    n_cells <- n
  sample(i, n_cells)
})

sce <- sce[, unlist(cs_keep)]

##Deconvolution of ST data
#Take HCC-5A for example
GUJIN_ST_HCC5A <- readRDS("/GUJIN_ST_HCC5A.rds")

res1 <- SPOTlight(
  x = sce,
  y = GUJIN_ST_HCC5A@assays$SCT@counts,
  groups = as.character(sce$THY1.status),
  mgs = mgs_df,
  hvg = hvg,
  weight_id = "mean.AUC",
  group_id = "cluster",
  gene_id = "gene")

mat <- res1$mat

write.csv(mat, file = "/mat_HCC5A.csv")



###==========================================================================================================
###Part VII: Statistical analysis of the spatial co-localization pattern of THY1 with ADGRE5, TGFB1, and CCL5
###==========================================================================================================

##Spatial co-localization pattern within the same spot of ST chips
#Take HCC-1L for example
GUJIN_ST_HCC1L <- readRDS("/GUJIN_ST_HCC1L.rds")

expr_HCC1L <- as.data.frame(GUJIN_ST_HCC1L@assays$SCT@data)

mole.expr_HCC1L <- expr_HCC1L[row.names(expr_HCC1L) %in% c('THY1', 'ADGRE5', 'TGFB1', 'CCL5'),]

mole.expr_HCC1L <- as.data.frame(t(mole.expr_HCC1L))

mole.expr_HCC1L <- cbind(cellName = row.names(mole.expr_HCC1L), mole.expr_HCC1L)

#Calculate the co-expression significance
cor.test(mole.expr_HCC1L$THY1, mole.expr_HCC1L$ADGRE5, method = 'spearman')

cor.test(mole.expr_HCC1L$THY1, mole.expr_HCC1L$TGFB1, method = 'spearman')

cor.test(mole.expr_HCC1L$THY1, mole.expr_HCC1L$CCL5, method = 'spearman')


##Spatial co-localization pattern between adjacent spots of ST chips
#Define adjacent spots
findSeqPts <- function(allDotMats, oneDotMat){
  dotRow <- oneDotMat$row
  dotCol <- oneDotMat$col
  upMat <- allDotMats[allDotMats$row == dotRow-2 & allDotMats$col == dotCol,]
  dwMat <- allDotMats[allDotMats$row == dotRow+2 & allDotMats$col == dotCol,]
  rupMat <- allDotMats[allDotMats$row == dotRow-1 & allDotMats$col == dotCol+1,]
  rdwMat <- allDotMats[allDotMats$row == dotRow+1 & allDotMats$col == dotCol+1,]
  lupMat <- allDotMats[allDotMats$row == dotRow-1 & allDotMats$col == dotCol-1,]
  ldwMat <- allDotMats[allDotMats$row == dotRow+1 & allDotMats$col == dotCol-1,]
  nearMats <- rbind(upMat, dwMat, rupMat, rdwMat, lupMat, ldwMat, oneDotMat)
  seqPts <- sum(nearMats$ifSeq)
  return(seqPts)
}

geneMeanExpr <- function(allDotMats, oneDotMat, geneName){
  dotRow <- oneDotMat$row
  dotCol <- oneDotMat$col
  upName <- allDotMats[allDotMats$row == dotRow-2 & allDotMats$col == dotCol,'cellName']
  dwName <- allDotMats[allDotMats$row == dotRow+2 & allDotMats$col == dotCol,'cellName']
  rupName <- allDotMats[allDotMats$row == dotRow-1 & allDotMats$col == dotCol+1,'cellName']
  rdwName <- allDotMats[allDotMats$row == dotRow+1 & allDotMats$col == dotCol+1,'cellName']
  lupName <- allDotMats[allDotMats$row == dotRow-1 & allDotMats$col == dotCol-1,'cellName']
  ldwName <- allDotMats[allDotMats$row == dotRow+1 & allDotMats$col == dotCol-1,'cellName']
  nearMats <- allDotMats[allDotMats$cellName %in% c(upName, dwName, rupName, rdwName, lupName, ldwName, oneDotMat$cellName),]
  meanExpr <- sum(nearMats[,geneName])/oneDotMat$near
  return(meanExpr)
}

#Load spatial information
dot_pos_HCC1L <- read.csv('/tissue_positions_list.csv', col.names = c('cellName', 'ifSeq', 'row', 'col', 'unknown-1', 'unknown-2'))
dot_pos_HCC1L$near <- NA
near.expr_HCC1L <- mole.expr_HCC1L
near.expr_HCC1L$near <- NA
near.expr_HCC1L$row <- NA
near.expr_HCC1L$col <- NA

for(r in 1:nrow(dot_pos_HCC1L)){
  if(dot_pos_HCC1L[r,]$ifSeq == 1){
    dot_pos_HCC1L[r,]$near <- findSeqPts(dot_pos_HCC1L,dot_pos_HCC1L[r,])
  }
}

for(r in 1:nrow(near.expr_HCC1L)){
  seqDotName <- near.expr_HCC1L[r,'cellName']
  selectedDotPos <- dot_pos_HCC1L[dot_pos_HCC1L$cellName == seqDotName,]
  near.expr_HCC1L[r,]$near <- selectedDotPos$near
  near.expr_HCC1L[r,]$row <- selectedDotPos$row
  near.expr_HCC1L[r,]$col <- selectedDotPos$col
}

near.expr_HCC1L$THY1.near <- NA
for(r in 1:nrow(near.expr_HCC1L)){
  near.expr_HCC1L[r, 'THY1.near'] <- geneMeanExpr(near.expr_HCC1L, near.expr_HCC1L[r,], 'THY1')
}

near.expr_HCC1L$ADGRE5.near <- NA
for(r in 1:nrow(near.expr_HCC1L)){
  near.expr_HCC1L[r, 'ADGRE5.near'] <- geneMeanExpr(near.expr_HCC1L, near.expr_HCC1L[r,], 'ADGRE5')
}

near.expr_HCC1L$TGFB1.near <- NA
for(r in 1:nrow(near.expr_HCC1L)){
  near.expr_HCC1L[r, 'TGFB1.near'] <- geneMeanExpr(near.expr_HCC1L, near.expr_HCC1L[r,], 'TGFB1')
}

near.expr_HCC1L$CCL5.near <- NA
for(r in 1:nrow(near.expr_HCC1L)){
  near.expr_HCC1L[r, 'CCL5.near'] <- geneMeanExpr(near.expr_HCC1L, near.expr_HCC1L[r,], 'CCL5')
}

#Calculate the co-expression significance
cor.test(near.expr_HCC1L$THY1.near, near.expr_HCC1L$ADGRE5.near, method = 'spearman')

cor.test(near.expr_HCC1L$THY1.near, near.expr_HCC1L$TGFB1.near,method = 'spearman')

cor.test(near.expr_HCC1L$THY1.near, near.expr_HCC1L$CCL5.near,method = 'spearman')
