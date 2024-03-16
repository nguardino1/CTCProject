##scRNA-seq Analysis

library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)
library(SoupX)
library(ggplot2)
library(knitr)

Blood_scRNAseq.data <- Read10X(data.dir = "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/10x_Files")
# Initialize the Seurat object with the raw (non-normalized data)
Blood4614 <- CreateSeuratObject(counts = Blood_scRNAseq.data, project = "Blood4614", min.cells = 3, min.features = 200)
Blood4614
# QC The percentage of reads that map to the mitochondrial genome.
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Blood4614[["percent.mt"]] <- PercentageFeatureSet(Blood4614, pattern = "^mt-")
head(Blood4614@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(Blood4614, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(Blood4614, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Blood4614, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Based on violin and feature plots, remove cells with mitochondria% >= 5, or #genes >= 5000
Blood4614 <- subset(Blood4614, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
#check the new distribution
VlnPlot(Blood4614, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Normalize the feature expression measurements for each cell by the total expression,
#   multiplied by a scale factor (10,000 by default), and log-transform the result.
# Normalized values are stored in pbmc[["RNA"]]@data.
Blood4614 <- NormalizeData(Blood4614, normalization.method = "LogNormalize", scale.factor = 10000)
# Calculate a subset of features that exhibit high cell-to-cell variation in the dataset
#   (i.e, they are highly expressed in some cells, and lowly expressed in others).
# Focusing on these genes in downstream analysis helps to highlight biological signals
Blood4614 <- FindVariableFeatures(Blood4614, selection.method = "vst", nfeatures = 2000)
# Linear transformation (scaling):
# 1. Shifts the expression of each gene, so that the mean expression across cells is 0
# 2. Scales the expression of each gene, so that the variance across cells is 1
all.genes <- rownames(Blood4614)
Blood4614 <- ScaleData(Blood4614, features = all.genes)
# Linear dimension reduction (PCA)
Blood4614 <- RunPCA(Blood4614, features = VariableFeatures(object = Blood4614))
ElbowPlot(Blood4614)
#identify the dimensionality of the dataset (how many PCA clusters are
#needed to explain the majority of the variance of the dataset.
#We chose 10 here, but encourage users to err on the higher side
#when choosing this parameter.
#For example, performing downstream analyses with only 5 PCs does
#signifcanltly and adversely affect results.
#use 10 PCs to find nearest neighbors during clustering
Blood4614 <- FindNeighbors(Blood4614, dims = 1:11)   #alter the dimsension when finding nearest neighbors
Blood4614 <- FindClusters(Blood4614, resolution = 0.4)
#Look at cluster IDs of the first 5 cells
head(Idents(Blood4614), 5)
# UMAP dimension reduction plot
Blood4614 <- RunUMAP(Blood4614, dims = 1:11)
DimPlot(Blood4614, reduction = "umap", label = TRUE)
#We will save this Seurat object so that it can be read in later without having to run all these steps again
saveRDS(Blood4614, file = "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/Blood4614.rds")

Blood4614.markers <- FindAllMarkers(Blood4614, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100.markers <- Blood4614.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)
top100.markers
write.csv(top100.markers, "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/Blood4614_top100_markers.csv")

Blood4614_nonplatelet=subset(Blood4614, idents=c(3, 4, 6, 7, 9, 10))
DimPlot(Blood4614_nonplatelet, label = TRUE, pt.size=0.5)


#Subsetted without platelets
Blood4614_nonplateletnonRBC=subset(Blood4614_nonplatelet, idents=c(2, 4, 6, 7, 5, 0))

# Normalized values are stored in pbmc[["RNA"]]@data.
Blood4614_nonplatelet <- NormalizeData(Blood4614_nonplatelet, normalization.method = "LogNormalize", scale.factor = 10000)
# Calculate a subset of features that exhibit high cell-to-cell variation in the dataset
#   (i.e, they are highly expressed in some cells, and lowly expressed in others).
# Focusing on these genes in downstream analysis helps to highlight biological signals
Blood4614_nonplatelet <- FindVariableFeatures(Blood4614_nonplatelet, selection.method = "vst", nfeatures = 2000)
# Linear transformation (scaling):
# 1. Shifts the expression of each gene, so that the mean expression across cells is 0
# 2. Scales the expression of each gene, so that the variance across cells is 1
all.genes <- rownames(Blood4614_nonplatelet)
Blood4614_nonplatelet <- ScaleData(Blood4614_nonplatelet, features = all.genes)
# Linear dimension reduction (PCA)
Blood4614_nonplatelet <- RunPCA(Blood4614_nonplatelet, features = VariableFeatures(object = Blood4614_nonplatelet))
ElbowPlot(Blood4614_nonplatelet)
#identify the dimensionality of the dataset (how many PCA clusters are
#needed to explain the majority of the variance of the dataset.
#We chose 10 here, but encourage users to err on the higher side
#when choosing this parameter.
#For example, performing downstream analyses with only 5 PCs does
#signifcanltly and adversely affect results.
#use 10 PCs to find nearest neighbors during clustering
Blood4614_nonplatelet <- FindNeighbors(Blood4614_nonplatelet, dims = 1:10)   #alter the dimsension when finding nearest neighbors
Blood4614_nonplatelet <- FindClusters(Blood4614_nonplatelet, resolution = 0.3)
#Look at cluster IDs of the first 5 cells
head(Idents(Blood4614_nonplatelet), 5)
# UMAP dimension reduction plot
Blood4614_nonplatelet <- RunUMAP(Blood4614_nonplatelet, dims = 1:10)
#We will save this Seurat object so that it can be read in later without having to run all these steps again
saveRDS(Blood4614_nonplateletnonRBC, file = "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/Blood4614_nonplateletnonRBC.rds")
Atlas1=readRDS("/Users/nguardin/Documents/RStudio/Intestine/RDSFiles/Atlas1.Rds")
FeaturePlot(Atlas1, reduction = "umap", features = c("Cdx2", "Cdx1"), sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

Blood4614_nonplatelet.markers <- FindAllMarkers(Blood4614_nonplatelet, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100.markers <- Blood4614_nonplatelet.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)
top100.markers
write.csv(top100.markers, "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/Blood4614_nonplatelet_top100_markers.csv")

#Subsetted without RBCs
Blood4614_nonplateletnonRBC=subset(Blood4614_nonplatelet, idents=c(2, 4, 6, 7, 5, 0))
DimPlot(Blood4614_nonplateletnonRBC, label = TRUE, pt.size=0.5)
# Normalized values are stored in pbmc[["RNA"]]@data.
Blood4614_nonplateletnonRBC <- NormalizeData(Blood4614_nonplateletnonRBC, normalization.method = "LogNormalize", scale.factor = 10000)
# Calculate a subset of features that exhibit high cell-to-cell variation in the dataset
#   (i.e, they are highly expressed in some cells, and lowly expressed in others).
# Focusing on these genes in downstream analysis helps to highlight biological signals
Blood4614_nonplateletnonRBC <- FindVariableFeatures(Blood4614_nonplateletnonRBC, selection.method = "vst", nfeatures = 2000)
# Linear transformation (scaling):
# 1. Shifts the expression of each gene, so that the mean expression across cells is 0
# 2. Scales the expression of each gene, so that the variance across cells is 1
all.genes <- rownames(Blood4614_nonplateletnonRBC)
Blood4614_nonplateletnonRBC <- ScaleData(Blood4614_nonplateletnonRBC, features = all.genes)
# Linear dimension reduction (PCA)
Blood4614_nonplateletnonRBC <- RunPCA(Blood4614_nonplateletnonRBC, features = VariableFeatures(object = Blood4614_nonplateletnonRBC))
ElbowPlot(Blood4614_nonplateletnonRBC)
#identify the dimensionality of the dataset (how many PCA clusters are
#needed to explain the majority of the variance of the dataset.
#We chose 10 here, but encourage users to err on the higher side
#when choosing this parameter.
#For example, performing downstream analyses with only 5 PCs does
#signifcanltly and adversely affect results.
#use 10 PCs to find nearest neighbors during clustering
Blood4614_nonplateletnonRBC <- FindNeighbors(Blood4614_nonplateletnonRBC, dims = 1:10)   #alter the dimsension when finding nearest neighbors
Blood4614_nonplateletnonRBC <- FindClusters(Blood4614_nonplateletnonRBC, resolution = 0.3)
#Look at cluster IDs of the first 5 cells
head(Idents(Blood4614_nonplateletnonRBC), 5)
# UMAP dimension reduction plot
Blood4614_nonplateletnonRBC <- RunUMAP(Blood4614_nonplateletnonRBC, dims = 1:10)
DimPlot(Blood4614_nonplateletnonRBC, reduction = "umap", label = TRUE, pt.size = 1.5)
#We will save this Seurat object so that it can be read in later without having to run all these steps again
saveRDS(Blood4614_nonplatelet, file = "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/Blood4614_nonplateletnonRBC.rds")
Atlas1=readRDS("/Users/nguardin/Documents/RStudio/Intestine/RDSFiles/Atlas1.Rds")
FeaturePlot(Atlas1, reduction = "umap", features = c("Cdx2", "Cdx1"), sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

Blood4614_nonplateletnonRBC.markers <- FindAllMarkers(Blood4614_nonplateletnonRBC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100.markers <- Blood4614_nonplateletnonRBC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)
top100.markers
write.csv(top100.markers, "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/Blood4614_nonplateletnonRBC_top100_markers.csv")

FeaturePlot(Blood4614_nonplateletnonRBC, reduction = "umap", features = c("Pecam1", "Krt19","Itgam", "Krt14", "Ptprc", "Vim", "Epcam", "Gp1ba", "Gypa"), sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)
FeaturePlot(Blood4614_nonplateletnonRBC, reduction = "umap", features = c("Krt19", "Ptprc","Cd34", "Col1a2", "Ly6d"), sort.cell = TRUE, min.cutoff = 'q10', label = TRUE)

MouseBCC=readRDS("C:/Users/nguardin/Documents/RStudio/MouseBCC/RDSFiles/P0andP1mBCC.Rds")
FeaturePlot(MouseBCC, reduction = "umap", features = c("Krt19", "Dcn","Pdpn", "Pdgfra", "Thy1", "Cd34", "Epcam", "Gp1ba", "Gypa"))


new.cluster.ids <- c("Immune 1", "Immune 2", "CTCs", "Immune 3", "Immune 4", "Immune 5", "Immune 6")
names(new.cluster.ids) <- levels(Blood4614_nonplateletnonRBC)
Blood4614_nonplateletnonRBC <- RenameIdents(Blood4614_nonplateletnonRBC, new.cluster.ids)
DimPlot(Blood4614_nonplateletnonRBC, reduction = "umap", label = TRUE, pt.size = 1.0)    


##Integrating DH_P0_mBCC with CTC only##
DH_P0_mBCC_Epithelial_Fibroblast=subset(DH_P0_mBCC, idents=c(0,1,2,3,5,6,7,8,9,10))
#Anchoring data, this will take some time
samp.anchors <- FindIntegrationAnchors(object.list = list(DH_P0_mBCC, Blood4614_nonplateletnonRBC), dims = 1:20)
normal_BCC_andbloodcombined <- IntegrateData(anchorset = samp.anchors, dims = 1:20)
#At this point you must switch the DefaultAssay to "integrated", note that you will have to change the DefaultAssay later, this is an important thing to keep track of!
DefaultAssay(normal_BCC_andbloodcombined) <- "integrated"
#Generate your UMAP, note that the number of dims used and the resolution will depend on
normal_BCC_andbloodcombined <- FindVariableFeatures(normal_BCC_andbloodcombined, selection.method = "vst", nfeatures = 2000)
normal_BCC_andbloodcombined <- ScaleData(normal_BCC_andbloodcombined, verbose = FALSE)
normal_BCC_andbloodcombined <- RunPCA(normal_BCC_andbloodcombined, npcs = 30, verbose = FALSE)
ElbowPlot(normal_BCC_andbloodcombined, ndims=30)
normal_BCC_andbloodcombined <- RunUMAP(normal_BCC_andbloodcombined, reduction = "pca", dims = 1:16)
normal_BCC_andbloodcombined <- FindNeighbors(normal_BCC_andbloodcombined, reduction = "pca", dims = 1:16)
normal_BCC_andbloodcombined <- FindClusters(normal_BCC_andbloodcombined, resolution = 0.3)
DimPlot(normal_BCC_andbloodcombined, reduction = "umap", group.by = "orig.ident", pt.size=1.5)
DimPlot(normal_BCC_andbloodcombined, reduction = "umap", label = TRUE, pt.size=1.5)
DefaultAssay(normal_BCC_andbloodcombined) <- "RNA"
#general markers:
FeaturePlot(normal_BCC_andCTCs.combined, features = c("Gli1", "CD34", "Dcn", "Krt19", "Pdgfra", "Ptprc", "Pecam1", "Pmel", "Cd3d"), cols=c("gray", "red"))

saveRDS(normal_BCC_andbloodcombined, file = "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/normal_BCC_andbloodcombinbed.rds")


normal_BCC_andblood.combined <- FindAllMarkers(normal_BCC_andblood.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top100.markers <- normal_BCC_andblood.combined %>%
  group_by(cluster) %>%
  slice_max(n = 100, order_by = avg_log2FC)
top100.markers
write.csv(top100.markers, "C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/normal_BCC_andBloodintegrated1_top100_markers.csv")

new.cluster.ids <- c("Tumor Epithelial 1", "Tumor Epithelial 2", "Fibroblast 1", "Fibroblast 2", "Immune 1", "Tumor Epithelial 3", "Tumor Epithelial 4",
                     "Immune 2", "Tumor Epithelial 5", "Tumor EMT", "Endothelial", "Immune 3", "Muscle", "RBC", "Fibroblast 3")
names(new.cluster.ids) <- levels(normal_BCC_andbloodcombined)
normal_BCC_andbloodcombined <- RenameIdents(normal_BCC_andbloodcombined, new.cluster.ids)
DimPlot(normal_BCC_andbloodcombined, reduction = "umap", label = TRUE, pt.size = 1.0)

