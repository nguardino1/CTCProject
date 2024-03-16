
##COPYKAT
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(copykat)
library(SeuratObject)


DH_P0_mBCC=readRDS("C:/Users/nguardin/Documents/RStudio/Mouse Blood/4164_P0_Blood/MouseBCC_andCTC_andNormalSkin.rds")



mBCC_exp.rawdata <- as.matrix(DH_P0_mBCC@assays$RNA@counts)
mBCC_copykat.test=copykat(rawmat=mBCC_exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name="test", distance="euclidean", norm.cell.names="",output.seg="FLASE", plot.genes="TRUE", genome="mm10",n.cores=4)
mBCC_pred.test <- data.frame(mBCC_copykat.test$prediction)
mBCC_CNA.test <- data.frame(mBCC_copykat.test$CNAmat)
write.csv(mBCC_pred.test, "/oak/stanford/groups/oro/arj58/Nick/copyKat/DH_P0_mBCC_6_29_copyKAT_predictions.csv")
write.csv(mBCC_CNA.test, "/oak/stanford/groups/oro/arj58/Nick/copyKat/DH_P0_mBCC_copyKAT_CNA.csv")


#post-copykat rds:
DH_P0_mBCC_copykat=readRDS("/Users/anna/data/R_data_and_figures/NG_inferCNV/test_copykat_clustering_results.rds")
DimPlot(DH_P0_mBCC_copykat)
DH_P0_mBCC_copykat_pred=read.table(file = "/Users/anna/data/R_data_and_figures/NG_inferCNV/test_copykat_prediction2.txt", header = T, row.names=1)
DH_P0_mBCC_copykat_manual <-Seurat::AddMetaData(DH_P0_mBCC, metadata = DH_P0_mBCC_copykat_pred)
DimPlot(DH_P0_mBCC_copykat_manual, group.by="copykat.pred")

