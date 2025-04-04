library(Seurat)
library(clustree)
library(ggraph)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(cowplot)



N1A.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
N1A <- CreateSeuratObject(counts = N1A.data, project = "N1A", min.cells = 3, min.features = 200)
N1A@meta.data$human<-"Normal"
N1A@meta.data$sample<-"N1A"




N2A.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
N2A <- CreateSeuratObject(counts = N2A.data, project = "N2A", min.cells = 3, min.features = 200)
N2A@meta.data$human<-"Normal"
N2A@meta.data$sample<-"N2A"





N3A.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
N3A <- CreateSeuratObject(counts = N3A.data, project = "N3A", min.cells = 3, min.features = 200)
N3A@meta.data$human<-"Normal"
N3A@meta.data$sample<-"N3A"






N4A.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
N4A <- CreateSeuratObject(counts = N4A.data, project = "N4A", min.cells = 3, min.features = 200)
N4A@meta.data$human<-"Normal"
N4A@meta.data$sample<-"N4A"







M1.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
M1 <- CreateSeuratObject(counts = M1.data, project = "M1", min.cells = 3, min.features = 200)
M1@meta.data$human<-"Myxomatous"
M1@meta.data$sample<-"M1"


M2.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
M2 <- CreateSeuratObject(counts = M2.data, project = "M2", min.cells = 3, min.features = 200)
M2@meta.data$human<-"Myxomatous"
M2@meta.data$sample<-"M2"


M3.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
M3 <- CreateSeuratObject(counts = M3.data, project = "M3", min.cells = 3, min.features = 200)
M3@meta.data$human<-"Myxomatous"
M3@meta.data$sample<-"M3"


M4.data<- Read10X(data.dir = "filtered_feature_bc_matrix/")
M4 <- CreateSeuratObject(counts = M4.data, project = "M4", min.cells = 3, min.features = 200)
M4@meta.data$human<-"Myxomatous"
M4@meta.data$sample<-"M4"









##merge all samples，name it as MV，and preattach on each cells
MV <- merge(N1A, y = c(N1P,N2A,N3A,N4A,M1,M2,M3,M4), add.cell.ids = c("N1A","N2A","N3A","N4A","M1","M2","M3","M4"), project = "MV_all")



###QC
levels(MV)
MV[["percent.mt"]] <- PercentageFeatureSet(MV, pattern = "^MT-")


VlnPlot(MV, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MV2 <- subset(MV, subset = nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt <10)


###NORMALIZE
MV2_1 <- NormalizeData(MV2, normalization.method = "LogNormalize", scale.factor = 10000)

MV2_1 <- FindVariableFeatures(MV2_1, selection.method = "vst", nfeatures = 2000)




MV2_1 <- ScaleData(MV2_1)

MV2_1 <- RunPCA(MV2_1, features = VariableFeatures(object = MV2_1))


ElbowPlot(MV2_1,ndims=50)


###Harmony

MV = MV %>% RunHarmony("sample", plot_convergence = TRUE)


MV <- FindNeighbors(MV, reduction = "harmony", dims = 1:30)


MV <- FindClusters(MV, resolution = c(seq(0,1,0.1)))


MV <- RunUMAP(MV, reduction = "harmony", dims = 1:30)











###UMAP 
MV2_1 <- RunUMAP(MV2_1, dims = 1:30)



pdf(file = "UMAP.pdf", width=5.3, height=4)
DimPlot(MV2_1, reduction = "umap")
dev.off()





