library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)
#Load data
WT.data<-Read10X(data.dir="DATA/sample1")
CD44h.data<-Read10X(data.dir="DATA/sample2")
#Seurat object initalization
WT <- CreateSeuratObject(counts = WT.data, project = "AT2 WT", min.cells = 3, min.features = 200)
CD44h <- CreateSeuratObject(counts = CD44h.data, project = "CD44 High", min.cells = 3, min.features = 200)

WT@meta.data[,"protocol"] <- "WT"
CD44h@meta.data[,"protocol"] <- "CD44 high"

#Brief Normalization
WT <- NormalizeData(WT)
CD44h <- NormalizeData(CD44h)

#subsetting by relative CD44 expression.
WTCD44low <- subset(x=WT,subset = Cd44 <=0.5)
CD44highsubset <- subset(x=CD44h,subset = Cd44 >= 1)
WTCD44low@meta.data[,"protocol"] <- "CD44 low"
CD44highsubset@meta.data[,"protocol"] <- "CD44 high"


#merging Dataset into one matrix
data.combine <- merge(WTCD44low,y = CD44highsubset, add.cell.ids = c("CD44 low","CD44 high"), project = "CD44 low vs CD44 high")

#Harmony
#Split objects base on protocol
data.list <- SplitObject(data.combine, split.by = "protocol")
all.genes <-rownames(data.combine)
data.list <-  lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
#Highlight features
features <- SelectIntegrationFeatures(object.list = data.list)
#Anchoring
gene.anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = features)
gene.combined <- IntegrateData(anchorset = gene.anchors,features.to.integrate = all.genes)

DefaultAssay(gene.combined) <- "integrated"
# Run the standard workflow for visualization and clustering

gene.combined <- ScaleData(gene.combined, verbose = FALSE,features = all.genes)
gene.combined <- RunPCA(gene.combined, npcs = 30, verbose = FALSE)
gene.combined <- FindNeighbors(gene.combined, reduction = "pca", dims = 1:30)
gene.combined <- FindClusters(gene.combined, resolution = 0.5)
gene.combined <- RunTSNE(gene.combined, reduction = "pca", dims = 1:30)
gene.combined <- RunUMAP(gene.combined, reduction = "pca", dims = 1:30)

#Visualization, both tsne and umap are available 
p1 <- DimPlot(gene.combined, reduction = "tsne", group.by = "protocol", pt.size = 1)
p2 <- DimPlot(gene.combined, reduction = "tsne", label = TRUE, repel = TRUE,pt.size = 1)+NoLegend()
p3 <- DimPlot(gene.combined, reduction = "tsne", split.by = "protocol", pt.size = 1)
cellcount <- table(gene.combined@meta.data$seurat_clusters, gene.combined@meta.data$orig.ident)
t1 <- tableGrob(cellcount)
p1 + ggtitle("CD44 low vs CD44 High")
p2 + ggtitle("CD44 low vs CD44 High") + t1
p3 + ggtitle("CD44 low vs CD44 High")

#Feature plot/Violin plot examples
FeaturePlot(gene.combined, features = c("Tdtomato","Zsgreen1","Sftpc", "Cdh1"), split.by = "protocol")
VlnPlot(gene.combined, features = c("Tdtomato","Zsgreen","Sftpc", "Cdh1","Thbs1", "Ifitm3"),pt.size = 0,split.by = "protocol",split.plot = TRUE)