library(Seurat)
library(patchwork)
library(ggplot2)

col_p = c("#F3E5C2","#EFB89A","#EF8F50", "#F5731D","#FA4602","#FE5100")
#WT.data<-Read10X(data.dir="WT")
#CD44h.data<-Read10X(data.dir="CD44hi")
#create Seurat objects
WT <- CreateSeuratObject(counts = WT.data, project = "AT2 WT", min.cells = 3, min.features = 200)
CD44h <- CreateSeuratObject(counts = CD44h.data, project = "CD44 High", min.cells = 3, min.features = 200)
WT[["percent.mt"]] <- PercentageFeatureSet(WT, pattern = "^mt-")
CD44h[["percent.mt"]] <- PercentageFeatureSet(CD44h, pattern = "^mt-")

#Qc
#VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(CD44h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#plot1 <- FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2


CD44h@meta.data[,"protocol"] <- "CD44 high"
WT@meta.data[,"protocol"] <- "WT"
CD44 <- merge(CD44h,y = WT, add.cell.ids = c("CD44hi","WT"), project = "CD44hi vs WT")
CD44<- subset(CD44, subset = nFeature_RNA > 2000 & nFeature_RNA < 5500 & nCount_RNA < 30000 & percent.mt <7.5)

s.genes<-cc.genes.updated.2019$s.genes
g2m.genes<-cc.genes.updated.2019$g2m.genes
m.s.genes<-str_to_title(s.genes)
m.g2m.genes<-str_to_title(g2m.genes)
CD44<-CellCycleScoring(CD44, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = F)
CD44$CC.Difference<-CD44$S.Score-CD44$G2M.Score
CD44<-NormalizeData(CD44)
CD44<-ScaleData(CD44, verbose = T, vars.to.regress = "CC.Difference", features = rownames(CD44))
saveRDS(CD44, "CD44_ccscale.rds")
CD44<-readRDS("CD44_ccscale.rds")

CD44<-FindVariableFeatures(CD44, selection.method = "vst", nfeatures = 5000)
CD44<-RunPCA(CD44, npcs = 30, verbose = F)
ElbowPlot(CD44)
CD44<-FindNeighbors(CD44, reduction = "pca", dims = 1:10)
CD44<-FindClusters(CD44, resolution = 0.1)
CD44<-RunUMAP(CD44, reduction = "pca", dims = 1:10, min.dist = 0.5, n.neighbors = 40, spread = 3, local.connectivity = 18)
DimPlot(CD44, reduction = "umap", label = T)
FeaturePlot(CD44, features = c("Foxj1","Ptprc","Mki67", "Sox2", "Muc5b", "Scgb1a1", "Hopx","Pdpn", "Tdtomato"),cols=col_p, min.cutoff = "q10", order = T)
saveRDS(CD44, "CD44.rds")
CD44<- readRDS("CD44.rds")

CD44_2 <- RenameIdents(CD44, '0'="AEpiC",'1'="Ciliated", '2'="Proliferating",'3'="Immune" )     
AEpiC<-subset(CD44_2, idents = c('AEpiC', 'Proliferating'))
AEpiC<-NormalizeData(AEpiC, scale.factor = 10000)
AEpiC<-ScaleData(AEpiC, vars.to.regress = "CC.Difference", features = rownames(AEpiC))
AEpiC<-FindVariableFeatures(AEpiC, selection.method = "vst", nfeatures = 5000)
AEpiC<-RunPCA(AEpiC)
ElbowPlot(AEpiC)
AEpiC<-FindNeighbors(AEpiC, dims = 1:15)
AEpiC<-FindClusters(AEpiC, resolution = 0.8)
AEpiC<-RunUMAP(AEpiC, dims = 1:15, min.dist = 1, n.neighbors = 30, spread = 3, local.connectivity = 18)
DimPlot(AEpiC, reduction = "umap", label = T)
DimPlot(AEpiC, reduction = "umap", group.by = "protocol")
FeaturePlot(AEpiC, features = c("Foxj1","Ptprc","Mki67", "Sox2", "Muc5b", "Scgb1a1", "Hopx","Pdpn", "Tdtomato", "Cd44"),cols=col_p, min.cutoff = "q10", order = T)

saveRDS(AEpiC, "AEpiCccdifference.rds")
AEpiC <- readRDS("AEpiCccdifference.rds")
AT2<-subset(AEpiC, idents = c('0','1','2','3','4','5','6','7','8','9','11'))
AT2<-NormalizeData(AT2)
AT2<-ScaleData(AT2, vars.to.regress = "CC.Difference", features = rownames(AT2))
AT2<-FindVariableFeatures(AT2)
AT2<-RunPCA(AT2)
ElbowPlot(AT2, ndims = 30)
AT2<-FindNeighbors(AT2, dims = 1:15)
AT2<-FindClusters(AT2, resolution = 0.3)
AT2<-RunUMAP(AT2, dims = 1:15, min.dist = 1, n.neighbors = 40, spread = 3, local.connectivity = 20)
DimPlot(AT2, reduction = "umap", label = T)
DimPlot(AT2, reduction = "umap", split.by = "protocol")
FeaturePlot(AT2, features = c("Foxj1","Ptprc","Mki67", "Sox2", "Muc5b", "Scgb1a1", "Hopx","Pdpn", "Tdtomato", "Cd44"), min.cutoff = "q10", order = T)
FeaturePlot(AT2, features = c("Sftpc", "Etv5"), min.cutoff = "q10", order = T)
FeaturePlot(AT2, features = c("Tlr4"),cols=col_p, min.cutoff = "q10", order = T,split.by = "protocol")
VlnPlot(AT2, features = "Cd44", pt.size = 0)

DotPlot(AT2, features = c("Cd44","Thbs1", "Nfkb2","Tnfrsf9","Ccl20","Mki67","Top2a","Krt19","Tlr2","Tlr3","Tlr4"))
DotPlot(AT2, features = c("Thbs1", "Nfkb2","Tnfrsf9","Ccl20","Mki67","Top2a"))
DotPlot(AT2, features = c("Krt19", "Cd44","Tlr2","Tlr3","Tlr4","Top2a"))
FeaturePlot(AT2, features = "Cd44",cols = col_p)
DotPlot(AT2, features = "Cd44")
cluster4_6.markers <- FindMarkers(AT2, ident.1 = 4, ident.2 = 6, min.pct = 0.25)
cluster4.markers<-FindMarkers(AT2, ident.1 = 4, min.pct = 0.25)

AT2_2 <- RenameIdents(CD44, '0'="WT AT2",'1'="WT AT2", '2'="CD44h AT2",'3'="CD44h AT2",'4'="Proliferating 1",'5'="CD44h AT2",'6'="Proliferating 2") 
VlnPlot(AT2_2, features = c("Thbs1", "Nfkb2","Tnfrsf9","Ccl20","Mki67","Top2a"))
VlnPlot(AT2_2, features = c("Thbs1", "Nfkb2","Tnfrsf9"), pt.size = 0)
VlnPlot(AT2_2,features = c("Krt8", "Lgals3","Krt19", "Sfn") )
saveRDS(AT2, "AT2.rds")












