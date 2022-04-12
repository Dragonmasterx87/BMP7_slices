# LOAD LIBRARIES ####
# Restart Rstudio or R
# control (3) vs BMP-7 (3) (stage 2)
install.packages("devtools")
install.packages("usethis")
library(usethis)
library(devtools)
install.packages("pkgbuild")
library(pkgbuild)
install.packages("Matrix")
install.packages("ggridges")
install.packages("cowplot")
install.packages('ggrepel')
install.packages("R.utils")
install.packages("gridExtra")
install.packages("Seurat")
install.packages("plotly")
install.packages("clustree")
install.packages('multtest')
install.packages("ggplot2")
install.packages("ggraph")
install.packages('Seurat')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")

# Install multtest for Seurat
BiocManager::install("multtest")
library(ggplot2)
library(cowplot)
library(Matrix)
library(ggridges)
library(ggrepel)
library(dplyr)
library(Seurat)
library(monocle)
library(plotly)
library(ggraph)
library(clustree)
library(Seurat)
library(ggraph)
library(Matrix)
library(VGAM)
library(stats4)
library(splines)
library(DDRTree)
library(irlba)
library(BiocGenerics)
library(Biobase)
library(monocle)
packageVersion("Seurat")
packageVersion("monocle")
library(Seurat)
library(cowplot)
library(patchwork)
library(SeuratObject)
library(Seurat)
library(sctransform)
library(dplyr)
library(RColorBrewer)
library(ggthemes)
library(ggplot2)
library(cowplot)
library(data.table)

ctrl1st2.data <-Read10X(data.dir = "D:/30-649959914/01_analysis/cellranger_count/ST2-CTRL/filtered_feature_bc_matrix")

ctrl2st2.data <-Read10X(data.dir = "D:/30-649959914/01_analysis/cellranger_count/ST-2-CTRL-1/filtered_feature_bc_matrix")

                               
ctrl3st2.data <-Read10X(data.dir = "D:/30-649959914/01_analysis/cellranger_count/ST-2-CTRL-2/filtered_feature_bc_matrix")
                         
control.list = list() # First create an empty list to hold the Seurat objects
control.list[[1]] = CreateSeuratObject(counts = ctrl1st2.data, 
                                       project = "ctrl1st2"
)

control.list[[2]] = CreateSeuratObject(counts = ctrl2st2.data, 
                                       project = "ctrl2st2")

control.list[[3]] = CreateSeuratObject(counts = ctrl3st2.data, 
                                       project = "ctrl3st2"
)


controldiabetes <- merge(x=control.list[[1]], y=c(control.list[[2]],control.list[[3]]), add.cell.ids = c("ctrl1st2","ctrl2st2","ctrl3st2"), project="CSHL")

controldiabetes[["percent.mt"]] <- PercentageFeatureSet(object = controldiabetes, pattern = "^MT-")
                               
 p1 <- VlnPlot(object = controldiabetes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
                            
 plot(p1)
                              
 p2 <- FeatureScatter(object = controldiabetes, feature1 = "nCount_RNA", feature2 = "percent.mt")
                               
plot(p2)
controldiabetes <- subset(x =  controldiabetes, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
 p3 <- VlnPlot(object = controldiabetes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
                               
plot(p3)

p4 <- FeatureScatter(object = controldiabetes, feature1 = "nCount_RNA", feature2 = "percent.mt")
                              
plot(p4)                               


BMP1st2.data <-Read10X(data.dir = "D:/30-649959914/01_analysis/cellranger_count/ST2-BMP-7/filtered_feature_bc_matrix")
BMP2st2.data <-Read10X(data.dir = "D:/30-649959914/01_analysis/cellranger_count/ST-2-BMP-7/filtered_feature_bc_matrix")
BMP3st2.data <-Read10X(data.dir = "D:/30-649959914/01_analysis/cellranger_count/ST--2-bmp-7/filtered_feature_bc_matrix")

BMP.list = list() # First create an empty list to hold the Seurat objects
BMP.list[[1]] = CreateSeuratObject(counts = BMP1st2.data, 
                                   project = "BMP1st2"
)

BMP.list[[2]] = CreateSeuratObject(counts = BMP2st2.data, 
                                   project = "BMP2st2")

BMP.list[[3]] = CreateSeuratObject(counts = BMP3st2.data, 
                                   project = "BMP3st2"
)


BMPdiabetes <- merge(x=BMP.list[[1]], y=c(BMP.list[[2]],BMP.list[[3]]), add.cell.ids = c("BMP1st2","BMP2st2","BMP3st2"), project="CSHL")

BMPdiabetes[["percent.mt"]] <- PercentageFeatureSet(object = BMPdiabetes, pattern = "^MT-")

p5 <- VlnPlot(object = BMPdiabetes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot(p5)

p6 <- FeatureScatter(object = BMPdiabetes, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot(p6)
BMPdiabetes <- subset(x =  BMPdiabetes, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
p7 <- VlnPlot(object = BMPdiabetes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot(p7)

p8 <- FeatureScatter(object = BMPdiabetes, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot(p8)                               

options(future.globals.maxSize= 891289600000000)
memory.limit()
memory.limit(size=21474836480)


controldiabetes$BMPdiabetes <- "CONTROL" # Assign new cell annotations to a new "identity class" in the meta data
BMPdiabetes$BMPdiabetes <- "BMP7"
# Make the list of 2 Seurat objects
diabetesst2comp <- c(controldiabetes,BMPdiabetes)

diabetesst2comp <- lapply(X = diabetesst2comp, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = diabetesst2comp)
diabetesst2comp <- lapply(X = diabetesst2comp, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

Diabetes.anchors <- FindIntegrationAnchors(object.list = diabetesst2comp, anchor.features = features, reduction = "rpca")
Diabetes.combined <- IntegrateData(anchorset = Diabetes.anchors)
DefaultAssay(Diabetes.combined) <- "integrated"
Diabetes.combined <- ScaleData(Diabetes.combined, verbose = FALSE)
Diabetes.combined <- RunPCA(Diabetes.combined, npcs = 30, verbose = FALSE)
Diabetes.combined <- RunUMAP(Diabetes.combined, reduction = "pca", dims = 1:30)
Diabetes.combined <- FindNeighbors(Diabetes.combined, reduction = "pca", dims = 1:30)
Diabetes.combined <- FindClusters(Diabetes.combined, resolution = 0.5)

p11 <- DimPlot(Diabetes.combined, reduction = "umap", raster=FALSE)
p12 <- DimPlot(Diabetes.combined, reduction = "umap",  label = TRUE,
              repel = TRUE, raster=FALSE)

p11 + p12

install.packages("ctv")
install.packages("iterators")
install.packages("doParallel")

ibrary(parallel)
library(iterators)
library(foreach)
library(doParallel)

cl = makeCluster(split)
init = clusterEvalQ(cl, { library(MASS); NULL })
results = parLapplyLB(cl
                      ,rep(eachStart, split)
                      ,function(nstart) kmeans(Boston, 4, nstart=nstart))
withinss = sapply(results, function(result) result$tot.withinss)
result = results[[which.min(withinss)]]
stopCluster(cl)
result$tot.withinss
split = detectCores()
eachStart = 25
# set up iterators
iters = iter(rep(eachStart, split))
# set up combine function
comb = function(res1, res2) {
  if(res1$tot.withinss < res2$tot.withinss) res1 else res2
}

cl = makeCluster(split)
registerDoParallel(cl)
result = foreach(nstart=iters, .combine="comb", .packages="MASS") %dopar%
  kmeans(Boston, 4, nstart=nstart)
stopCluster(cl)
result$tot.withinss
                  
DefaultAssay(Diabetes.combined) <- "RNA"

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

feature.plot.ad <- c("SPP1", "CFTR",
                     "AQP1", "ALDH1A3",
                     "KRT19", "CRP",
                     "DEFB1", "CEACAM6",
                     "MMP7", "TSPAN8",
                     "ONECUT2", "LITAF",
                     "SOX4", "DAB2",
                     "CREB5", "HLA-DQB1",
                     "WWTR1", "PPARGC1A",
                     "PKHD1", "NFIB",
                     "PNLIP", "REG1B",
                     "PRSS1", "ALB",
                     "CPA2", "CTRB2",
                     "CEL", "PLA2G1B",
                     "CELA3A", "GATA4",
                     "MECOM", "NR5A2",
                     "ZFP36L1", "CEBPD",
                     "CREB3L1", "XBP1",
                     "LGR4", "NUPR1")
DotPlot(object = Diabetes.combined, 
        features = feature.plot.ad, 
        dot.scale = 10, 
        col.min = 0,
        col.max = 2,
) + scale_colour_gradient2(low = "#0200ad", mid = "#ffe272", high = "#ff0000")
FeaturePlot(Diabetes.combined, 
            features = c("ID1", "ID2"), 
            blend = TRUE, 
            pt.size = 1, 
            order = TRUE, 
            blend.threshold = 0.1, 
            label.size = 6, 
            sort.cell = TRUE, 
)



DefaultAssay(object = Diabetes.combined)
DefaultAssay(object = Diabetes.combined) <- "RNA"
Diabetes.combined <- NormalizeData(Diabetes.combined)
Diabetes.combined <- SCTransform(Diabetes.combined, assay = "RNA", new.assay.name = "SCT", verbose = TRUE, return.only.var.genes = TRUE)
DefaultAssay(object = Diabetes.combined)
DefaultAssay(object = Diabetes.combined) <- "SCT"
top10.genes.vst <- head(x = VariableFeatures(object = Diabetes.combined), 10)
p106 <- VariableFeaturePlot(object = Diabetes.combined, assay = 'SCT', selection.method = c('sct'))
p106
LabelPoints(plot = p106, points = top10.genes.vst, repel = TRUE)


head(x = colnames(x = Diabetes.combined))
tail(x = colnames(x = Diabetes.combined))

unique(x = sapply(X = strsplit(x = colnames(x = Diabetes.combined), split = "_"), FUN = "[", 1))
table(Diabetes.combined$orig.ident)
# LINEAR DIMENSIONALITY REDUCTION ####
DefaultAssay(object = Diabetes.combined)
DefaultAssay(object = Diabetes.combined) <- "integrated"
Diabetes.combined <- RunPCA(object = Diabetes.combined, features = VariableFeatures(object = Diabetes.combined))
print(x = Diabetes.combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(object = Diabetes.combined, dims = 1:2, reduction = "pca")
DimPlot(object = Diabetes.combined, reduction = "pca")
DimHeatmap(object = Diabetes.combined, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = Diabetes.combined, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(object = Diabetes.combined)
Diabetes.combined <- FindNeighbors(object = Diabetes.combined, dims = 1:20)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.1)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.2)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.3)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.4)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.5)
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.6)
clustree(Diabetes.combined, prefix = "integrated_snn_res.")
Diabetes.combined <- FindClusters(object = Diabetes.combined, resolution = 0.4)
table(Diabetes.combined$seurat_clusters)
table(Idents(Diabetes.combined), Diabetes.combined$orig.ident)
Diabetes.combined <- RunUMAP(Diabetes.combined, dims = 1:30)
p1010 <- DimPlot(object = Diabetes.combined, group.by = c("orig.ident", "integrated_snn_res.0.4"), combine = FALSE, pt.size = 1, raster = FALSE)
p1010 <- lapply(X = p1010, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 3))))
CombinePlots(p1010)

FeaturePlot(Diabetes.combined, features = c("ID1", "ID2", "ID3", "ID4"), split.by = "BMPdiabetes", max.cutoff = 4,
            cols = c("grey", "red"))
FeaturePlot(Diabetes.combined, features = c("MALAT1", "PLCG2", "REG1B", "ID3", "PRSS2"), split.by = "BMPdiabetes", min.cutoff = "q5")
plots <- VlnPlot(Diabetes.combined, features = c("ID1", "ID2", "ID3", "ID4"), split.by = "BMPdiabetes", 
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)    

                
DefaultAssay(object = Diabetes.combined)
DefaultAssay(object = Diabetes.combined) <- "SCT"
Diabetes.markers <- FindAllMarkers(object = Diabetes.combined, 
                                   features = VariableFeatures(Diabetes.combined, assay = 'SCT'), 
                                   only.pos = TRUE, 
                                   min.pct = 0.1, 
                                   logfc.threshold = 0.41, 
                                   assay = 'SCT',
                                   slot = c('data'))
library(magrittr)
library(dplyr)
Diabetes.markers %>% group_by(cluster) %>% top_n(n = 2)

write.csv(Diabetes.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/db.markerscombined.csv')


DefaultAssay(object = Diabetes.combined) <- "SCT"
top30.nomes <- Diabetes.markers %>% group_by(cluster) %>% top_n(n = 15)
DoHeatmap(object = Diabetes.combined, 
          features = top30.nomes$gene, 
          disp.min = -2.5, 
          disp.max = 2.5,
          label = FALSE) + scale_fill_gradientn(colors = colorRampPalette(c("#0200ad", 
                                                                            "#fbfcbd", 
                                                                            "#ff0000"))(256))
-----------------------------------------------------------------------------------------------------------------------------
  
  # FOR DE ANALYSIS- 
DefaultAssay(object = Diabetes.combined) <- "RNA"  
#object_RNA <- NormalizeData(object, verbose = FALSE)
Diabetes.combined <- NormalizeData(Diabetes.combined)
cluster2.markers <- FindMarkers(Diabetes.combined, ident.1 = 2, min.pct = 0.25)
write.csv(cluster2.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster2.csv')

cluster0.markers <- FindMarkers(Diabetes.combined, ident.1 = 0, min.pct = 0.25)
write.csv(cluster0.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster0.csv')

cluster1.markers <- FindMarkers(Diabetes.combined, ident.1 = 1, min.pct = 0.25)
write.csv(cluster1.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster1.csv')

cluster3.markers <- FindMarkers(Diabetes.combined, ident.1 = 3, min.pct = 0.25)
write.csv(cluster3.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster3.csv')

cluster4.markers <- FindMarkers(Diabetes.combined, ident.1 = 4, min.pct = 0.25)
write.csv(cluster4.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster4.csv')

cluster5.markers <- FindMarkers(Diabetes.combined, ident.1 = 5, min.pct = 0.25)
write.csv(cluster5.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster5.csv')

cluster6.markers <- FindMarkers(Diabetes.combined, ident.1 = 6, min.pct = 0.25)
write.csv(cluster6.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster6.csv')

cluster7.markers <- FindMarkers(Diabetes.combined, ident.1 = 7, min.pct = 0.25)
write.csv(cluster7.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster7.csv')

cluster8.markers <- FindMarkers(Diabetes.combined, ident.1 = 8, min.pct = 0.25)
write.csv(cluster8.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster8.csv')

cluster9.markers <- FindMarkers(Diabetes.combined, ident.1 = 9, min.pct = 0.25)
write.csv(cluster9.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster9.csv')

cluster10.markers <- FindMarkers(Diabetes.combined, ident.1 = 10, min.pct = 0.25)
write.csv(cluster10.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster10.csv')

cluster11.markers <- FindMarkers(Diabetes.combined, ident.1 = 11, min.pct = 0.25)
write.csv(cluster11.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster11.csv')

cluster12.markers <- FindMarkers(Diabetes.combined, ident.1 = 12, min.pct = 0.25)
write.csv(cluster12.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster12.csv')

cluster13.markers <- FindMarkers(Diabetes.combined, ident.1 = 13, min.pct = 0.25)
write.csv(cluster13.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster13.csv')

cluster14.markers <- FindMarkers(Diabetes.combined, ident.1 = 14, min.pct = 0.25)
write.csv(cluster14.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster14.csv')

cluster15.markers <- FindMarkers(Diabetes.combined, ident.1 = 15, min.pct = 0.25)
write.csv(cluster15.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster15.csv')

cluster16.markers <- FindMarkers(Diabetes.combined, ident.1 = 16, min.pct = 0.25)
write.csv(cluster16.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster16.csv')

cluster17.markers <- FindMarkers(Diabetes.combined, ident.1 = 17, min.pct = 0.25)
write.csv(cluster17.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster17.csv')

cluster18.markers <- FindMarkers(Diabetes.combined, ident.1 = 18, min.pct = 0.25)
write.csv(cluster18.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster18.csv')

cluster19.markers <- FindMarkers(Diabetes.combined, ident.1 = 19, min.pct = 0.25)
write.csv(cluster19.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/db.markerscombined_cluster19.csv')


diabetescombinedall.markers <- FindAllMarkers(Diabetes.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
diabetescombinedall.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(diabetescombinedall.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/diabetescombinedall.markers.csv')


diabetescombinedall1.markers <- FindAllMarkers(Diabetes.combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
diabetescombinedall1.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(diabetescombinedall1.markers, 'C:/Users/mad1188/Desktop/Diabetes_new/control vs BMP-7/MARKERS ANALYSIS/diabetescombinedall1.markers.csv')

cluster0_conserved_markers <- FindConservedMarkers(Diabetes.combined,
                                                   ident.1 = 0,
                                                   grouping.var = "BMP7",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)

cluster1_conserved_markers <- FindConservedMarkers(Diabetes.combined,
                                                   ident.1 = 1,
                                                   grouping.var = "BMP7",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
Diabetes.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Diabetes.combined, features = top10$gene) + NoLegend()

Diabetes.combined1 <- RenameIdents(Diabetes.combined, `0` = "Acinar", `1` = "A", `2` = "B",
                                `3` = "C", `4` = "D", `5` = "Ductal", `6` = "E", `7` = "F", `8` = "G", `9` = "H",
                                `10` = "I", `11` = "J", `12` = "K", `13` = "L", `14` = "Stellate", `15` = "M",`16` = "N",`17` = "Beta",`18` = "Vascular",`19` = "O")
DimPlot(Diabetes.combined1, label = TRUE, raster = FALSE) 
