# LOAD LIBRARIES ####
# Restart Rstudio or R
# Run the following code once you have Seurat installed

install.packages('clustree')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

suppressWarnings(
  {
    library(ggplot2)
    library(cowplot)
    library(Matrix)
    library(ggridges)
    library(ggrepel)
    library(dplyr)
    library(Seurat)
    library(monocle3)
    library(plotly)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(EnhancedVolcano)
  }
)

# CONFIRM CORRECT INSTALL ####
# Confirm package version of Seurat and Monocle
packageVersion("Seurat")
packageVersion("monocle3")

# Set global environment parameter
options(future.globals.maxSize = 8000 * 1024^2)

# OBJECT SETUP AND NORMALIZATION ####
# STEP 1: Load 10X data ####
ctrl1_1.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST2-CTRL\filtered_feature_bc_matrix)")
ctrl1_2.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST-2-CTRL-1\filtered_feature_bc_matrix)")
ctrl1_3.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST-2-CTRL-2\filtered_feature_bc_matrix)")
bmp71_1.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST2-BMP-7\filtered_feature_bc_matrix)")
bmp71_2.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST-2-BMP-7\filtered_feature_bc_matrix)")
bmp71_3.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST--2-bmp-7\filtered_feature_bc_matrix)")

ctrl2_1.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST-3-CTRL-1\filtered_feature_bc_matrix)")
ctrl2_2.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST-3-CTRL-2\filtered_feature_bc_matrix)")
ctrl2_3.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST-3-CTRL-3\filtered_feature_bc_matrix)")
bmp72_1.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST3-BMP-7-1\filtered_feature_bc_matrix)")
bmp72_2.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST3-BMP-7-2\filtered_feature_bc_matrix)")
bmp72_3.data <- Read10X(data.dir = r"(C:\Users\mqadir\Box\Fahd Shared with JDB\Sequencing Data\cellranger_count\ST3-BMP-7-3\filtered_feature_bc_matrix)")

# STEP 2: Create Seurat objects ####
ctrl1_1 <- CreateSeuratObject(counts = ctrl1_1.data, min.features = 500)
ctrl1_2 <- CreateSeuratObject(counts = ctrl1_2.data, min.features = 500)
ctrl1_3 <- CreateSeuratObject(counts = ctrl1_3.data, min.features = 500)
bmp71_1 <- CreateSeuratObject(counts = bmp71_1.data, min.features = 500)
bmp71_2 <- CreateSeuratObject(counts = bmp71_2.data, min.features = 500)
bmp71_3 <- CreateSeuratObject(counts = bmp71_3.data, min.features = 500)

ctrl2_1 <- CreateSeuratObject(counts = ctrl2_1.data, min.features = 500)
ctrl2_2 <- CreateSeuratObject(counts = ctrl2_2.data, min.features = 500)
ctrl2_3 <- CreateSeuratObject(counts = ctrl2_3.data, min.features = 500)
bmp72_1 <- CreateSeuratObject(counts = bmp72_1.data, min.features = 500)
bmp72_2 <- CreateSeuratObject(counts = bmp72_2.data, min.features = 500)
bmp72_3 <- CreateSeuratObject(counts = bmp72_3.data, min.features = 500)

# Sample specific Metadata addition
ctrl1_1$sample <- "ctrl1_1"
ctrl1_2$sample <- "ctrl1_2"
ctrl1_3$sample <- "ctrl1_3"
bmp71_1$sample <- "bmp71_1"
bmp71_2$sample <- "bmp71_2"
bmp71_3$sample <- "bmp71_3"

ctrl2_1$sample <- "ctrl2_1"
ctrl2_2$sample <- "ctrl2_2"
ctrl2_3$sample <- "ctrl2_3"
bmp72_1$sample <- "bmp72_1"
bmp72_2$sample <- "bmp72_2"
bmp72_3$sample <- "bmp72_3"

# Sex specific Metadata addition
ctrl1_1$sample <- "NA"
ctrl1_2$sample <- "NA"
ctrl1_3$sample <- "NA"
bmp71_1$sample <- "NA"
bmp71_2$sample <- "NA"
bmp71_3$sample <- "NA"

ctrl2_1$sample <- "NA"
ctrl2_2$sample <- "NA"
ctrl2_3$sample <- "NA"
bmp72_1$sample <- "NA"
bmp72_2$sample <- "NA"
bmp72_3$sample <- "NA"

# Treatment specific Metadata addition
ctrl1_1$treatment <- "ctrl_1"
ctrl1_2$treatment <- "ctrl_1"
ctrl1_3$treatment <- "ctrl_1"
bmp71_1$treatment <- "bmp7_1"
bmp71_2$treatment <- "bmp7_1"
bmp71_3$treatment <- "bmp7_1"

ctrl2_1$treatment <- "ctrl_2"
ctrl2_2$treatment <- "ctrl_2"
ctrl2_3$treatment <- "ctrl_2"
bmp72_1$treatment <- "bmp7_2"
bmp72_2$treatment <- "bmp7_2"
bmp72_3$treatment <- "bmp7_2"

# STEP 3: Thresholding ####
# The operator can add columns to object metadata. This is a great place to stash QC stats
ctrl1_1[["percent.mt"]] <- PercentageFeatureSet(object = ctrl1_1, pattern = "^MT-")
ctrl1_2[["percent.mt"]] <- PercentageFeatureSet(object = ctrl1_2, pattern = "^MT-")
ctrl1_3[["percent.mt"]] <- PercentageFeatureSet(object = ctrl1_3, pattern = "^MT-")
bmp71_1[["percent.mt"]] <- PercentageFeatureSet(object = bmp71_1, pattern = "^MT-")
bmp71_2[["percent.mt"]] <- PercentageFeatureSet(object = bmp71_2, pattern = "^MT-")
bmp71_3[["percent.mt"]] <- PercentageFeatureSet(object = bmp71_3, pattern = "^MT-")

ctrl2_1[["percent.mt"]] <- PercentageFeatureSet(object = ctrl2_1, pattern = "^MT-")
ctrl2_2[["percent.mt"]] <- PercentageFeatureSet(object = ctrl2_2, pattern = "^MT-")
ctrl2_3[["percent.mt"]] <- PercentageFeatureSet(object = ctrl2_3, pattern = "^MT-")
bmp72_1[["percent.mt"]] <- PercentageFeatureSet(object = bmp72_1, pattern = "^MT-")
bmp72_2[["percent.mt"]] <- PercentageFeatureSet(object = bmp72_2, pattern = "^MT-")
bmp72_3[["percent.mt"]] <- PercentageFeatureSet(object = bmp72_3, pattern = "^MT-")

# QC information before thresholding
summary(head(ctrl1_1@meta.data))
summary(head(ctrl1_2@meta.data))
summary(head(ctrl1_3@meta.data))
summary(head(bmp71_1@meta.data))
summary(head(bmp71_2@meta.data))
summary(head(bmp71_3@meta.data))

summary(head(ctrl2_1@meta.data))
summary(head(ctrl2_2@meta.data))
summary(head(ctrl2_3@meta.data))
summary(head(bmp72_1@meta.data))
summary(head(bmp72_2@meta.data))
summary(head(bmp72_3@meta.data))

# Visualize QC metrics as a violin plot
VlnPlot(object = ctrl1_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl1_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl1_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp71_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp71_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp71_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(object = ctrl2_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl2_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl2_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp72_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp72_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp72_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# RNA based cell thresholding
ctrl1_1 <- subset(x = ctrl1_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
ctrl1_2 <- subset(x = ctrl1_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
ctrl1_3 <- subset(x = ctrl1_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
bmp71_1 <- subset(x = bmp71_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
bmp71_2 <- subset(x = bmp71_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
bmp71_3 <- subset(x = bmp71_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

ctrl2_1 <- subset(x = ctrl2_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
ctrl2_2 <- subset(x = ctrl2_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
ctrl2_3 <- subset(x = ctrl2_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
bmp72_1 <- subset(x = bmp72_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
bmp72_2 <- subset(x = bmp72_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)
bmp72_3 <- subset(x = bmp72_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# QC information after thresholding
summary(head(ctrl1_1@meta.data))
summary(head(ctrl1_2@meta.data))
summary(head(ctrl1_3@meta.data))
summary(head(bmp71_1@meta.data))
summary(head(bmp71_2@meta.data))
summary(head(bmp71_3@meta.data))

summary(head(ctrl2_1@meta.data))
summary(head(ctrl2_2@meta.data))
summary(head(ctrl2_3@meta.data))
summary(head(bmp72_1@meta.data))
summary(head(bmp72_2@meta.data))
summary(head(bmp72_3@meta.data))

# Visualize QC metrics post thresholding as a violin plot
VlnPlot(object = ctrl1_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl1_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl1_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp71_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp71_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp71_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

VlnPlot(object = ctrl2_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl2_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = ctrl2_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp72_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp72_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(object = bmp72_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#VlnPlot(object = pancreas.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, cols = c("red", "blue"))
# Step 4: Add cell IDs ####
# Add cell IDs
ctrl1_1 <- RenameCells(ctrl1_1, add.cell.id = "ctrl1_1")
ctrl1_2 <- RenameCells(ctrl1_2, add.cell.id = "ctrl1_2")
ctrl1_3 <- RenameCells(ctrl1_3, add.cell.id = "ctrl1_3")
bmp71_1 <- RenameCells(bmp71_1, add.cell.id = "bmp71_1")
bmp71_2 <- RenameCells(bmp71_2, add.cell.id = "bmp71_2")
bmp71_3 <- RenameCells(bmp71_3, add.cell.id = "bmp71_3")

ctrl2_1 <- RenameCells(ctrl2_1, add.cell.id = "ctrl2_1")
ctrl2_2 <- RenameCells(ctrl2_2, add.cell.id = "ctrl2_2")
ctrl2_3 <- RenameCells(ctrl2_3, add.cell.id = "ctrl2_3")
bmp72_1 <- RenameCells(bmp72_1, add.cell.id = "bmp72_1")
bmp72_2 <- RenameCells(bmp72_2, add.cell.id = "bmp72_2")
bmp72_3 <- RenameCells(bmp72_3, add.cell.id = "bmp72_3")

# Step 5: Merge Datasets
# Based on comment to Issue #4753 https://github.com/satijalab/seurat/issues/4753
# We use RPCA to yield conserved mapping and set Tx as control reference samples
# Merge bmp7 datasets
bmp7.list <- list("ctrl1_1" = ctrl1_1, "ctrl1_2" = ctrl1_2, "ctrl1_3" = ctrl1_3,
                  "bmp71_1" = bmp71_1, "bmp71_2" = bmp71_2, "bmp71_3" = bmp71_3,
                  "ctrl2_1" = ctrl2_1, "ctrl2_2" = ctrl2_2, "ctrl2_3" = ctrl2_3,
                  "bmp72_1" = bmp72_1, "bmp72_2" = bmp72_2, "bmp72_3" = bmp72_3)

# Step 6: Data normalization
#Normalise data
bmp7.list <- lapply(X = bmp7.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

# Step 7: Feature selection
# Select features for downstream integration
bmp7.features <- SelectIntegrationFeatures(object.list = bmp7.list)
bmp7.list <- lapply(X = bmp7.list, FUN = function(x) {
  x <- ScaleData(x, features = bmp7.features, verbose = FALSE)
  x <- RunPCA(x, features = bmp7.features, verbose = FALSE)
})

# Step 8: Anchor identification and data integration
# Identify anchors and integrate dataset
bmp7.list[c(1, 2, 3)] # check that you are correctly picking up control datasets for refrence integration
bmp7.anchors <- FindIntegrationAnchors(object.list = bmp7.list, reference = c(1,2,3), # Takes 18min 13sec to run when using cca as reduction
                                           reduction = "rpca", dims = 1:30, verbose = TRUE)
bmp7.integrated <- IntegrateData(anchorset = bmp7.anchors, dims = 1:30, verbose = TRUE)

# Step 9: Linear dimensionality assessment
# Look at your default assay
DefaultAssay(object = bmp7.integrated)

# Change default assay to integrated, to view dimensionality
DefaultAssay(object = bmp7.integrated) <- "integrated"

# Scaling this is weird, but as done in https://satijalab.org/seurat/articles/integration_large_datasets.html
bmp7.integrated <- ScaleData(bmp7.integrated, verbose = FALSE)

# Dimensionality assessment using PCA analysis
bmp7.integrated <- RunPCA(bmp7.integrated, features = VariableFeatures(object = bmp7.integrated))

# Examine data dimensionality
ElbowPlot(bmp7.integrated)

# Step 9a: CLUSTER ANALYSIS ####
# Clustering, we run multiple permutations to allow clustree to analyze optimal clustering resolution.
bmp7.integrated <- FindNeighbors(object = bmp7.integrated, dims = 1:30)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0.1)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0.2)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0.3)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0.4)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0.5)
bmp7.integrated <- FindClusters(object = bmp7.integrated, resolution = 0.6)

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(bmp7.integrated, prefix = "integrated_snn_res.")

# Based of clustree assessment choose res = 0.5
bmp7.integrated <- FindClusters(bmp7.integrated, resolution = 0.5)

# Alternatively build a cluster tree
DefaultAssay(object = bmp7.integrated) <- "integrated"
bmp7.integrated=BuildClusterTree(bmp7.integrated, slot = "scale.data")
PlotClusterTree(bmp7.integrated)

# Step 10: non-linear dimensionality assessment ####
# Run PCA and UMAP calculations
bmp7.integrated <- RunUMAP(bmp7.integrated, dims = 1:30)

# Change default assay to integrated, to view dimensionality
Idents(bmp7.integrated) <- "treatment"
Idents(bmp7.integrated) <- "sex"
Idents(bmp7.integrated) <- "sample"
Idents(bmp7.integrated) <- "seurat_clusters"
#Idents(pancreas.integrated) <- "integrated_snn_res.0.3"
DimPlot(bmp7.integrated, reduction = "umap", label = FALSE)

#Visualize gene expression
DefaultAssay(object = bmp7.integrated) <- "RNA"
DefaultAssay(object = bmp7.integrated)
FeaturePlot(object = bmp7.integrated,
            features = c("INS", "GCG", "SST", "PPY", "GHRL",
                         "KRT19", "CPA1",
                         "COL1A1", "VWF", "SOX10",
                         "TPSAB1", "SDS", "TRAC"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            #max.cutoff = 1,
            slot = 'counts',
            order = TRUE)

FeaturePlot(object = bmp7.integrated,
            features = c("SST"),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            #min.cutoff = 0,
            max.cutoff = 100,
            slot = 'counts',
            order = TRUE)
