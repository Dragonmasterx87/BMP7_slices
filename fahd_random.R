install.packages("corrplot")
library(corrplot)

M = cor(mtcars)
mtcars


pancreas.integrated <- readRDS(r"(C:\Users\mqadir\Box\Fahd shared with investigators scRNAseq data\Seurat objects\alk3n3.integrated.nomes.rds)")
alk3n3.integrated.nomes <- UpdateSeuratObject(alk3n3.integrated.nomes)

UMAPPlot(alk3n3.integrated.nomes, reduction = "umap")

pancreas.integrated <- UpdateSeuratObject(pancreas.integrated)
UMAPPlot(pancreas.integrated, reduction = "umap")
