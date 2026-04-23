library(celda)
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("decontX")
library(celda)
sce_t = as.SingleCellExperiment(dat_tmde1)
sce_p = as.SingleCellExperiment(dat_pmde1)
qsave(sce_t, file="sce-til-dat-t-obj-for-decont-doubletfinder.qs")
qsave(sce_p, file="sce-pbmc-dat-p-obj-for-decont-doubletfinder.qs")
Idents(sce_t) = "gen_cell_types2"
sce_t = decontX(sce_t, z=sce_t$gen_cell_types2)
sce_t
DimPlot(sce_t, reduction = "decontX_clusters")
dat_tmde1$decontX_contamination <- sce_t$decontX_contamination
dat_tmde1[["decontX"]] <- CreateAssayObject(counts = round(decontXcounts(sce_t)))
cat("decontX complete\n")
cat("Median contamination score:", median(dat_tmde1$decontX_contamination), "\n")
p_cont <- FeaturePlot(dat_tmde1,
features = "decontX_contamination",
cols     = c("lightgrey", "coral")) +
ggtitle("DecontX Contamination Score")
p_cont
p_cont <- FeaturePlot(dat_tmde1,
features = "decontX_contamination",
cols     = c("grey90", "blueviolet")) +
ggtitle("DecontX Contamination Score")
p_cont
Idents(dat_tmde1) = "gen_cell_types2"
DimPlot(dat_tmde1)
til_umap = DimPlot(dat_tmde1)
til_umap + p_cont
head(dat_tmde1@meta.data)
summary(dat_tmde1@meta.data$decontX_contamination)
DefaultAssay(dat_tmde1) <- "decontX"
p_after <- FeaturePlot(dat_tmde1,
features = c("MS4A1", "CD79A"),
ncol     = 2) &
ggtitle("After decontX")
print(p_cont)
print(p_before / p_after)
print(p_cont / p_after)
print(p_cont / p_after, ncol=3)
qsave(dat_tmde1, file="dat-tmde1-til-post-decontx-pre-ms4a1-removal.qs")
DefaultAssay(dat_tmde1) <- "decontX"
ms4a1_counts <- GetAssayData(dat_tmde1,
assay = "decontX",
slot  = "counts")["MS4A1", ]
cd79a_counts <- GetAssayData(dat_tmde1,
assay = "decontX",
slot  = "counts")["CD79A", ]
# Flag true B cells: co-expressing both markers
is_bcell <- ms4a1_counts > 0 & cd79a_counts > 0
cat("\nCells flagged as B cells (MS4A1+ CD79A+):", sum(is_bcell), "\n")
cat("Percentage of total: ", round(mean(is_bcell) * 100, 2), "%\n")
# Keep non-B cells
cells_keep <- names(is_bcell)[!is_bcell]
dat_tmde1_clean <- subset(dat_tmde1, cells = cells_keep)
cat("\nCells before removal:", ncol(dat_tmde1),       "\n")
cat("Cells after removal: ", ncol(dat_tmde1_clean), "\n")
cat("Total removed:       ", ncol(dat_tmde1) - ncol(dat_tmde1_clean), "\n")
qsave(dat_tmde1_clean, file="dat-tmde1-clean-til-post-ms4a1-cd79a-removal.qs")
DefaultAssay(dat_tmde1_clean) <- "RNA"
dat_tmde1_clean <- NormalizeData(dat_tmde1_clean)
dat_tmde1_clean <- FindVariableFeatures(dat_tmde1_clean)
dat_tmde1_clean <- ScaleData(dat_tmde1_clean)
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean,
group.by = "gen_cell_types2",
label    = TRUE) +
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
p_final
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean,
group.by = "gen_cell_types2",
label    = TRUE) +
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
p_final + p_cont + p_umap
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean,
group.by = "gen_cell_types2",
label    = TRUE) +
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
p_final + p_cont + p_umap
p_final + p_cont +p_umap +til_umap
til_umap + p_cont + p_umap + p_final
qsave(dat_tmde1_clean, file="dat-tmde1-clean-post-conserv-removal.qs")
is_bcell <- ms4a1_counts > 0
cat("\nCells flagged as B cells (MS4A1+ CD79A+):", sum(is_bcell), "\n")
cat("Percentage of total: ", round(mean(is_bcell) * 100, 2), "%\n")
# Keep non-B cells
cells_keep <- names(is_bcell)[!is_bcell]
dat_tmde1_clean2 <- subset(dat_tmde1, cells = cells_keep)
cat("\nCells before removal:", ncol(dat_tmde1_clean),       "\n")
cat("Cells after removal: ", ncol(dat_tmde1_clean2), "\n")
cat("Total removed:       ", ncol(dat_tmde1_clean) - ncol(dat_tmde1_clean2), "\n")
qsave(dat_tmde1_clean2, file="dat-tmde1-clean2-til-post-agress-ms4a1-removal\.qs")
# Keep non-B cells
cells_keep <- names(is_bcell)[!is_bcell]
dat_tmde1_clean2 <- subset(dat_tmde1_clean, cells = cells_keep)
cat("\nCells before removal:", ncol(dat_tmde1_clean),       "\n")
cat("Cells after removal: ", ncol(dat_tmde1_clean2), "\n")
cat("Total removed:       ", ncol(dat_tmde1_clean) - ncol(dat_tmde1_clean2), "\n")
qsave(dat_tmde1_clean2, file="dat-tmde1-clean2-til-post-agress-ms4a1-removal.qs")
DefaultAssay(dat_tmde1_clean2) <- "RNA"
dat_tmde1_clean2 <- NormalizeData(dat_tmde1_clean2)
dat_tmde1_clean2 <- FindVariableFeatures(dat_tmde1_clean2)
dat_tmde1_clean2 <- ScaleData(dat_tmde1_clean2)
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean2,
group.by = "gen_cell_types2",
label    = TRUE,     legend=FALSE) +
ggtitle("UMAP after B cell removal")
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean2,
group.by = "gen_cell_types2",
label    = TRUE) +NoLegend
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean2,
group.by = "gen_cell_types2",
label    = TRUE) +NoLegend()
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
til_umap + p_cont + p_umap + p_final
til_umap = DimPlot(dat_tmde1, label = TRUE)+NoLegend()
til_umap + p_cont + p_umap + p_final
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle(features)
features
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle('features')
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("UMAP after B cell removal")
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2)
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4)
p_final
p_prev <- FeaturePlot(dat_tmde1,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4)
p_final+ p_prev
p_prev
p_final
p_prev + p_final
p_prev + p_final
p_final
p_mid = FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4)
p_mid
p_prev
p_mid
p_final

library(celda)
if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("decontX")
library(celda)
sce_t = as.SingleCellExperiment(dat_tmde1)
sce_p = as.SingleCellExperiment(dat_pmde1)
qsave(sce_t, file="sce-til-dat-t-obj-for-decont-doubletfinder.qs")
qsave(sce_p, file="sce-pbmc-dat-p-obj-for-decont-doubletfinder.qs")
Idents(sce_t) = "gen_cell_types2"
sce_t = decontX(sce_t, z=sce_t$gen_cell_types2)
sce_t
DimPlot(sce_t, reduction = "decontX_clusters")
dat_tmde1$decontX_contamination <- sce_t$decontX_contamination
dat_tmde1[["decontX"]] <- CreateAssayObject(counts = round(decontXcounts(sce_t)))
cat("decontX complete\n")
cat("Median contamination score:", median(dat_tmde1$decontX_contamination), "\n")
p_cont <- FeaturePlot(dat_tmde1,
features = "decontX_contamination",
cols     = c("lightgrey", "coral")) +
ggtitle("DecontX Contamination Score")
p_cont
p_cont <- FeaturePlot(dat_tmde1,
features = "decontX_contamination",
cols     = c("grey90", "blueviolet")) +
ggtitle("DecontX Contamination Score")
p_cont
Idents(dat_tmde1) = "gen_cell_types2"
DimPlot(dat_tmde1)
til_umap = DimPlot(dat_tmde1)
til_umap + p_cont
head(dat_tmde1@meta.data)
summary(dat_tmde1@meta.data$decontX_contamination)
DefaultAssay(dat_tmde1) <- "decontX"
p_after <- FeaturePlot(dat_tmde1,
features = c("MS4A1", "CD79A"),
ncol     = 2) &
ggtitle("After decontX")
print(p_cont)
print(p_before / p_after)
print(p_cont / p_after)
print(p_cont / p_after, ncol=3)
qsave(dat_tmde1, file="dat-tmde1-til-post-decontx-pre-ms4a1-removal.qs")
DefaultAssay(dat_tmde1) <- "decontX"
ms4a1_counts <- GetAssayData(dat_tmde1,
assay = "decontX",
slot  = "counts")["MS4A1", ]
cd79a_counts <- GetAssayData(dat_tmde1,
assay = "decontX",
slot  = "counts")["CD79A", ]
# Flag true B cells: co-expressing both markers
is_bcell <- ms4a1_counts > 0 & cd79a_counts > 0
cat("\nCells flagged as B cells (MS4A1+ CD79A+):", sum(is_bcell), "\n")
cat("Percentage of total: ", round(mean(is_bcell) * 100, 2), "%\n")
# Keep non-B cells
cells_keep <- names(is_bcell)[!is_bcell]
dat_tmde1_clean <- subset(dat_tmde1, cells = cells_keep)
cat("\nCells before removal:", ncol(dat_tmde1),       "\n")
cat("Cells after removal: ", ncol(dat_tmde1_clean), "\n")
cat("Total removed:       ", ncol(dat_tmde1) - ncol(dat_tmde1_clean), "\n")
qsave(dat_tmde1_clean, file="dat-tmde1-clean-til-post-ms4a1-cd79a-removal.qs")
DefaultAssay(dat_tmde1_clean) <- "RNA"
dat_tmde1_clean <- NormalizeData(dat_tmde1_clean)
dat_tmde1_clean <- FindVariableFeatures(dat_tmde1_clean)
dat_tmde1_clean <- ScaleData(dat_tmde1_clean)
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean,
group.by = "gen_cell_types2",
label    = TRUE) +
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
p_final
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean,
group.by = "gen_cell_types2",
label    = TRUE) +
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
p_final + p_cont + p_umap
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean,
group.by = "gen_cell_types2",
label    = TRUE) +
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
p_final + p_cont + p_umap
p_final + p_cont +p_umap +til_umap
til_umap + p_cont + p_umap + p_final
qsave(dat_tmde1_clean, file="dat-tmde1-clean-post-conserv-removal.qs")
is_bcell <- ms4a1_counts > 0
cat("\nCells flagged as B cells (MS4A1+ CD79A+):", sum(is_bcell), "\n")
cat("Percentage of total: ", round(mean(is_bcell) * 100, 2), "%\n")
# Keep non-B cells
cells_keep <- names(is_bcell)[!is_bcell]
dat_tmde1_clean2 <- subset(dat_tmde1, cells = cells_keep)
cat("\nCells before removal:", ncol(dat_tmde1_clean),       "\n")
cat("Cells after removal: ", ncol(dat_tmde1_clean2), "\n")
cat("Total removed:       ", ncol(dat_tmde1_clean) - ncol(dat_tmde1_clean2), "\n")
qsave(dat_tmde1_clean2, file="dat-tmde1-clean2-til-post-agress-ms4a1-removal\.qs")
# Keep non-B cells
cells_keep <- names(is_bcell)[!is_bcell]
dat_tmde1_clean2 <- subset(dat_tmde1_clean, cells = cells_keep)
cat("\nCells before removal:", ncol(dat_tmde1_clean),       "\n")
cat("Cells after removal: ", ncol(dat_tmde1_clean2), "\n")
cat("Total removed:       ", ncol(dat_tmde1_clean) - ncol(dat_tmde1_clean2), "\n")
qsave(dat_tmde1_clean2, file="dat-tmde1-clean2-til-post-agress-ms4a1-removal.qs")
DefaultAssay(dat_tmde1_clean2) <- "RNA"
dat_tmde1_clean2 <- NormalizeData(dat_tmde1_clean2)
dat_tmde1_clean2 <- FindVariableFeatures(dat_tmde1_clean2)
dat_tmde1_clean2 <- ScaleData(dat_tmde1_clean2)
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean2,
group.by = "gen_cell_types2",
label    = TRUE,     legend=FALSE) +
ggtitle("UMAP after B cell removal")
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean2,
group.by = "gen_cell_types2",
label    = TRUE) +NoLegend
DefaultAssay(dat_tmde1_clean) <- "RNA"
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("After B cell removal")
p_umap <- DimPlot(dat_tmde1_clean2,
group.by = "gen_cell_types2",
label    = TRUE) +NoLegend()
ggtitle("UMAP after B cell removal")
print(p_final)
print(p_umap)
til_umap + p_cont + p_umap + p_final
til_umap = DimPlot(dat_tmde1, label = TRUE)+NoLegend()
til_umap + p_cont + p_umap + p_final
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle(features)
features
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle('features')
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2) &
ggtitle("UMAP after B cell removal")
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 2)
p_final
p_final <- FeaturePlot(dat_tmde1_clean2,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4)
p_final
p_prev <- FeaturePlot(dat_tmde1,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4)
p_final+ p_prev
p_prev
p_final
p_prev + p_final
p_prev + p_final
p_final
p_mid = FeaturePlot(dat_tmde1_clean,                                               features = c("MS4A1", "CD79A", "IGHG1", "CD27"), ncol=4)
p_mid = FeaturePlot(dat_tmde1_clean,
features = c("MS4A1", "CD79A", "IGHG1", "CD27"),
ncol     = 4)
p_mid
p_prev
p_mid
p_final