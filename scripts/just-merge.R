
library(qs)
library(SingleCellExperiment)
library(Seurat)
prj = "HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise"
base = file.path('/ix/rbao/Projects', prj)
data = file.path(base,'results','jwork')
scripts= file.path(base,'scripts')
results = file.path(base,'results','jwork')
out.dir <- file.path(results, "scDblFinder_tmp1")

# Check all files are there first
saved.files <- list.files(out.dir, pattern = "_singlets.qs", full.names = TRUE)
cat("Found", length(saved.files), "sample files:\n")
print(basename(saved.files))

# Load and merge
cat("Loading samples from disk...\n")
sample_list_clean <- lapply(saved.files, function(f) {
  cat("Loading:", basename(f), "\n")
  qread(f)
})

cat("Merging...\n")
dat_tc4_clean <- merge(sample_list_clean[[1]], y = sample_list_clean[-1])
cat("Final merged cell count:", ncol(dat_tc4_clean), "\n")

dat_tc4_clean <- JoinLayers(dat_tc4_clean)
qsave(dat_tc4_clean,
      file = "dat-til-redo-scdblfinder-doublets-removed-pre-normd-042426-1.qs")