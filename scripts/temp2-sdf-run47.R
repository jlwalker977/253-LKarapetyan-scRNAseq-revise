library(qs)
library(Seurat)
library(scDblFinder)
library(SingleCellExperiment)
dat_tc4 = qread("dat-til-redo_decontX-042426.qs")
sample_list <- SplitObject(dat_tc4, split.by = "Patient.Number")
sample_names <- names(sample_list)

# ── Process and save each sample to disk individually ────────────────────
prj = "HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise"
base = file.path('/ix/rbao/Projects', prj)
data = file.path(base,'results','jwork')
scripts= file.path(base,'scripts')
results = file.path(base,'results','jwork')
out.dir <- file.path(results, "scDblFinder_tmp2")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

for (nm in sample_names) {
  out.file <- file.path(out.dir, paste0(nm, "_singlets.qs"))
  
  # Skip if already done — safe to re-run after interruption
  if (file.exists(out.file)) {
    cat("Already processed, skipping:", nm, "\n\n")
    next
  }
  
  s <- sample_list[[nm]]
  cat("Processing sample:", nm, "| cells:", ncol(s), "\n")
  
  DefaultAssay(s) <- "RNA"
  s <- JoinLayers(s)
  
  sce <- as.SingleCellExperiment(s)
  
  sce <- tryCatch({
    scDblFinder(sce)
  }, error = function(e) {
    message("scDblFinder failed for ", nm, ": ", e$message,
            " — skipping doublet removal.")
    return(sce)
  })
  
  s$scDblFinder.class <- if ("scDblFinder.class" %in% colnames(colData(sce))) {
    sce$scDblFinder.class
  } else {
    "singlet"
  }
  s$scDblFinder.score <- if ("scDblFinder.score" %in% colnames(colData(sce))) {
    sce$scDblFinder.score
  } else {
    NA_real_
  }
  
  cat("Doublets found:", sum(s$scDblFinder.class == "doublet"), "\n")
  
  s_clean <- subset(s, subset = scDblFinder.class == "singlet")
  cat("Cells before:", ncol(s), "| after:", ncol(s_clean),
      "| removed:", ncol(s) - ncol(s_clean), "\n\n")
  
  qsave(s_clean, file = out.file)
  
  # Free memory before next sample
  rm(s, s_clean, sce)
  gc()
}

# ── Merge from disk rather than from RAM ─────────────────────────────────
saved.files <- file.path(out.dir, paste0(sample_names, "_singlets.qs"))

cat("Loading and merging samples...\n")
sample_list_clean <- lapply(saved.files, qread)

dat_tc4_clean <- merge(sample_list_clean[[1]], y = sample_list_clean[-1])
cat("Final merged cell count:", ncol(dat_tc4_clean), "\n")

dat_tc4_clean <- JoinLayers(dat_tc4_clean)
qsave(dat_tc4_clean,
      file = "dat-til-redo-scdblfinder-doublets-removed-pre-normd-042426-2.qs")