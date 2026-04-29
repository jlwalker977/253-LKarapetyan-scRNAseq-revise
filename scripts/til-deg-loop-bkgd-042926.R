cat("starting lib loading")
shhh <- suppressPackageStartupMessages
shhh(library(Matrix))
shhh(library(DESeq2))
shhh(library(celldex))
shhh(library(hdf5r))
shhh(library(dplyr))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(tidyr))
shhh(library(ggplot2))
shhh(library(patchwork))
shhh(library(Seurat))
shhh(library(SingleR))
shhh(library(ggrepel))
shhh(library(qs))
shhh(library(tidytext))
shhh(library(pheatmap))
shhh(library(ComplexHeatmap))
shhh(library(ggpubr))
shhh(library(rstatix))
shhh(library(SingleR))
shhh(library(edgeR))
shhh(library(slingshot))
shhh(library(RColorBrewer))
shhh(library(tradeSeq))
shhh(library(SingleCellExperiment))
shhh(library(tidyverse))
shhh(library(celda))

cat("loading functions")
### SAVE FIGS, VOLCANO FUNCTIONS
save.figs <- function(p1,out.dir,base.fn,width,height){
  fn=paste0(base.fn,'.png')
  ggsave(file.path(out.dir,fn),
         plot = p1,
         width = width,
         height = height,
         units = 'in')
  fn=paste0(base.fn,'.pdf')
  ggsave(file.path(out.dir,fn),
         plot = p1,
         width = width,
         height = height,
         units = 'in')
}

volcano.plot <- function(de,p_thresh,log_thresh,title){
  max_abs_x <- max(abs(de$log2FC))
  suppressWarnings({
    p1 <- ggplot(data=de, aes(x=log2FC, y=-log10(p.adj+1e-300), col=diffexpressed, label=delabel)) +
      geom_point() + 
      theme_minimal() +
      geom_text_repel() +
      scale_color_manual(values=c("blue", "black", "red")) +
      geom_vline(xintercept=c(-log_thresh, log_thresh), col="red") +
      geom_hline(yintercept=-log10(p_thresh), col="red") +
      ggtitle(title) +  xlim(-1.05 * max_abs_x, 1.05 * max_abs_x)
  })
  return(p1)
}

filter_uncharacterized_genes <- function(de) {
  
  remove_patterns <- c(
    "^ENSG", "^ENS",
    "^LOC", "^LINC",
    "-AS[0-9]*$", "-OT[0-9]*", "-DT$", "-IT[0-9]*",
    "^AC[0-9]+\\.[0-9]+", "^AL[0-9]+\\.[0-9]+",
    "^AP[0-9]+\\.[0-9]+", "^AF[0-9]+\\.[0-9]+",
    "^AJ[0-9]+\\.[0-9]+", "^BX[0-9]+\\.[0-9]+",
    "^CR[0-9]+\\.[0-9]+", "^CT[0-9]+\\.[0-9]+",
    "^CU[0-9]+\\.[0-9]+", "^FP[0-9]+\\.[0-9]+",
    "^Z[0-9]+\\.[0-9]+",
    "^[A-Z]{2}[0-9]+\\.[0-9]+",
    "^KIAA[0-9]+",
    "^HIST[0-9]", "^H[1-4][A-Z]", "^H[1-4]-",
    "^MT-", "^MTRNR", "^MTRN",
    "^FAM[0-9]+",
    "^C[0-9]+orf[0-9]+"
  )
  
  combined_pattern <- paste(remove_patterns, collapse = "|")
  
  genes <- if ("Gene" %in% colnames(de)) de$Gene else rownames(de)
  keep  <- !grepl(combined_pattern, genes, perl = TRUE)
  de[keep, ]
}

### BEST HEATMAP FUNCTION
# Annotation color generator
ann.colors <- function(annotation_col) {
  color_map <- list(
    PFS.6month      = c("NoPrg6M" = "coral",     "Prg6M"        = "navy"),
    timepoint       = c("Baseline" = "grey90",    "Week4"        = "hotpink1", "Week12.or.16" = "purple"),
    treatment       = c("combo"    = "magenta",   "nivo"         = "cyan",     "rela"         = "grey60"),
    gen_cell_types2 = c("CD3plus"  = "yellowgreen","Macrophages" = "violet",   "CD4"          = "#FDC086")
  )
  annotation_col <- as.data.frame(annotation_col)
  result <- list()
  for (col in colnames(annotation_col)) {
    if (!col %in% names(color_map)) next
    actual_levels <- unique(as.character(annotation_col[[col]]))
    actual_levels <- actual_levels[!is.na(actual_levels)]
    defined       <- color_map[[col]]
    result[[col]] <- defined[names(defined) %in% actual_levels]
  }
  result
}

# pheatmap wrapper
pheatmap.wrapper <- function(lognorm.dat, annotation_col, ann_colors, title,
                             cluster_rows = TRUE,
                             cluster_cols = FALSE) {
  annotation_col <- as.data.frame(annotation_col)
  ann_colors     <- ann_colors[names(ann_colors) %in% colnames(annotation_col)]
  
  for (col in colnames(annotation_col)) {
    vals <- as.character(annotation_col[[col]])
    lvls <- if (col %in% names(ann_colors)) {
      names(ann_colors[[col]])[names(ann_colors[[col]]) %in% unique(vals)]
    } else {
      sort(unique(vals))
    }
    annotation_col[[col]] <- factor(vals, levels = lvls)
  }
  
  pheatmap::pheatmap(
    lognorm.dat,
    main              = title,
    scale             = "none",
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    cluster_rows      = cluster_rows,
    cluster_cols      = cluster_cols,
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    show_colnames     = FALSE,
    show_rownames     = TRUE
  )
}

# Manual row-scaling
row.scale <- function(mat) {
  t(apply(mat, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }))
}

# Main heatmap function
heatmap.plot <- function(de, sub, title, ds.pt = 100,
                         ndeg_top_bottom = 15,
                         p.thresh        = 0.05,
                         plot.tag        = c('PFS.6month', 'timepoint', 'treatment', 'gen_cell_types2')) {
  
  if (!"avg_log2FC" %in% colnames(de) && "log2FC" %in% colnames(de))
    de$avg_log2FC <- de$log2FC
  if (!"p_val_adj" %in% colnames(de) && "p.adj" %in% colnames(de))
    de$p_val_adj <- de$p.adj
  
  top      <- de %>% filter(p_val_adj < p.thresh, avg_log2FC >  0) %>% top_n(ndeg_top_bottom, wt =  avg_log2FC)
  bottom   <- de %>% filter(p_val_adj < p.thresh, avg_log2FC <  0) %>% top_n(ndeg_top_bottom, wt = -avg_log2FC)
  use.genes <- unique(rbind(top, bottom)$Gene)
  
  if (length(use.genes) <= 1) {
    warning("Fewer than 2 DEGs passed filters for: ", title)
    return(NA)
  }
  
  Idents(sub) <- 'gen_cell_types2'
  suppressWarnings({ sub <- subset(sub, downsample = ds.pt) })
  
  use.assay <- names(which.max(sapply(Assays(sub), function(a) {
    sum(use.genes %in% rownames(sub[[a]]))
  })))
  DefaultAssay(sub) <- use.assay
  
  annotation_col <- as.data.frame(
    sub@meta.data[, intersect(plot.tag, colnames(sub@meta.data)), drop = FALSE]
  )
  annotation_col[] <- lapply(annotation_col, as.character)
  
  available.genes <- intersect(use.genes, rownames(sub[[use.assay]]))
  if (length(available.genes) <= 1) {
    warning("Fewer than 2 available genes for: ", title)
    return(NA)
  }
  
  lognorm.dat <- as.matrix(
    GetAssayData(sub, assay = use.assay, slot = "data")[available.genes,
                                                        rownames(annotation_col),
                                                        drop = FALSE]
  )
  
  nonzero.rows <- apply(lognorm.dat, 1, function(x) any(x != 0))
  if (sum(nonzero.rows) <= 1) {
    warning("Fewer than 2 non-zero genes for: ", title)
    return(NA)
  }
  lognorm.dat <- lognorm.dat[nonzero.rows, , drop = FALSE]
  lognorm.dat <- row.scale(lognorm.dat)
  
  ann_colors <- ann.colors(annotation_col)
  pheatmap.wrapper(lognorm.dat, annotation_col, ann_colors, title)
}
cat("TIL setup for loop")
dat_t = qread("dat-til-redo-scdblfinder-doublets-removed-042426-light.qs")
#TIL claude fix 4
#TIL LOOP SETUP
prj = "HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise"
base = file.path('/ix/rbao/Projects', prj)
data = file.path(base,'results','jwork')
scripts= file.path(base,'scripts')
results = file.path(base,'results','jwork')
cell.types <- unique(dat_t@meta.data$gen_cell_types2)
##TIL LOOP SETUP 2
use.tags <- c( 'PFS.6month','treatment', 'timepoint') # 'treatment', 'timepoint',
test <- 'MAST'
log.name <- "avg_log2FC"
p.name <- 'p_val_adj'
p = 0.001
log_threshes = c(1, 1.5)
save.filtered <- TRUE
width = 8
height = 6
v=11
suppressWarnings({
  out.gene.lists <- file.path(results,sprintf('DEG_contrasts_TIL_wbgf_%s_v%d',test,v),'csv')
  dir.create(out.gene.lists,recursive = TRUE)
  filt.out.gene.lists <- file.path(results,sprintf('DEG_contrasts_TIL_wbgf_%s_v%d',test,v),'filtered_csv')
  dir.create(filt.out.gene.lists,recursive = TRUE)
  out.de.volcano <- file.path(results,sprintf('DEG_contrasts_TIL_wbgf_%s_v%d',test,v),'volcano')
  dir.create(out.de.volcano,recursive=TRUE)
  out.de.heatmap <- file.path(results,sprintf('DEG_contrasts_TIL_wbgf_%s_v%d',test,v),'heatmap')
  dir.create(out.de.heatmap,recursive=TRUE)
})
plot.volc <- TRUE
plot.heat <- TRUE
min.pct <- 0.03
min.cells.feat = 20
tissue.types <- c('TIL')
treatments <- c('rela', 'nivo', 'combo')
timepoints <- c('Baseline','Week4')
PFS.6month <- c('Prg6M','NoPrg6M')

cat("starting TIL deg loop")
#TIL LOOP
system.time({
  for (log_thresh in log_threshes){
    for (time in timepoints){
      for (tag in use.tags){
        #--- Define idents ---
        if (tag == 'PFS.6month'){
          idents <- c('NoPrg6M','Prg6M')
        } else if (tag == 'timepoint'){
          idents <- c('Week4','Baseline')
        } else if (tag == 'treatment'){
          idents <- c('combo','rela','nivo')
        } else {
          idents <- c('yes','no')
        }
        # --- Define pairwise contrasts ---
        if (tag == 'treatment'){
          contrast.pairs <- list(c('combo','rela'), c('combo','nivo'), c('rela','nivo'))
        } else {
          contrast.pairs <- list(idents[1:2])
        }
        # --- Timepoint subsetting: skip redundant iterations when tag == 'timepoint' ---
        if (tag == 'timepoint'){
          if (time != timepoints[1]) next
          time.sub <- dat_t
        } else {
          time.sub <- subset(dat_t, subset = timepoint == time)
        }
        comp.sub <- subset(time.sub, subset = !!sym(tag) %in% idents)
        for (tissue in tissue.types){
          tissue.sub <- subset(comp.sub, subset = sample_type == tissue)
          for (trt in treatments){
            # --- Treatment subsetting: skip redundant iterations when tag == 'treatment' ---
            if (tag == 'treatment'){
              if (trt != treatments[1]) next
              treat.sub <- tissue.sub
            } else {
              treat.sub <- subset(tissue.sub, subset = treatment == trt)
            }
            for (cell in cell.types){
              if (ncol(treat.sub) == 0){
                message('Skipping: no cells for trt=', trt, ' tissue=', tissue, ' tag=', tag)
                next
              }
              if (cell != 'all.cells'){
                sub <- suppressWarnings(subset(treat.sub, subset = gen_cell_types2 == cell))
              } else {
                sub <- treat.sub
              }
              if (ncol(sub) == 0){
                message('Skipping: no cells for cell=', cell)
                next
              }
              DefaultAssay(sub) <- "RNA"
              Idents(sub) <- tag
              # --- Title ---
              if (tag == 'PFS.6month'){
                title = sprintf('NoPrg6m vs Prg6M, %s DEGs', cell)
              } else if (tag == 'timepoint'){
                title = sprintf('Week4 vs Baseline, %s DEGs', cell)
              } else if (tag == 'treatment'){
                title = sprintf('Treatment comparisons, %s DEGs', cell)
              } else {
                title = sprintf('yes vs no, %s DEGs', cell)
              }
              title = paste(tissue, title)
              for (pair in contrast.pairs){
                ident.1 <- pair[1]
                ident.2 <- pair[2]
                # --- Guard: both idents must exist in subset ---
                available.idents <- unique(as.character(sub[[tag]][,1]))
                if (!ident.1 %in% available.idents || !ident.2 %in% available.idents){
                  message('Skipping: missing ident. Need: ', ident.1, ' and ', ident.2,
                          '. Found: ', paste(available.idents, collapse=', '))
                  next
                }
                if (tag == 'treatment'){
                  title = paste(tissue, sprintf('%s vs %s, %s DEGs', ident.1, ident.2, cell))
                }
                # --- For timepoint tag, use 'all' as time label in filename ---
                time.label <- if (tag == 'timepoint') 'alltime' else time
                trt.label  <- if (tag == 'treatment') 'alltrt'  else trt
                contrast.str <- paste0(time.label,'.',trt.label,'.',tissue,'.',
                                       tag,'.',ident.1,'.vs.',ident.2)
                base.fn <- sprintf('%s.%s.DEG.p%1.3f.log.%1.1f.%s.%d.minpct.v%d',
                                   contrast.str, cell, p, log_thresh, test, min.pct*100, v)
                base.fn <- stringr::str_replace_all(base.fn, ' ', '.')
                fn <- paste0(base.fn, '.csv')
                if (file.exists(file.path(out.gene.lists, fn))){
                  print(paste('Loading', fn))
                  de <- read.csv(file.path(out.gene.lists, fn), row.names=1)
                  de <- filter_uncharacterized_genes(de)              # ← ADDED
                } else {
                  print(paste('Generating', fn))
                  markers <- tryCatch({
                    FindMarkers(object = sub,
                                ident.1 = ident.1,
                                ident.2 = ident.2,
                                base = 2,
                                min.pct = min.pct,
                                min.cells.feature = min.cells.feat,
                                min.cells.group = 50,
                                test.use = 'MAST')
                  }, error = function(e){
                    message('FindMarkers failed for ', contrast.str, ': ', e$message)
                    return(NULL)
                  })
                  if (is.null(markers) || nrow(markers) == 0){
                    message('No markers found for ', contrast.str, ', skipping')
                    next
                  }
                  de <- markers
                  de$Gene <- row.names(markers)
                  de$log2FC <- de[[log.name]]
                  max_fc <- 100
                  de <- de %>%
                    mutate(
                      log2FC = case_when(
                        log2FC == Inf    ~ max_fc,
                        log2FC > max_fc  ~ max_fc,
                        log2FC < -max_fc ~ -max_fc,
                        log2FC == -Inf   ~ -max_fc,
                        TRUE             ~ log2FC
                      )
                    ) %>%
                    filter(!is.na(p_val_adj))
                  de$p.adj <- de[[p.name]]
                  de$diffexpressed <- "NO"
                  de$diffexpressed[de$log2FC > log_thresh  & de$p.adj < p] <- "UP"
                  de$diffexpressed[de$log2FC < -log_thresh & de$p.adj < p] <- "DOWN"
                  de$diffexpressed <- as.character(de$diffexpressed)
                  de$delabel <- NA
                  de$delabel[de$diffexpressed != "NO"] <- de$Gene[de$diffexpressed != "NO"]
                  de <- de %>% arrange(-log2FC, p_val_adj)
                  de$contrast <- contrast.str
                  de <- filter_uncharacterized_genes(de)              # ← ADDED
                  write.csv(de, file.path(out.gene.lists, fn))
                }
                if (save.filtered){
                  filt.de <- de %>% filter(diffexpressed %in% c('UP','DOWN'))
                  fn.filt <- paste0('filt.', base.fn, '.csv')
                  write.csv(filt.de, file.path(filt.out.gene.lists, fn.filt))
                }
                if (plot.volc){
                  p1 <- volcano.plot(de, p, log_thresh, title)
                  save.figs(p1, out.de.volcano, base.fn, width=width, height=height)
                }
                if (plot.heat){
                  if (cell == 'all.cells'){ ds.pt <- 10;  ndeg <- 10
                  } else {                  ds.pt <- 200; ndeg <- 15 }
                  p2 <- tryCatch({
                    heatmap.plot(de, sub, title, plot.tag = tag,
                                 ndeg_top_bottom = ndeg, p.thresh = p,
                                 ds.pt = ds.pt)
                  }, error = function(e) {
                    warning("heatmap.plot error for ", base.fn, ": ", e$message)
                    NA
                  })
                  if (is.list(p2)) {
                    save.figs(p2, out.de.heatmap, base.fn, width=10, height=5)
                  } else {
                    print(paste('Too few genes for heatmap:', base.fn))
                  }
                }
              } # end pair loop
            } # end cell loop
          } # end trt loop
        } # end tissue loop
      } # end tag loop
    } # end time loop
  } # end log_thresh loop
}) # end system.time
cat(" done! hooray!")


