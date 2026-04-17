```{r}
#TIL claude fix 4
#TIL LOOP SETUP
prj = "HCC-CBS-253-Hillman-LKarapetyan-scRNAseq-revise"
base = file.path('/ix/rbao/Projects', prj)
data = file.path(base,'results','jwork')
scripts= file.path(base,'scripts')
results = file.path(base,'results','jwork')
cell.types <- unique(dat_tmde1@meta.data$gen_cell_types2)
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
v=7
suppressWarnings({
    out.gene.lists <- file.path(results,sprintf('DEG_contrasts_TIL_%s_v%d',test,v),'csv')
    dir.create(out.gene.lists,recursive = TRUE)
    filt.out.gene.lists <- file.path(results,sprintf('DEG_contrasts_TIL_%s_v%d',test,v),'filtered_csv')
    dir.create(filt.out.gene.lists,recursive = TRUE)
    out.de.volcano <- file.path(results,sprintf('DEG_contrasts_TIL_%s_v%d',test,v),'volcano')
    dir.create(out.de.volcano,recursive=TRUE)
    out.de.heatmap <- file.path(results,sprintf('DEG_contrasts_TIL_%s_v%d',test,v),'heatmap')
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
```
```{r}
#TIL LOOP
system.time({
for (log_thresh in log_threshes){
for (time in timepoints){
for (tag in use.tags){

    # --- Define idents ---
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
      if (time != timepoints[1]) next  # run only once across time loop
      time.sub <- dat_tmde1                # keep all timepoints for comparison
    } else {
      time.sub <- subset(dat_tmde1, subset = timepoint == time)
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

  # FIX: check class — pheatmap returns a list, NA returns logical
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
} # end tag loop  ← moved inside time loop so tag is defined before time.sub logic
} # end time loop
} # end log_thresh loop
}) # end system.time
```

```{r}
# Example for Seurat v4/v5, excluding pseudo, histones, mito, LOC, and ENSEMBL
genes_to_exclude <- c("^MT-", "^mt-", "pseudo", "\\w+P\\d+$", "^LOC", "^HIST", "^ENS\\w*G\\d{11}")
pattern <- paste(genes_to_exclude, collapse = "|")
filtered_genes <- grep(pattern, all_genes, value = TRUE, invert = TRUE)

(^(MT-|mt-)|^LOC|^FAM|^HIST|^pFAM|^ENS\w*G\d{11}|\w+P\d+$|pseudo)

```
