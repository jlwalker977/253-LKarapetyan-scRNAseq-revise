That's a clean QC validation question. Since you have multiple versions of the object at different cleaning stages, you just want to compare `%` cells expressing IG genes across those snapshots:

```r
check_ig_prevalence <- function(seurat_obj, label,
                                 group.by = 'gen_cell_types2') {
  
  ig.patterns <- c("^IGLV", "^IGHV", "^IGKV")
  ig.genes <- grep(paste(ig.patterns, collapse = "|"),
                   rownames(seurat_obj), value = TRUE)
  
  if (length(ig.genes) == 0) {
    message("No IG  genes found in object: ", label)
    return(NULL)
  }
  
  counts <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
  
  # Per-cell: is ANY IG gene detected?
  ig.detected <- colSums(counts[ig.genes, ] > 0) > 0
  
  meta <- seurat_obj@meta.data
  meta$ig.detected <- ig.detected
  meta$object      <- label
  
  # Summary by cell type
  summary <- meta %>%
    group_by(object, .data[[group.by]]) %>%
    summarise(
      n.cells       = n(),
      n.ig.positive = sum(ig.detected),
      pct.ig        = round(100 * mean(ig.detected), 2),
      .groups = 'drop'
    )
  
  return(summary)
}

# --- Run on each cleaning stage ---
# Swap in whatever your intermediate objects are named
results <- bind_rows(
  check_ig_prevalence(dat_raw,    label = "1_raw"),
  check_ig_prevalence(dat_decontx, label = "2_decontX"),
  check_ig_prevalence(dat_tc4c,   label = "3_decontX_dblFinder")
)

# --- Table view ---
results %>% arrange(gen_cell_types2, object)

# --- Plot ---
ggplot(results, aes(x = object, y = pct.ig, 
                    color = gen_cell_types2, 
                    group = gen_cell_types2)) +
  geom_line() +
  geom_point(aes(size = n.cells)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "IG gene prevalence across cleaning steps",
       x     = "Cleaning stage",
       y     = "% cells with any IGLV/IGHV/IGKV detected",
       size  = "N cells",
       color = "Cell type")
```

A few notes on interpreting the output:

- You'd expect `pct.ig` to drop most sharply after **decontX** since ambient RNA from plasma/B cells is the main source of IG contamination in T cell clusters
- **scDblFinder** would help if any T cell/B cell doublets slipped through, but the effect should be smaller
- If `pct.ig` stays flat across stages in a given cell type, that's worth investigating — it could mean those cells are genuine IG expressors or that decontX parameters need tuning
- The `size = n.cells` aesthetic helps flag cases where a percentage drop is just due to losing cells rather than cleaning