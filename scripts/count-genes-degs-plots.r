library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# =============================================================================
# PART 1: Process raw combined CSVs → annotated DEG table + gene counts summary
# (equivalent to the Python processing script)
# =============================================================================

process_deg_csv <- function(input_csv, tissue) {
  
  df <- read.csv(input_csv, stringsAsFactors = FALSE)
  
  # Extract cell type from source_file path
  df$cell.type <- str_match(df$source_file, 
                             "NoPrg6M\\.vs\\.Prg6M\\.(.+?)\\.DEG")[, 2]
  df$cell.type <- str_replace_all(df$cell.type, "\\.", " ")
  
  # Build cell.contrasts column
  df$cell.contrasts <- paste0(df$contrast, ".", 
                               str_replace_all(df$cell.type, " ", "."))
  
  # Summary column: direction of regulation
  # Contrast is NoPrg6M vs Prg6M → positive FC = up in non-progressor
  df$summary <- ifelse(df$avg_log2FC > 0,
                       "PFS.6month.up.in.non.progressor",
                       "PFS.6month.up.in.progressor")
  
  # Select final columns
  all_de <- df %>%
    select(Gene, p_val_adj, avg_log2FC, summary,
           pct.1, pct.2, contrast, cell.type, cell.contrasts, diffexpressed)
  
  # Gene count table: how many times each gene appears
  gene_counts <- all_de %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts)
  
  # Directionality breakdown per gene
  dir_counts <- all_de %>%
    group_by(Gene, summary) %>%
    reframe(n = n()) %>%
    pivot_wider(names_from = summary, values_from = n, values_fill = 0)
  
  # Merge counts + directionality
  gene_summary <- gene_counts %>%
    left_join(dir_counts, by = "Gene")
  
  # Save outputs
  write.csv(all_de,      sprintf("all_pfs_%s_DEG_annotated.csv",      tissue), row.names = FALSE)
  write.csv(gene_summary, sprintf("gene_counts_summary_%s.csv",        tissue), row.names = FALSE)
  
  cat(sprintf("%s: annotated DEG table (%d rows), gene counts (%d genes) saved.\n",
              tissue, nrow(all_de), nrow(gene_summary)))
  
  return(list(all_de = all_de, gene_summary = gene_summary))
}

# Run for both tissues
pbmc_results <- process_deg_csv("all-pfs-PBMC-combined-data.csv", "PBMC")
til_results  <- process_deg_csv("all-pfs-TIL-combined-data.csv",  "TIL")


# =============================================================================
# PART 2: Recurring DEG stacked bar plot
# (equivalent to the Python plotting script)
# =============================================================================

# --- Artifact filter ---
is_artifact <- function(gene) {
  histone_prefixes <- c("HIST", "H1F", "H2AF", "H3F", "H4F",
                         "H1-", "H2A", "H2B", "H3-", "H4-", "H1.")
  
  grepl("^ENSG", gene) |                                  # unannotated Ensembl IDs
  grepl("^MT-", gene) |                                   # mitochondrial genes
  grepl("^MTRN", gene) |                                  # mitochondrial RNA genes
  sapply(gene, function(g)                                # histone genes
    any(startsWith(g, histone_prefixes))) |
  grepl("^AC[0-9]+", gene) |                             # uncharacterized AC loci
  grepl("^LINC", gene) |                                  # lncRNAs
  grepl("-AS[0-9]+$", gene)                               # antisense transcripts
}

# --- Plot function ---
make_recurring_deg_plot <- function(label, all_de, contrast_levels, out_prefix) {
  
  # Filter artifacts
  all_de <- all_de %>% filter(!is_artifact(Gene))
  
  # Top 50 genes by count, ascending for coord_flip effect
  gg_tbl <- all_de %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts) %>%
    slice_head(n = 50) %>%
    arrange(counts)
  
  top_genes <- gg_tbl$Gene
  
  # Filter all_de to top genes, set factor level order
  temp <- all_de %>%
    filter(Gene %in% top_genes) %>%
    mutate(Gene = factor(Gene, levels = top_genes),
           cell.contrasts = factor(cell.contrasts, levels = contrast_levels))
  
  # Colors: Blues (Baseline) + YlOrRd (Week4), 9 each = 18 total
  cols <- c(brewer.pal(9, "Blues")[2:9],   # skip lightest
            brewer.pal(9, "YlOrRd")[2:9])
  # Pad to 18 if needed
  cols <- c(brewer.pal(9, "Blues")  %>% colorRampPalette()(9),
            brewer.pal(9, "YlOrRd") %>% colorRampPalette()(9))

  # Plot
  p <- ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = cols,
      labels = str_remove(contrast_levels,
                          paste0(label, "\\.PFS\\.6month\\.NoPrg6M\\.vs\\.Prg6M\\.")) %>%
               str_replace_all("\\.", " ")
    ) +
    scale_y_continuous(limits = c(0, 19)) +
    coord_flip() +
    ylab("Number of times found") +
    xlab("") +
    ggtitle(
      sprintf("Top 50 Recurring DEGs — %s", label),
      subtitle = "(ENS, MT, Histone, AC, LINC, -AS genes excluded)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      plot.subtitle   = element_text(size = 9, color = "grey40")
    ) +
    guides(fill = guide_legend(title = "Contrast", ncol = 1))
  
  # Save
  ggsave(sprintf("%s.png", out_prefix), p, width = 9,  height = 14, dpi = 150)
  ggsave(sprintf("%s.pdf", out_prefix), p, width = 9,  height = 14)
  cat(sprintf("%s plot saved.\n", label))
  
  return(p)
}

# --- Contrast levels ---
pbmc_contrasts <- c(
  "Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg"
)

# TIL contrasts: same but swap PBMC → TIL
til_contrasts <- str_replace(pbmc_contrasts, "PBMC", "TIL")

# --- Run plots ---
p_pbmc <- make_recurring_deg_plot(
  label          = "PBMC",
  all_de         = pbmc_results$all_de,
  contrast_levels = pbmc_contrasts,
  out_prefix     = "recurring_DEG_PBMC_final2"
)

p_til <- make_recurring_deg_plot(
  label          = "TIL",
  all_de         = til_results$all_de,
  contrast_levels = til_contrasts,
  out_prefix     = "recurring_DEG_TIL_final2"
)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# --- Artifact filter ---
is_artifact <- function(gene) {
  histone_prefixes <- c("HIST", "H1F", "H2AF", "H3F", "H4F",
                        "H1-", "H2A", "H2B", "H3-", "H4-", "H1.")
  grepl("^ENSG", gene)                                          |
  grepl("^MT-",  gene)                                          |
  grepl("^MTRN", gene)                                          |
  sapply(gene, function(g) any(startsWith(g, histone_prefixes))) |
  grepl("^AC[0-9]+", gene)                                      |
  grepl("^LINC",     gene)                                      |
  grepl("-AS[0-9]+$", gene)
}

# --- Plot function ---
make_recurring_deg_plot <- function(label, all_de, contrast_levels, out_prefix) {
  
  # Filter artifacts
  all_de <- all_de %>% filter(!is_artifact(Gene))
  
  # Top 50 genes by count, ascending for coord_flip effect
  gg_tbl <- all_de %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts) %>%
    slice_head(n = 50) %>%
    arrange(counts)
  
  top_genes <- gg_tbl$Gene
  
  # Filter to top genes and set factor levels
  temp <- all_de %>%
    filter(Gene %in% top_genes) %>%
    mutate(
      Gene          = factor(Gene, levels = top_genes),
      cell.contrasts = factor(cell.contrasts, levels = contrast_levels)
    )
  
  # --- Dominant direction per gene ---
  dominant_dir <- temp %>%
    group_by(Gene, summary) %>%
    reframe(n = n()) %>%
    arrange(-n) %>%
    distinct(Gene, .keep_all = TRUE) %>%       # keep top direction per gene
    select(Gene, dominant = summary)
  
  # Merge direction back onto gene count table for annotation
  gg_tbl <- gg_tbl %>%
    left_join(dominant_dir, by = "Gene") %>%
    mutate(
      dir_label = ifelse(grepl("non.progressor", dominant), "▲ Non-prog", "▲ Prog"),
      dir_color = ifelse(grepl("non.progressor", dominant), "#4575b4",    "#d73027")
    )
  
  # Colors: Blues (Baseline) + YlOrRd (Week4), 9 shades each
  cols <- c(colorRampPalette(brewer.pal(9, "Blues"))(9),
            colorRampPalette(brewer.pal(9, "YlOrRd"))(9))
  
  # Shortened legend labels
  short_labels <- contrast_levels %>%
    str_remove(paste0(label, "\\.PFS\\.6month\\.NoPrg6M\\.vs\\.Prg6M\\.")) %>%
    str_replace_all("\\.", " ")
  
  # Max count for x-axis limit
  max_count <- all_de %>%
    filter(Gene %in% top_genes) %>%
    group_by(Gene) %>%
    reframe(n = n()) %>%
    pull(n) %>%
    max()
  
  # --- Base plot ---
  p <- ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cols, labels = short_labels) +
    scale_y_continuous(limits = c(0, max_count + 3)) +
    coord_flip() +
    ylab("Number of times found") +
    xlab("") +
    ggtitle(
      sprintf("Top 50 Recurring DEGs — %s", label),
      subtitle = "(ENS, MT, Histone, AC, LINC, -AS genes excluded)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.title    = element_text(size = 8),
      legend.text     = element_text(size = 7),
      legend.key.size = unit(0.4, "cm"),
      plot.subtitle   = element_text(size = 9, color = "grey40")
    ) +
    guides(fill = guide_legend(title = "Contrast", ncol = 1))
  
  # --- Direction annotations ---
  # geom_text with per-gene color requires a separate column and scale_color_identity
  p <- p +
    geom_text(
      data        = gg_tbl,
      aes(x       = Gene,
          y       = counts + 0.2,       # just past end of bar
          label   = dir_label,
          color   = dir_color),
      inherit.aes = FALSE,
      hjust       = 0,
      size        = 2.5,
      fontface    = "bold"
    ) +
    scale_color_identity(
      guide  = guide_legend(title = "Dominant direction", order = 2),
      labels = c("#4575b4" = "▲ Up in Non-progressor",
                 "#d73027" = "▲ Up in Progressor"),
      breaks = c("#4575b4", "#d73027")
    )
  
  # Save
  ggsave(sprintf("%s.png", out_prefix), p, width = 10, height = 14, dpi = 150)
  ggsave(sprintf("%s.pdf", out_prefix), p, width = 10, height = 14)
  cat(sprintf("%s plot saved.\n", label))
  
  return(p)
}

# --- Contrast levels ---
pbmc_contrasts <- c(
  "Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Baseline.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Baseline.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Baseline.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Week4.combo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Week4.nivo.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg",
  "Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.CD8plus.T.cell",
  "Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Macrophages",
  "Week4.rela.PBMC.PFS.6month.NoPrg6M.vs.Prg6M.Treg"
)

til_contrasts <- str_replace(pbmc_contrasts, "PBMC", "TIL")

# --- Run ---
p_pbmc <- make_recurring_deg_plot(
  label           = "PBMC",
  all_de          = pbmc_results$all_de,   # from Part 1 processing script
  contrast_levels = pbmc_contrasts,
  out_prefix      = "recurring_DEG_PBMC_combined_direction"
)

p_til <- make_recurring_deg_plot(
  label           = "TIL",
  all_de          = til_results$all_de,    # from Part 1 processing script
  contrast_levels = til_contrasts,
  out_prefix      = "recurring_DEG_TIL_combined_direction"
)