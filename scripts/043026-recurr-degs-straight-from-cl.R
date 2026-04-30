library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(purrr)

# =============================================================================
# SHARED UTILITIES
# =============================================================================

# --- Artifact filter ---
is_artifact <- function(gene) {
  histone_prefixes <- c("HIST","H1F","H2AF","H3F","H4F","H1-","H2A","H2B","H3-","H4-","H1.")
  grepl("^ENSG",      gene) |
    grepl("^MT-",       gene) |
    grepl("^MTRN",      gene) |
    sapply(gene, function(g) any(startsWith(g, histone_prefixes))) |
    grepl("^AC[0-9]+",  gene) |
    grepl("^LINC",      gene) |
    grepl("-AS[0-9]+$", gene)
}

# --- Extract cell type from source_file ---
files.til.v11 = list.files(path = "DEG_contrasts_TIL_wbgf_MAST_v11/filtered_csv/", pattern="*\\.csv$", full.names = TRUE)
combined_data.til.v11 <- files.til.v11 %>%
  set_names() %>% # Uses filenames as identifiers
  map_dfr(read_csv, .id = "source_file") # Adds 'source_file' column
files.pbmc.v11 = list.files(path = "DEG_contrasts_PBMC_wbgf_MAST_v11/filtered_csv/", pattern="*\\.csv$", full.names = TRUE)
combined_data.pbmc.v11 <- files.pbmc.v11 %>%
  set_names() %>% # Uses filenames as identifiers
  map_dfr(read_csv, .id = "source_file") # Adds 'source_file' column
write.csv(combined_data.pbmc.v11, file="allcats-pbmc-v11.csv")
write.csv(combined_data.til.v11, file="allcats-til-v11.csv")

extract_cell_type <- function(src) {
  m <- str_match(src, "\\.([^.]+)\\.DEG\\.")
  ifelse(is.na(m[,2]), NA, m[,2])
}

# --- Direction summary label from contrast + fold change ---
make_summary_label <- function(contrast, fc) {
  case_when(
    grepl("NoPrg6M\\.vs\\.Prg6M", contrast) ~
      ifelse(fc > 0, "up.in.non.progressor", "up.in.progressor"),
    grepl("treatment", contrast) ~ {
      m <- str_match(contrast, "treatment\\.(\\w+)\\.vs\\.(\\w+)")
      ifelse(fc > 0,
             paste0("up.in.", m[,2]),
             paste0("up.in.", m[,3]))
    },
    grepl("timepoint", contrast) ~ {
      m <- str_match(contrast, "timepoint\\.(\\w+)\\.vs\\.(\\w+)")
      ifelse(fc > 0,
             paste0("up.in.", m[,2]),
             paste0("up.in.", m[,3]))
    },
    TRUE ~ ifelse(fc > 0, "up.in.group1", "up.in.group2")
  )
}

# --- Process CSV → annotated DEG table ---
process_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$cell.type      <- extract_cell_type(df$source_file)
  df$cell.contrasts <- paste0(df$contrast, ".", df$cell.type)
  df$summary        <- make_summary_label(df$contrast, df$avg_log2FC)
  df
}

# --- Direction colour helper ---
dir_color <- function(summary_str) {
  ifelse(grepl("non\\.progressor|Week4", summary_str), "#4575b4", "#d73027")
}

dir_label <- function(summary_str) {
  summary_str %>%
    str_replace("up\\.in\\.", "\u25b2 ") %>%
    str_replace_all("\\.", " ") %>%
    str_replace("non progressor", "Non-prog") %>%
    str_replace("progressor",     "Prog")
}

# --- Build Blues + YlOrRd palette scaled to n contrasts ---
make_palette <- function(n) {
  half <- n %/% 2
  c(colorRampPalette(brewer.pal(9, "Blues"))(max(half, 1)),
    colorRampPalette(brewer.pal(9, "YlOrRd"))(max(n - half, 1)))
}


# =============================================================================
# VERSION 1: Combined cell types WITH artifact filter
# =============================================================================

make_combined_filtered <- function(title, df, out_prefix, n_top = 50) {
  
  df <- df %>% filter(!is_artifact(Gene))
  
  contrast_levels <- sort(unique(df$cell.contrasts))
  cols            <- make_palette(length(contrast_levels))
  short_labels    <- str_replace_all(contrast_levels, "\\.", " ")
  
  gg_tbl <- df %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts) %>%
    slice_head(n = n_top) %>%
    arrange(counts)
  
  top_genes <- gg_tbl$Gene
  
  temp <- df %>%
    filter(Gene %in% top_genes) %>%
    mutate(Gene          = factor(Gene, levels = top_genes),
           cell.contrasts = factor(cell.contrasts, levels = contrast_levels))
  
  dominant_dir <- temp %>%
    group_by(Gene, summary) %>%
    reframe(n = n()) %>%
    arrange(-n) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, dominant = summary)
  
  gg_tbl <- gg_tbl %>%
    left_join(dominant_dir, by = "Gene") %>%
    mutate(dlabel = dir_label(dominant),
           dcolor = dir_color(dominant))
  
  max_count <- max(gg_tbl$counts)
  
  p <- ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cols, labels = short_labels) +
    scale_y_continuous(limits = c(0, max_count + 4)) +
    coord_flip() +
    ylab("Number of times found") + xlab("") +
    ggtitle(title, subtitle = "(ENS, MT, Histone, AC, LINC, -AS excluded)") +
    theme_classic(base_size = 11) +
    theme(legend.title    = element_text(size = 8),
          legend.text     = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.subtitle   = element_text(size = 9, color = "grey40")) +
    guides(fill = guide_legend(title = "Contrast · Cell type", ncol = 1)) +
    geom_text(data        = gg_tbl,
              aes(x       = Gene,
                  y       = counts + 0.2,
                  label   = dlabel,
                  color   = dcolor),
              inherit.aes = FALSE, hjust = 0, size = 2.5, fontface = "bold") +
    scale_color_identity(
      guide  = guide_legend(title = "Dominant direction", order = 2),
      labels = c("#4575b4" = "\u25b2 Up in Non-prog / Week4",
                 "#d73027" = "\u25b2 Up in Prog / Baseline"),
      breaks = c("#4575b4", "#d73027")
    )
  
  ggsave(paste0(out_prefix, ".png"), p, width = 11, height = 14, dpi = 150)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 11, height = 14)
  message(out_prefix, " saved.")
  invisible(p)
}


# =============================================================================
# VERSION 2: Combined cell types WITHOUT artifact filter
# =============================================================================

make_combined_unfiltered <- function(title, df, out_prefix, n_top = 50) {
  
  # identical to VERSION 1 except no is_artifact() call and subtitle differs
  contrast_levels <- sort(unique(df$cell.contrasts))
  cols            <- make_palette(length(contrast_levels))
  short_labels    <- str_replace_all(contrast_levels, "\\.", " ")
  
  gg_tbl <- df %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts) %>%
    slice_head(n = n_top) %>%
    arrange(counts)
  
  top_genes <- gg_tbl$Gene
  
  temp <- df %>%
    filter(Gene %in% top_genes) %>%
    mutate(Gene           = factor(Gene, levels = top_genes),
           cell.contrasts = factor(cell.contrasts, levels = contrast_levels))
  
  dominant_dir <- temp %>%
    group_by(Gene, summary) %>%
    reframe(n = n()) %>%
    arrange(-n) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, dominant = summary)
  
  gg_tbl <- gg_tbl %>%
    left_join(dominant_dir, by = "Gene") %>%
    mutate(dlabel = dir_label(dominant),
           dcolor = dir_color(dominant))
  
  max_count <- max(gg_tbl$counts)
  
  p <- ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cols, labels = short_labels) +
    scale_y_continuous(limits = c(0, max_count + 4)) +
    coord_flip() +
    ylab("Number of times found") + xlab("") +
    ggtitle(title, subtitle = "(no artifact filter)") +
    theme_classic(base_size = 11) +
    theme(legend.title    = element_text(size = 8),
          legend.text     = element_text(size = 7),
          legend.key.size = unit(0.4, "cm"),
          plot.subtitle   = element_text(size = 9, color = "grey40")) +
    guides(fill = guide_legend(title = "Contrast · Cell type", ncol = 1)) +
    geom_text(data        = gg_tbl,
              aes(x       = Gene,
                  y       = counts + 0.2,
                  label   = dlabel,
                  color   = dcolor),
              inherit.aes = FALSE, hjust = 0, size = 2.5, fontface = "bold") +
    scale_color_identity(
      guide  = guide_legend(title = "Dominant direction", order = 2),
      labels = c("#4575b4" = "\u25b2 Up in Non-prog / Week4",
                 "#d73027" = "\u25b2 Up in Prog / Baseline"),
      breaks = c("#4575b4", "#d73027")
    )
  
  ggsave(paste0(out_prefix, ".png"), p, width = 11, height = 14, dpi = 150)
  ggsave(paste0(out_prefix, ".pdf"), p, width = 11, height = 14)
  message(out_prefix, " saved.")
  invisible(p)
}


# =============================================================================
# VERSION 3: Split by cell type WITH artifact filter
# =============================================================================

make_celltype_filtered <- function(title, df, out_prefix, n_top = 25) {
  
  df         <- df %>% filter(!is_artifact(Gene))
  cell_types <- sort(unique(na.omit(df$cell.type)))
  n_cells    <- length(cell_types)
  
  # Build one ggplot per cell type, then combine with patchwork
  library(patchwork)
  
  panels <- map(cell_types, function(cell) {
    
    sub <- df %>% filter(cell.type == cell)
    
    cell_contrasts <- sub %>%
      pull(cell.contrasts) %>% unique() %>%
      sort_by_timepoint()   # helper defined below
    
    cols         <- make_palette(length(cell_contrasts))
    short_labels <- str_remove(cell_contrasts, paste0("\\.", cell)) %>%
      str_replace_all("\\.", " ")
    
    gg_tbl <- sub %>%
      group_by(Gene) %>%
      reframe(counts = n()) %>%
      arrange(-counts) %>%
      slice_head(n = n_top) %>%
      arrange(counts)
    
    top_genes <- gg_tbl$Gene
    
    temp <- sub %>%
      filter(Gene %in% top_genes) %>%
      mutate(Gene           = factor(Gene, levels = top_genes),
             cell.contrasts = factor(cell.contrasts, levels = cell_contrasts))
    
    dominant_dir <- temp %>%
      group_by(Gene, summary) %>%
      reframe(n = n()) %>%
      arrange(-n) %>%
      distinct(Gene, .keep_all = TRUE) %>%
      select(Gene, dominant = summary)
    
    gg_tbl <- gg_tbl %>%
      left_join(dominant_dir, by = "Gene") %>%
      mutate(dlabel = dir_label(dominant),
             dcolor = dir_color(dominant))
    
    max_count <- max(gg_tbl$counts)
    
    ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = cols, labels = short_labels) +
      scale_y_continuous(limits = c(0, max_count + 3.5)) +
      coord_flip() +
      ylab("Number of times found") + xlab("") +
      ggtitle(cell) +
      theme_classic(base_size = 10) +
      theme(legend.title    = element_text(size = 7),
            legend.text     = element_text(size = 6.5),
            legend.key.size = unit(0.35, "cm")) +
      guides(fill = guide_legend(title = "Contrast", ncol = 1)) +
      geom_text(data        = gg_tbl,
                aes(x       = Gene,
                    y       = counts + 0.1,
                    label   = dlabel,
                    color   = dcolor),
                inherit.aes = FALSE, hjust = 0, size = 2.2, fontface = "bold") +
      scale_color_identity(guide = "none")
  })
  
  # Combine panels side by side
  combined <- wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = title,
      subtitle = "(ENS, MT, Histone, AC, LINC, -AS excluded — by cell type)",
      theme    = theme(plot.title    = element_text(size = 13, face = "bold"),
                       plot.subtitle = element_text(size =  9, color = "grey40"))
    )
  
  ggsave(paste0(out_prefix, ".png"), combined,
         width = 7 * n_cells, height = 14, dpi = 150)
  ggsave(paste0(out_prefix, ".pdf"), combined,
         width = 7 * n_cells, height = 14)
  message(out_prefix, " saved.")
  invisible(combined)
}

# Helper: sort contrasts Baseline first, then Week4, alphabetically within
sort_by_timepoint <- function(contrasts) {
  contrasts[order(
    ifelse(grepl("Baseline", contrasts), 0, 1),
    contrasts
  )]
}


# =============================================================================
# RUN ALL 8 FILES THROUGH ALL 3 VERSIONS
# =============================================================================

configs <- list(
  list(fname = "alltrt_pbmc_wgfilt_v6.csv",  title = "Treatment PBMC — with gene filter (v6)"),
  list(fname = "alltrt_pbmc_nogfilt_v6.csv", title = "Treatment PBMC — no gene filter (v6)"),
  list(fname = "alltime_pbmc_wgfilt_v6.csv", title = "Timepoint PBMC — with gene filter (v6)"),
  list(fname = "alltime_pbmc_nogfilt_v6.csv",title = "Timepoint PBMC — no gene filter (v6)"),
  list(fname = "allpfs_til_wgfilt_v6.csv",   title = "PFS TIL — with gene filter (v6)"),
  list(fname = "allpfs_til_nogfilt_v6.csv",  title = "PFS TIL — no gene filter (v6)"),
  list(fname = "allpfs_pbmc_wgfilt_v5.csv",  title = "PFS PBMC — with gene filter (v5)"),
  list(fname = "allpfs_pbmc_nogfilt_v6.csv", title = "PFS PBMC — no gene filter (v6)")
)

for (cfg in configs) {
  stem <- tools::file_path_sans_ext(cfg$fname)
  df   <- process_csv(cfg$fname)
  
  make_combined_filtered(
    title      = cfg$title,
    df         = df,
    out_prefix = paste0("recurring_DEG_", stem, "_artfilt")
  )
  
  make_combined_unfiltered(
    title      = cfg$title,
    df         = df,
    out_prefix = paste0("recurring_DEG_", stem, "_noartfilt")
  )
  
  make_celltype_filtered(
    title      = cfg$title,
    df         = df,
    out_prefix = paste0("recurring_DEG_", stem, "_bycelltype")
  )
}