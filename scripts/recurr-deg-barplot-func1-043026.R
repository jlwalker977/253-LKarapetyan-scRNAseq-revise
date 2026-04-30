library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# =============================================================================
# SHARED UTILITIES
# =============================================================================

# --- Artifact filter ---------------------------------------------------------
is_artifact <- function(gene) {
  histone_prefixes <- c("HIST","H1F","H2AF","H3F","H4F",
                        "H1-","H2A","H2B","H3-","H4-","H1.")
  grepl("^ENSG",      gene) |
    grepl("^MT-",       gene) |
    grepl("^MTRN",      gene) |
    sapply(gene, function(g) any(startsWith(g, histone_prefixes))) |
    grepl("^AC[0-9]+",  gene) |
    grepl("^LINC",      gene) |
    grepl("-AS[0-9]+$", gene)
}

# ------ make combined csv files for each sample type # ------
#files.til.v11 = list.files(path = "DEG_contrasts_TIL_wbgf_MAST_v11/filtered_csv", pattern="*\\.csv$", full.names = TRUE)
#combined_data.til.v11 <- files.til.v11 %>%
#  set_names() %>% # Uses filenames as identifiers
#  map_dfr(read_csv, .id = "source_file") # Adds 'source_file' column
#files.pbmc.v11 = list.files(path = "DEG_contrasts_PBMC_wbgf_MAST_v11/filtered_csv", pattern="*\\.csv$", full.names = TRUE)
#combined_data.pbmc.v11 <- files.pbmc.v11 %>%
#  set_names() %>% # Uses filenames as identifiers
#  map_dfr(read_csv, .id = "source_file") # Adds 'source_file' column
#write.csv(combined_data.pbmc.v11, file="allcats-pbmc-v11.csv")
#write.csv(combined_data.til.v11, file="allcats-til-v11.csv")

# --- Extract cell type from source_file path ---------------------------------
extract_cell_type <- function(src) {
  m <- str_match(src, "\\.([^.]+)\\.DEG\\.")
  ifelse(is.na(m[,2]), NA, m[,2])
}

# --- Classify contrast into comparison category ------------------------------
get_category <- function(contrast) {
  case_when(
    grepl("PFS\\.6month", contrast) ~ "PFS",
    grepl("treatment",    contrast) ~ "Treatment",
    grepl("timepoint",    contrast) ~ "Timepoint",
    TRUE                            ~ "Other"
  )
}

# --- Direction summary label -------------------------------------------------
make_summary <- function(contrast, fc) {
  case_when(
    grepl("NoPrg6M\\.vs\\.Prg6M", contrast) ~
      ifelse(fc > 0, "up.in.non.progressor", "up.in.progressor"),
    grepl("treatment|timepoint", contrast) ~ {
      m <- str_match(contrast, "(?:treatment|timepoint)\\.(\\w+)\\.vs\\.(\\w+)")
      ifelse(fc > 0, paste0("up.in.", m[,2]), paste0("up.in.", m[,3]))
    },
    TRUE ~ ifelse(fc > 0, "up.in.group1", "up.in.group2")
  )
}

# --- Direction colour + short label ------------------------------------------
dir_color <- function(d) {
  ifelse(grepl("non\\.progressor|Week4|combo", d), "#4575b4", "#d73027")
}

dir_label_short <- function(d) {
  d %>%
    str_replace("up\\.in\\.", "\u25b2 ") %>%
    str_replace_all("\\.", " ") %>%
    str_replace("non progressor", "Non-prog") %>%
    str_replace("progressor",     "Prog")
}

# --- Shorten contrast label for legend ---------------------------------------
shorten_label <- function(contrast, category) {
  contrast %>%
    str_replace("\\.(PBMC|TIL)\\.", ".") %>%
    str_replace("PFS\\.6month\\.", "") %>%
    str_replace("treatment\\.",    "") %>%
    str_replace("timepoint\\.",    "") %>%
    str_replace("NoPrg6M\\.vs\\.Prg6M", "NonProg vs Prog") %>%
    str_replace_all("\\.", " ")
}

# --- Process CSV → annotated DEG table --------------------------------------
process_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$cell.type      <- extract_cell_type(df$source_file)
  df$category       <- get_category(df$contrast)
  df$cell.contrasts <- paste0(df$contrast, ".", df$cell.type)
  df$summary        <- make_summary(df$contrast, df$avg_log2FC)
  df
}

# --- Assign colours per contrast based on category --------------------------
#
#   PFS / Treatment : Blues  = Baseline contrasts
#                     Purples = Week4 contrasts
#   Timepoint       : Blues   = combo arm
#                     Greens  = nivo arm
#                     Oranges = rela arm
#
assign_colors <- function(contrast_levels, category) {
  
  make_shades <- function(palette_name, n) {
    colorRampPalette(brewer.pal(9, palette_name))(n + 2)[3:(n + 2)]  # skip lightest 2
  }
  
  color_map <- setNames(character(length(contrast_levels)), contrast_levels)
  
  if (category %in% c("PFS", "Treatment")) {
    
    baseline <- sort(contrast_levels[grepl("^Baseline", contrast_levels)])
    week4    <- sort(contrast_levels[grepl("^Week4",    contrast_levels)])
    
    if (length(baseline) > 0)
      color_map[baseline] <- make_shades("Blues",   length(baseline))
    if (length(week4) > 0)
      color_map[week4]    <- make_shades("Purples", length(week4))
    
    ordered <- c(baseline, week4)
    
  } else if (category == "Timepoint") {
    
    combo <- sort(contrast_levels[grepl("combo", contrast_levels)])
    nivo  <- sort(contrast_levels[grepl("nivo",  contrast_levels)])
    rela  <- sort(contrast_levels[grepl("rela",  contrast_levels)])
    
    if (length(combo) > 0) color_map[combo] <- make_shades("Blues",   length(combo))
    if (length(nivo)  > 0) color_map[nivo]  <- make_shades("Greens",  length(nivo))
    if (length(rela)  > 0) color_map[rela]  <- make_shades("Oranges", length(rela))
    
    ordered <- c(combo, nivo, rela)
    
  } else {
    
    ordered           <- sort(contrast_levels)
    color_map[ordered] <- make_shades("Greys", length(ordered))
    
  }
  
  list(ordered = ordered, color_map = color_map)
}

# --- Build legend group list -------------------------------------------------
build_legend_groups <- function(ordered_contrasts, color_map, category) {
  # Returns: list of list(title, data.frame(label, color))
  groups <- list()
  
  if (category %in% c("PFS", "Treatment")) {
    for (tp in list(c("^Baseline", "\u2500\u2500 Baseline \u2500\u2500"),
                    c("^Week4",    "\u2500\u2500 Week 4 \u2500\u2500"))) {
      idx <- grepl(tp[1], ordered_contrasts)
      if (any(idx)) {
        groups[[length(groups) + 1]] <- list(
          title = tp[2],
          items = data.frame(
            label = shorten_label(ordered_contrasts[idx], category),
            color = color_map[ordered_contrasts[idx]],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  } else if (category == "Timepoint") {
    for (arm in list(c("combo", "\u2500\u2500 Combo \u2500\u2500"),
                     c("nivo",  "\u2500\u2500 Nivo \u2500\u2500"),
                     c("rela",  "\u2500\u2500 Rela \u2500\u2500"))) {
      idx <- grepl(arm[1], ordered_contrasts)
      if (any(idx)) {
        groups[[length(groups) + 1]] <- list(
          title = arm[2],
          items = data.frame(
            label = shorten_label(ordered_contrasts[idx], category),
            color = color_map[ordered_contrasts[idx]],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  groups
}


# =============================================================================
# PLOT FUNCTION
# =============================================================================

make_plot <- function(label, sub_df, category, out_png, n_top = 50) {
  
  # Filter artifacts
  sub_df <- sub_df %>% filter(!is_artifact(Gene))
  
  # Colour assignment
  contrast_levels <- unique(sub_df$cell.contrasts)
  color_info      <- assign_colors(contrast_levels, category)
  ordered         <- color_info$ordered
  color_map       <- color_info$color_map
  
  # Top n genes, ascending for coord_flip
  gg_tbl <- sub_df %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts) %>%
    slice_head(n = n_top) %>%
    arrange(counts)
  
  top_genes <- gg_tbl$Gene
  
  temp <- sub_df %>%
    filter(Gene %in% top_genes) %>%
    mutate(
      Gene           = factor(Gene, levels = top_genes),
      cell.contrasts = factor(cell.contrasts, levels = ordered)
    )
  
  # Dominant direction per gene
  dominant_dir <- temp %>%
    group_by(Gene, summary) %>%
    reframe(n = n()) %>%
    arrange(-n) %>%
    distinct(Gene, .keep_all = TRUE) %>%
    select(Gene, dominant = summary)
  
  gg_tbl <- gg_tbl %>%
    left_join(dominant_dir, by = "Gene") %>%
    mutate(
      dlabel = dir_label_short(dominant),
      dcolor = dir_color(dominant)
    )
  
  max_count <- max(gg_tbl$counts)
  
  # Base plot
  p <- ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(
      values = color_map[ordered],
      breaks = ordered,
      guide  = "none"           # legend drawn manually below
    ) +
    scale_y_continuous(limits = c(0, max_count + 5)) +
    coord_flip() +
    ylab("Number of times found") + xlab("") +
    ggtitle(
      sprintf("Top %d Recurring DEGs \u2014 %s | %s", n_top, label, category),
      subtitle = "ENS, MT, Histone, AC, LINC, -AS excluded"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size =  9, color = "grey40"),
      legend.position = "none"
    )
  
  # Direction annotations
  p <- p +
    geom_text(
      data        = gg_tbl,
      aes(x       = Gene,
          y       = counts + 0.2,
          label   = dlabel,
          color   = dcolor),
      inherit.aes = FALSE,
      hjust       = 0,
      size        = 2.4,
      fontface    = "bold"
    ) +
    scale_color_identity()
  
  # --- Manual legends via annotation_custom / cowplot -------------------------
  # Build a combined legend as a separate grob using ggplot, then stitch
  library(cowplot)
  library(grid)
  
  legend_groups <- build_legend_groups(ordered, color_map, category)
  
  # Contrast legend grob
  legend_df <- bind_rows(lapply(legend_groups, function(g) {
    bind_rows(
      data.frame(label = g$title, color = NA,  is_header = TRUE,  stringsAsFactors = FALSE),
      data.frame(label = g$items$label, color = g$items$color, is_header = FALSE, stringsAsFactors = FALSE)
    )
  }))
  legend_df$y <- rev(seq_len(nrow(legend_df)))
  
  leg_plot <- ggplot(legend_df) +
    geom_tile(data = subset(legend_df, !is_header),
              aes(x = 0.2, y = y, fill = color), width = 0.3, height = 0.8) +
    scale_fill_identity() +
    geom_text(data = subset(legend_df, !is_header),
              aes(x = 0.4, y = y, label = label),
              hjust = 0, size = 2.5) +
    geom_text(data = subset(legend_df,  is_header),
              aes(x = 0.05, y = y, label = label),
              hjust = 0, size = 2.6, fontface = "bold", color = "grey30") +
    # Direction legend at bottom
    annotate("tile",   x = 0.2, y = -1.5, width = 0.3, height = 0.7, fill = "#4575b4") +
    annotate("text",   x = 0.4, y = -1.5, label = "\u25b2 Up in Non-prog / Week4 / Combo",
             hjust = 0, size = 2.5, color = "#4575b4", fontface = "bold") +
    annotate("tile",   x = 0.2, y = -2.5, width = 0.3, height = 0.7, fill = "#d73027") +
    annotate("text",   x = 0.4, y = -2.5, label = "\u25b2 Up in Prog / Baseline / Nivo / Rela",
             hjust = 0, size = 2.5, color = "#d73027", fontface = "bold") +
    xlim(0, 4) +
    theme_void()
  
  # Combine main plot + legend side by side
  combined <- plot_grid(p, leg_plot,
                        nrow  = 1,
                        rel_widths = c(3, 1.1))
  
  ggsave(out_png, combined, width = 13, height = 16, dpi = 150)
  message(out_png, " saved.")
  invisible(combined)
}


# =============================================================================
# RUN
# =============================================================================

pbmc <- process_csv("allcats-pbmc-v11.csv")
til  <- process_csv("allcats-til-v11.csv")

for (tissue_label in c("PBMC", "TIL")) {
  df <- if (tissue_label == "PBMC") pbmc else til
  
  for (cat in c("PFS", "Treatment", "Timepoint")) {
    sub <- df %>% filter(category == cat)
    make_plot(
      label    = tissue_label,
      sub_df   = sub,
      category = cat,
      out_png  = sprintf("recurring_DEG_%s_%s_v11_gradients.png", tissue_label, cat)
    )
  }
}