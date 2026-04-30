library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# =============================================================================
# SHARED UTILITIES
# =============================================================================

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

extract_cell_type <- function(src) {
  m <- str_match(src, "\\.([^.]+)\\.DEG\\.")
  ifelse(is.na(m[,2]), NA, m[,2])
}

get_category <- function(contrast) {
  case_when(
    grepl("PFS\\.6month", contrast) ~ "PFS",
    grepl("treatment",    contrast) ~ "Treatment",
    grepl("timepoint",    contrast) ~ "Timepoint",
    TRUE                            ~ "Other"
  )
}

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

process_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$cell.type      <- extract_cell_type(df$source_file)
  df$category       <- get_category(df$contrast)
  df$cell.contrasts <- paste0(df$contrast, ".", df$cell.type)
  df$summary        <- make_summary(df$contrast, df$avg_log2FC)
  df
}

# --- Colour assignment (on base contrasts, i.e. without cell type suffix) ----
#   PFS / Treatment : Blues = Baseline, Purples = Week4
#   Timepoint       : Blues = combo, Greens = nivo, Oranges = rela
assign_colors <- function(base_contrasts, category) {
  
  make_shades <- function(palette_name, n) {
    colorRampPalette(brewer.pal(9, palette_name))(n + 2)[3:(n + 2)]
  }
  
  color_map <- setNames(character(length(base_contrasts)), base_contrasts)
  
  if (category %in% c("PFS", "Treatment")) {
    baseline <- sort(base_contrasts[grepl("^Baseline", base_contrasts)])
    week4    <- sort(base_contrasts[grepl("^Week4",    base_contrasts)])
    if (length(baseline) > 0) color_map[baseline] <- make_shades("Blues",   length(baseline))
    if (length(week4)    > 0) color_map[week4]    <- make_shades("Purples", length(week4))
    ordered <- c(baseline, week4)
    
  } else if (category == "Timepoint") {
    combo <- sort(base_contrasts[grepl("combo", base_contrasts)])
    nivo  <- sort(base_contrasts[grepl("nivo",  base_contrasts)])
    rela  <- sort(base_contrasts[grepl("rela",  base_contrasts)])
    if (length(combo) > 0) color_map[combo] <- make_shades("Blues",   length(combo))
    if (length(nivo)  > 0) color_map[nivo]  <- make_shades("Greens",  length(nivo))
    if (length(rela)  > 0) color_map[rela]  <- make_shades("Oranges", length(rela))
    ordered <- c(combo, nivo, rela)
    
  } else {
    ordered            <- sort(base_contrasts)
    color_map[ordered] <- make_shades("Greys", length(ordered))
  }
  
  list(ordered = ordered, color_map = color_map)
}

# --- Shorten contrast label for legend (strips cell type + tissue + category)-
shorten_contrast <- function(contrast) {
  contrast %>%
    str_remove("\\.(CD3plus|CD4|Macrophages)$") %>%
    str_replace("\\.(PBMC|TIL)\\.", ".") %>%
    str_replace("PFS\\.6month\\.", "") %>%
    str_replace("treatment\\.",    "") %>%
    str_replace("timepoint\\.",    "") %>%
    str_replace("NoPrg6M\\.vs\\.Prg6M", "NonProg vs Prog") %>%
    str_replace_all("\\.", " ")
}

# --- Build grouped legend entries --------------------------------------------
build_legend_groups <- function(ordered_base, color_map, category) {
  groups <- list()
  
  if (category %in% c("PFS", "Treatment")) {
    for (grp in list(c("^Baseline", "\u2500\u2500 Baseline \u2500\u2500"),
                     c("^Week4",    "\u2500\u2500 Week 4 \u2500\u2500"))) {
      idx <- grepl(grp[1], ordered_base)
      if (any(idx)) {
        groups[[length(groups) + 1]] <- list(
          title = grp[2],
          items = data.frame(
            label = shorten_contrast(ordered_base[idx]),
            color = color_map[ordered_base[idx]],
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
  } else if (category == "Timepoint") {
    for (grp in list(c("combo", "\u2500\u2500 Combo \u2500\u2500"),
                     c("nivo",  "\u2500\u2500 Nivo \u2500\u2500"),
                     c("rela",  "\u2500\u2500 Rela \u2500\u2500"))) {
      idx <- grepl(grp[1], ordered_base)
      if (any(idx)) {
        groups[[length(groups) + 1]] <- list(
          title = grp[2],
          items = data.frame(
            label = shorten_contrast(ordered_base[idx]),
            color = color_map[ordered_base[idx]],
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  groups
}


# =============================================================================
# SINGLE PANEL PLOT (one cell type)
# =============================================================================

make_panel <- function(cell_df, cell_type, ordered_full, color_map_full,
                       max_x, n_top = 25) {
  
  # Top n genes for this cell type
  gg_tbl <- cell_df %>%
    group_by(Gene) %>%
    reframe(counts = n()) %>%
    arrange(-counts) %>%
    slice_head(n = n_top) %>%
    arrange(counts)
  
  top_genes <- gg_tbl$Gene
  
  temp <- cell_df %>%
    filter(Gene %in% top_genes) %>%
    mutate(
      Gene           = factor(Gene, levels = top_genes),
      cell.contrasts = factor(cell.contrasts,
                              levels = ordered_full[ordered_full %in% unique(cell_df$cell.contrasts)])
    )
  
  # Dominant direction
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
  
  # Only keep contrasts present for this cell type
  present_contrasts <- intersect(ordered_full, unique(cell_df$cell.contrasts))
  fill_values       <- color_map_full[present_contrasts]
  fill_labels       <- shorten_contrast(present_contrasts)
  
  p <- ggplot(temp, aes(x = Gene, y = 1, fill = cell.contrasts)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = fill_values,
                      labels = fill_labels,
                      breaks = present_contrasts,
                      guide  = "none") +           # legend handled separately
    scale_y_continuous(limits = c(0, max_x)) +
    coord_flip() +
    ylab("Number of times found") + xlab("") +
    ggtitle(cell_type) +
    theme_classic(base_size = 10) +
    theme(plot.title      = element_text(size = 11, face = "bold"),
          axis.text.y     = element_text(size = 7.5),
          axis.title.x    = element_text(size = 9),
          legend.position = "none") +
    # Direction annotations
    geom_text(data        = gg_tbl,
              aes(x       = Gene,
                  y       = counts + 0.1,
                  label   = dlabel,
                  color   = dcolor),
              inherit.aes = FALSE,
              hjust       = 0,
              size        = 2.2,
              fontface    = "bold") +
    scale_color_identity()
  
  p
}


# =============================================================================
# LEGEND PLOT (drawn as a ggplot, combined via patchwork)
# =============================================================================

make_legend_plot <- function(legend_groups) {
  
  # Stack all items into one data frame with group headers
  legend_df <- bind_rows(lapply(legend_groups, function(g) {
    bind_rows(
      data.frame(label = g$title,        color = NA,           is_header = TRUE,
                 stringsAsFactors = FALSE),
      data.frame(label = g$items$label,  color = g$items$color, is_header = FALSE,
                 stringsAsFactors = FALSE)
    )
  }))
  
  # Add direction rows at the bottom
  legend_df <- bind_rows(
    legend_df,
    data.frame(label = "\u2500\u2500 Dominant direction \u2500\u2500",
               color = NA, is_header = TRUE,  stringsAsFactors = FALSE),
    data.frame(label = "\u25b2 Up in Non-prog / Week4 / Combo",
               color = "#4575b4", is_header = FALSE, stringsAsFactors = FALSE),
    data.frame(label = "\u25b2 Up in Prog / Baseline / Nivo / Rela",
               color = "#d73027", is_header = FALSE, stringsAsFactors = FALSE)
  )
  
  legend_df$y <- rev(seq_len(nrow(legend_df)))
  
  ggplot(legend_df) +
    # Colour swatches for non-headers
    geom_tile(data    = subset(legend_df, !is_header),
              aes(x   = 0.2, y = y, fill = color),
              width   = 0.28, height = 0.75) +
    scale_fill_identity() +
    # Labels for non-headers — colour matches swatch for direction rows
    geom_text(data        = subset(legend_df, !is_header),
              aes(x       = 0.38, y = y, label = label,
                  color   = ifelse(grepl("\u25b2", label), color, "grey20")),
              hjust = 0, size = 2.5) +
    scale_color_identity() +
    # Bold headers
    geom_text(data     = subset(legend_df, is_header),
              aes(x    = 0.05, y = y, label = label),
              hjust    = 0, size = 2.6, fontface = "bold", color = "grey30") +
    xlim(0, 4.5) +
    theme_void()
}


# =============================================================================
# MAIN PLOT FUNCTION — combines panels + legend via patchwork
# =============================================================================

make_celltype_plot <- function(label, sub_df, category, out_png, n_top = 25) {
  
  sub_df     <- sub_df %>% filter(!is_artifact(Gene))
  cell_types <- sort(unique(na.omit(sub_df$cell.type)))
  
  # Build colour map from base contrasts (without cell type suffix)
  all_contrasts  <- unique(sub_df$cell.contrasts)
  base_contrasts <- unique(str_remove(all_contrasts,
                                      "\\.(CD3plus|CD4|Macrophages)$"))
  color_info     <- assign_colors(base_contrasts, category)
  ordered_base   <- color_info$ordered
  color_map_base <- color_info$color_map
  
  # Expand to full contrast × cell type combos, preserving order
  ordered_full <- unlist(lapply(ordered_base, function(bc) {
    paste0(bc, ".", cell_types)[paste0(bc, ".", cell_types) %in% all_contrasts]
  }))
  color_map_full <- setNames(
    color_map_base[str_remove(ordered_full, "\\.(CD3plus|CD4|Macrophages)$")],
    ordered_full
  )
  
  # Shared x-axis limit across all panels
  max_x <- max(unlist(lapply(cell_types, function(ct) {
    sub_df %>% filter(cell.type == ct) %>%
      group_by(Gene) %>% reframe(n = n()) %>%
      arrange(-n) %>% slice_head(n = n_top) %>% pull(n) %>% sum()
  }))) + 3.5
  
  # Build one panel per cell type
  panels <- lapply(cell_types, function(ct) {
    make_panel(
      cell_df        = sub_df %>% filter(cell.type == ct),
      cell_type      = ct,
      ordered_full   = ordered_full,
      color_map_full = color_map_full,
      max_x          = max_x,
      n_top          = n_top
    )
  })
  
  # Build legend
  legend_groups <- build_legend_groups(ordered_base, color_map_base, category)
  leg_plot      <- make_legend_plot(legend_groups)
  
  # Combine: panels side by side, legend on right
  combined <- wrap_plots(c(panels, list(leg_plot)),
                         nrow       = 1,
                         widths     = c(rep(3, length(panels)), 1.2)) +
    plot_annotation(
      title    = sprintf("Top %d Recurring DEGs \u2014 %s | %s | by Cell Type",
                         n_top, label, category),
      subtitle = "ENS, MT, Histone, AC, LINC, -AS excluded",
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size =  9, color = "grey40")
      )
    )
  
  ggsave(out_png, combined,
         width  = 7 * length(cell_types) + 2,
         height = 15,
         dpi    = 150)
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
    make_celltype_plot(
      label   = tissue_label,
      sub_df  = sub,
      category = cat,
      out_png = sprintf("recurring_DEG_%s_%s_v11_bycelltype.png",
                        tissue_label, cat)
    )
  }
}