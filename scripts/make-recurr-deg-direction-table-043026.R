library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

# =============================================================================
# SHARED UTILITIES
# (same as plotting script — keep in one source file if running together)
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

process_csv <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$cell.type      <- extract_cell_type(df$source_file)
  df$category       <- get_category(df$contrast)
  df$cell.contrasts <- paste0(df$contrast, ".", df$cell.type)
  df$summary        <- make_summary(df$contrast, df$avg_log2FC)
  df
}


# =============================================================================
# DIRECTIONALITY TABLE FUNCTION
# =============================================================================

make_direction_table <- function(sub_df, n_top = 50) {
  
  # Apply artifact filter
  sub_df <- sub_df %>% filter(!is_artifact(Gene))
  
  # Top n genes by total count — same selection as barplot
  top_genes <- sub_df %>%
    group_by(Gene) %>%
    reframe(total_count = n()) %>%
    arrange(-total_count) %>%
    slice_head(n = n_top)
  
  # Count appearances per gene per direction
  dir_counts <- sub_df %>%
    filter(Gene %in% top_genes$Gene) %>%
    group_by(Gene, summary) %>%
    reframe(n = n()) %>%
    pivot_wider(names_from  = summary,
                values_from = n,
                values_fill = 0)
  
  # Dominant direction = column with highest count per gene
  summary_cols <- setdiff(colnames(dir_counts), "Gene")
  
  dir_counts <- dir_counts %>%
    rowwise() %>%
    mutate(dominant_direction = summary_cols[which.max(c_across(all_of(summary_cols)))]) %>%
    ungroup()
  
  # Merge with total counts, sort, add rank
  tbl <- top_genes %>%
    left_join(dir_counts, by = "Gene") %>%
    arrange(-total_count) %>%
    mutate(rank = row_number()) %>%
    select(rank, Gene, total_count, everything())
  
  tbl
}


# =============================================================================
# RUN — build all 6 tables and write to one Excel file
# =============================================================================

pbmc <- process_csv("allcats-pbmc-v11.csv")
til  <- process_csv("allcats-til-v11.csv")

# Collect all sheets
sheets <- list()

for (tissue_label in c("PBMC", "TIL")) {
  df <- if (tissue_label == "PBMC") pbmc else til
  
  for (cat in c("PFS", "Treatment", "Timepoint")) {
    sub <- df %>% filter(category == cat)
    tbl <- make_direction_table(sub, n_top = 50)
    sheet_name <- paste0(tissue_label, "_", cat)
    sheets[[sheet_name]] <- tbl
    message(sheet_name, ": ", nrow(tbl), " genes, ",
            ncol(tbl) - 4, " direction column(s)")
  }
}

# Write to Excel — one sheet per tissue/category
write.xlsx(sheets, file = "DEG_directionality_tables_allcats_v11.xlsx",
           rowNames = FALSE)

message("Saved: DEG_directionality_tables_allcats_v11.xlsx")

#--- optional, for separate csvs ---
#for (name in names(sheets)) {
#  write.csv(sheets[[name]],
#            file      = paste0("DEG_directionality_", name, "_v11.csv"),
#            row.names = FALSE)
#}