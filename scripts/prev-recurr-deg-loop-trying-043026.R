#make deg-list obj
library(dplyr)
library(readr)

# Set path to your folder
csv_folder <- "DEG_contrasts_PBMC_wbgf_MAST_v11/filtered_csv"

# Get all CSV file paths
csv_files <- list.files(path = csv_folder, pattern = "\\.csv$", full.names = TRUE)

# Load all into a named list (names = filename without extension)
deg_list <- lapply(csv_files, function(f) {
  df <- read_csv(f)
  df <- as.data.frame(df)
  rownames(df) <- df[[1]]   # set first column as rownames (gene names)
  df <- df[, -1]            # remove gene column
  return(df)
}) %>%
  setNames(tools::file_path_sans_ext(basename(csv_files)))


filter_uncharacterized_genes <- function(deg_list) {
  
  remove_patterns <- c(
    # Ensembl IDs
    "^ENSG",
    "^ENS",
    
    # LOC / LINC
    "^LOC",
    "^LINC",
    
    # Antisense / overlapping / divergent / intronic transcripts
    "-AS[0-9]*$",
    "-OT[0-9]*",
    "-DT$",
    "-IT[0-9]*",
    
    # ── Accession-style genes (the ones slipping through) ──────────────────
    "^AC[0-9]+\\.[0-9]+",     # AC233755.2, AC119396.1
    "^AL[0-9]+\\.[0-9]+",     # AL606867.1, AL034397.3
    "^AP[0-9]+\\.[0-9]+",     # AP003328.1, AP003774.3
    "^AF[0-9]+\\.[0-9]+",     # AF117879.1
    "^AJ[0-9]+\\.[0-9]+",     # AJ-style accessions
    "^BX[0-9]+\\.[0-9]+",     # BX-style accessions
    "^CR[0-9]+\\.[0-9]+",     # CR-style accessions
    "^CT[0-9]+\\.[0-9]+",     # CT-style accessions
    "^CU[0-9]+\\.[0-9]+",     # CU-style accessions
    "^FP[0-9]+\\.[0-9]+",     # FP-style accessions
    "^Z[0-9]+\\.[0-9]+",      # Z-style accessions
    
    # Broad catch — any pattern: LETTERS followed by digits dot digits
    "^[A-Z]{2}[0-9]+\\.[0-9]+",   # catches AC, AL, AP, AF, BX etc.
    
    # KIAA genes (uncharacterized)
    "^KIAA[0-9]+",
    
    # Histones
    "^HIST[0-9]",
    "^H[1-4][A-Z]",
    "^H[1-4]-",
    
    # Mitochondrial
    "^MT-",
    "^MTRNR",
    "^MTRN",
    
    # FAM genes (family with sequence similarity — often uncharacterized)
    "^FAM[0-9]+",
    
    # C#orf genes (chromosome open reading frames)
    "^C[0-9]+orf[0-9]+"
  )
  
  combined_pattern <- paste(remove_patterns, collapse = "|")
  
  filtered_list <- lapply(deg_list, function(df) {
    genes <- rownames(df)
    keep  <- !grepl(combined_pattern, genes, perl = TRUE)
    df[keep, ]
  })
  
  return(filtered_list)
}

deg_list_filtered <- filter_uncharacterized_genes(deg_list)

#BEST PARSER
extract_contrast_annotation <- function(contrast_names, cat) {
  
  df <- data.frame(row.names = contrast_names)
  
  # ── Timepoint:  ─────
  df$Timepoint <- case_when(
      str_detect(contrast_names, "\\.Baseline\\.") ~ "Baseline",
      str_detect(contrast_names, "\\.Week4\\.")    ~ "Week4",
      TRUE ~ "Unknown"
    )
  
  # ── CellType: shared across all ───────────────────────────────────────────
  df$CellType <- case_when(
    str_detect(contrast_names, "\\.CD3plus\\.")     ~ "CD3+",
    str_detect(contrast_names, "\\.CD4\\.")         ~ "CD4",
    str_detect(contrast_names, "\\.Macrophages\\.") ~ "Macrophages",
    TRUE ~ "Unknown"
  )
  
  # ── Category-specific ─────────────────────────────────────────────────────
  df$Comparison <- case_when(
      str_detect(contrast_names, "\\.combo\\.vs\\.nivo\\.") ~ "Combo vs Nivo",
      str_detect(contrast_names, "\\.combo\\.vs\\.rela\\.") ~ "Combo vs Rela",
      str_detect(contrast_names, "\\.rela\\.vs\\.nivo\\.")  ~ "Rela vs Nivo",
      TRUE ~ "Unknown"
    )
    
  #df$Treatment <- case_when(
  #    str_detect(contrast_names, "\\.combo\\.") ~ "Combo",
  #    str_detect(contrast_names, "\\.nivo\\.")  ~ "Nivo",
  #    str_detect(contrast_names, "\\.rela\\.")  ~ "Rela",
   #   TRUE ~ "Unknown"
   # )
    
  return(df)
}

#best palette
palette_list <- list(
  Direction = c(
    "Consistent UP"   = "#d73027",
    "Consistent DOWN" = "#4575b4",
    "Mixed"           = "#fee090"
  ),
  Timepoint = c(
    "Baseline" = "grey90",
    "Week4"    = "hotpink"
  ),
  CellType = c(
    "CD3+"        = "gold",
    "CD4"         = "violet",
    "Macrophages" = "darkcyan"
  ),
  Comparison = c(
    "Combo vs Nivo" = "blue",
    "Combo vs Rela" = "darkmagenta",
    "Rela vs Nivo"  = "grey10"
  ),
  Treatment = c(
    "Combo" = "magenta",
    "Nivo"  = "cyan",
    "Rela"  = "grey60"
  ),
  Response = c(
    "NoPrg6M vs Prg6M" = "coral"   # only one level — single colour is fine
  )
)

library(dplyr)
library(tibble)
library(purrr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# ─────────────────────────────────────────────
# STEP 1: Build log2FC Table from deg_list_filtered
# ─────────────────────────────────────────────

build_log2fc_table <- function(deg_list_filtered, 
                               log2fc_col = "avg_log2FC",
                               min_contrasts = 2) {
  
  # Extract gene + log2FC per contrast
  log2fc_table <- deg_list_filtered %>%
    imap(function(df, contrast_name) {
      df %>%
        rownames_to_column("gene") %>%
        dplyr::select(gene, all_of(log2fc_col)) %>%
        dplyr::rename(!!contrast_name := all_of(log2fc_col))
    }) %>%
    reduce(full_join, by = "gene") %>%
    column_to_rownames("gene")
  
  # Filter to recurrent genes only (present in >= min_contrasts)
  recurrent <- log2fc_table[rowSums(!is.na(log2fc_table)) >= min_contrasts, ]
  
  cat("Total genes before recurrence filter:", nrow(log2fc_table), "\n")
  cat("Recurrent genes (>=", min_contrasts, "contrasts):", nrow(recurrent), "\n")
  
  return(recurrent)
}

# Run it
log2fc_mat <- build_log2fc_table(
  deg_list_filtered,
  log2fc_col    = "avg_log2FC",  # change to "log2FoldChange" for DESeq2
  min_contrasts = 2
)

head(log2fc_mat)
