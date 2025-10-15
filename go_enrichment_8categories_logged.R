#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot)
  library(tidyverse)
})

# Load gene-level comparison data
df <- readr::read_tsv("gene_level_comparison.tsv")

# Define 8-category classification
df <- df %>% mutate(category_8 = case_when(
  state_A == "Up" & state_B == "Up" ~ "Both Up",
  state_A == "Down" & state_B == "Down" ~ "Both Down",
  state_A == "Up" & state_B == "Down" ~ "Opposite (A Up / B Down)",
  state_A == "Down" & state_B == "Up" ~ "Opposite (A Down / B Up)",
  state_A == "Up" & state_B == "NotSig" ~ "A only Up",
  state_A == "Down" & state_B == "NotSig" ~ "A only Down",
  state_A == "NotSig" & state_B == "Up" ~ "B only Up",
  state_A == "NotSig" & state_B == "Down" ~ "B only Down",
  TRUE ~ "Neither"
))

# Filter out 'Neither'
df_filtered <- df %>% filter(category_8 != "Neither")

# Prepare universe and gene sets
universe_genes <- unique(df_filtered$gene)
gene_sets <- df_filtered %>%
  group_by(category_8) %>%
  summarize(genes = list(unique(gene)), .groups = "drop")

# Log gene counts per category
cat("Gene counts per category:\n")
for (i in seq_len(nrow(gene_sets))) {
  cat(sprintf("  %-25s : %d genes\n", gene_sets$category_8[i], length(gene_sets$genes[[i]])))
}

# GO enrichment parameters
go_params <- list(
  OrgDb = org.At.tair.db,
  keyType = "TAIR",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Create output directory
out_dir <- "GO_enrichment_results"
dir.create(out_dir, showWarnings = FALSE)

# Run enrichment and save results
for (i in seq_len(nrow(gene_sets))) {
  cat_name <- gene_sets$category_8[i]
  genes <- gene_sets$genes[[i]]
  cat(sprintf("\nRunning GO enrichment for category: %s (%d genes)\n", cat_name, length(genes)))

  ego <- tryCatch({
    enrichGO(
      gene = genes,
      universe = universe_genes,
      OrgDb = go_params$OrgDb,
      keyType = go_params$keyType,
      ont = go_params$ont,
      pvalueCutoff = go_params$pvalueCutoff,
      qvalueCutoff = go_params$qvalueCutoff,
      pAdjustMethod = go_params$pAdjustMethod,
      minGSSize = go_params$minGSSize,
      maxGSSize = go_params$maxGSSize,
      readable = go_params$readable
    )
  }, error = function(e) {
    message("GO enrichment failed for ", cat_name, ": ", e$message)
    NULL
  })

  if (!is.null(ego)) {
    result_df <- as.data.frame(ego)
    if (nrow(result_df) > 0) {
      readr::write_tsv(result_df, file.path(out_dir, paste0("GO_", gsub("[^A-Za-z0-9]", "_", cat_name), ".tsv")))
      p <- dotplot(ego, showCategory = 20) +
        ggtitle(paste("GO Enrichment -", cat_name)) +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, paste0("GO_", gsub("[^A-Za-z0-9]", "_", cat_name), "_dotplot.png")), p, width = 8, height = 6)
      cat("  ✓ Results and dotplot saved\n")
    } else {
      cat("  ⚠ No enriched GO terms found for ", cat_name, "\n")
    }
  }
}
