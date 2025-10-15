#!/usr/bin/env Rscript

# ============================================================
# Heatmap of GO BP enrichment across 4 sets (A=Col, B=cslc)
# Files:
#   GO_Aonly_Up.tsv
#   GO_Bonly_Down_BP.tsv
#   GO_Both_Down_BP.tsv
#   GO_Both_Up_BP.tsv
# Output: go_bp_heatmap_output/ (matrix + PNG/PDF)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)  # dplyr, tidyr, readr, ggplot2, tibble, etc.
})

cfg <- list(
  files        = c(
                   "GO_Bonly_Down.tsv",
                   "GO_Aonly_Up.tsv",
                   "GO_Both_Down.tsv",
                   "GO_Both_Up.tsv"),
  set_names    = c(
                   "cslc-only Down",
                   "Col-only Up",
                   "Both Down",
                   "Both Up"),
  outdir       = "go_bp_heatmap_output",
  value_metric = "neglog10_padj",  # or "GeneRatio"
  top_per_set  = 10,
  padj_cutoff  = 0.1,   # set NA to disable filtering
  cap_neglog10 = 50,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  drop_all_na_rows = TRUE
)

if (!dir.exists(cfg$outdir)) dir.create(cfg$outdir, recursive = TRUE)

# ---------------------------
# Helpers (explicit namespaces)
# ---------------------------
detect_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) NA_character_ else hit[1]
}

parse_gene_ratio <- function(x) {
  x <- as.character(x)
  ifelse(grepl("/", x),
         vapply(strsplit(x, "/"), function(v) {
           suppressWarnings(as.numeric(v[1]) / as.numeric(v[2]))
         }, numeric(1)),
         suppressWarnings(as.numeric(x)))
}

read_enrich <- function(path, set_label) {
  ext <- tolower(tools::file_ext(path))
  df <- if (ext %in% c("csv")) readr::read_csv(path, show_col_types = FALSE) else readr::read_tsv(path, show_col_types = FALSE)

  id_col   <- detect_col(df, c("ID","GO.ID","GOID","TermID"))
  desc_col <- detect_col(df, c("Description","Term","Name","Term Name"))
  padj_col <- detect_col(df, c("p.adjust","qvalue","FDR","adj.P.Val","padj"))
  p_col    <- detect_col(df, c("pvalue","P.value","p_value","p"))
  gr_col   <- detect_col(df, c("GeneRatio","Gene_Ratio","gene_ratio"))
  cnt_col  <- detect_col(df, c("Count","count","Gene.Count"))

  if (is.na(id_col) || is.na(desc_col)) {
    stop("Could not find GO ID or Description columns in file: ", path,
         "\nColumns available: ", paste(names(df), collapse = ", "))
  }
  if (is.na(padj_col)) {
    if (!is.na(p_col)) {
      warning("No adjusted p-value column in ", path, ". Falling back to raw p-values for plotting.")
      padj_col <- p_col
    } else {
      stop("Neither adjusted p-value nor p-value found in: ", path)
    }
  }

  tibble::tibble(
    set        = set_label,
    ID         = df[[id_col]],
    Description= df[[desc_col]],
    p_adjust   = suppressWarnings(as.numeric(df[[padj_col]])),
    pvalue     = if (!is.na(p_col)) suppressWarnings(as.numeric(df[[p_col]])) else NA_real_,
    GeneRatio  = if (!is.na(gr_col)) as.character(df[[gr_col]]) else NA_character_,
    Count      = if (!is.na(cnt_col)) suppressWarnings(as.numeric(df[[cnt_col]])) else NA_real_
  ) %>%
    dplyr::mutate(
      neglog10_padj = dplyr::if_else(is.na(p_adjust) | p_adjust <= 0, NA_real_, -log10(p_adjust)),
      neglog10_padj = pmin(neglog10_padj, cfg$cap_neglog10),
      GeneRatio_num = parse_gene_ratio(GeneRatio)
    ) %>%
    { if (!is.null(cfg$padj_cutoff) && !is.na(cfg$padj_cutoff))
        dplyr::filter(., !is.na(p_adjust) & p_adjust <= cfg$padj_cutoff) else . }
}

safe_order <- function(M, by_rows = TRUE, dist_method = "euclidean", hclust_method = "complete") {
  X <- if (by_rows) M else t(M)
  Xi <- X
  Xi[is.na(Xi)] <- 0  # impute NA for ordering only
  nobj <- if (by_rows) nrow(Xi) else ncol(Xi)
  if (nobj <= 1) return(seq_len(nobj))
  d <- stats::dist(Xi, method = dist_method)
  hc <- stats::hclust(d, method = hclust_method)
  hc$order
}

# ---------------------------
# Read all 4 files
# ---------------------------
stopifnot(length(cfg$files) == 4, length(cfg$set_names) == 4)

enrich_list <- purrr::map2(cfg$files, cfg$set_names, ~ read_enrich(.x, .y))

# ---------------------------
# Select top terms per file & combine
# ---------------------------
metric_col <- if (cfg$value_metric == "GeneRatio") "GeneRatio_num" else "neglog10_padj"

tops <- enrich_list %>%
  purrr::map(~ .x %>% dplyr::arrange(dplyr::desc(.data[[metric_col]])) %>% dplyr::slice_head(n = cfg$top_per_set))

selected_ids <- unique(unlist(purrr::map(tops, ~ .x$ID)))

combined <- dplyr::bind_rows(enrich_list) %>%
  dplyr::filter(ID %in% selected_ids) %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(Description = Description[which.max(nchar(tidyr::replace_na(Description, "")))]) %>%
  dplyr::ungroup()

# Long -> wide matrix
mat_df <- combined %>%
  dplyr::select(set, ID, Description, value = dplyr::all_of(metric_col)) %>%
  dplyr::mutate(set = factor(set, levels = cfg$set_names)) %>%
  dplyr::distinct(set, ID, .keep_all = TRUE)

row_labels_tbl <- mat_df %>%
  dplyr::select(ID, Description) %>%
  dplyr::distinct() %>%
  dplyr::mutate(label = ifelse(is.na(Description) | Description == "", ID,
                        paste0(Description, " (", ID, ")")))

mat_wide <- mat_df %>%
  dplyr::select(ID, set, value) %>%
  tidyr::pivot_wider(names_from = set, values_from = value) %>%
  dplyr::mutate(row_order_tmp = apply(dplyr::across(dplyr::all_of(cfg$set_names)), 1, function(x) max(x, na.rm = TRUE))) %>%
  dplyr::arrange(dplyr::desc(row_order_tmp)) %>%
  dplyr::select(-row_order_tmp)

# Drop rows that are NA across all sets (optional)
if (isTRUE(cfg$drop_all_na_rows)) {
  mat_wide <- mat_wide %>%
    dplyr::mutate(nonempty = rowSums(dplyr::across(dplyr::all_of(cfg$set_names), ~ !is.na(.))) > 0) %>%
    dplyr::filter(nonempty) %>% dplyr::select(-nonempty)
}

if (nrow(mat_wide) == 0) {
  stop("No terms left to plot after filtering. Try lowering 'padj_cutoff' or increasing 'top_per_set'.")
}

mat <- mat_wide %>% tibble::column_to_rownames("ID") %>% as.matrix()
mat[is.infinite(mat)] <- NA

row_labels <- row_labels_tbl %>%
  dplyr::filter(ID %in% rownames(mat)) %>%
  dplyr::arrange(match(ID, rownames(mat))) %>%
  dplyr::pull(label)

# ---------------------------
# Compute safe ordering to avoid NA in hclust
# ---------------------------
row_ord <- if (isTRUE(cfg$cluster_rows)) safe_order(mat, by_rows = TRUE)  else seq_len(nrow(mat))
col_ord <- if (isTRUE(cfg$cluster_cols)) safe_order(mat, by_rows = FALSE) else seq_len(ncol(mat))

mat_plot <- mat[row_ord, col_ord, drop = FALSE]
row_labels_plot <- row_labels[row_ord]
colnames_ordered <- colnames(mat)[col_ord]

# ---------------------------
# Save combined matrix (ordered)
# ---------------------------
combined_out <- as.data.frame(mat_plot) %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::left_join(row_labels_tbl, by = "ID") %>%
  dplyr::relocate(label, .before = ID)

readr::write_tsv(combined_out, file.path(cfg$outdir, "go_bp_combined_matrix.tsv"))

# ---------------------------
# Heatmap plotting
# ---------------------------
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pal <- grDevices::colorRampPalette(c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6",
                                       "#4292c6","#2171b5","#08519c","#08306b"))(100)
  main_title <- paste0(
    if (cfg$value_metric == "GeneRatio") "GO BP enrichment (GeneRatio) — "
    else "GO BP enrichment (-log10 adj p) — ",
    "A = Col, B = cslc"
  )

  # Preview (on-screen)
  pheatmap::pheatmap(
    mat_plot,
    color = pal,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = FALSE,
    main = main_title,
    na_col = "grey90",
    legend = TRUE,
    border_color = NA,
    angle_col = 0,
    labels_row = row_labels_plot
  )

  # PNG with dynamic height
  grDevices::png(file.path(cfg$outdir, "go_bp_enrichment_heatmap2.png"),
                 width = 2000, height = max(1600, 20 + 18 * nrow(mat_plot)), res = 200)
  pheatmap::pheatmap(
    mat_plot,
    color = pal,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = FALSE,
    main = main_title,
    na_col = "grey90",
    legend = TRUE,
    border_color = NA,
    angle_col = 0,
    labels_row = row_labels_plot
  )
  grDevices::dev.off()

  # PDF with height scaled to number of terms
  grDevices::pdf(file.path(cfg$outdir, "go_bp_enrichment_heatmap2.pdf"),
                 width = 9, height = max(6, 0.3 * nrow(mat_plot)))
  pheatmap::pheatmap(
    mat_plot,
    color = pal,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = FALSE,
    main = main_title,
    na_col = "grey90",
    legend = TRUE,
    border_color = NA,
    angle_col = 0,
    labels_row = row_labels_plot
  )
  grDevices::dev.off()

} else {
  message("Package 'pheatmap' not installed; using ggplot2 tile heatmap.")

  tidy_mat <- mat_plot %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ID") %>%
    tidyr::pivot_longer(-ID, names_to = "set", values_to = "value") %>%
    dplyr::left_join(row_labels_tbl, by = "ID") %>%
    dplyr::mutate(
      set = factor(set, levels = colnames_ordered),
      label = factor(label, levels = rev(row_labels_plot))
    )

  p <- ggplot2::ggplot(tidy_mat, ggplot2::aes(x = set, y = label, fill = value)) +
    ggplot2::geom_tile(color = NA) +
    ggplot2::scale_fill_gradientn(
      colors = c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6",
                 "#4292c6","#2171b5","#08519c","#08306b"),
      na.value = "grey90"
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0(
        if (cfg$value_metric == "GeneRatio") "GO BP enrichment (GeneRatio) — "
        else "GO BP enrichment (-log10 adj p) — ",
        "A = Col, B = cslc"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank()
    )

  ggplot2::ggsave(file.path(cfg$outdir, "go_bp_enrichment_heatmap2.png"), p,
         width = 9, height = max(6, 0.3 * nrow(mat_plot)), dpi = 300)
  ggplot2::ggsave(file.path(cfg$outdir, "go_bp_enrichment_heatmap2.pdf"), p,
         width = 9, height = max(6, 0.3 * nrow(mat_plot)))
}

message("Done. Outputs in: ", normalizePath(cfg$outdir))