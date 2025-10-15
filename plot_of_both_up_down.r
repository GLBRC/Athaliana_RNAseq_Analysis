# ===============================================================
# Summary plot: Both Up vs Both Down (Δlog2FC = cslc - col)
# - Violin + boxplot per group
# - Mean ± 95% CI overlay
# - Annotations: within-group (Δ vs 0) and between-group p-values
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# Ensure output dir
if (!exists("mag_outdir") || is.null(mag_outdir)) {
  mag_outdir <- if (exists("config") && !is.null(config$outdir)) {
    file.path(config$outdir, "magnitude_tests_both_updown")
  } else {
    "magnitude_tests_both_updown"
  }
}
if (!dir.exists(mag_outdir)) dir.create(mag_outdir, recursive = TRUE)

stopifnot(all(c("gene","log2FC_A","log2FC_B","state_A","state_B") %in% names(both)))

# Subsets (if not already defined)
if (!exists("both_up") || !exists("both_down")) {
  both_up   <- both  %>% filter(state_A == "Up",   state_B == "Up")
  both_down <- both  %>% filter(state_A == "Down", state_B == "Down")
}

# Combined Δ table
df_delta <- bind_rows(
  both_up   %>% transmute(group = "Both Up",   delta = log2FC_B - log2FC_A),
  both_down %>% transmute(group = "Both Down", delta = log2FC_B - log2FC_A)
) %>%
  filter(is.finite(delta)) %>%
  mutate(group = factor(group, levels = c("Both Up","Both Down")))

if (nrow(df_delta) == 0) stop("No data available for Both Up / Both Down to plot.")

# Summary stats per group (mean ± 95% CI)
summary_grp <- df_delta %>%
  group_by(group) %>%
  summarise(
    n   = n(),
    mean = mean(delta),
    sd   = sd(delta),
    se   = sd / sqrt(n),
    ci_lo = mean + qt(0.025, df = pmax(n - 1, 1)) * se,
    ci_hi = mean + qt(0.975, df = pmax(n - 1, 1)) * se,
    .groups = "drop"
  )

# Within-group tests (Δ vs 0)
p_within <- df_delta %>%
  group_by(group) %>%
  reframe(
    t_p = tryCatch(t.test(delta, mu = 0)$p.value,      error = function(e) NA_real_),
    w_p = tryCatch(wilcox.test(delta, mu = 0, exact=FALSE)$p.value, error = function(e) NA_real_)
  )

summary_grp <- left_join(summary_grp, p_within, by = "group")

# Between-group tests (Both Up vs Both Down)
tt_between <- tryCatch(t.test(delta ~ group, data = df_delta), error = function(e) NULL)
wt_between <- tryCatch(wilcox.test(delta ~ group, data = df_delta, exact = FALSE), error = function(e) NULL)
p_between_t <- if (!is.null(tt_between)) tt_between$p.value else NA_real_
p_between_w <- if (!is.null(wt_between)) wt_between$p.value else NA_real_

# Annotation positions
y_range <- range(df_delta$delta, na.rm = TRUE)
y_span  <- diff(y_range)
ann_y   <- y_range[2] + 0.05 * y_span

ann_within <- summary_grp %>%
  mutate(
    label = sprintf("t p=%.3g\nW p=%.3g", t_p, w_p),
    y = ann_y
  )

# ---- FIXED: add x = group when inherit.aes = FALSE ----
p_summary <- ggplot(df_delta, aes(x = group, y = delta, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9) +
  geom_hline(yintercept = 0, color = "red3", linetype = "dashed") +
  geom_point(
    data = summary_grp,
    aes(x = group, y = mean),
    color = "black", size = 2,
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = summary_grp,
    aes(x = group, ymin = ci_lo, ymax = ci_hi),
    width = 0.12, color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ann_within,
    aes(x = group, y = y, label = label),
    size = 3.2, vjust = 0
  ) +
  scale_fill_manual(values = c("Both Up" = "#e41a1c", "Both Down" = "#377eb8")) +
  labs(
    title    = expression(Delta~"log2FC (cslc - col): Both Up vs Both Down"),
    subtitle = sprintf("Between-group: t p=%.3g, Wilcoxon p=%.3g", p_between_t, p_between_w),
    x = NULL,
    y = expression(Delta~"log2FC (cslc - col)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 10)
  )

ggsave(file.path(mag_outdir, "summary_both_up_vs_both_down.png"), p_summary, width = 6.8, height = 5.4, dpi = 300)
ggsave(file.path(mag_outdir, "summary_both_up_vs_both_down.pdf"), p_summary, width = 6.8, height = 5.4)

# Save numeric summary table as well
readr::write_tsv(
  summary_grp %>%
    transmute(group, n, mean_delta = mean, ci95_lo = ci_lo, ci95_hi = ci_hi,
              t_pvalue = t_p, wilcox_pvalue = w_p),
  file.path(mag_outdir, "summary_both_up_vs_both_down_stats.tsv")
)

message("Summary plot saved to: ", normalizePath(file.path(mag_outdir, "summary_both_up_vs_both_down.png")))