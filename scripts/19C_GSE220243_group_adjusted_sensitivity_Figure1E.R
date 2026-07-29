# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 19C_GSE220243_group_adjusted_sensitivity_Figure1E.R
# PURPOSE :
#   Test whether donor-level pathway–MMP13 relationships remain
#   after accounting for disease group and donor sensitivity.
#
# STATISTICAL DESIGN:
#   - One row = one independent donor.
#   - Partial Spearman correlation is calculated by correlating
#     residualized ranks after adjustment for Normal/OA group.
#   - P values are estimated by permuting MMP13 within group,
#     preserving the Normal/OA structure.
#   - Leave-one-donor-out (LODO) ranges assess sensitivity to
#     individual donors.
#
# IMPORTANT:
#   - This script does not repeat single-cell analysis.
#   - It reads the corrected 12-donor table from Script 19B.
#   - Figure 1E remains a candidate until these results are read.
# ============================================================

rm(list = ls())
gc()

set.seed(20260723)

# ------------------------------------------------------------
# 1. Project root
# ------------------------------------------------------------

project_root <- paste0(
  "D:/PhD FILE Ymelda/",
  "Publikasi Protein buatan/",
  "OA_MMP13_Target_Discovery"
)

if (!dir.exists(project_root)) {
  stop(
    "Project folder tidak ditemukan:\n",
    project_root
  )
}

setwd(project_root)

cat("\nWorking directory:\n")
print(getwd())

# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

required_packages <- c(
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "patchwork",
  "ggrepel",
  "openxlsx",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

still_missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(still_missing) > 0) {
  stop(
    "Paket berikut belum tersedia: ",
    paste(
      still_missing,
      collapse = ", "
    )
  )
}

library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggrepel)
library(openxlsx)
library(scales)

# ------------------------------------------------------------
# 3. Input and output paths
# ------------------------------------------------------------

input_file <- file.path(
  "results",
  "tables",
  "GSE220243_corrected_whole_chondrocyte_donor_table.csv"
)

main_figure_folder <- file.path(
  "results",
  "figures",
  "Figure 1"
)

supplementary_figure_folder <- file.path(
  "results",
  "figures",
  "Figure 1 Suplementary",
  "S19_group_adjusted_sensitivity"
)

table_folder <- file.path(
  "results",
  "tables"
)

dir.create(
  main_figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  supplementary_figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Corrected donor table tidak ditemukan:\n",
    input_file,
    "\nJalankan Script 19B terlebih dahulu."
  )
}

# ------------------------------------------------------------
# 4. Load and validate the 12-donor table
# ------------------------------------------------------------

donor_table <- utils::read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "sample_label",
  "donor",
  "geo_accession",
  "group",
  "MMP13_CPM",
  "log2_MMP13_CPM",
  "NFkB_module_score",
  "Inflammatory_module_score",
  "p38_MAPK_module_score",
  "ECM_degradation_module_score"
)

missing_columns <- setdiff(
  required_columns,
  colnames(donor_table)
)

if (length(missing_columns) > 0) {
  stop(
    "Kolom berikut tidak ditemukan: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}

donor_table$group <- factor(
  donor_table$group,
  levels = c(
    "Normal",
    "OA"
  )
)

if (
  nrow(donor_table) != 12 ||
  dplyr::n_distinct(
    donor_table$donor
  ) != 12
) {
  stop(
    "Input harus berisi tepat 12 donor independen."
  )
}

score_columns <- c(
  "NFkB_module_score",
  "p38_MAPK_module_score",
  "Inflammatory_module_score",
  "ECM_degradation_module_score"
)

pathway_labels <- c(
  "NF-κB",
  "p38 MAPK",
  "Inflammatory response",
  "ECM degradation"
)

# ------------------------------------------------------------
# 5. Statistical helper functions
# ------------------------------------------------------------

safe_spearman <- function(x, y) {

  keep <- is.finite(x) &
    is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 5 ||
    length(unique(x)) < 3 ||
    length(unique(y)) < 2
  ) {
    return(
      list(
        n = length(x),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  current_test <- suppressWarnings(
    stats::cor.test(
      x = x,
      y = y,
      method = "spearman",
      exact = FALSE
    )
  )

  list(
    n = length(x),
    rho = unname(
      current_test$estimate
    ),
    p_value = current_test$p.value
  )
}

partial_spearman_group <- function(
  x,
  y,
  group,
  permutations = 10000,
  seed = 20260723
) {

  keep <- is.finite(x) &
    is.finite(y) &
    !is.na(group)

  x <- x[keep]
  y <- y[keep]
  group <- droplevels(
    factor(
      group[keep]
    )
  )

  if (
    length(x) < 8 ||
    length(unique(x)) < 3 ||
    length(unique(y)) < 2 ||
    nlevels(group) < 2
  ) {
    return(
      list(
        n = length(x),
        partial_rho = NA_real_,
        permutation_p = NA_real_
      )
    )
  }

  rank_x <- rank(
    x,
    ties.method = "average"
  )

  rank_y <- rank(
    y,
    ties.method = "average"
  )

  residual_x <- stats::residuals(
    stats::lm(
      rank_x ~ group
    )
  )

  residual_y <- stats::residuals(
    stats::lm(
      rank_y ~ group
    )
  )

  observed_rho <- stats::cor(
    residual_x,
    residual_y,
    method = "pearson",
    use = "complete.obs"
  )

  set.seed(seed)

  permuted_rho <- replicate(
    permutations,
    {
      permuted_y <- y

      for (
        current_group in
        levels(group)
      ) {

        current_index <- which(
          group == current_group
        )

        permuted_y[current_index] <-
          sample(
            y[current_index],
            size = length(current_index),
            replace = FALSE
          )
      }

      rank_permuted_y <- rank(
        permuted_y,
        ties.method = "average"
      )

      residual_permuted_y <- stats::residuals(
        stats::lm(
          rank_permuted_y ~ group
        )
      )

      stats::cor(
        residual_x,
        residual_permuted_y,
        method = "pearson",
        use = "complete.obs"
      )
    }
  )

  permuted_rho <- permuted_rho[
    is.finite(permuted_rho)
  ]

  permutation_p <- (
    1 +
      sum(
        abs(permuted_rho) >=
          abs(observed_rho)
      )
  ) /
    (
      1 +
        length(permuted_rho)
    )

  list(
    n = length(x),
    partial_rho = observed_rho,
    permutation_p = permutation_p
  )
}

leave_one_donor_out <- function(
  x,
  y,
  donor
) {

  rho_values <- vapply(
    seq_along(donor),
    function(current_index) {

      keep <- seq_along(donor) !=
        current_index

      current_result <- safe_spearman(
        x = x[keep],
        y = y[keep]
      )

      current_result$rho
    },
    FUN.VALUE = numeric(1)
  )

  tibble::tibble(
    omitted_donor = donor,
    LODO_rho = rho_values
  )
}

# ------------------------------------------------------------
# 6. Raw, group-adjusted, within-group, and LODO analyses
# ------------------------------------------------------------

analysis_list <- list()
lodo_list <- list()
within_group_list <- list()

for (
  current_index in
  seq_along(
    score_columns
  )
) {

  current_score <- score_columns[
    current_index
  ]

  current_pathway <- pathway_labels[
    current_index
  ]

  x_values <- donor_table[[
    current_score
  ]]

  y_values <- donor_table$
    log2_MMP13_CPM

  raw_result <- safe_spearman(
    x = x_values,
    y = y_values
  )

  adjusted_result <- partial_spearman_group(
    x = x_values,
    y = y_values,
    group = donor_table$group,
    permutations = 10000,
    seed = 20260723 +
      current_index
  )

  current_lodo <- leave_one_donor_out(
    x = x_values,
    y = y_values,
    donor = donor_table$sample_label
  ) %>%
    dplyr::mutate(
      pathway = current_pathway,
      score_column = current_score
    )

  lodo_list[[
    current_score
  ]] <- current_lodo

  analysis_list[[
    current_score
  ]] <- tibble::tibble(

    pathway =
      current_pathway,

    score_column =
      current_score,

    n_donors =
      raw_result$n,

    raw_spearman_rho =
      raw_result$rho,

    raw_p_value =
      raw_result$p_value,

    group_adjusted_partial_rho =
      adjusted_result$partial_rho,

    stratified_permutation_p =
      adjusted_result$permutation_p,

    LODO_min_rho =
      min(
        current_lodo$LODO_rho,
        na.rm = TRUE
      ),

    LODO_max_rho =
      max(
        current_lodo$LODO_rho,
        na.rm = TRUE
      ),

    LODO_median_rho =
      stats::median(
        current_lodo$LODO_rho,
        na.rm = TRUE
      )
  )

  for (
    current_group in
    levels(
      donor_table$group
    )
  ) {

    group_index <- donor_table$group ==
      current_group

    group_result <- safe_spearman(
      x = x_values[group_index],
      y = y_values[group_index]
    )

    within_group_list[[
      paste(
        current_score,
        current_group,
        sep = "__"
      )
    ]] <- tibble::tibble(

      pathway =
        current_pathway,

      group =
        current_group,

      n_donors =
        group_result$n,

      spearman_rho =
        group_result$rho,

      p_value =
        group_result$p_value
    )
  }
}

sensitivity_summary <- dplyr::bind_rows(
  analysis_list
) %>%
  dplyr::mutate(

    raw_FDR =
      stats::p.adjust(
        raw_p_value,
        method = "BH"
      ),

    adjusted_FDR =
      stats::p.adjust(
        stratified_permutation_p,
        method = "BH"
      ),

    adjusted_evidence =
      dplyr::case_when(
        adjusted_FDR < 0.05 ~
          "FDR < 0.05",

        stratified_permutation_p <
          0.05 ~
          "Nominal p < 0.05",

        TRUE ~
          "Not significant"
      ),

    pathway = factor(
      pathway,
      levels = rev(
        pathway_labels
      )
    )
  )

lodo_results <- dplyr::bind_rows(
  lodo_list
)

within_group_results <- dplyr::bind_rows(
  within_group_list
) %>%
  dplyr::mutate(
    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )

cat("\nGroup-adjusted sensitivity summary:\n")
print(
  sensitivity_summary,
  n = Inf,
  width = Inf
)

cat("\nWithin-group descriptive correlations:\n")
print(
  within_group_results,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 7. Select the strongest group-adjusted relationship
# ------------------------------------------------------------

best_pathway_row <- sensitivity_summary %>%
  dplyr::filter(
    is.finite(
      stratified_permutation_p
    )
  ) %>%
  dplyr::arrange(
    adjusted_FDR,
    stratified_permutation_p,
    dplyr::desc(
      abs(
        group_adjusted_partial_rho
      )
    )
  ) %>%
  dplyr::slice(1)

if (nrow(best_pathway_row) == 0) {
  stop(
    "Tidak ada group-adjusted pathway result."
  )
}

best_score_column <- as.character(
  best_pathway_row$score_column
)

best_pathway_name <- as.character(
  best_pathway_row$pathway
)

best_raw_rho <- best_pathway_row$
  raw_spearman_rho

best_raw_p <- best_pathway_row$
  raw_p_value

best_partial_rho <- best_pathway_row$
  group_adjusted_partial_rho

best_permutation_p <- best_pathway_row$
  stratified_permutation_p

best_adjusted_fdr <- best_pathway_row$
  adjusted_FDR

# ------------------------------------------------------------
# 8. Figure 1E candidate — sensitivity summary
# ------------------------------------------------------------

evidence_palette <- c(
  "FDR < 0.05" = "#B2182B",
  "Nominal p < 0.05" = "#EF8A62",
  "Not significant" = "#6E6E6E"
)

sensitivity_plot <- ggplot(
  sensitivity_summary,
  aes(
    x = raw_spearman_rho,
    y = pathway,
    color = adjusted_evidence
  )
) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.45,
    linetype = "dashed",
    color = "grey55"
  ) +
  geom_errorbarh(
    aes(
      xmin = LODO_min_rho,
      xmax = LODO_max_rho
    ),
    height = 0.16,
    linewidth = 0.70
  ) +
  geom_point(
    size = 3.2
  ) +
  geom_point(
    aes(
      x = group_adjusted_partial_rho
    ),
    shape = 21,
    size = 2.9,
    stroke = 0.75,
    fill = "white"
  ) +
  scale_color_manual(
    values = evidence_palette,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(
      -1,
      1
    ),
    breaks = seq(
      -1,
      1,
      by = 0.5
    )
  ) +
  labs(
    title =
      "Donor-level pathway–MMP13 sensitivity analysis",

    subtitle =
      paste(
        "Filled points: raw Spearman ρ;",
        "open points: group-adjusted partial ρ;",
        "bars: leave-one-donor-out range"
      ),

    x =
      "Correlation with log2(MMP13 CPM + 1)",

    y =
      NULL,

    color =
      "Group-adjusted evidence"
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11.5
      ),

    plot.subtitle =
      element_text(
        size = 8.8
      ),

    legend.position =
      "bottom",

    legend.title =
      element_blank(),

    legend.text =
      element_text(
        size = 8
      )
  )

# ------------------------------------------------------------
# 9. Figure 1E candidate — raw donor scatter
#
# No regression line is drawn because half of the donors have
# zero MMP13 and the sample size is only 12.
# ------------------------------------------------------------

scatter_table <- donor_table %>%
  dplyr::mutate(
    pathway_score =
      .data[[
        best_score_column
      ]]
  )

scatter_annotation <- paste0(
  "Raw Spearman ρ = ",
  sprintf(
    "%.2f",
    best_raw_rho
  ),
  "; P = ",
  scales::pvalue(
    best_raw_p,
    accuracy = 0.001
  ),
  "\nGroup-adjusted partial ρ = ",
  sprintf(
    "%.2f",
    best_partial_rho
  ),
  "; stratified P = ",
  scales::pvalue(
    best_permutation_p,
    accuracy = 0.001
  ),
  "\nAdjusted FDR = ",
  scales::pvalue(
    best_adjusted_fdr,
    accuracy = 0.001
  )
)

best_pathway_scatter <- ggplot(
  scatter_table,
  aes(
    x = pathway_score,
    y = log2_MMP13_CPM,
    color = group,
    shape = group
  )
) +
  geom_point(
    size = 3.1,
    alpha = 0.92
  ) +
  ggrepel::geom_text_repel(
    aes(
      label = sample_label
    ),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.28,
    point.padding = 0.16,
    min.segment.length = 0,
    segment.color = "grey55",
    seed = 20260723
  ) +
  scale_color_manual(
    values = c(
      "Normal" = "#D55E00",
      "OA" = "#009E9E"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Normal" = 16,
      "OA" = 17
    )
  ) +
  annotate(
    geom = "label",
    x = -Inf,
    y = Inf,
    hjust = -0.04,
    vjust = 1.08,
    size = 2.8,
    label = scatter_annotation,
    label.size = 0.20,
    fill = "white"
  ) +
  labs(
    title =
      paste0(
        best_pathway_name,
        " and MMP13"
      ),

    subtitle =
      "Each point represents one independent cartilage donor",

    x =
      paste0(
        "Weighted mean ",
        best_pathway_name,
        " module score"
      ),

    y =
      expression(
        log[2] *
          "(MMP13 CPM + 1)"
      ),

    color =
      "Group",

    shape =
      "Group"
  ) +
  theme_classic(
    base_size = 10.5
  ) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11.5
      ),

    plot.subtitle =
      element_text(
        size = 9
      ),

    legend.position =
      "bottom"
  )

figure_1e_sensitivity_candidate <- (
  sensitivity_plot |
    best_pathway_scatter
) +
  patchwork::plot_layout(
    widths = c(
      1,
      1
    )
  )

# ------------------------------------------------------------
# 10. Supplementary within-group correlation heatmap
# ------------------------------------------------------------

within_group_heatmap <- within_group_results %>%
  dplyr::mutate(
    pathway = factor(
      pathway,
      levels = pathway_labels
    ),
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = pathway,
      y = group,
      fill = spearman_rho
    )
  ) +
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(
      label = dplyr::case_when(
        is.na(
          spearman_rho
        ) ~
          "NA",

        TRUE ~
          sprintf(
            "%.2f",
            spearman_rho
          )
      )
    ),
    size = 3.2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    na.value = "grey90",
    name = "Spearman ρ"
  ) +
  labs(
    title =
      "Within-group donor-level pathway–MMP13 correlations",

    subtitle =
      "Descriptive only; each group contains six donors",

    x =
      "Pathway program",

    y =
      NULL
  ) +
  theme_bw(
    base_size = 10.5
  ) +
  theme(
    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(
        angle = 35,
        hjust = 1
      )
  )

# ------------------------------------------------------------
# 11. Save figures
# ------------------------------------------------------------

main_output_base <- file.path(
  main_figure_folder,
  paste0(
    "Figure1E_CANDIDATE_group_adjusted_",
    "sensitivity_pathway_MMP13"
  )
)

ggsave(
  filename = paste0(
    main_output_base,
    ".pdf"
  ),
  plot = figure_1e_sensitivity_candidate,
  width = 11.5,
  height = 5.4,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = paste0(
    main_output_base,
    ".tiff"
  ),
  plot = figure_1e_sensitivity_candidate,
  width = 11.5,
  height = 5.4,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

ggsave(
  filename = paste0(
    main_output_base,
    ".png"
  ),
  plot = figure_1e_sensitivity_candidate,
  width = 11.5,
  height = 5.4,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary_figure_folder,
    "FigureS19A_within_group_correlation_heatmap.pdf"
  ),
  plot = within_group_heatmap,
  width = 8.5,
  height = 4.5,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary_figure_folder,
    "FigureS19A_within_group_correlation_heatmap.tiff"
  ),
  plot = within_group_heatmap,
  width = 8.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

# ------------------------------------------------------------
# 12. Save tables and caption draft
# ------------------------------------------------------------

utils::write.csv(
  sensitivity_summary,
  file = file.path(
    table_folder,
    "GSE220243_group_adjusted_pathway_MMP13_sensitivity.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  lodo_results,
  file = file.path(
    table_folder,
    "GSE220243_pathway_MMP13_leave_one_donor_out.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  within_group_results,
  file = file.path(
    table_folder,
    "GSE220243_pathway_MMP13_within_group_correlations.csv"
  ),
  row.names = FALSE
)

summary_workbook <- file.path(
  table_folder,
  "GSE220243_group_adjusted_Figure1E_sensitivity_summary.xlsx"
)

openxlsx::write.xlsx(
  list(
    Donor_data =
      as.data.frame(
        donor_table
      ),

    Sensitivity_summary =
      as.data.frame(
        sensitivity_summary
      ),

    Leave_one_donor_out =
      as.data.frame(
        lodo_results
      ),

    Within_group =
      as.data.frame(
        within_group_results
      )
  ),
  file = summary_workbook,
  overwrite = TRUE
)

caption_text <- paste(
  "Figure 1E. Sensitivity analysis of donor-level pathway–MMP13",
  "relationships. Filled points show unadjusted Spearman correlations",
  "across 12 independent cartilage donors, open points show partial",
  "Spearman correlations after adjustment for Normal/OA group, and",
  "horizontal bars show the range obtained after sequentially omitting",
  "one donor. The donor scatter displays the pathway with the strongest",
  "group-adjusted evidence. Because MMP13 was undetected in half of the",
  "donors, all associations were interpreted as exploratory."
)

writeLines(
  caption_text,
  con = file.path(
    table_folder,
    "GSE220243_Figure1E_group_adjusted_caption_draft.txt"
  )
)

# ------------------------------------------------------------
# 13. Final summary
# ------------------------------------------------------------

message("")
message("================================================")
message("SCRIPT 19C SELESAI")
message("================================================")

message(
  "Independent donors             : ",
  nrow(
    donor_table
  )
)

message(
  "Group-adjusted method          : ",
  "partial Spearman + within-group permutation"
)

message(
  "LODO sensitivity              : COMPLETE"
)

message(
  "Best group-adjusted pathway   : ",
  best_pathway_name
)

message(
  "Raw Spearman rho              : ",
  sprintf(
    "%.3f",
    best_raw_rho
  )
)

message(
  "Group-adjusted partial rho     : ",
  sprintf(
    "%.3f",
    best_partial_rho
  )
)

message(
  "Stratified permutation P      : ",
  signif(
    best_permutation_p,
    4
  )
)

message(
  "Adjusted FDR                  : ",
  signif(
    best_adjusted_fdr,
    4
  )
)

message(
  "Figure 1E candidate TIFF      : ",
  paste0(
    main_output_base,
    ".tiff"
  )
)

message(
  "Figure 1E status              : CANDIDATE; inspect biological robustness"
)

message("================================================")
