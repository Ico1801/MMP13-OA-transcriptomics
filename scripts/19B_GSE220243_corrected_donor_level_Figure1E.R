# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 19B_GSE220243_corrected_donor_level_Figure1E.R
# PURPOSE :
#   Correct Figure 1E inference so that DONOR, not donor-state
#   rows or individual cells, is the independent biological unit.
#
# WORKFLOW:
#   1. Collapse all reliable chondrocyte states to one row/donor.
#   2. Calculate cell-number-weighted pathway scores per donor.
#   3. Calculate pseudobulk MMP13 CPM per donor.
#   4. Test donor-level Spearman associations across 12 donors.
#   5. Retain state-specific analyses only when each donor
#      contributes at most one value within that state.
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
  stop("Project folder tidak ditemukan:\n", project_root)
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
    paste(still_missing, collapse = ", ")
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
# 3. Input and output
# ------------------------------------------------------------

input_file <- file.path(
  "results",
  "tables",
  "GSE220243_NFkB_p38_donor_state_table.csv"
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
  "S18_corrected_donor_level"
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
    "Donor-state table tidak ditemukan:\n",
    input_file,
    "\nJalankan Script 19A terlebih dahulu."
  )
}

# ------------------------------------------------------------
# 4. Load and validate donor-state table
# ------------------------------------------------------------

message("Memuat donor-state table...")

donor_state_table <- utils::read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "sample_label",
  "donor",
  "geo_accession",
  "group",
  "cell_state_short",
  "cell_state",
  "n_cells",
  "total_RNA_UMI",
  "MMP13_total_UMI",
  "MMP13_positive_cells",
  "MMP13_CPM",
  "log2_MMP13_CPM",
  "reliable_state_sample",
  "NFkB_module_score_mean",
  "Inflammatory_module_score_mean",
  "p38_MAPK_module_score_mean",
  "ECM_degradation_module_score_mean"
)

missing_columns <- setdiff(
  required_columns,
  colnames(donor_state_table)
)

if (length(missing_columns) > 0) {
  stop(
    "Kolom berikut tidak ditemukan: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (is.character(donor_state_table$reliable_state_sample)) {
  donor_state_table$reliable_state_sample <-
    tolower(donor_state_table$reliable_state_sample) %in%
    c("true", "t", "1", "yes")
}

donor_state_table$group <- factor(
  donor_state_table$group,
  levels = c("Normal", "OA")
)

state_order <- c(
  "RegC",
  "EC",
  "RepC",
  "PreFC",
  "FC",
  "MTC-like",
  "HomC",
  "HTC",
  "Vascular/immune"
)

chondrocyte_states <- setdiff(
  state_order,
  "Vascular/immune"
)

donor_state_table$cell_state_short <- factor(
  donor_state_table$cell_state_short,
  levels = state_order
)

state_score_columns <- c(
  "NFkB_module_score_mean",
  "Inflammatory_module_score_mean",
  "p38_MAPK_module_score_mean",
  "ECM_degradation_module_score_mean"
)

donor_score_columns <- c(
  "NFkB_module_score",
  "Inflammatory_module_score",
  "p38_MAPK_module_score",
  "ECM_degradation_module_score"
)

pathway_labels <- c(
  "NF-κB",
  "Inflammatory response",
  "p38 MAPK",
  "ECM degradation"
)

donor_pathway_key <- tibble::tibble(
  score_column = donor_score_columns,
  pathway = pathway_labels
)

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------

safe_weighted_mean <- function(values, weights) {

  keep <- is.finite(values) &
    is.finite(weights) &
    weights > 0

  if (!any(keep)) {
    return(NA_real_)
  }

  stats::weighted.mean(
    x = values[keep],
    w = weights[keep],
    na.rm = TRUE
  )
}

safe_spearman <- function(x, y) {

  keep <- is.finite(x) &
    is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 6 ||
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
    rho = unname(current_test$estimate),
    p_value = current_test$p.value
  )
}

bootstrap_spearman_ci <- function(
  x,
  y,
  iterations = 2000,
  seed = 20260723
) {

  keep <- is.finite(x) &
    is.finite(y)

  x <- x[keep]
  y <- y[keep]

  if (
    length(x) < 6 ||
    length(unique(x)) < 3 ||
    length(unique(y)) < 2
  ) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  set.seed(seed)

  bootstrap_rho <- replicate(
    iterations,
    {
      sampled_index <- sample.int(
        n = length(x),
        size = length(x),
        replace = TRUE
      )

      suppressWarnings(
        stats::cor(
          x = x[sampled_index],
          y = y[sampled_index],
          method = "spearman",
          use = "complete.obs"
        )
      )
    }
  )

  bootstrap_rho <- bootstrap_rho[
    is.finite(bootstrap_rho)
  ]

  if (length(bootstrap_rho) < 100) {
    return(
      c(
        lower = NA_real_,
        upper = NA_real_
      )
    )
  }

  current_quantiles <- stats::quantile(
    bootstrap_rho,
    probabilities = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )

  c(
    lower = current_quantiles[1],
    upper = current_quantiles[2]
  )
}

# ------------------------------------------------------------
# 6. Collapse repeated donor-state rows to one row per donor
# ------------------------------------------------------------

reliable_chondrocyte_table <- donor_state_table %>%
  dplyr::filter(
    reliable_state_sample,
    cell_state_short %in%
      chondrocyte_states
  )

donor_level_table <- reliable_chondrocyte_table %>%
  dplyr::group_by(
    sample_label,
    donor,
    geo_accession,
    group
  ) %>%
  dplyr::summarise(

    chondrocyte_cells =
      sum(n_cells, na.rm = TRUE),

    total_RNA_UMI =
      sum(total_RNA_UMI, na.rm = TRUE),

    MMP13_total_UMI =
      sum(MMP13_total_UMI, na.rm = TRUE),

    MMP13_positive_cells =
      sum(MMP13_positive_cells, na.rm = TRUE),

    NFkB_module_score =
      safe_weighted_mean(
        NFkB_module_score_mean,
        n_cells
      ),

    Inflammatory_module_score =
      safe_weighted_mean(
        Inflammatory_module_score_mean,
        n_cells
      ),

    p38_MAPK_module_score =
      safe_weighted_mean(
        p38_MAPK_module_score_mean,
        n_cells
      ),

    ECM_degradation_module_score =
      safe_weighted_mean(
        ECM_degradation_module_score_mean,
        n_cells
      ),

    represented_states =
      dplyr::n_distinct(cell_state_short),

    .groups = "drop"
  ) %>%
  dplyr::mutate(

    MMP13_CPM =
      ifelse(
        total_RNA_UMI > 0,
        1e6 *
          MMP13_total_UMI /
          total_RNA_UMI,
        NA_real_
      ),

    log2_MMP13_CPM =
      log2(MMP13_CPM + 1)
  ) %>%
  dplyr::arrange(
    group,
    sample_label
  )

if (
  nrow(donor_level_table) != 12 ||
  dplyr::n_distinct(donor_level_table$donor) != 12
) {
  warning(
    "Expected 12 independent donors, found ",
    nrow(donor_level_table),
    " rows and ",
    dplyr::n_distinct(donor_level_table$donor),
    " unique donors."
  )
}

cat("\nCorrected donor-level table:\n")
print(
  donor_level_table,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 7. Correct global donor-level correlations
# ------------------------------------------------------------

global_correlation_list <- list()

for (current_index in seq_along(donor_score_columns)) {

  current_score <- donor_score_columns[current_index]

  current_result <- safe_spearman(
    x = donor_level_table[[current_score]],
    y = donor_level_table$log2_MMP13_CPM
  )

  current_ci <- bootstrap_spearman_ci(
    x = donor_level_table[[current_score]],
    y = donor_level_table$log2_MMP13_CPM,
    iterations = 2000,
    seed = 20260723 + current_index
  )

  global_correlation_list[[current_score]] <-
    tibble::tibble(

      score_column =
        current_score,

      pathway =
        pathway_labels[current_index],

      n_donors =
        current_result$n,

      spearman_rho =
        current_result$rho,

      bootstrap_CI_lower =
        unname(current_ci["lower"]),

      bootstrap_CI_upper =
        unname(current_ci["upper"]),

      p_value =
        current_result$p_value
    )
}

global_correlations <- dplyr::bind_rows(
  global_correlation_list
) %>%
  dplyr::mutate(

    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      ),

    evidence =
      dplyr::case_when(
        FDR < 0.05 ~
          "FDR < 0.05",

        p_value < 0.05 ~
          "Nominal p < 0.05",

        TRUE ~
          "Not significant"
      ),

    pathway = factor(
      pathway,
      levels = rev(
        c(
          "NF-κB",
          "p38 MAPK",
          "Inflammatory response",
          "ECM degradation"
        )
      )
    )
  )

cat("\nCorrected donor-level correlations:\n")
print(
  global_correlations,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 8. State-specific donor-level correlations
# ------------------------------------------------------------

state_correlation_list <- list()

for (current_state in chondrocyte_states) {

  current_state_table <- donor_state_table %>%
    dplyr::filter(
      reliable_state_sample,
      cell_state_short ==
        current_state
    )

  for (current_index in seq_along(state_score_columns)) {

    current_score <- state_score_columns[current_index]

    current_nonzero <- sum(
      current_state_table$MMP13_total_UMI > 0,
      na.rm = TRUE
    )

    current_result <- safe_spearman(
      x = current_state_table[[current_score]],
      y = current_state_table$log2_MMP13_CPM
    )

    if (current_nonzero < 3) {
      current_result$rho <- NA_real_
      current_result$p_value <- NA_real_
    }

    state_correlation_list[[
      paste(
        current_state,
        current_score,
        sep = "__"
      )
    ]] <- tibble::tibble(

      cell_state_short =
        current_state,

      score_column =
        current_score,

      pathway =
        pathway_labels[current_index],

      n_donors =
        current_result$n,

      nonzero_MMP13_donors =
        current_nonzero,

      spearman_rho =
        current_result$rho,

      p_value =
        current_result$p_value
    )
  }
}

state_correlations <- dplyr::bind_rows(
  state_correlation_list
) %>%
  dplyr::mutate(

    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      ),

    significance_label =
      dplyr::case_when(
        FDR < 0.05 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),

    cell_state_short = factor(
      cell_state_short,
      levels = rev(chondrocyte_states)
    ),

    pathway = factor(
      pathway,
      levels = c(
        "NF-κB",
        "p38 MAPK",
        "Inflammatory response",
        "ECM degradation"
      )
    )
  )

# ------------------------------------------------------------
# 9. Select strongest corrected donor-level relationship
# ------------------------------------------------------------

best_pathway_row <- global_correlations %>%
  dplyr::filter(
    is.finite(p_value)
  ) %>%
  dplyr::arrange(
    FDR,
    p_value,
    dplyr::desc(abs(spearman_rho))
  ) %>%
  dplyr::slice(1)

if (nrow(best_pathway_row) == 0) {
  stop(
    "Tidak ada pathway correlation yang dapat dipilih."
  )
}

best_score_column <- as.character(
  best_pathway_row$score_column
)

best_pathway_name <- as.character(
  best_pathway_row$pathway
)

best_rho <- best_pathway_row$spearman_rho
best_p <- best_pathway_row$p_value
best_fdr <- best_pathway_row$FDR

# ------------------------------------------------------------
# 10. Figure 1E candidate — correlation summary
# ------------------------------------------------------------

evidence_palette <- c(
  "FDR < 0.05" = "#B2182B",
  "Nominal p < 0.05" = "#EF8A62",
  "Not significant" = "#6E6E6E"
)

correlation_summary_plot <- ggplot(
  global_correlations,
  aes(
    x = spearman_rho,
    y = pathway,
    color = evidence
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
      xmin = bootstrap_CI_lower,
      xmax = bootstrap_CI_upper
    ),
    height = 0.16,
    linewidth = 0.65
  ) +
  geom_point(
    size = 3.1
  ) +
  scale_color_manual(
    values = evidence_palette,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.5)
  ) +
  labs(
    title =
      "Donor-level pathway associations with MMP13",

    subtitle =
      "One observation per donor; bootstrap 95% confidence intervals",

    x =
      "Spearman correlation with log2(MMP13 CPM + 1)",

    y =
      NULL,

    color =
      NULL
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
      element_text(size = 9),

    legend.position =
      "bottom",

    legend.text =
      element_text(size = 8)
  )

# ------------------------------------------------------------
# 11. Figure 1E candidate — strongest pathway scatter
# ------------------------------------------------------------

best_scatter_table <- donor_level_table %>%
  dplyr::mutate(
    pathway_score =
      .data[[best_score_column]]
  )

x_annotation <- min(
  best_scatter_table$pathway_score,
  na.rm = TRUE
)

y_annotation <- max(
  best_scatter_table$log2_MMP13_CPM,
  na.rm = TRUE
)

best_pathway_scatter <- ggplot(
  best_scatter_table,
  aes(
    x = pathway_score,
    y = log2_MMP13_CPM,
    color = group,
    shape = group
  )
) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.60,
    linetype = "dashed",
    color = "grey35"
  ) +
  geom_point(
    size = 3.0,
    alpha = 0.90
  ) +
  ggrepel::geom_text_repel(
    aes(label = sample_label),
    size = 2.7,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
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
    x = x_annotation,
    y = y_annotation,
    hjust = 0,
    vjust = 1,
    size = 3.0,
    label = paste0(
      "Spearman ρ = ",
      sprintf("%.2f", best_rho),
      "\nP = ",
      scales::pvalue(
        best_p,
        accuracy = 0.001
      ),
      "; FDR = ",
      scales::pvalue(
        best_fdr,
        accuracy = 0.001
      )
    ),
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
      element_text(size = 9),

    legend.position =
      "bottom"
  )

figure_1e_candidate <- (
  correlation_summary_plot |
    best_pathway_scatter
) +
  patchwork::plot_layout(
    widths = c(0.95, 1.05)
  )

# ------------------------------------------------------------
# 12. Supplementary state-specific heatmap
# ------------------------------------------------------------

state_correlation_heatmap <- ggplot(
  state_correlations,
  aes(
    x = pathway,
    y = cell_state_short,
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
        is.na(spearman_rho) ~
          paste0(
            "NA\nn+ = ",
            nonzero_MMP13_donors
          ),

        TRUE ~
          paste0(
            sprintf("%.2f", spearman_rho),
            significance_label
          )
      )
    ),
    size = 2.8
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey90",
    name = "Spearman ρ"
  ) +
  labs(
    title =
      "State-specific donor-level pathway–MMP13 relationships",

    subtitle =
      paste(
        "* nominal P < 0.05; ** FDR < 0.05;",
        "NA indicates fewer than three MMP13-positive donors"
      ),

    x =
      "Pathway program",

    y =
      "Annotated chondrocyte state"
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
# 13. Supplementary whole-chondrocyte pathway activity
# ------------------------------------------------------------

donor_long_table <- donor_level_table %>%
  dplyr::select(
    sample_label,
    donor,
    group,
    dplyr::all_of(donor_score_columns)
  ) %>%
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(donor_score_columns),

    names_to =
      "score_column",

    values_to =
      "module_score"
  ) %>%
  dplyr::left_join(
    donor_pathway_key,
    by = "score_column"
  )

donor_group_plot <- ggplot(
  donor_long_table,
  aes(
    x = group,
    y = module_score,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.60,
    outlier.shape = NA,
    alpha = 0.70
  ) +
  geom_point(
    position = position_jitter(
      width = 0.07,
      height = 0
    ),
    size = 2.0,
    alpha = 0.85
  ) +
  facet_wrap(
    ~ pathway,
    scales = "free_y",
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#D55E00",
      "OA" = "#009E9E"
    )
  ) +
  labs(
    title =
      "Whole-chondrocyte donor-level pathway activity",

    subtitle =
      "Each point represents one independent donor",

    x =
      NULL,

    y =
      "Cell-number-weighted mean module score",

    fill =
      "Group"
  ) +
  theme_bw(
    base_size = 10.5
  ) +
  theme(
    panel.grid.minor =
      element_blank(),

    legend.position =
      "top"
  )

# ------------------------------------------------------------
# 14. Save figures
# ------------------------------------------------------------

main_output_base <- file.path(
  main_figure_folder,
  paste0(
    "Figure1E_CANDIDATE_corrected_",
    "donor_level_pathway_MMP13"
  )
)

ggsave(
  filename = paste0(
    main_output_base,
    ".pdf"
  ),
  plot = figure_1e_candidate,
  width = 11.5,
  height = 5.2,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = paste0(
    main_output_base,
    ".tiff"
  ),
  plot = figure_1e_candidate,
  width = 11.5,
  height = 5.2,
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
  plot = figure_1e_candidate,
  width = 11.5,
  height = 5.2,
  units = "in",
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary_figure_folder,
    "FigureS18A_state_specific_correlation_heatmap.pdf"
  ),
  plot = state_correlation_heatmap,
  width = 9.0,
  height = 6.5,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary_figure_folder,
    "FigureS18A_state_specific_correlation_heatmap.tiff"
  ),
  plot = state_correlation_heatmap,
  width = 9.0,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary_figure_folder,
    "FigureS18B_whole_chondrocyte_pathway_activity.pdf"
  ),
  plot = donor_group_plot,
  width = 9.5,
  height = 7.0,
  units = "in",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    supplementary_figure_folder,
    "FigureS18B_whole_chondrocyte_pathway_activity.tiff"
  ),
  plot = donor_group_plot,
  width = 9.5,
  height = 7.0,
  units = "in",
  dpi = 600,
  compression = "lzw",
  limitsize = FALSE
)

# ------------------------------------------------------------
# 15. Save tables
# ------------------------------------------------------------

utils::write.csv(
  donor_level_table,
  file = file.path(
    table_folder,
    "GSE220243_corrected_whole_chondrocyte_donor_table.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  global_correlations,
  file = file.path(
    table_folder,
    "GSE220243_corrected_donor_level_pathway_MMP13_correlations.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  state_correlations,
  file = file.path(
    table_folder,
    "GSE220243_corrected_state_specific_pathway_MMP13_correlations.csv"
  ),
  row.names = FALSE
)

summary_workbook <- file.path(
  table_folder,
  "GSE220243_corrected_donor_level_Figure1E_summary.xlsx"
)

openxlsx::write.xlsx(
  list(
    Donor_level =
      as.data.frame(donor_level_table),

    Global_correlations =
      as.data.frame(global_correlations),

    State_correlations =
      as.data.frame(state_correlations)
  ),
  file = summary_workbook,
  overwrite = TRUE
)

caption_text <- paste(
  "Figure 1E. Donor-level relationships between cartilage pathway",
  "activity and MMP13 expression. Left, Spearman correlations between",
  "cell-number-weighted pathway module scores and pseudobulk MMP13",
  "expression across 12 independent cartilage donors; error bars",
  "represent bootstrap 95% confidence intervals. Right, donor-level",
  "relationship for the pathway with the strongest statistical",
  "support. Multiple cell states from the same donor were collapsed",
  "before inference to avoid treating repeated donor-state observations",
  "as independent biological replicates."
)

writeLines(
  caption_text,
  con = file.path(
    table_folder,
    "GSE220243_Figure1E_corrected_caption_draft.txt"
  )
)

# ------------------------------------------------------------
# 16. Final summary
# ------------------------------------------------------------

message("")
message("================================================")
message("SCRIPT 19B SELESAI")
message("================================================")

message(
  "Independent biological units : ",
  nrow(donor_level_table),
  " donors"
)

message(
  "Repeated donor-state rows     : COLLAPSED"
)

message(
  "Best-supported pathway        : ",
  best_pathway_name
)

message(
  "Spearman rho                  : ",
  sprintf("%.3f", best_rho)
)

message(
  "P value                       : ",
  signif(best_p, 4)
)

message(
  "FDR                           : ",
  signif(best_fdr, 4)
)

message(
  "Figure 1E candidate TIFF      : ",
  paste0(
    main_output_base,
    ".tiff"
  )
)

message(
  "Supplementary figures         : ",
  supplementary_figure_folder
)

message(
  "Figure 1E status              : CANDIDATE; inspect before manuscript use"
)

message("================================================")
