# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 12_GSE169077_RMA_QC.R
# PURPOSE : Quality control of RMA-processed GSE169077
#
# DATASET:
#   5 normal RNA pools
#   6 late-stage OA RNA pools
#
# IMPORTANT:
#   Statistical unit = RNA pool
#   The 5 individuals contributing to each pool are not
#   treated as separate biological replicates.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Memuat paket
# ------------------------------------------------------------

required_packages <- c(
  "limma",
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "openxlsx"
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
  stop(
    "Paket berikut belum tersedia: ",
    paste(
      missing_packages,
      collapse = ", "
    )
  )
}

library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(openxlsx)

message(
  "Semua paket QC GSE169077 berhasil dimuat."
)

# ------------------------------------------------------------
# 2. Memastikan folder tersedia
# ------------------------------------------------------------

required_folders <- c(
  "data_processed",
  "results",
  "results/figures",
  "results/tables"
)

for (current_folder in required_folders) {
  dir.create(
    current_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------
# 3. Memuat data preparation
# ------------------------------------------------------------

expression_rma <- readRDS(
  "data_processed/GSE169077_expression_RMA_11_pools.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE169077_metadata_11_pools.rds"
)

design_group <- readRDS(
  "data_processed/GSE169077_design_group.rds"
)

feature_annotation <- readRDS(
  "data_processed/GSE169077_feature_annotation_original.rds"
)

mmp13_annotation <- readRDS(
  "data_processed/GSE169077_MMP13_probe_annotation.rds"
)

cat("\nDimensi RMA expression matrix:\n")
print(
  dim(
    expression_rma
  )
)

cat("\nDimensi metadata:\n")
print(
  dim(
    sample_metadata
  )
)

cat("\nDimensi design matrix:\n")
print(
  dim(
    design_group
  )
)

cat("\nJumlah anotasi MMP13:\n")
print(
  nrow(
    mmp13_annotation
  )
)

# ------------------------------------------------------------
# 4. Validasi data
# ------------------------------------------------------------

stopifnot(
  is.matrix(
    expression_rma
  )
)

stopifnot(
  nrow(expression_rma) == 22283
)

stopifnot(
  ncol(expression_rma) == 11
)

stopifnot(
  nrow(sample_metadata) == 11
)

stopifnot(
  identical(
    colnames(expression_rma),
    sample_metadata$geo_accession
  )
)

stopifnot(
  identical(
    rownames(design_group),
    sample_metadata$geo_accession
  )
)

stopifnot(
  nrow(design_group) ==
    ncol(expression_rma)
)

stopifnot(
  qr(design_group)$rank ==
    ncol(design_group)
)

stopifnot(
  "groupOA" %in%
    colnames(design_group)
)

stopifnot(
  !anyNA(expression_rma)
)

stopifnot(
  all(
    is.finite(
      expression_rma
    )
  )
)

stopifnot(
  nrow(mmp13_annotation) >= 1
)

message(
  "Expression matrix, metadata, design, dan anotasi berhasil divalidasi."
)

# ------------------------------------------------------------
# 5. Menyiapkan metadata QC
# ------------------------------------------------------------

sample_metadata_qc <- sample_metadata %>%
  dplyr::mutate(
    
    group = factor(
      as.character(group),
      levels = c(
        "Normal",
        "OA"
      )
    ),
    
    short_label = dplyr::coalesce(
      as.character(pool_label),
      as.character(title_original),
      as.character(geo_accession)
    ),
    
    display_group = factor(
      dplyr::if_else(
        group == "Normal",
        "Normal cartilage",
        "Late-stage OA cartilage"
      ),
      levels = c(
        "Normal cartilage",
        "Late-stage OA cartilage"
      )
    )
  )

cat("\nDistribusi group:\n")

print(
  table(
    sample_metadata_qc$group,
    useNA = "ifany"
  )
)

stopifnot(
  sum(
    sample_metadata_qc$group ==
      "Normal"
  ) == 5
)

stopifnot(
  sum(
    sample_metadata_qc$group ==
      "OA"
  ) == 6
)

stopifnot(
  !anyNA(
    sample_metadata_qc$group
  )
)

cat("\nMetadata QC:\n")

print(
  sample_metadata_qc %>%
    dplyr::select(
      geo_accession,
      short_label,
      group,
      pooled_sample,
      individuals_per_pool,
      statistical_unit
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 6. Menyelaraskan feature annotation
# ------------------------------------------------------------

annotation_index <- match(
  rownames(expression_rma),
  feature_annotation$probe_id
)

if (anyNA(annotation_index)) {
  stop(
    "Ada probeset RMA yang tidak ditemukan ",
    "dalam feature annotation."
  )
}

annotation_aligned <- feature_annotation[
  annotation_index,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(expression_rma),
    as.character(
      annotation_aligned$probe_id
    )
  )
)

gene_symbol <- trimws(
  as.character(
    annotation_aligned[[
      "Gene Symbol"
    ]]
  )
)

keep_annotated <- (
  !is.na(gene_symbol) &
    gene_symbol != ""
)

expression_biological <- expression_rma[
  keep_annotated,
  ,
  drop = FALSE
]

annotation_biological <- annotation_aligned[
  keep_annotated,
  ,
  drop = FALSE
]

cat("\nJumlah seluruh probeset:\n")
print(
  nrow(
    expression_rma
  )
)

cat("\nJumlah probeset dengan gene symbol:\n")
print(
  nrow(
    expression_biological
  )
)

# ------------------------------------------------------------
# 7. Ringkasan distribusi setiap pool
# ------------------------------------------------------------

sample_expression_summary <- tibble::tibble(
  
  geo_accession =
    colnames(expression_rma),
  
  minimum =
    apply(
      expression_rma,
      2,
      min
    ),
  
  q1 =
    apply(
      expression_rma,
      2,
      stats::quantile,
      probs = 0.25
    ),
  
  median =
    apply(
      expression_rma,
      2,
      stats::median
    ),
  
  mean =
    colMeans(
      expression_rma
    ),
  
  q3 =
    apply(
      expression_rma,
      2,
      stats::quantile,
      probs = 0.75
    ),
  
  maximum =
    apply(
      expression_rma,
      2,
      max
    )
) %>%
  dplyr::left_join(
    sample_metadata_qc %>%
      dplyr::select(
        geo_accession,
        short_label,
        group,
        display_group
      ),
    by = "geo_accession"
  )

cat("\nRingkasan distribusi expression:\n")

print(
  sample_expression_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 8. Boxplot distribusi RMA
# ------------------------------------------------------------

sample_plot_colors <- ifelse(
  sample_metadata_qc$group == "Normal",
  "#F8766D",
  "#00BFC4"
)

pdf(
  file =
    "results/figures/GSE169077_RMA_expression_boxplot.pdf",
  width = 8,
  height = 5.5
)

boxplot(
  expression_rma,
  names =
    sample_metadata_qc$short_label,
  col =
    sample_plot_colors,
  outline = FALSE,
  las = 2,
  cex.axis = 0.9,
  main =
    "RMA expression distributions in GSE169077",
  xlab =
    "RNA pool",
  ylab =
    expression(
      RMA~log[2]~expression
    )
)

legend(
  "topright",
  legend = c(
    "Normal",
    "OA"
  ),
  fill = c(
    "#F8766D",
    "#00BFC4"
  ),
  border = NA,
  bty = "n"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE169077_RMA_expression_boxplot.tiff",
  width = 8,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

boxplot(
  expression_rma,
  names =
    sample_metadata_qc$short_label,
  col =
    sample_plot_colors,
  outline = FALSE,
  las = 2,
  cex.axis = 0.9,
  main =
    "RMA expression distributions in GSE169077",
  xlab =
    "RNA pool",
  ylab =
    expression(
      RMA~log[2]~expression
    )
)

legend(
  "topright",
  legend = c(
    "Normal",
    "OA"
  ),
  fill = c(
    "#F8766D",
    "#00BFC4"
  ),
  border = NA,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 9. Relative log expression
# ------------------------------------------------------------

probe_median_expression <- apply(
  expression_rma,
  MARGIN = 1,
  FUN = stats::median
)

relative_log_expression <- sweep(
  expression_rma,
  MARGIN = 1,
  STATS = probe_median_expression,
  FUN = "-"
)

rle_summary <- tibble::tibble(
  
  geo_accession =
    colnames(relative_log_expression),
  
  RLE_median =
    apply(
      relative_log_expression,
      2,
      stats::median
    ),
  
  RLE_q1 =
    apply(
      relative_log_expression,
      2,
      stats::quantile,
      probs = 0.25
    ),
  
  RLE_q3 =
    apply(
      relative_log_expression,
      2,
      stats::quantile,
      probs = 0.75
    ),
  
  RLE_IQR =
    apply(
      relative_log_expression,
      2,
      stats::IQR
    )
) %>%
  dplyr::left_join(
    sample_metadata_qc %>%
      dplyr::select(
        geo_accession,
        short_label,
        group,
        display_group
      ),
    by = "geo_accession"
  )

cat("\nRelative log-expression summary:\n")

print(
  rle_summary,
  n = Inf,
  width = Inf
)

pdf(
  file =
    "results/figures/GSE169077_RLE_boxplot.pdf",
  width = 8,
  height = 5.5
)

boxplot(
  relative_log_expression,
  names =
    sample_metadata_qc$short_label,
  col =
    sample_plot_colors,
  outline = FALSE,
  las = 2,
  cex.axis = 0.9,
  main =
    "Relative log expression in GSE169077",
  xlab =
    "RNA pool",
  ylab =
    "Relative log expression"
)

abline(
  h = 0,
  lty = 2,
  lwd = 1
)

legend(
  "topright",
  legend = c(
    "Normal",
    "OA"
  ),
  fill = c(
    "#F8766D",
    "#00BFC4"
  ),
  border = NA,
  bty = "n"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE169077_RLE_boxplot.tiff",
  width = 8,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

boxplot(
  relative_log_expression,
  names =
    sample_metadata_qc$short_label,
  col =
    sample_plot_colors,
  outline = FALSE,
  las = 2,
  cex.axis = 0.9,
  main =
    "Relative log expression in GSE169077",
  xlab =
    "RNA pool",
  ylab =
    "Relative log expression"
)

abline(
  h = 0,
  lty = 2,
  lwd = 1
)

legend(
  "topright",
  legend = c(
    "Normal",
    "OA"
  ),
  fill = c(
    "#F8766D",
    "#00BFC4"
  ),
  border = NA,
  bty = "n"
)

dev.off()

# ------------------------------------------------------------
# 10. Memilih probeset variabel untuk PCA
# ------------------------------------------------------------

probe_variance <- apply(
  expression_biological,
  MARGIN = 1,
  FUN = stats::var
)

probe_variance[
  !is.finite(
    probe_variance
  )
] <- 0

number_top_variable <- min(
  5000,
  length(
    probe_variance
  )
)

top_variable_index <- order(
  probe_variance,
  decreasing = TRUE
)[
  seq_len(
    number_top_variable
  )
]

expression_top_variable <- expression_biological[
  top_variable_index,
  ,
  drop = FALSE
]

cat("\nJumlah probeset untuk PCA:\n")
print(
  nrow(
    expression_top_variable
  )
)

# ------------------------------------------------------------
# 11. PCA
# ------------------------------------------------------------

pca_result <- prcomp(
  t(
    expression_top_variable
  ),
  center = TRUE,
  scale. = FALSE
)

pca_variance_explained <- (
  pca_result$sdev^2 /
    sum(
      pca_result$sdev^2
    )
) * 100

pca_table <- sample_metadata_qc %>%
  dplyr::mutate(
    
    PC1 =
      pca_result$x[, 1],
    
    PC2 =
      pca_result$x[, 2]
  )

cat("\nVariance explained PCA:\n")

print(
  round(
    pca_variance_explained[
      seq_len(
        min(
          5,
          length(
            pca_variance_explained
          )
        )
      )
    ],
    digits = 2
  )
)

pca_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = display_group
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  geom_text(
    aes(
      label = short_label
    ),
    vjust = -0.8,
    size = 3.4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "PCA of RMA-processed GSE169077",
    
    subtitle =
      "Top 5,000 variable annotated probesets",
    
    x = paste0(
      "PC1 (",
      round(
        pca_variance_explained[1],
        1
      ),
      "%)"
    ),
    
    y = paste0(
      "PC2 (",
      round(
        pca_variance_explained[2],
        1
      ),
      "%)"
    ),
    
    color =
      "Cartilage group"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(
  pca_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_PCA_RMA.pdf",
  plot =
    pca_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE169077_PCA_RMA.tiff",
  plot =
    pca_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 12. MDS
# ------------------------------------------------------------

mds_result <- limma::plotMDS(
  expression_biological,
  top = min(
    500,
    nrow(
      expression_biological
    )
  ),
  gene.selection = "pairwise",
  plot = FALSE
)

mds_table <- sample_metadata_qc %>%
  dplyr::mutate(
    
    MDS_dimension_1 =
      as.numeric(
        mds_result$x
      ),
    
    MDS_dimension_2 =
      as.numeric(
        mds_result$y
      )
  )

mds_plot <- ggplot(
  mds_table,
  aes(
    x = MDS_dimension_1,
    y = MDS_dimension_2,
    color = display_group
  )
) +
  geom_point(
    size = 4,
    alpha = 0.9
  ) +
  geom_text(
    aes(
      label = short_label
    ),
    vjust = -0.8,
    size = 3.4,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "MDS plot of GSE169077",
    
    subtitle =
      "Top 500 pairwise variable annotated probesets",
    
    x =
      "Leading log-fold-change dimension 1",
    
    y =
      "Leading log-fold-change dimension 2",
    
    color =
      "Cartilage group"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(
  mds_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_MDS_RMA.pdf",
  plot =
    mds_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE169077_MDS_RMA.tiff",
  plot =
    mds_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 13. Sample-to-sample correlation
# ------------------------------------------------------------

correlation_matrix <- stats::cor(
  expression_biological,
  method = "pearson"
)

rownames(
  correlation_matrix
) <- sample_metadata_qc$short_label

colnames(
  correlation_matrix
) <- sample_metadata_qc$short_label

correlation_distance <- 1 -
  correlation_matrix

diag(
  correlation_distance
) <- 0

sample_clustering <- stats::hclust(
  stats::as.dist(
    correlation_distance
  ),
  method = "average"
)

clustered_labels <- rownames(
  correlation_matrix
)[
  sample_clustering$order
]

correlation_long <- as.data.frame(
  as.table(
    correlation_matrix
  ),
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble()

colnames(
  correlation_long
) <- c(
  "Sample_1",
  "Sample_2",
  "Correlation"
)

correlation_long <- correlation_long %>%
  dplyr::mutate(
    
    Sample_1 = factor(
      Sample_1,
      levels =
        clustered_labels
    ),
    
    Sample_2 = factor(
      Sample_2,
      levels =
        rev(
          clustered_labels
        )
    )
  )

correlation_plot <- ggplot(
  correlation_long,
  aes(
    x = Sample_1,
    y = Sample_2,
    fill = Correlation
  )
) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "steelblue"
  ) +
  labs(
    title =
      "Sample-to-sample correlation in GSE169077",
    
    subtitle =
      "RMA-processed RNA pools",
    
    x = NULL,
    y = NULL,
    
    fill =
      "Pearson\ncorrelation"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1
      ),
    
    panel.grid =
      element_blank()
  )

print(
  correlation_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_sample_correlation.pdf",
  plot =
    correlation_plot,
  width = 7,
  height = 6
)

ggsave(
  filename =
    "results/figures/GSE169077_sample_correlation.tiff",
  plot =
    correlation_plot,
  width = 7,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 14. Hierarchical clustering dendrogram
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE169077_sample_dendrogram.pdf",
  width = 8,
  height = 5.5
)

plot(
  sample_clustering,
  main =
    "Hierarchical clustering of GSE169077 RNA pools",
  xlab =
    "RNA pool",
  sub = "",
  ylab =
    "1 - Pearson correlation",
  cex = 0.9
)

dev.off()

tiff(
  filename =
    "results/figures/GSE169077_sample_dendrogram.tiff",
  width = 8,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

plot(
  sample_clustering,
  main =
    "Hierarchical clustering of GSE169077 RNA pools",
  xlab =
    "RNA pool",
  sub = "",
  ylab =
    "1 - Pearson correlation",
  cex = 0.9
)

dev.off()

# ------------------------------------------------------------
# 15. Sample connectivity
# ------------------------------------------------------------

correlation_without_diagonal <-
  correlation_matrix

diag(
  correlation_without_diagonal
) <- NA_real_

mean_correlation_to_others <- rowMeans(
  correlation_without_diagonal,
  na.rm = TRUE
)

connectivity_z_score <- as.numeric(
  scale(
    mean_correlation_to_others
  )
)

names(
  connectivity_z_score
) <- names(
  mean_correlation_to_others
)

sample_connectivity <- sample_metadata_qc %>%
  dplyr::mutate(
    
    mean_correlation_to_others =
      mean_correlation_to_others[
        short_label
      ],
    
    connectivity_z_score =
      connectivity_z_score[
        short_label
      ],
    
    candidate_sample_outlier =
      connectivity_z_score < -3
  ) %>%
  dplyr::arrange(
    connectivity_z_score
  )

cat("\nSample connectivity:\n")

print(
  sample_connectivity %>%
    dplyr::select(
      short_label,
      geo_accession,
      group,
      mean_correlation_to_others,
      connectivity_z_score,
      candidate_sample_outlier
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

connectivity_plot <- ggplot(
  sample_connectivity,
  aes(
    x = reorder(
      short_label,
      connectivity_z_score
    ),
    y = connectivity_z_score,
    color = display_group
  )
) +
  geom_hline(
    yintercept = -3,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_point(
    size = 3.5
  ) +
  coord_flip() +
  labs(
    title =
      "Sample connectivity in GSE169077",
    
    subtitle =
      "Dashed line indicates the exploratory z-score threshold of -3",
    
    x =
      "RNA pool",
    
    y =
      "Connectivity z-score",
    
    color =
      "Cartilage group"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(
  connectivity_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_sample_connectivity.pdf",
  plot =
    connectivity_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE169077_sample_connectivity.tiff",
  plot =
    connectivity_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 16. MMP13 expression
# ------------------------------------------------------------

mmp13_probe_ids <- intersect(
  unique(
    as.character(
      mmp13_annotation$probe_id
    )
  ),
  rownames(
    expression_rma
  )
)

cat("\nProbe MMP13 tersedia:\n")
print(
  mmp13_probe_ids
)

if (length(mmp13_probe_ids) == 0) {
  stop(
    "Probe MMP13 tidak ditemukan dalam expression matrix."
  )
}

mmp13_expression_long <- as.data.frame(
  t(
    expression_rma[
      mmp13_probe_ids,
      ,
      drop = FALSE
    ]
  ),
  check.names = FALSE
) %>%
  tibble::rownames_to_column(
    "geo_accession"
  ) %>%
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        mmp13_probe_ids
      ),
    names_to =
      "probe_id",
    values_to =
      "MMP13_expression"
  ) %>%
  dplyr::left_join(
    sample_metadata_qc %>%
      dplyr::select(
        geo_accession,
        short_label,
        group,
        display_group
      ),
    by = "geo_accession"
  ) %>%
  tibble::as_tibble()

cat("\nEkspresi MMP13:\n")

print(
  mmp13_expression_long,
  n = Inf,
  width = Inf
)

mmp13_group_summary <- mmp13_expression_long %>%
  dplyr::group_by(
    probe_id,
    group
  ) %>%
  dplyr::summarise(
    
    number_of_pools =
      dplyr::n(),
    
    mean_expression =
      mean(
        MMP13_expression
      ),
    
    sd_expression =
      stats::sd(
        MMP13_expression
      ),
    
    median_expression =
      stats::median(
        MMP13_expression
      ),
    
    .groups = "drop"
  )

mmp13_descriptive_effect <- mmp13_group_summary %>%
  dplyr::select(
    probe_id,
    group,
    mean_expression
  ) %>%
  tidyr::pivot_wider(
    names_from =
      group,
    values_from =
      mean_expression
  ) %>%
  dplyr::mutate(
    
    mean_difference_OA_minus_Normal =
      OA - Normal,
    
    approximate_fold_change =
      2^(
        mean_difference_OA_minus_Normal
      )
  )

cat("\nRingkasan MMP13 berdasarkan group:\n")

print(
  mmp13_group_summary,
  n = Inf,
  width = Inf
)

cat("\nEfek deskriptif MMP13:\n")

print(
  mmp13_descriptive_effect,
  n = Inf,
  width = Inf
)

mmp13_plot <- ggplot(
  mmp13_expression_long,
  aes(
    x = display_group,
    y = MMP13_expression,
    color = display_group
  )
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.20,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.08,
    size = 3.5,
    alpha = 0.9
  ) +
  geom_text(
    aes(
      label = short_label
    ),
    vjust = -0.8,
    size = 3.2,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ probe_id,
    scales = "free_y"
  ) +
  labs(
    title =
      "MMP13 expression in GSE169077",
    
    subtitle =
      "RMA-processed pooled cartilage samples",
    
    x = NULL,
    
    y =
      expression(
        RMA~log[2]~expression
      ),
    
    color =
      "Cartilage group"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(
  mmp13_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_MMP13_RMA_expression.pdf",
  plot =
    mmp13_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE169077_MMP13_RMA_expression.tiff",
  plot =
    mmp13_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 17. Design diagnostics
# ------------------------------------------------------------

design_diagnostics <- tibble::tibble(
  
  model =
    "OA versus Normal RNA pools",
  
  formula =
    "~ group",
  
  number_of_statistical_units =
    nrow(
      design_group
    ),
  
  number_of_parameters =
    ncol(
      design_group
    ),
  
  rank =
    qr(
      design_group
    )$rank,
  
  full_rank =
    qr(
      design_group
    )$rank ==
    ncol(
      design_group
    ),
  
  residual_degrees_of_freedom =
    nrow(
      design_group
    ) -
    qr(
      design_group
    )$rank,
  
  condition_number =
    kappa(
      design_group
    )
)

cat("\nDesign diagnostics:\n")

print(
  design_diagnostics,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 18. QC summary
# ------------------------------------------------------------

qc_summary <- tibble::tibble(
  
  metric = c(
    "Number of RNA pools",
    "Normal pools",
    "OA pools",
    "Individuals contributing per pool",
    "Statistical unit",
    "Total RMA probesets",
    "Annotated biological probesets",
    "Top variable probesets used in PCA",
    "PC1 variance explained",
    "PC2 variance explained",
    "Minimum RLE median",
    "Maximum RLE median",
    "Candidate sample outliers",
    "Design residual degrees of freedom",
    "MMP13 probes available",
    "MMP13 descriptive OA-Normal difference"
  ),
  
  value = c(
    ncol(
      expression_rma
    ),
    
    sum(
      sample_metadata_qc$group ==
        "Normal"
    ),
    
    sum(
      sample_metadata_qc$group ==
        "OA"
    ),
    
    5,
    
    "RNA pool",
    
    nrow(
      expression_rma
    ),
    
    nrow(
      expression_biological
    ),
    
    nrow(
      expression_top_variable
    ),
    
    round(
      pca_variance_explained[1],
      4
    ),
    
    round(
      pca_variance_explained[2],
      4
    ),
    
    round(
      min(
        rle_summary$RLE_median
      ),
      6
    ),
    
    round(
      max(
        rle_summary$RLE_median
      ),
      6
    ),
    
    sum(
      sample_connectivity$
        candidate_sample_outlier
    ),
    
    design_diagnostics$
      residual_degrees_of_freedom,
    
    length(
      mmp13_probe_ids
    ),
    
    round(
      mmp13_descriptive_effect$
        mean_difference_OA_minus_Normal[1],
      6
    )
  )
)

cat("\nQC summary:\n")

print(
  qc_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 19. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  sample_metadata_qc,
  file =
    "data_processed/GSE169077_metadata_RMA_QC.rds"
)

saveRDS(
  expression_biological,
  file =
    "data_processed/GSE169077_expression_RMA_annotated_biological.rds"
)

saveRDS(
  annotation_biological,
  file =
    "data_processed/GSE169077_annotation_biological.rds"
)

saveRDS(
  relative_log_expression,
  file =
    "data_processed/GSE169077_relative_log_expression.rds"
)

saveRDS(
  pca_result,
  file =
    "data_processed/GSE169077_PCA_result.rds"
)

saveRDS(
  mds_result,
  file =
    "data_processed/GSE169077_MDS_result.rds"
)

saveRDS(
  correlation_matrix,
  file =
    "data_processed/GSE169077_sample_correlation_matrix.rds"
)

saveRDS(
  sample_connectivity,
  file =
    "data_processed/GSE169077_sample_connectivity.rds"
)

saveRDS(
  mmp13_expression_long,
  file =
    "data_processed/GSE169077_MMP13_expression_long.rds"
)

saveRDS(
  mmp13_descriptive_effect,
  file =
    "data_processed/GSE169077_MMP13_descriptive_effect.rds"
)

# ------------------------------------------------------------
# 20. Menyimpan Excel workbook
# ------------------------------------------------------------

correlation_matrix_export <- correlation_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "sample"
  )

openxlsx::write.xlsx(
  list(
    
    QC_Summary =
      as.data.frame(
        qc_summary
      ),
    
    Sample_Metadata =
      as.data.frame(
        sample_metadata_qc
      ),
    
    Expression_Summary =
      as.data.frame(
        sample_expression_summary
      ),
    
    RLE_Summary =
      as.data.frame(
        rle_summary
      ),
    
    PCA_Coordinates =
      as.data.frame(
        pca_table
      ),
    
    MDS_Coordinates =
      as.data.frame(
        mds_table
      ),
    
    Sample_Connectivity =
      as.data.frame(
        sample_connectivity
      ),
    
    Sample_Correlation =
      as.data.frame(
        correlation_matrix_export
      ),
    
    MMP13_Expression =
      as.data.frame(
        mmp13_expression_long
      ),
    
    MMP13_Group_Summary =
      as.data.frame(
        mmp13_group_summary
      ),
    
    MMP13_Descriptive_Effect =
      as.data.frame(
        mmp13_descriptive_effect
      ),
    
    Design_Diagnostics =
      as.data.frame(
        design_diagnostics
      )
  ),
  
  file =
    "results/tables/GSE169077_RMA_QC_results.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 21. Session information
# ------------------------------------------------------------

sink(
  "results/tables/GSE169077_RMA_QC_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 22. Pesan akhir
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE169077 RMA QC SELESAI")
message("Jumlah RNA pools            : ",
        ncol(expression_rma))
message("Normal pools                : ",
        sum(sample_metadata_qc$group == "Normal"))
message("OA pools                    : ",
        sum(sample_metadata_qc$group == "OA"))
message("Annotated probesets         : ",
        nrow(expression_biological))
message("PC1 variance                : ",
        round(pca_variance_explained[1], 2),
        "%")
message("PC2 variance                : ",
        round(pca_variance_explained[2], 2),
        "%")
message("Candidate sample outliers   : ",
        sum(
          sample_connectivity$
            candidate_sample_outlier
        ))
message("MMP13 probes available      : ",
        length(mmp13_probe_ids))
message("MMP13 OA-Normal difference  : ",
        round(
          mmp13_descriptive_effect$
            mean_difference_OA_minus_Normal[1],
          4
        ))
message("Approximate MMP13 fold change: ",
        round(
          mmp13_descriptive_effect$
            approximate_fold_change[1],
          3
        ))
message("Design residual df          : ",
        design_diagnostics$
          residual_degrees_of_freedom)
message("============================================")