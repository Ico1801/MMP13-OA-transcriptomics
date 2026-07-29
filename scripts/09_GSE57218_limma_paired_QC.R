# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 09_GSE57218_limma_paired_QC.R
# PURPOSE : Quality control untuk paired microarray analysis
#           GSE57218
# DESIGN  : OA-affected vs preserved cartilage dari pasien
#           yang sama
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
    paste(missing_packages, collapse = ", ")
  )
}

library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(openxlsx)

message("Semua paket QC GSE57218 berhasil dimuat.")

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
# 3. Memuat data hasil preparation
# ------------------------------------------------------------

expression_primary <- readRDS(
  "data_processed/GSE57218_expression_primary_33_pairs.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE57218_metadata_primary_33_pairs.rds"
)

feature_annotation <- readRDS(
  "data_processed/GSE57218_feature_annotation_original.rds"
)

mmp13_annotation <- readRDS(
  "data_processed/GSE57218_MMP13_probe_annotation.rds"
)

pair_summary <- readRDS(
  "data_processed/GSE57218_pair_summary.rds"
)

cat("\nDimensi expression matrix:\n")
print(dim(expression_primary))

cat("\nDimensi metadata:\n")
print(dim(sample_metadata))

cat("\nDimensi feature annotation:\n")
print(dim(feature_annotation))

cat("\nJumlah pasangan:\n")
print(nrow(pair_summary))

# ------------------------------------------------------------
# 4. Validasi expression matrix dan metadata
# ------------------------------------------------------------

stopifnot(
  nrow(expression_primary) == 48777
)

stopifnot(
  ncol(expression_primary) == 66
)

stopifnot(
  nrow(sample_metadata) == 66
)

stopifnot(
  ncol(expression_primary) ==
    nrow(sample_metadata)
)

stopifnot(
  identical(
    colnames(expression_primary),
    sample_metadata$expression_column
  )
)

stopifnot(
  !anyNA(expression_primary)
)

stopifnot(
  all(is.finite(expression_primary))
)

stopifnot(
  nrow(pair_summary) == 33
)

stopifnot(
  all(pair_summary$complete_pair)
)

message(
  "Expression matrix dan metadata GSE57218 berhasil divalidasi."
)

# ------------------------------------------------------------
# 5. Menyiapkan metadata QC
# ------------------------------------------------------------

sample_metadata_qc <- sample_metadata %>%
  dplyr::mutate(
    
    pair_id = factor(
      pair_id,
      levels = unique(pair_id)
    ),
    
    condition = factor(
      condition,
      levels = c(
        "Preserved",
        "OA"
      )
    ),
    
    age_years = as.numeric(
      `age (yrs):ch1`
    ),
    
    sex = factor(
      as.character(`Sex:ch1`)
    ),
    
    display_condition = factor(
      dplyr::if_else(
        condition == "Preserved",
        "Preserved cartilage",
        "OA-affected cartilage"
      ),
      levels = c(
        "Preserved cartilage",
        "OA-affected cartilage"
      )
    ),
    
    short_label = paste0(
      sprintf(
        "P%02d",
        match(
          pair_id,
          levels(pair_id)
        )
      ),
      dplyr::if_else(
        condition == "Preserved",
        "_P",
        "_OA"
      )
    )
  )

stopifnot(
  !anyNA(sample_metadata_qc$age_years)
)

stopifnot(
  !anyNA(sample_metadata_qc$sex)
)

cat("\nMetadata QC:\n")

print(
  sample_metadata_qc %>%
    dplyr::select(
      expression_column,
      pair_id,
      pair_number,
      condition,
      age_years,
      sex,
      short_label
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 6. Memastikan umur dan jenis kelamin sama dalam setiap pair
# ------------------------------------------------------------

pair_covariate_check <- sample_metadata_qc %>%
  dplyr::group_by(
    pair_id
  ) %>%
  dplyr::summarise(
    
    number_of_samples =
      dplyr::n(),
    
    unique_age_values =
      dplyr::n_distinct(age_years),
    
    unique_sex_values =
      dplyr::n_distinct(sex),
    
    pair_covariates_consistent =
      number_of_samples == 2 &
      unique_age_values == 1 &
      unique_sex_values == 1,
    
    .groups = "drop"
  )

cat("\nKonsistensi covariate dalam pasangan:\n")

print(
  pair_covariate_check,
  n = Inf,
  width = Inf
)

if (
  any(
    !pair_covariate_check$pair_covariates_consistent
  )
) {
  stop(
    "Ada umur atau jenis kelamin yang tidak konsisten ",
    "di dalam pasangan."
  )
}

# ------------------------------------------------------------
# 7. Menyelaraskan anotasi dengan expression matrix
# ------------------------------------------------------------

annotation_index <- match(
  rownames(expression_primary),
  feature_annotation$probe_id
)

if (anyNA(annotation_index)) {
  stop(
    "Ada probe expression matrix yang tidak ditemukan ",
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
    rownames(expression_primary),
    as.character(
      annotation_aligned$probe_id
    )
  )
)

gene_symbol <- trimws(
  as.character(
    annotation_aligned$Symbol
  )
)

keep_annotated <- (
  !is.na(gene_symbol) &
    gene_symbol != ""
)

expression_biological <- expression_primary[
  keep_annotated,
  ,
  drop = FALSE
]

annotation_biological <- annotation_aligned[
  keep_annotated,
  ,
  drop = FALSE
]

cat("\nJumlah seluruh probe:\n")
print(nrow(expression_primary))

cat("\nJumlah probe dengan gene symbol:\n")
print(nrow(expression_biological))

# ------------------------------------------------------------
# 8. Pemeriksaan distribusi antar-array
# ------------------------------------------------------------

sample_quantiles <- apply(
  expression_primary,
  MARGIN = 2,
  FUN = stats::quantile,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    1
  ),
  na.rm = TRUE
)

quantile_range_across_samples <- apply(
  sample_quantiles,
  MARGIN = 1,
  FUN = function(x) {
    max(x) - min(x)
  }
)

distribution_summary <- tibble::tibble(
  
  quantile =
    rownames(sample_quantiles),
  
  minimum_across_samples =
    apply(
      sample_quantiles,
      1,
      min
    ),
  
  maximum_across_samples =
    apply(
      sample_quantiles,
      1,
      max
    ),
  
  range_across_samples =
    quantile_range_across_samples
)

cat("\nRentang quantile antarsampel:\n")

print(
  distribution_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 9. Boxplot distribusi expression
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE57218_expression_boxplot.pdf",
  width = 14,
  height = 6
)

boxplot(
  expression_primary,
  names =
    sample_metadata_qc$short_label,
  outline = FALSE,
  las = 2,
  cex.axis = 0.55,
  main =
    "Processed expression distributions in GSE57218",
  xlab =
    "Paired cartilage samples",
  ylab =
    "Processed VST-like expression"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE57218_expression_boxplot.tiff",
  width = 14,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

boxplot(
  expression_primary,
  names =
    sample_metadata_qc$short_label,
  outline = FALSE,
  las = 2,
  cex.axis = 0.55,
  main =
    "Processed expression distributions in GSE57218",
  xlab =
    "Paired cartilage samples",
  ylab =
    "Processed VST-like expression"
)

dev.off()

# ------------------------------------------------------------
# 10. Memilih probe variabel untuk PCA
# ------------------------------------------------------------

probe_variance <- apply(
  expression_biological,
  MARGIN = 1,
  FUN = stats::var
)

probe_variance[
  !is.finite(probe_variance)
] <- 0

number_top_variable <- min(
  5000,
  length(probe_variance)
)

top_variable_index <- order(
  probe_variance,
  decreasing = TRUE
)[
  seq_len(number_top_variable)
]

expression_top_variable <-
  expression_biological[
    top_variable_index,
    ,
    drop = FALSE
  ]

cat("\nJumlah probe untuk PCA:\n")
print(nrow(expression_top_variable))

# ------------------------------------------------------------
# 11. PCA
# ------------------------------------------------------------

pca_result <- prcomp(
  t(expression_top_variable),
  center = TRUE,
  scale. = FALSE
)

pca_variance_explained <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pca_table <- sample_metadata_qc %>%
  dplyr::mutate(
    PC1 =
      pca_result$x[, 1],
    
    PC2 =
      pca_result$x[, 2]
  )

pca_pair_segments <- pca_table %>%
  dplyr::select(
    pair_id,
    condition,
    PC1,
    PC2
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = c(
      PC1,
      PC2
    )
  )

cat("\nVariance explained PCA:\n")

print(
  round(
    pca_variance_explained[1:5],
    digits = 2
  )
)

pca_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = display_condition,
    shape = sex
  )
) +
  geom_segment(
    data = pca_pair_segments,
    aes(
      x = PC1_Preserved,
      y = PC2_Preserved,
      xend = PC1_OA,
      yend = PC2_OA
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    alpha = 0.35
  ) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = short_label),
    vjust = -0.7,
    size = 2.2,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "PCA of paired GSE57218 cartilage samples",
    
    subtitle =
      "Lines connect preserved and OA-affected cartilage from the same patient",
    
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
      "Cartilage condition",
    
    shape =
      "Sex"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(pca_plot)

ggsave(
  filename =
    "results/figures/GSE57218_PCA_paired.pdf",
  plot = pca_plot,
  width = 8,
  height = 6
)

ggsave(
  filename =
    "results/figures/GSE57218_PCA_paired.tiff",
  plot = pca_plot,
  width = 8,
  height = 6,
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
    nrow(expression_biological)
  ),
  gene.selection = "pairwise",
  plot = FALSE
)

mds_table <- sample_metadata_qc %>%
  dplyr::mutate(
    MDS_dimension_1 =
      as.numeric(mds_result$x),
    
    MDS_dimension_2 =
      as.numeric(mds_result$y)
  )

mds_pair_segments <- mds_table %>%
  dplyr::select(
    pair_id,
    condition,
    MDS_dimension_1,
    MDS_dimension_2
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = c(
      MDS_dimension_1,
      MDS_dimension_2
    )
  )

mds_plot <- ggplot(
  mds_table,
  aes(
    x = MDS_dimension_1,
    y = MDS_dimension_2,
    color = display_condition,
    shape = sex
  )
) +
  geom_segment(
    data = mds_pair_segments,
    aes(
      x =
        MDS_dimension_1_Preserved,
      
      y =
        MDS_dimension_2_Preserved,
      
      xend =
        MDS_dimension_1_OA,
      
      yend =
        MDS_dimension_2_OA
    ),
    inherit.aes = FALSE,
    linewidth = 0.35,
    alpha = 0.35
  ) +
  geom_point(
    size = 3,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = short_label),
    vjust = -0.7,
    size = 2.2,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "MDS plot of paired GSE57218 samples",
    
    subtitle =
      "Top 500 pairwise variable annotated probes",
    
    x =
      "Leading log-fold-change dimension 1",
    
    y =
      "Leading log-fold-change dimension 2",
    
    color =
      "Cartilage condition",
    
    shape =
      "Sex"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(mds_plot)

ggsave(
  filename =
    "results/figures/GSE57218_MDS_paired.pdf",
  plot = mds_plot,
  width = 8,
  height = 6
)

ggsave(
  filename =
    "results/figures/GSE57218_MDS_paired.tiff",
  plot = mds_plot,
  width = 8,
  height = 6,
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

rownames(correlation_matrix) <-
  sample_metadata_qc$short_label

colnames(correlation_matrix) <-
  sample_metadata_qc$short_label

correlation_distance <-
  1 - correlation_matrix

diag(correlation_distance) <- 0

sample_clustering <- hclust(
  as.dist(correlation_distance),
  method = "average"
)

clustered_labels <- rownames(
  correlation_matrix
)[
  sample_clustering$order
]

correlation_long <- as.data.frame(
  as.table(correlation_matrix),
  stringsAsFactors = FALSE
) %>%
  tibble::as_tibble()

colnames(correlation_long) <- c(
  "Sample_1",
  "Sample_2",
  "Correlation"
)

correlation_long <- correlation_long %>%
  dplyr::mutate(
    
    Sample_1 = factor(
      Sample_1,
      levels = clustered_labels
    ),
    
    Sample_2 = factor(
      Sample_2,
      levels = rev(
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
  labs(
    title =
      "Sample-to-sample correlation in GSE57218",
    
    subtitle =
      "Primary paired cartilage samples",
    
    x = NULL,
    y = NULL,
    
    fill =
      "Pearson\ncorrelation"
  ) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        size = 5
      ),
    
    axis.text.y =
      element_text(
        size = 5
      ),
    
    panel.grid =
      element_blank()
  )

print(correlation_plot)

ggsave(
  filename =
    "results/figures/GSE57218_sample_correlation.pdf",
  plot = correlation_plot,
  width = 10,
  height = 9
)

ggsave(
  filename =
    "results/figures/GSE57218_sample_correlation.tiff",
  plot = correlation_plot,
  width = 10,
  height = 9,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 14. Hierarchical clustering dendrogram
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE57218_sample_dendrogram.pdf",
  width = 14,
  height = 6
)

plot(
  sample_clustering,
  main =
    "Hierarchical clustering of GSE57218 samples",
  xlab =
    "Sample",
  sub = "",
  ylab =
    "1 - Pearson correlation",
  cex = 0.55
)

dev.off()

tiff(
  filename =
    "results/figures/GSE57218_sample_dendrogram.tiff",
  width = 14,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

plot(
  sample_clustering,
  main =
    "Hierarchical clustering of GSE57218 samples",
  xlab =
    "Sample",
  sub = "",
  ylab =
    "1 - Pearson correlation",
  cex = 0.55
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

names(connectivity_z_score) <-
  names(
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
      pair_id,
      condition,
      mean_correlation_to_others,
      connectivity_z_score,
      candidate_sample_outlier
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 16. Pair-level diagnostics
# ------------------------------------------------------------

pair_ids <- levels(
  sample_metadata_qc$pair_id
)

pair_diagnostics <- lapply(
  pair_ids,
  function(current_pair) {
    
    preserved_column <-
      sample_metadata_qc$expression_column[
        sample_metadata_qc$pair_id ==
          current_pair &
          sample_metadata_qc$condition ==
          "Preserved"
      ]
    
    oa_column <-
      sample_metadata_qc$expression_column[
        sample_metadata_qc$pair_id ==
          current_pair &
          sample_metadata_qc$condition ==
          "OA"
      ]
    
    preserved_expression <-
      expression_biological[
        ,
        preserved_column
      ]
    
    oa_expression <-
      expression_biological[
        ,
        oa_column
      ]
    
    paired_difference <-
      oa_expression -
      preserved_expression
    
    tibble::tibble(
      pair_id =
        current_pair,
      
      within_pair_correlation =
        stats::cor(
          preserved_expression,
          oa_expression,
          method = "pearson"
        ),
      
      mean_absolute_change =
        mean(
          abs(
            paired_difference
          )
        ),
      
      rms_change =
        sqrt(
          mean(
            paired_difference^2
          )
        )
    )
  }
) %>%
  dplyr::bind_rows()

pair_metadata_unique <- sample_metadata_qc %>%
  dplyr::distinct(
    pair_id,
    pair_number,
    age_years,
    sex
  )

pair_diagnostics <- pair_diagnostics %>%
  dplyr::left_join(
    pair_metadata_unique,
    by = "pair_id"
  ) %>%
  dplyr::mutate(
    
    correlation_z_score =
      as.numeric(
        scale(
          within_pair_correlation
        )
      ),
    
    rms_change_z_score =
      as.numeric(
        scale(
          rms_change
        )
      ),
    
    candidate_pair_outlier =
      correlation_z_score < -3 |
      rms_change_z_score > 3
  ) %>%
  dplyr::arrange(
    correlation_z_score
  )

cat("\nPair-level diagnostics:\n")

print(
  pair_diagnostics,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 17. Plot within-pair correlation
# ------------------------------------------------------------

pair_correlation_plot <- ggplot(
  pair_diagnostics,
  aes(
    x = reorder(
      pair_id,
      within_pair_correlation
    ),
    y = within_pair_correlation
  )
) +
  geom_point(
    size = 2.8
  ) +
  geom_hline(
    yintercept =
      mean(
        pair_diagnostics$within_pair_correlation
      ),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  coord_flip() +
  labs(
    title =
      "Within-patient expression correlation",
    
    subtitle =
      "Preserved versus OA-affected cartilage",
    
    x =
      "Patient pair",
    
    y =
      "Pearson correlation"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(pair_correlation_plot)

ggsave(
  filename =
    "results/figures/GSE57218_within_pair_correlation.pdf",
  plot = pair_correlation_plot,
  width = 7,
  height = 8
)

ggsave(
  filename =
    "results/figures/GSE57218_within_pair_correlation.tiff",
  plot = pair_correlation_plot,
  width = 7,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 18. Plot pair-level RMS change
# ------------------------------------------------------------

pair_rms_plot <- ggplot(
  pair_diagnostics,
  aes(
    x = reorder(
      pair_id,
      rms_change
    ),
    y = rms_change
  )
) +
  geom_point(
    size = 2.8
  ) +
  geom_hline(
    yintercept =
      mean(
        pair_diagnostics$rms_change
      ),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  coord_flip() +
  labs(
    title =
      "Magnitude of paired transcriptomic change",
    
    subtitle =
      "Root-mean-square difference between OA and preserved cartilage",
    
    x =
      "Patient pair",
    
    y =
      "RMS expression difference"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(pair_rms_plot)

ggsave(
  filename =
    "results/figures/GSE57218_pair_RMS_change.pdf",
  plot = pair_rms_plot,
  width = 7,
  height = 8
)

ggsave(
  filename =
    "results/figures/GSE57218_pair_RMS_change.tiff",
  plot = pair_rms_plot,
  width = 7,
  height = 8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 19. Membuat paired design matrix
# ------------------------------------------------------------

design_paired <- model.matrix(
  ~ pair_id + condition,
  data = sample_metadata_qc
)

rownames(design_paired) <-
  sample_metadata_qc$expression_column

design_unpaired_sensitivity <- model.matrix(
  ~ condition,
  data = sample_metadata_qc
)

rownames(
  design_unpaired_sensitivity
) <- sample_metadata_qc$expression_column

cat("\nNama koefisien paired design:\n")
print(colnames(design_paired))

cat("\nDimensi paired design:\n")
print(dim(design_paired))

cat("\nRank paired design:\n")
print(qr(design_paired)$rank)

cat("\nCondition number paired design:\n")
print(kappa(design_paired))

design_diagnostics <- tibble::tibble(
  
  model = c(
    "Primary paired model",
    "Unpaired sensitivity model"
  ),
  
  formula = c(
    "~ pair_id + condition",
    "~ condition"
  ),
  
  number_of_samples = c(
    nrow(design_paired),
    nrow(
      design_unpaired_sensitivity
    )
  ),
  
  number_of_parameters = c(
    ncol(design_paired),
    ncol(
      design_unpaired_sensitivity
    )
  ),
  
  rank = c(
    qr(design_paired)$rank,
    qr(
      design_unpaired_sensitivity
    )$rank
  ),
  
  full_rank = c(
    qr(design_paired)$rank ==
      ncol(design_paired),
    
    qr(
      design_unpaired_sensitivity
    )$rank ==
      ncol(
        design_unpaired_sensitivity
      )
  ),
  
  residual_degrees_of_freedom = c(
    nrow(design_paired) -
      qr(design_paired)$rank,
    
    nrow(
      design_unpaired_sensitivity
    ) -
      qr(
        design_unpaired_sensitivity
      )$rank
  ),
  
  condition_number = c(
    kappa(design_paired),
    
    kappa(
      design_unpaired_sensitivity
    )
  )
)

cat("\nDesign diagnostics:\n")

print(
  design_diagnostics,
  n = Inf,
  width = Inf
)

if (
  qr(design_paired)$rank !=
  ncol(design_paired)
) {
  stop(
    "Paired design matrix tidak full rank."
  )
}

if (
  !"conditionOA" %in%
  colnames(design_paired)
) {
  stop(
    "Koefisien conditionOA tidak ditemukan."
  )
}

# ------------------------------------------------------------
# 20. Ekspresi MMP13
# ------------------------------------------------------------

stopifnot(
  nrow(mmp13_annotation) == 1
)

mmp13_probe_id <- as.character(
  mmp13_annotation$probe_id
)

stopifnot(
  mmp13_probe_id %in%
    rownames(expression_primary)
)

mmp13_expression <- sample_metadata_qc %>%
  dplyr::mutate(
    MMP13_expression =
      as.numeric(
        expression_primary[
          mmp13_probe_id,
          expression_column
        ]
      )
  )

mmp13_pair_table <- mmp13_expression %>%
  dplyr::select(
    pair_id,
    pair_number,
    age_years,
    sex,
    condition,
    MMP13_expression
  ) %>%
  tidyr::pivot_wider(
    names_from =
      condition,
    
    values_from =
      MMP13_expression
  ) %>%
  dplyr::mutate(
    
    difference_OA_minus_Preserved =
      OA - Preserved,
    
    direction =
      dplyr::case_when(
        difference_OA_minus_Preserved > 0 ~
          "Higher in OA",
        
        difference_OA_minus_Preserved < 0 ~
          "Lower in OA",
        
        TRUE ~
          "No difference"
      )
  ) %>%
  dplyr::arrange(
    pair_number
  )

cat("\nPaired MMP13 differences:\n")

print(
  mmp13_pair_table,
  n = Inf,
  width = Inf
)

mmp13_difference_summary <- tibble::tibble(
  
  number_of_pairs =
    nrow(mmp13_pair_table),
  
  mean_preserved =
    mean(
      mmp13_pair_table$Preserved
    ),
  
  mean_OA =
    mean(
      mmp13_pair_table$OA
    ),
  
  mean_paired_difference =
    mean(
      mmp13_pair_table$
        difference_OA_minus_Preserved
    ),
  
  sd_paired_difference =
    sd(
      mmp13_pair_table$
        difference_OA_minus_Preserved
    ),
  
  median_paired_difference =
    median(
      mmp13_pair_table$
        difference_OA_minus_Preserved
    ),
  
  pairs_higher_in_OA =
    sum(
      mmp13_pair_table$
        difference_OA_minus_Preserved > 0
    ),
  
  pairs_lower_in_OA =
    sum(
      mmp13_pair_table$
        difference_OA_minus_Preserved < 0
    ),
  
  pairs_equal =
    sum(
      mmp13_pair_table$
        difference_OA_minus_Preserved == 0
    )
)

cat("\nRingkasan paired MMP13:\n")

print(
  mmp13_difference_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 21. Paired MMP13 expression plot
# ------------------------------------------------------------

mmp13_paired_plot <- ggplot(
  mmp13_expression,
  aes(
    x = display_condition,
    y = MMP13_expression,
    group = pair_id
  )
) +
  geom_line(
    alpha = 0.45,
    linewidth = 0.55
  ) +
  geom_point(
    aes(
      color = display_condition,
      shape = sex
    ),
    size = 2.8,
    alpha = 0.9
  ) +
  labs(
    title =
      "Paired MMP13 expression in GSE57218",
    
    subtitle = paste0(
      "Probe ",
      mmp13_probe_id,
      "; each line represents one patient"
    ),
    
    x = NULL,
    
    y =
      "Processed VST-like expression",
    
    color =
      "Cartilage condition",
    
    shape =
      "Sex"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(mmp13_paired_plot)

ggsave(
  filename =
    "results/figures/GSE57218_MMP13_paired_expression.pdf",
  plot = mmp13_paired_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE57218_MMP13_paired_expression.tiff",
  plot = mmp13_paired_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 22. MMP13 paired-difference distribution
# ------------------------------------------------------------

mmp13_difference_plot <- ggplot(
  mmp13_pair_table,
  aes(
    x = difference_OA_minus_Preserved
  )
) +
  geom_histogram(
    bins = 12,
    boundary = 0,
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    xintercept =
      mean(
        mmp13_pair_table$
          difference_OA_minus_Preserved
      ),
    linetype = "dotted",
    linewidth = 0.7
  ) +
  labs(
    title =
      "Distribution of paired MMP13 differences",
    
    subtitle =
      "Positive values indicate higher expression in OA-affected cartilage",
    
    x =
      "OA minus preserved MMP13 expression",
    
    y =
      "Number of patient pairs"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(mmp13_difference_plot)

ggsave(
  filename =
    "results/figures/GSE57218_MMP13_paired_differences.pdf",
  plot = mmp13_difference_plot,
  width = 7,
  height = 5
)

ggsave(
  filename =
    "results/figures/GSE57218_MMP13_paired_differences.tiff",
  plot = mmp13_difference_plot,
  width = 7,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 23. Membuat ringkasan QC
# ------------------------------------------------------------

qc_summary <- tibble::tibble(
  
  metric = c(
    "Number of primary samples",
    "Number of complete patient pairs",
    "Preserved samples",
    "OA-affected samples",
    "Total probes",
    "Annotated biological probes",
    "Top variable probes used in PCA",
    "PC1 variance explained",
    "PC2 variance explained",
    "Candidate sample outliers",
    "Candidate pair outliers",
    "Paired design parameters",
    "Paired design residual df",
    "MMP13 probe"
  ),
  
  value = c(
    ncol(expression_primary),
    
    dplyr::n_distinct(
      sample_metadata_qc$pair_id
    ),
    
    sum(
      sample_metadata_qc$condition ==
        "Preserved"
    ),
    
    sum(
      sample_metadata_qc$condition ==
        "OA"
    ),
    
    nrow(expression_primary),
    
    nrow(expression_biological),
    
    nrow(expression_top_variable),
    
    round(
      pca_variance_explained[1],
      3
    ),
    
    round(
      pca_variance_explained[2],
      3
    ),
    
    sum(
      sample_connectivity$
        candidate_sample_outlier
    ),
    
    sum(
      pair_diagnostics$
        candidate_pair_outlier
    ),
    
    ncol(design_paired),
    
    nrow(design_paired) -
      qr(design_paired)$rank,
    
    mmp13_probe_id
  )
)

print(
  qc_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 24. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  sample_metadata_qc,
  file =
    "data_processed/GSE57218_metadata_paired_QC.rds"
)

saveRDS(
  expression_biological,
  file =
    "data_processed/GSE57218_expression_annotated_biological.rds"
)

saveRDS(
  design_paired,
  file =
    "data_processed/GSE57218_design_paired.rds"
)

saveRDS(
  design_unpaired_sensitivity,
  file =
    "data_processed/GSE57218_design_unpaired_sensitivity.rds"
)

saveRDS(
  pca_result,
  file =
    "data_processed/GSE57218_PCA_result.rds"
)

saveRDS(
  mds_result,
  file =
    "data_processed/GSE57218_MDS_result.rds"
)

saveRDS(
  correlation_matrix,
  file =
    "data_processed/GSE57218_sample_correlation_matrix.rds"
)

saveRDS(
  pair_diagnostics,
  file =
    "data_processed/GSE57218_pair_diagnostics.rds"
)

saveRDS(
  mmp13_pair_table,
  file =
    "data_processed/GSE57218_MMP13_pair_table.rds"
)

# ------------------------------------------------------------
# 25. Menyimpan Excel workbook
# ------------------------------------------------------------

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
    
    Distribution_Summary =
      as.data.frame(
        distribution_summary
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
    
    Pair_Diagnostics =
      as.data.frame(
        pair_diagnostics
      ),
    
    Pair_Covariate_Check =
      as.data.frame(
        pair_covariate_check
      ),
    
    Design_Diagnostics =
      as.data.frame(
        design_diagnostics
      ),
    
    MMP13_Sample_Expression =
      as.data.frame(
        mmp13_expression
      ),
    
    MMP13_Paired_Differences =
      as.data.frame(
        mmp13_pair_table
      ),
    
    MMP13_Difference_Summary =
      as.data.frame(
        mmp13_difference_summary
      )
  ),
  
  file =
    "results/tables/GSE57218_limma_paired_QC_results.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 26. Menyimpan session information
# ------------------------------------------------------------

sink(
  "results/tables/GSE57218_QC_sessionInfo.txt"
)

print(sessionInfo())

sink()

# ------------------------------------------------------------
# 27. Pesan akhir
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE57218 PAIRED QC SELESAI")
message("Jumlah sampel             : ",
        ncol(expression_primary))
message("Jumlah pasangan           : ",
        dplyr::n_distinct(
          sample_metadata_qc$pair_id
        ))
message("Probe biologis teranotasi : ",
        nrow(expression_biological))
message("PC1 variance              : ",
        round(pca_variance_explained[1], 2),
        "%")
message("PC2 variance              : ",
        round(pca_variance_explained[2], 2),
        "%")
message("Candidate sample outlier  : ",
        sum(
          sample_connectivity$
            candidate_sample_outlier
        ))
message("Candidate pair outlier    : ",
        sum(
          pair_diagnostics$
            candidate_pair_outlier
        ))
message("Paired design residual df : ",
        nrow(design_paired) -
          qr(design_paired)$rank)
message("MMP13 mean paired change  : ",
        round(
          mmp13_difference_summary$
            mean_paired_difference,
          4
        ))
message("MMP13 higher in OA pairs  : ",
        mmp13_difference_summary$
          pairs_higher_in_OA)
message("MMP13 lower in OA pairs   : ",
        mmp13_difference_summary$
          pairs_lower_in_OA)
message("============================================")