# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 06_GSE117999_limma_QC.R
# PURPOSE : Quality control processed microarray GSE117999
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
library(openxlsx)

message("Semua paket QC GSE117999 berhasil dimuat.")

# ------------------------------------------------------------
# 2. Memastikan folder tersedia
# ------------------------------------------------------------

dir.create(
  "results/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "data_processed",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 3. Memuat data
# ------------------------------------------------------------

expression_matrix <- readRDS(
  "data_processed/GSE117999_expression_included_20_samples.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE117999_metadata_included_20_samples.rds"
)

feature_annotation <- readRDS(
  "data_processed/GSE117999_feature_annotation_original.rds"
)

mmp13_annotation <- readRDS(
  "data_processed/GSE117999_MMP13_probe_annotation.rds"
)

cat("\nDimensi expression matrix:\n")
print(dim(expression_matrix))

cat("\nDimensi metadata:\n")
print(dim(sample_metadata))

cat("\nDimensi feature annotation:\n")
print(dim(feature_annotation))

# ------------------------------------------------------------
# 4. Validasi data
# ------------------------------------------------------------

stopifnot(
  ncol(expression_matrix) ==
    nrow(sample_metadata)
)

stopifnot(
  identical(
    colnames(expression_matrix),
    sample_metadata$expression_column
  )
)

stopifnot(
  !anyNA(expression_matrix)
)

stopifnot(
  all(
    is.finite(expression_matrix)
  )
)

stopifnot(
  nrow(mmp13_annotation) == 1
)

stopifnot(
  mmp13_annotation$probe_id %in%
    rownames(expression_matrix)
)

message("Expression matrix dan metadata berhasil divalidasi.")

# ------------------------------------------------------------
# 5. Menyiapkan metadata QC
# ------------------------------------------------------------

sample_metadata_qc <- sample_metadata %>%
  dplyr::mutate(
    age_years = as.numeric(
      `age (years):ch1`
    ),
    
    bmi = as.numeric(
      `bmi (kg/m2):ch1`
    ),
    
    sex = factor(
      as.character(`Sex:ch1`)
    ),
    
    group = factor(
      group,
      levels = c("Control", "OA")
    ),
    
    display_group = factor(
      dplyr::if_else(
        group == "Control",
        "Non-OA (APM)",
        "OA"
      ),
      levels = c(
        "Non-OA (APM)",
        "OA"
      )
    )
  )

stopifnot(
  !anyNA(sample_metadata_qc$age_years)
)

stopifnot(
  !anyNA(sample_metadata_qc$bmi)
)

stopifnot(
  !anyNA(sample_metadata_qc$sex)
)

cat("\nMetadata QC:\n")

print(
  sample_metadata_qc %>%
    dplyr::select(
      expression_column,
      subject_code,
      plot_label,
      group,
      age_years,
      bmi,
      sex
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 6. Memeriksa kesamaan distribusi antar-array
# ------------------------------------------------------------

quantile_probabilities <- c(
  0,
  0.25,
  0.50,
  0.75,
  1
)

sample_quantiles <- apply(
  expression_matrix,
  MARGIN = 2,
  FUN = stats::quantile,
  probs = quantile_probabilities,
  na.rm = TRUE
)

quantile_range_across_samples <- apply(
  sample_quantiles,
  MARGIN = 1,
  FUN = function(x) {
    max(x) - min(x)
  }
)

cat("\nRentang quantile antarsampel:\n")
print(quantile_range_across_samples)

distribution_summary <- data.frame(
  quantile = rownames(sample_quantiles),
  minimum_across_samples =
    apply(sample_quantiles, 1, min),
  maximum_across_samples =
    apply(sample_quantiles, 1, max),
  range_across_samples =
    quantile_range_across_samples
)

print(distribution_summary)

# ------------------------------------------------------------
# 7. Menyelaraskan anotasi dengan expression matrix
# ------------------------------------------------------------

annotation_match <- match(
  rownames(expression_matrix),
  feature_annotation$probe_id
)

if (anyNA(annotation_match)) {
  stop(
    "Ada probe expression matrix yang tidak ditemukan ",
    "pada feature annotation."
  )
}

annotation_aligned <- feature_annotation[
  annotation_match,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(expression_matrix),
    annotation_aligned$probe_id
  )
)

control_type <- toupper(
  trimws(
    as.character(
      annotation_aligned$CONTROL_TYPE
    )
  )
)

gene_symbol <- trimws(
  as.character(
    annotation_aligned$GENE_SYMBOL
  )
)

keep_non_control <- (
  is.na(control_type) |
    control_type == "" |
    control_type == "FALSE" |
    control_type == "0"
)

keep_annotated <- (
  !is.na(gene_symbol) &
    gene_symbol != ""
)

keep_biological_annotated <-
  keep_non_control &
  keep_annotated

cat("\nJumlah semua probe:\n")
print(nrow(expression_matrix))

cat("\nJumlah probe non-control:\n")
print(sum(keep_non_control))

cat("\nJumlah probe dengan gene symbol:\n")
print(sum(keep_annotated))

cat("\nJumlah probe biologis dan teranotasi:\n")
print(sum(keep_biological_annotated))

# ------------------------------------------------------------
# 8. Matriks untuk QC multivariat
# ------------------------------------------------------------

expression_biological <- expression_matrix[
  keep_biological_annotated,
  ,
  drop = FALSE
]

probe_variance <- apply(
  expression_biological,
  MARGIN = 1,
  FUN = stats::var
)

probe_variance[!is.finite(probe_variance)] <- 0

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
# 9. Boxplot distribusi expression
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE117999_expression_boxplot.pdf",
  width = 9,
  height = 5.5
)

boxplot(
  expression_matrix,
  names = sample_metadata_qc$plot_label,
  outline = FALSE,
  las = 2,
  main =
    "Processed expression distributions in GSE117999",
  xlab = "Sample",
  ylab = "Processed log2 expression"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE117999_expression_boxplot.tiff",
  width = 9,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

boxplot(
  expression_matrix,
  names = sample_metadata_qc$plot_label,
  outline = FALSE,
  las = 2,
  main =
    "Processed expression distributions in GSE117999",
  xlab = "Sample",
  ylab = "Processed log2 expression"
)

dev.off()

# ------------------------------------------------------------
# 10. PCA
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
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  )

cat("\nVariance explained PCA:\n")

print(
  round(
    pca_variance_explained[1:5],
    2
  )
)

pca_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    color = display_group,
    shape = sex
  )
) +
  geom_point(
    size = 3.4,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = plot_label),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of GSE117999",
    subtitle =
      "Top 5,000 variable annotated biological probes",
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
    color = "Disease group",
    shape = "Sex"
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
    "results/figures/GSE117999_PCA.pdf",
  plot = pca_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE117999_PCA.tiff",
  plot = pca_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 11. MDS
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

mds_plot <- ggplot(
  mds_table,
  aes(
    x = MDS_dimension_1,
    y = MDS_dimension_2,
    color = display_group,
    shape = sex
  )
) +
  geom_point(
    size = 3.4,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = plot_label),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "MDS plot of GSE117999",
    subtitle =
      "Top 500 pairwise variable biological probes",
    x = "Leading log-fold-change dimension 1",
    y = "Leading log-fold-change dimension 2",
    color = "Disease group",
    shape = "Sex"
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
    "results/figures/GSE117999_MDS.pdf",
  plot = mds_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE117999_MDS.tiff",
  plot = mds_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 12. Sample-to-sample correlation
# ------------------------------------------------------------

correlation_matrix <- stats::cor(
  expression_biological,
  method = "pearson"
)

rownames(correlation_matrix) <-
  sample_metadata_qc$plot_label

colnames(correlation_matrix) <-
  sample_metadata_qc$plot_label

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
  as.table(
    correlation_matrix
  ),
  stringsAsFactors = FALSE
)

colnames(correlation_long) <- c(
  "Sample_1",
  "Sample_2",
  "Correlation"
)

correlation_long$Sample_1 <- factor(
  correlation_long$Sample_1,
  levels = clustered_labels
)

correlation_long$Sample_2 <- factor(
  correlation_long$Sample_2,
  levels = rev(clustered_labels)
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
      "Sample-to-sample expression correlation",
    x = NULL,
    y = NULL,
    fill = "Pearson\ncorrelation"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1
      ),
    panel.grid =
      element_blank()
  )

print(correlation_plot)

ggsave(
  filename =
    "results/figures/GSE117999_sample_correlation.pdf",
  plot = correlation_plot,
  width = 7,
  height = 6.5
)

ggsave(
  filename =
    "results/figures/GSE117999_sample_correlation.tiff",
  plot = correlation_plot,
  width = 7,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 13. Hierarchical clustering dendrogram
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE117999_sample_dendrogram.pdf",
  width = 9,
  height = 5.5
)

plot(
  sample_clustering,
  labels = sample_metadata_qc$plot_label[
    match(
      sample_clustering$labels,
      sample_metadata_qc$plot_label
    )
  ],
  main =
    "Hierarchical clustering of GSE117999 samples",
  xlab = "Sample",
  sub = "",
  ylab = "1 - Pearson correlation"
)

dev.off()

# ------------------------------------------------------------
# 14. Sample connectivity
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

sample_connectivity <- sample_metadata_qc %>%
  dplyr::mutate(
    mean_correlation_to_others =
      mean_correlation_to_others[
        plot_label
      ],
    
    connectivity_z_score =
      connectivity_z_score[
        match(
          plot_label,
          names(
            mean_correlation_to_others
          )
        )
      ],
    
    candidate_outlier =
      connectivity_z_score < -3
  )

cat("\nSample connectivity:\n")

print(
  sample_connectivity %>%
    dplyr::select(
      plot_label,
      subject_code,
      group,
      mean_correlation_to_others,
      connectivity_z_score,
      candidate_outlier
    ) %>%
    dplyr::arrange(
      connectivity_z_score
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 15. Ringkasan umur, BMI, dan jenis kelamin
# ------------------------------------------------------------

covariate_summary <- sample_metadata_qc %>%
  dplyr::group_by(
    display_group
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    mean_age =
      mean(age_years),
    
    sd_age =
      sd(age_years),
    
    median_age =
      median(age_years),
    
    mean_bmi =
      mean(bmi),
    
    sd_bmi =
      sd(bmi),
    
    median_bmi =
      median(bmi),
    
    female =
      sum(sex == "Female"),
    
    male =
      sum(sex == "Male"),
    
    .groups = "drop"
  )

cat("\nRingkasan covariate berdasarkan kelompok:\n")

print(
  covariate_summary,
  n = Inf,
  width = Inf
)

cat("\nTabel sex × group:\n")

sex_group_table <- table(
  sample_metadata_qc$sex,
  sample_metadata_qc$display_group
)

print(sex_group_table)

group_binary <- as.numeric(
  sample_metadata_qc$group == "OA"
)

covariate_group_correlations <- data.frame(
  covariate = c(
    "Age",
    "BMI"
  ),
  
  correlation_with_OA_group = c(
    cor(
      sample_metadata_qc$age_years,
      group_binary
    ),
    
    cor(
      sample_metadata_qc$bmi,
      group_binary
    )
  )
)

cat("\nKorelasi covariate dengan kelompok OA:\n")
print(covariate_group_correlations)

# ------------------------------------------------------------
# 16. Memeriksa kelayakan design matrix
# ------------------------------------------------------------

design_primary <- model.matrix(
  ~ group,
  data = sample_metadata_qc
)

rownames(design_primary) <-
  sample_metadata_qc$expression_column

sample_metadata_qc <- sample_metadata_qc %>%
  dplyr::mutate(
    age_z = as.numeric(
      scale(age_years)
    ),
    
    bmi_z = as.numeric(
      scale(bmi)
    )
  )

design_adjusted <- model.matrix(
  ~ age_z + bmi_z + sex + group,
  data = sample_metadata_qc
)

rownames(design_adjusted) <-
  sample_metadata_qc$expression_column

design_diagnostics <- data.frame(
  model = c(
    "Primary: group only",
    "Sensitivity: age + BMI + sex + group"
  ),
  
  number_of_rows = c(
    nrow(design_primary),
    nrow(design_adjusted)
  ),
  
  number_of_columns = c(
    ncol(design_primary),
    ncol(design_adjusted)
  ),
  
  rank = c(
    qr(design_primary)$rank,
    qr(design_adjusted)$rank
  ),
  
  full_rank = c(
    qr(design_primary)$rank ==
      ncol(design_primary),
    
    qr(design_adjusted)$rank ==
      ncol(design_adjusted)
  ),
  
  condition_number = c(
    kappa(design_primary),
    kappa(design_adjusted)
  )
)

cat("\nDesign diagnostics:\n")
print(design_diagnostics)

if (
  qr(design_primary)$rank !=
  ncol(design_primary)
) {
  stop(
    "Primary design matrix tidak full rank."
  )
}

if (
  qr(design_adjusted)$rank !=
  ncol(design_adjusted)
) {
  warning(
    "Adjusted design matrix tidak full rank."
  )
}

# ------------------------------------------------------------
# 17. Ekspresi MMP13
# ------------------------------------------------------------

mmp13_probe_id <- as.character(
  mmp13_annotation$probe_id
)

mmp13_expression <- sample_metadata_qc %>%
  dplyr::mutate(
    MMP13_log2_expression =
      as.numeric(
        expression_matrix[
          mmp13_probe_id,
          expression_column
        ]
      )
  )

cat("\nEkspresi MMP13 seluruh sampel:\n")

print(
  mmp13_expression %>%
    dplyr::select(
      expression_column,
      subject_code,
      plot_label,
      display_group,
      age_years,
      bmi,
      sex,
      MMP13_log2_expression
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

mmp13_group_summary <- mmp13_expression %>%
  dplyr::group_by(
    display_group
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    
    mean =
      mean(
        MMP13_log2_expression
      ),
    
    sd =
      sd(
        MMP13_log2_expression
      ),
    
    median =
      median(
        MMP13_log2_expression
      ),
    
    minimum =
      min(
        MMP13_log2_expression
      ),
    
    maximum =
      max(
        MMP13_log2_expression
      ),
    
    .groups = "drop"
  )

cat("\nRingkasan MMP13:\n")

print(
  mmp13_group_summary,
  n = Inf,
  width = Inf
)

mmp13_plot <- ggplot(
  mmp13_expression,
  aes(
    x = display_group,
    y = MMP13_log2_expression,
    color = display_group,
    shape = sex
  )
) +
  geom_boxplot(
    aes(fill = display_group),
    width = 0.55,
    alpha = 0.18,
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 3,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = plot_label),
    position =
      position_jitter(
        width = 0.10,
        height = 0
      ),
    vjust = -0.8,
    size = 2.8,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "MMP13 expression in GSE117999",
    subtitle =
      "Probe A_33_P3221203",
    x = NULL,
    y =
      expression(
        "Processed " *
          log[2] *
          " expression"
      ),
    color = "Disease group",
    shape = "Sex"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank(),
    legend.position =
      "top"
  )

print(mmp13_plot)

ggsave(
  filename =
    "results/figures/GSE117999_MMP13_expression.pdf",
  plot = mmp13_plot,
  width = 5.8,
  height = 5
)

ggsave(
  filename =
    "results/figures/GSE117999_MMP13_expression.tiff",
  plot = mmp13_plot,
  width = 5.8,
  height = 5,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 18. Menyimpan hasil
# ------------------------------------------------------------

saveRDS(
  sample_metadata_qc,
  file =
    "data_processed/GSE117999_metadata_QC.rds"
)

saveRDS(
  expression_biological,
  file =
    "data_processed/GSE117999_expression_biological_annotated.rds"
)

saveRDS(
  design_primary,
  file =
    "data_processed/GSE117999_design_primary.rds"
)

saveRDS(
  design_adjusted,
  file =
    "data_processed/GSE117999_design_adjusted_sensitivity.rds"
)

saveRDS(
  pca_result,
  file =
    "data_processed/GSE117999_PCA_result.rds"
)

saveRDS(
  mds_result,
  file =
    "data_processed/GSE117999_MDS_result.rds"
)

saveRDS(
  correlation_matrix,
  file =
    "data_processed/GSE117999_sample_correlation_matrix.rds"
)

# ------------------------------------------------------------
# 19. Excel workbook
# ------------------------------------------------------------

qc_summary <- data.frame(
  metric = c(
    "Number of samples",
    "Number of non-OA APM samples",
    "Number of OA samples",
    "Total probes",
    "Biological annotated probes",
    "Top variable probes used in PCA",
    "PC1 variance explained",
    "PC2 variance explained",
    "Candidate connectivity outliers",
    "MMP13 probe"
  ),
  
  value = c(
    ncol(expression_matrix),
    
    sum(
      sample_metadata_qc$group ==
        "Control"
    ),
    
    sum(
      sample_metadata_qc$group ==
        "OA"
    ),
    
    nrow(expression_matrix),
    
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
      sample_connectivity$candidate_outlier
    ),
    
    mmp13_probe_id
  )
)

openxlsx::write.xlsx(
  list(
    QC_Summary =
      qc_summary,
    
    Sample_Metadata =
      as.data.frame(
        sample_metadata_qc
      ),
    
    Distribution_Summary =
      distribution_summary,
    
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
    
    Covariate_Summary =
      as.data.frame(
        covariate_summary
      ),
    
    Covariate_Correlations =
      covariate_group_correlations,
    
    Design_Diagnostics =
      design_diagnostics,
    
    MMP13_Expression =
      as.data.frame(
        mmp13_expression
      ),
    
    MMP13_Group_Summary =
      as.data.frame(
        mmp13_group_summary
      )
  ),
  
  file =
    "results/tables/GSE117999_limma_QC_results.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 20. Session information
# ------------------------------------------------------------

sink(
  "results/tables/GSE117999_QC_sessionInfo.txt"
)

print(sessionInfo())

sink()

# ------------------------------------------------------------
# 21. Pesan akhir
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE117999 LIMMA QC SELESAI")
message("Jumlah sampel             : ",
        ncol(expression_matrix))
message("Probe biologis teranotasi : ",
        nrow(expression_biological))
message("PC1 variance              : ",
        round(pca_variance_explained[1], 2),
        "%")
message("PC2 variance              : ",
        round(pca_variance_explained[2], 2),
        "%")
message("Connectivity outlier      : ",
        sum(sample_connectivity$candidate_outlier))
message("MMP13 probe               : ",
        mmp13_probe_id)
message("============================================")