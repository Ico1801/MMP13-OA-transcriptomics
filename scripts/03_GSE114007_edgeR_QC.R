# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 03_GSE114007_edgeR_QC.R
# PURPOSE : Filtering, TMM normalization, MDS dan PCA QC
# DATASET : GSE114007
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Memuat paket
# ------------------------------------------------------------

required_packages <- c(
  "edgeR",
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

library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(openxlsx)

message("Semua paket untuk QC berhasil dimuat.")
# ------------------------------------------------------------
# 2. Memuat data
# ------------------------------------------------------------

count_matrix <- readRDS(
  "data_processed/GSE114007_raw_count_matrix.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE114007_sample_metadata_with_platform.rds"
)

cat("\nDimensi count matrix:\n")
print(dim(count_matrix))

cat("\nDimensi metadata:\n")
print(dim(sample_metadata))

cat("\nDistribusi kelompok:\n")
print(table(sample_metadata$group))

cat("\nDistribusi platform:\n")
print(table(sample_metadata$platform_id))
# ------------------------------------------------------------
# 3. Validasi count matrix dan metadata
# ------------------------------------------------------------

stopifnot(
  ncol(count_matrix) == nrow(sample_metadata)
)

stopifnot(
  identical(
    colnames(count_matrix),
    sample_metadata$sample_id
  )
)

stopifnot(
  !anyNA(count_matrix)
)

stopifnot(
  !any(count_matrix < 0)
)

stopifnot(
  all(
    count_matrix == round(count_matrix)
  )
)

stopifnot(
  !anyNA(sample_metadata$group)
)

stopifnot(
  !anyNA(sample_metadata$platform_id)
)

message("Count matrix dan metadata sudah sesuai.")
# ------------------------------------------------------------
# 4. Membuat label pendek sampel
# ------------------------------------------------------------

sample_metadata <- sample_metadata %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    group_number = dplyr::row_number(),
    plot_label = dplyr::if_else(
      group == "Control",
      paste0("C", group_number),
      paste0("OA", group_number)
    )
  ) %>%
  dplyr::ungroup()

print(
  sample_metadata %>%
    dplyr::select(
      sample_id,
      plot_label,
      group,
      platform_id
    )
)
# ------------------------------------------------------------
# 5. Menghitung raw library size
# ------------------------------------------------------------

raw_library_size <- colSums(count_matrix)

sample_qc <- sample_metadata %>%
  dplyr::mutate(
    raw_library_size = as.numeric(
      raw_library_size[sample_id]
    ),
    raw_library_size_million =
      raw_library_size / 1e6
  )

cat("\nRingkasan raw library size:\n")
print(
  summary(sample_qc$raw_library_size_million)
)

cat("\nRaw library size setiap sampel:\n")
print(
  sample_qc %>%
    dplyr::select(
      plot_label,
      sample_id,
      group,
      platform_id,
      raw_library_size_million
    )
)
# ------------------------------------------------------------
# 6. Grafik raw library size
# ------------------------------------------------------------

library_size_plot <- ggplot(
  sample_qc,
  aes(
    x = reorder(plot_label, raw_library_size_million),
    y = raw_library_size_million,
    fill = group
  )
) +
  geom_col(
    width = 0.75
  ) +
  coord_flip() +
  facet_grid(
    platform_id ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title = "Raw library sizes of GSE114007 samples",
    subtitle = "Samples are separated according to sequencing platform",
    x = "Sample",
    y = "Raw library size (million reads)",
    fill = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "top"
  )

print(library_size_plot)
# ------------------------------------------------------------
# 7. Memastikan level referensi
# ------------------------------------------------------------

sample_metadata$group <- factor(
  sample_metadata$group,
  levels = c("Control", "OA")
)

sample_metadata$platform_id <- factor(
  sample_metadata$platform_id,
  levels = c("GPL11154", "GPL18573")
)

# ------------------------------------------------------------
# 8. Membuat design matrix
# ------------------------------------------------------------

design <- model.matrix(
  ~ platform_id + group,
  data = sample_metadata
)

rownames(design) <- sample_metadata$sample_id

cat("\nDesign matrix:\n")
print(design)

cat("\nNama koefisien:\n")
print(colnames(design))

cat("\nRank design matrix:\n")
print(qr(design)$rank)

cat("\nJumlah kolom design matrix:\n")
print(ncol(design))
# ------------------------------------------------------------
# 9. Membuat DGEList
# ------------------------------------------------------------

dge_original <- edgeR::DGEList(
  counts = count_matrix,
  group = sample_metadata$group
)

# Tambahkan informasi sampel
dge_original$samples$sample_id <-
  sample_metadata$sample_id

dge_original$samples$plot_label <-
  sample_metadata$plot_label

dge_original$samples$platform_id <-
  sample_metadata$platform_id

dge_original$samples$geo_accession <-
  sample_metadata$geo_accession

cat("\nObjek DGEList awal:\n")
print(dge_original)
# ------------------------------------------------------------
# 10. Menyaring gen berekspresi rendah
# ------------------------------------------------------------

keep_genes <- edgeR::filterByExpr(
  dge_original,
  design = design
)

cat("\nJumlah gen sebelum filtering:\n")
print(nrow(dge_original))

cat("\nJumlah gen yang dipertahankan:\n")
print(sum(keep_genes))

cat("\nJumlah gen yang dibuang:\n")
print(sum(!keep_genes))

cat("\nPersentase gen yang dipertahankan:\n")
print(
  round(
    100 * mean(keep_genes),
    digits = 2
  )
)
# ------------------------------------------------------------
# 11. Status filtering MMP13
# ------------------------------------------------------------

mmp13_in_original <-
  "MMP13" %in% rownames(dge_original)

mmp13_pass_filter <-
  if (mmp13_in_original) {
    keep_genes[
      match(
        "MMP13",
        rownames(dge_original)
      )
    ]
  } else {
    FALSE
  }

cat("\nMMP13 tersedia sebelum filtering:\n")
print(mmp13_in_original)

cat("\nMMP13 lolos filtering:\n")
print(mmp13_pass_filter)
# ------------------------------------------------------------
# 12. Normalisasi TMM
# ------------------------------------------------------------

dge_normalized <- edgeR::calcNormFactors(
  dge_filtered,
  method = "TMM"
)

cat("\nTMM normalization factors:\n")
print(
  dge_normalized$samples[
    ,
    c(
      "sample_id",
      "group",
      "platform_id",
      "lib.size",
      "norm.factors"
    )
  ]
)
# ------------------------------------------------------------
# 13. Menghitung normalized log2-CPM
# ------------------------------------------------------------

log_cpm <- edgeR::cpm(
  dge_normalized,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

cat("\nDimensi normalized log-CPM:\n")
print(dim(log_cpm))

cat("\nAda nilai non-finite:\n")
print(
  any(!is.finite(log_cpm))
)

cat("\nRentang normalized log-CPM:\n")
print(
  range(
    log_cpm,
    na.rm = TRUE
  )
)
# ------------------------------------------------------------
# 14. Menghitung MDS
# ------------------------------------------------------------

mds_result <- plotMDS(
  dge_normalized,
  top = 500,
  gene.selection = "pairwise",
  plot = FALSE
)

mds_table <- sample_metadata %>%
  dplyr::mutate(
    MDS_dimension_1 =
      as.numeric(mds_result$x),
    
    MDS_dimension_2 =
      as.numeric(mds_result$y)
  )

print(
  mds_table %>%
    dplyr::select(
      sample_id,
      plot_label,
      group,
      platform_id,
      MDS_dimension_1,
      MDS_dimension_2
    )
)
# ------------------------------------------------------------
# 15. PCA dari normalized log-CPM
# ------------------------------------------------------------

pca_result <- prcomp(
  t(log_cpm),
  center = TRUE,
  scale. = FALSE
)

variance_explained <- (
  pca_result$sdev^2 /
    sum(pca_result$sdev^2)
) * 100

pca_table <- sample_metadata %>%
  dplyr::mutate(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2]
  )

cat("\nVariance explained PC1 dan PC2:\n")
print(
  round(
    variance_explained[1:2],
    digits = 2
  )
)
# ------------------------------------------------------------
# 16. Platform-adjusted PCA untuk visualisasi saja
# ------------------------------------------------------------

group_design_for_visualization <- model.matrix(
  ~ group,
  data = sample_metadata
)

cat("\nDimensi design untuk visualisasi:\n")
print(dim(group_design_for_visualization))

cat("\nNama kolom design untuk visualisasi:\n")
print(colnames(group_design_for_visualization))

stopifnot(
  nrow(group_design_for_visualization) == ncol(log_cpm)
)

# Menghilangkan efek platform sambil mempertahankan efek group
log_cpm_platform_adjusted <- limma::removeBatchEffect(
  log_cpm,
  batch = sample_metadata$platform_id,
  design = group_design_for_visualization
)

cat("\nDimensi log-CPM setelah platform adjustment:\n")
print(dim(log_cpm_platform_adjusted))

cat("\nAda nilai non-finite setelah adjustment:\n")
print(any(!is.finite(log_cpm_platform_adjusted)))

# PCA setelah platform adjustment
pca_adjusted_result <- prcomp(
  t(log_cpm_platform_adjusted),
  center = TRUE,
  scale. = FALSE
)

adjusted_variance_explained <- (
  pca_adjusted_result$sdev^2 /
    sum(pca_adjusted_result$sdev^2)
) * 100

pca_adjusted_table <- sample_metadata %>%
  dplyr::mutate(
    PC1 = pca_adjusted_result$x[, 1],
    PC2 = pca_adjusted_result$x[, 2]
  )

cat("\nVariance explained setelah platform adjustment:\n")
print(
  round(
    adjusted_variance_explained[1:2],
    digits = 2
  )
)

pca_adjusted_plot <- ggplot(
  pca_adjusted_table,
  aes(
    x = PC1,
    y = PC2,
    color = group,
    shape = platform_id
  )
) +
  geom_point(
    size = 3.2,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = plot_label),
    vjust = -0.8,
    size = 2.8,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "Platform-adjusted PCA of GSE114007",
    subtitle =
      "Platform effect removed for visualization while preserving disease group",
    x = paste0(
      "PC1 (",
      round(adjusted_variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(adjusted_variance_explained[2], 1),
      "%)"
    ),
    color = "Group",
    shape = "Platform"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(pca_adjusted_plot)

ggsave(
  filename =
    "results/figures/GSE114007_PCA_platform_adjusted_visualization.pdf",
  plot = pca_adjusted_plot,
  width = 7,
  height = 5.5
)

ggsave(
  filename =
    "results/figures/GSE114007_PCA_platform_adjusted_visualization.tiff",
  plot = pca_adjusted_plot,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)
# ------------------------------------------------------------
# 17. Ekspresi normalized MMP13
# ------------------------------------------------------------

if ("MMP13" %in% rownames(log_cpm)) {
  
  mmp13_expression <- sample_metadata %>%
    dplyr::mutate(
      MMP13_log2_CPM = as.numeric(
        log_cpm["MMP13", sample_id]
      )
    )
  
  cat("\nEkspresi normalized MMP13:\n")
  
  print(
    mmp13_expression %>%
      dplyr::select(
        sample_id,
        plot_label,
        group,
        platform_id,
        MMP13_log2_CPM
      ),
    n = Inf
  )
  
  cat("\nRingkasan MMP13 berdasarkan kelompok:\n")
  
  mmp13_group_summary <- mmp13_expression %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_log2_CPM = mean(MMP13_log2_CPM),
      sd_log2_CPM = sd(MMP13_log2_CPM),
      median_log2_CPM = median(MMP13_log2_CPM),
      minimum = min(MMP13_log2_CPM),
      maximum = max(MMP13_log2_CPM),
      .groups = "drop"
    )
  
  print(mmp13_group_summary)
  
} else {
  
  stop("MMP13 tidak ditemukan dalam matriks log-CPM.")
}
# ------------------------------------------------------------
# 18. Ringkasan filtering
# ------------------------------------------------------------

filtering_summary <- data.frame(
  metric = c(
    "Genes before filtering",
    "Genes retained",
    "Genes removed",
    "Percentage retained",
    "MMP13 present before filtering",
    "MMP13 passed filtering"
  ),
  value = c(
    nrow(dge_original),
    sum(keep_genes),
    sum(!keep_genes),
    round(100 * mean(keep_genes), 2),
    mmp13_in_original,
    as.logical(mmp13_pass_filter)
  )
)

print(filtering_summary)

# ------------------------------------------------------------
# 19. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  dge_normalized,
  "data_processed/GSE114007_edgeR_TMM_normalized_DGEList.rds"
)

saveRDS(
  log_cpm,
  "data_processed/GSE114007_TMM_logCPM.rds"
)

saveRDS(
  log_cpm_platform_adjusted,
  "data_processed/GSE114007_platform_adjusted_logCPM_visualization.rds"
)

saveRDS(
  design,
  "data_processed/GSE114007_design_matrix.rds"
)

saveRDS(
  keep_genes,
  "data_processed/GSE114007_filterByExpr_keep_genes.rds"
)

saveRDS(
  pca_result,
  "data_processed/GSE114007_PCA_result.rds"
)

saveRDS(
  pca_adjusted_result,
  "data_processed/GSE114007_adjusted_PCA_result.rds"
)

saveRDS(
  mds_result,
  "data_processed/GSE114007_MDS_result.rds"
)