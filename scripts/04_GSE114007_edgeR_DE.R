# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 04_GSE114007_edgeR_DE.R
# PURPOSE : Differential-expression analysis dan ekstraksi
#           MMP13 effect size untuk meta-analysis
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

message("Semua paket differential expression berhasil dimuat.")
# ------------------------------------------------------------
# 2. Memuat objek hasil QC
# ------------------------------------------------------------

dge_normalized <- readRDS(
  "data_processed/GSE114007_edgeR_TMM_normalized_DGEList.rds"
)

design <- readRDS(
  "data_processed/GSE114007_design_matrix.rds"
)

log_cpm <- readRDS(
  "data_processed/GSE114007_TMM_logCPM.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE114007_sample_metadata_with_platform.rds"
)

cat("\nDimensi DGEList:\n")
print(dim(dge_normalized))

cat("\nDimensi design matrix:\n")
print(dim(design))

cat("\nDimensi log-CPM:\n")
print(dim(log_cpm))

cat("\nNama koefisien design:\n")
print(colnames(design))

cat("\nJumlah sampel setiap kelompok:\n")
print(table(sample_metadata$group))

cat("\nDistribusi group × platform:\n")
print(
  table(
    sample_metadata$platform_id,
    sample_metadata$group
  )
)
# ------------------------------------------------------------
# 3. Validasi final
# ------------------------------------------------------------

stopifnot(
  ncol(dge_normalized) == nrow(sample_metadata)
)

stopifnot(
  nrow(design) == ncol(dge_normalized)
)

stopifnot(
  identical(
    colnames(dge_normalized),
    sample_metadata$sample_id
  )
)

stopifnot(
  identical(
    rownames(design),
    sample_metadata$sample_id
  )
)

stopifnot(
  identical(
    rownames(log_cpm),
    rownames(dge_normalized)
  )
)

stopifnot(
  identical(
    colnames(log_cpm),
    colnames(dge_normalized)
  )
)

stopifnot(
  "groupOA" %in% colnames(design)
)

stopifnot(
  "MMP13" %in% rownames(dge_normalized)
)

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix tidak full rank.")
}

message("Semua validasi differential expression berhasil.")
# ------------------------------------------------------------
# 4. Estimasi dispersion untuk diagnostic
# ------------------------------------------------------------

dge_dispersion_diagnostic <- edgeR::estimateDisp(
  dge_normalized,
  design = design,
  robust = TRUE
)

cat("\nCommon dispersion diagnostic:\n")
print(
  dge_dispersion_diagnostic$common.dispersion
)

cat("\nBiological coefficient of variation:\n")
print(
  sqrt(
    dge_dispersion_diagnostic$common.dispersion
  )
)
# ------------------------------------------------------------
# 5. Robust quasi-likelihood model
# ------------------------------------------------------------

ql_fit <- edgeR::glmQLFit(
  dge_normalized,
  design = design,
  robust = TRUE
)

cat("\nClass QL fit:\n")
print(class(ql_fit))

cat("\nNama koefisien QL fit:\n")
print(colnames(ql_fit$coefficients))

cat("\nPrior degrees of freedom:\n")
print(summary(ql_fit$df.prior))
# ------------------------------------------------------------
# 6. Quasi-likelihood F-test: OA vs Control
# ------------------------------------------------------------

ql_test <- edgeR::glmQLFTest(
  ql_fit,
  coef = "groupOA"
)

cat("\nKoefisien yang diuji:\n")
print(ql_test$comparison)

de_table <- edgeR::topTags(
  ql_test,
  n = Inf,
  sort.by = "PValue"
)$table %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene")

cat("\nDimensi tabel DE:\n")
print(dim(de_table))

cat("\nSepuluh gen paling signifikan:\n")
print(head(de_table, 10))
# ------------------------------------------------------------
# 7. Rata-rata normalized expression per kelompok
# ------------------------------------------------------------

control_samples <- sample_metadata$sample_id[
  sample_metadata$group == "Control"
]

oa_samples <- sample_metadata$sample_id[
  sample_metadata$group == "OA"
]

mean_control <- rowMeans(
  log_cpm[, control_samples, drop = FALSE]
)

mean_oa <- rowMeans(
  log_cpm[, oa_samples, drop = FALSE]
)

expression_summary <- tibble::tibble(
  gene = rownames(log_cpm),
  
  mean_log2CPM_Control =
    as.numeric(mean_control),
  
  mean_log2CPM_OA =
    as.numeric(mean_oa),
  
  unadjusted_mean_difference =
    mean_log2CPM_OA -
    mean_log2CPM_Control
)

de_table <- de_table %>%
  dplyr::left_join(
    expression_summary,
    by = "gene"
  )
# ------------------------------------------------------------
# 8. Klasifikasi hasil
# ------------------------------------------------------------

de_table <- de_table %>%
  dplyr::mutate(
    statistical_status = dplyr::case_when(
      FDR < 0.05 & logFC > 0 ~ "Upregulated",
      FDR < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    ),
    
    volcano_status = dplyr::case_when(
      FDR < 0.05 & logFC >= 1 ~
        "Upregulated",
      
      FDR < 0.05 & logFC <= -1 ~
        "Downregulated",
      
      TRUE ~
        "Not significant"
    ),
    
    fold_change = 2^logFC,
    
    minus_log10_PValue =
      -log10(
        pmax(
          PValue,
          .Machine$double.xmin
        )
      ),
    
    signed_QL_statistic =
      sign(logFC) * sqrt(F)
  )
# ------------------------------------------------------------
# 9. Ringkasan DEG
# ------------------------------------------------------------

de_summary <- tibble::tibble(
  criterion = c(
    "Total genes tested",
    "FDR < 0.05",
    "Upregulated at FDR < 0.05",
    "Downregulated at FDR < 0.05",
    "FDR < 0.05 and abs(logFC) >= 1",
    "Upregulated with logFC >= 1",
    "Downregulated with logFC <= -1"
  ),
  
  number_of_genes = c(
    nrow(de_table),
    
    sum(
      de_table$FDR < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      de_table$FDR < 0.05 &
        de_table$logFC > 0,
      na.rm = TRUE
    ),
    
    sum(
      de_table$FDR < 0.05 &
        de_table$logFC < 0,
      na.rm = TRUE
    ),
    
    sum(
      de_table$FDR < 0.05 &
        abs(de_table$logFC) >= 1,
      na.rm = TRUE
    ),
    
    sum(
      de_table$FDR < 0.05 &
        de_table$logFC >= 1,
      na.rm = TRUE
    ),
    
    sum(
      de_table$FDR < 0.05 &
        de_table$logFC <= -1,
      na.rm = TRUE
    )
  )
)

print(de_summary)
# ------------------------------------------------------------
# 10. Hasil utama MMP13
# ------------------------------------------------------------

mmp13_edger_result <- de_table %>%
  dplyr::filter(
    gene == "MMP13"
  ) %>%
  dplyr::mutate(
    dataset = "GSE114007",
    
    comparison = "OA vs Control",
    
    adjusted_for =
      "Sequencing platform",
    
    primary_method =
      "edgeR robust quasi-likelihood"
  ) %>%
  dplyr::select(
    dataset,
    comparison,
    adjusted_for,
    primary_method,
    gene,
    logFC,
    fold_change,
    logCPM,
    F,
    PValue,
    FDR,
    mean_log2CPM_Control,
    mean_log2CPM_OA,
    statistical_status
  )

cat("\nHasil edgeR untuk MMP13:\n")
print(mmp13_edger_result)
# ------------------------------------------------------------
# 11. Menentukan gen yang akan diberi label
# ------------------------------------------------------------

top_up <- de_table %>%
  dplyr::filter(
    volcano_status == "Upregulated"
  ) %>%
  dplyr::arrange(PValue) %>%
  dplyr::slice_head(n = 5)

top_down <- de_table %>%
  dplyr::filter(
    volcano_status == "Downregulated"
  ) %>%
  dplyr::arrange(PValue) %>%
  dplyr::slice_head(n = 5)

label_genes <- dplyr::bind_rows(
  top_up,
  top_down,
  de_table %>%
    dplyr::filter(gene == "MMP13")
) %>%
  dplyr::distinct(gene, .keep_all = TRUE)
# ------------------------------------------------------------
# 12. MA plot
# ------------------------------------------------------------

ma_plot <- ggplot(
  de_table,
  aes(
    x = logCPM,
    y = logFC,
    color = volcano_status
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.2
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text(
    data = de_table %>%
      dplyr::filter(gene == "MMP13"),
    aes(label = gene),
    size = 3.2,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    title = "MA plot of GSE114007",
    subtitle =
      "OA versus Control, adjusted for sequencing platform",
    x = expression("Average " * log[2] * "-CPM"),
    y = expression(log[2] * " fold change"),
    color = "Status"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(ma_plot)

ggsave(
  filename =
    "results/figures/GSE114007_MA_plot.pdf",
  plot = ma_plot,
  width = 7,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE114007_MA_plot.tiff",
  plot = ma_plot,
  width = 7,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)
# ------------------------------------------------------------
# 13. Ranked gene list untuk GSEA
# ------------------------------------------------------------

gsea_ranked_table <- de_table %>%
  dplyr::select(
    gene,
    signed_QL_statistic,
    logFC,
    PValue,
    FDR
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      signed_QL_statistic
    )
  )

cat("\nSepuluh ranking GSEA tertinggi:\n")
print(head(gsea_ranked_table, 10))

cat("\nSepuluh ranking GSEA terendah:\n")
print(tail(gsea_ranked_table, 10))
# ------------------------------------------------------------
# 14. Interaction sensitivity analysis
# ------------------------------------------------------------

sample_metadata$group <- factor(
  sample_metadata$group,
  levels = c("Control", "OA")
)

sample_metadata$platform_id <- factor(
  sample_metadata$platform_id,
  levels = c("GPL11154", "GPL18573")
)

design_interaction <- model.matrix(
  ~ platform_id * group,
  data = sample_metadata
)

rownames(design_interaction) <-
  sample_metadata$sample_id

cat("\nNama koefisien interaction model:\n")
print(colnames(design_interaction))

cat("\nRank interaction model:\n")
print(qr(design_interaction)$rank)

cat("\nJumlah kolom interaction model:\n")
print(ncol(design_interaction))

stopifnot(
  qr(design_interaction)$rank ==
    ncol(design_interaction)
)
# ------------------------------------------------------------
# 15. limma-voom effect size untuk meta-analysis
# ------------------------------------------------------------

voom_fit <- edgeR::voomLmFit(
  counts = dge_normalized,
  design = design,
  sample.weights = FALSE,
  plot = FALSE,
  keep.EList = TRUE
)

voom_fit <- limma::eBayes(
  voom_fit,
  robust = TRUE
)

cat("\nKoefisien voom fit:\n")
print(colnames(voom_fit$coefficients))
# ------------------------------------------------------------
# 16. Menyimpan objek hasil analisis
# ------------------------------------------------------------

saveRDS(
  ql_fit,
  file =
    "data_processed/GSE114007_edgeR_QL_fit.rds"
)

saveRDS(
  ql_test,
  file =
    "data_processed/GSE114007_edgeR_QL_test_OA_vs_Control.rds"
)

saveRDS(
  de_table,
  file =
    "data_processed/GSE114007_edgeR_DE_table.rds"
)

saveRDS(
  interaction_test,
  file =
    "data_processed/GSE114007_platform_group_interaction_test.rds"
)

saveRDS(
  voom_fit,
  file =
    "data_processed/GSE114007_voomLmFit_meta_effect_model.rds"
)

saveRDS(
  gsea_ranked_table,
  file =
    "data_processed/GSE114007_GSEA_ranked_genes.rds"
)
# ------------------------------------------------------------
# 17. Tabel signifikan
# ------------------------------------------------------------

significant_fdr <- de_table %>%
  dplyr::filter(
    FDR < 0.05
  )

significant_fdr_lfc1 <- de_table %>%
  dplyr::filter(
    FDR < 0.05,
    abs(logFC) >= 1
  )

top_50_genes <- de_table %>%
  dplyr::slice_head(n = 50)