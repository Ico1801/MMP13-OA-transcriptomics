# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 07_GSE117999_limma_DE.R
# PURPOSE : Differential-expression analysis GSE117999
#           menggunakan limma
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

message("Semua paket differential expression berhasil dimuat.")

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
# 3. Memuat hasil QC
# ------------------------------------------------------------

expression_biological <- readRDS(
  "data_processed/GSE117999_expression_biological_annotated.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE117999_metadata_QC.rds"
)

feature_annotation <- readRDS(
  "data_processed/GSE117999_feature_annotation_original.rds"
)

design_unadjusted <- readRDS(
  "data_processed/GSE117999_design_primary.rds"
)

design_adjusted <- readRDS(
  "data_processed/GSE117999_design_adjusted_sensitivity.rds"
)

mmp13_annotation <- readRDS(
  "data_processed/GSE117999_MMP13_probe_annotation.rds"
)

cat("\nDimensi expression matrix:\n")
print(dim(expression_biological))

cat("\nDimensi unadjusted design:\n")
print(dim(design_unadjusted))

cat("\nDimensi adjusted design:\n")
print(dim(design_adjusted))

cat("\nKoefisien unadjusted design:\n")
print(colnames(design_unadjusted))

cat("\nKoefisien adjusted design:\n")
print(colnames(design_adjusted))

# ------------------------------------------------------------
# 4. Validasi data
# ------------------------------------------------------------

stopifnot(
  ncol(expression_biological) ==
    nrow(sample_metadata)
)

stopifnot(
  identical(
    colnames(expression_biological),
    sample_metadata$expression_column
  )
)

stopifnot(
  identical(
    rownames(design_unadjusted),
    sample_metadata$expression_column
  )
)

stopifnot(
  identical(
    rownames(design_adjusted),
    sample_metadata$expression_column
  )
)

stopifnot(
  "groupOA" %in%
    colnames(design_unadjusted)
)

stopifnot(
  "groupOA" %in%
    colnames(design_adjusted)
)

stopifnot(
  qr(design_unadjusted)$rank ==
    ncol(design_unadjusted)
)

stopifnot(
  qr(design_adjusted)$rank ==
    ncol(design_adjusted)
)

stopifnot(
  !anyNA(expression_biological)
)

stopifnot(
  all(is.finite(expression_biological))
)

stopifnot(
  nrow(mmp13_annotation) == 1
)

mmp13_probe_id <- as.character(
  mmp13_annotation$probe_id
)

stopifnot(
  mmp13_probe_id %in%
    rownames(expression_biological)
)

message("Semua data dan design matrix berhasil divalidasi.")

# ------------------------------------------------------------
# 5. Menyelaraskan anotasi probe
# ------------------------------------------------------------

annotation_index <- match(
  rownames(expression_biological),
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
    rownames(expression_biological),
    as.character(
      annotation_aligned$probe_id
    )
  )
)

probe_annotation <- tibble::tibble(
  probe_id =
    as.character(
      annotation_aligned$probe_id
    ),
  
  agilent_probe_name =
    as.character(
      annotation_aligned$NAME
    ),
  
  gene_symbol =
    trimws(
      as.character(
        annotation_aligned$GENE_SYMBOL
      )
    ),
  
  gene_name =
    as.character(
      annotation_aligned$GENE_NAME
    ),
  
  entrez_id =
    as.character(
      annotation_aligned$LOCUSLINK_ID
    ),
  
  refseq =
    as.character(
      annotation_aligned$REFSEQ
    )
)

stopifnot(
  !anyNA(probe_annotation$gene_symbol)
)

stopifnot(
  all(probe_annotation$gene_symbol != "")
)

# ------------------------------------------------------------
# 6. Model utama: adjusted
# ------------------------------------------------------------

fit_adjusted_raw <- limma::lmFit(
  expression_biological,
  design = design_adjusted
)

fit_adjusted <- limma::eBayes(
  fit_adjusted_raw,
  trend = TRUE,
  robust = TRUE
)

cat("\nKoefisien adjusted fit:\n")
print(colnames(fit_adjusted$coefficients))

cat("\nResidual degrees of freedom adjusted model:\n")
print(summary(fit_adjusted$df.residual))

# ------------------------------------------------------------
# 7. TREAT adjusted: menguji abs(log2FC) > 1
# ------------------------------------------------------------

fit_adjusted_treat <- limma::treat(
  fit_adjusted_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

# ------------------------------------------------------------
# 8. Model sensitivitas: group only
# ------------------------------------------------------------

fit_unadjusted_raw <- limma::lmFit(
  expression_biological,
  design = design_unadjusted
)

fit_unadjusted <- limma::eBayes(
  fit_unadjusted_raw,
  trend = TRUE,
  robust = TRUE
)

cat("\nKoefisien unadjusted fit:\n")
print(colnames(fit_unadjusted$coefficients))

# ------------------------------------------------------------
# 9. Fungsi mengambil tabel limma
# ------------------------------------------------------------

extract_probe_table <- function(
    fitted_model,
    coefficient_name
) {
  
  result <- limma::topTable(
    fitted_model,
    coef = coefficient_name,
    number = Inf,
    sort.by = "P",
    confint = TRUE
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "probe_id"
    ) %>%
    dplyr::rename(
      PValue = P.Value,
      FDR = adj.P.Val,
      CI_95_lower = CI.L,
      CI_95_upper = CI.R
    ) %>%
    dplyr::left_join(
      probe_annotation,
      by = "probe_id"
    )
  
  result
}

probe_de_adjusted <- extract_probe_table(
  fit_adjusted,
  "groupOA"
)

probe_de_unadjusted <- extract_probe_table(
  fit_unadjusted,
  "groupOA"
)

cat("\nDimensi adjusted probe-level table:\n")
print(dim(probe_de_adjusted))

cat("\nDimensi unadjusted probe-level table:\n")
print(dim(probe_de_unadjusted))

# ------------------------------------------------------------
# 10. Rata-rata ekspresi berdasarkan kelompok
# ------------------------------------------------------------

control_samples <- sample_metadata$expression_column[
  sample_metadata$group == "Control"
]

oa_samples <- sample_metadata$expression_column[
  sample_metadata$group == "OA"
]

probe_expression_summary <- tibble::tibble(
  probe_id =
    rownames(expression_biological),
  
  mean_log2_expression_nonOA_APM =
    rowMeans(
      expression_biological[
        ,
        control_samples,
        drop = FALSE
      ]
    ),
  
  mean_log2_expression_OA =
    rowMeans(
      expression_biological[
        ,
        oa_samples,
        drop = FALSE
      ]
    )
) %>%
  dplyr::mutate(
    unadjusted_mean_difference =
      mean_log2_expression_OA -
      mean_log2_expression_nonOA_APM
  )

probe_de_adjusted <- probe_de_adjusted %>%
  dplyr::left_join(
    probe_expression_summary,
    by = "probe_id"
  )

probe_de_unadjusted <- probe_de_unadjusted %>%
  dplyr::left_join(
    probe_expression_summary,
    by = "probe_id"
  )

# ------------------------------------------------------------
# 11. Klasifikasi probe-level adjusted result
# ------------------------------------------------------------

probe_de_adjusted <- probe_de_adjusted %>%
  dplyr::mutate(
    statistical_status =
      dplyr::case_when(
        FDR < 0.05 &
          logFC > 0 ~
          "Upregulated",
        
        FDR < 0.05 &
          logFC < 0 ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      ),
    
    volcano_status =
      dplyr::case_when(
        FDR < 0.05 &
          logFC >= 1 ~
          "Upregulated",
        
        FDR < 0.05 &
          logFC <= -1 ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      ),
    
    fold_change =
      2^logFC,
    
    minus_log10_PValue =
      -log10(
        pmax(
          PValue,
          .Machine$double.xmin
        )
      )
  )

# ------------------------------------------------------------
# 12. Gene-level expression untuk pathway dan visualization
# ------------------------------------------------------------

gene_ids <- probe_annotation$gene_symbol

expression_gene_level <- limma::avereps(
  expression_biological,
  ID = gene_ids
)

cat("\nDimensi gene-level expression matrix:\n")
print(dim(expression_gene_level))

stopifnot(
  "MMP13" %in%
    rownames(expression_gene_level)
)

# ------------------------------------------------------------
# 13. Gene-level adjusted model
# ------------------------------------------------------------

gene_fit_adjusted_raw <- limma::lmFit(
  expression_gene_level,
  design = design_adjusted
)

gene_fit_adjusted <- limma::eBayes(
  gene_fit_adjusted_raw,
  trend = TRUE,
  robust = TRUE
)

gene_fit_adjusted_treat <- limma::treat(
  gene_fit_adjusted_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

gene_de_adjusted <- limma::topTable(
  gene_fit_adjusted,
  coef = "groupOA",
  number = Inf,
  sort.by = "P",
  confint = TRUE
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "gene"
  ) %>%
  dplyr::rename(
    PValue = P.Value,
    FDR = adj.P.Val,
    CI_95_lower = CI.L,
    CI_95_upper = CI.R
  ) %>%
  dplyr::mutate(
    statistical_status =
      dplyr::case_when(
        FDR < 0.05 &
          logFC > 0 ~
          "Upregulated",
        
        FDR < 0.05 &
          logFC < 0 ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      ),
    
    volcano_status =
      dplyr::case_when(
        FDR < 0.05 &
          logFC >= 1 ~
          "Upregulated",
        
        FDR < 0.05 &
          logFC <= -1 ~
          "Downregulated",
        
        TRUE ~
          "Not significant"
      ),
    
    fold_change =
      2^logFC,
    
    minus_log10_PValue =
      -log10(
        pmax(
          PValue,
          .Machine$double.xmin
        )
      )
  )

gene_de_treat <- limma::topTreat(
  gene_fit_adjusted_treat,
  coef = "groupOA",
  number = Inf,
  sort.by = "P"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "gene"
  ) %>%
  dplyr::rename(
    PValue = P.Value,
    FDR = adj.P.Val
  )

# ------------------------------------------------------------
# 14. Gene-level unadjusted model
# ------------------------------------------------------------

gene_fit_unadjusted_raw <- limma::lmFit(
  expression_gene_level,
  design = design_unadjusted
)

gene_fit_unadjusted <- limma::eBayes(
  gene_fit_unadjusted_raw,
  trend = TRUE,
  robust = TRUE
)

gene_de_unadjusted <- limma::topTable(
  gene_fit_unadjusted,
  coef = "groupOA",
  number = Inf,
  sort.by = "P",
  confint = TRUE
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "gene"
  ) %>%
  dplyr::rename(
    PValue = P.Value,
    FDR = adj.P.Val,
    CI_95_lower = CI.L,
    CI_95_upper = CI.R
  )

# ------------------------------------------------------------
# 15. Ringkasan DEG adjusted
# ------------------------------------------------------------

de_summary <- tibble::tibble(
  criterion = c(
    "Total genes tested",
    "FDR < 0.05",
    "Upregulated at FDR < 0.05",
    "Downregulated at FDR < 0.05",
    "FDR < 0.05 and abs(logFC) >= 1",
    "Upregulated with logFC >= 1",
    "Downregulated with logFC <= -1",
    "TREAT FDR < 0.05 for abs(logFC) > 1"
  ),
  
  number_of_genes = c(
    nrow(gene_de_adjusted),
    
    sum(
      gene_de_adjusted$FDR < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      gene_de_adjusted$FDR < 0.05 &
        gene_de_adjusted$logFC > 0,
      na.rm = TRUE
    ),
    
    sum(
      gene_de_adjusted$FDR < 0.05 &
        gene_de_adjusted$logFC < 0,
      na.rm = TRUE
    ),
    
    sum(
      gene_de_adjusted$FDR < 0.05 &
        abs(
          gene_de_adjusted$logFC
        ) >= 1,
      na.rm = TRUE
    ),
    
    sum(
      gene_de_adjusted$FDR < 0.05 &
        gene_de_adjusted$logFC >= 1,
      na.rm = TRUE
    ),
    
    sum(
      gene_de_adjusted$FDR < 0.05 &
        gene_de_adjusted$logFC <= -1,
      na.rm = TRUE
    ),
    
    sum(
      gene_de_treat$FDR < 0.05,
      na.rm = TRUE
    )
  )
)

cat("\nRingkasan differential expression adjusted model:\n")
print(de_summary)

# ------------------------------------------------------------
# 16. Fungsi ekstraksi effect size dan standard error
# ------------------------------------------------------------

extract_limma_effect <- function(
    fitted_model,
    row_id,
    coefficient_name,
    model_name
) {
  
  row_index <- match(
    row_id,
    rownames(
      fitted_model$coefficients
    )
  )
  
  coefficient_index <- match(
    coefficient_name,
    colnames(
      fitted_model$coefficients
    )
  )
  
  if (is.na(row_index)) {
    stop(
      "Feature tidak ditemukan: ",
      row_id
    )
  }
  
  if (is.na(coefficient_index)) {
    stop(
      "Koefisien tidak ditemukan: ",
      coefficient_name
    )
  }
  
  effect <- as.numeric(
    fitted_model$coefficients[
      row_index,
      coefficient_index
    ]
  )
  
  standard_error <- as.numeric(
    fitted_model$stdev.unscaled[
      row_index,
      coefficient_index
    ] *
      sqrt(
        fitted_model$s2.post[
          row_index
        ]
      )
  )
  
  df_total <- if (
    length(
      fitted_model$df.total
    ) == 1
  ) {
    as.numeric(
      fitted_model$df.total
    )
  } else {
    as.numeric(
      fitted_model$df.total[
        row_index
      ]
    )
  }
  
  critical_value <- stats::qt(
    0.975,
    df = df_total
  )
  
  full_table <- limma::topTable(
    fitted_model,
    coef = coefficient_name,
    number = Inf,
    sort.by = "none"
  )
  
  current_statistics <- full_table[
    row_id,
    ,
    drop = FALSE
  ]
  
  tibble::tibble(
    model = model_name,
    
    feature_id = row_id,
    
    log2FC = effect,
    
    standard_error =
      standard_error,
    
    CI_95_lower =
      effect -
      critical_value *
      standard_error,
    
    CI_95_upper =
      effect +
      critical_value *
      standard_error,
    
    degrees_of_freedom =
      df_total,
    
    moderated_t =
      as.numeric(
        current_statistics$t
      ),
    
    PValue =
      as.numeric(
        current_statistics$P.Value
      ),
    
    FDR =
      as.numeric(
        current_statistics$adj.P.Val
      )
  )
}

# ------------------------------------------------------------
# 17. MMP13 adjusted dan unadjusted
# ------------------------------------------------------------

mmp13_adjusted <- extract_limma_effect(
  fitted_model =
    fit_adjusted,
  
  row_id =
    mmp13_probe_id,
  
  coefficient_name =
    "groupOA",
  
  model_name =
    "Adjusted for age, BMI, and sex"
)

mmp13_unadjusted <- extract_limma_effect(
  fitted_model =
    fit_unadjusted,
  
  row_id =
    mmp13_probe_id,
  
  coefficient_name =
    "groupOA",
  
  model_name =
    "Unadjusted group comparison"
)

mmp13_model_comparison <- dplyr::bind_rows(
  mmp13_unadjusted,
  mmp13_adjusted
) %>%
  dplyr::mutate(
    gene = "MMP13",
    
    probe_id =
      mmp13_probe_id,
    
    agilent_probe =
      as.character(
        mmp13_annotation$NAME
      ),
    
    n_nonOA_APM =
      sum(
        sample_metadata$group ==
          "Control"
      ),
    
    n_OA =
      sum(
        sample_metadata$group ==
          "OA"
      ),
    
    mean_nonOA_APM =
      mean(
        expression_biological[
          mmp13_probe_id,
          control_samples
        ]
      ),
    
    mean_OA =
      mean(
        expression_biological[
          mmp13_probe_id,
          oa_samples
        ]
      )
  )

cat("\nPerbandingan hasil MMP13:\n")

print(
  mmp13_model_comparison,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 18. Forest-input untuk GSE117999
# ------------------------------------------------------------

mmp13_forest_input <- mmp13_adjusted %>%
  dplyr::transmute(
    dataset =
      "GSE117999",
    
    tissue =
      "Human knee cartilage",
    
    comparison =
      "End-stage OA vs non-OA APM",
    
    n_control =
      sum(
        sample_metadata$group ==
          "Control"
      ),
    
    n_OA =
      sum(
        sample_metadata$group ==
          "OA"
      ),
    
    effect_measure =
      "Adjusted log2 expression difference",
    
    effect_method =
      "limma robust empirical Bayes",
    
    covariate_adjustment =
      "Age, BMI, and sex",
    
    log2FC =
      log2FC,
    
    standard_error =
      standard_error,
    
    CI_95_lower =
      CI_95_lower,
    
    CI_95_upper =
      CI_95_upper,
    
    degrees_of_freedom =
      degrees_of_freedom,
    
    PValue =
      PValue,
    
    FDR =
      FDR
  )

cat("\nMMP13 forest input:\n")

print(
  mmp13_forest_input,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 19. MMP13 effect comparison plot
# ------------------------------------------------------------

mmp13_model_comparison$model <- factor(
  mmp13_model_comparison$model,
  levels = c(
    "Unadjusted group comparison",
    "Adjusted for age, BMI, and sex"
  )
)

mmp13_effect_plot <- ggplot(
  mmp13_model_comparison,
  aes(
    x = model,
    y = log2FC
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = CI_95_lower,
      ymax = CI_95_upper
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3.2
  ) +
  coord_flip() +
  labs(
    title =
      "MMP13 effect estimates in GSE117999",
    
    subtitle =
      "OA versus non-OA APM cartilage",
    
    x = NULL,
    
    y =
      expression(
        log[2] *
          " expression difference"
      )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(mmp13_effect_plot)

ggsave(
  filename =
    "results/figures/GSE117999_MMP13_adjusted_vs_unadjusted.pdf",
  
  plot =
    mmp13_effect_plot,
  
  width = 7,
  height = 4
)

ggsave(
  filename =
    "results/figures/GSE117999_MMP13_adjusted_vs_unadjusted.tiff",
  
  plot =
    mmp13_effect_plot,
  
  width = 7,
  height = 4,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 20. Volcano plot adjusted gene-level results
# ------------------------------------------------------------

top_up <- gene_de_adjusted %>%
  dplyr::filter(
    volcano_status ==
      "Upregulated"
  ) %>%
  dplyr::slice_min(
    order_by = PValue,
    n = 5,
    with_ties = FALSE
  )

top_down <- gene_de_adjusted %>%
  dplyr::filter(
    volcano_status ==
      "Downregulated"
  ) %>%
  dplyr::slice_min(
    order_by = PValue,
    n = 5,
    with_ties = FALSE
  )

label_genes <- dplyr::bind_rows(
  top_up,
  top_down,
  gene_de_adjusted %>%
    dplyr::filter(
      gene == "MMP13"
    )
) %>%
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )

volcano_plot <- ggplot(
  gene_de_adjusted,
  aes(
    x = logFC,
    y = minus_log10_PValue,
    color = volcano_status
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.2
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept =
      -log10(0.05),
    
    linetype = "dotted",
    linewidth = 0.45
  ) +
  geom_text(
    data = label_genes,
    aes(label = gene),
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "Differential expression in GSE117999",
    
    subtitle =
      "OA versus non-OA APM, adjusted for age, BMI, and sex",
    
    x =
      expression(
        log[2] *
          " fold change"
      ),
    
    y =
      expression(
        -log[10] *
          " P-value"
      ),
    
    color =
      "Status"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(volcano_plot)

ggsave(
  filename =
    "results/figures/GSE117999_volcano_adjusted.pdf",
  
  plot =
    volcano_plot,
  
  width = 7,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE117999_volcano_adjusted.tiff",
  
  plot =
    volcano_plot,
  
  width = 7,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 21. MA plot adjusted
# ------------------------------------------------------------

ma_plot <- ggplot(
  gene_de_adjusted,
  aes(
    x = AveExpr,
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
    data =
      gene_de_adjusted %>%
      dplyr::filter(
        gene == "MMP13"
      ),
    
    aes(label = gene),
    size = 3.2,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    title =
      "MA plot of GSE117999",
    
    subtitle =
      "Adjusted for age, BMI, and sex",
    
    x =
      "Average processed expression",
    
    y =
      expression(
        log[2] *
          " fold change"
      ),
    
    color =
      "Status"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

print(ma_plot)

ggsave(
  filename =
    "results/figures/GSE117999_MA_adjusted.pdf",
  
  plot =
    ma_plot,
  
  width = 7,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE117999_MA_adjusted.tiff",
  
  plot =
    ma_plot,
  
  width = 7,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 22. Mean-variance trend diagnostic
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE117999_mean_variance_trend.pdf",
  
  width = 7,
  height = 5.5
)

limma::plotSA(
  fit_adjusted,
  main =
    "Mean-variance trend in GSE117999"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE117999_mean_variance_trend.tiff",
  
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

limma::plotSA(
  fit_adjusted,
  main =
    "Mean-variance trend in GSE117999"
)

dev.off()

# ------------------------------------------------------------
# 23. Ranked gene list untuk pathway analysis
# ------------------------------------------------------------

gsea_ranked_table <- gene_de_adjusted %>%
  dplyr::select(
    gene,
    moderated_t = t,
    logFC,
    PValue,
    FDR
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      moderated_t
    )
  )

cat("\nSepuluh gene ranking tertinggi:\n")
print(head(gsea_ranked_table, 10))

cat("\nSepuluh gene ranking terendah:\n")
print(tail(gsea_ranked_table, 10))

# ------------------------------------------------------------
# 24. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  fit_adjusted,
  file =
    "data_processed/GSE117999_limma_fit_adjusted.rds"
)

saveRDS(
  fit_unadjusted,
  file =
    "data_processed/GSE117999_limma_fit_unadjusted.rds"
)

saveRDS(
  probe_de_adjusted,
  file =
    "data_processed/GSE117999_probe_DE_adjusted.rds"
)

saveRDS(
  probe_de_unadjusted,
  file =
    "data_processed/GSE117999_probe_DE_unadjusted.rds"
)

saveRDS(
  gene_de_adjusted,
  file =
    "data_processed/GSE117999_gene_DE_adjusted.rds"
)

saveRDS(
  gene_de_unadjusted,
  file =
    "data_processed/GSE117999_gene_DE_unadjusted.rds"
)

saveRDS(
  gene_de_treat,
  file =
    "data_processed/GSE117999_gene_TREAT_lfc1.rds"
)

saveRDS(
  gsea_ranked_table,
  file =
    "data_processed/GSE117999_GSEA_ranked_genes.rds"
)

saveRDS(
  mmp13_forest_input,
  file =
    "data_processed/GSE117999_MMP13_forest_input.rds"
)

# ------------------------------------------------------------
# 25. Simpan full probe-level CSV
# ------------------------------------------------------------

utils::write.csv(
  probe_de_adjusted,
  file =
    "results/tables/GSE117999_probe_DE_adjusted.csv",
  row.names = FALSE
)

utils::write.csv(
  probe_de_unadjusted,
  file =
    "results/tables/GSE117999_probe_DE_unadjusted.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 26. Menyiapkan Excel workbook
# ------------------------------------------------------------

significant_adjusted <- gene_de_adjusted %>%
  dplyr::filter(
    FDR < 0.05
  )

significant_adjusted_lfc1 <- gene_de_adjusted %>%
  dplyr::filter(
    FDR < 0.05,
    abs(logFC) >= 1
  )

top_100_probe_adjusted <- probe_de_adjusted %>%
  dplyr::slice_head(n = 100)

openxlsx::write.xlsx(
  list(
    DE_Summary =
      as.data.frame(
        de_summary
      ),
    
    MMP13_Model_Comparison =
      as.data.frame(
        mmp13_model_comparison
      ),
    
    MMP13_Forest_Input =
      as.data.frame(
        mmp13_forest_input
      ),
    
    Gene_DE_Adjusted =
      as.data.frame(
        gene_de_adjusted
      ),
    
    Gene_DE_Unadjusted =
      as.data.frame(
        gene_de_unadjusted
      ),
    
    Gene_TREAT_LFC1 =
      as.data.frame(
        gene_de_treat
      ),
    
    Significant_Adjusted =
      as.data.frame(
        significant_adjusted
      ),
    
    Significant_Adjusted_LFC1 =
      as.data.frame(
        significant_adjusted_lfc1
      ),
    
    GSEA_Ranked_Genes =
      as.data.frame(
        gsea_ranked_table
      ),
    
    Top_100_Probe_Adjusted =
      as.data.frame(
        top_100_probe_adjusted
      )
  ),
  
  file =
    "results/tables/GSE117999_limma_DE_results.xlsx",
  
  overwrite = TRUE
)

openxlsx::write.xlsx(
  as.data.frame(
    mmp13_forest_input
  ),
  
  file =
    "results/tables/GSE117999_MMP13_forest_input.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 27. Session information
# ------------------------------------------------------------

sink(
  "results/tables/GSE117999_DE_sessionInfo.txt"
)

print(sessionInfo())

sink()

# ------------------------------------------------------------
# 28. Pesan akhir
# ------------------------------------------------------------

message("")
message("================================================")
message("GSE117999 DIFFERENTIAL EXPRESSION SELESAI")
message("================================================")
message(
  "Jumlah gen diuji        : ",
  nrow(gene_de_adjusted)
)
message(
  "Adjusted FDR < 0.05     : ",
  sum(
    gene_de_adjusted$FDR < 0.05
  )
)
message(
  "Adjusted FDR & |LFC|>=1 : ",
  sum(
    gene_de_adjusted$FDR < 0.05 &
      abs(
        gene_de_adjusted$logFC
      ) >= 1
  )
)
message(
  "MMP13 unadjusted logFC  : ",
  round(
    mmp13_unadjusted$log2FC,
    4
  )
)
message(
  "MMP13 adjusted logFC    : ",
  round(
    mmp13_adjusted$log2FC,
    4
  )
)
message(
  "MMP13 adjusted SE       : ",
  round(
    mmp13_adjusted$standard_error,
    4
  )
)
message(
  "MMP13 adjusted 95% CI   : ",
  round(
    mmp13_adjusted$CI_95_lower,
    4
  ),
  " to ",
  round(
    mmp13_adjusted$CI_95_upper,
    4
  )
)
message(
  "MMP13 adjusted FDR      : ",
  format(
    mmp13_adjusted$FDR,
    scientific = TRUE,
    digits = 4
  )
)
message("================================================")