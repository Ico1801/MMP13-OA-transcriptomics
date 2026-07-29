# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 10_GSE57218_limma_paired_DE.R
# PURPOSE : Paired differential-expression analysis GSE57218
# DESIGN  : OA-affected versus preserved cartilage from the
#           same patient
#
# PRIMARY ANALYSIS:
#   33 complete patient pairs
#
# SENSITIVITY ANALYSIS:
#   32 pairs after excluding the QC-flagged RAAK_14 pair
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
    paste(missing_packages, collapse = ", ")
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
  "Semua paket differential expression GSE57218 berhasil dimuat."
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
# 3. Memuat hasil QC
# ------------------------------------------------------------

expression_biological <- readRDS(
  "data_processed/GSE57218_expression_annotated_biological.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE57218_metadata_paired_QC.rds"
)

feature_annotation <- readRDS(
  "data_processed/GSE57218_feature_annotation_original.rds"
)

design_primary <- readRDS(
  "data_processed/GSE57218_design_paired.rds"
)

mmp13_annotation <- readRDS(
  "data_processed/GSE57218_MMP13_probe_annotation.rds"
)

pair_diagnostics <- readRDS(
  "data_processed/GSE57218_pair_diagnostics.rds"
)

mmp13_pair_table <- readRDS(
  "data_processed/GSE57218_MMP13_pair_table.rds"
)

cat("\nDimensi expression matrix biologis:\n")
print(dim(expression_biological))

cat("\nDimensi metadata:\n")
print(dim(sample_metadata))

cat("\nDimensi primary design:\n")
print(dim(design_primary))

cat("\nKoefisien primary design:\n")
print(colnames(design_primary))

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
    rownames(design_primary),
    sample_metadata$expression_column
  )
)

stopifnot(
  nrow(design_primary) ==
    ncol(expression_biological)
)

stopifnot(
  qr(design_primary)$rank ==
    ncol(design_primary)
)

stopifnot(
  "conditionOA" %in%
    colnames(design_primary)
)

stopifnot(
  !anyNA(expression_biological)
)

stopifnot(
  all(
    is.finite(expression_biological)
  )
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

message(
  "Expression matrix, metadata, dan paired design berhasil divalidasi."
)

# ------------------------------------------------------------
# 5. Memastikan factor reference
# ------------------------------------------------------------

sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    
    pair_id = factor(
      as.character(pair_id),
      levels = unique(
        as.character(pair_id)
      )
    ),
    
    condition = factor(
      as.character(condition),
      levels = c(
        "Preserved",
        "OA"
      )
    )
  )

stopifnot(
  levels(sample_metadata$condition)[1] ==
    "Preserved"
)

cat("\nDistribusi kondisi:\n")
print(table(sample_metadata$condition))

cat("\nJumlah pasangan:\n")
print(
  dplyr::n_distinct(
    sample_metadata$pair_id
  )
)

# ------------------------------------------------------------
# 6. Menyelaraskan feature annotation
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
  
  gene_symbol =
    trimws(
      as.character(
        annotation_aligned$Symbol
      )
    ),
  
  illumina_gene =
    as.character(
      annotation_aligned$ILMN_Gene
    ),
  
  entrez_id =
    as.character(
      annotation_aligned$Entrez_Gene_ID
    ),
  
  refseq_id =
    as.character(
      annotation_aligned$RefSeq_ID
    ),
  
  accession =
    as.character(
      annotation_aligned$Accession
    ),
  
  gene_definition =
    as.character(
      annotation_aligned$Definition
    )
)

stopifnot(
  !anyNA(
    probe_annotation$gene_symbol
  )
)

stopifnot(
  all(
    probe_annotation$gene_symbol != ""
  )
)

# ------------------------------------------------------------
# 7. Fungsi ekstraksi topTable
# ------------------------------------------------------------

extract_limma_table <- function(
    fitted_model,
    coefficient_name
) {
  
  result_table <- limma::topTable(
    fitted_model,
    coef = coefficient_name,
    number = Inf,
    sort.by = "P",
    confint = TRUE
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(
      "feature_id"
    ) %>%
    dplyr::rename(
      PValue = P.Value,
      FDR = adj.P.Val,
      CI_95_lower = CI.L,
      CI_95_upper = CI.R
    )
  
  result_table
}

# ------------------------------------------------------------
# 8. PRIMARY probe-level paired model: 33 pairs
# ------------------------------------------------------------

probe_fit_primary_raw <- limma::lmFit(
  expression_biological,
  design = design_primary
)

probe_fit_primary <- limma::eBayes(
  probe_fit_primary_raw,
  trend = TRUE,
  robust = TRUE
)

cat("\nPrimary probe fit coefficients:\n")
print(
  colnames(
    probe_fit_primary$coefficients
  )
)

cat("\nPrimary residual degrees of freedom:\n")
print(
  summary(
    probe_fit_primary$df.residual
  )
)

# TREAT: menguji apakah absolute effect > 1 VST-like unit
probe_fit_primary_treat <- limma::treat(
  probe_fit_primary_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

probe_de_primary <- extract_limma_table(
  fitted_model =
    probe_fit_primary,
  
  coefficient_name =
    "conditionOA"
) %>%
  dplyr::rename(
    probe_id = feature_id
  ) %>%
  dplyr::left_join(
    probe_annotation,
    by = "probe_id"
  )

cat("\nDimensi primary probe DE table:\n")
print(dim(probe_de_primary))

# ------------------------------------------------------------
# 9. Rata-rata ekspresi probe berdasarkan kondisi
# ------------------------------------------------------------

preserved_samples <-
  sample_metadata$expression_column[
    sample_metadata$condition ==
      "Preserved"
  ]

oa_samples <-
  sample_metadata$expression_column[
    sample_metadata$condition ==
      "OA"
  ]

probe_expression_summary <- tibble::tibble(
  
  probe_id =
    rownames(
      expression_biological
    ),
  
  mean_preserved =
    rowMeans(
      expression_biological[
        ,
        preserved_samples,
        drop = FALSE
      ]
    ),
  
  mean_OA =
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
      mean_OA -
      mean_preserved
  )

probe_de_primary <- probe_de_primary %>%
  dplyr::left_join(
    probe_expression_summary,
    by = "probe_id"
  ) %>%
  dplyr::mutate(
    
    statistical_status =
      dplyr::case_when(
        
        FDR < 0.05 &
          logFC > 0 ~
          "Higher in OA",
        
        FDR < 0.05 &
          logFC < 0 ~
          "Lower in OA",
        
        TRUE ~
          "Not significant"
      ),
    
    volcano_status =
      dplyr::case_when(
        
        FDR < 0.05 &
          logFC >= 1 ~
          "Higher in OA",
        
        FDR < 0.05 &
          logFC <= -1 ~
          "Lower in OA",
        
        TRUE ~
          "Not significant"
      ),
    
    minus_log10_PValue =
      -log10(
        pmax(
          PValue,
          .Machine$double.xmin
        )
      )
  )

# ------------------------------------------------------------
# 10. Sensitivity dataset tanpa seluruh pair RAAK_14
# ------------------------------------------------------------

excluded_sensitivity_pair <-
  "RAAK_14"

if (
  !excluded_sensitivity_pair %in%
  as.character(
    sample_metadata$pair_id
  )
) {
  stop(
    "RAAK_14 tidak ditemukan pada metadata."
  )
}

keep_sensitivity_samples <-
  as.character(
    sample_metadata$pair_id
  ) != excluded_sensitivity_pair

sample_metadata_sensitivity <-
  sample_metadata[
    keep_sensitivity_samples,
    ,
    drop = FALSE
  ] %>%
  droplevels()

expression_sensitivity <-
  expression_biological[
    ,
    sample_metadata_sensitivity$expression_column,
    drop = FALSE
  ]

sample_metadata_sensitivity$pair_id <-
  factor(
    as.character(
      sample_metadata_sensitivity$pair_id
    ),
    levels = unique(
      as.character(
        sample_metadata_sensitivity$pair_id
      )
    )
  )

sample_metadata_sensitivity$condition <-
  factor(
    as.character(
      sample_metadata_sensitivity$condition
    ),
    levels = c(
      "Preserved",
      "OA"
    )
  )

design_sensitivity <- model.matrix(
  ~ pair_id + condition,
  data = sample_metadata_sensitivity
)

rownames(design_sensitivity) <-
  sample_metadata_sensitivity$expression_column

cat("\nSensitivity sample number:\n")
print(nrow(sample_metadata_sensitivity))

cat("\nSensitivity pair number:\n")
print(
  dplyr::n_distinct(
    sample_metadata_sensitivity$pair_id
  )
)

cat("\nSensitivity design dimensions:\n")
print(dim(design_sensitivity))

cat("\nSensitivity design rank:\n")
print(qr(design_sensitivity)$rank)

cat("\nSensitivity residual df:\n")
print(
  nrow(design_sensitivity) -
    qr(design_sensitivity)$rank
)

stopifnot(
  ncol(expression_sensitivity) == 64
)

stopifnot(
  dplyr::n_distinct(
    sample_metadata_sensitivity$pair_id
  ) == 32
)

stopifnot(
  qr(design_sensitivity)$rank ==
    ncol(design_sensitivity)
)

stopifnot(
  "conditionOA" %in%
    colnames(design_sensitivity)
)

# ------------------------------------------------------------
# 11. Sensitivity probe-level model
# ------------------------------------------------------------

probe_fit_sensitivity_raw <- limma::lmFit(
  expression_sensitivity,
  design = design_sensitivity
)

probe_fit_sensitivity <- limma::eBayes(
  probe_fit_sensitivity_raw,
  trend = TRUE,
  robust = TRUE
)

probe_fit_sensitivity_treat <- limma::treat(
  probe_fit_sensitivity_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

probe_de_sensitivity <- extract_limma_table(
  fitted_model =
    probe_fit_sensitivity,
  
  coefficient_name =
    "conditionOA"
) %>%
  dplyr::rename(
    probe_id = feature_id
  ) %>%
  dplyr::left_join(
    probe_annotation,
    by = "probe_id"
  ) %>%
  dplyr::mutate(
    
    statistical_status =
      dplyr::case_when(
        
        FDR < 0.05 &
          logFC > 0 ~
          "Higher in OA",
        
        FDR < 0.05 &
          logFC < 0 ~
          "Lower in OA",
        
        TRUE ~
          "Not significant"
      )
  )

# ------------------------------------------------------------
# 12. Membuat clean gene symbols
# ------------------------------------------------------------

gene_symbol_clean <-
  trimws(
    probe_annotation$gene_symbol
  )

# Probe dengan anotasi ambigu tidak digunakan saat collapse
# ke tingkat gene.
ambiguous_symbol <- stringr::str_detect(
  gene_symbol_clean,
  "///|//|;|,|\\||\\s"
)

keep_unambiguous_gene_symbol <- (
  !is.na(gene_symbol_clean) &
    gene_symbol_clean != "" &
    !ambiguous_symbol
)

cat("\nProbe dengan unambiguous gene symbol:\n")
print(
  sum(
    keep_unambiguous_gene_symbol
  )
)

cat("\nProbe dengan ambiguous gene symbol:\n")
print(
  sum(
    ambiguous_symbol
  )
)

# ------------------------------------------------------------
# 13. Collapse probe menjadi gene-level expression
# ------------------------------------------------------------

expression_gene_primary <- limma::avereps(
  expression_biological[
    keep_unambiguous_gene_symbol,
    ,
    drop = FALSE
  ],
  ID =
    gene_symbol_clean[
      keep_unambiguous_gene_symbol
    ]
)

cat("\nDimensi gene-level expression:\n")
print(dim(expression_gene_primary))

stopifnot(
  "MMP13" %in%
    rownames(expression_gene_primary)
)

expression_gene_sensitivity <-
  expression_gene_primary[
    ,
    sample_metadata_sensitivity$expression_column,
    drop = FALSE
  ]

# ------------------------------------------------------------
# 14. Primary gene-level model
# ------------------------------------------------------------

gene_fit_primary_raw <- limma::lmFit(
  expression_gene_primary,
  design = design_primary
)

gene_fit_primary <- limma::eBayes(
  gene_fit_primary_raw,
  trend = TRUE,
  robust = TRUE
)

gene_fit_primary_treat <- limma::treat(
  gene_fit_primary_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

gene_de_primary <- extract_limma_table(
  fitted_model =
    gene_fit_primary,
  
  coefficient_name =
    "conditionOA"
) %>%
  dplyr::rename(
    gene = feature_id
  ) %>%
  dplyr::mutate(
    
    statistical_status =
      dplyr::case_when(
        
        FDR < 0.05 &
          logFC > 0 ~
          "Higher in OA",
        
        FDR < 0.05 &
          logFC < 0 ~
          "Lower in OA",
        
        TRUE ~
          "Not significant"
      ),
    
    volcano_status =
      dplyr::case_when(
        
        FDR < 0.05 &
          logFC >= 1 ~
          "Higher in OA",
        
        FDR < 0.05 &
          logFC <= -1 ~
          "Lower in OA",
        
        TRUE ~
          "Not significant"
      ),
    
    minus_log10_PValue =
      -log10(
        pmax(
          PValue,
          .Machine$double.xmin
        )
      )
  )

gene_de_primary_treat <- limma::topTreat(
  gene_fit_primary_treat,
  coef = "conditionOA",
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
# 15. Sensitivity gene-level model
# ------------------------------------------------------------

gene_fit_sensitivity_raw <- limma::lmFit(
  expression_gene_sensitivity,
  design = design_sensitivity
)

gene_fit_sensitivity <- limma::eBayes(
  gene_fit_sensitivity_raw,
  trend = TRUE,
  robust = TRUE
)

gene_fit_sensitivity_treat <- limma::treat(
  gene_fit_sensitivity_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

gene_de_sensitivity <- extract_limma_table(
  fitted_model =
    gene_fit_sensitivity,
  
  coefficient_name =
    "conditionOA"
) %>%
  dplyr::rename(
    gene = feature_id
  ) %>%
  dplyr::mutate(
    
    statistical_status =
      dplyr::case_when(
        
        FDR < 0.05 &
          logFC > 0 ~
          "Higher in OA",
        
        FDR < 0.05 &
          logFC < 0 ~
          "Lower in OA",
        
        TRUE ~
          "Not significant"
      )
  )

gene_de_sensitivity_treat <- limma::topTreat(
  gene_fit_sensitivity_treat,
  coef = "conditionOA",
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
# 16. Ringkasan differential expression
# ------------------------------------------------------------

summarize_de_result <- function(
    de_table,
    treat_table,
    model_name
) {
  
  tibble::tibble(
    
    model =
      model_name,
    
    total_genes_tested =
      nrow(de_table),
    
    FDR_less_than_0_05 =
      sum(
        de_table$FDR < 0.05,
        na.rm = TRUE
      ),
    
    higher_in_OA_FDR =
      sum(
        de_table$FDR < 0.05 &
          de_table$logFC > 0,
        na.rm = TRUE
      ),
    
    lower_in_OA_FDR =
      sum(
        de_table$FDR < 0.05 &
          de_table$logFC < 0,
        na.rm = TRUE
      ),
    
    FDR_and_abs_effect_at_least_1 =
      sum(
        de_table$FDR < 0.05 &
          abs(de_table$logFC) >= 1,
        na.rm = TRUE
      ),
    
    higher_in_OA_effect_at_least_1 =
      sum(
        de_table$FDR < 0.05 &
          de_table$logFC >= 1,
        na.rm = TRUE
      ),
    
    lower_in_OA_effect_at_least_1 =
      sum(
        de_table$FDR < 0.05 &
          de_table$logFC <= -1,
        na.rm = TRUE
      ),
    
    TREAT_FDR_less_than_0_05 =
      sum(
        treat_table$FDR < 0.05,
        na.rm = TRUE
      )
  )
}

de_summary <- dplyr::bind_rows(
  
  summarize_de_result(
    de_table =
      gene_de_primary,
    
    treat_table =
      gene_de_primary_treat,
    
    model_name =
      "Primary paired model: 33 pairs"
  ),
  
  summarize_de_result(
    de_table =
      gene_de_sensitivity,
    
    treat_table =
      gene_de_sensitivity_treat,
    
    model_name =
      "Sensitivity paired model: 32 pairs, excluding RAAK_14"
  )
)

cat("\nDifferential-expression summary:\n")

print(
  de_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 17. Fungsi ekstraksi effect size dan moderated SE
# ------------------------------------------------------------

extract_limma_effect <- function(
    fitted_model,
    feature_id,
    coefficient_name,
    analysis_name,
    number_of_pairs
) {
  
  feature_index <- match(
    feature_id,
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
  
  if (is.na(feature_index)) {
    stop(
      "Feature tidak ditemukan: ",
      feature_id
    )
  }
  
  if (is.na(coefficient_index)) {
    stop(
      "Koefisien tidak ditemukan: ",
      coefficient_name
    )
  }
  
  effect_estimate <- as.numeric(
    fitted_model$coefficients[
      feature_index,
      coefficient_index
    ]
  )
  
  moderated_standard_error <- as.numeric(
    fitted_model$stdev.unscaled[
      feature_index,
      coefficient_index
    ] *
      sqrt(
        fitted_model$s2.post[
          feature_index
        ]
      )
  )
  
  df_total_vector <-
    fitted_model$df.total
  
  degrees_of_freedom <- if (
    length(df_total_vector) == 1
  ) {
    
    as.numeric(
      df_total_vector
    )
    
  } else {
    
    as.numeric(
      df_total_vector[
        feature_index
      ]
    )
  }
  
  critical_value <- stats::qt(
    p = 0.975,
    df = degrees_of_freedom
  )
  
  statistics_table <- limma::topTable(
    fitted_model,
    coef = coefficient_name,
    number = Inf,
    sort.by = "none"
  )
  
  feature_statistics <- statistics_table[
    feature_id,
    ,
    drop = FALSE
  ]
  
  tibble::tibble(
    
    analysis =
      analysis_name,
    
    number_of_pairs =
      number_of_pairs,
    
    feature_id =
      feature_id,
    
    effect_scale =
      "Processed VST-like expression difference",
    
    effect_direction =
      "OA-affected minus preserved",
    
    effect =
      effect_estimate,
    
    standard_error =
      moderated_standard_error,
    
    CI_95_lower =
      effect_estimate -
      critical_value *
      moderated_standard_error,
    
    CI_95_upper =
      effect_estimate +
      critical_value *
      moderated_standard_error,
    
    degrees_of_freedom =
      degrees_of_freedom,
    
    moderated_t =
      as.numeric(
        feature_statistics$t
      ),
    
    PValue =
      as.numeric(
        feature_statistics$P.Value
      ),
    
    FDR =
      as.numeric(
        feature_statistics$adj.P.Val
      )
  )
}

# ------------------------------------------------------------
# 18. MMP13 primary dan sensitivity estimates
# ------------------------------------------------------------

mmp13_primary_effect <- extract_limma_effect(
  
  fitted_model =
    probe_fit_primary,
  
  feature_id =
    mmp13_probe_id,
  
  coefficient_name =
    "conditionOA",
  
  analysis_name =
    "Primary: all 33 complete pairs",
  
  number_of_pairs =
    33
)

mmp13_sensitivity_effect <- extract_limma_effect(
  
  fitted_model =
    probe_fit_sensitivity,
  
  feature_id =
    mmp13_probe_id,
  
  coefficient_name =
    "conditionOA",
  
  analysis_name =
    "Sensitivity: 32 pairs excluding RAAK_14",
  
  number_of_pairs =
    32
)

# ------------------------------------------------------------
# 19. Direct paired descriptive estimates
# ------------------------------------------------------------

mmp13_pair_table_primary <-
  mmp13_pair_table %>%
  dplyr::mutate(
    pair_id =
      as.character(pair_id)
  )

mmp13_pair_table_sensitivity <-
  mmp13_pair_table_primary %>%
  dplyr::filter(
    pair_id !=
      excluded_sensitivity_pair
  )

mmp13_direct_summary <- dplyr::bind_rows(
  
  tibble::tibble(
    
    analysis =
      "Primary: all 33 complete pairs",
    
    number_of_pairs =
      nrow(
        mmp13_pair_table_primary
      ),
    
    mean_preserved =
      mean(
        mmp13_pair_table_primary$Preserved
      ),
    
    mean_OA =
      mean(
        mmp13_pair_table_primary$OA
      ),
    
    mean_paired_difference =
      mean(
        mmp13_pair_table_primary$
          difference_OA_minus_Preserved
      ),
    
    sd_paired_difference =
      sd(
        mmp13_pair_table_primary$
          difference_OA_minus_Preserved
      ),
    
    median_paired_difference =
      median(
        mmp13_pair_table_primary$
          difference_OA_minus_Preserved
      ),
    
    higher_in_OA_pairs =
      sum(
        mmp13_pair_table_primary$
          difference_OA_minus_Preserved > 0
      ),
    
    lower_in_OA_pairs =
      sum(
        mmp13_pair_table_primary$
          difference_OA_minus_Preserved < 0
      )
  ),
  
  tibble::tibble(
    
    analysis =
      "Sensitivity: 32 pairs excluding RAAK_14",
    
    number_of_pairs =
      nrow(
        mmp13_pair_table_sensitivity
      ),
    
    mean_preserved =
      mean(
        mmp13_pair_table_sensitivity$Preserved
      ),
    
    mean_OA =
      mean(
        mmp13_pair_table_sensitivity$OA
      ),
    
    mean_paired_difference =
      mean(
        mmp13_pair_table_sensitivity$
          difference_OA_minus_Preserved
      ),
    
    sd_paired_difference =
      sd(
        mmp13_pair_table_sensitivity$
          difference_OA_minus_Preserved
      ),
    
    median_paired_difference =
      median(
        mmp13_pair_table_sensitivity$
          difference_OA_minus_Preserved
      ),
    
    higher_in_OA_pairs =
      sum(
        mmp13_pair_table_sensitivity$
          difference_OA_minus_Preserved > 0
      ),
    
    lower_in_OA_pairs =
      sum(
        mmp13_pair_table_sensitivity$
          difference_OA_minus_Preserved < 0
      )
  )
)

# ------------------------------------------------------------
# 20. Gabungkan hasil MMP13
# ------------------------------------------------------------

mmp13_model_comparison <- dplyr::bind_rows(
  mmp13_primary_effect,
  mmp13_sensitivity_effect
) %>%
  dplyr::left_join(
    mmp13_direct_summary,
    by = c(
      "analysis",
      "number_of_pairs"
    )
  ) %>%
  dplyr::mutate(
    
    gene =
      "MMP13",
    
    probe_id =
      mmp13_probe_id,
    
    outlier_pair_excluded =
      dplyr::if_else(
        stringr::str_detect(
          analysis,
          "Sensitivity"
        ),
        excluded_sensitivity_pair,
        "None"
      ),
    
    coefficient_minus_direct_mean =
      effect -
      mean_paired_difference
  )

cat("\nMMP13 primary versus sensitivity results:\n")

print(
  mmp13_model_comparison,
  n = Inf,
  width = Inf
)

if (
  any(
    abs(
      mmp13_model_comparison$
      coefficient_minus_direct_mean
    ) > 1e-6
  )
) {
  warning(
    "Model coefficient berbeda dari direct mean paired difference. ",
    "Periksa sample ordering dan design matrix."
  )
}

# ------------------------------------------------------------
# 21. Candidate forest input
# ------------------------------------------------------------

# CATATAN:
# Effect ini masih berada pada VST-like scale.
# Jangan langsung digabungkan dengan log2FC dataset lain
# sebelum cross-study effect metric diharmonisasi.

mmp13_forest_candidate <- mmp13_primary_effect %>%
  dplyr::transmute(
    
    dataset =
      "GSE57218",
    
    tissue =
      "Human knee articular cartilage",
    
    comparison =
      "OA-affected versus preserved cartilage within the same OA joint",
    
    design =
      "Paired analysis",
    
    number_of_pairs =
      number_of_pairs,
    
    probe_id =
      mmp13_probe_id,
    
    gene =
      "MMP13",
    
    effect_measure =
      "Paired difference on processed VST-like scale",
    
    effect_method =
      "limma paired fixed-effect model with robust trend empirical Bayes",
    
    effect =
      effect,
    
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
      FDR,
    
    cross_study_pooling_status =
      "Requires effect-metric harmonization before pooling"
  )

cat("\nMMP13 candidate forest input:\n")

print(
  mmp13_forest_candidate,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 22. Gene-level primary versus sensitivity comparison
# ------------------------------------------------------------

gene_sensitivity_comparison <- gene_de_primary %>%
  dplyr::select(
    gene,
    primary_logFC = logFC,
    primary_PValue = PValue,
    primary_FDR = FDR
  ) %>%
  dplyr::inner_join(
    
    gene_de_sensitivity %>%
      dplyr::select(
        gene,
        sensitivity_logFC = logFC,
        sensitivity_PValue = PValue,
        sensitivity_FDR = FDR
      ),
    
    by = "gene"
  ) %>%
  dplyr::mutate(
    
    logFC_change_after_exclusion =
      sensitivity_logFC -
      primary_logFC,
    
    same_direction =
      sign(primary_logFC) ==
      sign(sensitivity_logFC)
  )

logFC_sensitivity_correlation <- stats::cor(
  gene_sensitivity_comparison$primary_logFC,
  gene_sensitivity_comparison$sensitivity_logFC,
  method = "pearson"
)

cat("\nCorrelation primary versus sensitivity gene effects:\n")
print(logFC_sensitivity_correlation)

# ------------------------------------------------------------
# 23. MMP13 primary versus sensitivity effect plot
# ------------------------------------------------------------

mmp13_plot_table <-
  mmp13_model_comparison %>%
  dplyr::mutate(
    analysis = factor(
      analysis,
      levels = c(
        "Primary: all 33 complete pairs",
        "Sensitivity: 32 pairs excluding RAAK_14"
      )
    )
  )

mmp13_effect_plot <- ggplot(
  mmp13_plot_table,
  aes(
    x = analysis,
    y = effect
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
    linewidth = 0.75
  ) +
  geom_point(
    size = 3.2
  ) +
  coord_flip() +
  labs(
    title =
      "Paired MMP13 effect estimates in GSE57218",
    
    subtitle =
      "OA-affected minus preserved cartilage",
    
    x = NULL,
    
    y =
      "Difference on processed VST-like expression scale"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(mmp13_effect_plot)

ggsave(
  filename =
    "results/figures/GSE57218_MMP13_primary_vs_sensitivity.pdf",
  
  plot =
    mmp13_effect_plot,
  
  width = 7.5,
  height = 4.2
)

ggsave(
  filename =
    "results/figures/GSE57218_MMP13_primary_vs_sensitivity.tiff",
  
  plot =
    mmp13_effect_plot,
  
  width = 7.5,
  height = 4.2,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 24. Volcano plot primary gene-level result
# ------------------------------------------------------------

top_higher_in_oa <- gene_de_primary %>%
  dplyr::filter(
    volcano_status ==
      "Higher in OA"
  ) %>%
  dplyr::slice_min(
    order_by = PValue,
    n = 5,
    with_ties = FALSE
  )

top_lower_in_oa <- gene_de_primary %>%
  dplyr::filter(
    volcano_status ==
      "Lower in OA"
  ) %>%
  dplyr::slice_min(
    order_by = PValue,
    n = 5,
    with_ties = FALSE
  )

volcano_label_genes <- dplyr::bind_rows(
  
  top_higher_in_oa,
  
  top_lower_in_oa,
  
  gene_de_primary %>%
    dplyr::filter(
      gene == "MMP13"
    )
) %>%
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )

volcano_plot <- ggplot(
  gene_de_primary,
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
    data =
      volcano_label_genes,
    
    aes(label = gene),
    
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "Paired differential expression in GSE57218",
    
    subtitle =
      "OA-affected versus preserved cartilage; 33 patient pairs",
    
    x =
      "Effect on processed VST-like expression scale",
    
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
    "results/figures/GSE57218_volcano_paired_primary.pdf",
  
  plot =
    volcano_plot,
  
  width = 7.2,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE57218_volcano_paired_primary.tiff",
  
  plot =
    volcano_plot,
  
  width = 7.2,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 25. MA plot primary
# ------------------------------------------------------------

ma_plot <- ggplot(
  gene_de_primary,
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
      gene_de_primary %>%
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
      "MA plot of paired GSE57218 analysis",
    
    subtitle =
      "OA-affected versus preserved cartilage",
    
    x =
      "Average processed expression",
    
    y =
      "OA minus preserved expression difference",
    
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
    "results/figures/GSE57218_MA_paired_primary.pdf",
  
  plot =
    ma_plot,
  
  width = 7.2,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE57218_MA_paired_primary.tiff",
  
  plot =
    ma_plot,
  
  width = 7.2,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 26. Mean-variance trend
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE57218_mean_variance_trend.pdf",
  width = 7,
  height = 5.5
)

limma::plotSA(
  gene_fit_primary,
  main =
    "Mean-variance trend in paired GSE57218 analysis"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE57218_mean_variance_trend.tiff",
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

limma::plotSA(
  gene_fit_primary,
  main =
    "Mean-variance trend in paired GSE57218 analysis"
)

dev.off()

# ------------------------------------------------------------
# 27. Primary versus sensitivity effect scatter plot
#     dengan warna perubahan dan highlight MMP13
# ------------------------------------------------------------

gene_sensitivity_comparison <- gene_sensitivity_comparison %>%
  dplyr::mutate(
    
    effect_shift =
      sensitivity_logFC -
      primary_logFC,
    
    absolute_effect_shift =
      abs(effect_shift),
    
    is_MMP13 =
      gene == "MMP13"
  )

# Data khusus MMP13
mmp13_scatter_data <-
  gene_sensitivity_comparison %>%
  dplyr::filter(
    gene == "MMP13"
  )

cat("\nPosisi MMP13 pada sensitivity plot:\n")

print(
  mmp13_scatter_data %>%
    dplyr::select(
      gene,
      primary_logFC,
      sensitivity_logFC,
      effect_shift,
      primary_PValue,
      primary_FDR,
      sensitivity_PValue,
      sensitivity_FDR
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# Membuat batas x dan y sama agar jarak dari garis diagonal
# dapat dibandingkan secara objektif
scatter_axis_limit <- max(
  abs(
    c(
      gene_sensitivity_comparison$primary_logFC,
      gene_sensitivity_comparison$sensitivity_logFC
    )
  ),
  na.rm = TRUE
) * 1.07

sensitivity_scatter_plot <- ggplot(
  gene_sensitivity_comparison,
  aes(
    x = primary_logFC,
    y = sensitivity_logFC
  )
) +
  
  # Semua gen, diwarna berdasarkan perubahan effect
  geom_point(
    aes(
      color = effect_shift
    ),
    alpha = 0.60,
    size = 1.35
  ) +
  
  # Garis kesamaan sempurna
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.65
  ) +
  
  # Garis nol vertikal dan horizontal
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.35
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.35
  ) +
  
  # Highlight MMP13
  geom_point(
    data = mmp13_scatter_data,
    aes(
      x = primary_logFC,
      y = sensitivity_logFC
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "#FFD54F",
    color = "black",
    stroke = 1.2,
    size = 2
  ) +
  
  # Label MMP13 agar tidak tertutup
  ggrepel::geom_label_repel(
    data = mmp13_scatter_data,
    aes(
      x = primary_logFC,
      y = sensitivity_logFC,
      label = paste0(
        "MMP13",
        "\nPrimary = ",
        round(primary_logFC, 3),
        "\nSensitivity = ",
        round(sensitivity_logFC, 3),
        "\nShift = ",
        round(effect_shift, 3)
      )
    ),
    inherit.aes = FALSE,
    size = 3.4,
    fontface = "bold",
    color = "black",
    fill = "white",
    label.size = 0.35,
    box.padding = 0.7,
    point.padding = 0.6,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = 0.45,
    nudge_x = -0.20,
    nudge_y = 0.18,
    seed = 123,
    show.legend = FALSE
  ) +
  
  # Gradient perubahan effect
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey80",
    high = "#B2182B",
    midpoint = 0,
    name = paste0(
      "Effect shift\n",
      "(32 pairs − 33 pairs)"
    )
  ) +
  
  # Menyamakan skala kedua sumbu
  coord_equal(
    xlim = c(
      -scatter_axis_limit,
      scatter_axis_limit
    ),
    ylim = c(
      -scatter_axis_limit,
      scatter_axis_limit
    ),
    expand = FALSE
  ) +
  
  labs(
    title =
      "Effect stability after excluding RAAK_14",
    
    subtitle = paste0(
      "Pearson correlation = ",
      round(
        logFC_sensitivity_correlation,
        4
      ),
      "; color indicates the change after pair exclusion"
    ),
    
    x =
      "Primary effect: all 33 pairs",
    
    y =
      "Sensitivity effect: 32 pairs"
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "right",
    
    plot.title =
      element_text(
        face = "bold"
      )
  )

print(
  sensitivity_scatter_plot
)

# Simpan dengan nama lama agar otomatis menggantikan grafik sebelumnya
ggsave(
  filename =
    "results/figures/GSE57218_primary_vs_sensitivity_effects.pdf",
  
  plot =
    sensitivity_scatter_plot,
  
  width = 7.5,
  height = 6.2
)

ggsave(
  filename =
    "results/figures/GSE57218_primary_vs_sensitivity_effects.tiff",
  
  plot =
    sensitivity_scatter_plot,
  
  width = 7.5,
  height = 6.2,
  dpi = 600,
  compression = "lzw"
)
# ------------------------------------------------------------
# 28. Ranked gene lists untuk pathway analysis
# ------------------------------------------------------------

gsea_ranked_primary <- gene_de_primary %>%
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

gsea_ranked_sensitivity <- gene_de_sensitivity %>%
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

cat("\nSepuluh ranked genes tertinggi:\n")
print(
  head(
    gsea_ranked_primary,
    10
  )
)

cat("\nSepuluh ranked genes terendah:\n")
print(
  tail(
    gsea_ranked_primary,
    10
  )
)

# ------------------------------------------------------------
# 29. Significant result tables
# ------------------------------------------------------------

significant_gene_primary <- gene_de_primary %>%
  dplyr::filter(
    FDR < 0.05
  )

significant_gene_primary_effect1 <- gene_de_primary %>%
  dplyr::filter(
    FDR < 0.05,
    abs(logFC) >= 1
  )

significant_gene_sensitivity <- gene_de_sensitivity %>%
  dplyr::filter(
    FDR < 0.05
  )

top_100_probe_primary <- probe_de_primary %>%
  dplyr::slice_head(n = 100)

top_100_probe_sensitivity <- probe_de_sensitivity %>%
  dplyr::slice_head(n = 100)

# ------------------------------------------------------------
# 30. Menyimpan R objects
# ------------------------------------------------------------

saveRDS(
  probe_fit_primary,
  file =
    "data_processed/GSE57218_probe_limma_fit_primary_33_pairs.rds"
)

saveRDS(
  probe_fit_sensitivity,
  file =
    "data_processed/GSE57218_probe_limma_fit_sensitivity_32_pairs.rds"
)

saveRDS(
  gene_fit_primary,
  file =
    "data_processed/GSE57218_gene_limma_fit_primary_33_pairs.rds"
)

saveRDS(
  gene_fit_sensitivity,
  file =
    "data_processed/GSE57218_gene_limma_fit_sensitivity_32_pairs.rds"
)

saveRDS(
  probe_de_primary,
  file =
    "data_processed/GSE57218_probe_DE_primary.rds"
)

saveRDS(
  probe_de_sensitivity,
  file =
    "data_processed/GSE57218_probe_DE_sensitivity.rds"
)

saveRDS(
  gene_de_primary,
  file =
    "data_processed/GSE57218_gene_DE_primary.rds"
)

saveRDS(
  gene_de_sensitivity,
  file =
    "data_processed/GSE57218_gene_DE_sensitivity.rds"
)

saveRDS(
  gene_de_primary_treat,
  file =
    "data_processed/GSE57218_gene_TREAT_primary_lfc1.rds"
)

saveRDS(
  gene_de_sensitivity_treat,
  file =
    "data_processed/GSE57218_gene_TREAT_sensitivity_lfc1.rds"
)

saveRDS(
  mmp13_model_comparison,
  file =
    "data_processed/GSE57218_MMP13_primary_sensitivity_results.rds"
)

saveRDS(
  mmp13_forest_candidate,
  file =
    "data_processed/GSE57218_MMP13_forest_candidate.rds"
)

saveRDS(
  gsea_ranked_primary,
  file =
    "data_processed/GSE57218_GSEA_ranked_genes_primary.rds"
)

saveRDS(
  gsea_ranked_sensitivity,
  file =
    "data_processed/GSE57218_GSEA_ranked_genes_sensitivity.rds"
)

saveRDS(
  design_sensitivity,
  file =
    "data_processed/GSE57218_design_sensitivity_without_RAAK14.rds"
)

# ------------------------------------------------------------
# 31. Simpan full probe tables sebagai CSV
# ------------------------------------------------------------

utils::write.csv(
  probe_de_primary,
  file =
    "results/tables/GSE57218_probe_DE_primary_33_pairs.csv",
  row.names = FALSE
)

utils::write.csv(
  probe_de_sensitivity,
  file =
    "results/tables/GSE57218_probe_DE_sensitivity_32_pairs.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 32. Menyimpan Excel workbook
# ------------------------------------------------------------

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
    
    MMP13_Direct_Summary =
      as.data.frame(
        mmp13_direct_summary
      ),
    
    MMP13_Forest_Candidate =
      as.data.frame(
        mmp13_forest_candidate
      ),
    
    Gene_DE_Primary =
      as.data.frame(
        gene_de_primary
      ),
    
    Gene_DE_Sensitivity =
      as.data.frame(
        gene_de_sensitivity
      ),
    
    Gene_TREAT_Primary =
      as.data.frame(
        gene_de_primary_treat
      ),
    
    Gene_TREAT_Sensitivity =
      as.data.frame(
        gene_de_sensitivity_treat
      ),
    
    Significant_Primary =
      as.data.frame(
        significant_gene_primary
      ),
    
    Significant_Primary_Effect1 =
      as.data.frame(
        significant_gene_primary_effect1
      ),
    
    Significant_Sensitivity =
      as.data.frame(
        significant_gene_sensitivity
      ),
    
    GSEA_Ranked_Primary =
      as.data.frame(
        gsea_ranked_primary
      ),
    
    GSEA_Ranked_Sensitivity =
      as.data.frame(
        gsea_ranked_sensitivity
      ),
    
    Gene_Effect_Comparison =
      as.data.frame(
        gene_sensitivity_comparison
      ),
    
    Top_100_Probe_Primary =
      as.data.frame(
        top_100_probe_primary
      ),
    
    Top_100_Probe_Sensitivity =
      as.data.frame(
        top_100_probe_sensitivity
      )
  ),
  
  file =
    "results/tables/GSE57218_limma_paired_DE_results.xlsx",
  
  overwrite = TRUE
)

openxlsx::write.xlsx(
  as.data.frame(
    mmp13_forest_candidate
  ),
  
  file =
    "results/tables/GSE57218_MMP13_forest_candidate.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 33. Session information
# ------------------------------------------------------------

sink(
  "results/tables/GSE57218_DE_sessionInfo.txt"
)

print(sessionInfo())

sink()

# ------------------------------------------------------------
# 34. Pesan akhir
# ------------------------------------------------------------

message("")
message("================================================")
message("GSE57218 PAIRED DIFFERENTIAL EXPRESSION SELESAI")
message("================================================")

message(
  "Jumlah gen primary diuji       : ",
  nrow(gene_de_primary)
)

message(
  "Primary FDR < 0.05             : ",
  sum(
    gene_de_primary$FDR < 0.05
  )
)

message(
  "Primary FDR & |effect| >= 1    : ",
  sum(
    gene_de_primary$FDR < 0.05 &
      abs(gene_de_primary$logFC) >= 1
  )
)

message(
  "Primary TREAT FDR < 0.05       : ",
  sum(
    gene_de_primary_treat$FDR < 0.05
  )
)

message(
  "Sensitivity FDR < 0.05         : ",
  sum(
    gene_de_sensitivity$FDR < 0.05
  )
)

message(
  "MMP13 primary effect           : ",
  round(
    mmp13_primary_effect$effect,
    4
  )
)

message(
  "MMP13 primary SE               : ",
  round(
    mmp13_primary_effect$standard_error,
    4
  )
)

message(
  "MMP13 primary 95% CI           : ",
  round(
    mmp13_primary_effect$CI_95_lower,
    4
  ),
  " to ",
  round(
    mmp13_primary_effect$CI_95_upper,
    4
  )
)

message(
  "MMP13 primary FDR              : ",
  format(
    mmp13_primary_effect$FDR,
    scientific = TRUE,
    digits = 4
  )
)

message(
  "MMP13 sensitivity effect       : ",
  round(
    mmp13_sensitivity_effect$effect,
    4
  )
)

message(
  "MMP13 sensitivity 95% CI       : ",
  round(
    mmp13_sensitivity_effect$CI_95_lower,
    4
  ),
  " to ",
  round(
    mmp13_sensitivity_effect$CI_95_upper,
    4
  )
)

message(
  "Primary-sensitivity correlation: ",
  round(
    logFC_sensitivity_correlation,
    4
  )
)

message(
  "Forest pooling status          : ",
  "Requires effect-metric harmonization"
)

message("================================================")