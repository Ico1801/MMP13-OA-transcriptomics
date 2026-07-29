# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 13_GSE169077_limma_DE.R
# PURPOSE : Differential-expression analysis of GSE169077
#
# COMPARISON:
#   Late-stage OA cartilage versus normal cartilage
#
# STATISTICAL UNIT:
#   RNA pool
#
# SAMPLE SIZE:
#   5 normal pools versus 6 OA pools
#
# PRIMARY MODEL:
#   ~ group
#
# ADDITIONAL ANALYSIS:
#   Leave-one-pool-out influence analysis for MMP13
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
  "Semua paket differential expression GSE169077 berhasil dimuat."
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

expression_rma <- readRDS(
  "data_processed/GSE169077_expression_RMA_11_pools.rds"
)

expression_biological <- readRDS(
  "data_processed/GSE169077_expression_RMA_annotated_biological.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE169077_metadata_RMA_QC.rds"
)

annotation_biological <- readRDS(
  "data_processed/GSE169077_annotation_biological.rds"
)

mmp13_annotation <- readRDS(
  "data_processed/GSE169077_MMP13_probe_annotation.rds"
)

sample_connectivity <- readRDS(
  "data_processed/GSE169077_sample_connectivity.rds"
)

cat("\nDimensi full RMA matrix:\n")
print(
  dim(
    expression_rma
  )
)

cat("\nDimensi annotated biological matrix:\n")
print(
  dim(
    expression_biological
  )
)

cat("\nDimensi metadata:\n")
print(
  dim(
    sample_metadata
  )
)

cat("\nDimensi annotation biological:\n")
print(
  dim(
    annotation_biological
  )
)

# ------------------------------------------------------------
# 4. Validasi data
# ------------------------------------------------------------

stopifnot(
  is.matrix(
    expression_biological
  )
)

stopifnot(
  ncol(expression_biological) == 11
)

stopifnot(
  nrow(sample_metadata) == 11
)

stopifnot(
  identical(
    colnames(expression_biological),
    sample_metadata$geo_accession
  )
)

stopifnot(
  !anyNA(expression_biological)
)

stopifnot(
  all(
    is.finite(
      expression_biological
    )
  )
)

stopifnot(
  sum(
    sample_metadata$group ==
      "Normal"
  ) == 5
)

stopifnot(
  sum(
    sample_metadata$group ==
      "OA"
  ) == 6
)

stopifnot(
  sum(
    sample_connectivity$
      candidate_sample_outlier
  ) == 0
)

stopifnot(
  nrow(mmp13_annotation) >= 1
)

message(
  "Expression matrix, metadata, dan hasil QC berhasil divalidasi."
)

# ------------------------------------------------------------
# 5. Menyiapkan metadata dan design matrix
# ------------------------------------------------------------

sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    
    group = factor(
      as.character(group),
      levels = c(
        "Normal",
        "OA"
      )
    ),
    
    short_label = as.character(
      short_label
    )
  )

design_primary <- model.matrix(
  ~ group,
  data = sample_metadata
)

rownames(design_primary) <-
  sample_metadata$geo_accession

cat("\nPrimary design matrix:\n")
print(
  design_primary
)

cat("\nPrimary design dimensions:\n")
print(
  dim(
    design_primary
  )
)

cat("\nPrimary design rank:\n")
print(
  qr(
    design_primary
  )$rank
)

cat("\nPrimary residual degrees of freedom:\n")
print(
  nrow(design_primary) -
    qr(design_primary)$rank
)

stopifnot(
  qr(design_primary)$rank ==
    ncol(design_primary)
)

stopifnot(
  "groupOA" %in%
    colnames(design_primary)
)

stopifnot(
  identical(
    rownames(design_primary),
    colnames(expression_biological)
  )
)

# ------------------------------------------------------------
# 6. Menyelaraskan annotation dengan expression matrix
# ------------------------------------------------------------

annotation_index <- match(
  rownames(expression_biological),
  annotation_biological$probe_id
)

if (anyNA(annotation_index)) {
  stop(
    "Ada probeset expression yang tidak ditemukan ",
    "dalam annotation biological."
  )
}

annotation_aligned <- annotation_biological[
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

# ------------------------------------------------------------
# 7. Fungsi mengambil kolom anotasi secara aman
# ------------------------------------------------------------

get_annotation_column <- function(
    annotation_table,
    candidate_names
) {
  
  available_column <- intersect(
    candidate_names,
    colnames(annotation_table)
  )
  
  if (length(available_column) == 0) {
    return(
      rep(
        NA_character_,
        nrow(annotation_table)
      )
    )
  }
  
  as.character(
    annotation_table[[available_column[1]]]
  )
}

probe_annotation <- tibble::tibble(
  
  probe_id =
    as.character(
      annotation_aligned$probe_id
    ),
  
  gene_symbol =
    trimws(
      get_annotation_column(
        annotation_aligned,
        c(
          "Gene Symbol",
          "GENE_SYMBOL",
          "Symbol",
          "gene_symbol"
        )
      )
    ),
  
  gene_title =
    get_annotation_column(
      annotation_aligned,
      c(
        "Gene Title",
        "GENE_TITLE",
        "Gene title",
        "Description"
      )
    ),
  
  entrez_id =
    get_annotation_column(
      annotation_aligned,
      c(
        "Gene ID",
        "ENTREZ_GENE_ID",
        "Entrez Gene",
        "Entrez_Gene_ID"
      )
    ),
  
  accession =
    get_annotation_column(
      annotation_aligned,
      c(
        "GB_ACC",
        "GeneBank Accession",
        "Accession"
      )
    )
)

cat("\nProbe annotation summary:\n")

print(
  tibble::tibble(
    metric = c(
      "Total probesets",
      "Non-empty gene symbols",
      "Unique gene symbols"
    ),
    value = c(
      nrow(probe_annotation),
      sum(
        !is.na(
          probe_annotation$gene_symbol
        ) &
          probe_annotation$gene_symbol != ""
      ),
      dplyr::n_distinct(
        probe_annotation$gene_symbol[
          !is.na(
            probe_annotation$gene_symbol
          ) &
            probe_annotation$gene_symbol != ""
        ]
      )
    )
  ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 8. Fungsi mengambil topTable lengkap
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
    ) %>%
    tibble::as_tibble()
  
  result_table
}

# ------------------------------------------------------------
# 9. Primary probe-level limma model
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

probe_fit_primary_treat <- limma::treat(
  probe_fit_primary_raw,
  lfc = 1,
  trend = TRUE,
  robust = TRUE
)

cat("\nProbe-level coefficients:\n")
print(
  colnames(
    probe_fit_primary$coefficients
  )
)

cat("\nProbe-level residual df summary:\n")
print(
  summary(
    probe_fit_primary$df.residual
  )
)

probe_de_primary <- extract_limma_table(
  fitted_model =
    probe_fit_primary,
  
  coefficient_name =
    "groupOA"
) %>%
  dplyr::rename(
    probe_id = feature_id
  ) %>%
  dplyr::left_join(
    probe_annotation,
    by = "probe_id"
  )

probe_treat_primary <- limma::topTreat(
  probe_fit_primary_treat,
  coef = "groupOA",
  number = Inf,
  sort.by = "P"
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "probe_id"
  ) %>%
  dplyr::rename(
    PValue = P.Value,
    FDR = adj.P.Val
  ) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(
    probe_annotation,
    by = "probe_id"
  )

# ------------------------------------------------------------
# 10. Mean expression per group untuk setiap probeset
# ------------------------------------------------------------

normal_pool_ids <-
  sample_metadata$geo_accession[
    sample_metadata$group ==
      "Normal"
  ]

oa_pool_ids <-
  sample_metadata$geo_accession[
    sample_metadata$group ==
      "OA"
  ]

probe_expression_summary <- tibble::tibble(
  
  probe_id =
    rownames(
      expression_biological
    ),
  
  mean_normal =
    rowMeans(
      expression_biological[
        ,
        normal_pool_ids,
        drop = FALSE
      ]
    ),
  
  mean_OA =
    rowMeans(
      expression_biological[
        ,
        oa_pool_ids,
        drop = FALSE
      ]
    )
) %>%
  dplyr::mutate(
    
    direct_log2_difference =
      mean_OA -
      mean_normal,
    
    approximate_fold_change =
      2^direct_log2_difference
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
# 11. Menyiapkan gene symbol yang tidak ambigu
# ------------------------------------------------------------

gene_symbol_clean <- trimws(
  probe_annotation$gene_symbol
)

ambiguous_gene_symbol <- stringr::str_detect(
  gene_symbol_clean,
  "///|//|;|\\||,"
)

keep_unambiguous_gene_symbol <- (
  !is.na(gene_symbol_clean) &
    gene_symbol_clean != "" &
    !ambiguous_gene_symbol
)

cat("\nProbeset dengan unambiguous gene symbol:\n")
print(
  sum(
    keep_unambiguous_gene_symbol
  )
)

cat("\nProbeset dengan ambiguous gene symbol:\n")
print(
  sum(
    ambiguous_gene_symbol,
    na.rm = TRUE
  )
)

# ------------------------------------------------------------
# 12. Collapse probesets ke gene-level expression
# ------------------------------------------------------------

expression_gene <- limma::avereps(
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

cat("\nDimensi gene-level expression matrix:\n")
print(
  dim(
    expression_gene
  )
)

stopifnot(
  "MMP13" %in%
    rownames(expression_gene)
)

# ------------------------------------------------------------
# 13. Gene-level limma model
# ------------------------------------------------------------

gene_fit_primary_raw <- limma::lmFit(
  expression_gene,
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
    "groupOA"
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
      ),
    
    approximate_fold_change =
      2^logFC
  )

gene_treat_primary <- limma::topTreat(
  gene_fit_primary_treat,
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
  ) %>%
  tibble::as_tibble()

# ------------------------------------------------------------
# 14. Differential-expression summary
# ------------------------------------------------------------

de_summary <- tibble::tibble(
  
  analysis_level = c(
    "Probe level",
    "Gene level"
  ),
  
  total_features_tested = c(
    nrow(probe_de_primary),
    nrow(gene_de_primary)
  ),
  
  FDR_less_than_0_05 = c(
    sum(
      probe_de_primary$FDR < 0.05,
      na.rm = TRUE
    ),
    sum(
      gene_de_primary$FDR < 0.05,
      na.rm = TRUE
    )
  ),
  
  higher_in_OA_FDR = c(
    sum(
      probe_de_primary$FDR < 0.05 &
        probe_de_primary$logFC > 0,
      na.rm = TRUE
    ),
    sum(
      gene_de_primary$FDR < 0.05 &
        gene_de_primary$logFC > 0,
      na.rm = TRUE
    )
  ),
  
  lower_in_OA_FDR = c(
    sum(
      probe_de_primary$FDR < 0.05 &
        probe_de_primary$logFC < 0,
      na.rm = TRUE
    ),
    sum(
      gene_de_primary$FDR < 0.05 &
        gene_de_primary$logFC < 0,
      na.rm = TRUE
    )
  ),
  
  FDR_and_abs_log2FC_at_least_1 = c(
    sum(
      probe_de_primary$FDR < 0.05 &
        abs(probe_de_primary$logFC) >= 1,
      na.rm = TRUE
    ),
    sum(
      gene_de_primary$FDR < 0.05 &
        abs(gene_de_primary$logFC) >= 1,
      na.rm = TRUE
    )
  ),
  
  TREAT_FDR_less_than_0_05 = c(
    sum(
      probe_treat_primary$FDR < 0.05,
      na.rm = TRUE
    ),
    sum(
      gene_treat_primary$FDR < 0.05,
      na.rm = TRUE
    )
  )
)

cat("\nDifferential-expression summary:\n")

print(
  de_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 15. Fungsi ekstraksi effect estimate dan moderated SE
# ------------------------------------------------------------

extract_limma_effect <- function(
    fitted_model,
    feature_id,
    coefficient_name,
    analysis_name,
    number_of_units
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
    0.975,
    df =
      degrees_of_freedom
  )
  
  statistics_table <- limma::topTable(
    fitted_model,
    coef = coefficient_name,
    number = Inf,
    sort.by = "none"
  )
  
  feature_statistics <-
    statistics_table[
      feature_id,
      ,
      drop = FALSE
    ]
  
  tibble::tibble(
    
    analysis =
      analysis_name,
    
    number_of_statistical_units =
      number_of_units,
    
    feature_id =
      feature_id,
    
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
# 16. MMP13 probe identification
# ------------------------------------------------------------

mmp13_probe_ids <- intersect(
  unique(
    as.character(
      mmp13_annotation$probe_id
    )
  ),
  rownames(
    expression_biological
  )
)

cat("\nMMP13 probes in expression matrix:\n")
print(
  mmp13_probe_ids
)

if (length(mmp13_probe_ids) != 1) {
  stop(
    "Diharapkan tepat satu probe MMP13, ditemukan: ",
    length(mmp13_probe_ids)
  )
}

mmp13_probe_id <-
  mmp13_probe_ids[1]

# ------------------------------------------------------------
# 17. Primary MMP13 effect
# ------------------------------------------------------------

mmp13_primary_effect <- extract_limma_effect(
  
  fitted_model =
    probe_fit_primary,
  
  feature_id =
    mmp13_probe_id,
  
  coefficient_name =
    "groupOA",
  
  analysis_name =
    "Primary model: all 11 RNA pools",
  
  number_of_units =
    11
) %>%
  dplyr::mutate(
    
    gene =
      "MMP13",
    
    probe_id =
      mmp13_probe_id,
    
    effect_direction =
      "OA minus Normal",
    
    effect_measure =
      "log2 expression difference",
    
    approximate_fold_change =
      2^effect,
    
    normal_pools =
      5L,
    
    OA_pools =
      6L,
    
    statistical_unit =
      "RNA pool"
  )

cat("\nPrimary MMP13 effect:\n")

print(
  mmp13_primary_effect,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 18. Memastikan model effect sama dengan mean difference
# ------------------------------------------------------------

mmp13_expression <- tibble::tibble(
  
  geo_accession =
    sample_metadata$geo_accession,
  
  short_label =
    sample_metadata$short_label,
  
  group =
    sample_metadata$group,
  
  MMP13_expression =
    as.numeric(
      expression_biological[
        mmp13_probe_id,
        sample_metadata$geo_accession
      ]
    )
)

mmp13_direct_summary <- mmp13_expression %>%
  dplyr::group_by(
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
    
    .groups =
      "drop"
  )

mmp13_direct_effect <- mmp13_direct_summary %>%
  dplyr::select(
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
    
    direct_log2_difference =
      OA -
      Normal,
    
    direct_fold_change =
      2^direct_log2_difference,
    
    coefficient_minus_direct_difference =
      mmp13_primary_effect$effect -
      direct_log2_difference
  )

cat("\nMMP13 direct group summary:\n")

print(
  mmp13_direct_summary,
  n = Inf,
  width = Inf
)

cat("\nMMP13 direct effect check:\n")

print(
  mmp13_direct_effect,
  n = Inf,
  width = Inf
)

if (
  abs(
    mmp13_direct_effect$
    coefficient_minus_direct_difference
  ) > 1e-8
) {
  warning(
    "Model coefficient berbeda dari direct mean difference."
  )
}

# ------------------------------------------------------------
# 19. Leave-one-pool-out influence analysis
# ------------------------------------------------------------

message(
  "Menjalankan leave-one-pool-out influence analysis..."
)

leave_one_out_results <- lapply(
  seq_len(
    nrow(sample_metadata)
  ),
  function(excluded_index) {
    
    # Menggunakan nama berbeda agar tidak berbenturan
    # dengan nama kolom hasil
    excluded_sample_row <-
      sample_metadata[
        excluded_index,
        ,
        drop = FALSE
      ]
    
    excluded_short_label <-
      as.character(
        excluded_sample_row$short_label
      )
    
    excluded_geo_accession <-
      as.character(
        excluded_sample_row$geo_accession
      )
    
    excluded_group_name <-
      as.character(
        excluded_sample_row$group
      )
    
    metadata_loo <-
      sample_metadata[
        -excluded_index,
        ,
        drop = FALSE
      ] %>%
      droplevels()
    
    metadata_loo$group <- factor(
      as.character(
        metadata_loo$group
      ),
      levels = c(
        "Normal",
        "OA"
      )
    )
    
    expression_loo <-
      expression_biological[
        ,
        metadata_loo$geo_accession,
        drop = FALSE
      ]
    
    design_loo <- model.matrix(
      ~ group,
      data = metadata_loo
    )
    
    rownames(design_loo) <-
      metadata_loo$geo_accession
    
    if (
      qr(design_loo)$rank !=
      ncol(design_loo)
    ) {
      stop(
        "Leave-one-out design tidak full rank setelah mengeluarkan ",
        excluded_short_label
      )
    }
    
    fit_loo_raw <- limma::lmFit(
      expression_loo,
      design = design_loo
    )
    
    fit_loo <- limma::eBayes(
      fit_loo_raw,
      trend = TRUE,
      robust = TRUE
    )
    
    current_effect <- extract_limma_effect(
      
      fitted_model =
        fit_loo,
      
      feature_id =
        mmp13_probe_id,
      
      coefficient_name =
        "groupOA",
      
      analysis_name =
        paste0(
          "Excluding ",
          excluded_short_label
        ),
      
      number_of_units =
        nrow(metadata_loo)
    )
    
    current_effect %>%
      dplyr::mutate(
        
        excluded_pool =
          excluded_short_label,
        
        excluded_geo_accession =
          excluded_geo_accession,
        
        excluded_group =
          excluded_group_name,
        
        remaining_normal_pools =
          sum(
            metadata_loo$group ==
              "Normal"
          ),
        
        remaining_OA_pools =
          sum(
            metadata_loo$group ==
              "OA"
          )
      )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    
    effect_shift_from_primary =
      effect -
      mmp13_primary_effect$effect,
    
    absolute_effect_shift =
      abs(
        effect_shift_from_primary
      ),
    
    direction_consistent =
      sign(effect) ==
      sign(
        mmp13_primary_effect$effect
      )
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      absolute_effect_shift
    )
  )

cat("\nMMP13 leave-one-pool-out results:\n")

print(
  leave_one_out_results %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

mmp13_influence_summary <- tibble::tibble(
  
  primary_effect =
    mmp13_primary_effect$effect,
  
  minimum_leave_one_out_effect =
    min(
      leave_one_out_results$effect
    ),
  
  maximum_leave_one_out_effect =
    max(
      leave_one_out_results$effect
    ),
  
  maximum_absolute_shift =
    max(
      leave_one_out_results$
        absolute_effect_shift
    ),
  
  direction_consistent_in_all_models =
    all(
      leave_one_out_results$
        direction_consistent
    ),
  
  minimum_leave_one_out_PValue =
    min(
      leave_one_out_results$PValue
    ),
  
  maximum_leave_one_out_PValue =
    max(
      leave_one_out_results$PValue
    )
)

cat("\nMMP13 influence summary:\n")

print(
  mmp13_influence_summary %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 20. Candidate forest input
# ------------------------------------------------------------

mmp13_forest_candidate <- mmp13_primary_effect %>%
  dplyr::transmute(
    
    dataset =
      "GSE169077",
    
    platform =
      "GPL96 Affymetrix HG-U133A",
    
    tissue =
      "Human articular cartilage",
    
    comparison =
      "Late-stage OA versus normal cartilage",
    
    design =
      "Unpaired pool-level comparison",
    
    normal_pools =
      normal_pools,
    
    OA_pools =
      OA_pools,
    
    individuals_per_pool =
      5L,
    
    statistical_unit =
      statistical_unit,
    
    gene =
      gene,
    
    probe_id =
      probe_id,
    
    effect_measure =
      "log2 fold change from RMA expression",
    
    effect_direction =
      effect_direction,
    
    effect =
      effect,
    
    standard_error =
      standard_error,
    
    CI_95_lower =
      CI_95_lower,
    
    CI_95_upper =
      CI_95_upper,
    
    approximate_fold_change =
      approximate_fold_change,
    
    degrees_of_freedom =
      degrees_of_freedom,
    
    PValue =
      PValue,
    
    FDR =
      FDR,
    
    analysis_method =
      paste(
        "limma with trend and robust empirical Bayes;",
        "RNA pools treated as statistical units"
      ),
    
    pooling_note =
      paste(
        "Pool-level uncertainty;",
        "do not treat contributing individuals as independent replicates"
      )
  )

cat("\nMMP13 candidate forest input:\n")

print(
  mmp13_forest_candidate,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 21. MMP13 effect plot
# ------------------------------------------------------------

mmp13_effect_plot_table <- mmp13_primary_effect %>%
  dplyr::mutate(
    display_analysis =
      "GSE169077: 5 normal vs 6 OA pools"
  )

mmp13_effect_plot <- ggplot(
  mmp13_effect_plot_table,
  aes(
    x = display_analysis,
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
    width = 0.12,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5
  ) +
  coord_flip() +
  labs(
    title =
      "MMP13 effect estimate in GSE169077",
    
    subtitle =
      "Late-stage OA minus normal cartilage",
    
    x = NULL,
    
    y =
      expression(
        log[2]~fold~change
      )
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(
  mmp13_effect_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_MMP13_effect_estimate.pdf",
  plot =
    mmp13_effect_plot,
  width = 7,
  height = 4.2
)

ggsave(
  filename =
    "results/figures/GSE169077_MMP13_effect_estimate.tiff",
  plot =
    mmp13_effect_plot,
  width = 7,
  height = 4.2,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 22. Leave-one-pool-out plot
# ------------------------------------------------------------

leave_one_out_plot_table <- leave_one_out_results %>%
  dplyr::mutate(
    
    excluded_pool = factor(
      excluded_pool,
      levels =
        excluded_pool[
          order(effect)
        ]
    ),
    
    excluded_label = paste0(
      excluded_pool,
      " (",
      excluded_group,
      ")"
    )
  )

leave_one_out_plot <- ggplot(
  leave_one_out_plot_table,
  aes(
    x = reorder(
      excluded_label,
      effect
    ),
    y = effect
  )
) +
  geom_hline(
    yintercept =
      mmp13_primary_effect$effect,
    
    linetype =
      "solid",
    
    linewidth =
      0.6
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_errorbar(
    aes(
      ymin = CI_95_lower,
      ymax = CI_95_upper
    ),
    width = 0.15,
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      shape = excluded_group
    ),
    size = 3
  ) +
  coord_flip() +
  labs(
    title =
      "Leave-one-pool-out influence on MMP13",
    
    subtitle =
      "Solid line indicates the primary estimate using all 11 RNA pools",
    
    x =
      "Excluded RNA pool",
    
    y =
      expression(
        MMP13~log[2]~fold~change
      ),
    
    shape =
      "Excluded group"
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
  leave_one_out_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_MMP13_leave_one_pool_out.pdf",
  plot =
    leave_one_out_plot,
  width = 7.5,
  height = 6
)

ggsave(
  filename =
    "results/figures/GSE169077_MMP13_leave_one_pool_out.tiff",
  plot =
    leave_one_out_plot,
  width = 7.5,
  height = 6,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 23. Volcano plot gene-level results
# ------------------------------------------------------------

top_higher_genes <- gene_de_primary %>%
  dplyr::filter(
    statistical_status ==
      "Higher in OA"
  ) %>%
  dplyr::slice_min(
    order_by =
      PValue,
    
    n = 5,
    
    with_ties =
      FALSE
  )

top_lower_genes <- gene_de_primary %>%
  dplyr::filter(
    statistical_status ==
      "Lower in OA"
  ) %>%
  dplyr::slice_min(
    order_by =
      PValue,
    
    n = 5,
    
    with_ties =
      FALSE
  )

volcano_label_genes <- dplyr::bind_rows(
  
  top_higher_genes,
  
  top_lower_genes,
  
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
    color = statistical_status
  )
) +
  geom_point(
    alpha = 0.60,
    size = 1.3
  ) +
  geom_vline(
    xintercept = c(
      -1,
      1
    ),
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_text(
    data =
      volcano_label_genes,
    
    aes(
      label = gene
    ),
    
    size = 3,
    vjust = -0.7,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title =
      "Differential expression in GSE169077",
    
    subtitle =
      "Late-stage OA versus normal cartilage; colors indicate FDR < 0.05",
    
    x =
      expression(
        log[2]~fold~change
      ),
    
    y =
      expression(
        -log[10]~P-value
      ),
    
    color =
      "Statistical status"
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
  volcano_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_volcano_limma.pdf",
  plot =
    volcano_plot,
  width = 7.2,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE169077_volcano_limma.tiff",
  plot =
    volcano_plot,
  width = 7.2,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 24. MA plot
# ------------------------------------------------------------

ma_label_mmp13 <- gene_de_primary %>%
  dplyr::filter(
    gene == "MMP13"
  )

ma_plot <- ggplot(
  gene_de_primary,
  aes(
    x = AveExpr,
    y = logFC,
    color = statistical_status
  )
) +
  geom_point(
    alpha = 0.60,
    size = 1.3
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = c(
      -1,
      1
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_text(
    data =
      ma_label_mmp13,
    
    aes(
      label = gene
    ),
    
    size = 3.2,
    vjust = -0.8,
    show.legend = FALSE
  ) +
  labs(
    title =
      "MA plot of GSE169077",
    
    subtitle =
      "Late-stage OA versus normal cartilage",
    
    x =
      "Average RMA expression",
    
    y =
      expression(
        log[2]~fold~change
      ),
    
    color =
      "Statistical status"
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
  ma_plot
)

ggsave(
  filename =
    "results/figures/GSE169077_MA_limma.pdf",
  plot =
    ma_plot,
  width = 7.2,
  height = 5.8
)

ggsave(
  filename =
    "results/figures/GSE169077_MA_limma.tiff",
  plot =
    ma_plot,
  width = 7.2,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)

# ------------------------------------------------------------
# 25. Mean-variance trend
# ------------------------------------------------------------

pdf(
  file =
    "results/figures/GSE169077_mean_variance_trend.pdf",
  width = 7,
  height = 5.5
)

limma::plotSA(
  gene_fit_primary,
  main =
    "Mean-variance trend in GSE169077"
)

dev.off()

tiff(
  filename =
    "results/figures/GSE169077_mean_variance_trend.tiff",
  width = 7,
  height = 5.5,
  units = "in",
  res = 600,
  compression = "lzw"
)

limma::plotSA(
  gene_fit_primary,
  main =
    "Mean-variance trend in GSE169077"
)

dev.off()

# ------------------------------------------------------------
# 26. Ranked gene list untuk pathway enrichment
# ------------------------------------------------------------

gsea_ranked_genes <- gene_de_primary %>%
  dplyr::select(
    gene,
    moderated_t = t,
    logFC,
    AveExpr,
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
    gsea_ranked_genes,
    10
  ) %>%
    tibble::as_tibble(),
  n = 10,
  width = Inf
)

cat("\nSepuluh ranked genes terendah:\n")

print(
  tail(
    gsea_ranked_genes,
    10
  ) %>%
    tibble::as_tibble(),
  n = 10,
  width = Inf
)

# ------------------------------------------------------------
# 27. Significant result tables
# ------------------------------------------------------------

significant_gene_FDR <- gene_de_primary %>%
  dplyr::filter(
    FDR < 0.05
  )

significant_gene_FDR_logFC1 <- gene_de_primary %>%
  dplyr::filter(
    FDR < 0.05,
    abs(logFC) >= 1
  )

significant_gene_TREAT <- gene_treat_primary %>%
  dplyr::filter(
    FDR < 0.05
  )

top_100_probe_results <- probe_de_primary %>%
  dplyr::slice_head(
    n = 100
  )

top_100_gene_results <- gene_de_primary %>%
  dplyr::slice_head(
    n = 100
  )

# ------------------------------------------------------------
# 28. Analysis summary
# ------------------------------------------------------------

analysis_summary <- tibble::tibble(
  
  metric = c(
    "Dataset",
    "Comparison",
    "Statistical unit",
    "Normal pools",
    "OA pools",
    "Total probesets tested",
    "Total genes tested",
    "Gene-level FDR < 0.05",
    "Gene-level higher in OA",
    "Gene-level lower in OA",
    "Gene-level FDR and abs(log2FC) >= 1",
    "Gene-level TREAT FDR < 0.05",
    "MMP13 probe",
    "MMP13 log2FC",
    "MMP13 standard error",
    "MMP13 95% CI lower",
    "MMP13 95% CI upper",
    "MMP13 approximate fold change",
    "MMP13 P-value",
    "MMP13 FDR",
    "Maximum MMP13 leave-one-out shift",
    "MMP13 direction consistent in all LOO models"
  ),
  
  value = c(
    "GSE169077",
    "Late-stage OA versus normal cartilage",
    "RNA pool",
    5,
    6,
    nrow(probe_de_primary),
    nrow(gene_de_primary),
    sum(
      gene_de_primary$FDR < 0.05
    ),
    sum(
      gene_de_primary$FDR < 0.05 &
        gene_de_primary$logFC > 0
    ),
    sum(
      gene_de_primary$FDR < 0.05 &
        gene_de_primary$logFC < 0
    ),
    sum(
      gene_de_primary$FDR < 0.05 &
        abs(gene_de_primary$logFC) >= 1
    ),
    sum(
      gene_treat_primary$FDR < 0.05
    ),
    mmp13_probe_id,
    mmp13_primary_effect$effect,
    mmp13_primary_effect$standard_error,
    mmp13_primary_effect$CI_95_lower,
    mmp13_primary_effect$CI_95_upper,
    mmp13_primary_effect$approximate_fold_change,
    mmp13_primary_effect$PValue,
    mmp13_primary_effect$FDR,
    mmp13_influence_summary$
      maximum_absolute_shift,
    mmp13_influence_summary$
      direction_consistent_in_all_models
  )
)

cat("\nAnalysis summary:\n")

print(
  analysis_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 29. Menyimpan R objects
# ------------------------------------------------------------

saveRDS(
  design_primary,
  file =
    "data_processed/GSE169077_design_limma_primary.rds"
)

saveRDS(
  probe_fit_primary,
  file =
    "data_processed/GSE169077_probe_limma_fit.rds"
)

saveRDS(
  gene_fit_primary,
  file =
    "data_processed/GSE169077_gene_limma_fit.rds"
)

saveRDS(
  probe_de_primary,
  file =
    "data_processed/GSE169077_probe_DE_primary.rds"
)

saveRDS(
  gene_de_primary,
  file =
    "data_processed/GSE169077_gene_DE_primary.rds"
)

saveRDS(
  probe_treat_primary,
  file =
    "data_processed/GSE169077_probe_TREAT_lfc1.rds"
)

saveRDS(
  gene_treat_primary,
  file =
    "data_processed/GSE169077_gene_TREAT_lfc1.rds"
)

saveRDS(
  expression_gene,
  file =
    "data_processed/GSE169077_expression_gene_level.rds"
)

saveRDS(
  mmp13_primary_effect,
  file =
    "data_processed/GSE169077_MMP13_primary_effect.rds"
)

saveRDS(
  leave_one_out_results,
  file =
    "data_processed/GSE169077_MMP13_leave_one_out_results.rds"
)

saveRDS(
  mmp13_forest_candidate,
  file =
    "data_processed/GSE169077_MMP13_forest_candidate.rds"
)

saveRDS(
  gsea_ranked_genes,
  file =
    "data_processed/GSE169077_GSEA_ranked_genes.rds"
)

# ------------------------------------------------------------
# 30. Menyimpan full tables sebagai CSV
# ------------------------------------------------------------

utils::write.csv(
  probe_de_primary,
  file =
    "results/tables/GSE169077_probe_DE_primary.csv",
  row.names = FALSE
)

utils::write.csv(
  gene_de_primary,
  file =
    "results/tables/GSE169077_gene_DE_primary.csv",
  row.names = FALSE
)

utils::write.csv(
  gsea_ranked_genes,
  file =
    "results/tables/GSE169077_GSEA_ranked_genes.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 31. Menyimpan Excel workbook
# ------------------------------------------------------------

openxlsx::write.xlsx(
  list(
    
    Analysis_Summary =
      as.data.frame(
        analysis_summary
      ),
    
    DE_Summary =
      as.data.frame(
        de_summary
      ),
    
    MMP13_Primary_Effect =
      as.data.frame(
        mmp13_primary_effect
      ),
    
    MMP13_Direct_Summary =
      as.data.frame(
        mmp13_direct_summary
      ),
    
    MMP13_Direct_Check =
      as.data.frame(
        mmp13_direct_effect
      ),
    
    MMP13_Leave_One_Out =
      as.data.frame(
        leave_one_out_results
      ),
    
    MMP13_Influence_Summary =
      as.data.frame(
        mmp13_influence_summary
      ),
    
    MMP13_Forest_Candidate =
      as.data.frame(
        mmp13_forest_candidate
      ),
    
    Significant_Genes_FDR =
      as.data.frame(
        significant_gene_FDR
      ),
    
    Significant_Genes_FDR_LogFC1 =
      as.data.frame(
        significant_gene_FDR_logFC1
      ),
    
    Significant_Genes_TREAT =
      as.data.frame(
        significant_gene_TREAT
      ),
    
    Gene_TREAT_All =
      as.data.frame(
        gene_treat_primary
      ),
    
    GSEA_Ranked_Genes =
      as.data.frame(
        gsea_ranked_genes
      ),
    
    Top_100_Genes =
      as.data.frame(
        top_100_gene_results
      ),
    
    Top_100_Probes =
      as.data.frame(
        top_100_probe_results
      )
  ),
  
  file =
    "results/tables/GSE169077_limma_DE_results.xlsx",
  
  overwrite = TRUE
)

openxlsx::write.xlsx(
  as.data.frame(
    mmp13_forest_candidate
  ),
  
  file =
    "results/tables/GSE169077_MMP13_forest_candidate.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 32. Menyimpan session information
# ------------------------------------------------------------

sink(
  "results/tables/GSE169077_limma_DE_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 33. Pesan akhir
# ------------------------------------------------------------

message("")
message("================================================")
message("GSE169077 LIMMA DIFFERENTIAL EXPRESSION SELESAI")
message("================================================")

message(
  "Total probesets tested          : ",
  nrow(probe_de_primary)
)

message(
  "Total genes tested              : ",
  nrow(gene_de_primary)
)

message(
  "Gene-level FDR < 0.05           : ",
  sum(
    gene_de_primary$FDR < 0.05
  )
)

message(
  "Higher in OA                    : ",
  sum(
    gene_de_primary$FDR < 0.05 &
      gene_de_primary$logFC > 0
  )
)

message(
  "Lower in OA                     : ",
  sum(
    gene_de_primary$FDR < 0.05 &
      gene_de_primary$logFC < 0
  )
)

message(
  "FDR and |log2FC| >= 1           : ",
  sum(
    gene_de_primary$FDR < 0.05 &
      abs(gene_de_primary$logFC) >= 1
  )
)

message(
  "TREAT FDR < 0.05                : ",
  sum(
    gene_treat_primary$FDR < 0.05
  )
)

message(
  "MMP13 probe                     : ",
  mmp13_probe_id
)

message(
  "MMP13 log2FC                    : ",
  round(
    mmp13_primary_effect$effect,
    4
  )
)

message(
  "MMP13 standard error            : ",
  round(
    mmp13_primary_effect$standard_error,
    4
  )
)

message(
  "MMP13 95% CI                    : ",
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
  "MMP13 approximate fold change   : ",
  round(
    mmp13_primary_effect$
      approximate_fold_change,
    3
  )
)

message(
  "MMP13 P-value                   : ",
  format(
    mmp13_primary_effect$PValue,
    scientific = TRUE,
    digits = 4
  )
)

message(
  "MMP13 FDR                       : ",
  format(
    mmp13_primary_effect$FDR,
    scientific = TRUE,
    digits = 4
  )
)

message(
  "Maximum leave-one-out shift     : ",
  round(
    mmp13_influence_summary$
      maximum_absolute_shift,
    4
  )
)

message(
  "Direction consistent in all LOO : ",
  mmp13_influence_summary$
    direction_consistent_in_all_models
)

message(
  "Statistical unit                : RNA pool"
)

message("================================================")