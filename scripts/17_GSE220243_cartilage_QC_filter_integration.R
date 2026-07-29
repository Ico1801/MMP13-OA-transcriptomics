# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 17_GSE220243_cartilage_QC_filter_integration.R
# PURPOSE :
#   1. Donor-specific adaptive QC filtering
#   2. Per-donor doublet detection with scDblFinder
#   3. SCTransform v2 normalization
#   4. SCT-RPCA donor integration
#   5. PCA, UMAP, clustering, and integration diagnostics
#
# INPUT:
#   data_processed/
#   GSE220243_cartilage_raw_seurat_list.rds
#
# OUTPUT:
#   data_processed/
#   GSE220243_cartilage_integrated_seurat.rds
#
# IMPORTANT:
#   - Integrated representation is used for visualization
#     and clustering.
#   - RNA expression is retained for biological expression
#     analysis.
#   - Differential expression will not be performed using
#     the integrated assay.
# ============================================================

rm(list = ls())
gc()

set.seed(20260722)

# ------------------------------------------------------------
# 1. Package installation and loading
# ------------------------------------------------------------

cran_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "sctransform",
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "patchwork",
  "openxlsx"
)

bioconductor_packages <- c(
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scDblFinder",
  "BiocParallel",
  "glmGamPoi"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_cran) > 0) {
  
  message(
    "Menginstal paket CRAN: ",
    paste(
      missing_cran,
      collapse = ", "
    )
  )
  
  install.packages(
    missing_cran,
    repos = "https://cloud.r-project.org"
  )
}

missing_bioconductor <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioconductor) > 0) {
  
  if (
    !requireNamespace(
      "BiocManager",
      quietly = TRUE
    )
  ) {
    
    install.packages(
      "BiocManager",
      repos = "https://cloud.r-project.org"
    )
  }
  
  BiocManager::install(
    missing_bioconductor,
    ask = FALSE,
    update = FALSE
  )
}

required_packages <- c(
  cran_packages,
  bioconductor_packages
)

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

library(Seurat)
library(SeuratObject)
library(Matrix)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(patchwork)
library(openxlsx)

message(
  "Semua paket QC dan integration berhasil dimuat."
)

# ------------------------------------------------------------
# 2. Package versions
# ------------------------------------------------------------

package_versions <- tibble::tibble(
  
  package = c(
    "Seurat",
    "SeuratObject",
    "sctransform",
    "scDblFinder",
    "SingleCellExperiment",
    "glmGamPoi"
  ),
  
  version = c(
    as.character(
      packageVersion("Seurat")
    ),
    
    as.character(
      packageVersion("SeuratObject")
    ),
    
    as.character(
      packageVersion("sctransform")
    ),
    
    as.character(
      packageVersion("scDblFinder")
    ),
    
    as.character(
      packageVersion("SingleCellExperiment")
    ),
    
    as.character(
      packageVersion("glmGamPoi")
    )
  )
)

cat("\nVersi paket:\n")

print(
  package_versions,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 3. Input and output folders
# ------------------------------------------------------------

input_file <-
  paste0(
    "data_processed/",
    "GSE220243_cartilage_raw_seurat_list.rds"
  )

processed_folder <-
  "data_processed"

table_folder <-
  "results/tables"

supplementary_figure_folder <-
  paste0(
    "results/figures/",
    "Figure 1 Suplementary/",
    "S13"
  )

required_folders <- c(
  processed_folder,
  table_folder,
  supplementary_figure_folder
)

for (current_folder in required_folders) {
  
  dir.create(
    current_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (!file.exists(input_file)) {
  
  stop(
    "Input Seurat list tidak ditemukan:\n",
    input_file,
    "\nJalankan script 16 terlebih dahulu."
  )
}

# ------------------------------------------------------------
# 4. Load raw donor-level Seurat objects
# ------------------------------------------------------------

raw_seurat_list <- readRDS(
  input_file
)

if (!is.list(raw_seurat_list)) {
  
  stop(
    "Input bukan list of Seurat objects."
  )
}

cat("\nJumlah objek donor:\n")

print(
  length(
    raw_seurat_list
  )
)

if (
  length(
    raw_seurat_list
  ) != 12
) {
  
  stop(
    "Jumlah donor bukan 12. Ditemukan: ",
    length(
      raw_seurat_list
    )
  )
}

valid_seurat_objects <- vapply(
  
  raw_seurat_list,
  
  inherits,
  
  logical(1),
  
  what =
    "Seurat"
)

if (!all(valid_seurat_objects)) {
  
  stop(
    "Satu atau lebih elemen bukan objek Seurat."
  )
}

# ------------------------------------------------------------
# 5. Validate donor metadata and reorder objects
# ------------------------------------------------------------

required_metadata_columns <- c(
  "sample_label",
  "group",
  "donor",
  "geo_accession",
  "percent.mt"
)

sample_labels_from_objects <- vapply(
  
  raw_seurat_list,
  
  function(current_object) {
    
    current_metadata <-
      current_object[[]]
    
    missing_metadata <- setdiff(
      required_metadata_columns,
      colnames(
        current_metadata
      )
    )
    
    if (length(missing_metadata) > 0) {
      
      stop(
        "Metadata berikut tidak ditemukan: ",
        paste(
          missing_metadata,
          collapse = ", "
        )
      )
    }
    
    current_labels <- unique(
      as.character(
        current_metadata$sample_label
      )
    )
    
    if (length(current_labels) != 1) {
      
      stop(
        "Satu objek berisi lebih dari satu sample_label."
      )
    }
    
    current_labels[1]
  },
  
  FUN.VALUE =
    character(1)
)

if (
  anyDuplicated(
    sample_labels_from_objects
  ) > 0
) {
  
  stop(
    "Terdapat sample_label duplikat."
  )
}

names(raw_seurat_list) <-
  sample_labels_from_objects

expected_sample_order <- c(
  "Normal1",
  "Normal2",
  "Normal3",
  "Normal4",
  "Normal5",
  "Normal6",
  "OA1",
  "OA2",
  "OA3",
  "OA4",
  "OA5",
  "OA6"
)

missing_expected_samples <- setdiff(
  expected_sample_order,
  names(
    raw_seurat_list
  )
)

if (length(missing_expected_samples) > 0) {
  
  stop(
    "Sampel berikut tidak ditemukan: ",
    paste(
      missing_expected_samples,
      collapse = ", "
    )
  )
}

raw_seurat_list <-
  raw_seurat_list[
    expected_sample_order
  ]

# ------------------------------------------------------------
# 6. Build sample design table
# ------------------------------------------------------------

sample_design <- dplyr::bind_rows(
  
  lapply(
    
    names(
      raw_seurat_list
    ),
    
    function(current_sample) {
      
      current_object <-
        raw_seurat_list[[current_sample]]
      
      current_metadata <-
        current_object[[]]
      
      tibble::tibble(
        
        sample_label =
          current_sample,
        
        group =
          as.character(
            unique(
              current_metadata$group
            )[1]
          ),
        
        donor =
          as.character(
            unique(
              current_metadata$donor
            )[1]
          ),
        
        geo_accession =
          as.character(
            unique(
              current_metadata$geo_accession
            )[1]
          ),
        
        raw_cells =
          ncol(
            current_object
          ),
        
        genes =
          nrow(
            current_object
          )
      )
    }
  )
) %>%
  
  dplyr::mutate(
    
    sample_label = factor(
      sample_label,
      levels =
        expected_sample_order
    ),
    
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    )
  ) %>%
  
  dplyr::arrange(
    sample_label
  )

cat("\nSample design:\n")

print(
  sample_design,
  n = Inf,
  width = Inf
)

stopifnot(
  sum(sample_design$group == "Normal") == 6
)

stopifnot(
  sum(sample_design$group == "OA") == 6
)

# ------------------------------------------------------------
# 7. Robust spread function
# ------------------------------------------------------------

calculate_robust_spread <- function(
    input_values
) {
  
  input_values <- input_values[
    is.finite(
      input_values
    )
  ]
  
  if (length(input_values) < 2) {
    
    return(
      1
    )
  }
  
  current_spread <- stats::mad(
    input_values,
    constant = 1.4826,
    na.rm = TRUE
  )
  
  if (
    !is.finite(current_spread) ||
    current_spread <= 0
  ) {
    
    current_spread <-
      stats::IQR(
        input_values,
        na.rm = TRUE
      ) /
      1.349
  }
  
  if (
    !is.finite(current_spread) ||
    current_spread <= 0
  ) {
    
    current_spread <-
      stats::sd(
        input_values,
        na.rm = TRUE
      )
  }
  
  if (
    !is.finite(current_spread) ||
    current_spread <= 0
  ) {
    
    current_spread <-
      1
  }
  
  current_spread
}

# ------------------------------------------------------------
# 8. Function to calculate donor-specific QC thresholds
#
# Lower feature threshold:
#   median - 3 robust SD, minimum 200 genes
#
# Upper feature/count thresholds:
#   median + 4 robust SD
#
# Mitochondrial threshold:
#   median + 3 robust SD
#   with a permissive range from 10% to 30%
# ------------------------------------------------------------

calculate_qc_thresholds <- function(
    current_metadata,
    current_sample,
    current_group
) {
  
  feature_values <-
    current_metadata$nFeature_RNA
  
  count_values <-
    current_metadata$nCount_RNA
  
  mitochondrial_values <-
    current_metadata$percent.mt
  
  feature_median <- stats::median(
    feature_values,
    na.rm = TRUE
  )
  
  feature_spread <- calculate_robust_spread(
    feature_values
  )
  
  count_median <- stats::median(
    count_values,
    na.rm = TRUE
  )
  
  count_spread <- calculate_robust_spread(
    count_values
  )
  
  mt_median <- stats::median(
    mitochondrial_values,
    na.rm = TRUE
  )
  
  mt_spread <- calculate_robust_spread(
    mitochondrial_values
  )
  
  lower_feature_threshold <- max(
    200,
    floor(
      feature_median -
        3 *
        feature_spread
    )
  )
  
  upper_feature_threshold <- ceiling(
    feature_median +
      4 *
      feature_spread
  )
  
  upper_feature_threshold <- max(
    upper_feature_threshold,
    lower_feature_threshold + 1
  )
  
  upper_count_threshold <- ceiling(
    count_median +
      4 *
      count_spread
  )
  
  raw_mt_threshold <-
    mt_median +
    3 *
    mt_spread
  
  upper_mt_threshold <- min(
    30,
    max(
      10,
      raw_mt_threshold
    )
  )
  
  tibble::tibble(
    
    sample_label =
      current_sample,
    
    group =
      current_group,
    
    cells_before_QC =
      nrow(
        current_metadata
      ),
    
    median_nFeature_RNA =
      feature_median,
    
    robust_spread_nFeature =
      feature_spread,
    
    lower_nFeature_RNA =
      lower_feature_threshold,
    
    upper_nFeature_RNA =
      upper_feature_threshold,
    
    median_nCount_RNA =
      count_median,
    
    robust_spread_nCount =
      count_spread,
    
    upper_nCount_RNA =
      upper_count_threshold,
    
    median_percent_mt =
      mt_median,
    
    robust_spread_percent_mt =
      mt_spread,
    
    raw_adaptive_mt_threshold =
      raw_mt_threshold,
    
    upper_percent_mt =
      upper_mt_threshold
  )
}

# ------------------------------------------------------------
# 9. Apply adaptive QC separately to each donor
# ------------------------------------------------------------

basic_qc_list <- list()
qc_threshold_list <- list()
qc_flag_list <- list()
qc_failure_summary_list <- list()

for (
  current_sample in
  names(
    raw_seurat_list
  )
) {
  
  message("")
  message(
    "============================================"
  )
  
  message(
    "Adaptive QC: ",
    current_sample
  )
  
  message(
    "============================================"
  )
  
  current_object <-
    raw_seurat_list[[current_sample]]
  
  Seurat::DefaultAssay(
    current_object
  ) <- "RNA"
  
  current_metadata <-
    current_object[[]] %>%
    
    tibble::rownames_to_column(
      "cell_id"
    )
  
  current_group <- as.character(
    unique(
      current_metadata$group
    )[1]
  )
  
  current_thresholds <-
    calculate_qc_thresholds(
      
      current_metadata =
        current_metadata,
      
      current_sample =
        current_sample,
      
      current_group =
        current_group
    )
  
  qc_threshold_list[[current_sample]] <-
    current_thresholds
  
  current_flags <- current_metadata %>%
    
    dplyr::mutate(
      
      sample_label =
        current_sample,
      
      pass_lower_features =
        nFeature_RNA >=
        current_thresholds$
        lower_nFeature_RNA,
      
      pass_upper_features =
        nFeature_RNA <=
        current_thresholds$
        upper_nFeature_RNA,
      
      pass_upper_counts =
        nCount_RNA <=
        current_thresholds$
        upper_nCount_RNA,
      
      pass_mitochondrial =
        percent.mt <=
        current_thresholds$
        upper_percent_mt,
      
      pass_basic_qc =
        pass_lower_features &
        pass_upper_features &
        pass_upper_counts &
        pass_mitochondrial,
      
      failed_criteria_number =
        rowSums(
          cbind(
            !pass_lower_features,
            !pass_upper_features,
            !pass_upper_counts,
            !pass_mitochondrial
          )
        )
    )
  
  current_object$pass_lower_features <-
    current_flags$
    pass_lower_features[
      match(
        colnames(
          current_object
        ),
        current_flags$cell_id
      )
    ]
  
  current_object$pass_upper_features <-
    current_flags$
    pass_upper_features[
      match(
        colnames(
          current_object
        ),
        current_flags$cell_id
      )
    ]
  
  current_object$pass_upper_counts <-
    current_flags$
    pass_upper_counts[
      match(
        colnames(
          current_object
        ),
        current_flags$cell_id
      )
    ]
  
  current_object$pass_mitochondrial <-
    current_flags$
    pass_mitochondrial[
      match(
        colnames(
          current_object
        ),
        current_flags$cell_id
      )
    ]
  
  current_object$pass_basic_qc <-
    current_flags$
    pass_basic_qc[
      match(
        colnames(
          current_object
        ),
        current_flags$cell_id
      )
    ]
  
  retained_cell_ids <- current_flags %>%
    
    dplyr::filter(
      pass_basic_qc
    ) %>%
    
    dplyr::pull(
      cell_id
    )
  
  current_filtered_object <- subset(
    
    x =
      current_object,
    
    cells =
      retained_cell_ids
  )
  
  retention_percentage <-
    100 *
    ncol(
      current_filtered_object
    ) /
    ncol(
      current_object
    )
  
  message(
    "Cells before adaptive QC : ",
    ncol(
      current_object
    )
  )
  
  message(
    "Cells after adaptive QC  : ",
    ncol(
      current_filtered_object
    )
  )
  
  message(
    "Retention                : ",
    round(
      retention_percentage,
      2
    ),
    "%"
  )
  
  if (
    ncol(
      current_filtered_object
    ) < 200
  ) {
    
    stop(
      current_sample,
      " memiliki kurang dari 200 sel setelah basic QC."
    )
  }
  
  if (
    retention_percentage < 50
  ) {
    
    warning(
      current_sample,
      ": kurang dari 50% sel dipertahankan."
    )
  }
  
  qc_failure_summary_list[[current_sample]] <-
    tibble::tibble(
      
      sample_label =
        current_sample,
      
      group =
        current_group,
      
      cells_before =
        ncol(
          current_object
        ),
      
      failed_low_features =
        sum(
          !current_flags$
            pass_lower_features
        ),
      
      failed_high_features =
        sum(
          !current_flags$
            pass_upper_features
        ),
      
      failed_high_counts =
        sum(
          !current_flags$
            pass_upper_counts
        ),
      
      failed_high_mitochondrial =
        sum(
          !current_flags$
            pass_mitochondrial
        ),
      
      cells_after_basic_QC =
        ncol(
          current_filtered_object
        ),
      
      basic_QC_retention_percent =
        retention_percentage
    )
  
  qc_flag_list[[current_sample]] <-
    current_flags
  
  basic_qc_list[[current_sample]] <-
    current_filtered_object
  
  rm(
    current_object,
    current_filtered_object,
    current_metadata,
    current_flags
  )
  
  gc()
}

qc_threshold_table <- dplyr::bind_rows(
  qc_threshold_list
)

all_cell_qc_flags <- dplyr::bind_rows(
  qc_flag_list
)

qc_failure_summary <- dplyr::bind_rows(
  qc_failure_summary_list
)

cat("\nDonor-specific QC thresholds:\n")

print(
  qc_threshold_table,
  n = Inf,
  width = Inf
)

cat("\nBasic QC failure summary:\n")

print(
  qc_failure_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 10. Save objects after adaptive basic QC
# ------------------------------------------------------------

saveRDS(
  
  basic_qc_list,
  
  file = file.path(
    processed_folder,
    "GSE220243_cartilage_basic_QC_seurat_list.rds"
  ),
  
  compress =
    "gzip"
)

message(
  "Basic-QC Seurat list berhasil disimpan."
)

rm(
  raw_seurat_list
)

gc()

# ------------------------------------------------------------
# 11. Per-donor scDblFinder doublet detection
# ------------------------------------------------------------

singlet_seurat_list <- list()
doublet_annotated_list <- list()
doublet_summary_list <- list()
doublet_score_list <- list()

for (
  current_index in
  seq_along(
    basic_qc_list
  )
) {
  
  current_sample <- names(
    basic_qc_list
  )[current_index]
  
  message("")
  message(
    "============================================"
  )
  
  message(
    "scDblFinder: ",
    current_sample
  )
  
  message(
    "============================================"
  )
  
  current_object <-
    basic_qc_list[[current_sample]]
  
  Seurat::DefaultAssay(
    current_object
  ) <- "RNA"
  
  current_sce <- Seurat::as.SingleCellExperiment(
    
    x =
      current_object,
    
    assay =
      "RNA"
  )
  
  stopifnot(
    identical(
      colnames(
        current_sce
      ),
      colnames(
        current_object
      )
    )
  )
  
  current_sce$sample_label <-
    current_sample
  
  set.seed(
    20260722 +
      current_index
  )
  
  current_sce <- tryCatch(
    
    scDblFinder::scDblFinder(
      
      sce =
        current_sce,
      
      clusters =
        TRUE,
      
      samples =
        NULL,
      
      BPPARAM =
        BiocParallel::SerialParam(
          progressbar = TRUE
        ),
      
      verbose =
        TRUE
    ),
    
    error = function(error_object) {
      
      stop(
        "scDblFinder gagal untuk ",
        current_sample,
        ":\n",
        conditionMessage(
          error_object
        )
      )
    }
  )
  
  current_doublet_metadata <-
    as.data.frame(
      SummarizedExperiment::colData(
        current_sce
      )
    )
  
  required_doublet_columns <- c(
    "scDblFinder.score",
    "scDblFinder.class"
  )
  
  missing_doublet_columns <- setdiff(
    required_doublet_columns,
    colnames(
      current_doublet_metadata
    )
  )
  
  if (length(missing_doublet_columns) > 0) {
    
    stop(
      "Kolom scDblFinder tidak ditemukan untuk ",
      current_sample,
      ": ",
      paste(
        missing_doublet_columns,
        collapse = ", "
      )
    )
  }
  
  stopifnot(
    all(
      colnames(
        current_object
      ) %in%
        rownames(
          current_doublet_metadata
        )
    )
  )
  
  current_object$scDblFinder.score <-
    current_doublet_metadata[
      colnames(
        current_object
      ),
      "scDblFinder.score"
    ]
  
  current_object$scDblFinder.class <-
    as.character(
      current_doublet_metadata[
        colnames(
          current_object
        ),
        "scDblFinder.class"
      ]
    )
  
  current_doublet_scores <- tibble::tibble(
    
    cell_id =
      colnames(
        current_object
      ),
    
    sample_label =
      current_sample,
    
    group =
      as.character(
        current_object$group
      ),
    
    scDblFinder_score =
      current_object$
      scDblFinder.score,
    
    scDblFinder_class =
      current_object$
      scDblFinder.class
  )
  
  current_singlet_cells <- current_doublet_scores %>%
    
    dplyr::filter(
      scDblFinder_class ==
        "singlet"
    ) %>%
    
    dplyr::pull(
      cell_id
    )
  
  current_singlet_object <- subset(
    
    x =
      current_object,
    
    cells =
      current_singlet_cells
  )
  
  predicted_doublets <- sum(
    current_object$
      scDblFinder.class ==
      "doublet"
  )
  
  predicted_singlets <- sum(
    current_object$
      scDblFinder.class ==
      "singlet"
  )
  
  predicted_doublet_rate <-
    100 *
    predicted_doublets /
    ncol(
      current_object
    )
  
  message(
    "Cells before doublet removal : ",
    ncol(
      current_object
    )
  )
  
  message(
    "Predicted doublets           : ",
    predicted_doublets
  )
  
  message(
    "Predicted doublet rate       : ",
    round(
      predicted_doublet_rate,
      2
    ),
    "%"
  )
  
  message(
    "Retained singlets            : ",
    ncol(
      current_singlet_object
    )
  )
  
  if (
    ncol(
      current_singlet_object
    ) < 150
  ) {
    
    stop(
      current_sample,
      " memiliki kurang dari 150 singlets."
    )
  }
  
  doublet_summary_list[[current_sample]] <-
    tibble::tibble(
      
      sample_label =
        current_sample,
      
      group =
        as.character(
          unique(
            current_object$group
          )[1]
        ),
      
      cells_after_basic_QC =
        ncol(
          current_object
        ),
      
      predicted_singlets =
        predicted_singlets,
      
      predicted_doublets =
        predicted_doublets,
      
      predicted_doublet_percent =
        predicted_doublet_rate,
      
      cells_after_doublet_removal =
        ncol(
          current_singlet_object
        )
    )
  
  doublet_score_list[[current_sample]] <-
    current_doublet_scores
  
  doublet_annotated_list[[current_sample]] <-
    current_object
  
  singlet_seurat_list[[current_sample]] <-
    current_singlet_object
  
  rm(
    current_object,
    current_singlet_object,
    current_sce,
    current_doublet_metadata
  )
  
  gc()
}

doublet_summary <- dplyr::bind_rows(
  doublet_summary_list
)

all_doublet_scores <- dplyr::bind_rows(
  doublet_score_list
)

cat("\nDoublet summary:\n")

print(
  doublet_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 12. Save doublet-annotated and singlet lists
# ------------------------------------------------------------

saveRDS(
  
  doublet_annotated_list,
  
  file = file.path(
    processed_folder,
    paste0(
      "GSE220243_cartilage_",
      "doublet_annotated_seurat_list.rds"
    )
  ),
  
  compress =
    "gzip"
)

saveRDS(
  
  singlet_seurat_list,
  
  file = file.path(
    processed_folder,
    "GSE220243_cartilage_filtered_singlet_list.rds"
  ),
  
  compress =
    "gzip"
)

message(
  "Filtered singlet Seurat list berhasil disimpan."
)

rm(
  basic_qc_list,
  doublet_annotated_list
)

gc()

# ------------------------------------------------------------
# 13. Final retention summary
# ------------------------------------------------------------

retention_summary <- sample_design %>%
  
  dplyr::select(
    sample_label,
    group,
    raw_cells
  ) %>%
  
  dplyr::left_join(
    
    qc_failure_summary %>%
      
      dplyr::select(
        sample_label,
        cells_after_basic_QC
      ),
    
    by =
      "sample_label"
  ) %>%
  
  dplyr::left_join(
    
    doublet_summary %>%
      
      dplyr::select(
        sample_label,
        predicted_doublets,
        cells_after_doublet_removal
      ),
    
    by =
      "sample_label"
  ) %>%
  
  dplyr::mutate(
    
    retained_after_basic_QC_percent =
      100 *
      cells_after_basic_QC /
      raw_cells,
    
    retained_final_percent =
      100 *
      cells_after_doublet_removal /
      raw_cells
  )

cat("\nFinal retention summary:\n")

print(
  retention_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 14. MMP13 detection after QC and doublet removal
# ------------------------------------------------------------

mmp13_post_QC_list <- list()

for (
  current_sample in
  names(
    singlet_seurat_list
  )
) {
  
  current_object <-
    singlet_seurat_list[[current_sample]]
  
  current_group <- as.character(
    unique(
      current_object$group
    )[1]
  )
  
  mmp13_present <-
    "MMP13" %in%
    rownames(
      current_object[["RNA"]]
    )
  
  if (mmp13_present) {
    
    mmp13_count_matrix <-
      SeuratObject::LayerData(
        
        object =
          current_object,
        
        assay =
          "RNA",
        
        layer =
          "counts",
        
        features =
          "MMP13"
      )
    
    mmp13_counts <- as.numeric(
      mmp13_count_matrix[
        1,
        ,
        drop = TRUE
      ]
    )
    
    positive_cells <- sum(
      mmp13_counts > 0
    )
    
    total_counts <- sum(
      mmp13_counts
    )
    
  } else {
    
    positive_cells <- 0
    total_counts <- 0
  }
  
  detection_percent <-
    100 *
    positive_cells /
    ncol(
      current_object
    )
  
  mmp13_post_QC_list[[current_sample]] <-
    tibble::tibble(
      
      sample_label =
        current_sample,
      
      group =
        current_group,
      
      filtered_singlet_cells =
        ncol(
          current_object
        ),
      
      MMP13_present =
        mmp13_present,
      
      MMP13_positive_cells =
        positive_cells,
      
      MMP13_detection_percent =
        detection_percent,
      
      MMP13_total_counts =
        total_counts
    )
}

mmp13_post_QC_summary <- dplyr::bind_rows(
  mmp13_post_QC_list
)

cat("\nMMP13 detection after QC:\n")

print(
  mmp13_post_QC_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 15. SCTransform v2 normalization per donor
# ------------------------------------------------------------

sct_seurat_list <- list()

for (
  current_index in
  seq_along(
    singlet_seurat_list
  )
) {
  
  current_sample <- names(
    singlet_seurat_list
  )[current_index]
  
  message("")
  message(
    "============================================"
  )
  
  message(
    "SCTransform v2: ",
    current_sample
  )
  
  message(
    "============================================"
  )
  
  current_object <-
    singlet_seurat_list[[current_sample]]
  
  Seurat::DefaultAssay(
    current_object
  ) <- "RNA"
  
  set.seed(
    20260722 +
      100 +
      current_index
  )
  
  current_object <- Seurat::SCTransform(
    
    object =
      current_object,
    
    assay =
      "RNA",
    
    new.assay.name =
      "SCT",
    
    vars.to.regress =
      "percent.mt",
    
    vst.flavor =
      "v2",
    
    variable.features.n =
      3000,
    
    conserve.memory =
      TRUE,
    
    return.only.var.genes =
      TRUE,
    
    seed.use =
      20260722 +
      current_index,
    
    verbose =
      TRUE
  )
  
  sct_seurat_list[[current_sample]] <-
    current_object
  
  rm(
    current_object
  )
  
  gc()
}

message(
  "SCTransform selesai untuk seluruh donor."
)

# ------------------------------------------------------------
# 16. Select SCT integration features
# ------------------------------------------------------------

integration_features <- Seurat::SelectIntegrationFeatures(
  
  object.list =
    sct_seurat_list,
  
  nfeatures =
    3000,
  
  assay =
    rep(
      "SCT",
      length(
        sct_seurat_list
      )
    ),
  
  verbose =
    TRUE
)

cat("\nJumlah integration features:\n")

print(
  length(
    integration_features
  )
)

cat("\nMMP13 termasuk integration features:\n")

print(
  "MMP13" %in%
    integration_features
)

# ------------------------------------------------------------
# 17. Prepare SCT residuals for integration
# ------------------------------------------------------------

sct_seurat_list <- Seurat::PrepSCTIntegration(
  
  object.list =
    sct_seurat_list,
  
  anchor.features =
    integration_features,
  
  assay =
    rep(
      "SCT",
      length(
        sct_seurat_list
      )
    ),
  
  verbose =
    TRUE
)

# ------------------------------------------------------------
# 18. Determine PCA and anchor parameters
# ------------------------------------------------------------

cells_per_sct_object <- vapply(
  
  sct_seurat_list,
  
  ncol,
  
  integer(1)
)

minimum_cells_in_donor <- min(
  cells_per_sct_object
)

preintegration_npcs <- min(
  50,
  minimum_cells_in_donor - 1,
  length(
    integration_features
  ) - 1
)

if (preintegration_npcs < 20) {
  
  stop(
    "Terlalu sedikit sel/features untuk PCA integration."
  )
}

anchor_dimension_number <- min(
  30,
  preintegration_npcs
)

anchor_dimensions <- seq_len(
  anchor_dimension_number
)

k_filter_use <- min(
  200,
  minimum_cells_in_donor - 1
)

k_score_use <- min(
  30,
  minimum_cells_in_donor - 1
)

k_weight_use <- min(
  100,
  minimum_cells_in_donor - 1
)

integration_parameter_summary <- tibble::tibble(
  
  parameter = c(
    "Minimum donor cell number",
    "Integration features",
    "Pre-integration PCA components",
    "Anchor dimensions",
    "k.filter",
    "k.score",
    "k.weight"
  ),
  
  value = c(
    minimum_cells_in_donor,
    length(
      integration_features
    ),
    preintegration_npcs,
    anchor_dimension_number,
    k_filter_use,
    k_score_use,
    k_weight_use
  )
)

cat("\nIntegration parameters:\n")

print(
  integration_parameter_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 19. Run donor-level PCA for RPCA integration
# ------------------------------------------------------------

for (
  current_index in
  seq_along(
    sct_seurat_list
  )
) {
  
  current_sample <- names(
    sct_seurat_list
  )[current_index]
  
  message(
    "Running donor PCA: ",
    current_sample
  )
  
  sct_seurat_list[[current_index]] <-
    Seurat::RunPCA(
      
      object =
        sct_seurat_list[[current_index]],
      
      assay =
        "SCT",
      
      features =
        integration_features,
      
      npcs =
        preintegration_npcs,
      
      verbose =
        FALSE,
      
      seed.use =
        20260722 +
        200 +
        current_index
    )
}

# ------------------------------------------------------------
# 20. Select balanced reference donors
#
# One Normal and one OA donor with the largest retained
# cell number are used as balanced integration references.
# ------------------------------------------------------------

reference_candidate_table <- tibble::tibble(
  
  list_index =
    seq_along(
      sct_seurat_list
    ),
  
  sample_label =
    names(
      sct_seurat_list
    ),
  
  group =
    vapply(
      
      sct_seurat_list,
      
      function(current_object) {
        
        as.character(
          unique(
            current_object$group
          )[1]
        )
      },
      
      FUN.VALUE =
        character(1)
    ),
  
  cells =
    cells_per_sct_object
)

normal_reference_index <- reference_candidate_table %>%
  
  dplyr::filter(
    group ==
      "Normal"
  ) %>%
  
  dplyr::slice_max(
    order_by =
      cells,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  dplyr::pull(
    list_index
  )

oa_reference_index <- reference_candidate_table %>%
  
  dplyr::filter(
    group ==
      "OA"
  ) %>%
  
  dplyr::slice_max(
    order_by =
      cells,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  dplyr::pull(
    list_index
  )

reference_indices <- c(
  normal_reference_index,
  oa_reference_index
)

reference_donor_table <- reference_candidate_table %>%
  
  dplyr::filter(
    list_index %in%
      reference_indices
  ) %>%
  
  dplyr::mutate(
    selected_as_reference =
      TRUE
  )

cat("\nBalanced reference donors:\n")

print(
  reference_donor_table,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 21. Find SCT-RPCA integration anchors
# ------------------------------------------------------------

message(
  "Mencari SCT-RPCA integration anchors..."
)

integration_anchors <- Seurat::FindIntegrationAnchors(
  
  object.list =
    sct_seurat_list,
  
  assay =
    rep(
      "SCT",
      length(
        sct_seurat_list
      )
    ),
  
  reference =
    reference_indices,
  
  normalization.method =
    "SCT",
  
  anchor.features =
    integration_features,
  
  reduction =
    "rpca",
  
  dims =
    anchor_dimensions,
  
  k.anchor =
    5,
  
  k.filter =
    k_filter_use,
  
  k.score =
    k_score_use,
  
  verbose =
    TRUE
)

message(
  "Integration anchors berhasil dibuat."
)

# ------------------------------------------------------------
# 22. Integrate SCT data
# ------------------------------------------------------------

message(
  "Mengintegrasikan 12 cartilage donors..."
)

integrated_object <- Seurat::IntegrateData(
  
  anchorset =
    integration_anchors,
  
  new.assay.name =
    "integrated",
  
  normalization.method =
    "SCT",
  
  dims =
    anchor_dimensions,
  
  k.weight =
    k_weight_use,
  
  preserve.order =
    TRUE,
  
  verbose =
    TRUE
)

message(
  "SCT-RPCA integration selesai."
)

# ------------------------------------------------------------
# 23. Final PCA on integrated assay
# ------------------------------------------------------------

Seurat::DefaultAssay(
  integrated_object
) <- "integrated"

final_npcs <- min(
  50,
  ncol(
    integrated_object
  ) - 1,
  nrow(
    integrated_object[["integrated"]]
  ) - 1
)

set.seed(
  20260722
)

integrated_object <- Seurat::RunPCA(
  
  object =
    integrated_object,
  
  assay =
    "integrated",
  
  npcs =
    final_npcs,
  
  verbose =
    TRUE,
  
  seed.use =
    20260722
)

# ------------------------------------------------------------
# 24. Select number of PCs using cumulative variance
# ------------------------------------------------------------

pca_standard_deviation <- Seurat::Stdev(
  
  object =
    integrated_object,
  
  reduction =
    "pca"
)

pca_variance_percent <-
  100 *
  pca_standard_deviation^2 /
  sum(
    pca_standard_deviation^2
  )

pca_cumulative_variance <- cumsum(
  pca_variance_percent
)

first_pc_reaching_80_percent <- which(
  pca_cumulative_variance >=
    80
)[1]

if (
  length(
    first_pc_reaching_80_percent
  ) == 0 ||
  is.na(
    first_pc_reaching_80_percent
  )
) {
  
  first_pc_reaching_80_percent <-
    min(
      30,
      length(
        pca_standard_deviation
      )
    )
}

pcs_for_clustering <- min(
  30,
  max(
    15,
    first_pc_reaching_80_percent
  ),
  length(
    pca_standard_deviation
  )
)

pca_variance_table <- tibble::tibble(
  
  PC =
    seq_along(
      pca_standard_deviation
    ),
  
  standard_deviation =
    pca_standard_deviation,
  
  variance_percent =
    pca_variance_percent,
  
  cumulative_variance_percent =
    pca_cumulative_variance,
  
  selected_for_clustering =
    PC <=
    pcs_for_clustering
)

cat("\nPCs selected for clustering:\n")

print(
  pcs_for_clustering
)

cat("\nPCA variance table:\n")

print(
  pca_variance_table,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 25. Neighbors, clustering, and UMAP
# ------------------------------------------------------------

selected_dimensions <- seq_len(
  pcs_for_clustering
)

set.seed(
  20260722
)

integrated_object <- Seurat::FindNeighbors(
  
  object =
    integrated_object,
  
  reduction =
    "pca",
  
  dims =
    selected_dimensions,
  
  verbose =
    TRUE
)

cluster_resolutions <- c(
  0.2,
  0.4,
  0.6,
  0.8
)

integrated_object <- Seurat::FindClusters(
  
  object =
    integrated_object,
  
  resolution =
    cluster_resolutions,
  
  algorithm =
    1,
  
  random.seed =
    20260722,
  
  verbose =
    TRUE
)

integrated_object <- Seurat::RunUMAP(
  
  object =
    integrated_object,
  
  reduction =
    "pca",
  
  dims =
    selected_dimensions,
  
  reduction.name =
    "umap",
  
  reduction.key =
    "UMAP_",
  
  seed.use =
    20260722,
  
  verbose =
    TRUE
)

# ------------------------------------------------------------
# 26. Select resolution 0.6 as working cluster identity
# ------------------------------------------------------------

metadata_column_names <- colnames(
  integrated_object[[]]
)

cluster_column_candidates <- grep(
  
  pattern =
    "res\\.0\\.6$",
  
  x =
    metadata_column_names,
  
  value =
    TRUE
)

if (
  length(
    cluster_column_candidates
  ) == 0
) {
  
  cluster_column_candidates <- grep(
    
    pattern =
      "res.0.6",
    
    x =
      metadata_column_names,
    
    value =
      TRUE,
    
    fixed =
      TRUE
  )
}

if (
  length(
    cluster_column_candidates
  ) == 0
) {
  
  stop(
    "Cluster column untuk resolution 0.6 tidak ditemukan."
  )
}

working_cluster_column <-
  cluster_column_candidates[1]

integrated_object$analysis_cluster <-
  factor(
    integrated_object[[]][
      ,
      working_cluster_column
    ]
  )

Seurat::Idents(
  integrated_object
) <- "analysis_cluster"

cat("\nWorking cluster column:\n")

print(
  working_cluster_column
)

cat("\nCluster distribution:\n")

print(
  table(
    integrated_object$
      analysis_cluster
  )
)

# ------------------------------------------------------------
# 27. Join RNA layers for later expression analysis
# ------------------------------------------------------------

rna_assay_class <- class(
  integrated_object[["RNA"]]
)

cat("\nRNA assay class before joining layers:\n")

print(
  rna_assay_class
)

if (
  inherits(
    integrated_object[["RNA"]],
    "Assay5"
  )
) {
  
  message(
    "Joining RNA layers..."
  )
  
  integrated_object <- SeuratObject::JoinLayers(
    
    object =
      integrated_object,
    
    assay =
      "RNA"
  )
}

cat("\nRNA layers after joining:\n")

print(
  SeuratObject::Layers(
    
    object =
      integrated_object,
    
    assay =
      "RNA"
  )
)

# ------------------------------------------------------------
# 28. Normalize RNA assay for biological visualization
# ------------------------------------------------------------

Seurat::DefaultAssay(
  integrated_object
) <- "RNA"

integrated_object <- Seurat::NormalizeData(
  
  object =
    integrated_object,
  
  assay =
    "RNA",
  
  normalization.method =
    "LogNormalize",
  
  scale.factor =
    10000,
  
  verbose =
    FALSE
)

Seurat::Idents(
  integrated_object
) <- "analysis_cluster"

# ------------------------------------------------------------
# 29. Integrated cell metadata
# ------------------------------------------------------------

integrated_cell_metadata <- integrated_object[[]] %>%
  
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  
  dplyr::mutate(
    
    sample_label = factor(
      sample_label,
      levels =
        expected_sample_order
    ),
    
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    )
  )

cat("\nIntegrated object dimensions:\n")

print(
  dim(
    integrated_object
  )
)

cat("\nFinal cells per donor:\n")

print(
  table(
    integrated_cell_metadata$
      sample_label
  )
)

cat("\nFinal cells per disease group:\n")

print(
  table(
    integrated_cell_metadata$
      group
  )
)

# ------------------------------------------------------------
# 30. Cluster composition tables
# ------------------------------------------------------------

cluster_sample_composition <-
  integrated_cell_metadata %>%
  
  dplyr::count(
    analysis_cluster,
    sample_label,
    group,
    name =
      "cell_number"
  ) %>%
  
  dplyr::group_by(
    sample_label
  ) %>%
  
  dplyr::mutate(
    
    proportion_within_sample =
      cell_number /
      sum(
        cell_number
      )
  ) %>%
  
  dplyr::ungroup()

cluster_group_composition <-
  integrated_cell_metadata %>%
  
  dplyr::count(
    analysis_cluster,
    group,
    name =
      "cell_number"
  ) %>%
  
  dplyr::group_by(
    group
  ) %>%
  
  dplyr::mutate(
    
    proportion_within_group =
      cell_number /
      sum(
        cell_number
      )
  ) %>%
  
  dplyr::ungroup()

cat("\nCluster by donor composition:\n")

print(
  cluster_sample_composition,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 31. Retention plot
# ------------------------------------------------------------

retention_plot_data <- retention_summary %>%
  
  dplyr::select(
    sample_label,
    group,
    raw_cells,
    cells_after_basic_QC,
    cells_after_doublet_removal
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      raw_cells,
      cells_after_basic_QC,
      cells_after_doublet_removal
    ),
    
    names_to =
      "processing_stage",
    
    values_to =
      "cell_number"
  ) %>%
  
  dplyr::mutate(
    
    processing_stage = factor(
      
      processing_stage,
      
      levels = c(
        "raw_cells",
        "cells_after_basic_QC",
        "cells_after_doublet_removal"
      ),
      
      labels = c(
        "Initial object",
        "After adaptive QC",
        "Final singlets"
      )
    ),
    
    sample_label = factor(
      sample_label,
      levels =
        expected_sample_order
    )
  )

retention_plot <- ggplot(
  
  retention_plot_data,
  
  aes(
    x = sample_label,
    y = cell_number,
    fill = processing_stage
  )
) +
  
  geom_col(
    
    position =
      position_dodge(
        width = 0.8
      ),
    
    width =
      0.72
  ) +
  
  labs(
    
    title =
      "GSE220243 cell retention during quality control",
    
    subtitle =
      paste(
        "Adaptive donor-specific filtering followed by",
        "per-donor scDblFinder classification"
      ),
    
    x =
      "Cartilage donor",
    
    y =
      "Number of cells",
    
    fill =
      "Processing stage"
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    
    legend.position =
      "top"
  )

print(
  retention_plot
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13A_GSE220243_QC_cell_retention.pdf"
  ),
  
  plot =
    retention_plot,
  
  width =
    9.5,
  
  height =
    5.8
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13A_GSE220243_QC_cell_retention.tiff"
  ),
  
  plot =
    retention_plot,
  
  width =
    9.5,
  
  height =
    5.8,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 32. PCA variance plot
# ------------------------------------------------------------

pca_variance_plot <- ggplot(
  
  pca_variance_table,
  
  aes(
    x = PC,
    y = variance_percent
  )
) +
  
  geom_point(
    size = 1.5
  ) +
  
  geom_line(
    linewidth = 0.45
  ) +
  
  geom_vline(
    
    xintercept =
      pcs_for_clustering,
    
    linetype =
      "dashed",
    
    linewidth =
      0.6
  ) +
  
  annotate(
    
    geom =
      "text",
    
    x =
      pcs_for_clustering,
    
    y =
      max(
        pca_variance_table$
          variance_percent
      ) *
      0.85,
    
    label =
      paste0(
        "Selected PCs = ",
        pcs_for_clustering
      ),
    
    hjust =
      -0.1,
    
    size =
      3.4
  ) +
  
  labs(
    
    title =
      "Integrated PCA variance profile",
    
    subtitle =
      "Dashed line indicates dimensions used for UMAP and clustering",
    
    x =
      "Principal component",
    
    y =
      "Variance explained (%)"
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    panel.grid.minor =
      element_blank()
  )

print(
  pca_variance_plot
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13B_GSE220243_integrated_PCA_variance.pdf"
  ),
  
  plot =
    pca_variance_plot,
  
  width =
    7,
  
  height =
    5.2
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13B_GSE220243_integrated_PCA_variance.tiff"
  ),
  
  plot =
    pca_variance_plot,
  
  width =
    7,
  
  height =
    5.2,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 33. Integrated UMAP plots
# ------------------------------------------------------------

Seurat::DefaultAssay(
  integrated_object
) <- "integrated"

umap_by_cluster <- Seurat::DimPlot(
  
  object =
    integrated_object,
  
  reduction =
    "umap",
  
  group.by =
    "analysis_cluster",
  
  label =
    TRUE,
  
  repel =
    TRUE,
  
  raster =
    TRUE,
  
  pt.size =
    0.15
) +
  
  labs(
    title =
      "Integrated cartilage cell clusters"
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    legend.position =
      "none"
  )

umap_by_group <- Seurat::DimPlot(
  
  object =
    integrated_object,
  
  reduction =
    "umap",
  
  group.by =
    "group",
  
  raster =
    TRUE,
  
  pt.size =
    0.15
) +
  
  labs(
    title =
      "Cells by cartilage group"
  ) +
  
  theme_bw(
    base_size = 10
  )

umap_by_sample <- Seurat::DimPlot(
  
  object =
    integrated_object,
  
  reduction =
    "umap",
  
  group.by =
    "sample_label",
  
  raster =
    TRUE,
  
  pt.size =
    0.15
) +
  
  labs(
    title =
      "Cells by donor"
  ) +
  
  theme_bw(
    base_size = 10
  )

integrated_umap_panel <-
  umap_by_cluster +
  umap_by_group +
  umap_by_sample +
  patchwork::plot_layout(
    ncol = 3,
    guides = "collect"
  ) +
  
  patchwork::plot_annotation(
    
    title =
      "GSE220243 cartilage SCT–RPCA integration",
    
    subtitle =
      paste(
        "Integrated embedding is used for visualization",
        "and clustering, not differential-expression testing"
      )
  )

print(
  integrated_umap_panel
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13C_GSE220243_integrated_UMAPs.pdf"
  ),
  
  plot =
    integrated_umap_panel,
  
  width =
    15,
  
  height =
    5.5
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13C_GSE220243_integrated_UMAPs.tiff"
  ),
  
  plot =
    integrated_umap_panel,
  
  width =
    15,
  
  height =
    5.5,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 34. Cluster composition plot
# ------------------------------------------------------------

cluster_composition_plot <- ggplot(
  
  cluster_sample_composition,
  
  aes(
    x = sample_label,
    y = proportion_within_sample,
    fill = analysis_cluster
  )
) +
  
  geom_col(
    width = 0.78
  ) +
  
  scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 1
      )
  ) +
  
  labs(
    
    title =
      "Cluster composition across cartilage donors",
    
    subtitle =
      "Each column represents the proportion of cells within one donor",
    
    x =
      "Cartilage donor",
    
    y =
      "Proportion of retained cells",
    
    fill =
      "Cluster"
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    
    legend.position =
      "right"
  )

print(
  cluster_composition_plot
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13D_GSE220243_cluster_composition.pdf"
  ),
  
  plot =
    cluster_composition_plot,
  
  width =
    10,
  
  height =
    6
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13D_GSE220243_cluster_composition.tiff"
  ),
  
  plot =
    cluster_composition_plot,
  
  width =
    10,
  
  height =
    6,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 35. Exploratory MMP13 feature plot
#
# This is not yet final Figure 1D because cell types have
# not been annotated.
# ------------------------------------------------------------

Seurat::DefaultAssay(
  integrated_object
) <- "RNA"

if (
  "MMP13" %in%
  rownames(
    integrated_object[["RNA"]]
  )
) {
  
  mmp13_feature_plot <- Seurat::FeaturePlot(
    
    object =
      integrated_object,
    
    features =
      "MMP13",
    
    reduction =
      "umap",
    
    order =
      TRUE,
    
    raster =
      TRUE,
    
    pt.size =
      0.18
  ) +
    
    labs(
      
      title =
        "Exploratory MMP13 localization",
      
      subtitle =
        paste(
          "Log-normalized RNA expression before",
          "cell-type annotation"
        )
    ) +
    
    theme_bw(
      base_size = 10
    )
  
  print(
    mmp13_feature_plot
  )
  
  ggsave(
    
    filename = file.path(
      supplementary_figure_folder,
      "FigureS13E_GSE220243_exploratory_MMP13_UMAP.pdf"
    ),
    
    plot =
      mmp13_feature_plot,
    
    width =
      7,
    
    height =
      5.8
  )
  
  ggsave(
    
    filename = file.path(
      supplementary_figure_folder,
      "FigureS13E_GSE220243_exploratory_MMP13_UMAP.tiff"
    ),
    
    plot =
      mmp13_feature_plot,
    
    width =
      7,
    
    height =
      5.8,
    
    dpi =
      600,
    
    compression =
      "lzw"
  )
}

# ------------------------------------------------------------
# 36. Post-QC MMP13 donor detection plot
# ------------------------------------------------------------

mmp13_post_QC_summary <- mmp13_post_QC_summary %>%
  
  dplyr::mutate(
    
    sample_label = factor(
      sample_label,
      levels =
        expected_sample_order
    ),
    
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    )
  )

mmp13_detection_plot <- ggplot(
  
  mmp13_post_QC_summary,
  
  aes(
    x = sample_label,
    y = MMP13_detection_percent,
    fill = group
  )
) +
  
  geom_col(
    width = 0.75
  ) +
  
  geom_text(
    
    aes(
      label =
        sprintf(
          "%.2f%%",
          MMP13_detection_percent
        )
    ),
    
    vjust =
      -0.3,
    
    size =
      3
  ) +
  
  labs(
    
    title =
      "MMP13 detection after QC and doublet removal",
    
    subtitle =
      paste(
        "Percentage of retained singlet cells",
        "with at least one MMP13 count"
      ),
    
    x =
      "Cartilage donor",
    
    y =
      "MMP13-positive singlets (%)",
    
    fill =
      "Group"
  ) +
  
  theme_bw(
    base_size = 10
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1
      ),
    
    legend.position =
      "top"
  ) +
  
  expand_limits(
    
    y =
      max(
        mmp13_post_QC_summary$
          MMP13_detection_percent,
        na.rm = TRUE
      ) *
      1.15
  )

print(
  mmp13_detection_plot
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13F_GSE220243_postQC_MMP13_detection.pdf"
  ),
  
  plot =
    mmp13_detection_plot,
  
  width =
    9,
  
  height =
    5.5
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS13F_GSE220243_postQC_MMP13_detection.tiff"
  ),
  
  plot =
    mmp13_detection_plot,
  
  width =
    9,
  
  height =
    5.5,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 37. Analysis notes
# ------------------------------------------------------------

analysis_notes <- tibble::tibble(
  
  topic = c(
    "Basic cell QC",
    "Doublet detection",
    "Normalization",
    "Integration",
    "Integration references",
    "Clustering representation",
    "Expression representation",
    "Differential expression",
    "Figure 1D status"
  ),
  
  decision = c(
    
    paste(
      "Donor-specific median/MAD thresholds were applied",
      "to detected genes, UMI counts, and mitochondrial",
      "transcript percentage."
    ),
    
    paste(
      "scDblFinder was run independently for every",
      "cartilage donor after basic QC."
    ),
    
    paste(
      "SCTransform v2 was performed separately for",
      "each donor with percent.mt regressed."
    ),
    
    paste(
      "Donors were integrated using SCT-RPCA anchors."
    ),
    
    paste(
      "One Normal and one OA donor with the largest",
      "retained cell numbers were used as balanced",
      "integration references."
    ),
    
    paste(
      "The integrated assay and PCA/UMAP are used",
      "for visualization and clustering."
    ),
    
    paste(
      "Joined and log-normalized RNA counts are retained",
      "for expression visualization."
    ),
    
    paste(
      "The integrated assay must not be used for",
      "differential-expression testing. Donor-aware",
      "pseudobulk analysis will be used later."
    ),
    
    paste(
      "Cell-type annotation has not yet been performed;",
      "the MMP13 UMAP is exploratory and is not the",
      "final Figure 1D."
    )
  )
)

# ------------------------------------------------------------
# 38. Save integration features and supporting objects
# ------------------------------------------------------------

saveRDS(
  
  integration_features,
  
  file = file.path(
    processed_folder,
    "GSE220243_SCT_integration_features.rds"
  )
)

saveRDS(
  
  integration_anchors,
  
  file = file.path(
    processed_folder,
    "GSE220243_SCT_RPCA_integration_anchors.rds"
  ),
  
  compress =
    "gzip"
)

# ------------------------------------------------------------
# 39. Save final integrated object
# ------------------------------------------------------------

Seurat::DefaultAssay(
  integrated_object
) <- "RNA"

Seurat::Idents(
  integrated_object
) <- "analysis_cluster"

final_integrated_file <- file.path(
  
  processed_folder,
  
  "GSE220243_cartilage_integrated_seurat.rds"
)

message(
  "Menyimpan integrated Seurat object..."
)

saveRDS(
  
  integrated_object,
  
  file =
    final_integrated_file,
  
  compress =
    "gzip"
)

message(
  "Integrated Seurat object berhasil disimpan."
)

# ------------------------------------------------------------
# 40. Save CSV tables
# ------------------------------------------------------------

utils::write.csv(
  
  qc_threshold_table,
  
  file = file.path(
    table_folder,
    "GSE220243_donor_specific_QC_thresholds.csv"
  ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  qc_failure_summary,
  
  file = file.path(
    table_folder,
    "GSE220243_basic_QC_failure_summary.csv"
  ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  doublet_summary,
  
  file = file.path(
    table_folder,
    "GSE220243_scDblFinder_summary.csv"
  ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  retention_summary,
  
  file = file.path(
    table_folder,
    "GSE220243_final_cell_retention_summary.csv"
  ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  mmp13_post_QC_summary,
  
  file = file.path(
    table_folder,
    "GSE220243_postQC_MMP13_detection.csv"
  ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  pca_variance_table,
  
  file = file.path(
    table_folder,
    "GSE220243_integrated_PCA_variance.csv"
  ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  cluster_sample_composition,
  
  file = file.path(
    table_folder,
    "GSE220243_cluster_sample_composition.csv"
  ),
  
  row.names =
    FALSE
)

# ------------------------------------------------------------
# 41. Save large cell-level tables as compressed CSV
# ------------------------------------------------------------

qc_flag_connection <- gzfile(
  
  file.path(
    table_folder,
    "GSE220243_cell_level_QC_flags.csv.gz"
  ),
  
  open =
    "wt"
)

utils::write.csv(
  
  all_cell_qc_flags,
  
  file =
    qc_flag_connection,
  
  row.names =
    FALSE
)

close(
  qc_flag_connection
)

doublet_connection <- gzfile(
  
  file.path(
    table_folder,
    "GSE220243_cell_doublet_scores.csv.gz"
  ),
  
  open =
    "wt"
)

utils::write.csv(
  
  all_doublet_scores,
  
  file =
    doublet_connection,
  
  row.names =
    FALSE
)

close(
  doublet_connection
)

metadata_connection <- gzfile(
  
  file.path(
    table_folder,
    "GSE220243_integrated_cell_metadata.csv.gz"
  ),
  
  open =
    "wt"
)

utils::write.csv(
  
  integrated_cell_metadata,
  
  file =
    metadata_connection,
  
  row.names =
    FALSE
)

close(
  metadata_connection
)

# ------------------------------------------------------------
# 42. Save Excel workbook
# ------------------------------------------------------------

openxlsx::write.xlsx(
  
  list(
    
    Sample_Design =
      as.data.frame(
        sample_design
      ),
    
    QC_Thresholds =
      as.data.frame(
        qc_threshold_table
      ),
    
    QC_Failure_Summary =
      as.data.frame(
        qc_failure_summary
      ),
    
    Doublet_Summary =
      as.data.frame(
        doublet_summary
      ),
    
    Cell_Retention =
      as.data.frame(
        retention_summary
      ),
    
    MMP13_PostQC =
      as.data.frame(
        mmp13_post_QC_summary
      ),
    
    Integration_Parameters =
      as.data.frame(
        integration_parameter_summary
      ),
    
    Reference_Donors =
      as.data.frame(
        reference_donor_table
      ),
    
    PCA_Variance =
      as.data.frame(
        pca_variance_table
      ),
    
    Cluster_Group_Composition =
      as.data.frame(
        cluster_group_composition
      ),
    
    Analysis_Notes =
      as.data.frame(
        analysis_notes
      ),
    
    Package_Versions =
      as.data.frame(
        package_versions
      )
  ),
  
  file = file.path(
    table_folder,
    paste0(
      "GSE220243_QC_doublet_",
      "SCT_RPCA_integration_results.xlsx"
    )
  ),
  
  overwrite =
    TRUE
)

# ------------------------------------------------------------
# 43. Save session information
# ------------------------------------------------------------

sink(
  
  file.path(
    table_folder,
    paste0(
      "GSE220243_QC_integration_",
      "sessionInfo.txt"
    )
  )
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 44. Final summary
# ------------------------------------------------------------

total_initial_cells <- sum(
  retention_summary$
    raw_cells
)

total_after_basic_QC <- sum(
  retention_summary$
    cells_after_basic_QC
)

total_predicted_doublets <- sum(
  retention_summary$
    predicted_doublets
)

total_final_singlets <- sum(
  retention_summary$
    cells_after_doublet_removal
)

normal_final_cells <- sum(
  integrated_object$group ==
    "Normal"
)

oa_final_cells <- sum(
  integrated_object$group ==
    "OA"
)

total_clusters <- dplyr::n_distinct(
  integrated_object$
    analysis_cluster
)

message("")
message("================================================")
message("GSE220243 QC, DOUBLET, AND INTEGRATION SELESAI")
message("================================================")

message(
  "Initial cells                  : ",
  format(
    total_initial_cells,
    big.mark = ","
  )
)

message(
  "Cells after adaptive QC        : ",
  format(
    total_after_basic_QC,
    big.mark = ","
  )
)

message(
  "Predicted doublets removed     : ",
  format(
    total_predicted_doublets,
    big.mark = ","
  )
)

message(
  "Final retained singlets        : ",
  format(
    total_final_singlets,
    big.mark = ","
  )
)

message(
  "Final overall retention        : ",
  round(
    100 *
      total_final_singlets /
      total_initial_cells,
    2
  ),
  "%"
)

message(
  "Final Normal cells             : ",
  format(
    normal_final_cells,
    big.mark = ","
  )
)

message(
  "Final OA cells                 : ",
  format(
    oa_final_cells,
    big.mark = ","
  )
)

message(
  "Integration features           : ",
  format(
    length(
      integration_features
    ),
    big.mark = ","
  )
)

message(
  "PCs used                       : ",
  pcs_for_clustering
)

message(
  "Working cluster resolution     : 0.6"
)

message(
  "Number of working clusters     : ",
  total_clusters
)

message(
  "Final integrated object        : ",
  final_integrated_file
)

message(
  "Supplementary figures          : ",
  supplementary_figure_folder
)

message(
  "Cell-type annotation performed : FALSE"
)

message(
  "Final Figure 1D completed      : FALSE"
)

message("================================================")