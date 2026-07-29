# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 16_download_prepare_GSE220243_cartilage_scRNA.R
# PURPOSE : Download and prepare raw cartilage scRNA-seq data
# DATASET : GSE220243
#
# INCLUDED:
#   6 normal cartilage donors
#   6 OA cartilage donors
#
# EXCLUDED AT THIS STAGE:
#   Meniscus samples
#
# IMPORTANT:
#   This script performs preparation and exploratory QC only.
#   It does NOT yet filter cells, integrate donors, or assign
#   cell types.
# ============================================================

rm(list = ls())
gc()
# ------------------------------------------------------------
# PROJECT ROOT
# ------------------------------------------------------------

project_root <- paste0(
  "D:/PhD FILE Ymelda/",
  "Publikasi Protein buatan/",
  "OA_MMP13_Target_Discovery"
)

if (!dir.exists(project_root)) {
  stop(
    "Project folder tidak ditemukan:\n",
    project_root
  )
}

setwd(project_root)

cat("\nWorking directory aktif:\n")
print(getwd())

stopifnot(
  normalizePath(
    getwd(),
    winslash = "/",
    mustWork = TRUE
  ) ==
    normalizePath(
      project_root,
      winslash = "/",
      mustWork = TRUE
    )
)

message(
  "Working directory sudah diarahkan ke project root."
)
# ------------------------------------------------------------
# 1. Package installation and loading
# ------------------------------------------------------------

cran_packages <- c(
  "Seurat",
  "Matrix",
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "patchwork",
  "openxlsx"
)

bioconductor_packages <- c(
  "GEOquery"
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
library(Matrix)
library(GEOquery)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(patchwork)
library(openxlsx)

message(
  "Semua paket single-cell berhasil dimuat."
)

# ------------------------------------------------------------
# 2. Package versions
# ------------------------------------------------------------

package_versions <- tibble::tibble(
  
  package = c(
    "Seurat",
    "SeuratObject",
    "Matrix",
    "GEOquery"
  ),
  
  version = c(
    as.character(
      packageVersion("Seurat")
    ),
    
    as.character(
      packageVersion("SeuratObject")
    ),
    
    as.character(
      packageVersion("Matrix")
    ),
    
    as.character(
      packageVersion("GEOquery")
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
# 3. Download configuration
# ------------------------------------------------------------

options(
  timeout = 7200
)

options(
  download.file.method.GEOquery =
    "libcurl"
)

raw_data_folder <-
  "data_raw/GSE220243_cartilage"

processed_data_folder <-
  "data_processed"

table_folder <-
  "results/tables"

supplementary_figure_folder <-
  paste0(
    "results/figures/",
    "Figure 1 Suplementary/",
    "S12"
  )

required_folders <- c(
  raw_data_folder,
  processed_data_folder,
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

# ------------------------------------------------------------
# 4. Cartilage sample metadata
# ------------------------------------------------------------

sample_table <- tibble::tribble(
  
  ~geo_accession,
  ~sample_label,
  ~group,
  ~donor,
  ~tissue,
  
  "GSM6797148",
  "Normal1",
  "Normal",
  "Normal_Donor1",
  "Cartilage",
  
  "GSM6797149",
  "Normal2",
  "Normal",
  "Normal_Donor2",
  "Cartilage",
  
  "GSM6797150",
  "Normal3",
  "Normal",
  "Normal_Donor3",
  "Cartilage",
  
  "GSM6797151",
  "Normal4",
  "Normal",
  "Normal_Donor4",
  "Cartilage",
  
  "GSM6797152",
  "Normal5",
  "Normal",
  "Normal_Donor5",
  "Cartilage",
  
  "GSM6797153",
  "Normal6",
  "Normal",
  "Normal_Donor6",
  "Cartilage",
  
  "GSM6797154",
  "OA1",
  "OA",
  "OA_Donor1",
  "Cartilage",
  
  "GSM6797155",
  "OA2",
  "OA",
  "OA_Donor2",
  "Cartilage",
  
  "GSM6797156",
  "OA3",
  "OA",
  "OA_Donor3",
  "Cartilage",
  
  "GSM6797157",
  "OA4",
  "OA",
  "OA_Donor4",
  "Cartilage",
  
  "GSM6797158",
  "OA5",
  "OA",
  "OA_Donor5",
  "Cartilage",
  
  "GSM6797159",
  "OA6",
  "OA",
  "OA_Donor6",
  "Cartilage"
)

sample_table <- sample_table %>%
  dplyr::mutate(
    
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    ),
    
    sample_order =
      dplyr::row_number()
  )

cat("\nSampel kartilago yang akan dianalisis:\n")

print(
  sample_table,
  n = Inf,
  width = Inf
)

cat("\nJumlah sampel per kelompok:\n")

print(
  table(
    sample_table$group
  )
)

stopifnot(
  nrow(sample_table) == 12
)

stopifnot(
  sum(sample_table$group == "Normal") == 6
)

stopifnot(
  sum(sample_table$group == "OA") == 6
)

# ------------------------------------------------------------
# 5. Function to identify downloaded 10X files
# ------------------------------------------------------------

find_single_file <- function(
    search_directory,
    filename_pattern,
    file_description
) {
  
  matching_files <- list.files(
    search_directory,
    pattern = filename_pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  matching_files <- matching_files[
    file.info(matching_files)$isdir ==
      FALSE
  ]
  
  if (length(matching_files) == 0) {
    
    stop(
      "Tidak ditemukan ",
      file_description,
      " pada folder:\n",
      search_directory
    )
  }
  
  if (length(matching_files) > 1) {
    
    cat(
      "\nBeberapa kandidat ditemukan untuk ",
      file_description,
      ":\n",
      sep = ""
    )
    
    print(
      matching_files
    )
    
    stop(
      "Ditemukan lebih dari satu kandidat ",
      file_description,
      "."
    )
  }
  
  matching_files[1]
}

# ------------------------------------------------------------
# 6. Function to standardize 10X filenames
# ------------------------------------------------------------

prepare_standard_10x_directory <- function(
    gsm_id,
    gsm_directory
) {
  
  standard_directory <- file.path(
    gsm_directory,
    "filtered_feature_bc_matrix"
  )
  
  dir.create(
    standard_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  standard_matrix_file <- file.path(
    standard_directory,
    "matrix.mtx.gz"
  )
  
  standard_feature_file <- file.path(
    standard_directory,
    "features.tsv.gz"
  )
  
  standard_barcode_file <- file.path(
    standard_directory,
    "barcodes.tsv.gz"
  )
  
  standard_files_exist <- all(
    file.exists(
      c(
        standard_matrix_file,
        standard_feature_file,
        standard_barcode_file
      )
    )
  )
  
  if (!standard_files_exist) {
    
    original_matrix_file <- find_single_file(
      search_directory =
        gsm_directory,
      
      filename_pattern =
        "_matrix\\.mtx\\.gz$",
      
      file_description =
        paste0(
          gsm_id,
          " matrix.mtx.gz"
        )
    )
    
    original_feature_file <- find_single_file(
      search_directory =
        gsm_directory,
      
      filename_pattern =
        "_features\\.tsv\\.gz$",
      
      file_description =
        paste0(
          gsm_id,
          " features.tsv.gz"
        )
    )
    
    original_barcode_file <- find_single_file(
      search_directory =
        gsm_directory,
      
      filename_pattern =
        "_barcodes\\.tsv\\.gz$",
      
      file_description =
        paste0(
          gsm_id,
          " barcodes.tsv.gz"
        )
    )
    
    copy_matrix_success <- file.copy(
      from =
        original_matrix_file,
      
      to =
        standard_matrix_file,
      
      overwrite =
        TRUE
    )
    
    copy_feature_success <- file.copy(
      from =
        original_feature_file,
      
      to =
        standard_feature_file,
      
      overwrite =
        TRUE
    )
    
    copy_barcode_success <- file.copy(
      from =
        original_barcode_file,
      
      to =
        standard_barcode_file,
      
      overwrite =
        TRUE
    )
    
    if (
      !all(
        c(
          copy_matrix_success,
          copy_feature_success,
          copy_barcode_success
        )
      )
    ) {
      
      stop(
        "Gagal menyalin file 10X standar untuk ",
        gsm_id
      )
    }
  }
  
  final_files <- c(
    matrix =
      standard_matrix_file,
    
    features =
      standard_feature_file,
    
    barcodes =
      standard_barcode_file
  )
  
  if (
    !all(
      file.exists(
        final_files
      )
    )
  ) {
    
    stop(
      "File 10X standar belum lengkap untuk ",
      gsm_id
    )
  }
  
  list(
    
    directory =
      standard_directory,
    
    files =
      final_files
  )
}

# ------------------------------------------------------------
# 7. Download supplementary files per cartilage GSM
# ------------------------------------------------------------

download_records <- list()
standard_10x_directories <- list()

for (
  current_index in
  seq_len(
    nrow(sample_table)
  )
) {
  
  current_gsm <- as.character(
    sample_table$geo_accession[
      current_index
    ]
  )
  
  current_label <- as.character(
    sample_table$sample_label[
      current_index
    ]
  )
  
  message("")
  message(
    "============================================"
  )
  
  message(
    "Memproses ",
    current_gsm,
    " — ",
    current_label
  )
  
  message(
    "============================================"
  )
  
  current_gsm_directory <- file.path(
    raw_data_folder,
    current_gsm
  )
  
  expected_standard_files <- c(
    
    file.path(
      current_gsm_directory,
      "filtered_feature_bc_matrix",
      "matrix.mtx.gz"
    ),
    
    file.path(
      current_gsm_directory,
      "filtered_feature_bc_matrix",
      "features.tsv.gz"
    ),
    
    file.path(
      current_gsm_directory,
      "filtered_feature_bc_matrix",
      "barcodes.tsv.gz"
    )
  )
  
  if (
    all(
      file.exists(
        expected_standard_files
      )
    )
  ) {
    
    message(
      "File standar sudah tersedia. Download dilewati."
    )
    
  } else {
    
    message(
      "Memeriksa supplementary file GEO..."
    )
    
    supplementary_listing <- GEOquery::getGEOSuppFiles(
      
      GEO =
        current_gsm,
      
      makeDirectory =
        TRUE,
      
      baseDir =
        raw_data_folder,
      
      fetch_files =
        FALSE,
      
      filter_regex =
        "(barcodes|features|matrix).*(tsv|mtx)\\.gz$"
    )
    
    cat(
      "\nDaftar supplementary files ",
      current_gsm,
      ":\n",
      sep = ""
    )
    
    print(
      supplementary_listing
    )
    
    if (
      nrow(
        supplementary_listing
      ) < 3
    ) {
      
      stop(
        "Supplementary file 10X tidak lengkap untuk ",
        current_gsm
      )
    }
    
    message(
      "Mengunduh tiga file 10X..."
    )
    
    download_result <- tryCatch(
      
      GEOquery::getGEOSuppFiles(
        
        GEO =
          current_gsm,
        
        makeDirectory =
          TRUE,
        
        baseDir =
          raw_data_folder,
        
        fetch_files =
          TRUE,
        
        filter_regex =
          "(barcodes|features|matrix).*(tsv|mtx)\\.gz$"
      ),
      
      error = function(error_object) {
        
        stop(
          "Download gagal untuk ",
          current_gsm,
          ":\n",
          conditionMessage(
            error_object
          )
        )
      }
    )
  }
  
  standardized_result <-
    prepare_standard_10x_directory(
      
      gsm_id =
        current_gsm,
      
      gsm_directory =
        current_gsm_directory
    )
  
  standard_10x_directories[[current_gsm]] <-
    standardized_result$directory
  
  current_file_info <- file.info(
    standardized_result$files
  )
  
  download_records[[current_gsm]] <-
    tibble::tibble(
      
      geo_accession =
        current_gsm,
      
      sample_label =
        current_label,
      
      file_type =
        names(
          standardized_result$files
        ),
      
      file_path =
        unname(
          standardized_result$files
        ),
      
      file_size_bytes =
        current_file_info$size,
      
      file_size_MB =
        round(
          current_file_info$size /
            1024^2,
          3
        ),
      
      file_exists =
        file.exists(
          standardized_result$files
        )
    )
  
  message(
    "File 10X siap untuk ",
    current_gsm
  )
}

download_manifest <- dplyr::bind_rows(
  download_records
)

cat("\nDownload manifest:\n")

print(
  download_manifest,
  n = Inf,
  width = Inf
)

stopifnot(
  nrow(download_manifest) == 36
)

stopifnot(
  all(
    download_manifest$file_exists
  )
)

# ------------------------------------------------------------
# 8. Read each 10X matrix and create raw Seurat objects
# ------------------------------------------------------------

seurat_object_list <- list()
sample_qc_summaries <- list()
mmp13_detection_summaries <- list()

for (
  current_index in
  seq_len(
    nrow(sample_table)
  )
) {
  
  current_gsm <- as.character(
    sample_table$geo_accession[
      current_index
    ]
  )
  
  current_label <- as.character(
    sample_table$sample_label[
      current_index
    ]
  )
  
  current_group <- as.character(
    sample_table$group[
      current_index
    ]
  )
  
  current_donor <- as.character(
    sample_table$donor[
      current_index
    ]
  )
  
  current_10x_directory <-
    standard_10x_directories[[current_gsm]]
  
  message("")
  message(
    "Membaca matriks ",
    current_label,
    "..."
  )
  
  current_counts <- Seurat::Read10X(
    
    data.dir =
      current_10x_directory,
    
    gene.column =
      2,
    
    unique.features =
      TRUE,
    
    strip.suffix =
      FALSE
  )
  
  if (is.list(current_counts)) {
    
    if (
      "Gene Expression" %in%
      names(current_counts)
    ) {
      
      current_counts <-
        current_counts[["Gene Expression"]]
      
    } else {
      
      warning(
        current_label,
        ": Read10X menghasilkan beberapa assay. ",
        "Assay pertama digunakan."
      )
      
      current_counts <-
        current_counts[[1]]
    }
  }
  
  if (
    !inherits(
      current_counts,
      "sparseMatrix"
    )
  ) {
    
    current_counts <- methods::as(
      current_counts,
      "dgCMatrix"
    )
  }
  
  raw_gene_number <- nrow(
    current_counts
  )
  
  raw_barcode_number <- ncol(
    current_counts
  )
  
  message(
    current_label,
    ": ",
    raw_gene_number,
    " genes × ",
    raw_barcode_number,
    " barcodes"
  )
  
  current_object <- Seurat::CreateSeuratObject(
    
    counts =
      current_counts,
    
    project =
      "GSE220243_cartilage",
    
    min.cells =
      3,
    
    min.features =
      200
  )
  
  current_object <- SeuratObject::RenameCells(
    
    object =
      current_object,
    
    add.cell.id =
      current_label
  )
  
  current_object$geo_accession <-
    current_gsm
  
  current_object$sample_label <-
    current_label
  
  current_object$group <-
    current_group
  
  current_object$donor <-
    current_donor
  
  current_object$tissue <-
    "Cartilage"
  
  current_object$dataset <-
    "GSE220243"
  
  current_object$percent.mt <-
    Seurat::PercentageFeatureSet(
      
      object =
        current_object,
      
      pattern =
        "^MT-"
    )
  
  current_object$percent.ribo <-
    Seurat::PercentageFeatureSet(
      
      object =
        current_object,
      
      pattern =
        "^RP[SL]"
    )
  
  current_metadata <- current_object[[]]
  
  sample_qc_summaries[[current_gsm]] <-
    tibble::tibble(
      
      geo_accession =
        current_gsm,
      
      sample_label =
        current_label,
      
      group =
        current_group,
      
      donor =
        current_donor,
      
      raw_gene_number =
        raw_gene_number,
      
      raw_barcode_number =
        raw_barcode_number,
      
      cells_after_initial_creation =
        ncol(
          current_object
        ),
      
      genes_in_seurat_object =
        nrow(
          current_object
        ),
      
      median_nFeature_RNA =
        stats::median(
          current_metadata$nFeature_RNA
        ),
      
      mean_nFeature_RNA =
        mean(
          current_metadata$nFeature_RNA
        ),
      
      minimum_nFeature_RNA =
        min(
          current_metadata$nFeature_RNA
        ),
      
      maximum_nFeature_RNA =
        max(
          current_metadata$nFeature_RNA
        ),
      
      median_nCount_RNA =
        stats::median(
          current_metadata$nCount_RNA
        ),
      
      mean_nCount_RNA =
        mean(
          current_metadata$nCount_RNA
        ),
      
      minimum_nCount_RNA =
        min(
          current_metadata$nCount_RNA
        ),
      
      maximum_nCount_RNA =
        max(
          current_metadata$nCount_RNA
        ),
      
      median_percent_mt =
        stats::median(
          current_metadata$percent.mt
        ),
      
      mean_percent_mt =
        mean(
          current_metadata$percent.mt
        ),
      
      maximum_percent_mt =
        max(
          current_metadata$percent.mt
        ),
      
      median_percent_ribo =
        stats::median(
          current_metadata$percent.ribo
        )
    )
  
  if (
    "MMP13" %in%
    rownames(
      current_counts
    )
  ) {
    
    mmp13_counts <- as.numeric(
      current_counts[
        "MMP13",
        ,
        drop = TRUE
      ]
    )
    
    mmp13_positive_cells <- sum(
      mmp13_counts > 0
    )
    
    mmp13_detection_percentage <-
      100 *
      mmp13_positive_cells /
      length(
        mmp13_counts
      )
    
    mmp13_total_counts <- sum(
      mmp13_counts
    )
    
  } else {
    
    mmp13_positive_cells <-
      0
    
    mmp13_detection_percentage <-
      0
    
    mmp13_total_counts <-
      0
  }
  
  mmp13_detection_summaries[[current_gsm]] <-
    tibble::tibble(
      
      geo_accession =
        current_gsm,
      
      sample_label =
        current_label,
      
      group =
        current_group,
      
      donor =
        current_donor,
      
      MMP13_present_in_feature_matrix =
        "MMP13" %in%
        rownames(
          current_counts
        ),
      
      total_raw_barcodes =
        ncol(
          current_counts
        ),
      
      MMP13_positive_raw_barcodes =
        mmp13_positive_cells,
      
      MMP13_detection_percent =
        mmp13_detection_percentage,
      
      MMP13_total_raw_counts =
        mmp13_total_counts
    )
  
  seurat_object_list[[current_label]] <-
    current_object
  
  rm(
    current_counts,
    current_object
  )
  
  gc()
}

# ------------------------------------------------------------
# 9. Combine QC summaries
# ------------------------------------------------------------

sample_qc_summary <- dplyr::bind_rows(
  sample_qc_summaries
) %>%
  dplyr::left_join(
    
    sample_table %>%
      dplyr::select(
        geo_accession,
        sample_order
      ),
    
    by =
      "geo_accession"
  ) %>%
  dplyr::arrange(
    sample_order
  )

mmp13_detection_summary <- dplyr::bind_rows(
  mmp13_detection_summaries
) %>%
  dplyr::left_join(
    
    sample_table %>%
      dplyr::select(
        geo_accession,
        sample_order
      ),
    
    by =
      "geo_accession"
  ) %>%
  dplyr::arrange(
    sample_order
  )

cat("\nSample QC summary:\n")

print(
  sample_qc_summary,
  n = Inf,
  width = Inf
)

cat("\nMMP13 raw detection summary:\n")

print(
  mmp13_detection_summary,
  n = Inf,
  width = Inf
)

stopifnot(
  nrow(sample_qc_summary) == 12
)

stopifnot(
  nrow(mmp13_detection_summary) == 12
)

# ------------------------------------------------------------
# 10. Save the individual raw Seurat-object list
# ------------------------------------------------------------

saveRDS(
  
  seurat_object_list,
  
  file =
    file.path(
      processed_data_folder,
      "GSE220243_cartilage_raw_seurat_list.rds"
    ),
  
  compress =
    FALSE
)

message(
  "Raw Seurat-object list berhasil disimpan."
)

# ------------------------------------------------------------
# 11. Merge all donors for global QC visualization
# ------------------------------------------------------------

message(
  "Menggabungkan 12 raw Seurat objects..."
)

merged_raw_object <- merge(
  
  x =
    seurat_object_list[[1]],
  
  y =
    seurat_object_list[-1],
  
  project =
    "GSE220243_cartilage",
  
  merge.data =
    FALSE
)

merged_raw_object$group <- factor(
  
  merged_raw_object$group,
  
  levels = c(
    "Normal",
    "OA"
  )
)

merged_raw_object$sample_label <- factor(
  
  merged_raw_object$sample_label,
  
  levels =
    sample_table$sample_label
)

cat("\nMerged raw Seurat object:\n")

print(
  merged_raw_object
)

cat("\nJumlah sel per sampel:\n")

print(
  table(
    merged_raw_object$sample_label
  )
)

cat("\nJumlah sel per kelompok:\n")

print(
  table(
    merged_raw_object$group
  )
)

stopifnot(
  ncol(merged_raw_object) ==
    sum(
      sample_qc_summary$
        cells_after_initial_creation
    )
)

# ------------------------------------------------------------
# 12. Save merged raw object
# ------------------------------------------------------------

saveRDS(
  
  merged_raw_object,
  
  file =
    file.path(
      processed_data_folder,
      "GSE220243_cartilage_merged_raw_seurat.rds"
    ),
  
  compress =
    FALSE
)

message(
  "Merged raw Seurat object berhasil disimpan."
)

# ------------------------------------------------------------
# 13. Prepare cell-level QC metadata
# ------------------------------------------------------------

cell_qc_metadata <- merged_raw_object[[]] %>%
  
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  
  dplyr::mutate(
    
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    ),
    
    sample_label = factor(
      sample_label,
      levels =
        sample_table$sample_label
    )
  )

cat("\nDimensi cell-level QC metadata:\n")

print(
  dim(
    cell_qc_metadata
  )
)

# ------------------------------------------------------------
# 14. QC distribution figure
# ------------------------------------------------------------

qc_distribution_data <- cell_qc_metadata %>%
  
  dplyr::select(
    cell_id,
    sample_label,
    group,
    nFeature_RNA,
    nCount_RNA,
    percent.mt,
    percent.ribo
  ) %>%
  
  tidyr::pivot_longer(
    
    cols = c(
      nFeature_RNA,
      nCount_RNA,
      percent.mt,
      percent.ribo
    ),
    
    names_to =
      "QC_metric",
    
    values_to =
      "QC_value"
  ) %>%
  
  dplyr::mutate(
    
    QC_metric = factor(
      
      QC_metric,
      
      levels = c(
        "nFeature_RNA",
        "nCount_RNA",
        "percent.mt",
        "percent.ribo"
      ),
      
      labels = c(
        "Detected genes",
        "UMI counts",
        "Mitochondrial transcripts (%)",
        "Ribosomal transcripts (%)"
      )
    )
  )

qc_distribution_plot <- ggplot(
  
  qc_distribution_data,
  
  aes(
    x = sample_label,
    y = QC_value,
    fill = group
  )
) +
  
  geom_violin(
    scale = "width",
    trim = TRUE,
    linewidth = 0.25
  ) +
  
  facet_wrap(
    facets =
      vars(
        QC_metric
      ),
    ncol = 1,
    scales = "free_y"
  ) +
  
  labs(
    
    title =
      "GSE220243 cartilage single-cell QC distributions",
    
    subtitle =
      paste(
        "Raw objects after only minimal CreateSeuratObject",
        "requirements; no final cell filtering applied"
      ),
    
    x =
      "Cartilage donor",
    
    y =
      NULL,
    
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
        hjust = 1,
        vjust = 1
      ),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "top"
  )

print(
  qc_distribution_plot
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12A_GSE220243_raw_QC_distributions.pdf"
    ),
  
  plot =
    qc_distribution_plot,
  
  width =
    10,
  
  height =
    10
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12A_GSE220243_raw_QC_distributions.tiff"
    ),
  
  plot =
    qc_distribution_plot,
  
  width =
    10,
  
  height =
    10,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 15. nCount versus nFeature scatter plot
# ------------------------------------------------------------

qc_scatter_plot <- ggplot(
  
  cell_qc_metadata,
  
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    color = percent.mt
  )
) +
  
  geom_point(
    size = 0.35,
    alpha = 0.25
  ) +
  
  scale_x_log10() +
  
  scale_y_log10() +
  
  facet_wrap(
    facets =
      vars(
        sample_label
      ),
    ncol =
      4
  ) +
  
  labs(
    
    title =
      "Relationship between UMI depth and detected genes",
    
    subtitle =
      "Color represents mitochondrial transcript percentage",
    
    x =
      "UMI counts per cell, log10 scale",
    
    y =
      "Detected genes per cell, log10 scale",
    
    color =
      "Mitochondrial\ntranscripts (%)"
  ) +
  
  theme_bw(
    base_size = 9
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold"
      )
  )

print(
  qc_scatter_plot
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12B_GSE220243_count_feature_scatter.pdf"
    ),
  
  plot =
    qc_scatter_plot,
  
  width =
    10,
  
  height =
    8
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12B_GSE220243_count_feature_scatter.tiff"
    ),
  
  plot =
    qc_scatter_plot,
  
  width =
    10,
  
  height =
    8,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 16. Cell-number plot
# ------------------------------------------------------------

cell_number_plot_data <- sample_qc_summary %>%
  
  dplyr::mutate(
    
    sample_label = factor(
      sample_label,
      levels =
        sample_table$sample_label
    )
  )

cell_number_plot <- ggplot(
  
  cell_number_plot_data,
  
  aes(
    x = sample_label,
    y = cells_after_initial_creation,
    fill = group
  )
) +
  
  geom_col(
    width = 0.75
  ) +
  
  geom_text(
    
    aes(
      label =
        format(
          cells_after_initial_creation,
          big.mark = ","
        )
    ),
    
    vjust = -0.35,
    size = 3
  ) +
  
  labs(
    
    title =
      "Cell recovery across GSE220243 cartilage donors",
    
    subtitle =
      "Cell numbers before final QC filtering",
    
    x =
      "Cartilage donor",
    
    y =
      "Number of retained barcodes",
    
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
        cell_number_plot_data$
          cells_after_initial_creation
      ) *
      1.1
  )

print(
  cell_number_plot
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12C_GSE220243_cell_recovery.pdf"
    ),
  
  plot =
    cell_number_plot,
  
  width =
    9,
  
  height =
    5.5
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12C_GSE220243_cell_recovery.tiff"
    ),
  
  plot =
    cell_number_plot,
  
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
# 17. Raw MMP13 detection plot
# ------------------------------------------------------------

mmp13_detection_plot_data <-
  mmp13_detection_summary %>%
  
  dplyr::mutate(
    
    sample_label = factor(
      sample_label,
      levels =
        sample_table$sample_label
    )
  )

mmp13_detection_plot <- ggplot(
  
  mmp13_detection_plot_data,
  
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
    
    vjust = -0.35,
    size = 3
  ) +
  
  labs(
    
    title =
      "Exploratory raw MMP13 detection across cartilage donors",
    
    subtitle =
      paste(
        "Percentage of raw barcodes with at least one",
        "MMP13 count; final interpretation requires QC",
        "and cell-type annotation"
      ),
    
    x =
      "Cartilage donor",
    
    y =
      "MMP13-positive raw barcodes (%)",
    
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
        mmp13_detection_plot_data$
          MMP13_detection_percent,
        na.rm = TRUE
      ) *
      1.15
  )

print(
  mmp13_detection_plot
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12D_GSE220243_raw_MMP13_detection.pdf"
    ),
  
  plot =
    mmp13_detection_plot,
  
  width =
    9,
  
  height =
    5.5
)

ggsave(
  
  filename =
    file.path(
      supplementary_figure_folder,
      "FigureS12D_GSE220243_raw_MMP13_detection.tiff"
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
# 18. Group-level raw QC summary
# ------------------------------------------------------------

group_qc_summary <- sample_qc_summary %>%
  
  dplyr::group_by(
    group
  ) %>%
  
  dplyr::summarise(
    
    donors =
      dplyr::n(),
    
    total_cells_after_initial_creation =
      sum(
        cells_after_initial_creation
      ),
    
    median_cells_per_donor =
      stats::median(
        cells_after_initial_creation
      ),
    
    median_of_donor_median_features =
      stats::median(
        median_nFeature_RNA
      ),
    
    median_of_donor_median_counts =
      stats::median(
        median_nCount_RNA
      ),
    
    median_of_donor_median_percent_mt =
      stats::median(
        median_percent_mt
      ),
    
    .groups =
      "drop"
  )

cat("\nGroup-level raw QC summary:\n")

print(
  group_qc_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 19. Save tabular outputs
# ------------------------------------------------------------

utils::write.csv(
  
  sample_qc_summary,
  
  file =
    file.path(
      table_folder,
      "GSE220243_cartilage_raw_sample_QC_summary.csv"
    ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  mmp13_detection_summary,
  
  file =
    file.path(
      table_folder,
      "GSE220243_cartilage_raw_MMP13_detection.csv"
    ),
  
  row.names =
    FALSE
)

utils::write.csv(
  
  download_manifest,
  
  file =
    file.path(
      table_folder,
      "GSE220243_cartilage_download_manifest.csv"
    ),
  
  row.names =
    FALSE
)

# ------------------------------------------------------------
# 20. Save Excel workbook
# ------------------------------------------------------------

openxlsx::write.xlsx(
  
  list(
    
    Sample_Design =
      as.data.frame(
        sample_table
      ),
    
    Download_Manifest =
      as.data.frame(
        download_manifest
      ),
    
    Sample_QC =
      as.data.frame(
        sample_qc_summary
      ),
    
    Group_QC =
      as.data.frame(
        group_qc_summary
      ),
    
    Raw_MMP13_Detection =
      as.data.frame(
        mmp13_detection_summary
      ),
    
    Package_Versions =
      as.data.frame(
        package_versions
      )
  ),
  
  file =
    file.path(
      table_folder,
      "GSE220243_cartilage_raw_preparation_QC.xlsx"
    ),
  
  overwrite =
    TRUE
)

# ------------------------------------------------------------
# 21. Save session information
# ------------------------------------------------------------

sink(
  
  file.path(
    table_folder,
    "GSE220243_cartilage_raw_preparation_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 22. Final preparation summary
# ------------------------------------------------------------

total_cells <- ncol(
  merged_raw_object
)

total_genes <- nrow(
  merged_raw_object
)

total_normal_cells <- sum(
  merged_raw_object$group ==
    "Normal"
)

total_OA_cells <- sum(
  merged_raw_object$group ==
    "OA"
)

message("")
message("================================================")
message("GSE220243 CARTILAGE RAW PREPARATION SELESAI")
message("================================================")

message(
  "Cartilage donors analyzed      : ",
  nrow(sample_table)
)

message(
  "Normal donors                  : ",
  sum(sample_table$group == "Normal")
)

message(
  "OA donors                      : ",
  sum(sample_table$group == "OA")
)

message(
  "Merged genes                   : ",
  format(
    total_genes,
    big.mark = ","
  )
)

message(
  "Merged cells                   : ",
  format(
    total_cells,
    big.mark = ","
  )
)

message(
  "Normal cells                   : ",
  format(
    total_normal_cells,
    big.mark = ","
  )
)

message(
  "OA cells                       : ",
  format(
    total_OA_cells,
    big.mark = ","
  )
)

message(
  "Raw Seurat list                : ",
  file.path(
    processed_data_folder,
    "GSE220243_cartilage_raw_seurat_list.rds"
  )
)

message(
  "Merged raw Seurat object       : ",
  file.path(
    processed_data_folder,
    "GSE220243_cartilage_merged_raw_seurat.rds"
  )
)

message(
  "QC figures                     : ",
  supplementary_figure_folder
)

message(
  "Final cell filtering performed : FALSE"
)

message(
  "Cell-type annotation performed : FALSE"
)

message("================================================")