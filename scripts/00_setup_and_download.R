# ============================================================
# PROJECT: OA MMP13 Target Discovery
# SCRIPT : 00_setup_and_download.R
# PURPOSE: Install packages, create folders, download GEO data
# ============================================================

# Bersihkan environment
rm(list = ls())

# ------------------------------------------------------------
# 1. Instal BiocManager
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ------------------------------------------------------------
# 2. Paket CRAN
# ------------------------------------------------------------

cran_packages <- c(
  "tidyverse",
  "here",
  "janitor",
  "metafor",
  "openxlsx",
  "patchwork"
)

cran_missing <- cran_packages[
  !cran_packages %in% rownames(installed.packages())
]

if (length(cran_missing) > 0) {
  install.packages(cran_missing)
}

# ------------------------------------------------------------
# 3. Paket Bioconductor
# ------------------------------------------------------------

bioc_packages <- c(
  "GEOquery",
  "limma",
  "edgeR",
  "Biobase",
  "SummarizedExperiment",
  "sva",
  "clusterProfiler",
  "org.Hs.eg.db",
  "GSVA"
)

bioc_missing <- bioc_packages[
  !bioc_packages %in% rownames(installed.packages())
]

if (length(bioc_missing) > 0) {
  BiocManager::install(
    bioc_missing,
    ask = FALSE,
    update = FALSE
  )
}

# ------------------------------------------------------------
# 4. Muat paket
# ------------------------------------------------------------

library(GEOquery)
library(tidyverse)
library(here)
library(janitor)
library(Biobase)
library(SummarizedExperiment)
library(openxlsx)
library(clusterProfiler)
library(enrichplot)

message("Semua paket utama berhasil dimuat.")
# ------------------------------------------------------------
# 5. Membuat struktur folder project
# ------------------------------------------------------------

project_folders <- c(
  "data_raw",
  "data_processed",
  "metadata",
  "results",
  "results/tables",
  "results/figures",
  "scripts"
)

for (folder in project_folders) {
  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

message("Struktur folder berhasil dibuat.")

list.files()
# ------------------------------------------------------------
# ------------------------------------------------------------
# 6. Daftar dataset bulk transcriptomics
# ------------------------------------------------------------

bulk_ids <- c(
  "GSE114007",
  "GSE117999",
  "GSE57218",
  "GSE169077"
)

# ------------------------------------------------------------
# 7. Pastikan GEOquery tersedia
# ------------------------------------------------------------

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  stop(
    "Paket GEOquery belum terinstal. ",
    "Jalankan BiocManager::install('GEOquery') terlebih dahulu."
  )
}

message("GEOquery tersedia dan siap digunakan.")

# ------------------------------------------------------------
# 8. Fungsi download GEO
# ------------------------------------------------------------

download_geo_dataset <- function(gse_id) {
  
  message("")
  message("=========================================")
  message("Mengunduh dataset: ", gse_id)
  message("=========================================")
  
  result <- tryCatch(
    {
      GEOquery::getGEO(
        GEO = gse_id,
        GSEMatrix = TRUE,
        destdir = "data_raw"
      )
    },
    error = function(e) {
      message("DOWNLOAD GAGAL: ", gse_id)
      message("Penyebab: ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (!is.null(result)) {
    message(
      "DOWNLOAD BERHASIL: ",
      gse_id,
      " | jumlah platform/objek = ",
      length(result)
    )
  }
  
  return(result)
}

# ------------------------------------------------------------
# 9. Jalankan download
# ------------------------------------------------------------

bulk_geo <- lapply(
  bulk_ids,
  download_geo_dataset
)

names(bulk_geo) <- bulk_ids

# Catat dataset yang berhasil dan gagal
download_status <- tibble::tibble(
  dataset = bulk_ids,
  status = ifelse(
    vapply(bulk_geo, is.null, logical(1)),
    "Failed",
    "Successful"
  )
)

print(download_status)

# Pertahankan hanya dataset yang berhasil
bulk_geo <- bulk_geo[
  !vapply(bulk_geo, is.null, logical(1))
]

# Jangan menyimpan bila semuanya gagal
if (length(bulk_geo) == 0) {
  
  stop(
    "Tidak ada dataset yang berhasil diunduh. ",
    "Jangan lanjut ke tahap berikutnya."
  )
  
} else {
  
  saveRDS(
    bulk_geo,
    file = "data_processed/bulk_geo_original_objects.rds"
  )
  
  message(
    length(bulk_geo),
    " dari ",
    length(bulk_ids),
    " dataset berhasil diunduh dan disimpan."
  )
}
# ------------------------------------------------------------
# 9. Fungsi mengambil metadata
# ------------------------------------------------------------

get_sample_metadata <- function(object) {
  
  if (inherits(object, "ExpressionSet")) {
    metadata <- Biobase::pData(object)
    
  } else if (
    inherits(object, "SummarizedExperiment") ||
    inherits(object, "RangedSummarizedExperiment")
  ) {
    metadata <- as.data.frame(
      SummarizedExperiment::colData(object)
    )
    
  } else {
    stop(
      "Jenis objek tidak dikenali: ",
      paste(class(object), collapse = ", ")
    )
  }
  
  metadata
}

# ------------------------------------------------------------
# 10. Fungsi mengambil expression matrix
# ------------------------------------------------------------

get_expression_matrix <- function(object) {
  
  if (inherits(object, "ExpressionSet")) {
    expression_matrix <- Biobase::exprs(object)
    
  } else if (
    inherits(object, "SummarizedExperiment") ||
    inherits(object, "RangedSummarizedExperiment")
  ) {
    expression_matrix <- SummarizedExperiment::assay(object)
    
  } else {
    stop(
      "Jenis objek tidak dikenali: ",
      paste(class(object), collapse = ", ")
    )
  }
  
  expression_matrix
}

# ------------------------------------------------------------
# 11. Fungsi mengambil anotasi feature
# ------------------------------------------------------------

get_feature_annotation <- function(object) {
  
  if (inherits(object, "ExpressionSet")) {
    feature_data <- Biobase::fData(object)
    
  } else if (
    inherits(object, "SummarizedExperiment") ||
    inherits(object, "RangedSummarizedExperiment")
  ) {
    feature_data <- as.data.frame(
      SummarizedExperiment::rowData(object)
    )
    
  } else {
    stop(
      "Jenis objek tidak dikenali: ",
      paste(class(object), collapse = ", ")
    )
  }
  
  feature_data
}
# ------------------------------------------------------------
# 12. Ringkasan platform setiap dataset
# ------------------------------------------------------------

platform_summary <- tibble::tibble(
  dataset = names(bulk_geo),
  number_of_platforms = vapply(
    bulk_geo,
    length,
    integer(1)
  )
)

print(platform_summary)

openxlsx::write.xlsx(
  platform_summary,
  file = "metadata/platform_summary.xlsx",
  overwrite = TRUE
)
# ------------------------------------------------------------
# 13. Ringkasan ukuran setiap objek/platform
# ------------------------------------------------------------

dataset_dimension_summary <- list()

for (gse_id in names(bulk_geo)) {
  
  current_platforms <- bulk_geo[[gse_id]]
  
  for (platform_number in seq_along(current_platforms)) {
    
    current_object <- current_platforms[[platform_number]]
    current_expression <- get_expression_matrix(current_object)
    current_metadata <- get_sample_metadata(current_object)
    
    dataset_dimension_summary[[
      paste0(gse_id, "_platform_", platform_number)
    ]] <- tibble::tibble(
      dataset = gse_id,
      platform_number = platform_number,
      object_class = paste(
        class(current_object),
        collapse = "; "
      ),
      number_of_features = nrow(current_expression),
      number_of_samples_expression = ncol(current_expression),
      number_of_samples_metadata = nrow(current_metadata)
    )
  }
}

dataset_dimension_summary <- dplyr::bind_rows(
  dataset_dimension_summary
)

print(dataset_dimension_summary)

openxlsx::write.xlsx(
  dataset_dimension_summary,
  file = "metadata/dataset_dimension_summary.xlsx",
  overwrite = TRUE
)
# ------------------------------------------------------------
# 14. Fungsi memilih kolom metadata penting
# ------------------------------------------------------------

inspect_metadata <- function(object) {
  
  metadata <- get_sample_metadata(object) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id")
  
  important_columns <- grep(
    pattern = paste(
      c(
        "sample_id",
        "geo_accession",
        "title",
        "source_name",
        "characteristics",
        "description"
      ),
      collapse = "|"
    ),
    x = colnames(metadata),
    value = TRUE,
    ignore.case = TRUE
  )
  
  metadata %>%
    select(any_of(unique(important_columns)))
}