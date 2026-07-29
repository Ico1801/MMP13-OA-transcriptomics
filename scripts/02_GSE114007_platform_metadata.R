# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 02_GSE114007_platform_metadata.R
# PURPOSE : Menambahkan platform sequencing ke metadata
# ============================================================

rm(list = ls())

library(Biobase)
library(dplyr)
library(tibble)
library(purrr)
library(openxlsx)

# ------------------------------------------------------------
# 1. Memuat data yang sudah disiapkan
# ------------------------------------------------------------

count_matrix <- readRDS(
  "data_processed/GSE114007_raw_count_matrix.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE114007_sample_metadata.rds"
)

bulk_geo <- readRDS(
  "data_processed/bulk_geo_original_objects.rds"
)

gse114007_objects <- bulk_geo[["GSE114007"]]

cat("Jumlah objek/platform GSE114007:\n")
print(length(gse114007_objects))
# ------------------------------------------------------------
# 2. Mengambil sample ID dan platform dari metadata GEO
# ------------------------------------------------------------

platform_map <- purrr::map_dfr(
  gse114007_objects,
  function(current_object) {
    
    current_metadata <- Biobase::pData(
      current_object
    )
    
    # Mengambil platform ID
    if ("platform_id" %in% colnames(current_metadata)) {
      
      current_platform <- as.character(
        current_metadata$platform_id
      )
      
    } else {
      
      current_platform <- rep(
        Biobase::annotation(current_object),
        nrow(current_metadata)
      )
    }
    
    # Mengambil GEO accession
    if ("geo_accession" %in% colnames(current_metadata)) {
      
      current_geo_accession <- as.character(
        current_metadata$geo_accession
      )
      
    } else {
      
      current_geo_accession <- rownames(
        current_metadata
      )
    }
    
    # Mengambil nama sampel
    tibble::tibble(
      sample_id = as.character(
        current_metadata$title
      ),
      geo_accession = current_geo_accession,
      platform_id = current_platform
    )
  }
)

platform_map <- platform_map %>%
  dplyr::distinct(
    sample_id,
    .keep_all = TRUE
  )

print(platform_map)
# ------------------------------------------------------------
# 3. Menggabungkan platform dengan metadata sampel
# ------------------------------------------------------------

sample_metadata_platform <- sample_metadata %>%
  dplyr::left_join(
    platform_map,
    by = "sample_id"
  ) %>%
  dplyr::mutate(
    platform_name = dplyr::case_when(
      platform_id == "GPL11154" ~ "Illumina HiSeq 2000",
      platform_id == "GPL18573" ~ "Illumina NextSeq 500",
      TRUE ~ NA_character_
    ),
    
    group = factor(
      group,
      levels = c("Control", "OA")
    ),
    
    platform_id = factor(
      platform_id,
      levels = c("GPL11154", "GPL18573")
    )
  )

print(sample_metadata_platform)
# ------------------------------------------------------------
# 4. Quality control metadata
# ------------------------------------------------------------

cat("\nJumlah baris metadata:\n")
print(nrow(sample_metadata_platform))

cat("\nJumlah platform yang kosong:\n")
print(
  sum(is.na(sample_metadata_platform$platform_id))
)

cat("\nJumlah GEO accession yang kosong:\n")
print(
  sum(is.na(sample_metadata_platform$geo_accession))
)

cat("\nDistribusi kelompok berdasarkan platform:\n")
print(
  table(
    sample_metadata_platform$platform_id,
    sample_metadata_platform$group,
    useNA = "ifany"
  )
)

cat("\nUrutan metadata sesuai count matrix:\n")
print(
  identical(
    sample_metadata_platform$sample_id,
    colnames(count_matrix)
  )
)
# ------------------------------------------------------------
# 5. Validasi final sebelum menyimpan
# ------------------------------------------------------------

stopifnot(
  nrow(sample_metadata_platform) == ncol(count_matrix)
)

stopifnot(
  sum(is.na(sample_metadata_platform$platform_id)) == 0
)

stopifnot(
  sum(is.na(sample_metadata_platform$geo_accession)) == 0
)

stopifnot(
  identical(
    sample_metadata_platform$sample_id,
    colnames(count_matrix)
  )
)

message("Semua pemeriksaan metadata berhasil.")

# ------------------------------------------------------------
# 6. Menyimpan metadata final
# ------------------------------------------------------------

saveRDS(
  sample_metadata_platform,
  file =
    "data_processed/GSE114007_sample_metadata_with_platform.rds"
)

openxlsx::write.xlsx(
  as.data.frame(sample_metadata_platform),
  file =
    "metadata/GSE114007_sample_metadata_with_platform.xlsx",
  overwrite = TRUE
)

message(
  "Metadata GSE114007 dengan platform berhasil disimpan."
)