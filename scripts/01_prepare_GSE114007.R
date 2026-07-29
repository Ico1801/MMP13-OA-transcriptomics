# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 01_prepare_GSE114007.R
# PURPOSE : Download dan memeriksa raw counts GSE114007
# ============================================================

rm(list = ls())

library(GEOquery)
library(openxlsx)
library(tidyverse)

# ------------------------------------------------------------
# 1. Membuat folder bila belum tersedia
# ------------------------------------------------------------

dir.create(
  "data_raw",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "data_processed",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "metadata",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Mengunduh supplementary files GSE114007
# ------------------------------------------------------------

gse114007_download <- GEOquery::getGEOSuppFiles(
  GEO = "GSE114007",
  makeDirectory = TRUE,
  baseDir = "data_raw"
)

print(gse114007_download)
# ------------------------------------------------------------
# 3. Melihat file supplementary
# ------------------------------------------------------------

gse114007_folder <- file.path(
  "data_raw",
  "GSE114007"
)

downloaded_files <- list.files(
  path = gse114007_folder,
  full.names = TRUE
)

print(downloaded_files)
# ------------------------------------------------------------
# 4. Menemukan file raw counts
# ------------------------------------------------------------

raw_counts_file <- list.files(
  path = gse114007_folder,
  pattern = "GSE114007_raw_counts\\.xlsx$",
  full.names = TRUE
)

print(raw_counts_file)

if (length(raw_counts_file) == 0) {
  stop("File GSE114007_raw_counts.xlsx tidak ditemukan.")
}

if (length(raw_counts_file) > 1) {
  stop("Ditemukan lebih dari satu file raw counts.")
}

# ------------------------------------------------------------
# 5. Memeriksa nama sheet Excel
# ------------------------------------------------------------

sheet_names <- openxlsx::getSheetNames(
  raw_counts_file
)

print(sheet_names)
# ------------------------------------------------------------
# 6. Membaca raw-count matrix
# ------------------------------------------------------------

gse114007_raw <- openxlsx::read.xlsx(
  xlsxFile = raw_counts_file,
  sheet = 1,
  check.names = FALSE
)

# Periksa dimensi
print(dim(gse114007_raw))

# Periksa nama kolom pertama
print(
  colnames(gse114007_raw)[
    seq_len(min(15, ncol(gse114007_raw)))
  ]
)

# Melihat bagian awal tabel
View(
  gse114007_raw[
    1:min(20, nrow(gse114007_raw)),
    1:min(15, ncol(gse114007_raw))
  ]
)
# ------------------------------------------------------------
# 7. Menyimpan raw counts
# ------------------------------------------------------------

saveRDS(
  gse114007_raw,
  file = "data_processed/GSE114007_raw_counts_original.rds"
)

openxlsx::write.xlsx(
  data.frame(
    variable = c(
      "number_of_rows",
      "number_of_columns"
    ),
    value = c(
      nrow(gse114007_raw),
      ncol(gse114007_raw)
    )
  ),
  file = "metadata/GSE114007_raw_dimension.xlsx",
  overwrite = TRUE
)

message(
  "GSE114007 raw counts berhasil dibaca dan disimpan."
)
# ============================================================
# 6. Membaca kedua sheet: Normal dan OA
# ============================================================

normal_raw <- openxlsx::read.xlsx(
  xlsxFile = raw_counts_file,
  sheet = "Normal",
  check.names = FALSE
)

oa_raw <- openxlsx::read.xlsx(
  xlsxFile = raw_counts_file,
  sheet = "OA",
  check.names = FALSE
)

message("Sheet Normal dan OA berhasil dibaca.")

# Periksa dimensi kedua sheet
cat("\nDimensi sheet Normal:\n")
print(dim(normal_raw))

cat("\nDimensi sheet OA:\n")
print(dim(oa_raw))

# Periksa nama kolom
cat("\nNama kolom Normal:\n")
print(colnames(normal_raw))

cat("\nNama kolom OA:\n")
print(colnames(oa_raw))
# ============================================================
# 7. Memeriksa kolom gene symbol
# ============================================================

cat("\nKolom pertama sheet Normal:\n")
print(colnames(normal_raw)[1])

cat("\nKolom pertama sheet OA:\n")
print(colnames(oa_raw)[1])

if (colnames(normal_raw)[1] != "symbol") {
  stop("Kolom pertama sheet Normal bukan 'symbol'.")
}

if (colnames(oa_raw)[1] != "symbol") {
  stop("Kolom pertama sheet OA bukan 'symbol'.")
}

# Periksa apakah jumlah baris sama
cat("\nJumlah baris sama:\n")
print(nrow(normal_raw) == nrow(oa_raw))

# Periksa apakah urutan gene symbol sama persis
symbols_identical <- identical(
  as.character(normal_raw$symbol),
  as.character(oa_raw$symbol)
)

cat("\nUrutan gene symbol identik:\n")
print(symbols_identical)
# ============================================================
# 8. Pemeriksaan gene symbol kosong dan duplikat
# ============================================================

normal_symbols <- trimws(
  as.character(normal_raw$symbol)
)

oa_symbols <- trimws(
  as.character(oa_raw$symbol)
)

symbol_qc <- data.frame(
  check = c(
    "Normal missing symbols",
    "Normal blank symbols",
    "Normal duplicated symbols",
    "OA missing symbols",
    "OA blank symbols",
    "OA duplicated symbols"
  ),
  value = c(
    sum(is.na(normal_symbols)),
    sum(normal_symbols == "", na.rm = TRUE),
    sum(duplicated(normal_symbols)),
    sum(is.na(oa_symbols)),
    sum(oa_symbols == "", na.rm = TRUE),
    sum(duplicated(oa_symbols))
  )
)

print(symbol_qc)

openxlsx::write.xlsx(
  symbol_qc,
  file = "metadata/GSE114007_symbol_QC.xlsx",
  overwrite = TRUE
)
# ============================================================
# 9. Menggabungkan Normal dan OA
# ============================================================

if (!symbols_identical) {
  stop(
    "Urutan gene symbol Normal dan OA tidak sama. ",
    "Penggabungan dihentikan."
  )
}

gse114007_combined <- data.frame(
  symbol = normal_symbols,
  normal_raw[, -1, drop = FALSE],
  oa_raw[, -1, drop = FALSE],
  check.names = FALSE
)

cat("\nDimensi sebelum membersihkan simbol:\n")
print(dim(gse114007_combined))

cat("\nContoh nama sampel gabungan:\n")
print(colnames(gse114007_combined))
# ============================================================
# 10. Membuat sample metadata
# ============================================================

normal_sample_names <- colnames(normal_raw)[-1]
oa_sample_names <- colnames(oa_raw)[-1]

sample_metadata <- tibble::tibble(
  sample_id = c(
    normal_sample_names,
    oa_sample_names
  ),
  group = factor(
    c(
      rep("Control", length(normal_sample_names)),
      rep("OA", length(oa_sample_names))
    ),
    levels = c("Control", "OA")
  )
)

print(sample_metadata)

cat("\nJumlah sampel setiap kelompok:\n")
print(table(sample_metadata$group))

# Pastikan urutannya sama dengan count matrix
sample_order_identical <- identical(
  sample_metadata$sample_id,
  colnames(gse114007_combined)[-1]
)

cat("\nUrutan metadata dan count matrix sama:\n")
print(sample_order_identical)

if (!sample_order_identical) {
  stop(
    "Urutan sampel pada metadata dan count matrix tidak sama."
  )
}
# ============================================================
# 11. Memastikan semua expression values numerik
# ============================================================

gse114007_combined <- gse114007_combined %>%
  dplyr::mutate(
    dplyr::across(
      -symbol,
      as.numeric
    )
  )

count_matrix_temporary <- as.matrix(
  gse114007_combined[, -1, drop = FALSE]
)

cat("\nAda nilai NA dalam counts:\n")
print(anyNA(count_matrix_temporary))

cat("\nAda nilai negatif:\n")
print(any(count_matrix_temporary < 0, na.rm = TRUE))

cat("\nSemua nilai merupakan bilangan bulat:\n")
print(
  all(
    count_matrix_temporary ==
      round(count_matrix_temporary),
    na.rm = TRUE
  )
)

cat("\nRentang raw counts:\n")
print(range(count_matrix_temporary, na.rm = TRUE))
# ============================================================
# 12. Membersihkan gene symbol dan collapse duplikat
# ============================================================

gse114007_clean <- gse114007_combined %>%
  dplyr::filter(
    !is.na(symbol),
    symbol != ""
  ) %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

cat("\nDimensi sebelum collapse duplicate symbols:\n")
print(dim(gse114007_combined))

cat("\nDimensi setelah collapse duplicate symbols:\n")
print(dim(gse114007_clean))

cat("\nJumlah simbol duplikat setelah dibersihkan:\n")
print(sum(duplicated(gse114007_clean$symbol)))
# ============================================================
# 13. Membuat count matrix final
# ============================================================

count_matrix <- as.matrix(
  gse114007_clean[, -1, drop = FALSE]
)

rownames(count_matrix) <- gse114007_clean$symbol

storage.mode(count_matrix) <- "integer"

cat("\nDimensi count matrix final:\n")
print(dim(count_matrix))

cat("\nLima gene pertama:\n")
print(rownames(count_matrix)[1:5])

cat("\nJumlah sampel berdasarkan metadata:\n")
print(table(sample_metadata$group))

# Pastikan sampel sama dan berurutan
stopifnot(
  identical(
    colnames(count_matrix),
    sample_metadata$sample_id
  )
)
# ============================================================
# 14. Memeriksa target MMP13
# ============================================================

mmp13_present <- "MMP13" %in% rownames(count_matrix)

cat("\nApakah MMP13 tersedia:\n")
print(mmp13_present)

if (mmp13_present) {
  
  mmp13_counts <- data.frame(
    sample_id = colnames(count_matrix),
    group = sample_metadata$group,
    MMP13_raw_count = as.numeric(
      count_matrix["MMP13", ]
    )
  )
  
  print(mmp13_counts)
  
  openxlsx::write.xlsx(
    mmp13_counts,
    file = "metadata/GSE114007_MMP13_raw_counts.xlsx",
    overwrite = TRUE
  )
  
} else {
  
  warning(
    "MMP13 tidak ditemukan persis sebagai row name."
  )
}
# ============================================================
# 15. Menyimpan data final GSE114007
# ============================================================

saveRDS(
  count_matrix,
  file = "data_processed/GSE114007_raw_count_matrix.rds"
)

saveRDS(
  sample_metadata,
  file = "data_processed/GSE114007_sample_metadata.rds"
)

openxlsx::write.xlsx(
  sample_metadata,
  file = "metadata/GSE114007_sample_metadata.xlsx",
  overwrite = TRUE
)

gse114007_summary <- data.frame(
  variable = c(
    "Number of genes",
    "Number of control samples",
    "Number of OA samples",
    "Total samples",
    "MMP13 available"
  ),
  value = c(
    nrow(count_matrix),
    sum(sample_metadata$group == "Control"),
    sum(sample_metadata$group == "OA"),
    ncol(count_matrix),
    mmp13_present
  )
)

openxlsx::write.xlsx(
  gse114007_summary,
  file = "metadata/GSE114007_final_summary.xlsx",
  overwrite = TRUE
)

message(
  "GSE114007 Normal dan OA berhasil digabungkan."
)