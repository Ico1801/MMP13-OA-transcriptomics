# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 05_prepare_GSE117999.R
# PURPOSE : Menyiapkan expression matrix, metadata sampel,
#           eksklusi sampel, dan anotasi probe GSE117999
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Memuat paket
# ------------------------------------------------------------

required_packages <- c(
  "GEOquery",
  "Biobase",
  "limma",
  "dplyr",
  "tibble",
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

library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(tibble)
library(stringr)
library(openxlsx)

message("Semua paket GSE117999 berhasil dimuat.")

# ------------------------------------------------------------
# 2. Memastikan folder tersedia
# ------------------------------------------------------------

required_folders <- c(
  "data_processed",
  "metadata",
  "results",
  "results/tables",
  "results/figures"
)

for (folder in required_folders) {
  dir.create(
    folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------
# 3. Memuat objek GEO yang sudah diunduh
# ------------------------------------------------------------

bulk_geo_file <-
  "data_processed/bulk_geo_original_objects.rds"

if (!file.exists(bulk_geo_file)) {
  stop(
    "File tidak ditemukan: ",
    bulk_geo_file,
    "\nJalankan script 00_setup_and_download.R terlebih dahulu."
  )
}

bulk_geo <- readRDS(
  bulk_geo_file
)

if (!"GSE117999" %in% names(bulk_geo)) {
  stop(
    "GSE117999 tidak tersedia di bulk_geo_original_objects.rds."
  )
}

gse117999_list <- bulk_geo[["GSE117999"]]

cat("\nJumlah platform GSE117999:\n")
print(length(gse117999_list))

if (length(gse117999_list) != 1) {
  warning(
    "Jumlah platform bukan satu. ",
    "Periksa objek sebelum melanjutkan."
  )
}

gse117999_eset <- gse117999_list[[1]]

cat("\nClass objek:\n")
print(class(gse117999_eset))

cat("\nPlatform annotation:\n")
print(Biobase::annotation(gse117999_eset))

if (!inherits(gse117999_eset, "ExpressionSet")) {
  stop(
    "Objek GSE117999 bukan ExpressionSet."
  )
}

# ------------------------------------------------------------
# 4. Mengambil expression matrix, metadata, dan feature data
# ------------------------------------------------------------

expression_full <- Biobase::exprs(
  gse117999_eset
)

metadata_full <- Biobase::pData(
  gse117999_eset
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "expression_column"
  )

feature_annotation <- Biobase::fData(
  gse117999_eset
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "probe_id"
  )

cat("\nDimensi expression matrix lengkap:\n")
print(dim(expression_full))

cat("\nDimensi metadata lengkap:\n")
print(dim(metadata_full))

cat("\nDimensi feature annotation:\n")
print(dim(feature_annotation))

stopifnot(
  ncol(expression_full) == nrow(metadata_full)
)

stopifnot(
  identical(
    colnames(expression_full),
    metadata_full$expression_column
  )
)

message(
  "Expression matrix dan metadata memiliki urutan yang sesuai."
)

# ------------------------------------------------------------
# 5. Memeriksa metadata penting
# ------------------------------------------------------------

cat("\nNama seluruh kolom metadata:\n")
print(colnames(metadata_full))

important_metadata_columns <- grep(
  pattern = paste(
    c(
      "expression_column",
      "title",
      "geo_accession",
      "source_name",
      "characteristics",
      "description",
      "age",
      "bmi",
      "sex",
      "tissue"
    ),
    collapse = "|"
  ),
  x = colnames(metadata_full),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nKolom metadata penting:\n")
print(important_metadata_columns)

metadata_inspection <- metadata_full %>%
  dplyr::select(
    dplyr::any_of(
      unique(important_metadata_columns)
    )
  ) %>%
  tibble::as_tibble()

cat("\nMetadata penting:\n")

print(
  metadata_inspection,
  n = Inf,
  width = Inf
)

# Tampilkan seperti spreadsheet di RStudio
View(metadata_inspection)

# ------------------------------------------------------------
# 6. Menggabungkan teks metadata untuk ekstraksi informasi
# ------------------------------------------------------------

metadata_text_columns <- intersect(
  c(
    "title",
    "source_name_ch1",
    "characteristics_ch1",
    "characteristics_ch1.1",
    "characteristics_ch1.2",
    "characteristics_ch1.3",
    "description",
    "tissue:ch1"
  ),
  colnames(metadata_full)
)

if (length(metadata_text_columns) == 0) {
  stop(
    "Tidak ditemukan kolom metadata untuk identifikasi sampel."
  )
}

metadata_combined_text <- apply(
  metadata_full[
    ,
    metadata_text_columns,
    drop = FALSE
  ],
  MARGIN = 1,
  FUN = function(x) {
    paste(
      as.character(x),
      collapse = " | "
    )
  }
)

# ------------------------------------------------------------
# 7. Mengekstrak kode subjek dan kelompok
# ------------------------------------------------------------

metadata_upper <- stringr::str_to_upper(
  metadata_combined_text
)

metadata_lower <- stringr::str_to_lower(
  metadata_combined_text
)

# Cari kode seperti P4-001, P4_001, atau P4 001
subject_code_direct <- stringr::str_extract(
  metadata_upper,
  "P4[-_ ]?\\d{3}"
)

subject_code_direct <- subject_code_direct %>%
  stringr::str_replace_all(
    pattern = "[_ ]",
    replacement = "-"
  )

# Bila P4-xxx tidak ditemukan, ambil tiga angka pertama dari title
fallback_sample_code <- stringr::str_extract(
  as.character(metadata_full$title),
  "(?<!\\d)\\d{3}(?!\\d)"
)

fallback_subject_code <- dplyr::if_else(
  !is.na(fallback_sample_code),
  paste0(
    "P4-",
    fallback_sample_code
  ),
  NA_character_
)

subject_code_final <- dplyr::coalesce(
  subject_code_direct,
  fallback_subject_code
)

# Penentuan group
group_final <- dplyr::case_when(
  
  stringr::str_detect(
    metadata_lower,
    paste(
      c(
        "arthroscopic partial meniscectomy",
        "non[- ]?oa",
        "non osteoarthritic",
        "without osteoarthritis",
        "control cartilage",
        "normal cartilage"
      ),
      collapse = "|"
    )
  ) ~ "Control",
  
  stringr::str_detect(
    metadata_lower,
    paste(
      c(
        "end[- ]?stage osteoarthritis",
        "total knee arthroplasty",
        "osteoarthritis",
        "oa cartilage"
      ),
      collapse = "|"
    )
  ) ~ "OA",
  
  TRUE ~ NA_character_
)

sample_metadata <- metadata_full %>%
  dplyr::mutate(
    title_original = as.character(title),
    
    combined_metadata_text =
      metadata_combined_text,
    
    subject_code =
      subject_code_final,
    
    group =
      group_final
  )

cat("\nMetadata setelah ekstraksi kode dan kelompok:\n")

print(
  sample_metadata %>%
    dplyr::select(
      expression_column,
      geo_accession,
      title_original,
      subject_code,
      group
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 8. Validasi ekstraksi metadata
# ------------------------------------------------------------

cat("\nJumlah subject code yang kosong:\n")
print(
  sum(
    is.na(sample_metadata$subject_code)
  )
)

cat("\nJumlah group yang kosong:\n")
print(
  sum(
    is.na(sample_metadata$group)
  )
)

cat("\nJumlah sampel setiap group sebelum eksklusi:\n")
print(
  table(
    sample_metadata$group,
    useNA = "ifany"
  )
)

if (anyNA(sample_metadata$subject_code)) {
  
  cat(
    "\nSampel dengan subject code kosong:\n"
  )
  
  print(
    sample_metadata %>%
      dplyr::filter(
        is.na(subject_code)
      ) %>%
      dplyr::select(
        expression_column,
        geo_accession,
        title_original,
        combined_metadata_text
      ) %>%
      tibble::as_tibble(),
    n = Inf,
    width = Inf
  )
  
  stop(
    "Ada subject code yang belum berhasil diidentifikasi."
  )
}

if (anyNA(sample_metadata$group)) {
  
  cat(
    "\nSampel dengan group kosong:\n"
  )
  
  print(
    sample_metadata %>%
      dplyr::filter(
        is.na(group)
      ) %>%
      dplyr::select(
        expression_column,
        geo_accession,
        title_original,
        combined_metadata_text
      ) %>%
      tibble::as_tibble(),
    n = Inf,
    width = Inf
  )
  
  stop(
    "Ada kelompok sampel yang belum berhasil diidentifikasi."
  )
}

# ------------------------------------------------------------
# 9. Menandai sampel yang dikeluarkan
# ------------------------------------------------------------

excluded_subjects <- c(
  "P4-010",
  "P4-012",
  "P4-104",
  "P4-108"
)

sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    include_in_primary_analysis =
      !subject_code %in% excluded_subjects,
    
    exclusion_status =
      dplyr::if_else(
        include_in_primary_analysis,
        "Included",
        "Excluded according to study design"
      )
  )

excluded_samples_found <- sample_metadata %>%
  dplyr::filter(
    !include_in_primary_analysis
  )

cat("\nSampel yang dikeluarkan:\n")

print(
  excluded_samples_found %>%
    dplyr::select(
      expression_column,
      geo_accession,
      title_original,
      subject_code,
      group,
      exclusion_status
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

cat("\nJumlah sampel eksklusi yang ditemukan:\n")
print(nrow(excluded_samples_found))

missing_excluded_subjects <- setdiff(
  excluded_subjects,
  sample_metadata$subject_code
)

if (length(missing_excluded_subjects) > 0) {
  
  warning(
    "Subjek eksklusi berikut tidak ditemukan: ",
    paste(
      missing_excluded_subjects,
      collapse = ", "
    )
  )
}

if (nrow(excluded_samples_found) != 4) {
  
  stop(
    "Jumlah sampel yang dikeluarkan bukan empat. ",
    "Periksa subject_code dan metadata sebelum melanjutkan."
  )
}

# ------------------------------------------------------------
# 10. Membuat metadata final
# ------------------------------------------------------------

sample_metadata_included <- sample_metadata %>%
  dplyr::filter(
    include_in_primary_analysis
  ) %>%
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "Control",
        "OA"
      )
    )
  ) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(
    group_number =
      dplyr::row_number(),
    
    plot_label =
      dplyr::if_else(
        group == "Control",
        paste0(
          "C",
          group_number
        ),
        paste0(
          "OA",
          group_number
        )
      )
  ) %>%
  dplyr::ungroup()

cat("\nMetadata final setelah eksklusi:\n")

print(
  sample_metadata_included %>%
    dplyr::select(
      expression_column,
      geo_accession,
      subject_code,
      group,
      plot_label
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

cat("\nJumlah sampel setiap group setelah eksklusi:\n")

print(
  table(
    sample_metadata_included$group
  )
)

cat("\nTotal sampel setelah eksklusi:\n")
print(nrow(sample_metadata_included))

if (
  sum(
    sample_metadata_included$group ==
    "Control"
  ) != 10
) {
  stop(
    "Jumlah Control setelah eksklusi bukan 10."
  )
}

if (
  sum(
    sample_metadata_included$group ==
    "OA"
  ) != 10
) {
  stop(
    "Jumlah OA setelah eksklusi bukan 10."
  )
}

# ------------------------------------------------------------
# 11. Membuat expression matrix final
# ------------------------------------------------------------

included_sample_columns <-
  sample_metadata_included$expression_column

if (
  !all(
    included_sample_columns %in%
    colnames(expression_full)
  )
) {
  stop(
    "Ada sample ID metadata yang tidak ditemukan ",
    "pada expression matrix."
  )
}

expression_included <- expression_full[
  ,
  included_sample_columns,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(expression_included),
    sample_metadata_included$expression_column
  )
)

cat("\nDimensi expression matrix setelah eksklusi:\n")
print(dim(expression_included))

stopifnot(
  ncol(expression_included) == 20
)

# ------------------------------------------------------------
# 12. Pemeriksaan nilai expression matrix
# ------------------------------------------------------------

cat("\nAda NA pada expression matrix lengkap:\n")
print(anyNA(expression_full))

cat("\nAda NA setelah eksklusi:\n")
print(anyNA(expression_included))

cat("\nAda nilai non-finite setelah eksklusi:\n")
print(
  any(
    !is.finite(
      expression_included
    )
  )
)

cat("\nRentang nilai setelah eksklusi:\n")

expression_range <- range(
  expression_included,
  na.rm = TRUE
)

print(expression_range)

cat("\nQuantile expression values:\n")

expression_quantiles <- quantile(
  expression_included,
  probabilities = c(
    0,
    0.01,
    0.25,
    0.50,
    0.75,
    0.99,
    1
  ),
  na.rm = TRUE
)

print(expression_quantiles)

if (
  anyNA(expression_included) ||
  any(
    !is.finite(
      expression_included
    )
  )
) {
  stop(
    "Expression matrix mengandung NA atau nilai non-finite."
  )
}

# ------------------------------------------------------------
# 13. Menilai apakah data kemungkinan sudah log2
# ------------------------------------------------------------

likely_log2_scaled <- (
  expression_quantiles[
    "100%"
  ] < 50
)

cat("\nData kemungkinan sudah dalam skala log2:\n")
print(likely_log2_scaled)

if (!likely_log2_scaled) {
  
  warning(
    "Nilai ekspresi tampaknya belum dalam skala log2. ",
    "Jangan lanjut ke limma sebelum preprocessing diperiksa."
  )
}

# ------------------------------------------------------------
# 14. Ringkasan expression setiap sampel
# ------------------------------------------------------------

sample_expression_summary <- tibble::tibble(
  
  expression_column =
    colnames(
      expression_included
    ),
  
  minimum =
    apply(
      expression_included,
      2,
      min,
      na.rm = TRUE
    ),
  
  q1 =
    apply(
      expression_included,
      2,
      stats::quantile,
      probs = 0.25,
      na.rm = TRUE
    ),
  
  median =
    apply(
      expression_included,
      2,
      stats::median,
      na.rm = TRUE
    ),
  
  mean =
    colMeans(
      expression_included,
      na.rm = TRUE
    ),
  
  q3 =
    apply(
      expression_included,
      2,
      stats::quantile,
      probs = 0.75,
      na.rm = TRUE
    ),
  
  maximum =
    apply(
      expression_included,
      2,
      max,
      na.rm = TRUE
    )
) %>%
  dplyr::left_join(
    sample_metadata_included %>%
      dplyr::select(
        expression_column,
        geo_accession,
        subject_code,
        group,
        plot_label
      ),
    by = "expression_column"
  )

cat("\nRingkasan ekspresi per sampel:\n")

print(
  sample_expression_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 15. Memeriksa feature annotation
# ------------------------------------------------------------

cat("\nNama seluruh kolom feature annotation:\n")
print(colnames(feature_annotation))

candidate_annotation_columns <- grep(
  pattern = paste(
    c(
      "symbol",
      "gene",
      "entrez",
      "accession",
      "description",
      "systematic"
    ),
    collapse = "|"
  ),
  x = colnames(feature_annotation),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nKolom kandidat anotasi:\n")
print(candidate_annotation_columns)

if (length(candidate_annotation_columns) > 0) {
  
  annotation_inspection <- feature_annotation %>%
    dplyr::select(
      probe_id,
      dplyr::any_of(
        candidate_annotation_columns
      )
    ) %>%
    tibble::as_tibble()
  
} else {
  
  annotation_inspection <- feature_annotation %>%
    tibble::as_tibble()
}

cat("\nDua puluh baris pertama anotasi:\n")

print(
  head(
    annotation_inspection,
    20
  ),
  n = 20,
  width = Inf
)

View(
  annotation_inspection[
    seq_len(
      min(
        200,
        nrow(annotation_inspection)
      )
    ),
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 16. Mencari anotasi MMP13
# ------------------------------------------------------------

annotation_as_character <- feature_annotation %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::everything(),
      as.character
    )
  )

annotation_text <- do.call(
  paste,
  c(
    annotation_as_character,
    sep = " | "
  )
)

mmp13_index <- stringr::str_detect(
  stringr::str_to_upper(
    annotation_text
  ),
  "(^|[^A-Z0-9])MMP13([^A-Z0-9]|$)"
)

mmp13_annotation_rows <- feature_annotation[
  mmp13_index,
  ,
  drop = FALSE
] %>%
  tibble::as_tibble()

cat("\nJumlah probe yang mengandung MMP13:\n")
print(nrow(mmp13_annotation_rows))

cat("\nAnotasi probe MMP13:\n")

print(
  mmp13_annotation_rows,
  n = Inf,
  width = Inf
)

if (nrow(mmp13_annotation_rows) == 0) {
  warning(
    "Tidak ditemukan probe dengan anotasi MMP13."
  )
}

# ------------------------------------------------------------
# 17. Membuat tabel nilai ekspresi probe MMP13
# ------------------------------------------------------------

if (nrow(mmp13_annotation_rows) > 0) {
  
  mmp13_probe_ids <-
    mmp13_annotation_rows$probe_id
  
  mmp13_probe_ids_present <- intersect(
    mmp13_probe_ids,
    rownames(expression_included)
  )
  
  mmp13_probe_expression <- expression_included[
    mmp13_probe_ids_present,
    ,
    drop = FALSE
  ]
  
  mmp13_expression_long <- as.data.frame(
    t(
      mmp13_probe_expression
    ),
    check.names = FALSE
  ) %>%
    tibble::rownames_to_column(
      "expression_column"
    ) %>%
    dplyr::left_join(
      sample_metadata_included %>%
        dplyr::select(
          expression_column,
          geo_accession,
          subject_code,
          group,
          plot_label
        ),
      by = "expression_column"
    )
  
} else {
  
  mmp13_probe_expression <- matrix(
    nrow = 0,
    ncol = ncol(expression_included)
  )
  
  mmp13_expression_long <- data.frame()
}

# ------------------------------------------------------------
# 18. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  expression_full,
  file =
    "data_processed/GSE117999_expression_full_24_samples.rds"
)

saveRDS(
  expression_included,
  file =
    "data_processed/GSE117999_expression_included_20_samples.rds"
)

saveRDS(
  sample_metadata,
  file =
    "data_processed/GSE117999_metadata_full_24_samples.rds"
)

saveRDS(
  sample_metadata_included,
  file =
    "data_processed/GSE117999_metadata_included_20_samples.rds"
)

saveRDS(
  feature_annotation,
  file =
    "data_processed/GSE117999_feature_annotation_original.rds"
)

saveRDS(
  mmp13_annotation_rows,
  file =
    "data_processed/GSE117999_MMP13_probe_annotation.rds"
)

saveRDS(
  mmp13_probe_expression,
  file =
    "data_processed/GSE117999_MMP13_probe_expression.rds"
)

# ------------------------------------------------------------
# 19. Membuat ringkasan persiapan
# ------------------------------------------------------------

preparation_summary <- data.frame(
  metric = c(
    "Platform",
    "Number of probes",
    "Initial sample number",
    "Included sample number",
    "Excluded sample number",
    "Control sample number",
    "OA sample number",
    "Expression minimum",
    "Expression maximum",
    "Likely log2 scaled",
    "MMP13-related probe number"
  ),
  
  value = c(
    Biobase::annotation(
      gse117999_eset
    ),
    
    nrow(expression_included),
    
    ncol(expression_full),
    
    ncol(expression_included),
    
    nrow(excluded_samples_found),
    
    sum(
      sample_metadata_included$group ==
        "Control"
    ),
    
    sum(
      sample_metadata_included$group ==
        "OA"
    ),
    
    expression_range[1],
    
    expression_range[2],
    
    likely_log2_scaled,
    
    nrow(mmp13_annotation_rows)
  )
)

# ------------------------------------------------------------
# 20. Menyimpan laporan Excel
# ------------------------------------------------------------

excel_sheets <- list(
  
  Preparation_Summary =
    preparation_summary,
  
  Metadata_All_24 =
    as.data.frame(
      sample_metadata
    ),
  
  Metadata_Included_20 =
    as.data.frame(
      sample_metadata_included
    ),
  
  Excluded_Samples =
    as.data.frame(
      excluded_samples_found
    ),
  
  Expression_Summary =
    as.data.frame(
      sample_expression_summary
    ),
  
  MMP13_Annotation =
    as.data.frame(
      mmp13_annotation_rows
    ),
  
  Annotation_First_500 =
    as.data.frame(
      head(
        feature_annotation,
        500
      )
    )
)

if (nrow(mmp13_expression_long) > 0) {
  
  excel_sheets$MMP13_Probe_Expression <-
    as.data.frame(
      mmp13_expression_long
    )
}

openxlsx::write.xlsx(
  excel_sheets,
  file =
    "metadata/GSE117999_preparation_summary.xlsx",
  overwrite = TRUE
)

# ------------------------------------------------------------
# 21. Menyimpan session information
# ------------------------------------------------------------

sink(
  "metadata/GSE117999_preparation_sessionInfo.txt"
)

print(sessionInfo())

sink()

# ------------------------------------------------------------
# 22. Pesan akhir
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE117999 PREPARATION SELESAI")
message("Platform                 : ",
        Biobase::annotation(gse117999_eset))
message("Jumlah sampel awal       : ",
        ncol(expression_full))
message("Jumlah sampel dianalisis : ",
        ncol(expression_included))
message("Control                  : ",
        sum(sample_metadata_included$group == "Control"))
message("OA                       : ",
        sum(sample_metadata_included$group == "OA"))
message("Jumlah probe             : ",
        nrow(expression_included))
message("Rentang ekspresi         : ",
        round(expression_range[1], 3),
        " sampai ",
        round(expression_range[2], 3))
message("Kemungkinan sudah log2   : ",
        likely_log2_scaled)
message("Probe terkait MMP13      : ",
        nrow(mmp13_annotation_rows))
message("============================================")