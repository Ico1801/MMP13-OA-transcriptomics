# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 08_prepare_GSE57218.R
# PURPOSE : Menyiapkan expression matrix, metadata berpasangan,
#           healthy samples, dan anotasi probe GSE57218
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

message("Semua paket GSE57218 berhasil dimuat.")

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
# 3. Memuat objek GEO
# ------------------------------------------------------------

bulk_geo_file <-
  "data_processed/bulk_geo_original_objects.rds"

if (!file.exists(bulk_geo_file)) {
  stop(
    "File tidak ditemukan: ",
    bulk_geo_file
  )
}

bulk_geo <- readRDS(
  bulk_geo_file
)

if (!"GSE57218" %in% names(bulk_geo)) {
  stop(
    "GSE57218 tidak ditemukan dalam ",
    "bulk_geo_original_objects.rds."
  )
}

gse57218_list <- bulk_geo[["GSE57218"]]

cat("\nJumlah platform GSE57218:\n")
print(length(gse57218_list))

if (length(gse57218_list) != 1) {
  stop(
    "GSE57218 diharapkan hanya memiliki satu platform."
  )
}

gse57218_eset <- gse57218_list[[1]]

cat("\nClass objek:\n")
print(class(gse57218_eset))

cat("\nPlatform annotation:\n")
print(Biobase::annotation(gse57218_eset))

if (!inherits(gse57218_eset, "ExpressionSet")) {
  stop(
    "Objek GSE57218 bukan ExpressionSet."
  )
}

if (
  Biobase::annotation(gse57218_eset) !=
  "GPL6947"
) {
  warning(
    "Platform tidak sama dengan GPL6947."
  )
}

# ------------------------------------------------------------
# 4. Mengambil expression matrix dan metadata
# ------------------------------------------------------------

expression_full <- Biobase::exprs(
  gse57218_eset
)

metadata_full <- Biobase::pData(
  gse57218_eset
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "expression_column"
  )

feature_annotation_original <- Biobase::fData(
  gse57218_eset
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
print(dim(feature_annotation_original))

stopifnot(
  ncol(expression_full) ==
    nrow(metadata_full)
)

stopifnot(
  identical(
    colnames(expression_full),
    metadata_full$expression_column
  )
)

message(
  "Urutan expression matrix dan metadata sesuai."
)

# ------------------------------------------------------------
# 5. Memeriksa metadata
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
      "sex",
      "gender",
      "joint",
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
      unique(
        important_metadata_columns
      )
    )
  ) %>%
  tibble::as_tibble()

cat("\nMetadata penting GSE57218:\n")

print(
  metadata_inspection,
  n = Inf,
  width = Inf
)

View(metadata_inspection)

# ------------------------------------------------------------
# 6. Mengekstrak nomor RAAK dan kondisi cartilage
# ------------------------------------------------------------

title_original <- as.character(
  metadata_full$title
)

title_upper <- stringr::str_to_upper(
  title_original
)

ra_match <- stringr::str_match(
  title_upper,
  "RAAK[_-]?(\\d+)"
)

pair_number <- suppressWarnings(
  as.integer(
    ra_match[, 2]
  )
)

pair_id <- dplyr::if_else(
  !is.na(pair_number),
  paste0(
    "RAAK_",
    pair_number
  ),
  NA_character_
)

condition <- dplyr::case_when(
  
  stringr::str_detect(
    title_upper,
    "_PRESERVED$"
  ) ~ "Preserved",
  
  stringr::str_detect(
    title_upper,
    "_OA$"
  ) ~ "OA",
  
  stringr::str_detect(
    title_upper,
    "_HEALTHY$"
  ) ~ "Healthy",
  
  TRUE ~ NA_character_
)

sample_metadata <- metadata_full %>%
  dplyr::mutate(
    title_original =
      title_original,
    
    pair_number =
      pair_number,
    
    pair_id =
      pair_id,
    
    condition =
      condition,
    
    analysis_category =
      dplyr::case_when(
        condition %in%
          c(
            "Preserved",
            "OA"
          ) ~
          "Primary paired analysis",
        
        condition ==
          "Healthy" ~
          "Exploratory healthy reference",
        
        TRUE ~
          "Unclassified"
      )
  )

cat("\nMetadata setelah ekstraksi kondisi:\n")

print(
  sample_metadata %>%
    dplyr::select(
      expression_column,
      geo_accession,
      title_original,
      pair_number,
      pair_id,
      condition,
      analysis_category
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 7. Validasi klasifikasi sampel
# ------------------------------------------------------------

cat("\nJumlah pair ID yang kosong:\n")
print(
  sum(
    is.na(sample_metadata$pair_id)
  )
)

cat("\nJumlah condition yang kosong:\n")
print(
  sum(
    is.na(sample_metadata$condition)
  )
)

cat("\nDistribusi kondisi seluruh sampel:\n")
print(
  table(
    sample_metadata$condition,
    useNA = "ifany"
  )
)

if (anyNA(sample_metadata$pair_id)) {
  
  print(
    sample_metadata %>%
      dplyr::filter(
        is.na(pair_id)
      ) %>%
      dplyr::select(
        expression_column,
        geo_accession,
        title_original
      ) %>%
      tibble::as_tibble(),
    n = Inf,
    width = Inf
  )
  
  stop(
    "Ada RAAK pair ID yang tidak dapat diekstrak."
  )
}

if (anyNA(sample_metadata$condition)) {
  
  print(
    sample_metadata %>%
      dplyr::filter(
        is.na(condition)
      ) %>%
      dplyr::select(
        expression_column,
        geo_accession,
        title_original
      ) %>%
      tibble::as_tibble(),
    n = Inf,
    width = Inf
  )
  
  stop(
    "Ada kondisi sampel yang belum diklasifikasikan."
  )
}

expected_condition_counts <- c(
  Healthy = 7,
  OA = 33,
  Preserved = 33
)

observed_condition_counts <- table(
  sample_metadata$condition
)

for (
  current_condition in
  names(expected_condition_counts)
) {
  
  observed_number <- if (
    current_condition %in%
    names(observed_condition_counts)
  ) {
    as.integer(
      observed_condition_counts[
        current_condition
      ]
    )
  } else {
    0L
  }
  
  expected_number <-
    expected_condition_counts[
      current_condition
    ]
  
  if (observed_number != expected_number) {
    stop(
      "Jumlah sampel ",
      current_condition,
      " tidak sesuai. Ditemukan ",
      observed_number,
      ", diharapkan ",
      expected_number,
      "."
    )
  }
}

# ------------------------------------------------------------
# 8. Memeriksa kelengkapan pasangan
# ------------------------------------------------------------

pair_summary <- sample_metadata %>%
  dplyr::filter(
    condition %in%
      c(
        "Preserved",
        "OA"
      )
  ) %>%
  dplyr::group_by(
    pair_id,
    pair_number
  ) %>%
  dplyr::summarise(
    total_samples =
      dplyr::n(),
    
    preserved_samples =
      sum(
        condition ==
          "Preserved"
      ),
    
    oa_samples =
      sum(
        condition ==
          "OA"
      ),
    
    complete_pair =
      total_samples == 2 &
      preserved_samples == 1 &
      oa_samples == 1,
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    pair_number
  )

cat("\nRingkasan pasangan preserved–OA:\n")

print(
  pair_summary,
  n = Inf,
  width = Inf
)

cat("\nJumlah pasangan lengkap:\n")
print(
  sum(
    pair_summary$complete_pair
  )
)

cat("\nJumlah pasangan tidak lengkap:\n")
print(
  sum(
    !pair_summary$complete_pair
  )
)

if (
  any(
    !pair_summary$complete_pair
  )
) {
  
  print(
    pair_summary %>%
      dplyr::filter(
        !complete_pair
      ),
    n = Inf,
    width = Inf
  )
  
  stop(
    "Ada pasangan preserved–OA yang tidak lengkap."
  )
}

if (
  nrow(pair_summary) != 33
) {
  stop(
    "Jumlah pasangan lengkap bukan 33."
  )
}

# ------------------------------------------------------------
# 9. Metadata primary paired analysis
# ------------------------------------------------------------

complete_pair_ids <- pair_summary$pair_id[
  pair_summary$complete_pair
]

sample_metadata_primary <- sample_metadata %>%
  dplyr::filter(
    pair_id %in%
      complete_pair_ids,
    
    condition %in%
      c(
        "Preserved",
        "OA"
      )
  ) %>%
  dplyr::mutate(
    condition = factor(
      condition,
      levels = c(
        "Preserved",
        "OA"
      )
    )
  ) %>%
  dplyr::arrange(
    pair_number,
    condition
  ) %>%
  dplyr::mutate(
    pair_id = factor(
      pair_id,
      levels = unique(
        pair_id
      )
    ),
    
    plot_label =
      paste0(
        "P",
        sprintf(
          "%02d",
          match(
            pair_id,
            levels(pair_id)
          )
        ),
        "_",
        dplyr::if_else(
          condition ==
            "Preserved",
          "Pres",
          "OA"
        )
      )
  )

cat("\nMetadata primary paired analysis:\n")

print(
  sample_metadata_primary %>%
    dplyr::select(
      expression_column,
      geo_accession,
      title_original,
      pair_id,
      pair_number,
      condition,
      plot_label
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

cat("\nDistribusi primary condition:\n")

print(
  table(
    sample_metadata_primary$condition
  )
)

cat("\nJumlah primary samples:\n")
print(
  nrow(
    sample_metadata_primary
  )
)

cat("\nJumlah primary pairs:\n")
print(
  dplyr::n_distinct(
    sample_metadata_primary$pair_id
  )
)

stopifnot(
  nrow(sample_metadata_primary) ==
    66
)

stopifnot(
  dplyr::n_distinct(
    sample_metadata_primary$pair_id
  ) ==
    33
)

# ------------------------------------------------------------
# 10. Metadata healthy reference
# ------------------------------------------------------------

sample_metadata_healthy <- sample_metadata %>%
  dplyr::filter(
    condition ==
      "Healthy"
  ) %>%
  dplyr::arrange(
    pair_number
  ) %>%
  dplyr::mutate(
    plot_label =
      paste0(
        "H",
        dplyr::row_number()
      )
  )

cat("\nMetadata healthy reference:\n")

print(
  sample_metadata_healthy %>%
    dplyr::select(
      expression_column,
      geo_accession,
      title_original,
      pair_id,
      pair_number,
      condition,
      plot_label
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

stopifnot(
  nrow(sample_metadata_healthy) ==
    7
)

# ------------------------------------------------------------
# 11. Expression matrix primary paired
# ------------------------------------------------------------

primary_sample_columns <-
  sample_metadata_primary$expression_column

stopifnot(
  all(
    primary_sample_columns %in%
      colnames(expression_full)
  )
)

expression_primary <- expression_full[
  ,
  primary_sample_columns,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(expression_primary),
    sample_metadata_primary$expression_column
  )
)

cat("\nDimensi expression matrix primary:\n")
print(
  dim(
    expression_primary
  )
)

stopifnot(
  ncol(expression_primary) ==
    66
)

# ------------------------------------------------------------
# 12. Expression matrix healthy reference
# ------------------------------------------------------------

healthy_sample_columns <-
  sample_metadata_healthy$expression_column

expression_healthy <- expression_full[
  ,
  healthy_sample_columns,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(expression_healthy),
    sample_metadata_healthy$expression_column
  )
)

cat("\nDimensi expression matrix healthy:\n")
print(
  dim(
    expression_healthy
  )
)

stopifnot(
  ncol(expression_healthy) ==
    7
)

# ------------------------------------------------------------
# 13. Pemeriksaan nilai expression matrix
# ------------------------------------------------------------

cat("\nAda NA pada expression matrix lengkap:\n")
print(
  anyNA(
    expression_full
  )
)

cat("\nAda NA pada primary paired matrix:\n")
print(
  anyNA(
    expression_primary
  )
)

cat("\nAda nilai non-finite pada primary matrix:\n")
print(
  any(
    !is.finite(
      expression_primary
    )
  )
)

expression_range <- range(
  expression_primary,
  na.rm = TRUE
)

expression_quantiles <- quantile(
  expression_primary,
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

cat("\nRentang expression primary:\n")
print(
  expression_range
)

cat("\nQuantile expression primary:\n")
print(
  expression_quantiles
)

likely_processed_log_scale <- isTRUE(
  unname(
    expression_range[2] < 50
  )
)

cat(
  "\nData kemungkinan sudah processed/VST scale:\n"
)

print(
  likely_processed_log_scale
)

if (
  anyNA(expression_primary) ||
  any(
    !is.finite(
      expression_primary
    )
  )
) {
  stop(
    "Primary expression matrix memiliki NA ",
    "atau nilai non-finite."
  )
}

if (!likely_processed_log_scale) {
  warning(
    "Nilai expression mungkin belum berada ",
    "pada processed log-like scale."
  )
}

# ------------------------------------------------------------
# 14. Ringkasan distribusi setiap sampel
# ------------------------------------------------------------

sample_expression_summary <- tibble::tibble(
  expression_column =
    colnames(
      expression_primary
    ),
  
  minimum =
    apply(
      expression_primary,
      2,
      min,
      na.rm = TRUE
    ),
  
  q1 =
    apply(
      expression_primary,
      2,
      stats::quantile,
      probs = 0.25,
      na.rm = TRUE
    ),
  
  median =
    apply(
      expression_primary,
      2,
      stats::median,
      na.rm = TRUE
    ),
  
  mean =
    colMeans(
      expression_primary,
      na.rm = TRUE
    ),
  
  q3 =
    apply(
      expression_primary,
      2,
      stats::quantile,
      probs = 0.75,
      na.rm = TRUE
    ),
  
  maximum =
    apply(
      expression_primary,
      2,
      max,
      na.rm = TRUE
    )
) %>%
  dplyr::left_join(
    sample_metadata_primary %>%
      dplyr::select(
        expression_column,
        geo_accession,
        pair_id,
        pair_number,
        condition,
        plot_label
      ),
    by = "expression_column"
  )

cat("\nRingkasan ekspresi primary samples:\n")

print(
  sample_expression_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 15. Memeriksa feature annotation
# ------------------------------------------------------------

cat("\nNama seluruh kolom feature annotation:\n")

print(
  colnames(
    feature_annotation_original
  )
)

candidate_annotation_columns <- grep(
  pattern = paste(
    c(
      "symbol",
      "gene",
      "entrez",
      "refseq",
      "accession",
      "description",
      "definition",
      "transcript"
    ),
    collapse = "|"
  ),
  x = colnames(
    feature_annotation_original
  ),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nKolom kandidat anotasi:\n")

print(
  candidate_annotation_columns
)

if (
  length(
    candidate_annotation_columns
  ) > 0
) {
  
  annotation_inspection <-
    feature_annotation_original %>%
    dplyr::select(
      probe_id,
      dplyr::any_of(
        candidate_annotation_columns
      )
    ) %>%
    tibble::as_tibble()
  
} else {
  
  annotation_inspection <-
    feature_annotation_original %>%
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
        nrow(
          annotation_inspection
        )
      )
    ),
    ,
    drop = FALSE
  ]
)

# ------------------------------------------------------------
# 16. Fungsi mencari anotasi MMP13
# ------------------------------------------------------------

find_mmp13_rows <- function(
    annotation_table
) {
  
  annotation_character <-
    annotation_table %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::everything(),
        as.character
      )
    )
  
  annotation_text <- do.call(
    paste,
    c(
      annotation_character,
      sep = " | "
    )
  )
  
  mmp13_index <- stringr::str_detect(
    stringr::str_to_upper(
      annotation_text
    ),
    "(^|[^A-Z0-9])MMP13([^A-Z0-9]|$)"
  )
  
  annotation_table[
    mmp13_index,
    ,
    drop = FALSE
  ]
}

mmp13_annotation_rows <- find_mmp13_rows(
  feature_annotation_original
) %>%
  tibble::as_tibble()

annotation_source <- "ExpressionSet feature data"

# ------------------------------------------------------------
# 17. Fallback ke GPL6947 jika anotasi MMP13 belum ditemukan
# ------------------------------------------------------------

if (
  nrow(
    mmp13_annotation_rows
  ) == 0
) {
  
  message(
    "MMP13 belum ditemukan di fData. ",
    "Mencoba mengambil tabel GPL6947 dari GEO."
  )
  
  gpl6947_object <- GEOquery::getGEO(
    "GPL6947"
  )
  
  gpl6947_table <- GEOquery::Table(
    gpl6947_object
  ) %>%
    as.data.frame()
  
  if (
    !"ID" %in%
    colnames(
      gpl6947_table
    )
  ) {
    stop(
      "Kolom ID tidak ditemukan pada GPL6947."
    )
  }
  
  gpl6947_table <- gpl6947_table %>%
    dplyr::mutate(
      probe_id =
        as.character(
          ID
        )
    )
  
  mmp13_annotation_rows <-
    find_mmp13_rows(
      gpl6947_table
    ) %>%
    dplyr::filter(
      probe_id %in%
        rownames(
          expression_primary
        )
    ) %>%
    tibble::as_tibble()
  
  annotation_source <-
    "GEO GPL6947 platform table"
  
  saveRDS(
    gpl6947_table,
    file =
      "data_processed/GPL6947_platform_annotation.rds"
  )
}

cat("\nSumber anotasi MMP13:\n")
print(annotation_source)

cat("\nJumlah probe terkait MMP13:\n")
print(
  nrow(
    mmp13_annotation_rows
  )
)

cat("\nAnotasi probe MMP13:\n")

print(
  mmp13_annotation_rows,
  n = Inf,
  width = Inf
)

if (
  nrow(
    mmp13_annotation_rows
  ) == 0
) {
  stop(
    "Tidak ditemukan probe MMP13 pada GSE57218."
  )
}

# ------------------------------------------------------------
# 18. Memastikan kolom probe ID untuk MMP13
# ------------------------------------------------------------

if (
  !"probe_id" %in%
  colnames(
    mmp13_annotation_rows
  )
) {
  stop(
    "Kolom probe_id tidak tersedia ",
    "pada anotasi MMP13."
  )
}

mmp13_probe_ids <- unique(
  as.character(
    mmp13_annotation_rows$probe_id
  )
)

mmp13_probe_ids_present <- intersect(
  mmp13_probe_ids,
  rownames(
    expression_primary
  )
)

cat("\nProbe MMP13 yang tersedia dalam matriks:\n")

print(
  mmp13_probe_ids_present
)

if (
  length(
    mmp13_probe_ids_present
  ) == 0
) {
  stop(
    "Probe MMP13 tidak ditemukan ",
    "pada expression matrix."
  )
}

# ------------------------------------------------------------
# 19. Expression MMP13 primary paired samples
# ------------------------------------------------------------

mmp13_expression_matrix <-
  expression_primary[
    mmp13_probe_ids_present,
    ,
    drop = FALSE
  ]

mmp13_expression_long <- as.data.frame(
  t(
    mmp13_expression_matrix
  ),
  check.names = FALSE
) %>%
  tibble::rownames_to_column(
    "expression_column"
  ) %>%
  dplyr::left_join(
    sample_metadata_primary %>%
      dplyr::select(
        expression_column,
        geo_accession,
        pair_id,
        pair_number,
        condition,
        plot_label
      ),
    by = "expression_column"
  ) %>%
  tibble::as_tibble()

cat("\nEkspresi probe MMP13 primary samples:\n")

print(
  mmp13_expression_long,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 20. Membuat tabel pasangan eksplisit
# ------------------------------------------------------------

pair_mapping <- sample_metadata_primary %>%
  dplyr::select(
    pair_id,
    pair_number,
    condition,
    expression_column,
    geo_accession,
    title_original,
    plot_label
  ) %>%
  dplyr::arrange(
    pair_number,
    condition
  ) %>%
  tibble::as_tibble()

cat("\nMapping pasangan OA dan preserved:\n")

print(
  pair_mapping,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 21. Ringkasan persiapan
# ------------------------------------------------------------

preparation_summary <- data.frame(
  metric = c(
    "Dataset",
    "Platform",
    "Total number of probes",
    "Initial sample number",
    "Preserved samples",
    "OA-affected samples",
    "Healthy reference samples",
    "Complete preserved-OA pairs",
    "Primary paired sample number",
    "Expression minimum",
    "Expression maximum",
    "Likely processed VST scale",
    "MMP13-related probe number",
    "MMP13 annotation source"
  ),
  
  value = c(
    "GSE57218",
    
    Biobase::annotation(
      gse57218_eset
    ),
    
    nrow(
      expression_full
    ),
    
    ncol(
      expression_full
    ),
    
    sum(
      sample_metadata$condition ==
        "Preserved"
    ),
    
    sum(
      sample_metadata$condition ==
        "OA"
    ),
    
    sum(
      sample_metadata$condition ==
        "Healthy"
    ),
    
    sum(
      pair_summary$complete_pair
    ),
    
    ncol(
      expression_primary
    ),
    
    expression_range[1],
    
    expression_range[2],
    
    likely_processed_log_scale,
    
    length(
      mmp13_probe_ids_present
    ),
    
    annotation_source
  )
)

print(
  preparation_summary
)

# ------------------------------------------------------------
# 22. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  expression_full,
  file =
    "data_processed/GSE57218_expression_full_73_samples.rds"
)

saveRDS(
  expression_primary,
  file =
    "data_processed/GSE57218_expression_primary_33_pairs.rds"
)

saveRDS(
  expression_healthy,
  file =
    "data_processed/GSE57218_expression_healthy_7_samples.rds"
)

saveRDS(
  sample_metadata,
  file =
    "data_processed/GSE57218_metadata_full_73_samples.rds"
)

saveRDS(
  sample_metadata_primary,
  file =
    "data_processed/GSE57218_metadata_primary_33_pairs.rds"
)

saveRDS(
  sample_metadata_healthy,
  file =
    "data_processed/GSE57218_metadata_healthy_7_samples.rds"
)

saveRDS(
  pair_summary,
  file =
    "data_processed/GSE57218_pair_summary.rds"
)

saveRDS(
  pair_mapping,
  file =
    "data_processed/GSE57218_pair_mapping.rds"
)

saveRDS(
  feature_annotation_original,
  file =
    "data_processed/GSE57218_feature_annotation_original.rds"
)

saveRDS(
  mmp13_annotation_rows,
  file =
    "data_processed/GSE57218_MMP13_probe_annotation.rds"
)

saveRDS(
  mmp13_expression_matrix,
  file =
    "data_processed/GSE57218_MMP13_primary_expression.rds"
)

# ------------------------------------------------------------
# 23. Menyimpan laporan Excel
# ------------------------------------------------------------

excel_sheets <- list(
  Preparation_Summary =
    preparation_summary,
  
  Metadata_All_73 =
    as.data.frame(
      sample_metadata
    ),
  
  Metadata_Primary_66 =
    as.data.frame(
      sample_metadata_primary
    ),
  
  Metadata_Healthy_7 =
    as.data.frame(
      sample_metadata_healthy
    ),
  
  Pair_Summary =
    as.data.frame(
      pair_summary
    ),
  
  Pair_Mapping =
    as.data.frame(
      pair_mapping
    ),
  
  Expression_Summary =
    as.data.frame(
      sample_expression_summary
    ),
  
  MMP13_Annotation =
    as.data.frame(
      mmp13_annotation_rows
    ),
  
  MMP13_Primary_Expression =
    as.data.frame(
      mmp13_expression_long
    ),
  
  Annotation_First_500 =
    as.data.frame(
      head(
        feature_annotation_original,
        500
      )
    )
)

openxlsx::write.xlsx(
  excel_sheets,
  file =
    "metadata/GSE57218_preparation_summary.xlsx",
  overwrite = TRUE
)

# ------------------------------------------------------------
# 24. Menyimpan session information
# ------------------------------------------------------------

sink(
  "metadata/GSE57218_preparation_sessionInfo.txt"
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 25. Pesan akhir
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE57218 PREPARATION SELESAI")
message("Platform                   : ",
        Biobase::annotation(gse57218_eset))
message("Jumlah sampel awal         : ",
        ncol(expression_full))
message("Preserved                  : ",
        sum(sample_metadata$condition == "Preserved"))
message("OA affected                : ",
        sum(sample_metadata$condition == "OA"))
message("Healthy reference          : ",
        sum(sample_metadata$condition == "Healthy"))
message("Pasangan lengkap           : ",
        sum(pair_summary$complete_pair))
message("Sampel analisis utama      : ",
        ncol(expression_primary))
message("Jumlah probe               : ",
        nrow(expression_primary))
message("Rentang expression         : ",
        round(expression_range[1], 3),
        " sampai ",
        round(expression_range[2], 3))
message("Processed/VST-like scale   : ",
        likely_processed_log_scale)
message("Probe MMP13 tersedia       : ",
        length(mmp13_probe_ids_present))
message("============================================")