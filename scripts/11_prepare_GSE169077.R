# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 11_prepare_GSE169077.R
# PURPOSE : Menyiapkan metadata dan memproses raw CEL
#           GSE169077 menggunakan RMA
#
# DATASET:
#   GPL96 Affymetrix Human Genome U133A Array
#   5 normal pooled samples
#   6 OA pooled samples
#
# IMPORTANT:
#   Setiap array merupakan pooled RNA sample.
#   Unit statistik tetap pool, bukan individu penyusun pool.
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Memuat paket
# ------------------------------------------------------------

required_packages <- c(
  "GEOquery",
  "Biobase",
  "affy",
  "hgu133acdf",
  "limma",
  "dplyr",
  "tibble",
  "stringr",
  "openxlsx",
  "R.utils"
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

library(GEOquery)
library(Biobase)
library(affy)
library(hgu133acdf)
library(limma)
library(dplyr)
library(tibble)
library(stringr)
library(openxlsx)

message(
  "Semua paket preparation GSE169077 berhasil dimuat."
)

# ------------------------------------------------------------
# 2. Memastikan folder tersedia
# ------------------------------------------------------------

required_folders <- c(
  "data_raw",
  "data_raw/GSE169077",
  "data_raw/GSE169077/extracted",
  "data_raw/GSE169077/CEL",
  "data_processed",
  "metadata",
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
# 3. Memuat objek GEO yang sebelumnya diunduh
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

if (!"GSE169077" %in% names(bulk_geo)) {
  stop(
    "GSE169077 tidak ditemukan dalam ",
    "bulk_geo_original_objects.rds."
  )
}

gse169077_list <-
  bulk_geo[["GSE169077"]]

cat("\nJumlah platform GSE169077:\n")
print(length(gse169077_list))

if (length(gse169077_list) != 1) {
  stop(
    "GSE169077 diharapkan memiliki satu platform."
  )
}

gse169077_eset <-
  gse169077_list[[1]]

cat("\nClass objek GEO:\n")
print(class(gse169077_eset))

cat("\nPlatform GEO:\n")
print(
  Biobase::annotation(
    gse169077_eset
  )
)

stopifnot(
  inherits(
    gse169077_eset,
    "ExpressionSet"
  )
)

if (
  Biobase::annotation(
    gse169077_eset
  ) != "GPL96"
) {
  warning(
    "Platform bukan GPL96. ",
    "Periksa kembali objek GEO."
  )
}

# ------------------------------------------------------------
# 4. Mengambil processed matrix dan metadata GEO
# ------------------------------------------------------------

expression_geo_processed <-
  Biobase::exprs(
    gse169077_eset
  )

metadata_full <-
  Biobase::pData(
    gse169077_eset
  ) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "expression_column"
  )

feature_annotation_original <-
  Biobase::fData(
    gse169077_eset
  ) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(
    "probe_id"
  )

cat("\nDimensi processed GEO matrix:\n")
print(dim(expression_geo_processed))

cat("\nDimensi metadata:\n")
print(dim(metadata_full))

cat("\nDimensi feature annotation:\n")
print(dim(feature_annotation_original))

stopifnot(
  ncol(expression_geo_processed) ==
    nrow(metadata_full)
)

stopifnot(
  identical(
    colnames(expression_geo_processed),
    metadata_full$expression_column
  )
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
      "disease",
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

cat("\nMetadata penting GSE169077:\n")

print(
  metadata_inspection,
  n = Inf,
  width = Inf
)

View(metadata_inspection)

# ------------------------------------------------------------
# 6. Menggabungkan teks metadata
# ------------------------------------------------------------

metadata_text_columns <- intersect(
  c(
    "title",
    "source_name_ch1",
    "characteristics_ch1",
    "characteristics_ch1.1",
    "characteristics_ch1.2",
    "description"
  ),
  colnames(metadata_full)
)

if (length(metadata_text_columns) == 0) {
  stop(
    "Tidak ada kolom metadata yang dapat digunakan ",
    "untuk klasifikasi kelompok."
  )
}

metadata_combined_text <- apply(
  metadata_full[
    ,
    metadata_text_columns,
    drop = FALSE
  ],
  MARGIN = 1,
  FUN = function(current_row) {
    paste(
      as.character(current_row),
      collapse = " | "
    )
  }
)

metadata_lower <- stringr::str_to_lower(
  metadata_combined_text
)

title_upper <- stringr::str_to_upper(
  as.character(
    metadata_full$title
  )
)

# ------------------------------------------------------------
# 7. Membuat metadata analisis
# ------------------------------------------------------------

group_classification <- dplyr::case_when(
  
  stringr::str_detect(
    metadata_lower,
    paste(
      c(
        "late stage oa",
        "osteoarthritis",
        "disease state: oa"
      ),
      collapse = "|"
    )
  ) ~ "OA",
  
  stringr::str_detect(
    metadata_lower,
    paste(
      c(
        "normal cartilage",
        "disease state: normal",
        "non-osteoarthritis"
      ),
      collapse = "|"
    )
  ) ~ "Normal",
  
  stringr::str_detect(
    title_upper,
    "^OA"
  ) ~ "OA",
  
  stringr::str_detect(
    title_upper,
    "^N[0-9]+"
  ) ~ "Normal",
  
  TRUE ~ NA_character_
)

sample_metadata <- metadata_full %>%
  dplyr::mutate(
    
    title_original =
      as.character(title),
    
    combined_metadata_text =
      metadata_combined_text,
    
    group =
      group_classification,
    
    pooled_sample =
      TRUE,
    
    individuals_per_pool =
      5L,
    
    statistical_unit =
      "RNA pool",
    
    group = factor(
      group,
      levels = c(
        "Normal",
        "OA"
      )
    ),
    
    pool_label =
      title_original
  )

cat("\nMetadata analysis:\n")

print(
  sample_metadata %>%
    dplyr::select(
      expression_column,
      geo_accession,
      title_original,
      group,
      pooled_sample,
      individuals_per_pool,
      statistical_unit
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 8. Validasi kelompok
# ------------------------------------------------------------

cat("\nJumlah group yang kosong:\n")
print(
  sum(
    is.na(sample_metadata$group)
  )
)

cat("\nDistribusi group:\n")
print(
  table(
    sample_metadata$group,
    useNA = "ifany"
  )
)

if (anyNA(sample_metadata$group)) {
  
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
    "Ada sampel yang belum dapat diklasifikasikan."
  )
}

if (
  sum(
    sample_metadata$group ==
    "Normal"
  ) != 5
) {
  stop(
    "Jumlah normal pool bukan 5."
  )
}

if (
  sum(
    sample_metadata$group ==
    "OA"
  ) != 6
) {
  stop(
    "Jumlah OA pool bukan 6."
  )
}

if (
  anyDuplicated(
    sample_metadata$geo_accession
  )
) {
  stop(
    "Terdapat duplikasi GEO accession."
  )
}

# ------------------------------------------------------------
# 9. Ringkasan processed GEO matrix
# ------------------------------------------------------------

cat("\nRentang processed GEO matrix:\n")
print(
  range(
    expression_geo_processed,
    na.rm = TRUE
  )
)

cat("\nQuantile processed GEO matrix:\n")
print(
  quantile(
    expression_geo_processed,
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
)

# Processed GEO signal disimpan hanya sebagai pembanding.
# Analisis utama menggunakan raw CEL yang diproses ulang dengan RMA.

# ------------------------------------------------------------
# 10. Menemukan atau mengunduh raw CEL archive
# ------------------------------------------------------------

raw_data_directory <-
  "data_raw/GSE169077"

raw_tar_candidates <- list.files(
  path =
    raw_data_directory,
  pattern =
    "^GSE169077_RAW\\.tar$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(raw_tar_candidates) == 0) {
  
  message(
    "GSE169077_RAW.tar belum tersedia. ",
    "Mengunduh supplementary files dari GEO..."
  )
  
  GEOquery::getGEOSuppFiles(
    GEO = "GSE169077",
    makeDirectory = FALSE,
    baseDir = raw_data_directory
  )
  
  raw_tar_candidates <- list.files(
    path =
      raw_data_directory,
    pattern =
      "^GSE169077_RAW\\.tar$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

if (length(raw_tar_candidates) == 0) {
  stop(
    "GSE169077_RAW.tar tidak ditemukan ",
    "setelah proses download."
  )
}

raw_tar_file <-
  raw_tar_candidates[1]

cat("\nRaw CEL archive:\n")
print(raw_tar_file)

# ------------------------------------------------------------
# 11. Mengekstrak raw archive
# ------------------------------------------------------------

extracted_directory <-
  "data_raw/GSE169077/extracted"

cel_gz_files <- list.files(
  path =
    extracted_directory,
  pattern =
    "\\.CEL\\.gz$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(cel_gz_files) != 11) {
  
  message(
    "Mengekstrak GSE169077_RAW.tar..."
  )
  
  utils::untar(
    tarfile =
      raw_tar_file,
    exdir =
      extracted_directory
  )
  
  cel_gz_files <- list.files(
    path =
      extracted_directory,
    pattern =
      "\\.CEL\\.gz$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

cat("\nJumlah compressed CEL files:\n")
print(length(cel_gz_files))

if (length(cel_gz_files) != 11) {
  stop(
    "Jumlah compressed CEL files bukan 11."
  )
}

# ------------------------------------------------------------
# 12. Mengekstrak file .CEL.gz menjadi .CEL
# ------------------------------------------------------------

cel_directory <-
  "data_raw/GSE169077/CEL"

cel_files <- character(
  length(cel_gz_files)
)

for (
  current_index in
  seq_along(cel_gz_files)
) {
  
  current_gz_file <-
    cel_gz_files[current_index]
  
  output_cel_file <- file.path(
    cel_directory,
    stringr::str_remove(
      basename(current_gz_file),
      regex(
        "\\.gz$",
        ignore_case = TRUE
      )
    )
  )
  
  if (!file.exists(output_cel_file)) {
    
    R.utils::gunzip(
      filename =
        current_gz_file,
      destname =
        output_cel_file,
      overwrite =
        FALSE,
      remove =
        FALSE
    )
  }
  
  cel_files[current_index] <-
    output_cel_file
}

cat("\nJumlah uncompressed CEL files:\n")
print(length(cel_files))

stopifnot(
  length(cel_files) == 11
)

stopifnot(
  all(
    file.exists(cel_files)
  )
)

# ------------------------------------------------------------
# 13. Membuat manifest CEL
# ------------------------------------------------------------

cel_manifest <- tibble::tibble(
  
  cel_file =
    normalizePath(
      cel_files,
      winslash = "/",
      mustWork = TRUE
    ),
  
  cel_filename =
    basename(cel_files),
  
  geo_accession =
    stringr::str_extract(
      basename(cel_files),
      "GSM[0-9]+"
    )
)

cat("\nCEL file manifest:\n")

print(
  cel_manifest,
  n = Inf,
  width = Inf
)

if (anyNA(cel_manifest$geo_accession)) {
  stop(
    "Ada CEL file tanpa GSM accession pada nama file."
  )
}

if (
  anyDuplicated(
    cel_manifest$geo_accession
  )
) {
  stop(
    "Terdapat duplikasi GSM accession pada CEL manifest."
  )
}

missing_cel_accessions <- setdiff(
  sample_metadata$geo_accession,
  cel_manifest$geo_accession
)

extra_cel_accessions <- setdiff(
  cel_manifest$geo_accession,
  sample_metadata$geo_accession
)

cat("\nGSM metadata tanpa CEL:\n")
print(missing_cel_accessions)

cat("\nGSM CEL tanpa metadata:\n")
print(extra_cel_accessions)

if (
  length(missing_cel_accessions) > 0 ||
  length(extra_cel_accessions) > 0
) {
  stop(
    "Mapping CEL dan metadata tidak lengkap."
  )
}

# ------------------------------------------------------------
# 14. Mengurutkan CEL sesuai metadata
# ------------------------------------------------------------

cel_match <- match(
  sample_metadata$geo_accession,
  cel_manifest$geo_accession
)

stopifnot(
  !anyNA(cel_match)
)

cel_manifest_ordered <-
  cel_manifest[
    cel_match,
    ,
    drop = FALSE
  ]

stopifnot(
  identical(
    cel_manifest_ordered$geo_accession,
    sample_metadata$geo_accession
  )
)

# ------------------------------------------------------------
# 15. Membaca raw Affymetrix CEL files
# ------------------------------------------------------------

message(
  "Membaca 11 raw CEL files..."
)

affy_raw <- affy::ReadAffy(
  filenames =
    cel_manifest_ordered$cel_file
)

Biobase::sampleNames(
  affy_raw
) <- sample_metadata$geo_accession

cat("\nAffyBatch dimensions:\n")
print(dim(affy_raw))

cat("\nCDF name:\n")
print(
  affy::cdfName(
    affy_raw
  )
)

# ------------------------------------------------------------
# 16. RMA preprocessing
# ------------------------------------------------------------

message(
  "Menjalankan RMA background correction, ",
  "quantile normalization, dan summarization..."
)

rma_eset <- affy::rma(
  object =
    affy_raw,
  background =
    TRUE,
  normalize =
    TRUE,
  verbose =
    TRUE
)

expression_rma <- Biobase::exprs(
  rma_eset
)

colnames(expression_rma) <-
  sample_metadata$geo_accession

Biobase::sampleNames(
  rma_eset
) <- sample_metadata$geo_accession

cat("\nDimensi RMA expression matrix:\n")
print(dim(expression_rma))

cat("\nRentang RMA expression:\n")
print(
  range(
    expression_rma,
    na.rm = TRUE
  )
)

cat("\nQuantile RMA expression:\n")
print(
  quantile(
    expression_rma,
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
)

stopifnot(
  ncol(expression_rma) == 11
)

stopifnot(
  identical(
    colnames(expression_rma),
    sample_metadata$geo_accession
  )
)

stopifnot(
  !anyNA(expression_rma)
)

stopifnot(
  all(
    is.finite(expression_rma)
  )
)

# ------------------------------------------------------------
# 17. Ringkasan RMA setiap pool
# ------------------------------------------------------------

sample_expression_summary <- tibble::tibble(
  
  geo_accession =
    colnames(expression_rma),
  
  minimum =
    apply(
      expression_rma,
      2,
      min
    ),
  
  q1 =
    apply(
      expression_rma,
      2,
      stats::quantile,
      probs = 0.25
    ),
  
  median =
    apply(
      expression_rma,
      2,
      stats::median
    ),
  
  mean =
    colMeans(
      expression_rma
    ),
  
  q3 =
    apply(
      expression_rma,
      2,
      stats::quantile,
      probs = 0.75
    ),
  
  maximum =
    apply(
      expression_rma,
      2,
      max
    )
) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(
        geo_accession,
        title_original,
        group,
        pool_label
      ),
    by =
      "geo_accession"
  )

cat("\nRMA expression summary per pool:\n")

print(
  sample_expression_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 18. Memeriksa feature annotation
# ------------------------------------------------------------

cat("\nNama feature annotation columns:\n")
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
      "accession",
      "description",
      "title"
    ),
    collapse = "|"
  ),
  x =
    colnames(
      feature_annotation_original
    ),
  value = TRUE,
  ignore.case = TRUE
)

cat("\nKolom kandidat annotation:\n")
print(candidate_annotation_columns)

annotation_inspection <-
  feature_annotation_original %>%
  dplyr::select(
    probe_id,
    dplyr::any_of(
      candidate_annotation_columns
    )
  ) %>%
  tibble::as_tibble()

cat("\nDua puluh annotation rows pertama:\n")

print(
  head(
    annotation_inspection,
    20
  ),
  n = 20,
  width = Inf
)

# ------------------------------------------------------------
# 19. Fungsi mencari MMP13
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

annotation_source <-
  "ExpressionSet feature data"

# ------------------------------------------------------------
# 20. Fallback ke GPL96 platform annotation
# ------------------------------------------------------------

if (nrow(mmp13_annotation_rows) == 0) {
  
  message(
    "MMP13 tidak ditemukan pada fData. ",
    "Mengambil GPL96 platform table..."
  )
  
  gpl96_object <- GEOquery::getGEO(
    "GPL96"
  )
  
  gpl96_table <- GEOquery::Table(
    gpl96_object
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(
      probe_id =
        as.character(ID)
    )
  
  mmp13_annotation_rows <-
    find_mmp13_rows(
      gpl96_table
    ) %>%
    dplyr::filter(
      probe_id %in%
        rownames(expression_rma)
    ) %>%
    tibble::as_tibble()
  
  annotation_source <-
    "GEO GPL96 platform table"
  
  saveRDS(
    gpl96_table,
    file =
      "data_processed/GPL96_platform_annotation.rds"
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

if (nrow(mmp13_annotation_rows) == 0) {
  stop(
    "Tidak ditemukan probe MMP13."
  )
}

if (
  !"probe_id" %in%
  colnames(
    mmp13_annotation_rows
  )
) {
  stop(
    "Kolom probe_id tidak ditemukan ",
    "pada MMP13 annotation."
  )
}

# ------------------------------------------------------------
# 21. Memastikan MMP13 probe tersedia dalam RMA matrix
# ------------------------------------------------------------

mmp13_probe_ids <- unique(
  as.character(
    mmp13_annotation_rows$probe_id
  )
)

mmp13_probe_ids_present <- intersect(
  mmp13_probe_ids,
  rownames(expression_rma)
)

cat("\nProbe MMP13 dalam RMA matrix:\n")
print(mmp13_probe_ids_present)

if (length(mmp13_probe_ids_present) == 0) {
  stop(
    "Probe MMP13 tidak ditemukan ",
    "pada RMA expression matrix."
  )
}

# ------------------------------------------------------------
# 22. Ekspresi MMP13 per pool
# ------------------------------------------------------------

mmp13_expression_matrix <-
  expression_rma[
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
    "geo_accession"
  ) %>%
  dplyr::left_join(
    sample_metadata %>%
      dplyr::select(
        geo_accession,
        title_original,
        group,
        pool_label,
        individuals_per_pool,
        statistical_unit
      ),
    by =
      "geo_accession"
  ) %>%
  tibble::as_tibble()

cat("\nEkspresi MMP13 setiap RNA pool:\n")

print(
  mmp13_expression_long,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 23. Membuat design matrix awal
# ------------------------------------------------------------

design_group <- model.matrix(
  ~ group,
  data =
    sample_metadata
)

rownames(design_group) <-
  sample_metadata$geo_accession

cat("\nDesign matrix:\n")
print(design_group)

cat("\nDesign dimensions:\n")
print(dim(design_group))

cat("\nDesign rank:\n")
print(qr(design_group)$rank)

stopifnot(
  qr(design_group)$rank ==
    ncol(design_group)
)

stopifnot(
  "groupOA" %in%
    colnames(design_group)
)

# ------------------------------------------------------------
# 24. Membuat preparation summary
# ------------------------------------------------------------

preparation_summary <- tibble::tibble(
  
  metric = c(
    "Dataset",
    "Platform",
    "Raw CEL files",
    "RMA probesets",
    "Number of RNA pools",
    "Normal pools",
    "OA pools",
    "Individuals contributing to each pool",
    "Statistical unit",
    "RMA expression minimum",
    "RMA expression maximum",
    "MMP13 probe number",
    "MMP13 annotation source"
  ),
  
  value = c(
    "GSE169077",
    "GPL96",
    length(cel_files),
    nrow(expression_rma),
    ncol(expression_rma),
    sum(sample_metadata$group == "Normal"),
    sum(sample_metadata$group == "OA"),
    5,
    "RNA pool",
    min(expression_rma),
    max(expression_rma),
    length(mmp13_probe_ids_present),
    annotation_source
  )
)

cat("\nPreparation summary:\n")

print(
  preparation_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 25. Menyimpan objek R
# ------------------------------------------------------------

saveRDS(
  expression_geo_processed,
  file =
    "data_processed/GSE169077_expression_GEO_processed_original.rds"
)

saveRDS(
  expression_rma,
  file =
    "data_processed/GSE169077_expression_RMA_11_pools.rds"
)

saveRDS(
  rma_eset,
  file =
    "data_processed/GSE169077_RMA_ExpressionSet.rds"
)

saveRDS(
  sample_metadata,
  file =
    "data_processed/GSE169077_metadata_11_pools.rds"
)

saveRDS(
  design_group,
  file =
    "data_processed/GSE169077_design_group.rds"
)

saveRDS(
  feature_annotation_original,
  file =
    "data_processed/GSE169077_feature_annotation_original.rds"
)

saveRDS(
  mmp13_annotation_rows,
  file =
    "data_processed/GSE169077_MMP13_probe_annotation.rds"
)

saveRDS(
  mmp13_expression_matrix,
  file =
    "data_processed/GSE169077_MMP13_RMA_expression.rds"
)

saveRDS(
  cel_manifest_ordered,
  file =
    "data_processed/GSE169077_CEL_manifest.rds"
)

# ------------------------------------------------------------
# 26. Menyimpan Excel report
# ------------------------------------------------------------

openxlsx::write.xlsx(
  list(
    
    Preparation_Summary =
      as.data.frame(
        preparation_summary
      ),
    
    Sample_Metadata =
      as.data.frame(
        sample_metadata
      ),
    
    CEL_Manifest =
      as.data.frame(
        cel_manifest_ordered
      ),
    
    RMA_Expression_Summary =
      as.data.frame(
        sample_expression_summary
      ),
    
    MMP13_Annotation =
      as.data.frame(
        mmp13_annotation_rows
      ),
    
    MMP13_RMA_Expression =
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
  ),
  
  file =
    "metadata/GSE169077_preparation_summary.xlsx",
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 27. Session information
# ------------------------------------------------------------

sink(
  "metadata/GSE169077_preparation_sessionInfo.txt"
)

print(sessionInfo())

sink()

# ------------------------------------------------------------
# 28. Pesan akhir
# ------------------------------------------------------------

message("")
message("============================================")
message("GSE169077 PREPARATION SELESAI")
message("Platform                    : GPL96")
message("Raw CEL files               : ",
        length(cel_files))
message("RMA probesets               : ",
        nrow(expression_rma))
message("Total RNA pools             : ",
        ncol(expression_rma))
message("Normal pools                : ",
        sum(sample_metadata$group == "Normal"))
message("OA pools                    : ",
        sum(sample_metadata$group == "OA"))
message("Individuals per pool        : 5")
message("Statistical unit            : RNA pool")
message("RMA expression range        : ",
        round(min(expression_rma), 3),
        " sampai ",
        round(max(expression_rma), 3))
message("MMP13 probes available      : ",
        length(mmp13_probe_ids_present))
message("Design residual df          : ",
        nrow(design_group) -
          qr(design_group)$rank)
message("============================================")