# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 14_prepare_MMP13_meta_analysis.R
# PURPOSE : Harmonize MMP13 effects and create Figure 1B
#
# PRIMARY META-ANALYSIS:
#   GSE114007
#   GSE117999
#   GSE169077
#
# SUPPORTIVE PAIRED EVIDENCE, NOT POOLED:
#   GSE57218
#
# EFFECT DIRECTION:
#   OA minus reference cartilage
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Memuat dan memeriksa paket
# ------------------------------------------------------------

required_packages <- c(
  "metafor",
  "ggplot2",
  "dplyr",
  "tibble",
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
  
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

still_missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(still_missing_packages) > 0) {
  
  stop(
    "Paket berikut masih belum tersedia: ",
    paste(
      still_missing_packages,
      collapse = ", "
    )
  )
}

library(metafor)
library(ggplot2)
library(dplyr)
library(tibble)
library(openxlsx)

message(
  "Semua paket meta-analysis berhasil dimuat."
)

# ------------------------------------------------------------
# 2. Memastikan output folder tersedia
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
# 3. Exact effect files dari analisis sebelumnya
# ------------------------------------------------------------

effect_files <- c(
  
  GSE114007 =
    "results/tables/GSE114007_MMP13_forest_input.xlsx",
  
  GSE117999 =
    "data_processed/GSE117999_MMP13_forest_input.rds",
  
  GSE57218 =
    "data_processed/GSE57218_MMP13_forest_candidate.rds",
  
  GSE169077 =
    "data_processed/GSE169077_MMP13_forest_candidate.rds"
)

cat("\nEffect files yang akan digunakan:\n")
print(effect_files)

missing_effect_files <- effect_files[
  !file.exists(effect_files)
]

if (length(missing_effect_files) > 0) {
  
  stop(
    "File effect berikut tidak ditemukan:\n",
    paste(
      paste0(
        names(missing_effect_files),
        ": ",
        missing_effect_files
      ),
      collapse = "\n"
    ),
    "\n\nPeriksa kembali folder data_processed dan results/tables."
  )
}

message(
  "Seluruh effect files ditemukan."
)

# ------------------------------------------------------------
# 4. Fungsi membaca effect table
# ------------------------------------------------------------

read_effect_table <- function(
    file_path
) {
  
  file_extension <- tolower(
    tools::file_ext(
      file_path
    )
  )
  
  if (file_extension == "rds") {
    
    current_object <- readRDS(
      file_path
    )
    
    if (
      !is.data.frame(current_object) &&
      !is.matrix(current_object)
    ) {
      
      stop(
        "Objek RDS bukan tabel: ",
        file_path
      )
    }
    
    current_table <- as.data.frame(
      current_object,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else if (
    file_extension %in%
    c(
      "xlsx",
      "xlsm"
    )
  ) {
    
    current_table <- openxlsx::read.xlsx(
      file_path,
      sheet = 1,
      check.names = FALSE
    )
    
    current_table <- as.data.frame(
      current_table,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else {
    
    stop(
      "Ekstensi file tidak didukung: ",
      file_path
    )
  }
  
  if (nrow(current_table) == 0) {
    
    stop(
      "Tabel effect kosong: ",
      file_path
    )
  }
  
  current_table
}

# ------------------------------------------------------------
# 5. Fungsi mencari kolom secara aman
# ------------------------------------------------------------

find_first_column <- function(
    input_table,
    candidate_names,
    required = TRUE
) {
  
  available_names <- colnames(
    input_table
  )
  
  matched_names <- candidate_names[
    candidate_names %in%
      available_names
  ]
  
  if (length(matched_names) == 0) {
    
    if (required) {
      
      stop(
        "Tidak ditemukan kolom berikut: ",
        paste(
          candidate_names,
          collapse = ", "
        ),
        "\nKolom yang tersedia: ",
        paste(
          available_names,
          collapse = ", "
        )
      )
    }
    
    return(
      NA_character_
    )
  }
  
  matched_names[1]
}

# ------------------------------------------------------------
# 6. Fungsi mengambil nilai numerik
# ------------------------------------------------------------

extract_numeric_value <- function(
    input_table,
    candidate_names,
    required = TRUE
) {
  
  current_column <- find_first_column(
    input_table = input_table,
    candidate_names = candidate_names,
    required = required
  )
  
  if (is.na(current_column)) {
    
    return(
      NA_real_
    )
  }
  
  current_value <- suppressWarnings(
    as.numeric(
      as.character(
        input_table[[current_column]][1]
      )
    )
  )
  
  if (
    required &&
    (
      length(current_value) != 1 ||
      !is.finite(current_value)
    )
  ) {
    
    stop(
      "Nilai numerik tidak valid pada kolom: ",
      current_column
    )
  }
  
  current_value
}

# ------------------------------------------------------------
# 7. Fungsi mengambil nilai character
# ------------------------------------------------------------

extract_character_value <- function(
    input_table,
    candidate_names,
    default_value = NA_character_
) {
  
  current_column <- find_first_column(
    input_table = input_table,
    candidate_names = candidate_names,
    required = FALSE
  )
  
  if (is.na(current_column)) {
    
    return(
      default_value
    )
  }
  
  current_value <- as.character(
    input_table[[current_column]][1]
  )
  
  if (
    length(current_value) == 0 ||
    is.na(current_value) ||
    current_value == ""
  ) {
    
    return(
      default_value
    )
  }
  
  current_value
}

# ------------------------------------------------------------
# 8. Fungsi standardisasi satu effect row
# ------------------------------------------------------------

standardize_effect_row <- function(
    dataset_id,
    input_table,
    source_file,
    effect_candidates
) {
  
  if (nrow(input_table) != 1) {
    
    stop(
      dataset_id,
      " harus memiliki tepat satu baris effect, tetapi ditemukan ",
      nrow(input_table),
      " baris dalam file:\n",
      source_file
    )
  }
  
  effect_value <- extract_numeric_value(
    input_table,
    effect_candidates
  )
  
  standard_error_value <- extract_numeric_value(
    input_table,
    c(
      "standard_error",
      "standard error",
      "SE",
      "se",
      "stderr"
    )
  )
  
  lower_value <- extract_numeric_value(
    input_table,
    c(
      "CI_95_lower",
      "CI95_lower",
      "ci_lower",
      "lower"
    ),
    required = FALSE
  )
  
  upper_value <- extract_numeric_value(
    input_table,
    c(
      "CI_95_upper",
      "CI95_upper",
      "ci_upper",
      "upper"
    ),
    required = FALSE
  )
  
  if (!is.finite(lower_value)) {
    
    lower_value <-
      effect_value -
      1.96 *
      standard_error_value
  }
  
  if (!is.finite(upper_value)) {
    
    upper_value <-
      effect_value +
      1.96 *
      standard_error_value
  }
  
  p_value <- extract_numeric_value(
    input_table,
    c(
      "PValue",
      "P.Value",
      "P-value",
      "pvalue",
      "PValue_voom"
    ),
    required = FALSE
  )
  
  fdr_value <- extract_numeric_value(
    input_table,
    c(
      "FDR",
      "adj.P.Val",
      "FDR_voom",
      "padj"
    ),
    required = FALSE
  )
  
  tibble::tibble(
    
    dataset =
      dataset_id,
    
    source_file =
      source_file,
    
    source_effect_measure =
      extract_character_value(
        input_table,
        c(
          "effect_measure",
          "effect_scale"
        ),
        default_value =
          "Log2-scale expression difference"
      ),
    
    source_effect_method =
      extract_character_value(
        input_table,
        c(
          "effect_method",
          "analysis_method"
        ),
        default_value =
          NA_character_
      ),
    
    source_direction =
      extract_character_value(
        input_table,
        c(
          "effect_direction",
          "direction"
        ),
        default_value =
          "OA minus reference"
      ),
    
    source_probe =
      extract_character_value(
        input_table,
        c(
          "probe_id",
          "feature_id",
          "agilent_probe"
        ),
        default_value =
          "MMP13"
      ),
    
    effect =
      effect_value,
    
    standard_error =
      standard_error_value,
    
    CI_95_lower =
      lower_value,
    
    CI_95_upper =
      upper_value,
    
    PValue =
      p_value,
    
    FDR =
      fdr_value
  )
}

# ------------------------------------------------------------
# 9. Membaca empat effect tables
# ------------------------------------------------------------

gse114007_table <- read_effect_table(
  effect_files[["GSE114007"]]
)

gse117999_table <- read_effect_table(
  effect_files[["GSE117999"]]
)

gse57218_table <- read_effect_table(
  effect_files[["GSE57218"]]
)

gse169077_table <- read_effect_table(
  effect_files[["GSE169077"]]
)

cat("\nDimensi input tables:\n")

print(
  tibble::tibble(
    
    dataset =
      names(effect_files),
    
    rows = c(
      nrow(gse114007_table),
      nrow(gse117999_table),
      nrow(gse57218_table),
      nrow(gse169077_table)
    ),
    
    columns = c(
      ncol(gse114007_table),
      ncol(gse117999_table),
      ncol(gse57218_table),
      ncol(gse169077_table)
    )
  ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 10. Mengekstrak effect dari setiap dataset
# ------------------------------------------------------------

selected_effects <- dplyr::bind_rows(
  
  standardize_effect_row(
    
    dataset_id =
      "GSE114007",
    
    input_table =
      gse114007_table,
    
    source_file =
      effect_files[["GSE114007"]],
    
    effect_candidates = c(
      "log2FC",
      "effect"
    )
  ),
  
  standardize_effect_row(
    
    dataset_id =
      "GSE117999",
    
    input_table =
      gse117999_table,
    
    source_file =
      effect_files[["GSE117999"]],
    
    effect_candidates = c(
      "log2FC",
      "effect"
    )
  ),
  
  standardize_effect_row(
    
    dataset_id =
      "GSE57218",
    
    input_table =
      gse57218_table,
    
    source_file =
      effect_files[["GSE57218"]],
    
    effect_candidates = c(
      "effect",
      "log2FC"
    )
  ),
  
  standardize_effect_row(
    
    dataset_id =
      "GSE169077",
    
    input_table =
      gse169077_table,
    
    source_file =
      effect_files[["GSE169077"]],
    
    effect_candidates = c(
      "effect",
      "log2FC"
    )
  )
)

cat("\nSelected MMP13 effects:\n")

print(
  selected_effects,
  n = Inf,
  width = Inf
)

stopifnot(
  nrow(selected_effects) == 4
)

stopifnot(
  dplyr::n_distinct(
    selected_effects$dataset
  ) == 4
)

stopifnot(
  all(
    is.finite(
      selected_effects$effect
    )
  )
)

stopifnot(
  all(
    is.finite(
      selected_effects$standard_error
    )
  )
)

stopifnot(
  all(
    selected_effects$standard_error > 0
  )
)

# ------------------------------------------------------------
# 11. Study specification
# ------------------------------------------------------------

study_specification <- tibble::tibble(
  
  dataset = c(
    "GSE114007",
    "GSE117999",
    "GSE57218",
    "GSE169077"
  ),
  
  platform = c(
    "RNA-seq: GPL11154 and GPL18573",
    "Agilent microarray: GPL20844",
    "Illumina microarray: GPL6947",
    "Affymetrix HG-U133A: GPL96"
  ),
  
  tissue = c(
    "Human knee articular cartilage",
    "Human knee cartilage",
    "Human knee articular cartilage",
    "Human articular cartilage"
  ),
  
  comparison = c(
    "OA versus healthy control cartilage",
    "End-stage OA versus non-OA APM cartilage",
    paste(
      "OA-affected versus preserved cartilage",
      "within the same OA joint"
    ),
    "Late-stage OA versus normal cartilage"
  ),
  
  study_design = c(
    "Independent case-control",
    "Independent case-control with covariate adjustment",
    "Paired within-patient",
    "Independent pool-level comparison"
  ),
  
  reference_units = c(
    18L,
    10L,
    33L,
    5L
  ),
  
  OA_units = c(
    20L,
    10L,
    33L,
    6L
  ),
  
  statistical_unit = c(
    "Individual cartilage sample",
    "Individual cartilage sample",
    "Patient pair",
    "RNA pool"
  ),
  
  included_in_primary_meta = c(
    TRUE,
    TRUE,
    FALSE,
    TRUE
  ),
  
  evidence_role = c(
    "Primary independent cohort",
    "Primary independent cohort",
    "Supportive paired cohort",
    "Primary independent cohort"
  )
)

# ------------------------------------------------------------
# 12. Harmonisasi arah effect
# ------------------------------------------------------------

harmonized_effects <- selected_effects %>%
  
  dplyr::left_join(
    study_specification,
    by = "dataset"
  ) %>%
  
  dplyr::mutate(
    
    harmonized_direction =
      "OA minus reference cartilage",
    
    direction_adjustment =
      paste(
        "No sign reversal;",
        "source analyses already used OA minus reference"
      ),
    
    sampling_variance =
      standard_error^2,
    
    approximate_fold_change =
      2^effect
  ) %>%
  
  dplyr::arrange(
    match(
      dataset,
      study_specification$dataset
    )
  )

cat("\nHarmonized MMP13 effects:\n")

print(
  harmonized_effects %>%
    
    dplyr::select(
      dataset,
      comparison,
      study_design,
      statistical_unit,
      reference_units,
      OA_units,
      effect,
      standard_error,
      CI_95_lower,
      CI_95_upper,
      approximate_fold_change,
      PValue,
      FDR,
      included_in_primary_meta,
      source_file
    ),
  
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 13. Validasi broad range berdasarkan hasil sebelumnya
# ------------------------------------------------------------

validation_table <- harmonized_effects %>%
  
  dplyr::transmute(
    
    dataset,
    effect,
    standard_error,
    
    validation_pass =
      dplyr::case_when(
        
        dataset == "GSE114007" ~
          effect > 3.0 &
          effect < 4.0 &
          standard_error > 0.4 &
          standard_error < 0.9,
        
        dataset == "GSE117999" ~
          effect > -0.5 &
          effect < 1.0 &
          standard_error > 0.5 &
          standard_error < 1.5,
        
        dataset == "GSE57218" ~
          effect > -0.5 &
          effect < 0 &
          standard_error > 0.05 &
          standard_error < 0.2,
        
        dataset == "GSE169077" ~
          effect > 0.5 &
          effect < 0.9 &
          standard_error > 0.1 &
          standard_error < 0.3,
        
        TRUE ~
          FALSE
      )
  )

cat("\nEffect validation:\n")

print(
  validation_table,
  n = Inf,
  width = Inf
)

if (!all(validation_table$validation_pass)) {
  
  stop(
    "Satu atau lebih effect tidak cocok dengan hasil sebelumnya. ",
    "Periksa Effect_Validation dan effect input files."
  )
}

# ------------------------------------------------------------
# 14. Memisahkan primary dan supportive evidence
# ------------------------------------------------------------

primary_meta_data <- harmonized_effects %>%
  
  dplyr::filter(
    included_in_primary_meta
  ) %>%
  
  dplyr::arrange(
    match(
      dataset,
      c(
        "GSE114007",
        "GSE117999",
        "GSE169077"
      )
    )
  )

supportive_paired_data <- harmonized_effects %>%
  
  dplyr::filter(
    !included_in_primary_meta
  )

stopifnot(
  nrow(primary_meta_data) == 3
)

stopifnot(
  nrow(supportive_paired_data) == 1
)

stopifnot(
  supportive_paired_data$dataset ==
    "GSE57218"
)

message(
  "Primary meta-analysis berisi tiga independent cohorts."
)

message(
  "GSE57218 disimpan sebagai supportive paired evidence."
)

# ------------------------------------------------------------
# 15. Random-effects REML + Knapp-Hartung
# ------------------------------------------------------------

meta_random_REML_KH <- metafor::rma.uni(
  
  yi =
    effect,
  
  sei =
    standard_error,
  
  data =
    primary_meta_data,
  
  method =
    "REML",
  
  test =
    "knha",
  
  slab =
    dataset,
  
  level =
    95
)

cat("\nRandom-effects REML + Knapp-Hartung:\n")

print(
  meta_random_REML_KH
)

# ------------------------------------------------------------
# 16. Sensitivity models
# ------------------------------------------------------------

meta_common_effect <- metafor::rma.uni(
  
  yi =
    effect,
  
  sei =
    standard_error,
  
  data =
    primary_meta_data,
  
  method =
    "FE",
  
  test =
    "z",
  
  slab =
    dataset,
  
  level =
    95
)

meta_random_REML_z <- metafor::rma.uni(
  
  yi =
    effect,
  
  sei =
    standard_error,
  
  data =
    primary_meta_data,
  
  method =
    "REML",
  
  test =
    "z",
  
  slab =
    dataset,
  
  level =
    95
)

cat("\nCommon-effect sensitivity model:\n")

print(
  meta_common_effect
)

cat("\nRandom-effects REML z-test model:\n")

print(
  meta_random_REML_z
)

# ------------------------------------------------------------
# 17. Prediction interval
# ------------------------------------------------------------

meta_prediction <- predict(
  meta_random_REML_KH,
  level = 95
)

cat("\nRandom-effects prediction:\n")

print(
  meta_prediction
)

prediction_lower <- if (
  !is.null(
    meta_prediction$pi.lb
  )
) {
  
  as.numeric(
    meta_prediction$pi.lb
  )
  
} else {
  
  NA_real_
}

prediction_upper <- if (
  !is.null(
    meta_prediction$pi.ub
  )
) {
  
  as.numeric(
    meta_prediction$pi.ub
  )
  
} else {
  
  NA_real_
}

# ------------------------------------------------------------
# 18. Model summaries
# ------------------------------------------------------------

random_effects_summary <- tibble::tibble(
  
  model =
    "Random-effects REML with Knapp-Hartung",
  
  studies =
    as.integer(
      meta_random_REML_KH$k
    ),
  
  pooled_effect =
    as.numeric(
      meta_random_REML_KH$b[1]
    ),
  
  standard_error =
    as.numeric(
      meta_random_REML_KH$se
    ),
  
  CI_95_lower =
    as.numeric(
      meta_random_REML_KH$ci.lb
    ),
  
  CI_95_upper =
    as.numeric(
      meta_random_REML_KH$ci.ub
    ),
  
  approximate_fold_change =
    2^as.numeric(
      meta_random_REML_KH$b[1]
    ),
  
  degrees_of_freedom =
    as.numeric(
      meta_random_REML_KH$k -
        meta_random_REML_KH$p
    ),
  
  PValue =
    as.numeric(
      meta_random_REML_KH$pval
    ),
  
  tau_squared =
    as.numeric(
      meta_random_REML_KH$tau2
    ),
  
  I_squared_percent =
    as.numeric(
      meta_random_REML_KH$I2
    ),
  
  H_squared =
    as.numeric(
      meta_random_REML_KH$H2
    ),
  
  Cochran_Q =
    as.numeric(
      meta_random_REML_KH$QE
    ),
  
  Cochran_Q_PValue =
    as.numeric(
      meta_random_REML_KH$QEp
    ),
  
  prediction_interval_lower =
    prediction_lower,
  
  prediction_interval_upper =
    prediction_upper
)

common_effect_summary <- tibble::tibble(
  
  model =
    "Common-effect inverse-variance",
  
  studies =
    as.integer(
      meta_common_effect$k
    ),
  
  pooled_effect =
    as.numeric(
      meta_common_effect$b[1]
    ),
  
  standard_error =
    as.numeric(
      meta_common_effect$se
    ),
  
  CI_95_lower =
    as.numeric(
      meta_common_effect$ci.lb
    ),
  
  CI_95_upper =
    as.numeric(
      meta_common_effect$ci.ub
    ),
  
  approximate_fold_change =
    2^as.numeric(
      meta_common_effect$b[1]
    ),
  
  PValue =
    as.numeric(
      meta_common_effect$pval
    )
)

random_z_summary <- tibble::tibble(
  
  model =
    "Random-effects REML with z-test",
  
  studies =
    as.integer(
      meta_random_REML_z$k
    ),
  
  pooled_effect =
    as.numeric(
      meta_random_REML_z$b[1]
    ),
  
  standard_error =
    as.numeric(
      meta_random_REML_z$se
    ),
  
  CI_95_lower =
    as.numeric(
      meta_random_REML_z$ci.lb
    ),
  
  CI_95_upper =
    as.numeric(
      meta_random_REML_z$ci.ub
    ),
  
  PValue =
    as.numeric(
      meta_random_REML_z$pval
    )
)

cat("\nRandom-effects summary:\n")

print(
  random_effects_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 19. Study weights
# ------------------------------------------------------------

random_effect_weights <- as.numeric(
  stats::weights(
    meta_random_REML_KH
  )
)

common_effect_weights <- as.numeric(
  stats::weights(
    meta_common_effect
  )
)

cat("\nRandom-effects weights:\n")
print(random_effect_weights)

cat("\nCommon-effect weights:\n")
print(common_effect_weights)

stopifnot(
  length(random_effect_weights) ==
    nrow(primary_meta_data)
)

stopifnot(
  length(common_effect_weights) ==
    nrow(primary_meta_data)
)

stopifnot(
  abs(
    sum(random_effect_weights) -
      100
  ) < 0.01
)

stopifnot(
  abs(
    sum(common_effect_weights) -
      100
  ) < 0.01
)

study_weights <- primary_meta_data %>%
  
  dplyr::transmute(
    
    dataset =
      dataset,
    
    effect =
      effect,
    
    standard_error =
      standard_error,
    
    sampling_variance =
      sampling_variance,
    
    random_effect_weight_percent =
      random_effect_weights,
    
    common_effect_weight_percent =
      common_effect_weights
  )

cat("\nStudy weights:\n")

print(
  study_weights,
  n = Inf,
  width = Inf
)
# ------------------------------------------------------------
# 20. Leave-one-study-out sensitivity
# ------------------------------------------------------------

meta_leave_one_out <- tryCatch(
  
  {
    
    current_result <- metafor::leave1out(
      meta_random_REML_KH
    )
    
    as.data.frame(
      current_result
    ) %>%
      
      tibble::rownames_to_column(
        "omitted_study"
      ) %>%
      
      tibble::as_tibble()
  },
  
  error = function(error_object) {
    
    warning(
      "Leave-one-study-out gagal: ",
      conditionMessage(
        error_object
      )
    )
    
    tibble::tibble(
      
      note =
        conditionMessage(
          error_object
        )
    )
  }
)

cat("\nLeave-one-study-out results:\n")

print(
  meta_leave_one_out,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 21. Data untuk primary forest plot
# ------------------------------------------------------------

maximum_random_weight <- max(
  study_weights$
    random_effect_weight_percent,
  na.rm = TRUE
)

primary_study_rows <- primary_meta_data %>%
  
  dplyr::left_join(
    
    study_weights %>%
      
      dplyr::select(
        dataset,
        random_effect_weight_percent
      ),
    
    by = "dataset"
  ) %>%
  
  dplyr::mutate(
    
    study_display =
      dplyr::case_when(
        
        dataset ==
          "GSE114007" ~
          "GSE114007 — OA vs healthy control",
        
        dataset ==
          "GSE117999" ~
          "GSE117999 — OA vs non-OA APM",
        
        dataset ==
          "GSE169077" ~
          "GSE169077 — late-stage OA vs normal",
        
        TRUE ~
          dataset
      ),
    
    row_type =
      "Independent cohort",
    
    point_size =
      3 +
      2 *
      sqrt(
        random_effect_weight_percent /
          maximum_random_weight
      ),
    
    line_width =
      0.8,
    
    estimate_text =
      sprintf(
        "%.2f [%.2f, %.2f]",
        effect,
        CI_95_lower,
        CI_95_upper
      )
  ) %>%
  
  dplyr::transmute(
    
    study_display,
    
    estimate =
      effect,
    
    lower =
      CI_95_lower,
    
    upper =
      CI_95_upper,
    
    row_type,
    
    point_size,
    
    line_width,
    
    estimate_text
  )

primary_pooled_row <- tibble::tibble(
  
  study_display =
    "Random-effects pooled estimate",
  
  estimate =
    random_effects_summary$
    pooled_effect,
  
  lower =
    random_effects_summary$
    CI_95_lower,
  
  upper =
    random_effects_summary$
    CI_95_upper,
  
  row_type =
    "Pooled estimate",
  
  point_size =
    5,
  
  line_width =
    1.2,
  
  estimate_text =
    sprintf(
      "%.2f [%.2f, %.2f]",
      random_effects_summary$
        pooled_effect,
      random_effects_summary$
        CI_95_lower,
      random_effects_summary$
        CI_95_upper
    )
)

primary_forest_rows <- dplyr::bind_rows(
  primary_study_rows,
  primary_pooled_row
)

primary_order_top_to_bottom <- c(
  "GSE114007 — OA vs healthy control",
  "GSE117999 — OA vs non-OA APM",
  "GSE169077 — late-stage OA vs normal",
  "Random-effects pooled estimate"
)

primary_forest_rows <- primary_forest_rows %>%
  
  dplyr::mutate(
    
    study_display = factor(
      study_display,
      levels = rev(
        primary_order_top_to_bottom
      )
    ),
    
    row_type = factor(
      row_type,
      levels = c(
        "Independent cohort",
        "Pooled estimate"
      )
    )
  )

primary_plot_range <- range(
  c(
    primary_forest_rows$lower,
    primary_forest_rows$upper,
    0
  ),
  finite = TRUE
)

primary_range_width <- diff(
  primary_plot_range
)

if (
  !is.finite(primary_range_width) ||
  primary_range_width <= 0
) {
  
  primary_range_width <-
    1
}

primary_annotation_x <-
  max(
    primary_forest_rows$upper,
    na.rm = TRUE
  ) +
  0.10 *
  primary_range_width

primary_x_lower <-
  min(
    primary_forest_rows$lower,
    0,
    na.rm = TRUE
  ) -
  0.08 *
  primary_range_width

primary_x_upper <-
  primary_annotation_x +
  0.42 *
  primary_range_width

# ------------------------------------------------------------
# 22. Primary forest plot
# ------------------------------------------------------------

primary_forest_plot <- ggplot(
  
  primary_forest_rows,
  
  aes(
    y = study_display
  )
) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  geom_segment(
    
    aes(
      x = lower,
      xend = upper,
      yend = study_display,
      color = row_type,
      linewidth = line_width
    )
  ) +
  
  geom_point(
    
    aes(
      x = estimate,
      shape = row_type,
      color = row_type,
      size = point_size
    )
  ) +
  
  geom_text(
    
    aes(
      x = primary_annotation_x,
      label = estimate_text
    ),
    
    hjust = 0,
    size = 3.5
  ) +
  
  scale_color_manual(
    
    values = c(
      "Independent cohort" = "black",
      "Pooled estimate" = "#B2182B"
    )
  ) +
  
  scale_shape_manual(
    
    values = c(
      "Independent cohort" = 15,
      "Pooled estimate" = 18
    )
  ) +
  
  scale_size_identity(
    guide = "none"
  ) +
  
  scale_linewidth_identity(
    guide = "none"
  ) +
  
  coord_cartesian(
    
    xlim = c(
      primary_x_lower,
      primary_x_upper
    ),
    
    clip = "off"
  ) +
  
  labs(
    
    title =
      paste(
        "MMP13 meta-analysis across independent",
        "human cartilage cohorts"
      ),
    
    subtitle = paste0(
      "Exploratory random-effects estimate from three independent cohorts ",
      "(I\u00B2 = ",
      round(
        random_effects_summary$I_squared_percent,
        1
      ),
      "%).\n",
      "Paired within-joint evidence is shown separately."
    ),
    
    x =
      expression(
        MMP13~log[2]~
          expression~difference~
          "(OA minus reference)"
      ),
    
    y = NULL,
    
    color =
      "Evidence type",
    
    shape =
      "Evidence type"
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.y =
      element_blank(),
    
    legend.position =
      "top",
    
    axis.text.y =
      element_text(
        size = 9.5
      ),
    
    plot.margin =
      margin(
        10,
        90,
        10,
        10
      )
  )

print(
  primary_forest_plot
)

ggsave(
  
  filename =
    "results/figures/MMP13_primary_meta_analysis_forest.pdf",
  
  plot =
    primary_forest_plot,
  
  width =
    10.5,
  
  height =
    5.8
)

ggsave(
  
  filename =
    "results/figures/MMP13_primary_meta_analysis_forest.tiff",
  
  plot =
    primary_forest_plot,
  
  width =
    10.5,
  
  height =
    5.8,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 23. Data untuk Figure 1B
# ------------------------------------------------------------

# Nama section dibuat singkat agar facet strip tidak terlalu lebar
independent_section_name <-
  "Independent cohorts"

supportive_section_name <-
  "Paired evidence"

# ------------------------------------------------------------
# 23A. Baris tiga independent cohorts
# ------------------------------------------------------------

figure1b_independent_rows <- primary_study_rows %>%
  dplyr::mutate(
    
    section =
      independent_section_name,
    
    row_order =
      dplyr::case_when(
        
        study_display ==
          "GSE114007 — OA vs healthy control" ~
          1,
        
        study_display ==
          "GSE117999 — OA vs non-OA APM" ~
          2,
        
        study_display ==
          "GSE169077 — late-stage OA vs normal" ~
          3,
        
        TRUE ~
          99
      )
  )

# ------------------------------------------------------------
# 23B. Baris pooled estimate
# ------------------------------------------------------------

figure1b_pooled_row <- primary_pooled_row %>%
  dplyr::mutate(
    
    section =
      independent_section_name,
    
    row_order =
      4
  )

# ------------------------------------------------------------
# 23C. Baris supportive paired cohort GSE57218
# ------------------------------------------------------------

figure1b_supportive_row <- supportive_paired_data %>%
  dplyr::transmute(
    
    study_display =
      "GSE57218 — OA-affected vs preserved",
    
    estimate =
      effect,
    
    lower =
      CI_95_lower,
    
    upper =
      CI_95_upper,
    
    row_type =
      "Paired supportive cohort",
    
    point_size =
      4,
    
    line_width =
      0.8,
    
    estimate_text =
      sprintf(
        "%.2f [%.2f, %.2f]",
        effect,
        CI_95_lower,
        CI_95_upper
      ),
    
    section =
      supportive_section_name,
    
    row_order =
      5
  )

# ------------------------------------------------------------
# 23D. Menggabungkan seluruh evidence
# ------------------------------------------------------------

figure1b_rows <- dplyr::bind_rows(
  
  figure1b_independent_rows,
  
  figure1b_pooled_row,
  
  figure1b_supportive_row
)

# Urutan yang diinginkan dari atas ke bawah
figure1b_order_top_to_bottom <- c(
  
  "GSE114007 — OA vs healthy control",
  
  "GSE117999 — OA vs non-OA APM",
  
  "GSE169077 — late-stage OA vs normal",
  
  "Random-effects pooled estimate",
  
  "GSE57218 — OA-affected vs preserved"
)

figure1b_rows <- figure1b_rows %>%
  dplyr::mutate(
    
    study_display = factor(
      study_display,
      levels = rev(
        figure1b_order_top_to_bottom
      )
    ),
    
    row_type = factor(
      row_type,
      levels = c(
        "Independent cohort",
        "Pooled estimate",
        "Paired supportive cohort"
      )
    ),
    
    section = factor(
      section,
      levels = c(
        independent_section_name,
        supportive_section_name
      )
    )
  )

cat("\nData Figure 1B:\n")

print(
  figure1b_rows %>%
    dplyr::select(
      section,
      study_display,
      estimate,
      lower,
      upper,
      row_type,
      estimate_text
    ) %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 23E. Menentukan rentang sumbu dan posisi annotation
# ------------------------------------------------------------

figure1b_plot_range <- range(
  c(
    figure1b_rows$lower,
    figure1b_rows$upper,
    0
  ),
  finite = TRUE
)

figure1b_range_width <- diff(
  figure1b_plot_range
)

if (
  !is.finite(
    figure1b_range_width
  ) ||
  figure1b_range_width <= 0
) {
  
  figure1b_range_width <-
    1
}

# Posisi teks nilai effect di sebelah kanan CI
figure1b_annotation_x <-
  max(
    figure1b_rows$upper,
    na.rm = TRUE
  ) +
  0.10 *
  figure1b_range_width

# Batas kiri plot
figure1b_x_lower <-
  min(
    figure1b_rows$lower,
    0,
    na.rm = TRUE
  ) -
  0.08 *
  figure1b_range_width

# Batas kanan plot, termasuk ruang untuk annotation
figure1b_x_upper <-
  figure1b_annotation_x +
  0.42 *
  figure1b_range_width

cat("\nFigure 1B x-axis limits:\n")

print(
  c(
    lower =
      figure1b_x_lower,
    
    upper =
      figure1b_x_upper,
    
    annotation =
      figure1b_annotation_x
  )
)

# ------------------------------------------------------------
# 24. Figure 1B forest plot
# ------------------------------------------------------------

figure1b_plot <- ggplot(
  
  figure1b_rows,
  
  aes(
    y = study_display
  )
) +
  
  # Garis nol
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey25"
  ) +
  
  # Confidence intervals
  geom_segment(
    
    aes(
      x = lower,
      xend = upper,
      yend = study_display,
      color = row_type,
      linewidth = line_width
    ),
    
    lineend = "butt"
  ) +
  
  # Titik effect estimate
  geom_point(
    
    aes(
      x = estimate,
      shape = row_type,
      color = row_type,
      size = point_size
    )
  ) +
  
  # Nilai effect dan 95% CI
  geom_text(
    
    aes(
      x = figure1b_annotation_x,
      label = estimate_text
    ),
    
    hjust = 0,
    size = 4.1,
    color = "black"
  ) +
  
  # Memisahkan independent evidence dan paired evidence
  facet_grid(
    
    rows =
      vars(
        section
      ),
    
    scales =
      "free_y",
    
    space =
      "free_y"
  ) +
  
  # Warna evidence
  scale_color_manual(
    
    values = c(
      
      "Independent cohort" =
        "black",
      
      "Pooled estimate" =
        "#B2182B",
      
      "Paired supportive cohort" =
        "#2166AC"
    ),
    
    drop = FALSE
  ) +
  
  # Bentuk titik
  scale_shape_manual(
    
    values = c(
      
      "Independent cohort" =
        15,
      
      "Pooled estimate" =
        18,
      
      "Paired supportive cohort" =
        17
    ),
    
    drop = FALSE
  ) +
  
  # Ukuran titik sudah ditentukan di tabel
  scale_size_identity(
    guide = "none"
  ) +
  
  # Ketebalan CI sudah ditentukan di tabel
  scale_linewidth_identity(
    guide = "none"
  ) +
  
  # Menyediakan ruang untuk teks CI di sebelah kanan
  coord_cartesian(
    
    xlim = c(
      figure1b_x_lower,
      figure1b_x_upper
    ),
    
    clip = "off"
  ) +
  
  labs(
    
    title =
      paste(
        "Figure 1B | MMP13 expression effects",
        "across human cartilage cohorts"
      ),
    
    subtitle =
      paste0(
        "Exploratory random-effects estimate from three ",
        "independent cohorts; I² = ",
        round(
          random_effects_summary$
            I_squared_percent,
          1
        ),
        "%. Paired within-joint evidence is shown separately."
      ),
    
    x =
      expression(
        MMP13~log[2]~
          expression~difference~
          "(OA minus reference)"
      ),
    
    y =
      NULL,
    
    color =
      "Evidence type",
    
    shape =
      "Evidence type",
    
    caption =
      paste0(
        "Random-effects pooled estimate = ",
        sprintf(
          "%.2f",
          random_effects_summary$
            pooled_effect
        ),
        " log2 units (95% CI: ",
        sprintf(
          "%.2f",
          random_effects_summary$
            CI_95_lower
        ),
        " to ",
        sprintf(
          "%.2f",
          random_effects_summary$
            CI_95_upper
        ),
        "). The wide interval reflects substantial ",
        "between-cohort heterogeneity."
      )
  ) +
  
  theme_bw(
    base_size = 13
  ) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.y =
      element_blank(),
    
    # Facet strip lebih sempit
    strip.background =
      element_rect(
        fill = "grey95",
        color = "grey40"
      ),
    
    strip.text.y =
      element_text(
        angle = 90,
        face = "bold",
        size = 9.5,
        margin = margin(
          t = 4,
          r = 4,
          b = 4,
          l = 4
        )
      ),
    
    legend.position =
      "top",
    
    legend.title =
      element_text(
        face = "bold"
      ),
    
    axis.text.y =
      element_text(
        size = 11.5,
        angle = 25,
        hjust = 1,
        vjust = 0.5,
        lineheight = 0.90,
        color = "black"
      ),
    
    axis.title.x =
      element_text(
        size = 11.5,
        color = "black",
        margin = margin(t = 7)
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 15
      ),
    
    plot.subtitle =
      element_text(
        size = 11.5,
        lineheight = 1.05
      ),
    
    plot.caption =
      element_text(
        hjust = 0,
        size = 9.5,
        lineheight = 1.05,
        color = "grey30",
        margin = margin(
          t = 8
        )
      ),
    
    plot.margin =
      margin(
        t = 12,
        r = 95,
        b = 12,
        l = 6
      )
  )

print(
  figure1b_plot
)

# ------------------------------------------------------------
# 24A. Menyimpan Figure 1B sebagai PDF
# ------------------------------------------------------------

ggsave(
  
  filename =
    "results/figures/Figure1B_MMP13_cross_cohort_forest.pdf",
  
  plot =
    figure1b_plot,
  
  width =
    11,
  
  height =
    7.2
)

# ------------------------------------------------------------
# 24B. Menyimpan Figure 1B sebagai TIFF 600 dpi
# ------------------------------------------------------------

ggsave(
  
  filename =
    "results/figures/Figure1B_MMP13_cross_cohort_forest.tiff",
  
  plot =
    figure1b_plot,
  
  width =
    11,
  
  height =
    7.2,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

message(
  "Figure 1B berhasil diperbarui dan disimpan."
)
# ------------------------------------------------------------
# 25. Methodological notes
# ------------------------------------------------------------

methodological_notes <- tibble::tibble(
  
  topic = c(
    "Primary pooled datasets",
    "Supportive paired dataset",
    "Effect direction",
    "Effect scale",
    "Random-effects estimator",
    "Inference method",
    "Small number of studies",
    "Publication-bias tests"
  ),
  
  decision = c(
    
    paste(
      "GSE114007, GSE117999,",
      "and GSE169077"
    ),
    
    paste(
      "GSE57218 is shown separately because it compares",
      "OA-affected with preserved cartilage",
      "from the same OA joints"
    ),
    
    paste(
      "All effects are oriented as",
      "OA minus reference cartilage"
    ),
    
    paste(
      "Study-specific log2-scale expression differences;",
      "cross-platform pooled interpretation",
      "must remain cautious"
    ),
    
    "Restricted maximum likelihood (REML)",
    
    "Knapp-Hartung adjustment",
    
    paste(
      "Heterogeneity and prediction interval",
      "are interpreted cautiously because only",
      "three independent cohorts are pooled"
    ),
    
    paste(
      "Funnel plot and Egger regression are not performed",
      "because only three independent cohorts are available"
    )
  )
)

# ------------------------------------------------------------
# 26. Preparation summary
# ------------------------------------------------------------

preparation_summary <- tibble::tibble(
  
  metric = c(
    "Total MMP13 datasets",
    "Datasets in primary meta-analysis",
    "Supportive paired datasets",
    "Random-effects estimator",
    "Inference method",
    "Primary pooled log2 effect",
    "Primary pooled approximate fold change",
    "Primary pooled 95% CI lower",
    "Primary pooled 95% CI upper",
    "Primary pooled P-value",
    "Tau squared",
    "I squared percent",
    "Cochran Q",
    "Cochran Q P-value",
    "Prediction interval lower",
    "Prediction interval upper"
  ),
  
  value = c(
    nrow(
      harmonized_effects
    ),
    nrow(
      primary_meta_data
    ),
    nrow(
      supportive_paired_data
    ),
    "REML",
    "Knapp-Hartung",
    random_effects_summary$
      pooled_effect,
    random_effects_summary$
      approximate_fold_change,
    random_effects_summary$
      CI_95_lower,
    random_effects_summary$
      CI_95_upper,
    random_effects_summary$
      PValue,
    random_effects_summary$
      tau_squared,
    random_effects_summary$
      I_squared_percent,
    random_effects_summary$
      Cochran_Q,
    random_effects_summary$
      Cochran_Q_PValue,
    random_effects_summary$
      prediction_interval_lower,
    random_effects_summary$
      prediction_interval_upper
  )
)

cat("\nMeta-analysis preparation summary:\n")

print(
  preparation_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 27. Source audit
# ------------------------------------------------------------

source_audit <- harmonized_effects %>%
  
  dplyr::select(
    dataset,
    source_file,
    source_effect_measure,
    source_effect_method,
    source_direction,
    source_probe,
    effect,
    standard_error,
    CI_95_lower,
    CI_95_upper,
    PValue,
    FDR
  )

# ------------------------------------------------------------
# 28. Menyimpan R objects
# ------------------------------------------------------------

saveRDS(
  
  selected_effects,
  
  file =
    "data_processed/MMP13_selected_effects.rds"
)

saveRDS(
  
  source_audit,
  
  file =
    "data_processed/MMP13_effect_source_audit.rds"
)

saveRDS(
  
  harmonized_effects,
  
  file =
    "data_processed/MMP13_harmonized_effects.rds"
)

saveRDS(
  
  primary_meta_data,
  
  file =
    "data_processed/MMP13_primary_meta_input.rds"
)

saveRDS(
  
  supportive_paired_data,
  
  file =
    "data_processed/MMP13_supportive_paired_input.rds"
)

saveRDS(
  
  meta_random_REML_KH,
  
  file =
    "data_processed/MMP13_meta_random_REML_KH.rds"
)

saveRDS(
  
  meta_common_effect,
  
  file =
    "data_processed/MMP13_meta_common_effect.rds"
)

saveRDS(
  
  meta_random_REML_z,
  
  file =
    "data_processed/MMP13_meta_random_REML_z.rds"
)

saveRDS(
  
  meta_leave_one_out,
  
  file =
    "data_processed/MMP13_meta_leave_one_out.rds"
)

saveRDS(
  
  random_effects_summary,
  
  file =
    "data_processed/MMP13_meta_random_summary.rds"
)

# ------------------------------------------------------------
# 29. Menyimpan CSV
# ------------------------------------------------------------

utils::write.csv(
  
  harmonized_effects,
  
  file =
    "results/tables/MMP13_harmonized_effects.csv",
  
  row.names =
    FALSE
)

utils::write.csv(
  
  study_weights,
  
  file =
    "results/tables/MMP13_meta_study_weights.csv",
  
  row.names =
    FALSE
)

# ------------------------------------------------------------
# 30. Menyimpan Excel workbook
# ------------------------------------------------------------

openxlsx::write.xlsx(
  
  list(
    
    Preparation_Summary =
      as.data.frame(
        preparation_summary
      ),
    
    Harmonized_Effects =
      as.data.frame(
        harmonized_effects
      ),
    
    Effect_Validation =
      as.data.frame(
        validation_table
      ),
    
    Primary_Meta_Input =
      as.data.frame(
        primary_meta_data
      ),
    
    Supportive_Paired =
      as.data.frame(
        supportive_paired_data
      ),
    
    Random_Effects_Summary =
      as.data.frame(
        random_effects_summary
      ),
    
    Common_Effect_Summary =
      as.data.frame(
        common_effect_summary
      ),
    
    Random_Z_Sensitivity =
      as.data.frame(
        random_z_summary
      ),
    
    Study_Weights =
      as.data.frame(
        study_weights
      ),
    
    Leave_One_Out =
      as.data.frame(
        meta_leave_one_out
      ),
    
    Source_Audit =
      as.data.frame(
        source_audit
      ),
    
    Method_Notes =
      as.data.frame(
        methodological_notes
      )
  ),
  
  file =
    "results/tables/MMP13_meta_analysis_preparation.xlsx",
  
  overwrite =
    TRUE
)

# ------------------------------------------------------------
# 31. Menyimpan model output text
# ------------------------------------------------------------

sink(
  "results/tables/MMP13_meta_analysis_model_output.txt"
)

cat(
  "PRIMARY RANDOM-EFFECTS MODEL\n"
)

print(
  meta_random_REML_KH
)

cat(
  "\n\nCOMMON-EFFECT SENSITIVITY MODEL\n"
)

print(
  meta_common_effect
)

cat(
  "\n\nRANDOM-EFFECTS Z-TEST SENSITIVITY MODEL\n"
)

print(
  meta_random_REML_z
)

cat(
  "\n\nPREDICTION INTERVAL\n"
)

print(
  meta_prediction
)

cat(
  "\n\nLEAVE-ONE-STUDY-OUT\n"
)

print(
  meta_leave_one_out
)

sink()

# ------------------------------------------------------------
# 32. Session information
# ------------------------------------------------------------

sink(
  "results/tables/MMP13_meta_analysis_sessionInfo.txt"
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
message("MMP13 META-ANALYSIS PREPARATION SELESAI")
message("================================================")

message(
  "Total datasets                  : ",
  nrow(
    harmonized_effects
  )
)

message(
  "Primary independent cohorts     : ",
  nrow(
    primary_meta_data
  )
)

message(
  "Supportive paired cohorts       : ",
  nrow(
    supportive_paired_data
  )
)

message(
  "Primary pooled log2 effect      : ",
  round(
    random_effects_summary$
      pooled_effect,
    4
  )
)

message(
  "Primary pooled 95% CI           : ",
  round(
    random_effects_summary$
      CI_95_lower,
    4
  ),
  " to ",
  round(
    random_effects_summary$
      CI_95_upper,
    4
  )
)

message(
  "Approximate pooled fold change  : ",
  round(
    random_effects_summary$
      approximate_fold_change,
    3
  )
)

message(
  "Primary pooled P-value          : ",
  format(
    random_effects_summary$
      PValue,
    scientific = TRUE,
    digits = 4
  )
)

message(
  "Tau squared                     : ",
  round(
    random_effects_summary$
      tau_squared,
    4
  )
)

message(
  "I squared                       : ",
  round(
    random_effects_summary$
      I_squared_percent,
    2
  ),
  "%"
)

message(
  "Prediction interval             : ",
  round(
    random_effects_summary$
      prediction_interval_lower,
    4
  ),
  " to ",
  round(
    random_effects_summary$
      prediction_interval_upper,
    4
  )
)

message(
  "GSE57218 included in pooled model: FALSE"
)

message(
  "Figure 1B saved                : ",
  "results/figures/Figure1B_MMP13_cross_cohort_forest"
)

message("================================================")