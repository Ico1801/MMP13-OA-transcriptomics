# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 15_prepare_cross_cohort_pathway_enrichment.R
# PURPOSE : Cross-cohort preranked GSEA and Figure 1C
#
# PRIMARY INDEPENDENT COHORTS:
#   GSE114007
#   GSE117999
#   GSE169077
#
# SUPPORTIVE PAIRED COHORT:
#   GSE57218
#
# GENE-SET DATABASES:
#   MSigDB Hallmark
#   MSigDB Reactome
#
# FIGURE 1C:
#   Prespecified OA/MMP13-related pathway heatmap
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Package installation and loading
# ------------------------------------------------------------

cran_packages <- c(
  "msigdbr",
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "openxlsx"
)

bioconductor_packages <- c(
  "clusterProfiler"
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

library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(openxlsx)

message(
  "Semua paket pathway enrichment berhasil dimuat."
)

# ------------------------------------------------------------
# 2. Check current msigdbr interface
# ------------------------------------------------------------

msigdbr_arguments <- names(
  formals(
    msigdbr::msigdbr
  )
)

if (
  !"collection" %in%
  msigdbr_arguments
) {
  
  stop(
    "Versi msigdbr terlalu lama. ",
    "Jalankan install.packages('msigdbr') ",
    "kemudian restart RStudio."
  )
}

cat("\nVersi paket:\n")

print(
  tibble::tibble(
    
    package = c(
      "clusterProfiler",
      "msigdbr"
    ),
    
    version = c(
      as.character(
        packageVersion(
          "clusterProfiler"
        )
      ),
      as.character(
        packageVersion(
          "msigdbr"
        )
      )
    )
  ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 3. Output folders
# ------------------------------------------------------------

main_figure_folder <-
  "results/figures/Figure 1"

supplementary_figure_folder <-
  paste0(
    "results/figures/",
    "Figure 1 Suplementary/",
    "S10-S11"
  )

table_folder <-
  "results/tables"

processed_folder <-
  "data_processed"

required_folders <- c(
  main_figure_folder,
  supplementary_figure_folder,
  table_folder,
  processed_folder
)

for (current_folder in required_folders) {
  
  dir.create(
    current_folder,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------
# 4. Preferred ranked-gene files
# ------------------------------------------------------------

preferred_ranked_files <- c(
  
  GSE114007 =
    "data_processed/GSE114007_GSEA_ranked_genes.rds",
  
  GSE117999 =
    "data_processed/GSE117999_GSEA_ranked_genes.rds",
  
  GSE57218 =
    "data_processed/GSE57218_GSEA_ranked_genes.rds",
  
  GSE169077 =
    "data_processed/GSE169077_GSEA_ranked_genes.rds"
)

# ------------------------------------------------------------
# 5. Function to locate a ranked-gene file
# ------------------------------------------------------------

locate_ranked_file <- function(
    dataset_id,
    preferred_file
) {
  
  if (
    file.exists(
      preferred_file
    )
  ) {
    
    return(
      preferred_file
    )
  }
  
  search_pattern <- paste0(
    dataset_id,
    ".*(GSEA|RANKED).*(RDS|CSV|XLSX)$"
  )
  
  candidate_files <- unique(
    c(
      
      list.files(
        "data_processed",
        pattern = search_pattern,
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
      ),
      
      list.files(
        "results/tables",
        pattern = search_pattern,
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
      )
    )
  )
  
  candidate_files <- candidate_files[
    !stringr::str_detect(
      stringr::str_to_upper(
        basename(
          candidate_files
        )
      ),
      "SESSION|SUMMARY|META"
    )
  ]
  
  if (
    length(
      candidate_files
    ) == 0
  ) {
    
    stop(
      "Ranked-gene file tidak ditemukan untuk ",
      dataset_id,
      ".\nExpected file: ",
      preferred_file
    )
  }
  
  file_names_upper <- stringr::str_to_upper(
    basename(
      candidate_files
    )
  )
  
  selection_score <- rep(
    0,
    length(
      candidate_files
    )
  )
  
  selection_score <- selection_score +
    ifelse(
      stringr::str_detect(
        file_names_upper,
        "GSEA_RANKED_GENES"
      ),
      100,
      0
    )
  
  selection_score <- selection_score +
    ifelse(
      tools::file_ext(
        candidate_files
      ) %in%
        c(
          "rds",
          "RDS"
        ),
      20,
      0
    )
  
  selection_score <- selection_score +
    ifelse(
      stringr::str_detect(
        file_names_upper,
        "RANKED"
      ),
      10,
      0
    )
  
  selected_index <- order(
    selection_score,
    decreasing = TRUE
  )[1]
  
  cat(
    "\nCandidate ranked files untuk ",
    dataset_id,
    ":\n",
    sep = ""
  )
  
  print(
    candidate_files
  )
  
  message(
    "File yang dipilih untuk ",
    dataset_id,
    ": ",
    candidate_files[selected_index]
  )
  
  candidate_files[selected_index]
}

# ------------------------------------------------------------
# 6. Locate all four ranked-gene files
# ------------------------------------------------------------

ranked_files <- c(
  
  GSE114007 =
    locate_ranked_file(
      "GSE114007",
      preferred_ranked_files[["GSE114007"]]
    ),
  
  GSE117999 =
    locate_ranked_file(
      "GSE117999",
      preferred_ranked_files[["GSE117999"]]
    ),
  
  GSE57218 =
    locate_ranked_file(
      "GSE57218",
      preferred_ranked_files[["GSE57218"]]
    ),
  
  GSE169077 =
    locate_ranked_file(
      "GSE169077",
      preferred_ranked_files[["GSE169077"]]
    )
)

cat("\nRanked-gene files yang digunakan:\n")

print(
  ranked_files
)

stopifnot(
  all(
    file.exists(
      ranked_files
    )
  )
)

# ------------------------------------------------------------
# 7. Read ranked-gene table
# ------------------------------------------------------------

read_ranked_file <- function(
    file_path
) {
  
  file_extension <- tolower(
    tools::file_ext(
      file_path
    )
  )
  
  if (
    file_extension ==
    "rds"
  ) {
    
    current_object <- readRDS(
      file_path
    )
    
    if (
      is.numeric(current_object) &&
      !is.null(
        names(
          current_object
        )
      )
    ) {
      
      return(
        tibble::tibble(
          
          gene =
            names(
              current_object
            ),
          
          rank_metric =
            as.numeric(
              current_object
            )
        )
      )
    }
    
    if (
      is.data.frame(current_object) ||
      is.matrix(current_object)
    ) {
      
      return(
        as.data.frame(
          current_object,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      )
    }
    
    stop(
      "Objek RDS bukan ranked table/vector: ",
      file_path
    )
  }
  
  if (
    file_extension ==
    "csv"
  ) {
    
    return(
      utils::read.csv(
        file_path,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }
  
  if (
    file_extension ==
    "xlsx"
  ) {
    
    return(
      openxlsx::read.xlsx(
        file_path,
        sheet = 1,
        check.names = FALSE
      )
    )
  }
  
  stop(
    "Format ranked-gene file tidak didukung: ",
    file_path
  )
}

# ------------------------------------------------------------
# 8. Column-name utilities
# ------------------------------------------------------------

normalize_column_name <- function(x) {
  
  x <- tolower(
    as.character(x)
  )
  
  gsub(
    "[^a-z0-9]+",
    "",
    x
  )
}

find_first_column <- function(
    input_table,
    candidate_names,
    required = TRUE
) {
  
  available_names <- colnames(
    input_table
  )
  
  normalized_available <-
    normalize_column_name(
      available_names
    )
  
  normalized_candidates <-
    normalize_column_name(
      candidate_names
    )
  
  for (
    current_candidate in
    normalized_candidates
  ) {
    
    current_match <- which(
      normalized_available ==
        current_candidate
    )
    
    if (
      length(
        current_match
      ) > 0
    ) {
      
      return(
        available_names[
          current_match[1]
        ]
      )
    }
  }
  
  if (required) {
    
    stop(
      "Kolom tidak ditemukan. Candidate columns: ",
      paste(
        candidate_names,
        collapse = ", "
      ),
      "\nAvailable columns: ",
      paste(
        available_names,
        collapse = ", "
      )
    )
  }
  
  NA_character_
}

# ------------------------------------------------------------
# 9. Prepare a unique named ranked vector
# ------------------------------------------------------------

prepare_ranked_vector <- function(
    input_table,
    dataset_id,
    source_file
) {
  
  input_table <- as.data.frame(
    input_table,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  original_rows <- nrow(
    input_table
  )
  
  gene_column <- find_first_column(
    input_table,
    c(
      "gene",
      "gene_symbol",
      "symbol",
      "Gene",
      "GeneSymbol",
      "feature"
    )
  )
  
  metric_column <- find_first_column(
    input_table,
    c(
      "moderated_t",
      "signed_sqrt_QLF",
      "signed_sqrt_F",
      "signed_statistic",
      "rank_metric",
      "statistic",
      "stat",
      "t",
      "logFC"
    ),
    required = FALSE
  )
  
  ranking_method <- metric_column
  
  if (
    !is.na(
      metric_column
    )
  ) {
    
    rank_metric <- suppressWarnings(
      as.numeric(
        as.character(
          input_table[[metric_column]]
        )
      )
    )
    
  } else {
    
    logfc_column <- find_first_column(
      input_table,
      c(
        "logFC",
        "log2FC",
        "effect"
      )
    )
    
    pvalue_column <- find_first_column(
      input_table,
      c(
        "PValue",
        "P.Value",
        "pvalue",
        "p_value"
      )
    )
    
    logfc_value <- suppressWarnings(
      as.numeric(
        as.character(
          input_table[[logfc_column]]
        )
      )
    )
    
    pvalue_value <- suppressWarnings(
      as.numeric(
        as.character(
          input_table[[pvalue_column]]
        )
      )
    )
    
    rank_metric <-
      sign(
        logfc_value
      ) *
      -log10(
        pmax(
          pvalue_value,
          1e-300
        )
      )
    
    ranking_method <-
      paste0(
        "sign(",
        logfc_column,
        ") × -log10(",
        pvalue_column,
        ")"
      )
  }
  
  ranked_table <- tibble::tibble(
    
    gene =
      toupper(
        trimws(
          as.character(
            input_table[[gene_column]]
          )
        )
      ),
    
    rank_metric =
      rank_metric
  )
  
  valid_gene <- (
    !is.na(
      ranked_table$gene
    ) &
      ranked_table$gene != "" &
      is.finite(
        ranked_table$rank_metric
      )
  )
  
  ambiguous_gene <- stringr::str_detect(
    ranked_table$gene,
    "///|//|;|\\||,"
  )
  
  ranked_table <- ranked_table[
    valid_gene &
      !ambiguous_gene,
    ,
    drop = FALSE
  ]
  
  valid_rows <- nrow(
    ranked_table
  )
  
  ranked_table <- ranked_table %>%
    
    dplyr::arrange(
      dplyr::desc(
        abs(
          rank_metric
        )
      ),
      gene
    ) %>%
    
    dplyr::distinct(
      gene,
      .keep_all = TRUE
    )
  
  unique_genes <- nrow(
    ranked_table
  )
  
  number_of_tied_statistics <- sum(
    duplicated(
      ranked_table$rank_metric
    )
  )
  
  # Deterministic tiny tie-breaking adjustment.
  # Magnitude is too small to change biological interpretation.
  ranked_table <- ranked_table %>%
    
    dplyr::arrange(
      gene
    ) %>%
    
    dplyr::mutate(
      
      tie_break_offset =
        dplyr::row_number() *
        1e-12,
      
      rank_metric_final =
        rank_metric +
        tie_break_offset
    )
  
  ranked_vector <- ranked_table$
    rank_metric_final
  
  names(
    ranked_vector
  ) <- ranked_table$gene
  
  ranked_vector <- sort(
    ranked_vector,
    decreasing = TRUE
  )
  
  mmp13_rank_metric <- if (
    "MMP13" %in%
    ranked_table$gene
  ) {
    
    ranked_table$rank_metric[
      match(
        "MMP13",
        ranked_table$gene
      )
    ]
    
  } else {
    
    NA_real_
  }
  
  input_summary <- tibble::tibble(
    
    dataset =
      dataset_id,
    
    source_file =
      source_file,
    
    gene_column =
      gene_column,
    
    ranking_method =
      ranking_method,
    
    original_rows =
      original_rows,
    
    valid_rows =
      valid_rows,
    
    unique_genes =
      unique_genes,
    
    removed_or_duplicate_rows =
      original_rows -
      unique_genes,
    
    tied_statistics_before_adjustment =
      number_of_tied_statistics,
    
    minimum_rank_metric =
      min(
        ranked_table$rank_metric,
        na.rm = TRUE
      ),
    
    maximum_rank_metric =
      max(
        ranked_table$rank_metric,
        na.rm = TRUE
      ),
    
    MMP13_present =
      "MMP13" %in%
      ranked_table$gene,
    
    MMP13_rank_metric =
      mmp13_rank_metric
  )
  
  list(
    
    vector =
      ranked_vector,
    
    table =
      ranked_table,
    
    summary =
      input_summary
  )
}

# ------------------------------------------------------------
# 10. Prepare ranked vectors for all datasets
# ------------------------------------------------------------

ranked_inputs <- list()
ranked_input_summaries <- list()

for (
  current_dataset in
  names(
    ranked_files
  )
) {
  
  message(
    "Menyiapkan ranked vector: ",
    current_dataset
  )
  
  current_table <- read_ranked_file(
    ranked_files[[current_dataset]]
  )
  
  current_prepared <- prepare_ranked_vector(
    
    input_table =
      current_table,
    
    dataset_id =
      current_dataset,
    
    source_file =
      ranked_files[[current_dataset]]
  )
  
  ranked_inputs[[current_dataset]] <-
    current_prepared
  
  ranked_input_summaries[[current_dataset]] <-
    current_prepared$summary
  
  saveRDS(
    
    current_prepared$vector,
    
    file = file.path(
      processed_folder,
      paste0(
        current_dataset,
        "_standardized_GSEA_vector.rds"
      )
    )
  )
  
  utils::write.csv(
    
    current_prepared$table,
    
    file = file.path(
      table_folder,
      paste0(
        current_dataset,
        "_standardized_GSEA_ranking.csv"
      )
    ),
    
    row.names = FALSE
  )
}

ranked_input_summary <- dplyr::bind_rows(
  ranked_input_summaries
)

cat("\nRanked input summary:\n")

print(
  ranked_input_summary %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

if (
  any(
    !ranked_input_summary$
    MMP13_present
  )
) {
  
  warning(
    "MMP13 tidak ditemukan pada satu atau lebih ranked lists."
  )
}

# ------------------------------------------------------------
# 11. Inspect available MSigDB collections
# ------------------------------------------------------------

msigdb_collections <- msigdbr::msigdbr_collections(
  db_species = "HS"
)

cat("\nAvailable MSigDB collections:\n")

print(
  msigdb_collections %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

hallmark_available <- any(
  msigdb_collections$
    gs_collection ==
    "H"
)

reactome_available <- any(
  msigdb_collections$
    gs_collection ==
    "C2" &
    msigdb_collections$
    gs_subcollection ==
    "CP:REACTOME"
)

if (!hallmark_available) {
  
  stop(
    "MSigDB Hallmark collection tidak tersedia."
  )
}

if (!reactome_available) {
  
  stop(
    "MSigDB Reactome subcollection tidak tersedia."
  )
}

# ------------------------------------------------------------
# 12. Retrieve Hallmark and Reactome gene sets
# ------------------------------------------------------------

message(
  "Mengambil MSigDB Hallmark gene sets..."
)

hallmark_database <- msigdbr::msigdbr(
  
  db_species =
    "HS",
  
  species =
    "Homo sapiens",
  
  collection =
    "H"
)

message(
  "Mengambil MSigDB Reactome gene sets..."
)

reactome_database <- msigdbr::msigdbr(
  
  db_species =
    "HS",
  
  species =
    "Homo sapiens",
  
  collection =
    "C2",
  
  subcollection =
    "CP:REACTOME"
)

cat("\nDimensi Hallmark database:\n")
print(
  dim(
    hallmark_database
  )
)

cat("\nDimensi Reactome database:\n")
print(
  dim(
    reactome_database
  )
)

msigdb_version <- unique(
  c(
    as.character(
      hallmark_database$db_version
    ),
    as.character(
      reactome_database$db_version
    )
  )
)

cat("\nMSigDB version:\n")
print(
  msigdb_version
)

# ------------------------------------------------------------
# 13. Build TERM2GENE and TERM2NAME
# ------------------------------------------------------------

build_gene_set_tables <- function(
    gene_set_database,
    collection_label
) {
  
  term2gene <- gene_set_database %>%
    
    dplyr::transmute(
      
      term =
        as.character(
          gs_name
        ),
      
      gene =
        toupper(
          trimws(
            as.character(
              gene_symbol
            )
          )
        )
    ) %>%
    
    dplyr::filter(
      !is.na(term),
      term != "",
      !is.na(gene),
      gene != ""
    ) %>%
    
    dplyr::distinct()
  
  if (
    "gs_description" %in%
    colnames(
      gene_set_database
    )
  ) {
    
    term2name <- gene_set_database %>%
      
      dplyr::transmute(
        
        term =
          as.character(
            gs_name
          ),
        
        name =
          dplyr::coalesce(
            as.character(
              gs_description
            ),
            as.character(
              gs_name
            )
          )
      ) %>%
      
      dplyr::distinct(
        term,
        .keep_all = TRUE
      )
    
  } else {
    
    term2name <- gene_set_database %>%
      
      dplyr::transmute(
        
        term =
          as.character(
            gs_name
          ),
        
        name =
          as.character(
            gs_name
          )
      ) %>%
      
      dplyr::distinct()
  }
  
  catalog <- gene_set_database %>%
    
    dplyr::transmute(
      
      collection =
        collection_label,
      
      pathway_id =
        as.character(
          gs_name
        ),
      
      pathway_name =
        if (
          "gs_description" %in%
          colnames(
            gene_set_database
          )
        ) {
          
          dplyr::coalesce(
            as.character(
              gs_description
            ),
            as.character(
              gs_name
            )
          )
          
        } else {
          
          as.character(
            gs_name
          )
        }
    ) %>%
    
    dplyr::distinct(
      collection,
      pathway_id,
      .keep_all = TRUE
    )
  
  list(
    
    term2gene =
      term2gene,
    
    term2name =
      term2name,
    
    catalog =
      catalog
  )
}

hallmark_gene_sets <- build_gene_set_tables(
  hallmark_database,
  "Hallmark"
)

reactome_gene_sets <- build_gene_set_tables(
  reactome_database,
  "Reactome"
)

gene_set_catalog <- dplyr::bind_rows(
  
  hallmark_gene_sets$catalog,
  
  reactome_gene_sets$catalog
)

gene_set_summary <- tibble::tibble(
  
  collection = c(
    "Hallmark",
    "Reactome"
  ),
  
  pathways = c(
    dplyr::n_distinct(
      hallmark_gene_sets$
        term2gene$term
    ),
    dplyr::n_distinct(
      reactome_gene_sets$
        term2gene$term
    )
  ),
  
  unique_genes = c(
    dplyr::n_distinct(
      hallmark_gene_sets$
        term2gene$gene
    ),
    dplyr::n_distinct(
      reactome_gene_sets$
        term2gene$gene
    )
  )
)

cat("\nGene-set summary:\n")

print(
  gene_set_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 14. Function to run preranked GSEA
# ------------------------------------------------------------

run_preranked_gsea <- function(
    ranked_vector,
    dataset_id,
    collection_label,
    term2gene,
    term2name
) {
  
  message(
    "GSEA: ",
    dataset_id,
    " — ",
    collection_label
  )
  
  set.seed(
    20260722
  )
  
  gsea_object <- tryCatch(
    
    clusterProfiler::GSEA(
      
      geneList =
        ranked_vector,
      
      exponent =
        1,
      
      minGSSize =
        15,
      
      maxGSSize =
        500,
      
      eps =
        0,
      
      pvalueCutoff =
        1,
      
      pAdjustMethod =
        "BH",
      
      TERM2GENE =
        term2gene,
      
      TERM2NAME =
        term2name,
      
      verbose =
        FALSE,
      
      seed =
        TRUE,
      
      by =
        "fgsea"
    ),
    
    error = function(error_object) {
      
      stop(
        "GSEA gagal untuk ",
        dataset_id,
        " — ",
        collection_label,
        ":\n",
        conditionMessage(
          error_object
        )
      )
    }
  )
  
  result_table <- as.data.frame(
    gsea_object
  )
  
  if (
    nrow(
      result_table
    ) == 0
  ) {
    
    warning(
      "Tidak ada GSEA results untuk ",
      dataset_id,
      " — ",
      collection_label
    )
    
    return(
      list(
        object =
          gsea_object,
        
        table =
          tibble::tibble()
      )
    )
  }
  
  optional_columns <- c(
    "Description",
    "setSize",
    "enrichmentScore",
    "NES",
    "pvalue",
    "p.adjust",
    "qvalues",
    "rank",
    "leading_edge",
    "core_enrichment"
  )
  
  for (
    current_column in
    optional_columns
  ) {
    
    if (
      !current_column %in%
      colnames(
        result_table
      )
    ) {
      
      result_table[[current_column]] <-
        NA
    }
  }
  
  standardized_table <- result_table %>%
    
    tibble::rownames_to_column(
      "result_row"
    ) %>%
    
    tibble::as_tibble() %>%
    
    dplyr::transmute(
      
      dataset =
        dataset_id,
      
      collection =
        collection_label,
      
      pathway_id =
        as.character(
          ID
        ),
      
      pathway =
        dplyr::coalesce(
          as.character(
            Description
          ),
          as.character(
            ID
          )
        ),
      
      set_size =
        as.numeric(
          setSize
        ),
      
      enrichment_score =
        as.numeric(
          enrichmentScore
        ),
      
      NES =
        as.numeric(
          NES
        ),
      
      PValue =
        as.numeric(
          pvalue
        ),
      
      FDR =
        as.numeric(
          p.adjust
        ),
      
      qvalue =
        as.numeric(
          qvalues
        ),
      
      rank_at_maximum =
        as.numeric(
          rank
        ),
      
      leading_edge =
        as.character(
          leading_edge
        ),
      
      core_enrichment =
        as.character(
          core_enrichment
        )
    )
  
  list(
    
    object =
      gsea_object,
    
    table =
      standardized_table
  )
}

# ------------------------------------------------------------
# 15. Run GSEA for all datasets and both collections
# ------------------------------------------------------------

gsea_objects <- list()
gsea_result_tables <- list()

for (
  current_dataset in
  names(
    ranked_inputs
  )
) {
  
  current_vector <-
    ranked_inputs[[current_dataset]]$
    vector
  
  hallmark_result <- run_preranked_gsea(
    
    ranked_vector =
      current_vector,
    
    dataset_id =
      current_dataset,
    
    collection_label =
      "Hallmark",
    
    term2gene =
      hallmark_gene_sets$term2gene,
    
    term2name =
      hallmark_gene_sets$term2name
  )
  
  hallmark_key <- paste0(
    current_dataset,
    "_Hallmark"
  )
  
  gsea_objects[[hallmark_key]] <-
    hallmark_result$object
  
  gsea_result_tables[[hallmark_key]] <-
    hallmark_result$table
  
  reactome_result <- run_preranked_gsea(
    
    ranked_vector =
      current_vector,
    
    dataset_id =
      current_dataset,
    
    collection_label =
      "Reactome",
    
    term2gene =
      reactome_gene_sets$term2gene,
    
    term2name =
      reactome_gene_sets$term2name
  )
  
  reactome_key <- paste0(
    current_dataset,
    "_Reactome"
  )
  
  gsea_objects[[reactome_key]] <-
    reactome_result$object
  
  gsea_result_tables[[reactome_key]] <-
    reactome_result$table
}

all_gsea_results <- dplyr::bind_rows(
  gsea_result_tables
)

if (
  nrow(
    all_gsea_results
  ) == 0
) {
  
  stop(
    "Tidak ada hasil GSEA yang berhasil diperoleh."
  )
}

cat("\nJumlah GSEA results:\n")

print(
  all_gsea_results %>%
    
    dplyr::count(
      dataset,
      collection,
      name =
        "pathways_tested"
    ) %>%
    
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 16. Dataset specification
# ------------------------------------------------------------

dataset_specification <- tibble::tribble(
  
  ~dataset,
  ~dataset_label,
  ~evidence_group,
  ~dataset_order,
  
  "GSE114007",
  "GSE114007",
  "Independent cohorts",
  1,
  
  "GSE117999",
  "GSE117999",
  "Independent cohorts",
  2,
  
  "GSE169077",
  "GSE169077",
  "Independent cohorts",
  3,
  
  "GSE57218",
  "GSE57218\npaired",
  "Paired supportive cohort",
  4
)

primary_datasets <- c(
  "GSE114007",
  "GSE117999",
  "GSE169077"
)

# ------------------------------------------------------------
# 17. Cross-cohort pathway consensus summary
# ------------------------------------------------------------

pathway_consensus <- all_gsea_results %>%
  
  dplyr::filter(
    dataset %in%
      primary_datasets
  ) %>%
  
  dplyr::group_by(
    collection,
    pathway_id,
    pathway
  ) %>%
  
  dplyr::summarise(
    
    primary_cohorts_tested =
      dplyr::n_distinct(
        dataset
      ),
    
    primary_significant_FDR05 =
      sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
    
    primary_positive_NES =
      sum(
        NES > 0,
        na.rm = TRUE
      ),
    
    primary_negative_NES =
      sum(
        NES < 0,
        na.rm = TRUE
      ),
    
    mean_NES =
      mean(
        NES,
        na.rm = TRUE
      ),
    
    median_NES =
      stats::median(
        NES,
        na.rm = TRUE
      ),
    
    mean_absolute_NES =
      mean(
        abs(
          NES
        ),
        na.rm = TRUE
      ),
    
    minimum_FDR =
      min(
        FDR,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::mutate(
    
    direction_consistent =
      (
        primary_positive_NES ==
          primary_cohorts_tested
      ) |
      (
        primary_negative_NES ==
          primary_cohorts_tested
      ),
    
    consensus_direction =
      dplyr::case_when(
        
        primary_positive_NES ==
          primary_cohorts_tested ~
          "Positive in all tested cohorts",
        
        primary_negative_NES ==
          primary_cohorts_tested ~
          "Negative in all tested cohorts",
        
        TRUE ~
          "Mixed direction"
      )
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      primary_significant_FDR05
    ),
    dplyr::desc(
      direction_consistent
    ),
    dplyr::desc(
      mean_absolute_NES
    )
  )

cat("\nTop cross-cohort pathway consensus results:\n")

print(
  pathway_consensus %>%
    dplyr::slice_head(
      n = 20
    ) %>%
    tibble::as_tibble(),
  n = 20,
  width = Inf
)

# ------------------------------------------------------------
# 18. Prespecified OA/MMP13-related pathways
#
# These pathways are specified before inspecting Figure 1C.
# ------------------------------------------------------------

target_pathways <- tibble::tribble(
  
  ~target_order,
  ~theme,
  ~collection,
  ~preferred_id,
  ~fallback_pattern,
  ~display_name,
  
  1,
  "Inflammatory signaling",
  "Hallmark",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "TNFA.*SIGNALING.*NFKB",
  "TNFα–NF-κB signaling",
  
  2,
  "Inflammatory signaling",
  "Hallmark",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "INFLAMMATORY_RESPONSE",
  "Inflammatory response",
  
  3,
  "Inflammatory signaling",
  "Hallmark",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "IL6.*JAK.*STAT3",
  "IL-6–JAK–STAT3 signaling",
  
  4,
  "Inflammatory signaling",
  "Reactome",
  "REACTOME_SIGNALING_BY_INTERLEUKINS",
  "SIGNALING_BY_INTERLEUKINS",
  "Interleukin signaling",
  
  5,
  "MAPK signaling",
  "Reactome",
  "REACTOME_P38MAPK_EVENTS",
  "P38.*MAPK.*EVENT",
  "p38 MAPK events",
  
  6,
  "MAPK signaling",
  "Reactome",
  "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES",
  "MAPK_FAMILY_SIGNALING_CASCADES",
  "MAPK-family cascades",
  
  7,
  "Extracellular matrix",
  "Reactome",
  "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "EXTRACELLULAR_MATRIX_ORGANIZATION",
  "Extracellular-matrix organization",
  
  8,
  "Extracellular matrix",
  "Reactome",
  "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
  "Extracellular-matrix degradation",
  
  9,
  "Extracellular matrix",
  "Reactome",
  "REACTOME_COLLAGEN_DEGRADATION",
  "COLLAGEN_DEGRADATION",
  "Collagen degradation",
  
  10,
  "Extracellular matrix",
  "Reactome",
  "REACTOME_ACTIVATION_OF_MATRIX_METALLOPROTEINASES",
  "ACTIVATION_OF_MATRIX_METALLOPROTEINASES",
  "Matrix-metalloproteinase activation",
  
  11,
  "Cartilage remodeling",
  "Hallmark",
  "HALLMARK_TGF_BETA_SIGNALING",
  "TGF_BETA_SIGNALING",
  "TGF-β signaling",
  
  12,
  "Cellular stress",
  "Hallmark",
  "HALLMARK_HYPOXIA",
  "HALLMARK_HYPOXIA",
  "Hypoxia",
  
  13,
  "Cellular stress",
  "Hallmark",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_APOPTOSIS",
  "Apoptosis"
)

# ------------------------------------------------------------
# 19. Resolve target pathways against current MSigDB version
# ------------------------------------------------------------

resolve_one_target <- function(
    target_row,
    catalog
) {
  
  collection_value <-
    as.character(
      target_row$collection
    )
  
  preferred_value <-
    as.character(
      target_row$preferred_id
    )
  
  pattern_value <-
    as.character(
      target_row$fallback_pattern
    )
  
  available_catalog <- catalog %>%
    
    dplyr::filter(
      collection ==
        collection_value
    )
  
  exact_match <- available_catalog %>%
    
    dplyr::filter(
      pathway_id ==
        preferred_value
    )
  
  if (
    nrow(
      exact_match
    ) > 0
  ) {
    
    selected_pathway <-
      exact_match %>%
      dplyr::slice_head(
        n = 1
      )
    
    match_method <-
      "Exact preferred ID"
    
  } else {
    
    pattern_match <- available_catalog %>%
      
      dplyr::filter(
        stringr::str_detect(
          stringr::str_to_upper(
            pathway_id
          ),
          stringr::regex(
            pattern_value,
            ignore_case = TRUE
          )
        )
      ) %>%
      
      dplyr::arrange(
        nchar(
          pathway_id
        )
      )
    
    if (
      nrow(
        pattern_match
      ) > 0
    ) {
      
      selected_pathway <-
        pattern_match %>%
        dplyr::slice_head(
          n = 1
        )
      
      match_method <-
        "Fallback regular expression"
      
    } else {
      
      return(
        target_row %>%
          dplyr::mutate(
            
            pathway_id =
              NA_character_,
            
            pathway_name =
              NA_character_,
            
            match_method =
              "Not found"
          )
      )
    }
  }
  
  target_row %>%
    
    dplyr::mutate(
      
      pathway_id =
        selected_pathway$
        pathway_id[1],
      
      pathway_name =
        selected_pathway$
        pathway_name[1],
      
      match_method =
        match_method
    )
}

resolved_target_list <- lapply(
  seq_len(
    nrow(
      target_pathways
    )
  ),
  function(current_index) {
    
    resolve_one_target(
      
      target_row =
        target_pathways[
          current_index,
          ,
          drop = FALSE
        ],
      
      catalog =
        gene_set_catalog
    )
  }
)

resolved_target_pathways <- dplyr::bind_rows(
  resolved_target_list
)

cat("\nResolved target pathways:\n")

print(
  resolved_target_pathways %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

missing_target_pathways <- resolved_target_pathways %>%
  
  dplyr::filter(
    is.na(
      pathway_id
    )
  )

if (
  nrow(
    missing_target_pathways
  ) > 0
) {
  
  warning(
    "Beberapa target pathways tidak ditemukan pada MSigDB version ini."
  )
  
  print(
    missing_target_pathways %>%
      dplyr::select(
        collection,
        preferred_id,
        display_name
      ) %>%
      tibble::as_tibble(),
    n = Inf,
    width = Inf
  )
}

resolved_target_pathways <- resolved_target_pathways %>%
  
  dplyr::filter(
    !is.na(
      pathway_id
    )
  ) %>%
  
  dplyr::distinct(
    pathway_id,
    .keep_all = TRUE
  )

if (
  nrow(
    resolved_target_pathways
  ) < 8
) {
  
  stop(
    "Terlalu sedikit target pathways yang berhasil ditemukan: ",
    nrow(
      resolved_target_pathways
    )
  )
}

# ------------------------------------------------------------
# 20. Figure 1C data
# ------------------------------------------------------------

figure1c_data <- tidyr::crossing(
  
  resolved_target_pathways %>%
    
    dplyr::select(
      target_order,
      theme,
      collection,
      pathway_id,
      display_name
    ),
  
  dataset_specification %>%
    
    dplyr::select(
      dataset,
      dataset_label,
      evidence_group,
      dataset_order
    )
) %>%
  
  dplyr::left_join(
    
    all_gsea_results %>%
      
      dplyr::select(
        dataset,
        collection,
        pathway_id,
        NES,
        PValue,
        FDR,
        set_size,
        core_enrichment
      ),
    
    by = c(
      "dataset",
      "collection",
      "pathway_id"
    )
  ) %>%
  
  dplyr::mutate(
    
    significance_symbol =
      dplyr::case_when(
        
        is.na(FDR) ~
          "",
        
        FDR < 0.001 ~
          "***",
        
        FDR < 0.01 ~
          "**",
        
        FDR < 0.05 ~
          "*",
        
        TRUE ~
          ""
      ),
    
    tile_label =
      dplyr::case_when(
        
        is.na(NES) ~
          "NA",
        
        TRUE ~
          paste0(
            sprintf(
              "%.2f",
              NES
            ),
            significance_symbol
          )
      )
  )

figure1c_pathway_order <- resolved_target_pathways %>%
  
  dplyr::arrange(
    target_order
  ) %>%
  
  dplyr::pull(
    display_name
  )

figure1c_dataset_order <- dataset_specification %>%
  
  dplyr::arrange(
    dataset_order
  ) %>%
  
  dplyr::pull(
    dataset_label
  )

figure1c_data <- figure1c_data %>%
  
  dplyr::mutate(
    
    display_name = factor(
      display_name,
      levels = rev(
        figure1c_pathway_order
      )
    ),
    
    dataset_label = factor(
      dataset_label,
      levels =
        figure1c_dataset_order
    ),
    
    evidence_group = factor(
      evidence_group,
      levels = c(
        "Independent cohorts",
        "Paired supportive cohort"
      )
    )
  )

maximum_absolute_NES <- max(
  abs(
    figure1c_data$NES
  ),
  na.rm = TRUE
)

if (
  !is.finite(
    maximum_absolute_NES
  ) ||
  maximum_absolute_NES < 1.5
) {
  
  maximum_absolute_NES <-
    1.5
}

figure1c_data <- figure1c_data %>%
  
  dplyr::mutate(
    
    text_color =
      dplyr::if_else(
        
        !is.na(NES) &
          abs(NES) >=
          0.62 *
          maximum_absolute_NES,
        
        "white",
        
        "black"
      )
  )

cat("\nFigure 1C data:\n")

print(
  figure1c_data %>%
    
    dplyr::select(
      theme,
      display_name,
      dataset,
      NES,
      FDR,
      tile_label
    ) %>%
    
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 21. Figure 1C pathway heatmap
# ------------------------------------------------------------

figure1c_plot <- ggplot(
  
  figure1c_data,
  
  aes(
    x = dataset_label,
    y = display_name,
    fill = NES
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.7
  ) +
  
  geom_text(
    
    aes(
      label = tile_label,
      color = text_color
    ),
    
    size = 3.1,
    fontface = "bold"
  ) +
  
  facet_grid(
    
    cols =
      vars(
        evidence_group
      ),
    
    scales =
      "free_x",
    
    space =
      "free_x"
  ) +
  
  scale_fill_gradient2(
    
    low =
      "#2166AC",
    
    mid =
      "white",
    
    high =
      "#B2182B",
    
    midpoint =
      0,
    
    limits = c(
      -maximum_absolute_NES,
      maximum_absolute_NES
    ),
    
    na.value =
      "grey90",
    
    name =
      "Normalized\nenrichment\nscore"
  ) +
  
  scale_color_identity() +
  
  labs(
    
    title =
      paste(
        "Figure 1C | OA-related pathway activity",
        "across human cartilage cohorts"
      ),
    
    subtitle =
      paste0(
        "Preranked GSEA using MSigDB Hallmark and Reactome; ",
        "positive NES indicates enrichment toward OA. ",
        "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001."
      ),
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_bw(
    base_size = 10.5
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 9
      ),
    
    axis.text.y =
      element_text(
        size = 9
      ),
    
    strip.background =
      element_rect(
        fill = "grey95",
        color = "grey40"
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 9
      ),
    
    plot.title =
      element_text(
        face = "bold",
        size = 12
      ),
    
    plot.subtitle =
      element_text(
        size = 9
      ),
    
    legend.position =
      "right",
    
    plot.margin =
      margin(
        10,
        12,
        10,
        10
      )
  )

print(
  figure1c_plot
)

ggsave(
  
  filename = file.path(
    main_figure_folder,
    "Figure1C_cross_cohort_pathway_GSEA.pdf"
  ),
  
  plot =
    figure1c_plot,
  
  width =
    10.5,
  
  height =
    7.6
)

ggsave(
  
  filename = file.path(
    main_figure_folder,
    "Figure1C_cross_cohort_pathway_GSEA.tiff"
  ),
  
  plot =
    figure1c_plot,
  
  width =
    10.5,
  
  height =
    7.6,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 22. Supplementary Figure S10:
#     Complete Hallmark GSEA heatmap
# ------------------------------------------------------------

hallmark_pathway_order_table <- all_gsea_results %>%
  
  dplyr::filter(
    collection ==
      "Hallmark",
    dataset %in%
      primary_datasets
  ) %>%
  
  dplyr::group_by(
    pathway_id
  ) %>%
  
  dplyr::summarise(
    
    mean_NES =
      mean(
        NES,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      mean_NES
    )
  )

hallmark_pathways <- hallmark_pathway_order_table$
  pathway_id

hallmark_supplementary_data <- tidyr::crossing(
  
  pathway_id =
    hallmark_pathways,
  
  dataset_specification %>%
    
    dplyr::select(
      dataset,
      dataset_label,
      evidence_group,
      dataset_order
    )
) %>%
  
  dplyr::left_join(
    
    all_gsea_results %>%
      
      dplyr::filter(
        collection ==
          "Hallmark"
      ) %>%
      
      dplyr::select(
        dataset,
        pathway_id,
        NES,
        FDR
      ),
    
    by = c(
      "dataset",
      "pathway_id"
    )
  ) %>%
  
  dplyr::mutate(
    
    pathway_display =
      pathway_id %>%
      stringr::str_remove(
        "^HALLMARK_"
      ) %>%
      stringr::str_replace_all(
        "_",
        " "
      ) %>%
      stringr::str_to_title(),
    
    significance_symbol =
      dplyr::case_when(
        
        is.na(FDR) ~
          "",
        
        FDR < 0.001 ~
          "***",
        
        FDR < 0.01 ~
          "**",
        
        FDR < 0.05 ~
          "*",
        
        TRUE ~
          ""
      )
  )

hallmark_display_order <- hallmark_pathway_order_table$
  pathway_id %>%
  stringr::str_remove(
    "^HALLMARK_"
  ) %>%
  stringr::str_replace_all(
    "_",
    " "
  ) %>%
  stringr::str_to_title()

hallmark_supplementary_data <-
  hallmark_supplementary_data %>%
  
  dplyr::mutate(
    
    pathway_display = factor(
      pathway_display,
      levels = rev(
        hallmark_display_order
      )
    ),
    
    dataset_label = factor(
      dataset_label,
      levels =
        figure1c_dataset_order
    ),
    
    evidence_group = factor(
      evidence_group,
      levels = c(
        "Independent cohorts",
        "Paired supportive cohort"
      )
    )
  )

hallmark_maximum_NES <- max(
  abs(
    hallmark_supplementary_data$NES
  ),
  na.rm = TRUE
)

hallmark_plot <- ggplot(
  
  hallmark_supplementary_data,
  
  aes(
    x = dataset_label,
    y = pathway_display,
    fill = NES
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  
  geom_text(
    aes(
      label =
        significance_symbol
    ),
    size = 2.2,
    fontface = "bold"
  ) +
  
  facet_grid(
    
    cols =
      vars(
        evidence_group
      ),
    
    scales =
      "free_x",
    
    space =
      "free_x"
  ) +
  
  scale_fill_gradient2(
    
    low =
      "#2166AC",
    
    mid =
      "white",
    
    high =
      "#B2182B",
    
    midpoint =
      0,
    
    limits = c(
      -hallmark_maximum_NES,
      hallmark_maximum_NES
    ),
    
    na.value =
      "grey90",
    
    name =
      "NES"
  ) +
  
  labs(
    
    title =
      "Supplementary Figure S10 | Hallmark GSEA across cartilage cohorts",
    
    subtitle =
      "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001",
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_bw(
    base_size = 9
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 8
      ),
    
    axis.text.y =
      element_text(
        size = 6.8
      ),
    
    strip.background =
      element_rect(
        fill = "grey95"
      ),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      )
  )

print(
  hallmark_plot
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS10_Hallmark_GSEA_heatmap.pdf"
  ),
  
  plot =
    hallmark_plot,
  
  width =
    10.5,
  
  height =
    13
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS10_Hallmark_GSEA_heatmap.tiff"
  ),
  
  plot =
    hallmark_plot,
  
  width =
    10.5,
  
  height =
    13,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 23. Supplementary Figure S11:
#     Top cross-cohort Reactome pathways
# ------------------------------------------------------------

top_reactome_pathways <- pathway_consensus %>%
  
  dplyr::filter(
    collection ==
      "Reactome"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      primary_significant_FDR05
    ),
    dplyr::desc(
      direction_consistent
    ),
    dplyr::desc(
      mean_absolute_NES
    )
  ) %>%
  
  dplyr::slice_head(
    n = 24
  )

reactome_supplementary_data <- tidyr::crossing(
  
  top_reactome_pathways %>%
    
    dplyr::select(
      pathway_id,
      pathway
    ),
  
  dataset_specification %>%
    
    dplyr::select(
      dataset,
      dataset_label,
      evidence_group,
      dataset_order
    )
) %>%
  
  dplyr::left_join(
    
    all_gsea_results %>%
      
      dplyr::filter(
        collection ==
          "Reactome"
      ) %>%
      
      dplyr::select(
        dataset,
        pathway_id,
        NES,
        FDR
      ),
    
    by = c(
      "dataset",
      "pathway_id"
    )
  ) %>%
  
  dplyr::mutate(
    
    pathway_display =
      pathway_id %>%
      stringr::str_remove(
        "^REACTOME_"
      ) %>%
      stringr::str_replace_all(
        "_",
        " "
      ) %>%
      stringr::str_to_title(),
    
    significance_symbol =
      dplyr::case_when(
        
        is.na(FDR) ~
          "",
        
        FDR < 0.001 ~
          "***",
        
        FDR < 0.01 ~
          "**",
        
        FDR < 0.05 ~
          "*",
        
        TRUE ~
          ""
      ),
    
    tile_label =
      dplyr::case_when(
        
        is.na(NES) ~
          "NA",
        
        TRUE ~
          paste0(
            sprintf(
              "%.1f",
              NES
            ),
            significance_symbol
          )
      )
  )

reactome_display_order <- top_reactome_pathways$
  pathway_id %>%
  stringr::str_remove(
    "^REACTOME_"
  ) %>%
  stringr::str_replace_all(
    "_",
    " "
  ) %>%
  stringr::str_to_title()

reactome_supplementary_data <-
  reactome_supplementary_data %>%
  
  dplyr::mutate(
    
    pathway_display = factor(
      pathway_display,
      levels = rev(
        reactome_display_order
      )
    ),
    
    dataset_label = factor(
      dataset_label,
      levels =
        figure1c_dataset_order
    ),
    
    evidence_group = factor(
      evidence_group,
      levels = c(
        "Independent cohorts",
        "Paired supportive cohort"
      )
    )
  )

reactome_maximum_NES <- max(
  abs(
    reactome_supplementary_data$NES
  ),
  na.rm = TRUE
)

reactome_supplementary_data <-
  reactome_supplementary_data %>%
  
  dplyr::mutate(
    
    text_color =
      dplyr::if_else(
        
        !is.na(NES) &
          abs(NES) >=
          0.62 *
          reactome_maximum_NES,
        
        "white",
        
        "black"
      )
  )

reactome_plot <- ggplot(
  
  reactome_supplementary_data,
  
  aes(
    x = dataset_label,
    y = pathway_display,
    fill = NES
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  
  geom_text(
    
    aes(
      label = tile_label,
      color = text_color
    ),
    
    size = 2.6,
    fontface = "bold"
  ) +
  
  facet_grid(
    
    cols =
      vars(
        evidence_group
      ),
    
    scales =
      "free_x",
    
    space =
      "free_x"
  ) +
  
  scale_fill_gradient2(
    
    low =
      "#2166AC",
    
    mid =
      "white",
    
    high =
      "#B2182B",
    
    midpoint =
      0,
    
    limits = c(
      -reactome_maximum_NES,
      reactome_maximum_NES
    ),
    
    na.value =
      "grey90",
    
    name =
      "NES"
  ) +
  
  scale_color_identity() +
  
  labs(
    
    title =
      paste(
        "Supplementary Figure S11 | Top cross-cohort",
        "Reactome pathways"
      ),
    
    subtitle =
      paste0(
        "Selected by significance count, directional consistency, ",
        "and mean absolute NES across independent cohorts."
      ),
    
    x =
      NULL,
    
    y =
      NULL
  ) +
  
  theme_bw(
    base_size = 9
  ) +
  
  theme(
    
    panel.grid =
      element_blank(),
    
    axis.text.x =
      element_text(
        face = "bold",
        size = 8
      ),
    
    axis.text.y =
      element_text(
        size = 7.2
      ),
    
    strip.background =
      element_rect(
        fill = "grey95"
      ),
    
    strip.text =
      element_text(
        face = "bold"
      ),
    
    plot.title =
      element_text(
        face = "bold"
      )
  )

print(
  reactome_plot
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS11_Reactome_consensus_heatmap.pdf"
  ),
  
  plot =
    reactome_plot,
  
  width =
    10.5,
  
  height =
    9
)

ggsave(
  
  filename = file.path(
    supplementary_figure_folder,
    "FigureS11_Reactome_consensus_heatmap.tiff"
  ),
  
  plot =
    reactome_plot,
  
  width =
    10.5,
  
  height =
    9,
  
  dpi =
    600,
  
  compression =
    "lzw"
)

# ------------------------------------------------------------
# 24. Significant pathway summaries
# ------------------------------------------------------------

significant_gsea_results <- all_gsea_results %>%
  
  dplyr::filter(
    FDR < 0.05
  ) %>%
  
  dplyr::arrange(
    dataset,
    collection,
    FDR,
    dplyr::desc(
      abs(
        NES
      )
    )
  )

gsea_dataset_summary <- all_gsea_results %>%
  
  dplyr::group_by(
    dataset,
    collection
  ) %>%
  
  dplyr::summarise(
    
    pathways_tested =
      dplyr::n(),
    
    significant_FDR05 =
      sum(
        FDR < 0.05,
        na.rm = TRUE
      ),
    
    positive_significant =
      sum(
        FDR < 0.05 &
          NES > 0,
        na.rm = TRUE
      ),
    
    negative_significant =
      sum(
        FDR < 0.05 &
          NES < 0,
        na.rm = TRUE
      ),
    
    maximum_NES =
      max(
        NES,
        na.rm = TRUE
      ),
    
    minimum_NES =
      min(
        NES,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )

cat("\nGSEA dataset summary:\n")

print(
  gsea_dataset_summary %>%
    tibble::as_tibble(),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 25. Save R objects
# ------------------------------------------------------------

saveRDS(
  
  ranked_inputs,
  
  file = file.path(
    processed_folder,
    "cross_cohort_standardized_ranked_inputs.rds"
  )
)

saveRDS(
  
  gsea_objects,
  
  file = file.path(
    processed_folder,
    "cross_cohort_GSEA_objects.rds"
  )
)

saveRDS(
  
  all_gsea_results,
  
  file = file.path(
    processed_folder,
    "cross_cohort_GSEA_results.rds"
  )
)

saveRDS(
  
  pathway_consensus,
  
  file = file.path(
    processed_folder,
    "cross_cohort_pathway_consensus.rds"
  )
)

saveRDS(
  
  figure1c_data,
  
  file = file.path(
    processed_folder,
    "Figure1C_pathway_GSEA_data.rds"
  )
)

saveRDS(
  
  resolved_target_pathways,
  
  file = file.path(
    processed_folder,
    "Figure1C_resolved_target_pathways.rds"
  )
)

# ------------------------------------------------------------
# 26. Save CSV tables
# ------------------------------------------------------------

utils::write.csv(
  
  all_gsea_results,
  
  file = file.path(
    table_folder,
    "cross_cohort_all_GSEA_results.csv"
  ),
  
  row.names = FALSE
)

utils::write.csv(
  
  pathway_consensus,
  
  file = file.path(
    table_folder,
    "cross_cohort_pathway_consensus.csv"
  ),
  
  row.names = FALSE
)

utils::write.csv(
  
  figure1c_data,
  
  file = file.path(
    table_folder,
    "Figure1C_pathway_GSEA_data.csv"
  ),
  
  row.names = FALSE
)

# ------------------------------------------------------------
# 27. Save Excel workbook
# ------------------------------------------------------------

openxlsx::write.xlsx(
  
  list(
    
    Input_Summary =
      as.data.frame(
        ranked_input_summary
      ),
    
    GSEA_Summary =
      as.data.frame(
        gsea_dataset_summary
      ),
    
    Target_Resolution =
      as.data.frame(
        resolved_target_pathways
      ),
    
    Figure1C_Data =
      as.data.frame(
        figure1c_data
      ),
    
    Pathway_Consensus =
      as.data.frame(
        pathway_consensus
      ),
    
    Significant_GSEA =
      as.data.frame(
        significant_gsea_results
      ),
    
    Hallmark_Results =
      as.data.frame(
        all_gsea_results %>%
          dplyr::filter(
            collection ==
              "Hallmark"
          )
      ),
    
    Reactome_Results =
      as.data.frame(
        all_gsea_results %>%
          dplyr::filter(
            collection ==
              "Reactome"
          )
      ),
    
    MSigDB_Collections =
      as.data.frame(
        msigdb_collections
      ),
    
    Gene_Set_Summary =
      as.data.frame(
        gene_set_summary
      )
  ),
  
  file = file.path(
    table_folder,
    "cross_cohort_pathway_enrichment_results.xlsx"
  ),
  
  overwrite = TRUE
)

# ------------------------------------------------------------
# 28. Save session information
# ------------------------------------------------------------

sink(
  file.path(
    table_folder,
    "cross_cohort_pathway_enrichment_sessionInfo.txt"
  )
)

print(
  sessionInfo()
)

sink()

# ------------------------------------------------------------
# 29. Final messages
# ------------------------------------------------------------

message("")
message("================================================")
message("CROSS-COHORT PATHWAY ENRICHMENT SELESAI")
message("================================================")

message(
  "MSigDB version                  : ",
  paste(
    msigdb_version,
    collapse = ", "
  )
)

message(
  "Datasets analyzed              : ",
  length(
    ranked_inputs
  )
)

message(
  "Total GSEA result rows         : ",
  nrow(
    all_gsea_results
  )
)

message(
  "Resolved Figure 1C pathways    : ",
  nrow(
    resolved_target_pathways
  )
)

message(
  "Significant GSEA results       : ",
  nrow(
    significant_gsea_results
  )
)

message(
  "Figure 1C saved                : ",
  file.path(
    main_figure_folder,
    "Figure1C_cross_cohort_pathway_GSEA"
  )
)

message(
  "Supplementary S10 saved        : ",
  file.path(
    supplementary_figure_folder,
    "FigureS10_Hallmark_GSEA_heatmap"
  )
)

message(
  "Supplementary S11 saved        : ",
  file.path(
    supplementary_figure_folder,
    "FigureS11_Reactome_consensus_heatmap"
  )
)

message("================================================")