# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 19A_GSE220243_NFkB_p38_MMP13_exploratory.R
# PURPOSE :
#   Biology-first exploratory analysis for Figure 1E:
#
#   1. Retrieve curated human pathway gene sets from MSigDB.
#   2. Score NF-kB, p38 MAPK, inflammatory-response, and
#      extracellular-matrix-degradation programs.
#   3. Aggregate scores by donor and annotated cell state.
#   4. Quantify MMP13 from raw RNA counts at donor-state level.
#   5. Explore pathway-MMP13 relationships without treating
#      individual cells as independent biological replicates.
#
# IMPORTANT:
#   - Donor is the biological replicate.
#   - Integrated coordinates are NOT used for expression scoring.
#   - Module scores are calculated from log-normalized RNA data.
#   - MMP13 is summarized from raw RNA counts.
#   - Results are exploratory; Figure 1E is not finalized here.
# ============================================================

rm(list = ls())
gc()

set.seed(20260723)

# ------------------------------------------------------------
# 1. Project root
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

cat("\nWorking directory:\n")
print(getwd())

# ------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "msigdbr",
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "patchwork",
  "openxlsx",
  "scales"
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
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(patchwork)
library(openxlsx)
library(scales)

package_versions <- tibble::tibble(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package_name) {
      as.character(
        packageVersion(package_name)
      )
    },
    FUN.VALUE = character(1)
  )
)

cat("\nPackage versions:\n")
print(
  package_versions,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 3. Input and output paths
# ------------------------------------------------------------

input_file <- file.path(
  "data_processed",
  "GSE220243_cartilage_annotated_cell_states.rds"
)

processed_folder <- "data_processed"

table_folder <- file.path(
  "results",
  "tables"
)

figure_folder <- file.path(
  "results",
  "figures",
  "Figure 1 Suplementary",
  "S17_NFkB_p38_exploration"
)

dir.create(
  processed_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Annotated Seurat object tidak ditemukan:\n",
    input_file
  )
}

# ------------------------------------------------------------
# 4. Load and validate annotated object
# ------------------------------------------------------------

message("Memuat annotated Seurat object...")

cartilage_object <- readRDS(input_file)

if (!inherits(cartilage_object, "Seurat")) {
  stop("Input bukan objek Seurat.")
}

required_metadata <- c(
  "sample_label",
  "group",
  "donor",
  "geo_accession",
  "analysis_cluster",
  "cell_state",
  "cell_state_short"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(cartilage_object[[]])
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

DefaultAssay(cartilage_object) <- "RNA"

if (
  inherits(cartilage_object[["RNA"]], "Assay5") &&
  length(
    SeuratObject::Layers(
      cartilage_object[["RNA"]],
      search = "^counts"
    )
  ) > 1
) {
  message("Joining RNA layers...")
  cartilage_object <- SeuratObject::JoinLayers(
    cartilage_object,
    assay = "RNA"
  )
}

rna_layers <- SeuratObject::Layers(
  cartilage_object[["RNA"]]
)

if (!"counts" %in% rna_layers) {
  stop("RNA counts layer tidak ditemukan.")
}

if (!"data" %in% rna_layers) {
  message(
    "RNA data layer tidak ditemukan; menjalankan NormalizeData..."
  )
  
  cartilage_object <- NormalizeData(
    cartilage_object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
}

state_order <- c(
  "RegC",
  "EC",
  "RepC",
  "PreFC",
  "FC",
  "MTC-like",
  "HomC",
  "HTC",
  "Vascular/immune"
)

cartilage_object$cell_state_short <- factor(
  as.character(cartilage_object$cell_state_short),
  levels = state_order
)

cartilage_object$group <- factor(
  as.character(cartilage_object$group),
  levels = c("Normal", "OA")
)

cat("\nAnnotated object:\n")
print(cartilage_object)

# ------------------------------------------------------------
# 5. Retrieve curated MSigDB gene sets
#
# Selected pathways:
#   - HALLMARK_TNFA_SIGNALING_VIA_NFKB
#   - HALLMARK_INFLAMMATORY_RESPONSE
#   - REACTOME_P38MAPK_EVENTS
#   - REACTOME_ACTIVATED_TAK1_MEDIATES_P38_MAPK_ACTIVATION
#   - REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX
# ------------------------------------------------------------

message("Retrieving MSigDB 2026 human gene sets...")

hallmark_df <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "H"
)

reactome_df <- msigdbr::msigdbr(
  db_species = "HS",
  species = "Homo sapiens",
  collection = "C2",
  subcollection = "CP:REACTOME"
)

msigdb_version <- unique(
  c(
    hallmark_df$db_version,
    reactome_df$db_version
  )
)

target_gene_sets <- list(
  
  NFkB = hallmark_df %>%
    dplyr::filter(
      gs_name ==
        "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
    ) %>%
    dplyr::pull(gene_symbol) %>%
    unique(),
  
  Inflammatory = hallmark_df %>%
    dplyr::filter(
      gs_name ==
        "HALLMARK_INFLAMMATORY_RESPONSE"
    ) %>%
    dplyr::pull(gene_symbol) %>%
    unique(),
  
  p38_MAPK = unique(
    c(
      reactome_df %>%
        dplyr::filter(
          gs_name ==
            "REACTOME_P38MAPK_EVENTS"
        ) %>%
        dplyr::pull(gene_symbol),
      
      reactome_df %>%
        dplyr::filter(
          gs_name ==
            paste0(
              "REACTOME_ACTIVATED_TAK1_",
              "MEDIATES_P38_MAPK_ACTIVATION"
            )
        ) %>%
        dplyr::pull(gene_symbol)
    )
  ),
  
  ECM_degradation = reactome_df %>%
    dplyr::filter(
      gs_name ==
        paste0(
          "REACTOME_DEGRADATION_OF_",
          "THE_EXTRACELLULAR_MATRIX"
        )
    ) %>%
    dplyr::pull(gene_symbol) %>%
    unique()
)

empty_sets <- names(
  target_gene_sets
)[
  lengths(target_gene_sets) == 0
]

if (length(empty_sets) > 0) {
  stop(
    "Gene set berikut tidak ditemukan di MSigDB: ",
    paste(
      empty_sets,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 6. Match gene sets to the RNA assay
# ------------------------------------------------------------

available_features <- rownames(cartilage_object)

matched_gene_sets <- lapply(
  target_gene_sets,
  function(current_set) {
    intersect(
      unique(current_set),
      available_features
    )
  }
)

gene_set_qc <- tibble::tibble(
  pathway = names(target_gene_sets),
  requested_genes = lengths(target_gene_sets),
  matched_genes = lengths(matched_gene_sets),
  matched_percent =
    100 *
    matched_genes /
    requested_genes,
  matched_gene_symbols = vapply(
    matched_gene_sets,
    function(x) {
      paste(
        x,
        collapse = "; "
      )
    },
    FUN.VALUE = character(1)
  )
)

cat("\nGene-set matching summary:\n")
print(
  gene_set_qc,
  n = Inf,
  width = Inf
)

insufficient_sets <- gene_set_qc$pathway[
  gene_set_qc$matched_genes < 5
]

if (length(insufficient_sets) > 0) {
  stop(
    "Terlalu sedikit genes cocok untuk pathway: ",
    paste(
      insufficient_sets,
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# 7. Score pathway programs from log-normalized RNA expression
#
# AddModuleScore subtracts matched control-gene expression.
# A fixed seed is used for reproducibility.
# ------------------------------------------------------------

message("Calculating pathway module scores...")

score_column_map <- c()

for (current_pathway in names(matched_gene_sets)) {
  
  score_prefix <- paste0(
    current_pathway,
    "_Score"
  )
  
  cartilage_object <- Seurat::AddModuleScore(
    object = cartilage_object,
    features = list(
      matched_gene_sets[[current_pathway]]
    ),
    assay = "RNA",
    name = score_prefix,
    nbin = 24,
    ctrl = 50,
    seed = 20260723,
    search = FALSE,
    slot = "data"
  )
  
  generated_column <- paste0(
    score_prefix,
    "1"
  )
  
  final_column <- paste0(
    current_pathway,
    "_module_score"
  )
  
  cartilage_object[[final_column]] <-
    cartilage_object[[generated_column]][, 1]
  
  cartilage_object[[generated_column]] <- NULL
  
  score_column_map[
    current_pathway
  ] <- final_column
}

score_columns <- unname(
  score_column_map
)

# ------------------------------------------------------------
# 8. Quantify raw MMP13 and total RNA library size per cell
# ------------------------------------------------------------

mmp13_candidates <- available_features[
  toupper(available_features) == "MMP13"
]

if (length(mmp13_candidates) == 0) {
  stop("MMP13 tidak ditemukan pada RNA assay.")
}

mmp13_feature <- mmp13_candidates[1]

rna_counts <- SeuratObject::LayerData(
  cartilage_object,
  assay = "RNA",
  layer = "counts"
)

mmp13_raw_count <- as.numeric(
  rna_counts[
    mmp13_feature,
    ,
    drop = TRUE
  ]
)

cell_total_umi <- as.numeric(
  Matrix::colSums(rna_counts)
)

cartilage_object$MMP13_raw_count <-
  mmp13_raw_count

cartilage_object$MMP13_detected <-
  mmp13_raw_count > 0

cartilage_object$total_RNA_UMI <-
  cell_total_umi

# ------------------------------------------------------------
# 9. Cell-level score table
# ------------------------------------------------------------

cell_score_table <- cartilage_object[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  dplyr::select(
    cell_id,
    sample_label,
    donor,
    geo_accession,
    group,
    analysis_cluster,
    cell_state_short,
    cell_state,
    dplyr::all_of(score_columns),
    MMP13_raw_count,
    MMP13_detected,
    total_RNA_UMI
  )

saveRDS(
  cell_score_table,
  file = file.path(
    processed_folder,
    "GSE220243_NFkB_p38_cell_scores.rds"
  ),
  compress = "gzip"
)

# ------------------------------------------------------------
# 10. Donor-state aggregation
#
# Each row below represents one donor within one annotated state.
# This is the primary table for Figure 1E exploration.
# ------------------------------------------------------------

donor_state_table <- cell_score_table %>%
  dplyr::group_by(
    sample_label,
    donor,
    geo_accession,
    group,
    cell_state_short,
    cell_state
  ) %>%
  dplyr::summarise(
    
    n_cells = dplyr::n(),
    
    total_RNA_UMI =
      sum(
        total_RNA_UMI,
        na.rm = TRUE
      ),
    
    MMP13_total_UMI =
      sum(
        MMP13_raw_count,
        na.rm = TRUE
      ),
    
    MMP13_positive_cells =
      sum(
        MMP13_detected,
        na.rm = TRUE
      ),
    
    MMP13_detection_percent =
      100 *
      MMP13_positive_cells /
      n_cells,
    
    dplyr::across(
      dplyr::all_of(score_columns),
      list(
        mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        median = ~ stats::median(
          .x,
          na.rm = TRUE
        )
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    
    MMP13_CPM =
      ifelse(
        total_RNA_UMI > 0,
        1e6 *
          MMP13_total_UMI /
          total_RNA_UMI,
        NA_real_
      ),
    
    log2_MMP13_CPM =
      log2(
        MMP13_CPM + 1
      ),
    
    reliable_state_sample =
      n_cells >= 50
  )

cat("\nDonor-state table:\n")
print(
  donor_state_table,
  n = 25,
  width = Inf
)

# ------------------------------------------------------------
# 11. Group-level descriptive pathway summary
#
# Means are first calculated per donor-state, then summarized
# across donors. Cells are not treated as independent replicates.
# ------------------------------------------------------------

mean_score_columns <- paste0(
  score_columns,
  "_mean"
)

group_state_summary <- donor_state_table %>%
  dplyr::filter(
    reliable_state_sample
  ) %>%
  dplyr::group_by(
    group,
    cell_state_short,
    cell_state
  ) %>%
  dplyr::summarise(
    
    donors =
      dplyr::n_distinct(
        donor
      ),
    
    mean_cells =
      mean(
        n_cells
      ),
    
    mean_MMP13_CPM =
      mean(
        MMP13_CPM,
        na.rm = TRUE
      ),
    
    mean_log2_MMP13_CPM =
      mean(
        log2_MMP13_CPM,
        na.rm = TRUE
      ),
    
    dplyr::across(
      dplyr::all_of(
        mean_score_columns
      ),
      list(
        mean = ~ mean(
          .x,
          na.rm = TRUE
        ),
        sd = ~ stats::sd(
          .x,
          na.rm = TRUE
        )
      ),
      .names = "{.col}_{.fn}"
    ),
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 12. Exploratory OA-versus-Normal donor-level tests
#
# Wilcoxon tests are run per cell state and pathway.
# These tests are exploratory and not multiple-cell tests.
# ------------------------------------------------------------

oa_normal_test_list <- list()

for (current_state in state_order) {
  
  current_data <- donor_state_table %>%
    dplyr::filter(
      cell_state_short ==
        current_state,
      reliable_state_sample
    )
  
  for (current_score in mean_score_columns) {
    
    normal_values <- current_data %>%
      dplyr::filter(
        group == "Normal"
      ) %>%
      dplyr::pull(
        dplyr::all_of(
          current_score
        )
      )
    
    oa_values <- current_data %>%
      dplyr::filter(
        group == "OA"
      ) %>%
      dplyr::pull(
        dplyr::all_of(
          current_score
        )
      )
    
    normal_values <- normal_values[
      is.finite(normal_values)
    ]
    
    oa_values <- oa_values[
      is.finite(oa_values)
    ]
    
    if (
      length(normal_values) >= 3 &&
      length(oa_values) >= 3
    ) {
      
      current_test <- suppressWarnings(
        stats::wilcox.test(
          x = oa_values,
          y = normal_values,
          exact = FALSE,
          conf.int = FALSE
        )
      )
      
      p_value <- current_test$p.value
      
    } else {
      
      p_value <- NA_real_
    }
    
    oa_normal_test_list[[
      paste(
        current_state,
        current_score,
        sep = "__"
      )
    ]] <- tibble::tibble(
      
      cell_state_short =
        current_state,
      
      pathway_score_column =
        current_score,
      
      n_normal =
        length(normal_values),
      
      n_OA =
        length(oa_values),
      
      median_normal =
        ifelse(
          length(normal_values) > 0,
          stats::median(
            normal_values
          ),
          NA_real_
        ),
      
      median_OA =
        ifelse(
          length(oa_values) > 0,
          stats::median(
            oa_values
          ),
          NA_real_
        ),
      
      OA_minus_Normal =
        median_OA -
        median_normal,
      
      p_value =
        p_value
    )
  }
}

oa_normal_tests <- dplyr::bind_rows(
  oa_normal_test_list
) %>%
  dplyr::mutate(
    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )

# ------------------------------------------------------------
# 13. Exploratory pathway-MMP13 correlations
#
# Spearman correlations are calculated:
#   A. across all reliable chondrocyte donor-state samples;
#   B. separately within each annotated chondrocyte state.
#
# Vascular/immune cells are excluded from the main relationship.
# ------------------------------------------------------------

chondrocyte_states <- setdiff(
  state_order,
  "Vascular/immune"
)

correlation_input <- donor_state_table %>%
  dplyr::filter(
    reliable_state_sample,
    cell_state_short %in%
      chondrocyte_states
  )

correlation_list <- list()

for (current_score in mean_score_columns) {
  
  x_values <- correlation_input[[current_score]]
  
  y_values <- correlation_input$
    log2_MMP13_CPM
  
  keep_values <- is.finite(x_values) &
    is.finite(y_values)
  
  if (sum(keep_values) >= 6) {
    
    current_test <- suppressWarnings(
      stats::cor.test(
        x = x_values[keep_values],
        y = y_values[keep_values],
        method = "spearman",
        exact = FALSE
      )
    )
    
    rho_value <- unname(
      current_test$estimate
    )
    
    p_value <- current_test$p.value
    
  } else {
    
    rho_value <- NA_real_
    p_value <- NA_real_
  }
  
  correlation_list[[
    paste0(
      "ALL__",
      current_score
    )
  ]] <- tibble::tibble(
    
    analysis_scope =
      "All chondrocyte donor-state samples",
    
    cell_state_short =
      "All_chondrocytes",
    
    pathway_score_column =
      current_score,
    
    n_samples =
      sum(keep_values),
    
    spearman_rho =
      rho_value,
    
    p_value =
      p_value
  )
  
  for (current_state in chondrocyte_states) {
    
    current_state_data <- correlation_input %>%
      dplyr::filter(
        cell_state_short ==
          current_state
      )
    
    x_state <- current_state_data[[current_score]]
    
    y_state <- current_state_data$
      log2_MMP13_CPM
    
    keep_state <- is.finite(x_state) &
      is.finite(y_state)
    
    nonzero_mmp13 <- sum(
      y_state[keep_state] > 0
    )
    
    if (
      sum(keep_state) >= 6 &&
      nonzero_mmp13 >= 3
    ) {
      
      state_test <- suppressWarnings(
        stats::cor.test(
          x = x_state[keep_state],
          y = y_state[keep_state],
          method = "spearman",
          exact = FALSE
        )
      )
      
      state_rho <- unname(
        state_test$estimate
      )
      
      state_p <- state_test$p.value
      
    } else {
      
      state_rho <- NA_real_
      state_p <- NA_real_
    }
    
    correlation_list[[
      paste(
        current_state,
        current_score,
        sep = "__"
      )
    ]] <- tibble::tibble(
      
      analysis_scope =
        "Within annotated state",
      
      cell_state_short =
        current_state,
      
      pathway_score_column =
        current_score,
      
      n_samples =
        sum(keep_state),
      
      nonzero_MMP13_samples =
        nonzero_mmp13,
      
      spearman_rho =
        state_rho,
      
      p_value =
        state_p
    )
  }
}

pathway_mmp13_correlations <- dplyr::bind_rows(
  correlation_list
) %>%
  dplyr::mutate(
    FDR =
      stats::p.adjust(
        p_value,
        method = "BH"
      )
  )

cat("\nExploratory pathway-MMP13 correlations:\n")
print(
  pathway_mmp13_correlations,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 14. Figure S17A:
# Donor-level pathway activity heatmap by state and group
# ------------------------------------------------------------

heatmap_table <- donor_state_table %>%
  dplyr::filter(
    reliable_state_sample
  ) %>%
  dplyr::select(
    group,
    donor,
    cell_state_short,
    dplyr::all_of(
      mean_score_columns
    )
  ) %>%
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        mean_score_columns
      ),
    names_to =
      "pathway",
    values_to =
      "score"
  ) %>%
  dplyr::mutate(
    pathway = stringr::str_remove(
      pathway,
      "_module_score_mean$"
    )
  ) %>%
  dplyr::group_by(
    pathway
  ) %>%
  dplyr::mutate(
    pathway_z = as.numeric(
      scale(score)
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(
    group,
    cell_state_short,
    pathway
  ) %>%
  dplyr::summarise(
    donors =
      dplyr::n_distinct(
        donor
      ),
    mean_pathway_z =
      mean(
        pathway_z,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    state_group = paste(
      cell_state_short,
      group,
      sep = " | "
    ),
    state_group = factor(
      state_group,
      levels = as.vector(
        rbind(
          paste(
            state_order,
            "Normal",
            sep = " | "
          ),
          paste(
            state_order,
            "OA",
            sep = " | "
          )
        )
      )
    ),
    pathway = factor(
      pathway,
      levels = c(
        "NFkB",
        "p38_MAPK",
        "Inflammatory",
        "ECM_degradation"
      )
    )
  )

pathway_heatmap <- ggplot(
  heatmap_table,
  aes(
    x = pathway,
    y = state_group,
    fill = mean_pathway_z
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.30
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        mean_pathway_z
      )
    ),
    size = 2.7
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Mean donor-level\nz score"
  ) +
  labs(
    title =
      "Donor-level pathway activity across cartilage cell states",
    subtitle =
      "Scores were aggregated by donor before group summarization",
    x =
      "Pathway program",
    y =
      "Cell state and group"
  ) +
  theme_bw(
    base_size = 10.5
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  )

# ------------------------------------------------------------
# 15. Figure S17B:
# Exploratory donor-state relationship with MMP13
# ------------------------------------------------------------

scatter_table <- correlation_input %>%
  dplyr::select(
    sample_label,
    donor,
    group,
    cell_state_short,
    log2_MMP13_CPM,
    dplyr::all_of(
      mean_score_columns
    )
  ) %>%
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        mean_score_columns
      ),
    names_to =
      "pathway",
    values_to =
      "pathway_score"
  ) %>%
  dplyr::mutate(
    pathway = stringr::str_remove(
      pathway,
      "_module_score_mean$"
    ),
    pathway = factor(
      pathway,
      levels = c(
        "NFkB",
        "p38_MAPK",
        "Inflammatory",
        "ECM_degradation"
      )
    )
  )

state_palette <- c(
  "RegC" = "#E64B35",
  "EC" = "#E69F00",
  "RepC" = "#7A9A01",
  "PreFC" = "#009E73",
  "FC" = "#00A6A6",
  "MTC-like" = "#0072B2",
  "HomC" = "#56B4E9",
  "HTC" = "#BC3C82"
)

pathway_scatter <- ggplot(
  scatter_table,
  aes(
    x = pathway_score,
    y = log2_MMP13_CPM,
    color = cell_state_short,
    shape = group
  )
) +
  geom_point(
    size = 2.0,
    alpha = 0.80
  ) +
  geom_smooth(
    aes(
      group = 1
    ),
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "grey25",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  facet_wrap(
    ~ pathway,
    scales = "free_x",
    ncol = 2
  ) +
  scale_color_manual(
    values = state_palette,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "Normal" = 16,
      "OA" = 17
    )
  ) +
  labs(
    title =
      "Exploratory donor-state relationships between pathway activity and MMP13",
    subtitle =
      "Each point is one donor within one annotated chondrocyte state",
    x =
      "Mean pathway module score",
    y =
      expression(
        log[2] * "(MMP13 CPM + 1)"
      ),
    color =
      "Cell state",
    shape =
      "Group"
  ) +
  theme_bw(
    base_size = 10.5
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    legend.position =
      "right"
  )

# ------------------------------------------------------------
# 16. Figure S17C:
# FC and HTC focus because they are biologically relevant to MMP13
# ------------------------------------------------------------

focus_states <- c(
  "FC",
  "HTC"
)

focus_table <- donor_state_table %>%
  dplyr::filter(
    reliable_state_sample,
    cell_state_short %in%
      focus_states
  ) %>%
  dplyr::select(
    sample_label,
    donor,
    group,
    cell_state_short,
    log2_MMP13_CPM,
    dplyr::all_of(
      mean_score_columns
    )
  ) %>%
  tidyr::pivot_longer(
    cols =
      dplyr::all_of(
        mean_score_columns
      ),
    names_to =
      "pathway",
    values_to =
      "pathway_score"
  ) %>%
  dplyr::mutate(
    pathway = stringr::str_remove(
      pathway,
      "_module_score_mean$"
    ),
    pathway = factor(
      pathway,
      levels = c(
        "NFkB",
        "p38_MAPK",
        "Inflammatory",
        "ECM_degradation"
      )
    )
  )

focus_pathway_plot <- ggplot(
  focus_table,
  aes(
    x = group,
    y = pathway_score,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA,
    alpha = 0.70
  ) +
  geom_point(
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 1.9,
    alpha = 0.85
  ) +
  facet_grid(
    cell_state_short ~ pathway,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "Normal" = "#D55E00",
      "OA" = "#009E9E"
    )
  ) +
  labs(
    title =
      "Donor-level pathway activity in fibrocartilage and hypertrophic states",
    subtitle =
      "Descriptive comparison; n = 6 donors per disease group when available",
    x =
      NULL,
    y =
      "Mean pathway module score",
    fill =
      "Group"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    legend.position =
      "top",
    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      )
  )

# ------------------------------------------------------------
# 17. Save figures
# ------------------------------------------------------------

figure_list <- list(
  
  FigureS17A_donor_level_pathway_heatmap =
    pathway_heatmap,
  
  FigureS17B_pathway_MMP13_scatter =
    pathway_scatter,
  
  FigureS17C_FC_HTC_pathway_activity =
    focus_pathway_plot
)

figure_sizes <- list(
  c(8.5, 9.0),
  c(12.0, 8.0),
  c(12.0, 6.2)
)

for (current_index in seq_along(figure_list)) {
  
  current_name <- names(figure_list)[current_index]
  current_plot <- figure_list[[current_index]]
  current_size <- figure_sizes[[current_index]]
  
  ggsave(
    filename = file.path(
      figure_folder,
      paste0(
        current_name,
        ".pdf"
      )
    ),
    plot = current_plot,
    width = current_size[1],
    height = current_size[2],
    units = "in",
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(
      figure_folder,
      paste0(
        current_name,
        ".tiff"
      )
    ),
    plot = current_plot,
    width = current_size[1],
    height = current_size[2],
    units = "in",
    dpi = 600,
    compression = "lzw",
    limitsize = FALSE
  )
}

# ------------------------------------------------------------
# 18. Save tables
# ------------------------------------------------------------

utils::write.csv(
  gene_set_qc,
  file = file.path(
    table_folder,
    "GSE220243_NFkB_p38_gene_set_QC.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  donor_state_table,
  file = file.path(
    table_folder,
    "GSE220243_NFkB_p38_donor_state_table.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  group_state_summary,
  file = file.path(
    table_folder,
    "GSE220243_NFkB_p38_group_state_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  oa_normal_tests,
  file = file.path(
    table_folder,
    "GSE220243_NFkB_p38_OA_vs_Normal_tests.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  pathway_mmp13_correlations,
  file = file.path(
    table_folder,
    "GSE220243_NFkB_p38_MMP13_correlations.csv"
  ),
  row.names = FALSE
)

summary_workbook <- file.path(
  table_folder,
  "GSE220243_NFkB_p38_MMP13_exploratory_summary.xlsx"
)

openxlsx::write.xlsx(
  list(
    
    Gene_set_QC =
      as.data.frame(
        gene_set_qc
      ),
    
    Donor_state =
      as.data.frame(
        donor_state_table
      ),
    
    Group_state =
      as.data.frame(
        group_state_summary
      ),
    
    OA_vs_Normal_tests =
      as.data.frame(
        oa_normal_tests
      ),
    
    MMP13_correlations =
      as.data.frame(
        pathway_mmp13_correlations
      ),
    
    Package_versions =
      as.data.frame(
        package_versions
      )
  ),
  file =
    summary_workbook,
  overwrite =
    TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    table_folder,
    "GSE220243_NFkB_p38_sessionInfo.txt"
  )
)

# ------------------------------------------------------------
# 19. Final summary
# ------------------------------------------------------------

message("")
message("================================================")
message("SCRIPT 19A SELESAI")
message("================================================")

message(
  "MSigDB version             : ",
  paste(
    msigdb_version,
    collapse = ", "
  )
)

message(
  "Pathway programs scored    : ",
  paste(
    names(matched_gene_sets),
    collapse = ", "
  )
)

message(
  "Biological replicate       : DONOR"
)

message(
  "Donor-state rows           : ",
  nrow(
    donor_state_table
  )
)

message(
  "Exploratory figures        : ",
  figure_folder
)

message(
  "Summary workbook           : ",
  summary_workbook
)

message(
  "Figure 1E status           : NOT FINAL"
)

message(
  "Next step                  : inspect S17A-S17C and choose the defensible Figure 1E structure"
)

message("================================================")