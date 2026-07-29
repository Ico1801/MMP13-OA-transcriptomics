# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 17B_GSE220243_memory_efficient_RPCA_resume.R
# PURPOSE :
#   Resume after std::bad_alloc in IntegrateData()
#   Uses Seurat v5 low-dimensional RPCA integration
#
# IMPORTANT:
#   - Starts from the filtered singlet list saved by script 17.
#   - Does NOT repeat adaptive QC or scDblFinder.
#   - Uses LogNormalize + IntegrateLayers(RPCAIntegration).
#   - RNA counts/data are retained for expression analysis.
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
  stop("Project folder tidak ditemukan:\n", project_root)
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
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "patchwork",
  "openxlsx",
  "future"
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
    paste(still_missing, collapse = ", ")
  )
}

library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(openxlsx)
library(future)

future::plan("sequential")
options(future.globals.maxSize = 12 * 1024^3)
options(Seurat.object.assay.version = "v5")

package_versions <- tibble::tibble(
  package = required_packages,
  version = vapply(
    required_packages,
    function(x) as.character(packageVersion(x)),
    FUN.VALUE = character(1)
  )
)

cat("\nPackage versions:\n")
print(package_versions, n = Inf, width = Inf)

# ------------------------------------------------------------
# 3. Input and output paths
# ------------------------------------------------------------

processed_folder <- "data_processed"
table_folder <- "results/tables"

supplementary_figure_folder <- paste0(
  "results/figures/",
  "Figure 1 Suplementary/",
  "S13"
)

dir.create(processed_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(table_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(
  supplementary_figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

singlet_file <- file.path(
  processed_folder,
  "GSE220243_cartilage_filtered_singlet_list.rds"
)

raw_file <- file.path(
  processed_folder,
  "GSE220243_cartilage_raw_seurat_list.rds"
)

basic_qc_file <- file.path(
  processed_folder,
  "GSE220243_cartilage_basic_QC_seurat_list.rds"
)

if (!file.exists(singlet_file)) {
  stop(
    "Filtered singlet list tidak ditemukan:\n",
    singlet_file,
    "\nJangan jalankan script ini sebelum script 17 mencapai ",
    "'Filtered singlet Seurat list berhasil disimpan.'"
  )
}

# ------------------------------------------------------------
# 4. Load filtered singlet objects
# ------------------------------------------------------------

message("Memuat filtered singlet list...")

singlet_list <- readRDS(singlet_file)

if (!is.list(singlet_list)) {
  stop("Input filtered singlet bukan list.")
}

expected_sample_order <- c(
  "Normal1", "Normal2", "Normal3",
  "Normal4", "Normal5", "Normal6",
  "OA1", "OA2", "OA3",
  "OA4", "OA5", "OA6"
)

if (is.null(names(singlet_list))) {
  sample_names <- vapply(
    singlet_list,
    function(x) as.character(unique(x$sample_label)[1]),
    FUN.VALUE = character(1)
  )
  names(singlet_list) <- sample_names
}

missing_samples <- setdiff(
  expected_sample_order,
  names(singlet_list)
)

if (length(missing_samples) > 0) {
  stop(
    "Sampel berikut tidak ditemukan: ",
    paste(missing_samples, collapse = ", ")
  )
}

singlet_list <- singlet_list[expected_sample_order]

stopifnot(
  length(singlet_list) == 12,
  all(
    vapply(
      singlet_list,
      inherits,
      logical(1),
      what = "Seurat"
    )
  )
)

final_singlet_counts <- vapply(
  singlet_list,
  ncol,
  numeric(1)
)

message(
  "Filtered singlets loaded: ",
  format(sum(final_singlet_counts), big.mark = ",")
)

# ------------------------------------------------------------
# 5. Reconstruct retention summary without repeating QC
# ------------------------------------------------------------

raw_counts <- rep(NA_integer_, 12)
names(raw_counts) <- expected_sample_order

basic_qc_counts <- rep(NA_integer_, 12)
names(basic_qc_counts) <- expected_sample_order

if (file.exists(raw_file)) {
  message("Membaca raw list untuk jumlah sel awal...")
  raw_list_temp <- readRDS(raw_file)
  
  if (is.null(names(raw_list_temp))) {
    names(raw_list_temp) <- vapply(
      raw_list_temp,
      function(x) as.character(unique(x$sample_label)[1]),
      FUN.VALUE = character(1)
    )
  }
  
  raw_counts[names(raw_list_temp)] <- vapply(
    raw_list_temp,
    ncol,
    numeric(1)
  )
  
  rm(raw_list_temp)
  gc()
}

if (file.exists(basic_qc_file)) {
  message("Membaca basic-QC list untuk jumlah sel pasca-QC...")
  basic_list_temp <- readRDS(basic_qc_file)
  
  if (is.null(names(basic_list_temp))) {
    names(basic_list_temp) <- vapply(
      basic_list_temp,
      function(x) as.character(unique(x$sample_label)[1]),
      FUN.VALUE = character(1)
    )
  }
  
  basic_qc_counts[names(basic_list_temp)] <- vapply(
    basic_list_temp,
    ncol,
    numeric(1)
  )
  
  rm(basic_list_temp)
  gc()
}

group_vector <- c(
  rep("Normal", 6),
  rep("OA", 6)
)

retention_summary <- tibble::tibble(
  sample_label = expected_sample_order,
  group = factor(
    group_vector,
    levels = c("Normal", "OA")
  ),
  raw_cells = as.integer(
    raw_counts[expected_sample_order]
  ),
  cells_after_basic_QC = as.integer(
    basic_qc_counts[expected_sample_order]
  ),
  predicted_doublets = as.integer(
    basic_qc_counts[expected_sample_order] -
      final_singlet_counts[expected_sample_order]
  ),
  final_singlets = as.integer(
    final_singlet_counts[expected_sample_order]
  )
) %>%
  dplyr::mutate(
    predicted_doublet_percent =
      100 * predicted_doublets / cells_after_basic_QC,
    
    final_retention_percent =
      100 * final_singlets / raw_cells
  )

cat("\nReconstructed retention summary:\n")
print(retention_summary, n = Inf, width = Inf)

# ------------------------------------------------------------
# 6. Donor-level MMP13 detection before integration
# ------------------------------------------------------------

mmp13_donor_list <- list()

for (current_sample in expected_sample_order) {
  
  current_object <- singlet_list[[current_sample]]
  DefaultAssay(current_object) <- "RNA"
  
  current_group <- as.character(
    unique(current_object$group)[1]
  )
  
  current_donor <- as.character(
    unique(current_object$donor)[1]
  )
  
  current_geo <- as.character(
    unique(current_object$geo_accession)[1]
  )
  
  mmp13_candidates <- rownames(current_object)[
    toupper(rownames(current_object)) == "MMP13"
  ]
  
  if (length(mmp13_candidates) == 0) {
    
    positive_cells <- 0
    total_counts <- 0
    mmp13_present <- FALSE
    
  } else {
    
    mmp13_feature <- mmp13_candidates[1]
    
    current_counts <- SeuratObject::LayerData(
      current_object,
      assay = "RNA",
      layer = "counts"
    )
    
    mmp13_values <- as.numeric(
      current_counts[
        mmp13_feature,
        ,
        drop = TRUE
      ]
    )
    
    positive_cells <- sum(mmp13_values > 0)
    total_counts <- sum(mmp13_values)
    mmp13_present <- TRUE
  }
  
  mmp13_donor_list[[current_sample]] <- tibble::tibble(
    sample_label = current_sample,
    group = current_group,
    donor = current_donor,
    geo_accession = current_geo,
    final_singlets = ncol(current_object),
    MMP13_present = mmp13_present,
    MMP13_positive_cells = positive_cells,
    MMP13_detection_percent =
      100 * positive_cells / ncol(current_object),
    MMP13_total_counts = total_counts
  )
}

mmp13_donor_summary <- dplyr::bind_rows(
  mmp13_donor_list
)

cat("\nMMP13 detection among final singlets:\n")
print(mmp13_donor_summary, n = Inf, width = Inf)

# ------------------------------------------------------------
# 7. Merge 12 donors
# ------------------------------------------------------------

message("Menggabungkan 12 donor singlet...")

merged_object <- merge(
  x = singlet_list[[1]],
  y = singlet_list[-1],
  project = "GSE220243_cartilage_memory_efficient",
  merge.data = FALSE
)

rm(singlet_list)
gc()

DefaultAssay(merged_object) <- "RNA"

merged_object$group <- factor(
  merged_object$group,
  levels = c("Normal", "OA")
)

merged_object$sample_label <- factor(
  merged_object$sample_label,
  levels = expected_sample_order
)

cat("\nMerged object:\n")
print(merged_object)

cat("\nCells per donor:\n")
print(table(merged_object$sample_label))

stopifnot(
  ncol(merged_object) ==
    sum(retention_summary$final_singlets)
)

# ------------------------------------------------------------
# 8. Ensure donor-specific RNA layers
# ------------------------------------------------------------

rna_layers <- SeuratObject::Layers(
  merged_object[["RNA"]]
)

cat("\nRNA layers before preprocessing:\n")
print(rna_layers)

count_layer_number <- sum(
  grepl("^counts", rna_layers)
)

if (count_layer_number == 1) {
  
  message(
    "RNA counts masih satu layer; membaginya berdasarkan donor..."
  )
  
  merged_object[["RNA"]] <- split(
    merged_object[["RNA"]],
    f = merged_object$sample_label
  )
}

cat("\nRNA layers used for integration:\n")
print(
  SeuratObject::Layers(
    merged_object[["RNA"]]
  )
)

# ------------------------------------------------------------
# 9. Memory-efficient preprocessing
# ------------------------------------------------------------

message("Log-normalizing donor layers...")

merged_object <- NormalizeData(
  merged_object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

message("Selecting 2,000 variable features...")

merged_object <- FindVariableFeatures(
  merged_object,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

integration_features <- VariableFeatures(
  merged_object
)

if (length(integration_features) > 2000) {
  integration_features <- integration_features[1:2000]
}

if (length(integration_features) < 1000) {
  stop(
    "Variable features terlalu sedikit: ",
    length(integration_features)
  )
}

message(
  "Variable features retained: ",
  length(integration_features)
)

message("Scaling only integration features...")

merged_object <- ScaleData(
  merged_object,
  assay = "RNA",
  features = integration_features,
  block.size = 500,
  verbose = TRUE
)

message("Running unintegrated PCA...")

merged_object <- RunPCA(
  merged_object,
  assay = "RNA",
  features = integration_features,
  npcs = 30,
  reduction.name = "pca",
  reduction.key = "PC_",
  seed.use = 20260723,
  verbose = TRUE
)

# ------------------------------------------------------------
# 10. Unintegrated UMAP for diagnostic comparison
# ------------------------------------------------------------

message("Running unintegrated UMAP...")

merged_object <- RunUMAP(
  merged_object,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.unintegrated",
  reduction.key = "UMAPUNINT_",
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 20260723,
  verbose = TRUE
)

# ------------------------------------------------------------
# 11. Seurat v5 low-dimensional RPCA integration
# ------------------------------------------------------------

message(
  "Running memory-efficient Seurat v5 RPCA integration..."
)

# Layer order follows the donor merge order:
# 1 = Normal1 and 9 = OA3.
reference_indices <- c(1, 9)

merged_object <- Seurat::IntegrateLayers(
  object = merged_object,
  method = Seurat::RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  assay = "RNA",
  features = integration_features,
  normalization.method = "LogNormalize",
  dims = 1:30,
  reference = reference_indices,
  k.anchor = 5,
  k.filter = 200,
  k.weight = 50,
  preserve.order = FALSE,
  verbose = TRUE
)

if (
  !"integrated.rpca" %in%
  names(merged_object@reductions)
) {
  stop("Reduction integrated.rpca tidak berhasil dibuat.")
}

message("Low-dimensional RPCA integration selesai.")

# ------------------------------------------------------------
# 12. Integrated clustering and UMAP
# ------------------------------------------------------------

merged_object <- FindNeighbors(
  merged_object,
  reduction = "integrated.rpca",
  dims = 1:30,
  graph.name = c("rpca_nn", "rpca_snn"),
  verbose = TRUE
)

merged_object <- FindClusters(
  merged_object,
  graph.name = "rpca_snn",
  resolution = 0.4,
  cluster.name = "analysis_cluster",
  algorithm = 1,
  random.seed = 20260723,
  verbose = TRUE
)

merged_object <- RunUMAP(
  merged_object,
  reduction = "integrated.rpca",
  dims = 1:30,
  reduction.name = "umap.rpca",
  reduction.key = "UMAPRPCA_",
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 20260723,
  verbose = TRUE
)

Idents(merged_object) <- "analysis_cluster"

cat("\nIntegrated clusters:\n")
print(table(merged_object$analysis_cluster))

# ------------------------------------------------------------
# 13. Join RNA layers for biological visualization
# ------------------------------------------------------------

if (
  inherits(
    merged_object[["RNA"]],
    "Assay5"
  )
) {
  message("Joining RNA layers...")
  merged_object <- SeuratObject::JoinLayers(
    merged_object,
    assay = "RNA"
  )
}

DefaultAssay(merged_object) <- "RNA"

# Re-normalize after joining to guarantee one complete data layer.
merged_object <- NormalizeData(
  merged_object,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

cat("\nRNA layers after joining:\n")
print(
  SeuratObject::Layers(
    merged_object[["RNA"]]
  )
)

# ------------------------------------------------------------
# 14. Cell-level MMP13 data
# ------------------------------------------------------------

mmp13_candidates <- rownames(merged_object)[
  toupper(rownames(merged_object)) == "MMP13"
]

if (length(mmp13_candidates) == 0) {
  stop("MMP13 tidak ditemukan pada merged object.")
}

mmp13_feature <- mmp13_candidates[1]

rna_counts <- SeuratObject::LayerData(
  merged_object,
  assay = "RNA",
  layer = "counts"
)

rna_data <- SeuratObject::LayerData(
  merged_object,
  assay = "RNA",
  layer = "data"
)

mmp13_raw_counts <- as.numeric(
  rna_counts[
    mmp13_feature,
    ,
    drop = TRUE
  ]
)

mmp13_log_expression <- as.numeric(
  rna_data[
    mmp13_feature,
    ,
    drop = TRUE
  ]
)

mmp13_cell_table <- merged_object[[]] %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::mutate(
    MMP13_raw_count = mmp13_raw_counts,
    MMP13_detected = MMP13_raw_count > 0,
    MMP13_log_normalized = mmp13_log_expression
  )

mmp13_group_summary <- mmp13_cell_table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    donors = dplyr::n_distinct(donor),
    final_singlets = dplyr::n(),
    MMP13_positive_cells = sum(MMP13_detected),
    MMP13_detection_percent =
      100 * MMP13_positive_cells / final_singlets,
    MMP13_total_counts = sum(MMP13_raw_count),
    mean_MMP13_log_normalized =
      mean(MMP13_log_normalized),
    .groups = "drop"
  )

cat("\nMMP13 group summary:\n")
print(mmp13_group_summary, n = Inf, width = Inf)

# ------------------------------------------------------------
# 15. Supplementary figures
# ------------------------------------------------------------

retention_plot_data <- retention_summary %>%
  dplyr::select(
    sample_label,
    group,
    raw_cells,
    cells_after_basic_QC,
    final_singlets
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      raw_cells,
      cells_after_basic_QC,
      final_singlets
    ),
    names_to = "stage",
    values_to = "cells"
  ) %>%
  dplyr::mutate(
    stage = factor(
      stage,
      levels = c(
        "raw_cells",
        "cells_after_basic_QC",
        "final_singlets"
      ),
      labels = c(
        "Initial cells",
        "After adaptive QC",
        "Final singlets"
      )
    ),
    sample_label = factor(
      sample_label,
      levels = expected_sample_order
    )
  )

retention_plot <- ggplot(
  retention_plot_data,
  aes(
    x = sample_label,
    y = cells,
    fill = stage
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.72
  ) +
  facet_grid(
    ~ group,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Cell retention during quality control",
    x = "Cartilage donor",
    y = "Cells",
    fill = "Processing stage"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

doublet_plot <- ggplot(
  retention_summary,
  aes(
    x = factor(
      sample_label,
      levels = expected_sample_order
    ),
    y = predicted_doublet_percent,
    fill = group
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        predicted_doublet_percent
      )
    ),
    vjust = -0.3,
    size = 3
  ) +
  facet_grid(
    ~ group,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Predicted doublet rate by donor",
    x = "Cartilage donor",
    y = "Predicted doublets (%)",
    fill = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

umap_before <- DimPlot(
  merged_object,
  reduction = "umap.unintegrated",
  group.by = "sample_label",
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  labs(
    title = "Before integration",
    subtitle = "RNA PCA–UMAP"
  ) +
  theme(legend.position = "none")

umap_after <- DimPlot(
  merged_object,
  reduction = "umap.rpca",
  group.by = "sample_label",
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  labs(
    title = "After RPCA integration",
    subtitle = "Cells colored by donor"
  )

integration_comparison <- (
  umap_before +
    umap_after
) +
  patchwork::plot_annotation(
    title = paste(
      "Donor structure before and after",
      "memory-efficient RPCA integration"
    )
  )

umap_group <- DimPlot(
  merged_object,
  reduction = "umap.rpca",
  group.by = "group",
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  labs(
    title = "Integrated UMAP by disease group"
  )

umap_cluster <- DimPlot(
  merged_object,
  reduction = "umap.rpca",
  group.by = "analysis_cluster",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  labs(
    title = "Integrated transcriptional clusters"
  ) +
  theme(legend.position = "none")

integrated_umap_panel <- (
  umap_group +
    umap_cluster
) +
  patchwork::plot_annotation(
    title = "Integrated single-cell landscape of human cartilage"
  )

mmp13_detection_plot <- ggplot(
  mmp13_donor_summary,
  aes(
    x = factor(
      sample_label,
      levels = expected_sample_order
    ),
    y = MMP13_detection_percent,
    fill = group
  )
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f%%",
        MMP13_detection_percent
      )
    ),
    vjust = -0.3,
    size = 2.8
  ) +
  facet_grid(
    ~ group,
    scales = "free_x",
    space = "free_x"
  ) +
  labs(
    title = "Post-QC MMP13 detection across cartilage donors",
    x = "Cartilage donor",
    y = "MMP13-positive singlets (%)",
    fill = "Group"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

mmp13_feature_plot <- FeaturePlot(
  merged_object,
  features = mmp13_feature,
  reduction = "umap.rpca",
  min.cutoff = "q05",
  max.cutoff = "q99",
  raster = TRUE,
  raster.dpi = c(300, 300)
) +
  labs(
    title = "Preliminary MMP13 localization",
    subtitle = "Final interpretation requires cell-type annotation"
  )

plot_list <- list(
  FigureS13A_GSE220243_cell_retention =
    retention_plot,
  
  FigureS13B_GSE220243_doublet_rate =
    doublet_plot,
  
  FigureS13C_GSE220243_RPCA_integration_comparison =
    integration_comparison,
  
  FigureS13D_GSE220243_integrated_UMAP =
    integrated_umap_panel,
  
  FigureS13E_GSE220243_postQC_MMP13_detection =
    mmp13_detection_plot,
  
  FigureS13F_GSE220243_preliminary_MMP13_featureplot =
    mmp13_feature_plot
)

plot_dimensions <- list(
  c(10, 5.8),
  c(9, 5.5),
  c(13, 6),
  c(13, 6),
  c(9, 5.8),
  c(7, 6)
)

for (current_index in seq_along(plot_list)) {
  
  current_name <- names(plot_list)[current_index]
  current_plot <- plot_list[[current_index]]
  current_size <- plot_dimensions[[current_index]]
  
  ggsave(
    filename = file.path(
      supplementary_figure_folder,
      paste0(current_name, ".pdf")
    ),
    plot = current_plot,
    width = current_size[1],
    height = current_size[2]
  )
  
  ggsave(
    filename = file.path(
      supplementary_figure_folder,
      paste0(current_name, ".tiff")
    ),
    plot = current_plot,
    width = current_size[1],
    height = current_size[2],
    dpi = 600,
    compression = "lzw"
  )
}

# ------------------------------------------------------------
# 16. Save tables and final object
# ------------------------------------------------------------

utils::write.csv(
  retention_summary,
  file = file.path(
    table_folder,
    "GSE220243_memory_efficient_retention_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mmp13_donor_summary,
  file = file.path(
    table_folder,
    "GSE220243_memory_efficient_MMP13_donor_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mmp13_group_summary,
  file = file.path(
    table_folder,
    "GSE220243_memory_efficient_MMP13_group_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mmp13_cell_table,
  file = file.path(
    table_folder,
    "GSE220243_memory_efficient_MMP13_cell_level.csv"
  ),
  row.names = FALSE
)

openxlsx::write.xlsx(
  list(
    Retention =
      as.data.frame(retention_summary),
    
    MMP13_by_donor =
      as.data.frame(mmp13_donor_summary),
    
    MMP13_by_group =
      as.data.frame(mmp13_group_summary),
    
    Package_versions =
      as.data.frame(package_versions)
  ),
  file = file.path(
    table_folder,
    "GSE220243_memory_efficient_RPCA_summary.xlsx"
  ),
  overwrite = TRUE
)

final_object_file <- file.path(
  processed_folder,
  paste0(
    "GSE220243_cartilage_",
    "memory_efficient_RPCA_integrated.rds"
  )
)

saveRDS(
  merged_object,
  file = final_object_file,
  compress = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    table_folder,
    "GSE220243_memory_efficient_RPCA_sessionInfo.txt"
  )
)

# ------------------------------------------------------------
# 17. Final summary
# ------------------------------------------------------------

final_cells <- ncol(merged_object)
normal_cells <- sum(merged_object$group == "Normal")
oa_cells <- sum(merged_object$group == "OA")
cluster_number <- dplyr::n_distinct(
  merged_object$analysis_cluster
)

message("")
message("================================================")
message("GSE220243 MEMORY-EFFICIENT RPCA SELESAI")
message("================================================")

message(
  "Initial cells               : ",
  format(
    sum(retention_summary$raw_cells),
    big.mark = ","
  )
)

message(
  "After adaptive QC           : ",
  format(
    sum(retention_summary$cells_after_basic_QC),
    big.mark = ","
  )
)

message(
  "Predicted doublets removed  : ",
  format(
    sum(retention_summary$predicted_doublets),
    big.mark = ","
  )
)

message(
  "Final singlet cells         : ",
  format(final_cells, big.mark = ",")
)

message(
  "Normal singlet cells        : ",
  format(normal_cells, big.mark = ",")
)

message(
  "OA singlet cells            : ",
  format(oa_cells, big.mark = ",")
)

message(
  "Integrated clusters         : ",
  cluster_number
)

message(
  "Final integrated object     : ",
  final_object_file
)

message(
  "Supplementary figures       : ",
  supplementary_figure_folder
)

message(
  "Cell-type annotation        : FALSE"
)

message(
  "Figure 1D final             : FALSE"
)

message("================================================")