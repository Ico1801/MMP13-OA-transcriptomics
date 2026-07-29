# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 18A_GSE220243_cluster_marker_discovery.R
# PURPOSE :
#   1. Validate donor and group representation in each cluster
#   2. Identify positive cluster markers from the RNA assay
#   3. Visualize established cartilage cell-state markers
#   4. Localize MMP13-positive cells by cluster, donor, and group
#   5. Generate evidence needed for manual cell-type annotation
#
# IMPORTANT:
#   - Integrated RPCA coordinates are used only for clustering
#     and visualization.
#   - RNA expression is used for marker discovery.
#   - This script does NOT assign final cell-type names.
#   - OA-versus-Normal differential expression is NOT tested here.
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
  "ggplot2",
  "dplyr",
  "tibble",
  "tidyr",
  "stringr",
  "forcats",
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
    paste(
      still_missing,
      collapse = ", "
    )
  )
}

library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(forcats)
library(patchwork)
library(openxlsx)
library(future)

future::plan("sequential")
options(future.globals.maxSize = 12 * 1024^3)

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
  paste0(
    "GSE220243_cartilage_",
    "memory_efficient_RPCA_integrated.rds"
  )
)

table_folder <- "results/tables"

supplementary_figure_folder <- paste0(
  "results/figures/",
  "Figure 1 Suplementary/",
  "S14"
)

dir.create(
  table_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  supplementary_figure_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Integrated object tidak ditemukan:\n",
    input_file
  )
}

# ------------------------------------------------------------
# 4. Load and validate integrated object
# ------------------------------------------------------------

message("Memuat integrated Seurat object...")

cartilage_object <- readRDS(input_file)

if (!inherits(cartilage_object, "Seurat")) {
  stop("Input bukan objek Seurat.")
}

required_metadata <- c(
  "sample_label",
  "group",
  "donor",
  "geo_accession",
  "analysis_cluster"
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

if (
  !"umap.rpca" %in%
  names(cartilage_object@reductions)
) {
  stop(
    "Reduction 'umap.rpca' tidak ditemukan."
  )
}

if (
  !"integrated.rpca" %in%
  names(cartilage_object@reductions)
) {
  stop(
    "Reduction 'integrated.rpca' tidak ditemukan."
  )
}

DefaultAssay(cartilage_object) <- "RNA"

if (
  inherits(
    cartilage_object[["RNA"]],
    "Assay5"
  ) &&
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

cartilage_object$analysis_cluster <- factor(
  as.character(
    cartilage_object$analysis_cluster
  ),
  levels = sort(
    unique(
      as.integer(
        as.character(
          cartilage_object$analysis_cluster
        )
      )
    )
  )
)

Idents(cartilage_object) <- "analysis_cluster"

cluster_levels <- levels(
  cartilage_object$analysis_cluster
)

cat("\nIntegrated object:\n")
print(cartilage_object)

cat("\nCells per cluster:\n")
print(
  table(
    cartilage_object$analysis_cluster
  )
)

# ------------------------------------------------------------
# 5. Cluster composition by donor and group
# ------------------------------------------------------------

metadata_table <- cartilage_object[[]] %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::mutate(
    analysis_cluster =
      as.character(
        analysis_cluster
      ),
    sample_label =
      as.character(
        sample_label
      ),
    group =
      as.character(
        group
      )
  )

cluster_size_table <- metadata_table %>%
  dplyr::count(
    analysis_cluster,
    name = "cluster_cells"
  )

cluster_donor_composition <- metadata_table %>%
  dplyr::count(
    analysis_cluster,
    sample_label,
    group,
    name = "cells"
  ) %>%
  dplyr::left_join(
    cluster_size_table,
    by = "analysis_cluster"
  ) %>%
  dplyr::mutate(
    percent_within_cluster =
      100 * cells / cluster_cells
  ) %>%
  dplyr::arrange(
    as.integer(analysis_cluster),
    dplyr::desc(cells)
  )

cluster_group_composition <- metadata_table %>%
  dplyr::count(
    analysis_cluster,
    group,
    name = "cells"
  ) %>%
  dplyr::left_join(
    cluster_size_table,
    by = "analysis_cluster"
  ) %>%
  dplyr::mutate(
    percent_within_cluster =
      100 * cells / cluster_cells
  ) %>%
  dplyr::arrange(
    as.integer(analysis_cluster),
    group
  )

cluster_balance_summary <- cluster_donor_composition %>%
  dplyr::group_by(
    analysis_cluster
  ) %>%
  dplyr::summarise(
    cluster_cells =
      unique(cluster_cells),
    
    represented_donors =
      dplyr::n_distinct(
        sample_label[
          cells > 0
        ]
      ),
    
    maximum_single_donor_percent =
      max(
        percent_within_cluster
      ),
    
    dominant_donor =
      sample_label[
        which.max(
          percent_within_cluster
        )
      ],
    
    .groups = "drop"
  ) %>%
  dplyr::left_join(
    cluster_group_composition %>%
      dplyr::select(
        analysis_cluster,
        group,
        percent_within_cluster
      ) %>%
      tidyr::pivot_wider(
        names_from = group,
        values_from = percent_within_cluster,
        values_fill = 0,
        names_prefix = "percent_"
      ),
    by = "analysis_cluster"
  ) %>%
  dplyr::arrange(
    as.integer(
      analysis_cluster
    )
  )

cat("\nCluster balance summary:\n")
print(
  cluster_balance_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 6. Donor contribution heatmap
# ------------------------------------------------------------

donor_heatmap_data <- cluster_donor_composition %>%
  dplyr::mutate(
    analysis_cluster = factor(
      analysis_cluster,
      levels = cluster_levels
    ),
    sample_label = factor(
      sample_label,
      levels = c(
        "Normal1", "Normal2", "Normal3",
        "Normal4", "Normal5", "Normal6",
        "OA1", "OA2", "OA3",
        "OA4", "OA5", "OA6"
      )
    )
  )

donor_contribution_heatmap <- ggplot(
  donor_heatmap_data,
  aes(
    x = sample_label,
    y = analysis_cluster,
    fill = percent_within_cluster
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f",
        percent_within_cluster
      )
    ),
    size = 2.7
  ) +
  scale_fill_viridis_c(
    option = "C",
    name = "% of cluster"
  ) +
  labs(
    title =
      "Donor contribution to each integrated cluster",
    
    subtitle =
      "Values represent the percentage of each cluster contributed by a donor",
    
    x = "Cartilage donor",
    y = "Integrated cluster"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  )

group_composition_plot <- ggplot(
  cluster_group_composition %>%
    dplyr::mutate(
      analysis_cluster = factor(
        analysis_cluster,
        levels = cluster_levels
      )
    ),
  aes(
    x = analysis_cluster,
    y = percent_within_cluster,
    fill = group
  )
) +
  geom_col(
    width = 0.75
  ) +
  labs(
    title =
      "Normal and OA composition of integrated clusters",
    
    subtitle =
      "Descriptive composition only; this is not a differential abundance test",
    
    x = "Integrated cluster",
    y = "Cells within cluster (%)",
    fill = "Group"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

# ------------------------------------------------------------
# 7. Positive cluster-marker discovery
#
# Marker discovery is performed on log-normalized RNA data.
# Downsampling limits memory use and prevents very large
# clusters from completely dominating computation time.
# ------------------------------------------------------------

message("")
message("============================================")
message("Mencari positive RNA markers untuk 9 cluster...")
message("============================================")

set.seed(20260723)

cluster_markers <- Seurat::FindAllMarkers(
  object = cartilage_object,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  test.use = "wilcox",
  min.pct = 0.10,
  logfc.threshold = 0.25,
  max.cells.per.ident = 3000,
  random.seed = 20260723,
  return.thresh = 0.05,
  densify = FALSE,
  verbose = TRUE
)

if (nrow(cluster_markers) == 0) {
  stop(
    "FindAllMarkers tidak menghasilkan marker."
  )
}

cluster_markers <- cluster_markers %>%
  dplyr::mutate(
    cluster =
      as.character(
        cluster
      )
  ) %>%
  dplyr::arrange(
    as.integer(cluster),
    p_val_adj,
    dplyr::desc(
      avg_log2FC
    )
  )

top_markers_per_cluster <- cluster_markers %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 >= 0.10
  ) %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 20,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(
    as.integer(cluster),
    dplyr::desc(avg_log2FC)
  )

top10_markers_per_cluster <- top_markers_per_cluster %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::slice_head(
    n = 10
  ) %>%
  dplyr::ungroup()

cat("\nTop markers per cluster:\n")
print(
  top10_markers_per_cluster %>%
    dplyr::select(
      cluster,
      gene,
      avg_log2FC,
      pct.1,
      pct.2,
      p_val_adj
    ),
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 8. Curated cartilage marker panel
#
# These signatures are used as annotation evidence.
# Final labels must be based on several markers, not one gene.
# ------------------------------------------------------------

marker_categories <- list(
  
  Chondrocyte_core = c(
    "COL2A1",
    "ACAN",
    "SOX9",
    "COMP",
    "COL9A1",
    "COL9A2",
    "COL11A2",
    "MATN3"
  ),
  
  Regulatory_chondrocyte = c(
    "CHI3L1",
    "CHI3L2",
    "CLU",
    "CFH",
    "SERPINA3"
  ),
  
  Effector_chondrocyte = c(
    "FGFBP2",
    "CLEC3A",
    "CYTL1",
    "FRZB",
    "CHAD"
  ),
  
  Pre_fibrocartilage = c(
    "PRG4",
    "SEMA3A",
    "CRTAC1",
    "TNC"
  ),
  
  Fibrocartilage = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "COL5A1",
    "COL5A2",
    "S100A4",
    "IGFBP5",
    "MSMP",
    "DCN",
    "LUM"
  ),
  
  Prehypertrophic_hypertrophic = c(
    "COL10A1",
    "SPP1",
    "IBSP",
    "ALPL",
    "MMP13",
    "RUNX2",
    "IHH"
  ),
  
  Reparative_ECM = c(
    "COL2A1",
    "ACAN",
    "COMP",
    "FN1",
    "TGFBI",
    "CTHRC1",
    "COL11A1"
  ),
  
  Homeostatic_stress_response = c(
    "FOS",
    "JUN",
    "JUNB",
    "FOSB",
    "EGR1",
    "DUSP1",
    "HSPA1A",
    "HSPA1B"
  ),
  
  Proliferating = c(
    "MKI67",
    "TOP2A",
    "CENPF",
    "UBE2C",
    "TYMS",
    "STMN1"
  ),
  
  Metallothionein = c(
    "MT1E",
    "MT1F",
    "MT1G",
    "MT1H",
    "MT1X",
    "MT2A"
  ),
  
  Pathogenic_fibrotic_senescent = c(
    "FAP",
    "ZEB1",
    "TNC",
    "TGFBI",
    "SERPINE1",
    "COL1A1",
    "COL3A1",
    "MMP9",
    "CDKN1A",
    "CDKN2A"
  ),
  
  Immune = c(
    "PTPRC",
    "LST1",
    "TYROBP",
    "FCER1G",
    "CD68",
    "AIF1"
  ),
  
  Endothelial = c(
    "PECAM1",
    "VWF",
    "EMCN",
    "KDR",
    "EGFL7",
    "CLDN5"
  )
)

marker_annotation_table <- dplyr::bind_rows(
  lapply(
    names(marker_categories),
    function(current_category) {
      tibble::tibble(
        category =
          current_category,
        
        gene =
          marker_categories[[current_category]]
      )
    }
  )
) %>%
  dplyr::mutate(
    present_in_dataset =
      gene %in%
      rownames(
        cartilage_object
      )
  )

available_marker_table <- marker_annotation_table %>%
  dplyr::filter(
    present_in_dataset
  )

missing_marker_table <- marker_annotation_table %>%
  dplyr::filter(
    !present_in_dataset
  )

cat("\nCurated markers available in dataset:\n")
print(
  available_marker_table,
  n = Inf,
  width = Inf
)

if (nrow(missing_marker_table) > 0) {
  cat("\nCurated markers missing from dataset:\n")
  print(
    missing_marker_table,
    n = Inf,
    width = Inf
  )
}

# Keep a biologically informative but readable DotPlot.
dotplot_gene_order <- unique(
  c(
    "COL2A1", "ACAN", "SOX9", "COMP",
    "CHI3L1", "CHI3L2", "CLU",
    "FGFBP2", "CLEC3A", "CYTL1",
    "PRG4", "SEMA3A",
    "COL1A1", "COL1A2", "S100A4", "IGFBP5",
    "COL10A1", "SPP1", "IBSP", "MMP13",
    "FAP", "ZEB1", "TNC", "TGFBI",
    "MKI67", "TOP2A",
    "MT1G", "MT2A",
    "PTPRC", "LST1",
    "PECAM1", "VWF"
  )
)

dotplot_genes <- dotplot_gene_order[
  dotplot_gene_order %in%
    rownames(
      cartilage_object
    )
]

if (length(dotplot_genes) < 10) {
  stop(
    "Terlalu sedikit curated markers tersedia."
  )
}

marker_dotplot <- Seurat::DotPlot(
  object = cartilage_object,
  features = dotplot_genes,
  assay = "RNA",
  group.by = "analysis_cluster",
  scale = TRUE,
  dot.scale = 6
) +
  coord_flip() +
  labs(
    title =
      "Established cartilage and contaminating-cell markers",
    
    subtitle =
      "Cluster identities must be assigned from combined marker patterns",
    
    x = "Marker gene",
    y = "Integrated cluster",
    color = "Scaled average\nexpression",
    size = "Cells expressing\nmarker (%)"
  ) +
  theme_bw(
    base_size = 10
  ) +
  theme(
    panel.grid.major =
      element_line(
        linewidth = 0.2
      ),
    
    panel.grid.minor =
      element_blank(),
    
    axis.text.y =
      element_text(
        size = 8
      )
  )

# ------------------------------------------------------------
# 9. Heatmap of top cluster markers
# ------------------------------------------------------------

heatmap_marker_table <- cluster_markers %>%
  dplyr::filter(
    p_val_adj < 0.05,
    avg_log2FC > 0.25
  ) %>%
  dplyr::group_by(
    cluster
  ) %>%
  dplyr::slice_max(
    order_by = avg_log2FC,
    n = 5,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

heatmap_genes <- unique(
  heatmap_marker_table$gene
)

heatmap_genes <- heatmap_genes[
  heatmap_genes %in%
    rownames(
      cartilage_object
    )
]

if (length(heatmap_genes) >= 5) {
  
  message(
    "Scaling top cluster-marker genes for heatmap..."
  )
  
  cartilage_object <- ScaleData(
    cartilage_object,
    assay = "RNA",
    features = heatmap_genes,
    block.size = 500,
    verbose = FALSE
  )
  
  set.seed(20260723)
  
  heatmap_cells <- unlist(
    lapply(
      cluster_levels,
      function(current_cluster) {
        
        current_cells <- WhichCells(
          cartilage_object,
          idents = current_cluster
        )
        
        if (length(current_cells) > 120) {
          sample(
            current_cells,
            size = 120,
            replace = FALSE
          )
        } else {
          current_cells
        }
      }
    ),
    use.names = FALSE
  )
  
  top_marker_heatmap <- Seurat::DoHeatmap(
    object = cartilage_object,
    features = heatmap_genes,
    cells = heatmap_cells,
    assay = "RNA",
    group.by = "analysis_cluster",
    raster = TRUE,
    size = 3
  ) +
    labs(
      title =
        "Top positive RNA markers across integrated clusters",
      
      subtitle =
        "Up to 120 cells were sampled per cluster for visualization"
    )
  
} else {
  
  top_marker_heatmap <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label =
        "Insufficient significant markers for heatmap"
    ) +
    theme_void()
}

# ------------------------------------------------------------
# 10. MMP13 localization by cluster and condition
# ------------------------------------------------------------

mmp13_candidates <- rownames(
  cartilage_object
)[
  toupper(
    rownames(
      cartilage_object
    )
  ) == "MMP13"
]

if (length(mmp13_candidates) == 0) {
  stop(
    "MMP13 tidak ditemukan pada RNA assay."
  )
}

mmp13_feature <- mmp13_candidates[1]

rna_counts <- SeuratObject::LayerData(
  cartilage_object,
  assay = "RNA",
  layer = "counts"
)

mmp13_counts <- as.numeric(
  rna_counts[
    mmp13_feature,
    ,
    drop = TRUE
  ]
)

cartilage_object$MMP13_raw_count <-
  mmp13_counts

cartilage_object$MMP13_detected <-
  mmp13_counts > 0

cartilage_object$MMP13_status <- factor(
  ifelse(
    cartilage_object$MMP13_detected,
    "MMP13-positive",
    "MMP13-negative"
  ),
  levels = c(
    "MMP13-negative",
    "MMP13-positive"
  )
)

mmp13_cell_table <- cartilage_object[[]] %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  dplyr::mutate(
    analysis_cluster =
      as.character(
        analysis_cluster
      ),
    sample_label =
      as.character(
        sample_label
      ),
    group =
      as.character(
        group
      )
  )

mmp13_cluster_summary <- mmp13_cell_table %>%
  dplyr::group_by(
    analysis_cluster
  ) %>%
  dplyr::summarise(
    total_cells =
      dplyr::n(),
    
    MMP13_positive_cells =
      sum(
        MMP13_detected
      ),
    
    MMP13_detection_percent =
      100 *
      MMP13_positive_cells /
      total_cells,
    
    MMP13_total_counts =
      sum(
        MMP13_raw_count
      ),
    
    represented_positive_donors =
      dplyr::n_distinct(
        sample_label[
          MMP13_detected
        ]
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      MMP13_positive_cells
    )
  )

mmp13_cluster_group_summary <- mmp13_cell_table %>%
  dplyr::group_by(
    analysis_cluster,
    group
  ) %>%
  dplyr::summarise(
    total_cells =
      dplyr::n(),
    
    MMP13_positive_cells =
      sum(
        MMP13_detected
      ),
    
    MMP13_detection_percent =
      100 *
      MMP13_positive_cells /
      total_cells,
    
    MMP13_total_counts =
      sum(
        MMP13_raw_count
      ),
    
    represented_positive_donors =
      dplyr::n_distinct(
        sample_label[
          MMP13_detected
        ]
      ),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    as.integer(
      analysis_cluster
    ),
    group
  )

mmp13_positive_cell_table <- mmp13_cell_table %>%
  dplyr::filter(
    MMP13_detected
  ) %>%
  dplyr::select(
    cell_id,
    sample_label,
    donor,
    geo_accession,
    group,
    analysis_cluster,
    MMP13_raw_count,
    nCount_RNA,
    nFeature_RNA,
    percent.mt
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      MMP13_raw_count
    )
  )

cat("\nMMP13 localization by cluster:\n")
print(
  mmp13_cluster_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# 11. High-visibility MMP13 UMAP
# ------------------------------------------------------------

umap_coordinates <- as.data.frame(
  Embeddings(
    cartilage_object,
    reduction = "umap.rpca"
  )
) %>%
  tibble::rownames_to_column(
    "cell_id"
  ) %>%
  dplyr::left_join(
    mmp13_cell_table %>%
      dplyr::select(
        cell_id,
        MMP13_detected,
        MMP13_raw_count,
        analysis_cluster,
        sample_label,
        group
      ),
    by = "cell_id"
  )

umap_column_names <- colnames(
  Embeddings(
    cartilage_object,
    reduction = "umap.rpca"
  )
)

umap_x <- umap_column_names[1]
umap_y <- umap_column_names[2]

mmp13_umap <- ggplot() +
  geom_point(
    data =
      umap_coordinates %>%
      dplyr::filter(
        !MMP13_detected
      ),
    aes(
      x = .data[[umap_x]],
      y = .data[[umap_y]]
    ),
    size = 0.10,
    alpha = 0.35,
    color = "grey75"
  ) +
  geom_point(
    data =
      umap_coordinates %>%
      dplyr::filter(
        MMP13_detected
      ),
    aes(
      x = .data[[umap_x]],
      y = .data[[umap_y]],
      size = MMP13_raw_count
    ),
    alpha = 0.95,
    color = "red"
  ) +
  scale_size_continuous(
    range = c(0.9, 3.2),
    name = "MMP13\nraw count"
  ) +
  labs(
    title =
      "MMP13-positive cells on the integrated cartilage landscape",
    
    subtitle =
      paste0(
        sum(
          mmp13_cell_table$MMP13_detected
        ),
        " positive singlets among ",
        format(
          nrow(mmp13_cell_table),
          big.mark = ","
        ),
        " total cells"
      ),
    
    x = "UMAP-RPCA 1",
    y = "UMAP-RPCA 2"
  ) +
  coord_equal() +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid =
      element_blank(),
    
    legend.position =
      "right"
  )

mmp13_cluster_plot <- mmp13_cluster_group_summary %>%
  dplyr::mutate(
    analysis_cluster =
      factor(
        analysis_cluster,
        levels =
          cluster_levels
      )
  ) %>%
  ggplot(
    aes(
      x = analysis_cluster,
      y = MMP13_detection_percent,
      fill = group
    )
  ) +
  geom_col(
    position =
      position_dodge(
        width = 0.8
      ),
    width = 0.72
  ) +
  labs(
    title =
      "MMP13 detection by integrated cluster and disease group",
    
    subtitle =
      "Percentages are descriptive because expression is sparse and donor-dependent",
    
    x = "Integrated cluster",
    y = "MMP13-positive cells (%)",
    fill = "Group"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "top"
  )

mmp13_combined_panel <- (
  mmp13_umap +
    mmp13_cluster_plot
) +
  patchwork::plot_annotation(
    title =
      "Single-cell localization of MMP13 before cell-type annotation"
  )

# ------------------------------------------------------------
# 12. Save figures
# ------------------------------------------------------------

figure_list <- list(
  
  FigureS14A_cluster_donor_contribution =
    donor_contribution_heatmap,
  
  FigureS14B_cluster_group_composition =
    group_composition_plot,
  
  FigureS14C_curated_cartilage_marker_dotplot =
    marker_dotplot,
  
  FigureS14D_top_cluster_marker_heatmap =
    top_marker_heatmap,
  
  FigureS14E_MMP13_cluster_localization =
    mmp13_combined_panel
)

figure_dimensions <- list(
  c(11, 6.5),
  c(8, 5.5),
  c(10, 9),
  c(13, 8),
  c(14, 6.5)
)

for (
  current_index in
  seq_along(
    figure_list
  )
) {
  
  current_name <-
    names(
      figure_list
    )[current_index]
  
  current_plot <-
    figure_list[[current_index]]
  
  current_dimensions <-
    figure_dimensions[[current_index]]
  
  ggsave(
    filename = file.path(
      supplementary_figure_folder,
      paste0(
        current_name,
        ".pdf"
      )
    ),
    plot = current_plot,
    width = current_dimensions[1],
    height = current_dimensions[2],
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(
      supplementary_figure_folder,
      paste0(
        current_name,
        ".tiff"
      )
    ),
    plot = current_plot,
    width = current_dimensions[1],
    height = current_dimensions[2],
    dpi = 600,
    compression = "lzw",
    limitsize = FALSE
  )
}

# ------------------------------------------------------------
# 13. Save tables
# ------------------------------------------------------------

utils::write.csv(
  cluster_markers,
  file = file.path(
    table_folder,
    "GSE220243_cluster_markers_full.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  top_markers_per_cluster,
  file = file.path(
    table_folder,
    "GSE220243_top20_markers_per_cluster.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  cluster_balance_summary,
  file = file.path(
    table_folder,
    "GSE220243_cluster_balance_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mmp13_cluster_summary,
  file = file.path(
    table_folder,
    "GSE220243_MMP13_cluster_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mmp13_cluster_group_summary,
  file = file.path(
    table_folder,
    "GSE220243_MMP13_cluster_group_summary.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  mmp13_positive_cell_table,
  file = file.path(
    table_folder,
    "GSE220243_MMP13_positive_cells.csv"
  ),
  row.names = FALSE
)

summary_workbook <- file.path(
  table_folder,
  "GSE220243_cluster_marker_annotation_evidence.xlsx"
)

openxlsx::write.xlsx(
  list(
    
    Cluster_balance =
      as.data.frame(
        cluster_balance_summary
      ),
    
    Donor_composition =
      as.data.frame(
        cluster_donor_composition
      ),
    
    Group_composition =
      as.data.frame(
        cluster_group_composition
      ),
    
    Top20_markers =
      as.data.frame(
        top_markers_per_cluster
      ),
    
    Top10_markers =
      as.data.frame(
        top10_markers_per_cluster
      ),
    
    Marker_panel =
      as.data.frame(
        marker_annotation_table
      ),
    
    MMP13_by_cluster =
      as.data.frame(
        mmp13_cluster_summary
      ),
    
    MMP13_cluster_group =
      as.data.frame(
        mmp13_cluster_group_summary
      ),
    
    MMP13_positive_cells =
      as.data.frame(
        mmp13_positive_cell_table
      ),
    
    Package_versions =
      as.data.frame(
        package_versions
      )
  ),
  file = summary_workbook,
  overwrite = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    table_folder,
    "GSE220243_cluster_marker_discovery_sessionInfo.txt"
  )
)

# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

message("")
message("================================================")
message("SCRIPT 18A SELESAI")
message("================================================")

message(
  "Total cells                 : ",
  format(
    ncol(
      cartilage_object
    ),
    big.mark = ","
  )
)

message(
  "Integrated clusters         : ",
  length(
    cluster_levels
  )
)

message(
  "Positive marker rows        : ",
  format(
    nrow(
      cluster_markers
    ),
    big.mark = ","
  )
)

message(
  "MMP13-positive cells        : ",
  sum(
    cartilage_object$
      MMP13_detected
  )
)

message(
  "Marker evidence workbook    : ",
  summary_workbook
)

message(
  "Supplementary figures       : ",
  supplementary_figure_folder
)

message(
  "Final cell-type annotation  : FALSE"
)

message(
  "Next step                   : manual review and script 18B"
)

message("================================================")