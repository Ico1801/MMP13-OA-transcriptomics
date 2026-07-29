# ============================================================
# PROJECT : OA MMP13 Target Discovery
# SCRIPT  : 04b_GSE114007_platform_stratified_MMP13.R
# PURPOSE : Mengestimasi efek OA terhadap MMP13 secara terpisah
#           pada GPL11154 dan GPL18573
# ============================================================

rm(list = ls())

library(edgeR)
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(openxlsx)

# ------------------------------------------------------------
# 1. Memuat data
# ------------------------------------------------------------

dge_normalized <- readRDS(
  "data_processed/GSE114007_edgeR_TMM_normalized_DGEList.rds"
)

sample_metadata <- readRDS(
  "data_processed/GSE114007_sample_metadata_with_platform.rds"
)

sample_metadata$group <- factor(
  sample_metadata$group,
  levels = c("Control", "OA")
)

sample_metadata$platform_id <- factor(
  sample_metadata$platform_id,
  levels = c("GPL11154", "GPL18573")
)

stopifnot(
  identical(
    sample_metadata$sample_id,
    colnames(dge_normalized)
  )
)

stopifnot(
  "MMP13" %in% rownames(dge_normalized)
)

# ------------------------------------------------------------
# 2. Membentuk empat kombinasi platform dan kelompok
# ------------------------------------------------------------

sample_metadata <- sample_metadata %>%
  mutate(
    platform_group = paste(
      platform_id,
      group,
      sep = "_"
    )
  )

sample_metadata$platform_group <- factor(
  sample_metadata$platform_group,
  levels = c(
    "GPL11154_Control",
    "GPL11154_OA",
    "GPL18573_Control",
    "GPL18573_OA"
  )
)

cat("\nJumlah sampel setiap kombinasi:\n")

print(
  table(sample_metadata$platform_group)
)
# ------------------------------------------------------------
# 3. Design matrix cell-means
# ------------------------------------------------------------

design_stratified <- model.matrix(
  ~ 0 + platform_group,
  data = sample_metadata
)

colnames(design_stratified) <-
  levels(sample_metadata$platform_group)

rownames(design_stratified) <-
  sample_metadata$sample_id

cat("\nNama kolom design:\n")
print(colnames(design_stratified))

cat("\nRank design:\n")
print(qr(design_stratified)$rank)

stopifnot(
  qr(design_stratified)$rank ==
    ncol(design_stratified)
)

contrast_matrix <- limma::makeContrasts(
  
  GPL11154_OA_vs_Control =
    GPL11154_OA -
    GPL11154_Control,
  
  GPL18573_OA_vs_Control =
    GPL18573_OA -
    GPL18573_Control,
  
  levels = design_stratified
)

print(contrast_matrix)
# ------------------------------------------------------------
# 4. Voom model
# ------------------------------------------------------------

voom_stratified_fit <- edgeR::voomLmFit(
  counts = dge_normalized,
  design = design_stratified,
  sample.weights = FALSE,
  plot = FALSE,
  keep.EList = TRUE
)

voom_stratified_fit <- limma::contrasts.fit(
  voom_stratified_fit,
  contrasts = contrast_matrix
)

voom_stratified_fit <- limma::eBayes(
  voom_stratified_fit,
  robust = TRUE
)

voom_se_matrix <- sweep(
  voom_stratified_fit$stdev.unscaled,
  MARGIN = 1,
  STATS = sqrt(
    voom_stratified_fit$s2.post
  ),
  FUN = "*"
)

mmp13_index <- match(
  "MMP13",
  rownames(voom_stratified_fit$coefficients)
)

stopifnot(
  !is.na(mmp13_index)
)
# ------------------------------------------------------------
# 5. Fungsi mengekstrak effect size
# ------------------------------------------------------------

extract_mmp13_voom <- function(
    coefficient_name,
    platform_name,
    n_control,
    n_oa
) {
  
  coefficient_number <- match(
    coefficient_name,
    colnames(
      voom_stratified_fit$coefficients
    )
  )
  
  if (is.na(coefficient_number)) {
    stop(
      "Koefisien tidak ditemukan: ",
      coefficient_name
    )
  }
  
  effect <- as.numeric(
    voom_stratified_fit$coefficients[
      mmp13_index,
      coefficient_number
    ]
  )
  
  standard_error <- as.numeric(
    voom_se_matrix[
      mmp13_index,
      coefficient_number
    ]
  )
  
  df_vector <- voom_stratified_fit$df.total
  
  degrees_of_freedom <- if (
    length(df_vector) == 1
  ) {
    as.numeric(df_vector)
  } else {
    as.numeric(df_vector[mmp13_index])
  }
  
  critical_value <- qt(
    0.975,
    df = degrees_of_freedom
  )
  
  statistics_table <- limma::topTable(
    voom_stratified_fit,
    coef = coefficient_number,
    number = Inf,
    sort.by = "none"
  )
  
  mmp13_statistics <-
    statistics_table[
      "MMP13",
      ,
      drop = FALSE
    ]
  
  tibble(
    dataset = paste0(
      "GSE114007-",
      platform_name
    ),
    
    platform = platform_name,
    
    comparison = "OA vs Control",
    
    n_control = n_control,
    
    n_OA = n_oa,
    
    log2FC = effect,
    
    standard_error = standard_error,
    
    CI_95_lower =
      effect -
      critical_value *
      standard_error,
    
    CI_95_upper =
      effect +
      critical_value *
      standard_error,
    
    PValue_voom =
      mmp13_statistics$P.Value,
    
    FDR_voom =
      mmp13_statistics$adj.P.Val
  )
}
# ------------------------------------------------------------
# 6. edgeR QL tests per platform
# ------------------------------------------------------------

ql_stratified_fit <- edgeR::glmQLFit(
  dge_normalized,
  design = design_stratified,
  robust = TRUE
)

extract_mmp13_edger <- function(
    contrast_name
) {
  
  contrast_number <- match(
    contrast_name,
    colnames(contrast_matrix)
  )
  
  current_test <- edgeR::glmQLFTest(
    ql_stratified_fit,
    contrast =
      contrast_matrix[
        ,
        contrast_number
      ]
  )
  
  current_table <- edgeR::topTags(
    current_test,
    n = Inf,
    sort.by = "PValue"
  )$table
  
  tibble(
    contrast = contrast_name,
    
    edgeR_log2FC =
      current_table[
        "MMP13",
        "logFC"
      ],
    
    PValue_edgeR =
      current_table[
        "MMP13",
        "PValue"
      ],
    
    FDR_edgeR =
      current_table[
        "MMP13",
        "FDR"
      ]
  )
}

mmp13_edger_platform_results <- bind_rows(
  
  extract_mmp13_edger(
    "GPL11154_OA_vs_Control"
  ),
  
  extract_mmp13_edger(
    "GPL18573_OA_vs_Control"
  )
)

mmp13_platform_results <-
  mmp13_voom_platform_results %>%
  mutate(
    contrast = c(
      "GPL11154_OA_vs_Control",
      "GPL18573_OA_vs_Control"
    )
  ) %>%
  left_join(
    mmp13_edger_platform_results,
    by = "contrast"
  )

print(
  mmp13_platform_results
)
# ------------------------------------------------------------
# 7. Grafik effect size MMP13 per platform
# ------------------------------------------------------------

mmp13_platform_effect_plot <- ggplot(
  mmp13_platform_results,
  aes(
    x = platform,
    y = log2FC
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = CI_95_lower,
      ymax = CI_95_upper
    ),
    width = 0.15,
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  coord_flip() +
  labs(
    title =
      "Platform-specific MMP13 effect in GSE114007",
    
    subtitle =
      "OA versus Control with 95% confidence intervals",
    
    x = NULL,
    
    y =
      expression(
        "MMP13 adjusted " *
          log[2] *
          " fold change"
      )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor =
      element_blank()
  )

print(
  mmp13_platform_effect_plot
)

ggsave(
  filename =
    "results/figures/GSE114007_MMP13_platform_specific_effect.pdf",
  
  plot =
    mmp13_platform_effect_plot,
  
  width = 6.5,
  height = 4
)

ggsave(
  filename =
    "results/figures/GSE114007_MMP13_platform_specific_effect.tiff",
  
  plot =
    mmp13_platform_effect_plot,
  
  width = 6.5,
  height = 4,
  dpi = 600,
  compression = "lzw"
)