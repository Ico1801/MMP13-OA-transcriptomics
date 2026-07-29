# MMP13 Dysregulation in Human Osteoarthritic Cartilage

This repository contains the R scripts, analysis-ready data, metadata,
statistical outputs, figures, and software-environment records supporting
the manuscript:

> **Context-dependent MMP13 dysregulation in human osteoarthritic cartilage:
> cross-cohort transcriptomics and single-cell localization**

## Authors and affiliations

**Ymelda Agatha Christy Manurung¹**  
**Esra Yohanna Siburian²**

¹ Department of Chemical Engineering and Biotechnology,  
National Taipei University of Technology, Taipei, Taiwan

² Department of Chemical Engineering,  
National Taiwan University of Science and Technology, Taipei, Taiwan

## Study overview

This study integrates four publicly available human cartilage bulk-transcriptomic
cohorts and one single-cell RNA-sequencing cohort:

- GSE114007
- GSE117999
- GSE57218
- GSE169077
- GSE220243

The study examines whether MMP13 transcript dysregulation in human
osteoarthritic cartilage is consistent across cohorts, reference tissues,
study designs, and cartilage cell states.

The analyses include:

- cohort-specific preprocessing and quality control;
- differential-expression analysis;
- harmonization of cohort-specific MMP13 effect estimates;
- exploratory random-effects meta-analysis;
- Hallmark and Reactome gene-set enrichment analysis;
- single-cell quality control and RPCA integration;
- cartilage cell-state annotation;
- descriptive localization of MMP13-positive cells;
- donor-level pathway sensitivity analyses.

## Repository structure

```text
OA_MMP13_Target_Discovery/
│
├── data_raw/
│   └── Raw or downloaded source files used for local analysis.
│
├── data_processed/
│   └── Cleaned, normalized, filtered, or analysis-ready data objects.
│
├── metadata/
│   └── Sample information, cohort annotations, group definitions,
│       covariates, pairing information, and cell-state labels.
│
├── results/
│   └── Differential-expression results, quality-control outputs,
│       meta-analysis results, GSEA results, figures, tables,
│       donor-level summaries, and session-information files.
│
├── scripts/
│   └── R scripts used for data preparation, statistical analysis,
│       visualization, and generation of manuscript outputs.
│
└── OA_MMP13_Target_Discovery.Rproj
    └── RStudio project file defining the repository root.
