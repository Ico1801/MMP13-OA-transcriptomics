# MMP13-OA-transcriptomics
This repository contains the R scripts, processed analysis inputs,session information, and workflow documentation supporting the manuscript
"Context-dependent MMP13 dysregulation in human osteoarthritic cartilage: cross-cohort transcriptomics and single-cell localization"
## Authors

- Ymelda Agatha Christy Manurung
- Esra Siburian

Department of Chemical Engineering and Biotechnology,
National Taipei University of Technology, Taipei, Taiwan.

## Study overview

The study integrates four human cartilage bulk-transcriptomic cohorts
(GSE114007, GSE117999, GSE57218, and GSE169077) and one single-cell
RNA-sequencing cohort (GSE220243).

The analyses include:

- cohort-specific preprocessing and differential-expression analysis;
- harmonization of MMP13 effect estimates;
- exploratory random-effects meta-analysis;
- Hallmark and Reactome gene-set enrichment analysis;
- single-cell RPCA integration and cartilage cell-state annotation;
- descriptive localization of MMP13-positive cells;
- donor-level pathway sensitivity analyses.

## Data availability

All raw source data are publicly available from the NCBI Gene Expression
Omnibus under accession numbers GSE114007, GSE117999, GSE57218,
GSE169077, and GSE220243.

Raw GEO files are not redistributed in this repository.

## Repository structure

- `scripts/`: R scripts in recommended execution order.
- `metadata/`: cohort and sample metadata.
- `derived_data/`: analysis-ready derived results.
- `sessionInfo/`: R and package-version records.
- `figures/`: main and supplementary figure outputs.

## Software

Analyses were conducted in R 4.6.0. Exact package versions are reported
in the `sessionInfo/` directory.

## Reproducibility

Run the scripts numerically from `00_setup.R` onward. File paths should
be configured in `00_setup.R`.

## Citation

Please cite the archived Zenodo release associated with this repository.

## License

See the LICENSE file.

## Contact

Ymelda Agatha Christy Manurung  
Department of Chemical Engineering and Biotechnology  
National Taipei University of Technology
