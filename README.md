# Wildbird infectome figure code and data

This repository contains the scripts, tabular inputs, metadata, quality-control summaries, and derived outputs used to generate the main and supplementary figures of the wildbird infectome study.

## Directory map

| Figure set | Main analysis | Directory |
|---|---|---|
| Figure 1 | Sampling map, host metadata, and host phylogeny | `Figure1/` |
| Figure 4 | Novelty summaries, family counts, heatmap, and depth plots | `Figure4/` |
| Figure 5 | Tissue landscape, systemic infection heatmap, and t-SNE | `Figure5/` |
| Figure 6 | Geographic distribution, host-order network, host range, and wild-domestic transmission | `Figure6/` |
| Figure S2 | Quality-control summaries and taxonomic composition | `FigureS2/` |
| Figure S3 | Library-based rarefaction and diversity analysis | `FigureS3/` |
| Figure S4 | Brain-associated viral read-depth and gene-annotation plots | `FigureS4/` |
| Figure S5 | Sampling-effort rarefaction and matched t-SNE panels | `FigureS5/` |
| Metadata | Project, library, pathogen, and profile metadata | `metadata/` |

## Reproducibility

This repository contains the contents prepared in the `data_Figures` directory. Run each script from its corresponding figure directory. The scripts use relative paths based on the script directory where supported.

Examples:

```bash
cd Figure5/figure5C
Rscript Figure5C_tSNE.R

cd ../../FigureS3
Rscript FigureS3_rarefaction.R

cd ../FigureS4
Rscript FigureS4_brain_depth_with_annotations.R

cd ../FigureS5
Rscript FigureS5_rarefaction_tSNE.R
