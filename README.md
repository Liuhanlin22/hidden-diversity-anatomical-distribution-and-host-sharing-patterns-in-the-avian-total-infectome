# Wildbird infectome figure code and data

This directory contains the scripts, tabular inputs, quality-control summaries, and derived outputs used to generate the main and supplementary figures of the wildbird infectome study.

## Directory map

| Figure set | Main analysis | Primary script or scripts |
|---|---|---|
| Figure 1 | Sampling map, host metadata, and host phylogeny | `Figure1/` |
| Figure 4 | Novelty summaries, family counts, heatmap, and depth plots | `Figure4/` |
| Figure 5 | Tissue landscape, systemic infection heatmap, and t-SNE | `Figure5/` |
| Figure 6 | Geographic distribution, host-order network, host range, and wild-domestic transmission | `Figure6/` |
| Figure S2 | Quality-control summaries and taxonomic composition | `FigureS2/` |
| Figure S3 | Brain-associated viral read-depth and gene-annotation plots | `FigureS3/` |
| Figure S4 | Library-based rarefaction and diversity analysis | `FigureS4/` |
| Figure S5 | Sampling-effort rarefaction and matched t-SNE panels | `FigureS5/` |
| Metadata | Project, library, pathogen, and profile metadata | `metadata/` |

## Reproducibility

The repository root is this `data_Figures` directory. Run each script from its own figure directory, or execute it with `Rscript` using the script path. The current scripts use relative paths based on the script directory where supported.

Examples:

```bash
cd Figure5/figure5C
Rscript Figure5C_tSNE.R

cd ../../FigureS3
Rscript FigureS3_brain_depth_with_annotations.R

cd ../FigureS5
Rscript FigureS5_rarefaction_tSNE.R
```

The FigureS2B Python script can be run as follows:

```bash
cd FigureS2/FigureS2B
python 00_stacked_barplot.py
```

## GitHub Desktop upload

This directory can be uploaded as a local GitHub repository using GitHub Desktop:

1. Create an empty repository on GitHub. Do not initialize it with a README, `.gitignore`, or license.
2. In GitHub Desktop, select **File > Add local repository**.
3. Choose this `data_Figures` directory as the local repository. If GitHub Desktop reports that it is not a repository, choose **Create a repository** and use this directory as the local path.
4. Review the changed files, enter a commit message such as `Add reproducible figure code and metadata`, and select **Commit to main**.
5. Select **Publish repository** or **Push origin**. Use a public repository only when the sequence and metadata release has been approved.

The five final supplementary metadata files are stored in `metadata/`:

- `Supplement_Table1_project_info.xlsx`: project-level information.
- `Supplement_Table2_libs_hosts_qc_info.xlsx`: library, host, location, and QC information.
- `Supplement_Table3_library_info.xlsx`: library-level information.
- `Supplement_Table4_pathogen_annot_updated_ICTV.xlsx`: pathogen annotation, ICTV classification, host-range, reference, and phylogenetic information.
- `Supplement_Table5_pathogen_profiles.xlsx`: pathogen profile and abundance information.

The five metadata workbooks were checked for Chinese characters and none were found. Before publishing, verify that they contain no private identifiers or restricted information.

## R packages

Depending on the figure, the analyses use the following R packages:

`ggplot2`, `ggnewscale`, `gggenes`, `grid`, `gridExtra`, `igraph`, `maps`, `openxlsx`, `patchwork`, `pheatmap`, `readxl`, `Rtsne`, `scatterpie`, `tidyverse`, `vegan`, `permute`, and `circlize`.

## Python packages

The Python plotting script requires:

`pandas`, `matplotlib`, and `seaborn`.

## Input and output conventions

- Files ending in `_input`, `_annot`, `_metadata`, `_matrix`, or `_table` are analysis inputs or input annotations.
- Files ending in `_summary`, `_statistics`, `_coordinates`, `_audit`, or `_results` are derived QC or analysis outputs.
- PDF and PNG files are rendered figure outputs.
- Raw sequencing reads, BAM files, external databases, and upstream profiling tables are not duplicated in this directory. When an analysis depends on an upstream file outside this directory, the audit file records a relative source path or the source is described in the corresponding figure README.

## Statistical and ordination notes

Figure 5C uses t-SNE for visualization after `log10(RPM + 1)` transformation. PERMANOVA and PERMDISP are calculated from the Bray–Curtis distance matrix with the blocking scheme specified in the script. t-SNE does not have an NMDS stress statistic; the script therefore does not report stress. t-SNE is used as an exploratory visualization and is not treated as an independent biological validation.

Figure S5 combines library-based rarefaction with matched t-SNE panels for the wild-domestic, host-order, and geographic comparisons. The rarefaction and t-SNE panels use the corresponding processed input tables described in `FigureS5/`.

## File-language and portability audit

An audit of the files in this directory found no Chinese characters in file names, scripts, CSV/TSV/TXT/Markdown files, or Excel cell values. Machine-specific absolute paths were removed from the portable summaries where possible. The `tsne_sessionInfo.txt` file may retain the host operating system locale as part of the recorded R session information; this does not affect the analysis.

## Data availability note

This directory is prepared for code and supplementary metadata deposition. Viral contigs should be deposited separately in `sequences/` or in a public nucleotide archive, with the accession mapping recorded in a sequence manifest. Do not upload private credentials, local machine paths, raw BAM files, raw FASTQ files, or restricted databases.
