# Figure S5: rarefaction curves with matched tSNE panels

Run in R/RStudio from this directory:

```r
source("FigureS5_rarefaction_tSNE.R")
```

The script creates a six-panel figure:

- A–C: rarefaction curves for Wild/Domestic, host order, and Geography.
- D–F: the corresponding tSNE plots directly below A–C.

Required inputs:

- `FigureS5_virus_table_filtered.csv`
- `FigureS5_libs_annot.xlsx`
- `FigureS5_virus_table_selected_samples_updated.csv`
- `FigureS5_library_annotation.csv`

Main outputs:

- `FigureS5_rarefaction_tSNE_combined.pdf`
- `FigureS5_rarefaction_tSNE_combined.png`
- `FigureS5_tSNE_coordinates.csv`
- `FigureS5_tSNE_statistics.csv`
- renamed rarefaction summary tables beginning with `FigureS5_BCD_`

The cleanup script `cleanup_FigureS5_old_FigureS4_files.ps1` removes only the stale
`FigureS4*` files and `FigureS4D` subdirectory from this FigureS5 folder.
