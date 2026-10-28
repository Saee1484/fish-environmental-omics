# Published Differential-Expression Results

## Data source

The results summarized here come from Supplementary Table SI-1
in the published study associated with GEO accession GSE116687.
The supplementary workbook is stored at:
`data/raw/geo_supplementary/mmc2.xlsx`.

## Summary

| Tissue | Genes with published FDR < 0.05 |
|---|---:|
| Liver | 207 |
| Testis | 32 |

These counts were calculated from the published FDR columns
in SI-Table-1. Missing FDR values were excluded from the counts.

## Interpretation

The supplementary table contains more genes with FDR values
below 0.05 for liver than for testis. This is a descriptive
comparison of the published results, not proof that one tissue
is universally more sensitive to contaminant exposure.

## Important limitations

- These counts summarize the published supplementary table;
  they are not results from a new differential-expression analysis.
- The original study's gene-selection criteria and the FDR < 0.05
  threshold used for these counts are not interchangeable.
- The expression scale of the downloaded GEO matrices must be
  verified before selecting a new differential-expression method.
- GSE116687 does not include matched eDNA community measurements.
  It therefore cannot establish whether eDNA detects contaminant
  effects earlier than fish transcriptomics.

## Reproducibility

The descriptive counts are available in
`results/tables/published_FDR_summary.csv`.
The source workbook is retained in the raw supplementary-data folder.
