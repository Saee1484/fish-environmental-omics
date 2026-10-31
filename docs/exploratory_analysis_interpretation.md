# Exploratory Analysis and Interpretation — GSE116687

## 1. Dataset overview

This exploratory analysis uses the public GEO dataset **GSE116687**, which investigates the effects of benzo[a]pyrene (BaP) exposure in adult male Nile tilapia (*Oreochromis niloticus*).

Two tissues were analysed separately:

- **Liver:** 4 control samples and 4 BaP-exposed samples.
- **Testis:** 3 control samples and 3 BaP-exposed samples.

The analysis focused on expression-profile distributions, zero-expression patterns, sample-to-sample Spearman correlations, and principal component analysis (PCA).

The purpose was to assess data structure, sample similarity, and potential patterns before undertaking further differential expression analysis.

## 2. Methods

Exploratory analysis was performed separately for liver and testis.

The following analyses were conducted:

1. Sample-level expression overview.
2. Expression density distributions.
3. Percentage of zero-expression values.
4. Pairwise Spearman correlation between samples.
5. PCA to examine sample-level variation and clustering.

Figures were exported in PNG and PDF formats. Summary tables, PCA scores, variance explained, and correlation matrices were saved as CSV files.

## 3. Liver analysis

### 3.1 PCA findings

The first two principal components explained:

- PC1: 37.50% of the variance.
- PC2: 17.67% of the variance.
- Combined: 55.17% of the variance.

The PCA plot showed a possible separation tendency along PC1, with several control samples positioned toward the negative side and several BaP-exposed samples toward the positive side.

However, this pattern was not consistent across all replicates. In particular, 18L.CTRL appeared closer to the BaP samples along PC1, while 16L.BaP occupied a relatively distinct position along PC2.

These observations suggest sample-level variation and a possible treatment-associated pattern, but PCA alone cannot establish a BaP effect.

### 3.2 Spearman correlation findings

Pairwise Spearman correlations among liver samples ranged from approximately **0.858 to 0.940**, excluding self-correlations.

Notable observations:

- The highest off-diagonal correlation was between 14L.BaP and 15L.BaP (approximately 0.940).
- The lowest correlation was between 19L.CTRL and 14L.BaP (approximately 0.858).
- 18L.CTRL showed relatively strong correlations with BaP samples (approximately 0.907–0.933).
- 16L.BaP also showed strong correlations with the other samples (approximately 0.880–0.933).

Overall, the liver samples retained substantial expression-profile similarity. The correlation results do not independently identify a sample as an outlier.

## 4. Testis analysis

### 4.1 PCA findings

The first two principal components explained:

- PC1: 27.46% of the variance.
- PC2: 23.44% of the variance.
- Combined: 50.90% of the variance.

The PCA plot did not show a clear separation between control and BaP-exposed samples. Some samples occupied relatively distinct positions, including 14G.BaP and 18G.CTRL.

This indicates variation among samples, but the observed PCA pattern does not demonstrate a distinct treatment-associated separation.

### 4.2 Spearman correlation findings

Pairwise Spearman correlations among testis samples ranged from approximately **0.938 to 0.970**, excluding self-correlations.

The highest off-diagonal correlation was between 19G.CTRL and 15G.BaP (approximately 0.970). The lowest was between 18G.CTRL and 20G.CTRL (approximately 0.938).

All samples showed strong expression-profile similarity. The correlation matrix did not reveal an obvious sample with broadly reduced similarity to the others.

## 5. Comparative interpretation

The liver and testis datasets showed different exploratory patterns.

- Liver samples showed a possible treatment-associated tendency along PC1, although replicate-level variation remained evident.
- Testis samples did not show clear separation between control and BaP-exposed groups in the first two principal components.
- Testis samples had higher pairwise Spearman correlations overall than liver samples.

These differences may reflect tissue-specific expression patterns and variability. However, the exploratory results alone cannot establish tissue-specific biological responses to BaP.

## 6. Limitations

Several limitations should be considered:

1. The number of biological replicates is relatively small, particularly for testis (three samples per group).
2. PCA is an exploratory method and does not test statistical significance of treatment effects.
3. High sample correlations do not prove that treatment has no biological effect.
4. A sample's position in PCA does not, by itself, establish that it is an outlier.
5. The expression matrix scale and data-processing method must be verified before selecting a differential expression method. The matrices contain non-integer values, so they should not automatically be treated as raw integer counts.
6. No matched environmental DNA (eDNA) community measurements are included in this dataset.

## 7. Next steps

The next stages of the project should include:

1. Verify the expression data scale and original preprocessing method using GEO documentation and supplementary files.
2. Review the published differential expression results and compare their reported findings with the exploratory patterns.
3. Select an appropriate statistical workflow based on the confirmed data type.
4. Investigate biologically relevant contaminant-response pathways and gene sets.
5. Preserve all scripts, figures, tables, metadata, and analysis decisions for reproducibility.
6. In a future independent study, consider integrating fish transcriptomic responses with eDNA-based ecological community monitoring.

## 8. Overall conclusion

The exploratory analysis identified substantial sample-to-sample expression-profile similarity in both tissues, alongside differences in their PCA patterns. Liver showed a possible treatment-associated tendency that was not consistent across every replicate, whereas testis did not show clear separation between treatment groups.

These findings provide a foundation for further quality assessment and biologically informed analysis. They should be interpreted as exploratory observations rather than evidence of statistically significant BaP-associated changes.

---

*Project: fish-environmental-omics*  
*Dataset: GEO GSE116687*  
*Analysis stage: Exploratory analysis*
