# Fish Environmental Omics



## Project Overview



This project explores publicly available RNA-sequencing expression data from Nile tilapia (*Oreochromis niloticus*) exposed to benzo[a]pyrene (BaP).



The aim is to examine sample-level expression patterns and assess the data's suitability for further investigation of contaminant-associated molecular responses in fish.



## Research Questions



- How do liver and testis expression profiles vary across control and BaP-exposed fish?

- Do samples show clear patterns of similarity or separation?

- What additional analysis is needed to investigate contaminant-associated molecular responses?



## Dataset



- **GEO accession:** GSE116687

- **Species:** *Oreochromis niloticus* (Nile tilapia)

- **Exposure:** Benzo[a]pyrene (BaP)

- **Tissues:** Liver and testis

- **Data source:** [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116687)



## Current Analysis



The exploratory workflow includes:



- Expression-data diagnostics

- Expression distribution plots

- Zero-expression summaries

- Sample correlation heatmaps

- Principal component analysis (PCA)



## Repository Structure



- `02_exploratory_analysis.R` — exploratory analysis script

- `data/raw/geo_supplementary/` — downloaded supplementary expression matrices

- `results/figures/` — generated plots

- `results/tables/` — metadata and analysis tables



## Limitations



This repository currently focuses on exploratory analysis. The expression matrix scale and preprocessing must be verified before selecting an appropriate differential-expression method. No matched environmental DNA (eDNA) data are included in this dataset.



## Future Work



- Verify the expression data format and preprocessing.

- Evaluate appropriate differential-expression methods.

- Investigate biological pathways if supported by the data.

- Develop a future research framework linking fish molecular responses with environmental monitoring and eDNA.



## Reproducibility



The analysis script and supporting outputs are provided for review. Software requirements, data-processing details, and execution instructions will be documented as the workflow develops.



## Citation



Colli-Dula RC et al. (2018). See the [GSE116687 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116687) for the associated study and data.
