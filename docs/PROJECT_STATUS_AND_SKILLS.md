# Fish Environmental Omics

## Bioinformatics Skills Portfolio



**GitHub repository:** https://github.com/Saee1484/fish-environmental-omics



## 1. Project purpose



This repository is an independent bioinformatics skills portfolio developed using publicly available biological datasets.



Its purpose is to demonstrate practical abilities in:

- Retrieving and organising public biological datasets.

- Understanding experimental design and sample metadata.

- Handling gene-expression matrices and sequencing files.

- Performing sequencing quality control.

- Conducting exploratory data analysis using R.

- Interpreting biological results and analytical limitations.

- Maintaining reproducible workflows using Git and GitHub.



The portfolio includes fish transcriptomics analysis and preliminary sequencing quality-control work. eDNA metabarcoding will be explored separately as an additional bioinformatics skill.



This is a skills-development portfolio, not a continuation of the MSc thesis or a proposed laboratory project.



## 2. Transcriptomics dataset



**Dataset:** GSE116687  

**Organism:** Nile tilapia (*Oreochromis niloticus*)  

**Exposure:** Benzo[a]pyrene (BaP)  

**Tissues:** Liver and testis  

**Platform:** Illumina NextSeq 500  

**Reference genome reported in the study:** Orenil2.0



The dataset was selected to practise public gene-expression data handling and exploratory analysis of molecular responses to environmental contaminant exposure.



## 3. Work completed



### 3.1 Dataset retrieval and organisation



- Retrieved GEO supplementary expression matrices.

- Organised original supplementary files separately from analysis outputs.

- Prepared sample metadata for liver and testis.

- Retrieved SRA run information associated with the study.

- Created an SRA run manifest at `results/tables/GSE116687_SRA_run_manifest.csv`.



### 3.2 Exploratory expression analysis using R



Separate exploratory analyses were performed for liver and testis.



Methods included:

- Sample-level expression overview.

- Expression distribution and density plots.

- Zero-expression assessment.

- Sample-to-sample correlation analysis.

- Principal component analysis (PCA).



**Liver observations**

- Pairwise correlations ranged approximately from 0.858 to 0.940.

- PC1 explained 37.50% of variation.

- PC2 explained 17.67% of variation.

- PCA showed partial exposure-related grouping, but not definitive separation.



**Testis observations**

- Pairwise correlations ranged approximately from 0.938 to 0.970.

- PC1 explained 27.46% of variation.

- PC2 explained 23.44% of variation.

- PCA did not show clear control-versus-exposure separation.



These are exploratory observations. Correlation and PCA alone cannot establish whether BaP caused a biological effect.



### 3.3 Sequencing data handling and quality control



A preliminary pilot was performed using SRA run `SRR7474517`.



Tools used:

- SRA Toolkit

- FastQC



FastQC summary:

- Total reads: 6,682,257.

- Total bases: approximately 504.5 million.

- Mean read length: 75 bp.

- Overall GC content: approximately 47%.

- Per-base sequence quality: PASS.

- Adapter content: PASS.

- Per-base sequence content: FAIL.

- Per-sequence GC content: WARN.

- Sequence length distribution: WARN.



The QC report was saved under `results/fastqc/`.



This was a one-run pilot, not a complete quality assessment of all sequencing runs. The warnings were interpreted rather than treated as automatic reasons for trimming.



## 4. Analytical limitations and scientific reasoning



The GEO supplementary expression matrices contain decimal-valued measurements. Their exact measurement scale and processing provenance must be verified before selecting a differential-expression method.



Therefore:

- The matrices have not been assumed to contain raw integer counts.

- A new DESeq2 or edgeR analysis has not been claimed.

- Published differential-expression findings must be distinguished from independently generated results.

- The SRA pilot does not establish that all associated runs are suitable for independent remapping.



Recognising these limitations is an important part of responsible bioinformatics analysis.



## 5. Additional skill development: eDNA data handling



A separate public eDNA metabarcoding dataset will be used to practise the main stages of eDNA sequence processing.



Planned workflow:



1. Retrieve public eDNA sequencing data and sample metadata.

2. Understand sampling design, marker region, primers and sequencing layout.

3. Inspect raw FASTQ files using FastQC.

4. Perform primer and adapter removal using Cutadapt.

5. Apply quality filtering and sequence denoising using DADA2 or an appropriate alternative.

6. Identify and remove chimeric sequences.

7. Generate an amplicon sequence variant (ASV) table.

8. Assign taxonomy using a suitable reference database.

9. Assess negative and positive controls and potential contamination.

10. Explore taxonomic composition, richness, diversity and community differences.



The actual workflow will depend on the selected dataset, marker region and sequencing design. These steps are planned learning activities, not completed analyses.



## 6. Tools and computational skills



**Used in the current project**

- R and RStudio

- GEO and SRA data retrieval

- SRA Toolkit

- FastQC

- Exploratory data analysis and visualisation

- Git and GitHub

- Markdown documentation



**Planned for separate eDNA practice**

- Cutadapt

- DADA2

- Taxonomic reference databases

- ASV-based ecological analysis



## 7. Reproducibility practices



The project follows these principles:



- Preserve original raw data.

- Separate raw data, scripts, figures and tables.

- Record sample metadata and accession information.

- Document methods and software.

- Distinguish technical sequencing runs from biological replicates.

- Maintain version control using Git and GitHub.

- Report analytical limitations transparently.

- Distinguish completed work from planned work.



## 8. Current status



**Completed**

- GEO dataset retrieval and organisation.

- Liver and testis exploratory expression analysis.

- Exploratory figures and interpretation.

- Sample metadata and SRA run manifest.

- One-run SRA quality-control pilot.

- GitHub version control and repository updates.



**Further skills-development activities**

- Verify GEO expression matrix provenance.

- Investigate suitable count data for differential-expression analysis.

- Practise additional transcriptomics methods where data permit.

- Select an independent public eDNA metabarcoding dataset.

- Practise a reproducible eDNA processing workflow.



## 9. Portfolio objective



This repository demonstrates practical familiarity with biological data handling, sequencing quality assessment, R-based exploratory analysis, bioinformatics tools, scientific interpretation and reproducible computational workflows.



It is an independent portfolio built from public-data exercises. It does not claim completed eDNA research or imply that this particular dataset will be continued in another laboratory.

