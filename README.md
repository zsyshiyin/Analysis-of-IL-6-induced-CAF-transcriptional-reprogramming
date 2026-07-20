# Analysis-of-IL-6-induced-CAF-transcriptional-reprogramming
This study includes transcriptomic and secretomic analyses of fibroblasts, transcriptomic analysis of tumor cells treated with fibroblast-conditioned medium, and integrative analysis code, together with analysis code for public single-cell RNA-seq, spatial transcriptomics, and bulk RNA-seq data obtained from GEO, TCGA, and 10x Genomics Datasets.
Code for: Interleukin-6-reprogrammed cancer-associated fibroblasts accumulate at the invasive front and promote breast cancer invasion.

This repository contains the custom code and analysis scripts for the manuscript titled " Interleukin-6-reprogrammed cancer-associated fibroblasts accumulate at the invasive front and promote breast cancer invasion".

Abstract
This study includes transcriptomic and secretomic analyses of fibroblasts, transcriptomic analysis of tumor cells treated with fibroblast-conditioned medium, and integrative analysis code, together with analysis code for public single-cell RNA-seq, spatial transcriptomics, and bulk RNA-seq data obtained from GEO, TCGA, and 10x Genomics Datasets. These analyses aimed to investigate the mechanisms by which IL-6-reprogrammed cancer-associated fibroblasts (CAFs) contribute to breast cancer invasion.

Repository Structure
This repository is organized according to the figures presented in the manuscript. Each folder contains all the code and input files required to generate the corresponding figure.

*   `Figure3/`: Contains all scripts and input files for generating Figure 3.
*   `Figure4/`: Contains all scripts and input files for generating Figure 4.
*   `Figure5/`: Contains all scripts and input files for generating Figure 5.
*   `FigureS6/`: And so on for all supplementary figures.
*   `README.md`: This file, providing an overview of the project and instructions for reproduction.

Important Note on Raw Data
This repository does not include the original raw large-scale datasets due to file size limitations on GitHub.
Full raw data access:
RNA-seq data were deposited in the Gene Expression Omnibus (GEO): GSE337366 (fibroblasts),and GSE336649 (CM-treated tumor cells).
The mass spectrometry proteomics data were deposited in the ProteomeXchange Consortium via PRIDE (PXD080221).
Single-cell RNA sequencing (scRNA-seq) data from human TNBC tissue samples (GSE199515) and CAF/NF transcriptomic data from breast cancer patients (GSE296349) were obtained from the GEO database (https://www.ncbi.nlm.nih.gov/geo/).
Spatial transcriptomics data (Human Breast Cancer, Block A Section 1) were downloaded from the 10x Genomics Datasets portal (https://www.10xgenomics.com/datasets) (Space Ranger 1.1.0). 
The Cancer Genome Atlas (TCGA) BRCA cohort (112 normal, 1,150 tumors) was obtained from the Genomic Data Commons (GDC) Data Portal (https://portal.gdc.cancer.gov/). No further selection criteria were applied. Generic CAF scores were computed using xCell and MCPcounter.
Transcription factor enrichment analysis was performed using the ChIP-Atlas database (https://chip-atlas.org) (ChIP: TFs and others, All cell types, threshold for significance: 50).

Dependencies and Third-Party Software
1. Scripting Environment
All custom scripts were developed and tested in:
*   Perl: Version 5.12.30-v6
*   R: Version 4.5.2
Required packages are loaded within each script. Please refer to the individual scripts for specific package requirements.
2. Third-Party Software (GSEA)
This analysis uses the Gene Set Enrichment Analysis (GSEA) software, which is not included in this repository and must be downloaded separately.
*   Version Used: 4.4.0

How to Run the Code
1. Data Preparation
Please refer to the "Important Note on Raw Data" section above for information on accessing and preparing the raw datasets.
2. Run the Analysis
All scripts are organized by figure in folders named `Figure3/`, `Figure4/`, etc. The purpose and required parameters for each script are provided as comments within the script files themselves. Please refer to the header of each script for detailed instructions on its specific function and how to run it.

License
This project is licensed under the MIT License - see the `LICENSE` file for details.

Contact
For questions or issues regarding this code, please contact:
*   Corresponding Author: Ying-Hua Jin
*   Email: yhjin@jlu.edu.cn
