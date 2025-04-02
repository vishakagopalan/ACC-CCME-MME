Code and notebooks for generating figures in our analysis of cancer-conditioned and metastatic adrenocortical carcinoma microenvironments.

`rename_genes.R` : Functions to change gene names of all publicly available single-cell gene expression matrices to match gene names used in patient data generated from adrenocortical carcinoma (ACC) patients.

`deg_gsea_analysis.R` : MAST-based functions to compute DE genes across conditions for each cluster in combined CCME (cancer-conditioned microenvironment), MME (metastatic microenvironment) and HME (healthy microenvironment)).

`compute_gene_set_scores.R`, `gene_set_scoring.R` : Functions to use AUCell to compute gene set scores in single-cell RNA-seq datasets and normalization routines to compute gene set scores in bulk RNA-seq cohorts.

`process_data.R` : Code to process TCGA datasets.

`Compartment_Annotation.ipynb` : Separating single-cell RNA-seq data from HME, MME and CCME cohorts into individual compartments (myeloid, lymphoid, endothelial and fibroblasts).

`inferCNV, PMN, CCME classifier.ipynb` : Routines to separate malignant cells from non-malignant cells, comparing CCME vs HME signatures with pre-metastatic niche datasets and constructing and testing our CCME vs HME classifier. 

`Endothelial.ipynb`, `Lymphoid.ipynb`, `Myeloid.ipynb`, `Fibroblast.ipynb` : Code to sub-cluster and annotate cells in endothelial, myeloid, lymphoid and fibroblast compartments.

`Public Normal Data Processing.ipynb` : Code to process single-cell RNA-seq of HME datasets of lung and liver from publicly available cohorts.

`Patient Data Processing.ipynb` : Routines to process patient single-cell RNA-seq data from our ACC cohort. 

`CellChat.ipynb` : Cell-cell interaction analysis code using the Cellchat package in R.

Data associated with this code is available at https://zenodo.org/records/15119991 . Files will be made accessible to the public upon publication of our manuscript. For questions, please contact Etan Aber (etan.aber@nih.gov), Vishaka Gopalan (vishaka.gopalan@nih.gov), Sridhar Hannenhalli (sridhar.hannenhalli@nih.gov) and
Rosandra Kaplan (rosandra.kaplan@nih.gov). 
