# Predicting Liver Transplant Rejection
Team Members: Gonca Bulbul, Noah Mac, Winnie Lau

This is our Fall 2025 Erdős Institute Data Science Boot Camp project. 

## Introduction 
Liver transplantation is a life-saving therapy where the diseased liver is replaced with another person's healthy liver, but one of the main complications of liver transplantation is immune-mediated rejection. This occurrs when patient's immune system  attacks the transplanted liver.

The goal of our project is to develop a predictive model that can recognise patients with a high risk to develop rejection based on gene expression data. Genomic signatures can be identified by comparing differentially expressed genes in the rejected group compared to the non-rejected group, i.e. the control group.

Our project references the paper [Defining an NK Cell–enriched Rejection-like Phenotype in Liver Transplant Biopsies From the INTERLIVER Study](https://journals.lww.com/transplantjournal/fulltext/2025/08000/defining_an_nk_cell_enriched_rejection_like.24.aspx) where we accessed their data from the Gene Expression Omnibus under the accession number [GSE277334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277334).

See more details in [Checkpoint 1.](https://docs.google.com/document/d/10DiEbP9I1G3XYccKOn7K5ndFGJNbdZNPD-gfpPvCwfM/edit?usp=sharing)

## Data Analysis Workflow
# Preprocessing
Due to the large size of the raw microarray dataset containing 48744 features, preprocessing was needed to clean the data, detailed in [Checkpoint 2](https://docs.google.com/document/d/16MgFH9F_lhokG8RnbQstxtkFHi3-UE2dls6Zt_3UuN0/edit?usp=sharing) with the relevant notebooks found in [Checkpoint2_Microarray_Processing](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/tree/main/Checkpoint2_Microarray_Processing).

# Exploratory Data Analysis
The data was visualized with a correlation heatmap and PCA (...)

[Checkpoint 3](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Models/Checkpoint3_evaluation_plan.md#models)


## Modeling (and linking to specific files/folder in the github too)

***(maybe include diagram of modeling workflow here?)***

To further decrease the size of our dataset, the feature selection we used is 
[Checkpoint 3](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Models/Checkpoint3_evaluation_plan.md#models)

The final results for our models can be found in [Checkpoint 4](https://docs.google.com/document/d/1Y4EvtKCyqkBQrXOt7FIfvxMKX9Ov-_0U-zRjuhxCuEU/edit?usp=sharing)


We could add slide deck and link here at the end too once we upload it?
[Checkpoint 5](https://docs.google.com/document/d/1uGCoXnhvnfNjrd15Ap4mVyRvNdBr_XMwIifV_6-agKs/edit?usp=sharing)

