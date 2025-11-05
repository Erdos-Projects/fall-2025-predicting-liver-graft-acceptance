# Predicting Liver Transplant Rejection
Team Members: Gonca Bulbul, Noah Mac, Winnie Lau

This is our Fall 2025 Erdős Institute Data Science Boot Camp project. 

## Introduction 
Liver transplantation is a life-saving therapy where the diseased liver is replaced with another person's healthy liver, but one of the main complications of liver transplantation is immune-mediated rejection. This occurrs when patient's immune system  attacks the transplanted liver.

The goal of our project is to develop a predictive model that can recognise patients with a high risk to develop rejection based on gene expression data. Genomic signatures can be identified by comparing differentially expressed genes in the rejected group compared to the non-rejected group, i.e. the control group.

Our project references the paper [Defining an NK Cell–enriched Rejection-like Phenotype in Liver Transplant Biopsies From the INTERLIVER Study](https://journals.lww.com/transplantjournal/fulltext/2025/08000/defining_an_nk_cell_enriched_rejection_like.24.aspx) where we accessed their data from the Gene Expression Omnibus under the accession number [GSE277334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277334).

See more details in [Checkpoint 1](https://docs.google.com/document/d/10DiEbP9I1G3XYccKOn7K5ndFGJNbdZNPD-gfpPvCwfM/edit?usp=sharing).

## Data Analysis Workflow

### Preprocessing
Due to the large size of the raw microarray dataset containing 48744 features, preprocessing was needed to clean the data, detailed in [Checkpoint 2](https://docs.google.com/document/d/16MgFH9F_lhokG8RnbQstxtkFHi3-UE2dls6Zt_3UuN0/edit?usp=sharing). with the relevant notebooks found in [Checkpoint2_Microarray_Processing](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/tree/main/Checkpoint2_Microarray_Processing).

### Exploratory Data Analysis
The data was initially analyzed via t-test to identify genes that were significantly differentially expressed between rejection and non-rejection liver transplant biopsies. The data was then visualized using Principal Component Analysis (PCA) to illustrate the separation between the groups before and after a t-test, and a heatmap was used to reveal the distinct transcriptional profiles between rejection and non-rejection liver transplant biopsies. 

See more details in [EDA and statistical analysis_1104.ipynb](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/EDA%20and%20statistical%20analysis/EDA%20and%20statistical%20analysis_1104.ipynb).


## Modeling Workflow
Our dataset was split into a training set and testing set using [train_test_split.ipynb](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Data/train_test_split.ipynb). The files are detailed in [data_split_and_feature_selection.md](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Data/data_split_and_feature_selection.md).

To further decrease the size of our dataset, the feature selection we used was t-test, PCA, and Independent Component Analysis (ICA). The models we chose for classification were Logistic Regression, Random Forest, and XGBoost.

See more details in [Checkpoint 3](https://docs.google.com/document/d/1lXauiea1D0aFap8AOmuSvD5to4yxaVE1lUdqahYrsIk/edit?usp=sharing)

## Results
The training results for our models can be found in [Checkpoint 4.](https://docs.google.com/document/d/1Y4EvtKCyqkBQrXOt7FIfvxMKX9Ov-_0U-zRjuhxCuEU/edit?usp=sharing)

Our best model was identified to be Logistic Regression with 500 genes chosen from one component is ICA, where it had the highest accuracy and ROC-AUC score. This was fit to the test set using [final_model_test_results.ipynb](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Models/final_model_test_results.ipynb) which gave a final accuracy of 0.88889 and AUC = 0.961. Additional confusion matrix and ROC-AUC curve can be found in [Plots](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/tree/main/Models/Plots). 

