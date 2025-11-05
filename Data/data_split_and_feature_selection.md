# Data 

This folder contains [train_test_split.ipynb]() which splits our dataset into the training and testing dataset, labeled as 
- X_train.csv
- y_train.csv
- X_test.csv
- y_test.csv

The training files will be referenced by our models for training and cross-validation before the final model is selected and run on the test set.

# Feature Selection

In choosing features for our models, a t-test was conducted via [t-test_full_dataset.ipynb](), with the top genes selected detailed in 
- top_400_genes.csv
- top_1000_genes.csv

For Independent Component Analysis (ICA), GridSearchCV() was initially used to identify the top number of components to be used for the model in [logreg_ICA.ipynb](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Models/logreg_ICA.ipynb).

Additionally, we chose 10 components in reference to Figure 6 of [INTERLIVER study update_2024.pdf](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Literature/INTERLIVER%20study%20update_2024.pdf) where they identified 9-10 terms, or groups, of genes that represent various function. 

The top 500 genes were identified from Component 10 which showed a distinction between the two groups, as seen in our [pairplot](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Data/ICA_pairplot_10.png).

For PCA, GridSearchCV() was used to identify the optimal number of components to be used, which can be found in [logreg_PCA.ipynb](https://github.com/Erdos-Projects/fall-2025-predicting-liver-graft-acceptance/blob/main/Models/logreg_PCA.ipynb).

