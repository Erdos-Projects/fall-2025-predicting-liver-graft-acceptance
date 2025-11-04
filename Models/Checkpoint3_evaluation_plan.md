## **Modeling Approach**

Goal: to predict if a liver transplant will be rejected based on gene expression data. 
- Our dataset is composed of microarray data with gene expression values, where each sample is a liver biopsy from a patient. 
- We used a 80/20 train-test split for testing our models, with a 70/30 distribution of our non-rejection and rejection group, respectively, maintained in each dataset. 
- The dataset is split into X_train.csv, X_test.csv, y_train.csv, y_test.csv to prevent data leakage, and that only the training datasets are used for feature selection and model training. 
- For model cross validation, StratifiedKFold is used to ensure both groups are represented in each split due to the uneven distribution. 

### **Feature Selection**

Our dataset is of dimension (763, 18644) where the number of features far exceeds the number of samples. This makes our subsequent models prone to overfitting, so feature selection is needed to reduce the dimensions of the dataset. 

**Principal Component Analysis**
- One method to reduce the dimensions of our dataset by capturing the principal components of the dataset, i.e. the directions with maximum variance. 
- By creating a scree plot, we observe how the captured variance changes with each component, and choose the best number of components for our model. 
- A constraint is the mutual orthogonality of components, which is not always appropriate for biomedical data. PCA also only uses second-order statistics with the covariance between observed variables and may fail to extract other features that may require higher-order statistics. 
We use sklearn.decomposition with PCA. 

**t-test, p-value, and adjusted p-value**
- We have two independent groups for our samples, so a consideration is to run Welch’s t-test and compare the means between the two groups for each gene. 
- We then find the p-values and adjusted p-values with the Bonferroni correction due to the large number of genes in our dataset. Comparisons can then be made with the genes that are under the threshold of 0.05, and then taking the genes that are found to be most significant to be used for further analysis. 
- We use scipy.stats with ttest_ind and statsmodels with multipletests.

**Independent Component Analysis (ICA)**
- ICA is a method in signal processing used to separate mixed signals into additive subcomponents. It has been widely used for gesture extraction, gene clustering, and classification in microarray data analysis, making it a viable option for our dataset. 
- The goal of ICA is to find a linear representation of non-Gaussian data so that the components are statistically independent in order to capture the essential structure of the data. 
- In applying it to microarray data, it decomposes the dataset into components that are as statistically independent from the other as possible, potentially revealing the underlying gene expression patterns for gene clustering and classification.
- Theory:
  - ICA breaks down our original gene profile dataset $X$ into a matrix $A$ representing the mixing system with the latent variables and S containing the source signals of the gene signatures whose rows are mutually independent statistically independent coefficients, such that $X = A \cdot S$.
  - The matrices A and S are unknown, so we want to recover the original signals from the observations X using the linear transform 
				\[Y = WX\]
		where $Y$ is the estimation of the source signal $S$ and $W$ is the demixing matrix. 
- To implement this in Python, we use sklearn.decomposition with FastICA.


### **Models**

We have a binary classification defined as:
0: non-rejection
1: rejection
so we want to use classification models. Each model is primarily compared using the accuracy score, with the ROC-AUC score as an additional reference. 

**Logistic Regression**
- A simple supervised algorithm for classification, which fits the goal of our project. 
- Logistic regression can be run on both the entire dataset and a subset of the dataset after feature selection, PCA, and ICA has been used. 

**Random Forest**
- This is an ensemble learning method useful for classification by creating multiple decision trees to make better predictions during training. 
- Need to be wary of overfitting due to the large number of features in our dataset, where the model will fit to the additional noise.
- A large and complex model with many parameters that will have to be tuned. 

**XGBoost**
- XGBoost is a boosting model with the flexibility to work with large amounts of gene expression data by integrating multiple tree models to build a stronger learner model. 
- Also has potential to overfit to noise in large datasets with many features, along with being a large and complex model with many hyperparameters.
