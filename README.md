# SZClassifier
Classification of schizophrenia from blood methylation markers

## Classifying Schizophrenia Using DNA Methylation Markers

__Schizophrenia__ is a complex psychiatric disorder that affects about 1% of the world's population. It is characterized by a range of symptoms, such as hallucinations, delusions, cognitive impairment, and social withdrawal. The exact causes of schizophrenia are still unknown, but genetic and environmental factors are thought to play a role.

DNA methylation is a type of epigenetic modification that involves adding a methyl group to a cytosine base in the DNA sequence. DNA methylation can affect gene expression and influence various biological processes, such as development, aging, and disease. Previous studies have shown that DNA methylation patterns are altered in schizophrenia patients compared to healthy controls, suggesting that DNA methylation may be involved in the pathogenesis of schizophrenia.

This project aims to use DNA methylation markers to classify schizophrenia patients and healthy controls. We will use a dataset of DNA methylation profiles from the blood of schizophrenia patients and controls. The dataset contains approximately 450,000 CpG sites (regions of DNA where cytosine is followed by guanine) measured by the Illumina Infinium HumanMethylation450 BeadChip array.

We performed the following steps to analyze the data:
- __Data preprocessing__: We will clean and normalize the data, remove missing values and outliers, and select the most relevant CpG sites for classification.
- __Data analysis__: We will explore the data using descriptive statistics and visualization techniques, such as histograms, boxplots, scatterplots, and heatmaps.
- __Data modeling__: We will build and evaluate different machine learning models, such as logistic regression, support vector machines (SVM), k-nearest neighbors (KNN), decision trees, random forests, and neural networks. We will use cross-validation to optimize the hyperparameters and select the best model based on AUC, F1 score, and Matthews Correlation Coefficient (MCC).
- __Data interpretation__: We will interpret the results of our analysis and modeling, identify the most important CpG sites for classification, and discuss the biological implications of our findings. We will also compare our results with previous studies and suggest future directions for research.

By using DNA methylation markers to classify schizophrenia patients and healthy controls, we hope to gain more insights into the molecular mechanisms of schizophrenia and potentially identify new biomarkers for diagnosis and treatment.
