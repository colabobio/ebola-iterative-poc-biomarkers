# Ebola Prognostication Machine Learning Model Using Iterative Point-of-Care Biomarkers

This repository contains the data and code for the manuscript titled "Derivation and Internal Validation of a Mortality Prognostication Machine Learning Model in Ebola Virus Disease Using Iterative Point-of-Care Biomarkers" by Courtney Bearnot et al.

## Data

The data include a cohort study of patients 18 years and older admitted to International Medical Corps’ Mangina ETC with laboratory-confirmed EVD from December 2018 – January 2020. For each patient, measures of ALT, albumin, amylase, AST, BUN, calcium, CRP, creatinine, creatinine kinase (CK), glucose, potassium, sodium, and total bilirubin were obtiend daily during Point-of-Care (POC) serum biomarker testing between admission and discharge.

## Code

The code is organized in two R notebooks:

* 1-data-cleaning.Rmd: It cleans up the original data by removing unnecesary variables and aggregates biomarker results in treatment days one and two (D12), days three and four (D34) and days five and six (D56).
* 2-model-updating.Rmd: It creates four iterative prognostic logistic regression models using elastic net regularization. The base model used admission age and Ct as predictors. Ct and biomarkers collected on D12, D34, and D56, that were associated with mortality, were iteratively added to the preceding model to yield temporally dynamic mortality risk estimates.
