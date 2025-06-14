# Assessment of River Water Quality Using Statistical Analysis and Data Imputation in R

## Table of Contents

- [Project Overview](#project-overview)
- [Getting Started](#getting-started)
- [Workflow Steps](#workflow-steps)
- [Outputs](#outputs)
- [Conclusion](#conclusion)

---

# Project Overview
This R script performs a comprehensive water quality assessment by analyzing various parameters such as electrical conductivity, pH, hardness, calcium content, BOD, and COD. The aim is to measure these parameters at specific sampling points and assess the influence of effluent water on river quality. The workflow includes data cleaning, missing value imputation, statistical analysis (ANOVA, PCA, clustering), and visualization.

---

# Getting Started
Required Libraries:
The script uses several R packages: boot, MASS, rcompanion, ggplot2, gridExtra, and FSA for data analysis and plotting.
Data Files:
Primary data is read from "Water quality analysis (ANOVA).csv", with additional files such as "PCA analysis.csv" and "Cluster analysis with imputing data.csv" used for advanced analysis.
Working Directory:
The working directory must be set to the location where these CSV files are stored.
Workflow Steps
1. Data Cleaning and Preparation
Clears the R environment and closes open graphics devices.
Reads the main water quality dataset and prepares it for analysis by removing unused columns.
2. Handling Missing Values
For each parameter (e.g., Na, Ammonia, Calcium, Magnesium, Potassium):
Identifies correlated variables.
Fits a linear regression model to predict missing values based on related variables.
Imputes missing values using the resulting regression equation.
3. Data Transformation and Normality Check
For each parameter, the script:
Summarizes the data and checks for normality using Shapiro-Wilk tests.
Applies suitable transformations (log or Tukey) if data deviate from normality.
Visualizes distributions with density plots and histograms.
4. Visualization
Creates boxplots for each parameter grouped by treatment to visualize differences across sampling points.
5. Statistical Analysis
ANOVA and Post Hoc Testing:
Performs one-way ANOVA and Tukey's HSD post hoc tests for each parameter to identify significant differences across groups.
Principal Component Analysis (PCA):
Reduces dimensionality and finds patterns in the data.
Visualizes results with scree plots, biplots, and component plots.
Examines variable loadings and the proportion of explained variance.
Cluster Analysis:
Normalizes data and computes Euclidean distances.
Performs hierarchical clustering (complete and average linkage), generates dendrograms, and silhouette plots.
Conducts K-means clustering and visualizes cluster membership.
Outputs
Imputed Datasets:
Data with missing values filled using regression-based imputation.
Plots:
Density and histogram plots for normality checks.
Boxplots for group comparisons.
PCA scree plots, biplots, and component loading plots.
Dendrograms and silhouette plots for clustering.
K-means cluster plots.
Statistical Results:
ANOVA summaries, Tukey post-hoc results, and t-tests for each parameter.
PCA summaries and loadings.
Cluster membership tables and means.
Conclusion
This R code provides a robust workflow for water quality analysis, starting from data cleaning and imputation, proceeding through statistical analyses (ANOVA, PCA, clustering), and ending with comprehensive visualizations. The approach enables a detailed understanding of how effluent water affects river quality, highlighting both significant differences and underlying patterns across sampling locations.

