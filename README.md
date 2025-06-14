# Assessment of River Water Quality Using Statistical Analysis and Data Imputation in R

## Table of Contents

- [Project Overview](#project-overview)
- [Features](#features)
- [Getting Started](#getting-started)
- [Workflow Steps](#workflow-steps)
- [Outputs](#outputs)
- [Conclusion](#conclusion)

---

# Project Overview
This R script performs a comprehensive water quality assessment by analyzing various parameters such as electrical conductivity, pH, hardness, calcium content, BOD, and COD. The aim is to measure these parameters at specific sampling points and assess the influence of effluent water on river quality. The workflow includes data cleaning, missing value imputation, statistical analysis (ANOVA, PCA, clustering), and visualization.

---

# Features
- Handles missing data using regression-based imputation
- Checks and transforms data for normality
- Performs statistical tests (ANOVA, Tukey post-hoc, t-tests)
- Conducts Principal Component Analysis (PCA)
- Executes cluster analysis (hierarchical and k-means)
- Automated generation of publication-ready visualizations

---

# Data Description

1. **Primary Dataset:**

   - Water quality analysis (ANOVA).csv
   - Contains raw water quality measurements for different parameters and sampling points.
   - Supplementary Datasets:

2. **PCA analysis.csv**

   - Used for principal component analysis.
   - Cluster analysis with inputting data.csv
   - Used for cluster analysis after imputation.

---

# Getting Started

1. **Required Libraries:**
The script uses several R packages: boot, MASS, rcompanion, ggplot2, gridExtra, and FSA for data analysis and plotting.

2. **Data Files:**
Primary data is read from "Water quality analysis (ANOVA).csv", with additional files such as "PCA analysis.csv" and "Cluster analysis with imputing data.csv" used for advanced analysis.

3. **Working Directory:**
The working directory must be set to the location where these CSV files are stored.

---

# Workflow Steps

1. **Data Cleaning and Preparation**
   - Clears the R environment and closes open graphics devices.
   - Reads the main water quality dataset and prepares it for analysis by removing unused columns.

2. **Handling Missing Values**
   - For each parameter (e.g., Na, Ammonia, Calcium, Magnesium, Potassium):
   - Identifies correlated variables.
   - Fits a linear regression model to predict missing values based on related variables.
   - Imputes missing values using the resulting regression equation.

3. **Data Transformation and Normality Check**

For each parameter, the script:
   - Summarizes the data and checks for normality using Shapiro-Wilk tests.
   - Applies suitable transformations (log or Tukey) if data deviate from normality.
   - Visualizes distributions with density plots and histograms.

4. **Visualization**
Creates boxplots for each parameter grouped by treatment to visualize differences across sampling points.

5. **Statistical Analysis**

- **ANOVA and Post Hoc Testing:**
  - Performs one-way ANOVA and Tukey's HSD post hoc tests for each parameter to identify significant differences across groups.

- **Principal Component Analysis (PCA):**
  - Reduces dimensionality and finds patterns in the data.
  - Visualizes results with scree plots, biplots, and component plots.
  - Examines variable loadings and the proportion of explained variance.

- **Cluster Analysis:**
  - Normalizes data and computes Euclidean distances.
  - Performs hierarchical clustering (complete and average linkage), generates dendrograms, and silhouette plots.
  - Conducts K-means clustering and visualizes cluster membership.

---

# Outputs

1. **Imputed Datasets:**
   - Data with missing values filled using regression-based imputation.

2. **Plots:**
   - Density and histogram plots for normality checks.
   - Boxplots for group comparisons.
   - PCA scree plots, biplots, and component loading plots.
   - Dendrograms and silhouette plots for clustering.
   - K-means cluster plots.

3. **Statistical Results:**
   - ANOVA summaries, Tukey post-hoc results, and t-tests for each parameter.
   - PCA summaries and loadings.
   - Cluster membership tables and means.

# Conclusion
The R script systematically analyzed water quality parameters such as electrical conductivity (EC), pH, total hardness, calcium, magnesium, sodium, potassium, BOD, COD, nitrate, sulphate, phosphate, and dissolved oxygen across different sampling points and treatments (e.g., river, effluent, mixing zones). 

## Findings:

- The analysis likely revealed significant spatial differences in water quality parameters between the river, effluent, and mixing zones.
- Certain parameters (such as BOD, nitrate, sulphate, and EC) showed strong responses to effluent influence, as indicated by ANOVA and PCA loadings.
- The use of imputation and transformation ensured robust statistical testing.
- Multivariate analysis (PCA and clustering) highlighted that effluent discharge alters the chemical composition of the receiving water, resulting in distinct clusters for river, effluent, and mixed water samples.

