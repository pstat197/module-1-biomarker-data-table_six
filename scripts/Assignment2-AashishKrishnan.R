# Question 1: Why Log-Transform Protein Levels?
# ==============================================
library(tidyverse)
library(moments)

# Load data
data_raw <- read_csv("/Users/aashishkrishnan/Desktop/biomarker-raw.csv")

# Select 4 random proteins (exclude non-protein columns like ID, group)
set.seed(123)
protein_cols <- data_raw %>% select(where(is.numeric)) %>% names()
cat("Total numeric columns found:", length(protein_cols), "\n")

# Sample 4 proteins (or fewer if less than 4 available)
n_sample <- min(4, length(protein_cols))
sample_proteins <- sample(protein_cols, n_sample)
cat("Sampling", n_sample, "proteins for analysis\n\n")

# Visualize raw vs log-transformed
par(mfrow = c(2, 4))
for (protein in sample_proteins) {
  values <- data_raw[[protein]][data_raw[[protein]] > 0]
  
  # Raw distribution
  hist(values, main = paste("Raw:", substr(protein, 1, 15)),
       xlab = "Value", col = "steelblue", breaks = 30)
  
  # Log-transformed distribution
  hist(log(values), main = paste("Log:", substr(protein, 1, 15)),
       xlab = "Log(Value)", col = "coral", breaks = 30)
}


# Protein levels are log-transformed to address the severe positive skewness observed in raw concentration data. 
# As shown in Figure 1, raw protein values for Protein 4.1 exhibit a heavily right-skewed distribution with most observations clustered near zero and a long tail of extreme values. After log-transformation, the distribution becomes approximately normal and symmetric. This transformation is necessary because: (1) protein concentrations follow log-normal distributions due to multiplicative biological processes, (2) it stabilizes variance across different concentration ranges, and (3) it satisfies the normality assumptions required for parametric statistical tests such as t-tests and logistic regression used in subsequent variable selection procedures.