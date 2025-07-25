

library(pls)
library(caret)
library(dplyr)

hsi_wide <- read.csv("E:/DATA/HSI_Metrics.csv")

# Step 1: Prepare data
df <- hsi_wide %>%
  filter(!is.na(SpeciesID))  # remove rows with unknown species

# Step 2: Convert SpeciesID to factor
df$SpeciesID <- as.factor(df$SpeciesID)

# Step 3: Define predictors and response
# spectra_cols <- names(df)[!(names(df) %in% c("TreeID", "SpeciesID"))]
spectra_cols <- grep("^Mean_", names(df), value = TRUE)
X <- df[, spectra_cols]
Y <- df$SpeciesID

# Step 4: Convert Y (factor) to dummy matrix for PLS-DA
Y_dummy <- model.matrix(~ Y - 1)  # One-hot encode

# Step 5: Fit PLS model
set.seed(123)
plsda_model <- plsr(Y_dummy ~ as.matrix(X), ncomp = 35, validation = "CV", scale = TRUE)

# Step 6: Summary and variance explained
summary(plsda_model)
plot(RMSEP(plsda_model), legendpos = "topright")  # to choose number of components

# Step 7: Predict class labels
# Get predicted probabilities for each class
pred_probs <- predict(plsda_model, ncomp = 3, newdata = X)
# Choose class with highest probability
pred_class <- colnames(Y_dummy)[apply(pred_probs, 1, which.max)]
pred_class <- gsub("^Y", "", pred_class)  # clean up label names if needed

# Step 8: Evaluate model
conf_matrix <- confusionMatrix(as.factor(pred_class), Y)
print(conf_matrix)

# Optional: Plot scores (e.g., PLS components 1 and 2)
scores <- scores(plsda_model)
plot(scores[, 2], scores[, 3], col = Y, pch = 16,
     xlab = "PLS Component 1", ylab = "PLS Component 2",
     main = "PLS-DA Score Plot")
legend("topright", legend = levels(Y), col = 1:length(levels(Y)), pch = 16)

loadings <- plsda_model[["loadings"]]

top_loadings_abs <- function(loadings, comp_num, top_n = 5) {
  comp_loadings <- loadings[, comp_num]
  abs_loadings <- abs(comp_loadings)
  top_vars <- sort(abs_loadings, decreasing = TRUE)[1:top_n]
  return(data.frame(Variable = names(top_vars), AbsoluteLoading = top_vars))
}

# Example: Top 5 by absolute loading for component 1
top_loadings_abs(loadings, comp_num = 15, top_n = 5)
