

library(pdp)
library(CAST)

library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(tidyverse)
library(corrr)
library(igraph)
library(vegan)
library(beepr)

# Load data
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

colnames(samples_clean)

# Identify VI columns (exclude ID columns)
vi_cols <- names(combined_df)[!(names(combined_df) %in% c("TreeID", "SpeciesID"))]

Canopies_noVIs <- c("080736","090800","090830","090848","090863",
                    "100987","110563","110602","110618","110752","110819",
                    "110883","110913","120866","140806","150720","80917")

samples_clean <- combined_df[!combined_df$TreeID %in% Canopies_noVIs, ]
print(vi_cols)
vi_cols <- c("Mean_Intensity_LeafOn", "Max_PAD_30_35_on", "Mean_Height_LeafOn")
# Remove constant or all-NA columns
vi_data_clean <- samples_clean %>%
   select(all_of(vi_cols)) %>%
   select(where(~ sd(.x, na.rm = TRUE) > 0))



#######################################################################

help("ffs")
set.seed(10)
ffsmodel <- ffs((samples_clean)[,vi_cols],
                samples_clean$SpeciesID,
                method="rf", 
                tuneGrid=data.frame("mtry"=2),
                verbose=FALSE,
                ntree=25,
                trControl=trainControl(method="cv",
                                       number = 5,
                                       savePredictions = "final"))
beep()

ffsmodel$selectedvars

