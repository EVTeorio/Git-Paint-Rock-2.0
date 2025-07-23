

All_Metrics <- read.csv("E:/DATA/All_CanopyMetrics.csv")
print(All_Metrics)

values <- names(All_Metrics)[!(names(All_Metrics) %in% c("TreeID", "SpeciesID"))]

species_values <- All_Metrics %>%
  group_by(SpeciesID) %>%
  summarise(across(values, mean, na.rm = TRUE)) 

write.csv(species_values, "E:/DATA/All_SpeciesMetrics.csv")
