
library(dplyr)
library(readr)
library(tibble)

# ---- Step 1: Load the data ----
combined_df <- read_csv("E:/DATA/All_CanopyMetrics.csv")

# ---- Step 2: Identify VI columns (exclude TreeID and SpeciesID) ----
id_cols <- c("TreeID", "SpeciesID")
vi_cols <- setdiff(names(combined_df), id_cols)

# ---- Step 3: Calculate stats for each VI column ----
global_vi_stats <- combined_df %>%
  summarise(across(
    all_of(vi_cols),
    list(
      Mean = ~mean(.x, na.rm = TRUE),
      SD   = ~sd(.x, na.rm = TRUE),
      Min  = ~min(.x, na.rm = TRUE),
      Max  = ~max(.x, na.rm = TRUE)
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  separate(Metric, into = c("Variable", "Stat"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = Stat, values_from = Value)

species_vi_stats <- combined_df %>%
  group_by(SpeciesID) %>%
  summarise(across(
    all_of(vi_cols),
    list(
      Mean = ~mean(.x, na.rm = TRUE),
      SD   = ~sd(.x, na.rm = TRUE),
      Min  = ~min(.x, na.rm = TRUE),
      Max  = ~max(.x, na.rm = TRUE)
    )
  ))  # Wide format — one row per species

write_csv(species_vi_stats, "E:/DATA/VI_Stats_By_Species.csv")

# ---- Compute stats per Tree within Species ----
tree_within_species_stats <- combined_df %>%
  group_by(SpeciesID, TreeID) %>%
  summarise(across(
    all_of(vi_cols),
    list(
      Mean = ~mean(.x, na.rm = TRUE),
      SD   = ~sd(.x, na.rm = TRUE),
      Min  = ~min(.x, na.rm = TRUE),
      Max  = ~max(.x, na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  ungroup()


# ---- Step 4: View or save ----
print(head(global_vi_stats))
write_csv(global_vi_stats, "E:/DATA/VI_Stats_By_Variable.csv")




#############################################################################

library(ggplot2)


# ---- Long format for plotting ----
global_long <- global_vi_stats %>%
  pivot_longer(cols = c(Mean, SD, Min, Max), names_to = "Stat", values_to = "Value")

# ---- Bar plot of Mean per VI ----
ggplot(global_long %>% filter(Stat == "Mean"), aes(x = reorder(Variable, Value), y = Value)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Mean of Each Vegetation Index", x = "VI", y = "Mean Value") +
  theme_minimal()

###############################################################################

library(reshape2)

# ---- Compute species-wise VI mean ----
species_vi_means <- combined_df %>%
  group_by(SpeciesID) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm = TRUE))) %>%
  ungroup()

# ---- Reshape to long format for heatmap ----
vi_long <- species_vi_means %>%
  pivot_longer(-SpeciesID, names_to = "VI", values_to = "MeanValue")

# ---- Heatmap ----
ggplot(vi_long, aes(x = VI, y = SpeciesID, fill = MeanValue)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +
  labs(title = "Mean VI Values by Species", x = "Vegetation Index", y = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
