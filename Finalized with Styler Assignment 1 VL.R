##*************************************************************************
## SOFTWARE TOOLS 2025 - Assignment 1 Script
## Author: Vian Lelo
## Date: 2025-10-17
## Project: Cryptic Diversity Patterns in Deep-Sea Asteroids
##*************************************************************************

## ============================================================
## 0. Install Unique PACKAGES -----
## ============================================================
#this is a test change
# install packages:
install.packages("sf") # Simple Features for spatial data
install.packages("rnaturalearth") # World map polygons (Natural Earth)
install.packages("rnaturalearthdata") # Supporting data for rnaturalearth
install.packages("ggrepel") # Non-overlapping text labels for ggplot2

## ============================================================
## 1. LOAD PACKAGES AND SET GLOBAL THEMES -----
## ============================================================

# Justification:
#   Importing required packages for data wrangling (tidyverse),
#   spatial mapping (sf, rnaturalearth), visualization (ggplot2, viridis, ggrepel),
#   and string manipulation (Biostrings).
#   Conflicted package ensures correct function precedence.

library(conflicted)
library(ggrepel)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(viridis)
library(Biostrings)

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflicts_prefer(dplyr::rename)

theme_set(theme_light()) # consistent visualization theme across plots

# ==> Startup complete; environment ready

## ============================================================
## 2. PURPOSE AND OVERVIEW -----
## ============================================================

# Purpose:
#   Investigate cryptic diversity in
#   using DNA barcode data from BOLD Systems.
# Overview:
#   1. Import and clean metadata from BOLD output
#   2. Extract coordinates, BINs, and species
#   3. Classify species as "cryptic" vs "noncryptic"
#   4. Visualize spatial and regional diversity patterns

## ============================================================
## 3. OBTAINING AND CLEANING DATA FROM BOLD -----
## ============================================================

# set working directory
#setwd("../R/") # adjust as needed for your environment
# check working directory
getwd() # verify working directory

# Load dataset
dfastropec <- read_tsv("Astropectinidae_result.tsv")

# Quick exploration
class(dfastropec)
dim(dfastropec)
summary(dfastropec)

## ---- Rename columns for consistency ----
dfastropec <- dfastropec %>%
  rename(
    country_ocean = `country/ocean`,
    province_state = `province/state`
  )

## ---- Count barcodes by country ----
dfastropec %>%
  group_by(country_ocean) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

sum(!is.na(dfastropec$country_ocean)) # ==> Total records with country data

## ---- Extract latitude and longitude with less lines of code ----
dfastropec <- dfastropec %>% 
  mutate(coord_clean = str_remove_all(coord_clean, "\\[|\\]")) %>% 
  separate(coord_clean, into = c("longitude", "latitude"), sep = ",", convert = T)

# Confirm numeric conversion
summary(select(dfastropec, latitude, longitude))

## ---- Extract latitude and longitude from 'coord' column ----
#dfastropec <- dfastropec %>%
#  mutate(
#    coord_clean = str_remove_all(coord, "\\[|\\]"), # clean brackets
#    longitude = as.numeric(str_split_fixed(coord_clean, ",", 2)[, 1]),
#    latitude = as.numeric(str_split_fixed(coord_clean, ",", 2)[, 2])
#  )

# (>_<) In some BOLD exports, latitude/longitude can be flipped. Will verify later.

## ---- Group by province for a quick geographical sense ----
dfastro_countcountry <- dfastropec %>%
  group_by(province_state) %>%
  filter(!is.na(province_state)) %>%
  summarise(avgLat = mean(latitude, na.rm = TRUE)) %>%
  arrange(desc(avgLat))

dfastro_countcountry # ==> Table of provinces with mean latitude

## ============================================================
## 4. SELECT AND CLEAN ESSENTIAL COLUMNS -----
## ============================================================

dfastrostar <- dfastropec %>%
  select(species, bin_uri, latitude, longitude) %>%
  filter(!is.na(species), !is.na(bin_uri), !is.na(latitude), !is.na(longitude))

glimpse(dfastrostar)
# (*_*) Corresponds to Table S1: cleaned occurrence data

## ============================================================
## 5. GROUPING AND CRYPTIC DIVERSITY CLASSIFICATION -----
## ============================================================

# ---- Count BINs per species ----
species_bin_counts <- dfastrostar %>%
  group_by(species, bin_uri) %>%
  summarise(n_records = n(), .groups = "drop")

species_bin_summary <- species_bin_counts %>%
  group_by(species) %>%
  summarise(
    n_bins = n_distinct(bin_uri),
    total_records = sum(n_records),
    .groups = "drop"
  ) %>%
  arrange(desc(n_bins))

# ==> Species with multiple BINs are potential cryptic taxa

multi_bin_species <- species_bin_summary %>% filter(n_bins > 1)
head(multi_bin_species)

## ---- Join back to main dataset ----
dfastrostar <- dfastrostar %>%
  left_join(species_bin_summary, by = "species")

## ---- Classify cryptic status ----
dfastrostar <- dfastrostar %>%
  mutate(
    cryptic_status = if_else(n_bins > 1, "cryptic", "noncryptic"),
    cryptic_status = factor(cryptic_status, levels = c("noncryptic", "cryptic"))
  )

table(dfastrostar$cryptic_status) # ==> Distribution summary

# (*_*) Supports Figure 2A: Counts of cryptic vs noncryptic species

## ============================================================
## 5B. GROUP BY SPECIES AND BIN ID – RECORD-LEVEL REALIZATION -----
## ============================================================

# Justification:
# Although total record counts per BIN are informative, high sampling of a few BINs could bias cryptic diversity estimates if certain species or BINs dominate the dataset. Here, we explicitly calculate BIN richness per morphospecies to validate that cryptic diversity is not an artifact of uneven record sampling.

# (*_*) Analytical check to prevent sampling bias influencing diversity interpretation

species_bin_counts <- dfastrostar %>%
  group_by(species, bin_uri) %>%
  summarise(
    specimen_count = n(), # number of records for each species–BIN combination
    .groups = "drop"
  )

species_bin_summary_check <- species_bin_counts %>%
  group_by(species) %>%
  summarise(
    n_bins = n_distinct(bin_uri), # unique BINs per species
    total_specimens = sum(specimen_count),
    .groups = "drop"
  ) %>%
  arrange(desc(n_bins))

# View top results
head(species_bin_summary_check, 10)

# Interpretation:
#   This summary provides a secondary confirmation of BIN diversity per morphospecies. Species with >1 BIN but moderate specimen counts likely represent true cryptic taxa. Conversely, species with high record counts but single BINs suggest strong barcode stability and low cryptic diversity.

# (>_<) Dead-end caution:
#   Without this grouping check, cryptic diversity patterns could be confounded by record overrepresentation (e.g., hundreds of identical BIN entries for one species).

# ==> Supports the robustness of our classification and avoids sampling bias
# (*_*) Supplementary Table S2: BIN count validation per morphospecies


## ============================================================
## 5B-1. VISUALIZATION A – BIN COUNT VS SPECIMEN COUNT -----
## ============================================================

# Justification:
#   To visualize whether heavily sampled species (many specimens) also appear to have more BINs, which could suggest sampling bias rather than genuine cryptic diversity.

p_bias_scatter <- ggplot(
  species_bin_summary_check,
  aes(x = total_specimens, y = n_bins)
) +
  geom_point(aes(color = n_bins > 1),
    size = 3.5, alpha = 0.8
  ) +
  scale_color_viridis_d(
    name = "Cryptic (>1 BIN)",
    option = "C",
    direction = -1,
    labels = c("No", "Yes")
  ) +
  geom_smooth(
    method = "lm", se = TRUE,
    color = "gray30", linetype = "dashed"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Relationship between Sampling Effort and BIN Diversity",
    x = "Total Specimens per Species",
    y = "Number of BINs"
  ) 

p_bias_scatter
# ==> Scatterplot showing whether BIN richness scales with specimen number.
# (*_*) Figure S2A: bias assessment of BIN richness vs sampling effort.


## ============================================================
## 5B-2. VISUALIZATION B – DISTRIBUTION OF BIN COUNTS -----
## ============================================================

# Justification:
#   To examine how many species possess 1, 2, 3, ... BINs.
#   This shows whether most taxa are non-cryptic (1 BIN) or show multiple BINs.

p_bin_distribution <- ggplot(
  species_bin_summary_check,
  aes(x = factor(n_bins))
) +
  geom_bar(
    fill = viridis::viridis(6, option = "D")[3],
    color = "black", alpha = 0.9, width = 0.7
  ) +
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    vjust = -0.3, fontface = "bold", size = 5
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black")
  ) +
  labs(
    title = "Distribution of BIN Counts per Morpho-Species",
    x = "Number of BINs",
    y = "Number of Species"
  ) +
  expand_limits(y = max(table(species_bin_summary_check$n_bins)) * 1.15)

p_bin_distribution
# ==> Histogram reveals that most species have 1 BIN, confirming low cryptic inflation.
# (*_*) Figure S2B: overall distribution of BIN richness per morphospecies.

## ============================================================
## 6. VISUALIZATION 1 – CRYPTIC VS NON-CRYPTIC COUNTS -----
## ============================================================

count_df <- dfastrostar %>%
  count(cryptic_status)

p_overall <- ggplot(count_df, aes(x = cryptic_status, y = n, fill = cryptic_status)) +
  geom_col(width = 0.65, color = "black", alpha = 0.9) +
  geom_text(aes(label = n), vjust = -0.3, size = 6, fontface = "bold") +
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Overall Count of Cryptic vs Non-Cryptic Species",
    y = "Number of Records"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +   # ✅ centers the title
  expand_limits(y = max(count_df$n) * 1.15)

p_overall
# ==> Bar chart correctly displays relative frequencies
# (*_*) Figure 2A – summary counts of cryptic vs non-cryptic species

## ============================================================
## 7. CLASSIFY SPECIES BY REGION (TROPICAL VS TEMPERATE) -----
## ============================================================

# Justification:
#   Latitude < ±23.5° defines tropical regions (NASA Earth Observatory, 2024). Used for broad biogeographic comparison.

dfastrostar <- dfastrostar %>%
  mutate(
    region = if_else(abs(latitude) < 23.5, "Tropical", "Temperate"),
    region = factor(region, levels = c("Tropical", "Temperate"))
  )

table(dfastrostar$region, dfastrostar$cryptic_status)

## ============================================================
## 8. VISUALIZATION 2 – REGIONAL COMPARISON -----
## ============================================================

# ==> Tropical zones appear to host slightly higher cryptic diversity
# (*_*) Figure 2B – regional cryptic/noncryptic comparison

p_region <- dfastrostar %>%
  count(region, cryptic_status) %>%
  ggplot(aes(x = region, y = n, fill = cryptic_status)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, color = "black") +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.7),
    vjust = -0.3, size = 5, fontface = "bold"
  ) +
  scale_fill_viridis_d(option = "C", direction = -1, name = "Status") +
  theme_minimal(base_size = 16) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10), # padding
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    # Prevent clipping of text labels above panel
    plot.background = element_rect(fill = "white", color = NA),
    plot.title.position = "plot",
    panel.ontop = FALSE,
    clip = "off" # key fix: allows text above top of panel
  ) +
  labs(
    title = "Cryptic versus Non-Cryptic Species in Tropical versus Temperate Regions",
    x = "Region",
    y = "Number of Records"
  ) +
  # dynamically extend y limit by 10–15%
  expand_limits(y = max(dfastrostar %>% count(region, cryptic_status) %>% pull(n)) * 1.15)

p_region

## ============================================================
## 8B. VISUALIZATION 2 (SPECIES-LEVEL) – REGIONAL COMPARISON -----
## ============================================================

# Justification:
#   To avoid record-level bias, we again aggregate at the species level.
#   This ensures that regional comparisons (Tropical vs Temperate) are based on distinct species rather than heavily sampled BINs or individuals.

# (*_*) Figure 2B (species-level): regional cryptic diversity adjusted for sampling effort

species_region_summary <- dfastrostar %>%
  distinct(species, region, cryptic_status) %>% # one row per species-region
  count(region, cryptic_status)

p_region_species <- ggplot(
  species_region_summary,
  aes(x = region, y = n, fill = cryptic_status)
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.65, color = "black"
  ) +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.7),
    vjust = -0.3, size = 5, fontface = "bold"
  ) +
  scale_fill_viridis_d(option = "C", direction = -1, name = "Status") +
  theme_minimal(base_size = 16) +
  theme(
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    plot.title.position = "plot",
    panel.ontop = FALSE,
    clip = "off"
  ) +
  labs(
    title = "Cryptic versus Non-Cryptic Species by Region (Species-Level)",
    x = "Region",
    y = "Number of Species"
  ) +
  expand_limits(y = max(species_region_summary$n) * 1.15)

p_region_species
# ==> Each bar now represents the number of distinct species, not records
# (*_*) Figure 2B (revised): unbiased regional comparison of cryptic diversity

## ---- Statistical Test: Are cryptic proportions different by region? ----
region_table <- table(species_region_summary$region, species_region_summary$cryptic_status)
chisq_test <- chisq.test(region_table)
chisq_test

#Output readable summary
cat("Chi-squared test for independence between region and cryptic status: \n")
cat("Chi-squared = ", round(chisq_test$statistic, 3), "  df = ", chisq_test$parameter, "  p-value = ", signif(chisq_test$p.value, 3), "\n")

## ---- Add chi-squares test result as a caption ----
region_table <- table(species_region_summary$region, 
  species_region_summary$cryptic_status)
chisq_test <- chisq.test(region_table)

#Rebuild the plot with the caption
p_region_species <- p_region_species + labs(caption = paste0(
  "Chi-squared = ", round(chisq_test$statistic, 2),
  ", df = ", chisq_test$parameter,
  ", p = ", signif(chisq_test$p.value, 3)
)) + theme(plot.caption = element_text(hjust = 1, size = 12, face = 'italic'))

p_region_species

## ============================================================
## 9. FIX LATITUDE/LONGITUDE SWAP (IDENTIFIED ERROR) -----
## ============================================================

# (>_<) Dead-end issue discovered:
#   Latitude and longitude were reversed in raw data,
#   producing incorrect map projections (latitudes up to ±179°).

dfastrostar <- dfastrostar %>%
  rename(
    old_latitude = latitude,
    old_longitude = longitude
  ) %>%
  mutate(
    latitude = old_longitude,
    longitude = old_latitude
  ) %>%
  select(-old_latitude, -old_longitude)

summary(dfastrostar$latitude)
summary(dfastrostar$longitude)
# ==> Confirm latitudes now within ±90°, longitudes within ±180°

## ============================================================
## 10. CREATE SPECIES-LEVEL SUMMARY FOR MAPPING -----
## ============================================================

species_summary <- dfastrostar %>%
  group_by(species, cryptic_status) %>%
  summarise(
    mean_latitude = mean(latitude, na.rm = TRUE),
    mean_longitude = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(mean_latitude), !is.na(mean_longitude))

# (*_*) Dataset used for Figure 3: global map of species centroids

## ============================================================
## 11. VISUALIZATION 3 – GEOGRAPHIC MAP -----
## ============================================================

world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_point(
    data = species_summary,
    aes(x = mean_longitude, y = mean_latitude, color = cryptic_status),
    size = 3.2, alpha = 0.9
  ) +
  scale_color_viridis_d(option = "C", direction = -1, name = "Status") +
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Geographic Distribution of Cryptic versus Non-Cryptic Species",
    x = "Longitude", y = "Latitude"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# ==> Points represent mean species centroids colored by cryptic status
# (*_*) Figure 3 – spatial pattern of cryptic diversity

