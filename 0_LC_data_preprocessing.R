## Peatland - Metabolome Project - 2023 Samples 
## Metabolome LC
## Preprocessing


# Load Libraries ----------------------------------------------------------
suppressPackageStartupMessages({  
  library(tidyverse)
  library(readxl)
  source('functions_cdis_exploration_1.R')
  source('functions_cdis_norm_stats.R')
  }) 

#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(x){
  mtb_scaled <- apply(x, 2, AS_helper) 
  return(mtb_scaled)
}

#Log Transformation Functions:
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}

#Log Scaling:
log_transform <- function(x){
  x_nz <- x[ ,which(apply(x, 2, sum) != 0)]
  min.val <- min(abs(x_nz[x_nz!=0]))/10
  x_log_trans <- data.frame(apply(x_nz, 2, log_helper, min.val))
  return(x_log_trans)
}

calc_log2fc <- function(x, factor, groups){
  if(length(groups) != 2){
    stop('Must Be Exactly 2 Groups')
  }
  group1 <- x[x[factor] == groups[1], ][['area']]
  group2 <- x[x[factor] == groups[2], ][['area']]
  l2fc <- log2(mean(group2)/mean(group1))
  return(l2fc)
}

# Function to add spaces between chemical elements and quantities, but keep two-letter elements intact
add_spaces_to_formula <- function(formula) {
  # Add spaces between numbers and the next element
  spaced_formula <- gsub("([0-9])([A-Z])", "\\1 \\2", formula)
  # Add spaces between elements but avoid separating two-letter elements
  spaced_formula <- gsub("([A-Z][a-z]*)(?=[A-Z])", "\\1 ", spaced_formula, perl = TRUE)
  return(spaced_formula)
}

#HILIC -----
### input data:
cd_results_table <- readxl::read_xlsx("Input_Data/LC/SPRUCE_2023_Surface_HILIC_CD_Compounds_Table.xlsx")
# Select columns needed for downstream analysis
cd_results_table <- cd_results_table %>%
  arrange(desc(`Calc. MW`)) %>% 
  select(FeatureID, Name, Formula, `Calc. MW`, contains('Annotation source'), contains('Results'), contains('Area:'), contains('Gap Status:'), contains('Gap Fill Status:')) %>% 
  # Differentiate between features that share the same name using "peak#" at the end of the name
  group_by(Name) %>% 
  add_count(Name) %>% 
  # Create variable with names for plotting (useful in following scripts)
  mutate(name4plot = ifelse(is.na(Name), FeatureID, ifelse(n == 1, Name, paste0(Name, '-isomer', n():1)))) %>% 
  select(-n) %>% 
  ungroup() #%>%
  #rename(Formula = Final_formula) %>%
  #rename(Name = Final_name)
cd_results_table$Formula <- sapply(cd_results_table$Formula, add_spaces_to_formula)

write_csv(cd_results_table, "Tables/hilic_cd_features_table.csv")

# Import metadata and fix names
metadata <- read_csv("Input_Data/LC/SPRUCE_2023_Surface_metadata_HILIC.csv")

# Select only the useful columns and fix column names
metadata <- metadata %>%  
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw'), 
         SampleID = str_remove(SampleID, 'HILIC_Neg_'))
write_csv(metadata, "Tables/fixed_metadata_hilic.csv")

## Thermodynamics ----------------------------------------------------------
# Split formula column into elemental counts
cd_results_table <- separate_formula(cd_results_table)
# Calculate ratios and thermodynamic indices 
cd_results_table <- calc_ratios_n_idxs(cd_results_table)
# Calculate classes
cd_results_table <- calc_classes(cd_results_table)


## Gap-Filled Table - Needed for Statistical Analysis ----------------------
# Gather area under the curve (AUC) values per sample
compounds_table <- cd_results_table %>% 
  select(-contains('Gap Status:'), -contains('Gap Fill Status:')) %>% 
  pivot_longer(contains('Area:'), names_to = 'SampleID', values_to = 'AUC') %>% 
  filter(AUC > 0) %>%  
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw.*'),
         SampleID = str_remove(SampleID, 'HILIC_Neg_'))

# Save gap-filled table to be used in Statistical Analysis
write_csv(compounds_table, "Tables/gap_filled_compounds_table_HILIC.csv")

## Non-gap filled table ----------------------------------------------------
label = FALSE

# Gather "Gap Status" for filtering 
if(label == FALSE){
  gap_status <- cd_results_table %>% 
    select(FeatureID, contains('Gap Status:')) %>% 
    gather(contains('Gap Status:'), key = 'SampleID', value = 'gap_status')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'Gap Status: ')
  gap_status$SampleID <- str_remove(gap_status$SampleID, '.raw.*')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'HILIC_Neg_')
  
  ## Filtering
  
  compounds_table <- left_join(compounds_table, gap_status, by = c('FeatureID', 'SampleID')) %>% 
    filter(gap_status != 'Full gap')
  
  # Plotting types of gaps detected
  
  gap_fill_status <- cd_results_table%>% 
    select(FeatureID, Name, contains('Gap Fill Status:')) %>% 
    pivot_longer(contains('Gap Fill Status:'), names_to = 'SampleID', values_to = 'gap_fill') %>% 
    mutate(SampleID = str_remove(SampleID, 'Gap Fill Status: '),
           SampleID = str_remove(SampleID, '.raw.*'),
           SampleID = str_remove(SampleID, 'HILIC_Neg_')) %>% 
    mutate(gap_fill_type = case_when(gap_fill == 32 ~ 'Filled by spectrum noise',
                                     gap_fill == 128 ~ 'Filled by re-detected peak',
                                     gap_fill == 0 ~ 'No gap to fill',
                                     gap_fill == 16 ~ 'Filled by simulated peak',
                                     gap_fill == 8 ~ 'Filled by trace area',
                                     gap_fill == 64 ~ 'Filled by matching ion')) %>% 
    mutate(gap_fill_type = factor(gap_fill_type, levels = c('No gap to fill', 'Filled by re-detected peak', 
                                                            'Filled by simulated peak', 'Filled by matching ion',
                                                            'Filled by trace area', 'Filled by spectrum noise')))
}

non_gap_cmpds <- merge(gap_fill_status, gap_status, by = c("FeatureID", "SampleID"))

remove_full_gaps <- non_gap_cmpds %>%
  filter(gap_status != 'Full gap')
remove_real_gaps <- non_gap_cmpds%>%
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')

## Filtering

compounds_table_non_gap <- left_join(compounds_table, non_gap_cmpds, by = c('FeatureID', 'SampleID', "Name")) %>% 
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')


# Add metadata information
compounds_table_non_gap <- left_join(compounds_table_non_gap, metadata, by = 'SampleID')
write_csv(compounds_table_non_gap, "Tables/compounds_table_non_real_gap_filled_HILIC.csv")

## make zero matrix 
zero_auc_table <- compounds_table_non_gap %>% 
  select(FeatureID, SampleID, AUC) %>% 
  pivot_wider(names_from = 'SampleID', values_from = 'AUC', values_fill = 0) %>%  
  column_to_rownames(var = 'FeatureID')
write.csv(zero_auc_table, "Tables/zero_auc_table_HILIC.csv")


## Normalization ------
compounds_table <- read.csv("Tables/gap_filled_compounds_table_HILIC.csv")
# Create a new tibble with the AUC per each mass from each sample
auc_table <- compounds_table %>% 
  select(FeatureID, SampleID, AUC) %>% 
  pivot_wider(names_from = 'SampleID', values_from = 'AUC') %>%  
  column_to_rownames(var = 'FeatureID')
write.csv(auc_table, 'Tables/raw_auc_table_hilic.csv', row.names = TRUE)


# Obtaining non-transformed data for differential analysis
norm.matrix.mean <- mean.norm(auc_table)
# Save normalized data, non transformed data for differential analysis
write.csv(norm.matrix.mean, file = "Tables/HILIC_norm_mean.csv", row.names = TRUE)

norm.matrix.median <- median.norm(auc_table)
norm.matrix.max <- max.norm(auc_table)
write.csv(norm.matrix.max, file = "Tables/HILIC_norm_max.csv", row.names = TRUE)

## Pareto
HILIC_raw <- auc_table%>%
  t()%>%
  as.data.frame()

HILIC_clean <- HILIC_raw %>%
  log_transform() %>%
  pareto_scale() 
# Remove the "X" in front of column names if they exist
colnames(HILIC_clean) <- gsub("^X", "", colnames(HILIC_clean))

write.csv(HILIC_clean, file = "Tables/HILIC_norm_pareto.csv", row.names = TRUE)

HILIC_full <- HILIC_clean %>%
  rownames_to_column('SampleID') %>%
  left_join(metadata)
test_num <- sample(1:ncol(HILIC_raw), 300)

#Check the distribution - Clearly bad
boxplot(HILIC_raw[test_num])

#Much better
boxplot(HILIC_clean[test_num])


# RP -----
### input data:
cd_results_table <- readxl::read_xlsx("Input_Data/LC/SPRUCE_2023_Surface_RP_CD_Compounds_Table.xlsx")
# Select columns needed for downstream analysis
cd_results_table <- cd_results_table %>%
  arrange(desc(`Calc. MW`)) %>% 
  select(FeatureID, Name, Formula, `Calc. MW`, contains('Annotation source'), contains('Results'), contains('Area:'), contains('Gap Status:'), contains('Gap Fill Status:')) %>% 
  # Differentiate between features that share the same name using "peak#" at the end of the name
  group_by(Name) %>% 
  add_count(Name) %>% 
  # Create variable with names for plotting (useful in following scripts)
  mutate(name4plot = ifelse(is.na(Name), FeatureID, ifelse(n == 1, Name, paste0(Name, '-isomer', n():1)))) %>% 
  select(-n) %>% 
  ungroup() #%>%
#rename(Formula = Final_formula) %>%
#rename(Name = Final_name)
cd_results_table$Formula <- sapply(cd_results_table$Formula, add_spaces_to_formula)

write_csv(cd_results_table, "Tables/rp_cd_features_table.csv")

# Import metadata and fix names
metadata <- read_csv("Input_Data/LC/SPRUCE_2023_Surface_metadata_RP.csv")

# Select only the useful columns and fix column names
metadata <- metadata %>%  
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw'), 
         SampleID = str_remove(SampleID, 'RP_Pos_'))
write_csv(metadata, "Tables/fixed_metadata_rp.csv")

## Thermodynamics ----------------------------------------------------------
# Split formula column into elemental counts
cd_results_table <- separate_formula(cd_results_table)
# Calculate ratios and thermodynamic indices 
cd_results_table <- calc_ratios_n_idxs(cd_results_table)
# Calculate classes
cd_results_table <- calc_classes(cd_results_table)


## Gap-Filled Table - Needed for Statistical Analysis ----------------------
# Gather area under the curve (AUC) values per sample
compounds_table <- cd_results_table %>% 
  select(-contains('Gap Status:'), -contains('Gap Fill Status:')) %>% 
  pivot_longer(contains('Area:'), names_to = 'SampleID', values_to = 'AUC') %>% 
  filter(AUC > 0) %>%  
  mutate(SampleID = str_remove(SampleID, 'Area: '),
         SampleID = str_remove(SampleID, '.raw.*'),
         SampleID = str_remove(SampleID, 'RP_Pos_'))

# Save gap-filled table to be used in Statistical Analysis
write_csv(compounds_table, "Tables/gap_filled_compounds_table_RP.csv")



## Non-gap filled table ----------------------------------------------------
label = FALSE

# Gather "Gap Status" for filtering 
if(label == FALSE){
  gap_status <- cd_results_table %>% 
    select(FeatureID, contains('Gap Status:')) %>% 
    gather(contains('Gap Status:'), key = 'SampleID', value = 'gap_status')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'Gap Status: ')
  gap_status$SampleID <- str_remove(gap_status$SampleID, '.raw.*')
  gap_status$SampleID <- str_remove(gap_status$SampleID, 'RP_Pos_')
  
  ## Filtering
  
  compounds_table <- left_join(compounds_table, gap_status, by = c('FeatureID', 'SampleID')) %>% 
    filter(gap_status != 'Full gap')
  
  # Plotting types of gaps detected
  
  gap_fill_status <- cd_results_table%>% 
    select(FeatureID, Name, contains('Gap Fill Status:')) %>% 
    pivot_longer(contains('Gap Fill Status:'), names_to = 'SampleID', values_to = 'gap_fill') %>% 
    mutate(SampleID = str_remove(SampleID, 'Gap Fill Status: '),
           SampleID = str_remove(SampleID, '.raw.*'),
           SampleID = str_remove(SampleID, 'RP_Pos_')) %>% 
    mutate(gap_fill_type = case_when(gap_fill == 32 ~ 'Filled by spectrum noise',
                                     gap_fill == 128 ~ 'Filled by re-detected peak',
                                     gap_fill == 0 ~ 'No gap to fill',
                                     gap_fill == 16 ~ 'Filled by simulated peak',
                                     gap_fill == 8 ~ 'Filled by trace area',
                                     gap_fill == 64 ~ 'Filled by matching ion')) %>% 
    mutate(gap_fill_type = factor(gap_fill_type, levels = c('No gap to fill', 'Filled by re-detected peak', 
                                                            'Filled by simulated peak', 'Filled by matching ion',
                                                            'Filled by trace area', 'Filled by spectrum noise')))
}

non_gap_cmpds <- merge(gap_fill_status, gap_status, by = c("FeatureID", "SampleID"))

remove_full_gaps <- non_gap_cmpds %>%
  filter(gap_status != 'Full gap')
remove_real_gaps <- non_gap_cmpds%>%
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')

## Filtering

compounds_table_non_gap <- left_join(compounds_table, non_gap_cmpds, by = c('FeatureID', 'SampleID', "Name")) %>% 
  filter(gap_fill_type !='Filled by spectrum noise')%>%
  filter(gap_fill_type !='Filled by trace area')


# Add metadata information
compounds_table_non_gap <- left_join(compounds_table_non_gap, metadata, by = 'SampleID')
write_csv(compounds_table_non_gap, "Tables/compounds_table_non_real_gap_filled_RP.csv")

## make zero matrix 
zero_auc_table <- compounds_table_non_gap %>% 
  select(FeatureID, SampleID, AUC) %>% 
  pivot_wider(names_from = 'SampleID', values_from = 'AUC', values_fill = 0) %>%  
  column_to_rownames(var = 'FeatureID')
write.csv(zero_auc_table, "Tables/zero_auc_table_RP.csv")



## Normalization ------
compounds_table <- read.csv("Tables/gap_filled_compounds_table_rp.csv")
# Create a new tibble with the AUC per each mass from each sample
auc_table <- compounds_table %>% 
  select(FeatureID, SampleID, AUC) %>% 
  pivot_wider(names_from = 'SampleID', values_from = 'AUC') %>%  
  column_to_rownames(var = 'FeatureID')
write.csv(auc_table, 'Tables/raw_auc_table_RP.csv', row.names = TRUE)

# Obtaining non-transformed data for differential analysis
norm.matrix.mean <- mean.norm(auc_table)
# Save normalized data, non transformed data for differential analysis
write.csv(norm.matrix.mean, file = "Tables/RP_norm_mean.csv", row.names = TRUE)

## Pareto
RP_raw <- auc_table%>%
  t()%>%
  as.data.frame()

RP_clean <- RP_raw %>%
  log_transform() %>%
  pareto_scale() 
# Remove the "X" in front of column names if they exist
colnames(RP_clean) <- gsub("^X", "", colnames(RP_clean))

write.csv(RP_clean, file = "Tables/RP_norm_pareto.csv", row.names = TRUE)

RP_full <- RP_clean %>%
  rownames_to_column('SampleID') %>%
  left_join(metadata)
test_num <- sample(1:ncol(RP_raw), 300)

#Check the distribution - Clearly bad
boxplot(RP_raw[test_num])

#Much better
boxplot(RP_clean[test_num])


