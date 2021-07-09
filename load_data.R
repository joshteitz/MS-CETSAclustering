library(conflicted)
library(tidyverse)
library(readxl)
library(here)
library(drc)
library(progress)

conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("map", "purrr")

### Introduction ###################################################################################################
# This file contains functions for loading data from the "Supplementary Materials" of the paper
# "Thermal proximity coaggregation for system-wide profiling of protein complex dynamics in cells" by Tan et al..
# This paper is accessible from the following url: https://science.sciencemag.org/content/359/6380/1170.
#
# In order for these functions to work, make sure to:
#   1. Download "Tables S1 to S27" from https://science.sciencemag.org/content/suppl/2018/02/07/science.aan0346.DC1
#   2. Unzip folder and place it in the "data" directory.
#   3. Ensure the folder is named "TablesS1_to_S27" and contains two Excel files (.xlsx) named
#      "Tables_S1_to_S18.xlsx" and "Tables_S19_to_S27.xlsx".
#   ** If either Excel file is open, the code will not run.
#
# This file also contains a function called `param_mdata` that transforms the feature space of the melting data
# given as input. `convert_long` and `make_unit_var` are helper functions. 
####################################################################################################################

# Location of the two Excel spreadsheets.
LOC <- here("data", "TablesS1_to_S27")

# Load protein melting data from a cell line.
load_mdata <- function(cell_line = "K562") {
  file_sheet = switch(
    cell_line,
    "K562 lysate" = c("Tables_S1_to_S18.xlsx", "Table S6"),
    "K562" = c("Tables_S1_to_S18.xlsx", "Table S7"),
    "A375" = c("Tables_S19_to_S27.xlsx", "Table S19"),
    "HCT116" = c("Tables_S19_to_S27.xlsx", "Table S20"),
    "HEK293T" = c("Tables_S19_to_S27.xlsx", "Table S21"),
    "HL60" = c("Tables_S19_to_S27.xlsx", "Table S22"),
    "MCF7" = c("Tables_S19_to_S27.xlsx", "Table S23"),
    "mouse liver" = c("Tables_S19_to_S27.xlsx", "Table S24")
  )
  if (is.null(file_sheet)) {
    stop("Unrecognizable cell line! Please enter one of: \"K562 lysate\", \"K562\", \"A375\", \"HCT116\",
         \"HEK293T\", \"HL60\", \"MCF7\", or \"mouse liver\".")
  }
  mdata <- read_excel(paste0(LOC, "/", file_sheet[1]), sheet = file_sheet[2], skip = 2) %>%
    select(Protein = Accession, matches("^T\\d+$"))
  return(mdata)
}

# Load CORUM protein complex data
load_cdata <- function() {
  cdata <- read_excel(paste0(LOC, "/Tables_S1_to_S18.xlsx"), sheet = "Table S8", skip = 2) %>%
    select(Complex = Complex.id, Proteins = subunits..UniProt.IDs.) %>%
    mutate(Complex = Complex %>% as.integer) %>%
    mutate(Proteins = strsplit(Proteins, ","))
  cdata <- map2_dfr(cdata$Complex, cdata$Proteins, ~ tibble(Complex = .x, Protein = .y))
  return(cdata)
}

# Convert melting data to long format
convert_long <- function(mdata) {
  mdata_long <- mdata %>%
    pivot_longer(
      cols = matches("T\\d+$"), 
      names_to = "Temp", 
      values_to = "Sol", 
      names_prefix = "T", 
      names_transform = list(Temp = as.integer)
    )
  return(mdata_long)
}

# For each protein's melting data, fit a three parameter log-logistic function with upper limit equal to 1.
# Then extract the function's three parameters. Finally, scale each parameter to be unit variance.
param_mdata <- function(mdata) {
  pb <- progress_bar$new(total = nrow(mdata))
  mparams <- convert_long(mdata) %>%
    group_by(Protein) %>%
    nest() %>%
    mutate(Mod = map(data, ~ {pb$tick(); drm(Sol ~ Temp, data = .x, fct = LL.3u())})) %>%
    select(-data) %>%
    ungroup() %>%
    mutate(Param_b = map_dbl(Mod, ~ .x$coefficients[1])) %>%
    mutate(Param_c = map_dbl(Mod, ~ .x$coefficients[2])) %>%
    mutate(Param_e = map_dbl(Mod, ~ .x$coefficients[3])) %>%
    select(-Mod)
  return(mparams)
}

# Make specified columns unit variance.
make_unit_var <- function(mdata, numeric_cols = NULL) {
  scale1 <- function(x) {
    x_sd <- sd(x)
    x_mean <- mean(x)
    return((x - x_mean) / x_sd)
  }
  mdata_unit_var <- 
    mdata %>%
    mutate(across({{numeric_cols}}, scale1))
  return(mdata_unit_var)
}

