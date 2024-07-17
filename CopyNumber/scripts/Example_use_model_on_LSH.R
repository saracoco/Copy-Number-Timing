library(dplyr)
library(ggplot2)
library(patchwork)
library(loo)
library(bayesplot)
library(cmdstanr)
library(factoextra)
library(ppclust)


setwd("C:/Users/sarac/OneDrive/Documenti/DOCUMENTI/DESKTOP/0_TRIESTE_STAT/0_TRIESTE/CDS/progetto")
source("./CNTiming/R/simulate_functions.R")
source("./CNTiming/R/fitting_functions.R")
source("./CNTiming/R/plotting_functions.R")


# LOAD DATA #
setwd("C:/scratch/CDS_ORFEO/Timing_CDS")
UPN04_extra <- readRDS("Data/extra_cnloh/alpha_beta/UPN04/mutations.rds")
UPN05_extra <- readRDS("Data/extra_cnloh/alpha_beta/UPN05/mutations.rds")
UPN04_alpha_beta <- readRDS("Data/alpha_beta/UPN04/mutations.rds")
UPN05_alpha_beta <- readRDS("Data/alpha_beta/UPN05/mutations.rds")
# FILTERING #
UPN04_extra_NV = UPN04_extra %>% filter(timing_classification %in% c("alpha private", "beta"),PASS == TRUE)
UPN05_extra_NV = UPN05_extra %>% filter(timing_classification  %in% c("alpha private", "beta"),chr == "chr1",PASS == TRUE)
UPN04_alpha_beta_NV = UPN04_alpha_beta %>% ungroup %>% filter(timing_classification %in% c("alpha"),PASS == TRUE)
UPN05_alpha_beta_NV = UPN05_alpha_beta %>% filter(timing_classification %in% c("alpha", "beta"),PASS == TRUE)


data_lsh <- list(UPN04 = UPN04_extra_NV, UPN05 = UPN05_extra_NV, UPN04_LSH = UPN04_alpha_beta_NV, UPN05_LSH = UPN05_alpha_beta_NV)
#names <- c("UPN04","UPN04_LSH", "UPN05", "UPN05_LSH")
#names <- c("UPN04","UPN04_LSH")
names <- c("UPN05", "UPN05_LSH")




data <- dplyr::tibble()
karyo_all <- c()

for(i in 1:length(names)){
  
  data_single <- data_lsh[[names[i]]]

  data_single <- data_single %>% rename(karyo = segment.REL, 
                        DP = DP.REL,
                        NV = NV.REL) %>%
                  mutate(j = paste0(names[i]),
                        segment_id = i,
                        karyo = as.character(karyo))
  
  karyo <- data_single$karyo[1]
  
  data <- dplyr::bind_rows(data, data_single)
  karyo_all <- c(karyo_all, karyo)
}



##### INFERENCE MODEL SEELECTION MULTIPLE K ##################################################################################
results <- fit_model_selection(data, karyo=karyo_all, purity=0.98, INIT=FALSE)

results$model_selection_tibble





#### INFERENCE FO A SINGLE k ##################################################################################


input_data_1 <- prepare_input_data(data, karyo_all, K=1, purity=0.98)
res_1 <- fit_variational(input_data_1, max_attempts = 10, initialization = inits_chain1, INIT = FALSE, initial_iter = 10000)

p <- plotting(res_1,input_data_1,1)
ggsave(paste0("./plots/plot_inference",names[1],"_1.png"), width = 12, height = 16, device = "png", plot=p)



input_data_2 <- prepare_input_data(data, karyo_all, K=2, purity=0.98)
res_2 <- fit_variational(input_data_2, max_attempts = 10, initialization = inits_chain1, INIT = FALSE, initial_iter = 10000)

p <- plotting(res_2,input_data_2,2)
ggsave(paste0("./plots/plot_inference",names[2],"_2.png"), width = 12, height = 16, device = "png", plot=p)



