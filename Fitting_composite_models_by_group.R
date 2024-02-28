### fitting the composite group indicator models for the State of Canada's Birds

library(tidyverse)
library(cmdstanr)
library(readxl)
library(naturecounts)
library(patchwork)

source("functions/GAM_basis_function_mgcv.R")

# Quantile functions to generate 95% CI instead of stan's default 90
q2_5 <- function(x)c(q2_5 = quantile(x,probs = c(0.025),
                                     names = FALSE))
q97_5 <- function(x)c(q97_5 = quantile(x,probs = c(0.975),
                                       names = FALSE))




base_year <- 1970




# extract state of canada's birds -----------------------------------------
re_download <- FALSE

if(re_download){

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbTrendRank") -> rank_tbl
  saveRDS(rank_tbl,"data/SocbTrendRank.rds")

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbSpecies") -> sp_tbl
  saveRDS(sp_tbl,"data/SocbSpecies.rds")

  #Information on species' group#
  nc_query_table(username = "adam.smith",
                 table = "Groups") -> group_tbl
  saveRDS(group_tbl,"data/Groups.rds")

  #Data on all annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndices") -> indices_tbl
  saveRDS(indices_tbl,"data/TrendsIndices.rds")

  #Data on goal-based annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndicesGoals") -> goal_indices_tbl
  saveRDS(goal_indices_tbl,"data/TrendsIndicesGoals.rds")

  #Data on species trends #
  nc_query_table(username = "adam.smith",
                 table = "Trends") -> trends_tbl
  saveRDS(trends_tbl,"data/Trends.rds")

}else{
  sp_tbl <- readRDS("data/SocbSpecies.rds")
  group_tbl <- readRDS("data/Groups.rds")
  trend_tbl <- readRDS("data/Trends.rds")
  rank_tbl <- readRDS("data/SocbTrendRank.rds")
  indices_tbl <- readRDS("data/TrendsIndices.rds")
  goal_indices_tbl <- readRDS("data/TrendsIndicesGoals.rds")

}

# Generate species-level smooths ------------------------------------------

sp_simple <- sp_tbl %>%
  select(speciesCode,speciesID) %>% #,population,popType,popID) %>%
  rename(speciesId = speciesID) %>%
  distinct()
goal_indices_tbl<- goal_indices_tbl %>%
  left_join(.,sp_simple,
            by = "speciesId") %>%
  filter(year >= base_year) %>%
  distinct()

sp_gt_1_pop <- goal_indices_tbl %>%
  group_by(speciesId,speciesCode) %>%
  summarise(n_years = n()) %>%
  filter(n_years > 53)

## select indices for species with > 1 region
tst <- goal_indices_tbl %>%
  filter(speciesId %in% sp_gt_1_pop$speciesId)

# manually determine which region names to select for national models
tst1 <- tst %>%
  group_by(speciesId,speciesCode,areaCode) %>%
  summarise(n_years = n())

regions_sel <- c("Canada", "CAN", "USACAN")
# filter the indices to only the national-assessment regions
tst <- tst %>%
  filter(areaCode %in% regions_sel)

## test to make sure manual process worked
tst1 <- tst %>%
  group_by(speciesId,speciesCode,areaCode) %>%
  summarise(n_years = n())
if(max(tst1$n_years) > 53){
  stop("Some species have too many years of data")
}
## end test


# select indices for species with only one region
tst2 <- goal_indices_tbl %>%
  filter(!speciesId %in% sp_gt_1_pop$speciesId)

# combine indices for species with one region with national scale indices for species with >1 region
all_inds <- bind_rows(tst2,tst)


# final check that each species has only 1 time-series
n_yrs <- all_inds %>%
  group_by(speciesId,areaCode) %>%
  summarise(n_years = n())
if(max(n_yrs$n_years) > 53){
  stop("Some species have too many years of data")
}

miss_inds <- all_inds %>%
  filter(index <= 0 |
         indexLowerCI <= 0)

all_inds <- all_inds %>%
  filter(index > 0, indexLowerCI >= 0)

# Generate species-level smooths ------------------------------------------


#all_inds <- readRDS("data/Canadian_BBS_indices.rds")

#all_inds <- readRDS("socb_smoothed_indices.rds")

all_inds <- all_inds %>%
  mutate(ln_index = log(index),
         ln_lci = log(indexLowerCI),
         ln_uci = log(indexUpperCI),
         ln_index_sd = (ln_uci-ln_lci)/3.9) %>%
  group_by(speciesId) %>%
  mutate(yearn = year-(min(year)-1)) %>% # sets a yearn value specific to each species
  arrange(speciesId,year)

all_sp <- unique(all_inds$speciesId)
all_smoothed_indices <- NULL

for(species in all_sp){

  inds <- all_inds %>%
    filter(speciesId == species) %>%
    select(speciesId,year,
           ln_index,ln_index_sd,yearn) %>%
    distinct()

  years <- c(min(inds$year):max(inds$year))
  n_years <- as.integer(length(years))
  n_indices <- as.integer(nrow(inds))
  n_knots <- as.integer(round(n_indices/3))
  start_year <- min(inds$year)

  gam_data <- gam_basis(years,
                        nknots = n_knots,
                        sm_name = "year")
  stan_data <- list(
    n_years = n_years,
    n_indices = n_indices,
    n_knots_year = gam_data$nknots_year,
    year = inds$yearn,
    ln_index = inds$ln_index,
    ln_index_sd = inds$ln_index_sd,
    year_basis = gam_data$year_basis
  )
  ## fit model with cmdstanr
  file <- "models/GAM_smooth_model.stan"
  mod <- cmdstan_model(file)

  fit <- mod$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 0,
                    adapt_delta = 0.95,
                    show_messages = FALSE,
                    show_exceptions = FALSE)

  sum <- fit$summary(variables = NULL,
                     "mean",
                     "sd",
                     "ess_bulk",
                     "rhat",
                     q2_5 = q2_5,
                     q97_5 = q97_5)

  mx_rhat <- max(sum$rhat,na.rm = TRUE)
  if(mx_rhat > 1.05){
    fit <- mod$sample(data = stan_data,
                      parallel_chains = 4,
                      iter_warmup = 8000,
                      iter_sampling = 8000,
                      thin = 8,
                      refresh = 0)
    sum <- fit$summary(variables = NULL,
                       "mean",
                       "sd",
                       "ess_bulk",
                       "rhat",
                       q2_5 = q2_5,
                       q97_5 = q97_5)
  }
  smooth_inds <- sum %>%
    filter(grepl("mu_smooth",variable)) %>%
    mutate(yearn = as.integer(str_extract_all(variable,
                                              "[[:digit:]]{1,}",
                                              simplify = TRUE)),
           speciesId = species,
           smooth_ind = mean,
           smooth_ind_sd = sd,
           smooth_ind_lci = q2_5,
           smooth_ind_uci = q97_5) %>%
    select(speciesId,yearn,smooth_ind,smooth_ind_sd,
           smooth_ind_lci,smooth_ind_uci)

  diffs <- sum %>%
    filter(grepl("annual_diffs",variable)) %>%
    mutate(yearn = as.integer(str_extract_all(variable,
                                              "[[:digit:]]{1,}",
                                              simplify = TRUE)),
           annual_diff = mean,
           annual_diff_sd = sd) %>%
    select(yearn,annual_diff,annual_diff_sd)

  scaled_status <- sum %>%
    filter(grepl("scaled_status",variable)) %>%
    mutate(yearn = as.integer(str_extract_all(variable,
                                              "[[:digit:]]{1,}",
                                              simplify = TRUE)),
           scaled_status = mean,
           scaled_status_sd = sd) %>%
    select(yearn,scaled_status,scaled_status_sd)

  scaled_log_status <- sum %>%
    filter(grepl("scaled_log_status",variable)) %>%
    mutate(yearn = as.integer(str_extract_all(variable,
                                              "[[:digit:]]{1,}",
                                              simplify = TRUE)),
           scaled_log_status = (mean),
           scaled_log_status_sd = sd) %>%
    select(yearn,scaled_log_status,scaled_log_status_sd)


  smooth_inds <- smooth_inds %>%
    inner_join(.,diffs,
               by = c("yearn"))%>%
    inner_join(.,scaled_status,
               by = c("yearn"))%>%
    inner_join(.,scaled_log_status,
               by = c("yearn")) %>%
    mutate(year = yearn+(start_year-1))

  all_smoothed_indices <- bind_rows(all_smoothed_indices,smooth_inds)

  tst <- ggplot(data = smooth_inds,
                aes(x = year, y = annual_diff))+
    geom_line()
  tst

  print(paste(species,"complete",round(which(all_sp == species)/length(all_sp),2)))

}

saveRDS(all_smoothed_indices,"socb_smoothed_indices.rds")

# Prepare data ------------------------------------------------------------


all_smoothed_indices <- readRDS("socb_smoothed_indices.rds")




# Group-level models ------------------------------------------------------

my_custom_name_repair <- function(nms){
  nc <- tolower(str_replace_all( pattern = "[[:punct:]]+",replacement = "", nms))
  nc <- str_replace_all( pattern = "[\r\n]",replacement = "", nc)
  nc <- str_replace_all( pattern = "[[:blank:]]+",replacement = "_", nc)
}

species_groups <- read_xlsx("data/SoCB species guild associations Feb2024.xlsx",
                            .name_repair = my_custom_name_repair)

groups_to_fit <- species_groups %>%
  select(waterfowl:all_other_birds_tous_les_autres_oiseaux) %>%
  names()

# remove non-native and range expansion species
species_to_fit <- species_groups %>%
  rename(speciesId = speciesid) %>%
  mutate(speciesId = as.integer(speciesId)) %>%
  filter(!is.na(regularly_occurring_native_species_included_in_analyses),
         is.na(species_whose_range_expanded_into_canada_since_1970_excluded_from_indicators),
         !is.na(speciesId))



# Group loop --------------------------------------------------------------

pdf("figures/composite_summary_plots.pdf",
    width = 8.5,
    height = 11)
for(grp in groups_to_fit){

  species_sel <- species_to_fit[which(!is.na(species_to_fit[,grp])),c("speciesId",
                                                                      "english_name",grp)]

  if(nrow(species_sel) == 0){next}

  n_sp_w_data <- length(which(species_sel$speciesId %in% all_smoothed_indices$speciesId))

  print(paste("There are data for",n_sp_w_data,"of",nrow(species_sel),
              "in the",grp,"group"))

if(n_sp_w_data/nrow(species_sel) < 0.5){
  print(paste("Skipping",grp,"because only",round(n_sp_w_data/nrow(species_sel),2)*100,"% of species have data"))
    next}
  ### Drop the base-year values and other species
inds_all <- all_smoothed_indices %>%
  inner_join(.,species_sel,
             by = "speciesId") %>%
  mutate(species_ind = as.integer(factor(speciesId)))#

base_yr <- max(base_year,min(inds_all$year))

inds <- inds_all %>%
  group_by(speciesId,species_ind) %>%
  mutate(yearn2 = year-base_year) %>%
  filter(yearn2 > 0,
         scaled_log_status_sd > 0,
         scaled_status_sd > 0,
         annual_diff_sd > 0) %>%
  arrange(species_ind,yearn2)

## track the start and end years for each species
sp_y <- inds %>%
  group_by(species_ind,speciesId) %>%
  summarise(first_year = min(year),
            last_year = max(year),
            first_yearn2 = min(yearn2),
            last_yearn2 = max(yearn2),
            .groups = "drop")

# number of years and species
n_years <- max(inds$yearn2,na.rm = TRUE)
n_species <- max(inds$species_ind,na.rm = TRUE)


species <- matrix(data = NA,
                  nrow = n_years,
                  ncol = n_species)
ln_index <- matrix(data = 0,
                   nrow = n_years,
                   ncol = n_species)
ln_index_sd <- matrix(data = 0,
                      nrow = n_years,
                      ncol = n_species)
n_species_year <- vector("integer",length = n_years)
for(y in 1:n_years){
  tmp <- inds %>%
    filter(yearn2 == y) %>%
    arrange(species_ind)
  n_sp <- nrow(tmp)
  sp_inc <- as.integer(tmp$species_ind)
  species[y,1:n_sp] <- sp_inc
  if(n_sp < n_species){
    sp_miss <- c(1:n_species)[-sp_inc]
    species[y,c((n_sp+1):n_species)] <- sp_miss
  }
  n_species_year[y] <- n_sp
  for(s in sp_inc){
    ln_index[y,s] <- as.numeric(tmp[which(tmp$species_ind == s & tmp$yearn2 == y),"annual_diff"])
    ln_index_sd[y,s] <- as.numeric(tmp[which(tmp$species_ind == s & tmp$yearn2 == y),"annual_diff_sd"])
  }
}

stan_data2 <- list(n_years = n_years,
                   n_species = n_species,
                   n_species_year = n_species_year,
                   species = species,
                   ln_index = ln_index,
                   ln_index_sd = ln_index_sd)


file2 <- "models/State_of_Birds_Model_differences.stan"
mod2 <- cmdstan_model(file2)

fit2 <- mod2$sample(data = stan_data2,
                    parallel_chains = 4,
                    refresh = 0,
                    adapt_delta = 0.8)

sum2 <- fit2$summary(variables = NULL,
                     "mean",
                     "sd",
                     "ess_bulk",
                     "rhat",
                     q2_5 = q2_5,
                     q97_5 = q97_5)

mx_rhat2 <- max(sum2$rhat,na.rm = TRUE)
if(mx_rhat2 > 1.05){
  fit2 <- mod$sample(data = stan_data2,
                     parallel_chains = 4,
                     iter_warmup = 8000,
                     iter_sampling = 8000,
                     thin = 8,
                     refresh = 0)
  sum2 <- fit2$summary(variables = NULL,
                       "mean",
                       "sd",
                       "ess_bulk",
                       "rhat",
                       q2_5 = q2_5,
                       q97_5 = q97_5)
}

annual_status_difference <- sum2 %>%
  filter(grepl("annual_status",variable)) %>%
  mutate(yearn2 = as.integer(str_extract_all(variable,
                                             "[[:digit:]]{1,}",
                                             simplify = TRUE)),
         model = "difference",
         year = yearn2+(base_yr-1),
         percent_diff = (exp(mean)-1)*100,
         percent_diff_lci = (exp(q2_5)-1)*100,
         percent_diff_uci = (exp(q97_5)-1)*100)

sigma1 <- sum %>%
  filter(grepl("sigma",variable))
sigma2 <- sum2 %>%
  filter(grepl("sigma",variable))

  saveRDS(annual_status_difference,paste0("output/composite_fit_",grp,".rds"))
  saveRDS(inds_all,paste0("output/composite_data_",grp,".rds"))


# alternate original model ------------------------------------------------




  ## fill the data matrices and vectors
  for(y in 1:n_years){
    tmp <- inds %>%
      filter(yearn2 == y) %>%
      arrange(species_ind)
    n_sp <- nrow(tmp)
    sp_inc <- as.integer(tmp$species_ind)
    species[y,1:n_sp] <- sp_inc
    if(n_sp < n_species){
      sp_miss <- c(1:n_species)[-sp_inc]
      species[y,c((n_sp+1):n_species)] <- sp_miss
    }
    n_species_year[y] <- n_sp
    for(s in sp_inc){
      ln_index[y,s] <- as.numeric(tmp[which(tmp$species_ind == s & tmp$yearn2 == y),"scaled_status"])
      ln_index_sd[y,s] <- as.numeric(tmp[which(tmp$species_ind == s & tmp$yearn2 == y),"scaled_status_sd"])
    }
  }

  stan_data <- list(n_years = n_years,
                    n_species = n_species,
                    n_species_year = n_species_year,
                    species = species,
                    ln_index = ln_index,
                    ln_index_sd = ln_index_sd)


  ## fit Standard Stan model
  file <- "models/State_of_Birds_Model_standard.stan"
  mod <- cmdstan_model(file)

  fit <- mod$sample(data = stan_data,
                    parallel_chains = 4,
                    refresh = 0,
                    adapt_delta = 0.8)

  sum <- fit$summary(variables = NULL,
                     "mean",
                     "sd",
                     "ess_bulk",
                     "rhat",
                     q2_5 = q2_5,
                     q97_5 = q97_5)

  mx_rhat <- max(sum$rhat,na.rm = TRUE)
  if(mx_rhat > 1.05){ # if not converged run again with more iterations
    fit <- mod$sample(data = stan_data,
                      parallel_chains = 4,
                      iter_warmup = 8000,
                      iter_sampling = 8000,
                      thin = 8,
                      refresh = 0)
    sum <- fit$summary(variables = NULL,
                       "mean",
                       "sd",
                       "ess_bulk",
                       "rhat",
                       q2_5 = q2_5,
                       q97_5 = q97_5)
  }

  annual_status_standard <- sum %>%
    filter(grepl("annual_status",variable)) %>%
    mutate(yearn2 = as.integer(str_extract_all(variable,
                                               "[[:digit:]]{1,}",
                                               simplify = TRUE)),
           model = "standard",
           year = yearn2+(base_yr-1),
           percent_diff = (exp(mean)-1)*100,
           percent_diff_lci = (exp(q2_5)-1)*100,
           percent_diff_uci = (exp(q97_5)-1)*100)

  saveRDS(annual_status_standard,paste0("output/composite_fit_standard_",grp,".rds"))

  annual_status_combine <- bind_rows(annual_status_difference,
                                     annual_status_standard)


  brks_pch <- c(-95,-90,-75,-50,-40,-25,0,25,50,100,200,300,500)
  brks_log <- log((brks_pch/100)+1)
  brks_labs <- paste0(brks_pch,"%")



  inds_label <- inds_all %>%
    inner_join(.,sp_y,
               by = c("species_ind",
                      "speciesId",
                      "year" = "last_year"))

  tst <- ggplot(data = annual_status_combine,
                aes(x = year,y = mean))+
    # geom_line(data = inds_all,
    #           aes(x = year,y = scaled_status,
    #               group = species_ind),
    #           alpha = 0.2,
    #           inherit.aes = FALSE)+
    geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
                alpha = 0.5)+
    geom_line()+
    # geom_text(data = inds_label,
    #           aes(x = year,y = scaled_status,
    #               label = english_name),
    #           size = 2)+
    scale_y_continuous(breaks = brks_log,
                       labels = brks_labs)+
    labs(title = grp)+
    theme_bw()+
    facet_wrap(vars(model))

  tst2 <- ggplot(data = annual_status_difference,
                aes(x = year,y = mean))+
    geom_line(data = inds_all,
              aes(x = year,y = scaled_status,
                  group = species_ind),
              alpha = 0.2,
              inherit.aes = FALSE)+
    geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
                alpha = 0.5)+
    geom_line()+
    geom_text(data = inds_label,
              aes(x = year,y = scaled_status,
                  label = english_name),
              size = 2)+
    scale_y_continuous(breaks = brks_log,
                       labels = brks_labs)+
    labs(title = grp)+
    theme_bw()

  print(tst/
          tst2)



}


dev.off()



