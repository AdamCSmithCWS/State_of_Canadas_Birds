### fitting the composite group indicator models for the State of Canada's Birds
###
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
fit_alternate_socb_model <- FALSE #use TRUE if traditional SOCB model is desired
re_smooth <- TRUE



# extract state of canada's birds -----------------------------------------
re_download <- TRUE

specid_rename <- function(x){
  y <- vector("character",length(x))
  for(i in 1:length(x)){
    if(grepl("speciesId",x[i]) |
       grepl("species_id",x[i])){
      y[i] <- "speciesID"
    }else{
      y[i] <- x[i]
    }
  }
  return(y)
}

if(re_download){

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbTrendRank",
                 timeout = 120) -> rank_tbl
  saveRDS(rank_tbl,"data/SocbTrendRank.rds")

  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "SocbSpecies",
                 timeout = 120) -> sp_tbl
  saveRDS(sp_tbl,"data/SocbSpecies.rds")

  #Information on species' group#
  nc_query_table(username = "adam.smith",
                 table = "Groups",
                 timeout = 120) -> group_tbl
  saveRDS(group_tbl,"data/Groups.rds")

  #Data on all annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndices",
                 timeout = 240) -> indices_tbl

  saveRDS(indices_tbl,"data/TrendsIndices.rds")

  #Data on goal-based annual indices of abundance#
  nc_query_table(username = "adam.smith",
                 table = "TrendsIndicesGoals",
                 timeout = 120) -> goal_indices_tbl
  saveRDS(goal_indices_tbl,"data/TrendsIndicesGoals.rds")

  #Data on species trends #
  nc_query_table(username = "adam.smith",
                 table = "Trends",
                 timeout = 120) -> trends_tbl
  saveRDS(trends_tbl,"data/Trends.rds")

  species_names <- naturecounts::search_species() %>%
    rename_with(.,.fn = specid_rename)
  saveRDS(species_names,"data/species_names.rds")


}

  sp_tbl <- readRDS("data/SocbSpecies.rds") %>%
    rename_with(.,.fn = specid_rename)
  group_tbl <- readRDS("data/Groups.rds") %>%
    rename_with(.,.fn = specid_rename)
  trend_tbl <- readRDS("data/Trends.rds") %>%
    rename_with(.,.fn = specid_rename)
  rank_tbl <- readRDS("data/SocbTrendRank.rds") %>%
    rename_with(.,.fn = specid_rename)
  indices_tbl <- readRDS("data/TrendsIndices.rds") %>%
    rename_with(.,.fn = specid_rename)
  goal_indices_tbl <- readRDS("data/TrendsIndicesGoals.rds") %>%
    rename_with(.,.fn = specid_rename)
  species_names <- readRDS("data/species_names.rds") %>%
    rename_with(.,.fn = specid_rename)


if(re_smooth){
# Generate species-level smooths ------------------------------------------


  rank_tbl <- rank_tbl %>%
    select(trendID,goalTrend,popID,rank,speciesID,subspeciesID,trendID,
           resultsCode,popType,areaCode)%>%
    filter(goalTrend == "Y",
           popType == 1) %>%
    distinct()


sp_simple <- sp_tbl %>%
  select(speciesCode,speciesID,
         population,
         objective,
         popID,
         popType) %>% #,population,popType,popID) %>%
  filter(popType == 1) %>%
  left_join(.,species_names) %>%
  distinct()


traj_sel <- rank_tbl %>%
  inner_join(.,sp_simple)



goal_indices_tbl <- goal_indices_tbl %>%
  select(resultsCode,speciesID,
         year, index, indexUpperCI,indexLowerCI,
         areaCode)%>%
  distinct()








# trend_tbl <- trend_tbl %>%
#   select(resultsCode, speciesID,trendID, socb, rank)%>%
#   distinct()


goal_indices_tbl1 <- goal_indices_tbl %>%
  inner_join(.,traj_sel) %>%
  filter(year >= base_year) %>%
  distinct()


# gwte <- goal_indices_tbl %>%
#   filter(speciesCode == "gnwtea",
#          resultsCode == "WBPHS")
#
#
#
# tst <- ggplot(data = gwte,
#               aes(x = year,y = index))+
#   geom_pointrange(aes(ymin = indexLowerCI, ymax = indexUpperCI))+
#   scale_y_continuous(trans = "log10")+
#   geom_smooth(method = "lm")
# tst


sp_gt_1_pop <- goal_indices_tbl1 %>%
  group_by(speciesID,speciesCode) %>%
  summarise(n_years = n()) %>%
  filter(n_years > 53)

## select indices for species with > 1 region
tst <- goal_indices_tbl %>%
  filter(speciesID %in% sp_gt_1_pop$speciesID)

# manually determine which region names to select for national models
# tst1 <- tst %>%
#   group_by(speciesID,speciesCode,areaCode) %>%
#   summarise(n_years = n())
# #speciesID > 1000 &
# regions_sel <- tst1 %>%
#   filter(( areaCode %in% c("Canada", "CAN", "USACAN"))|
#          (speciesCode %in% c("snogoo") & areaCode %in% c("USACAN")) |
#            (speciesCode %in% c("gnwtea","rinduc") & areaCode %in% c("CANPRAIRIE_WBOREAL")) )
# # filter the indices to only the national-assessment regions
# tst <- tst %>%
#   inner_join(.,regions_sel,
#              by = c("speciesID","speciesCode","areaCode"))

## test to make sure manual process worked
# tst1 <- tst %>%
#   group_by(speciesID,speciesCode,areaCode) %>%
#   summarise(n_years = n())
# if(max(tst1$n_years) > 53){
#   stop("Some species have too many years of data")
# }
# ## end test


# select indices for species with only one region
# tst2 <- goal_indices_tbl %>%
#   filter(!speciesID %in% sp_gt_1_pop$speciesID)

# combine indices for species with one region with national scale indices for species with >1 region
all_inds <- goal_indices_tbl %>%
  inner_join(.,traj_sel) %>%
  filter(year >= base_year) %>%
  distinct()


# final check that each species has only 1 time-series
n_yrs <- all_inds %>%
  group_by(speciesID,areaCode) %>%
  summarise(n_years = n())
if(max(n_yrs$n_years) > 53){
  stop("Some species have too many years of data")
}

miss_inds <- all_inds %>%
  filter(index <= 0 |
         indexLowerCI <= 0 |
           is.infinite(indexUpperCI))
if(nrow(miss_inds) > 0){
  miss_sp <- unique(miss_inds$speciesCode)
  all_inds <- all_inds %>%
    rowwise() %>%
    mutate(ifelse(indexLowerCI == 0,max(0.1,index-indexUpperCI),indexLowerCI)) #fixes lower bound for a handful of
  warning(paste((paste(miss_sp,collapse = ", ")),"had missing index or CI information, that should be repaired in the next step"))
}

all_inds <- all_inds %>%
  filter(index > 0, indexLowerCI > 0)

saveRDS(all_inds,"data/all_socb_goal_indices.rds")
# Generate species-level smooths ------------------------------------------


#all_inds <- readRDS("data/Canadian_BBS_indices.rds")

#all_inds <- readRDS("data/all_socb_goal_indices.rds")

all_inds <- all_inds %>%
  mutate(ln_index = log(index),
         ln_lci = log(indexLowerCI),
         ln_uci = log(indexUpperCI),
         ln_index_sd = (ln_uci-ln_lci)/3.9) %>%
  group_by(speciesID) %>%
  mutate(yearn = year-(min(year)-1)) %>% # sets a yearn value specific to each species
  arrange(speciesID,year)

all_sp <- unique(all_inds$speciesID)
all_smoothed_indices <- NULL

for(species in all_sp){

  inds <- all_inds %>%
    filter(speciesID == species) %>%
    select(speciesID,year,
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
           speciesID = species,
           smooth_ind = mean,
           smooth_ind_sd = sd,
           smooth_ind_lci = q2_5,
           smooth_ind_uci = q97_5) %>%
    select(speciesID,yearn,smooth_ind,smooth_ind_sd,
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

}
# Prepare data ------------------------------------------------------------


all_smoothed_indices <- readRDS("socb_smoothed_indices.rds") %>%
  rename_with(.,
              .fn = specid_rename) %>%
  left_join(.,species_names,
            by = "speciesID")




# Group-level models ------------------------------------------------------

if(re_download){
# Species groups
species_groups <- nc_query_table(table = "SocbTrendGroups",
                                 username = "adam.smith") %>%
  rename_with(., .fn = specid_rename) %>%
  mutate(groupName = gsub("\r\n",
                          "",
                          x = groupName,
                          fixed = TRUE),
         groupNameFr = gsub("\r\n",
                          "",
                          x = groupNameFr,
                          fixed = TRUE))

saveRDS(species_groups,"data/species_groups.rds")
}else{
  species_groups <- readRDS("data/species_groups.rds")
}

species_non_native <- data.frame(english_name =
                                c("Mute Swan",
                                  "Mountain Quail",
                                  "California Quail",
                                  "Gray Partridge",
                                  "Ring-necked Pheasant",
                                  "Chukar",
                                  "Rock Pigeon",
                                  "Eurasian Collared-Dove",
                                  "Eurasian Skylark",
                                  "European Starling",
                                  "House Sparrow"))


ss_track <- NULL
for(j in 1:nrow(species_non_native)){
  ss <- naturecounts::search_species(as.character(species_non_native[j,"english_name"]))

  if(nrow(ss) > 1){
    ss_sel <- which(ss$english_name == as.character(species_non_native[j,"english_name"]))
  }else{
    ss_sel <- 1
  }
  species_non_native[j,"speciesID"] <- ss[ss_sel,"species_id"]
  species_non_native[j,"english_name_nc"] <- ss[ss_sel,"english_name"]

  ss_track <- bind_rows(ss_track,ss)
}


# # remove non-native and range expansion species

species_to_drop <- data.frame(english_name =
                                c("Wild Turkey",
"Anna's Hummingbird",
"Black-necked Stilt",
"Great Egret",
"White-faced Ibis",
"Red-bellied Woodpecker",
"Bushtit",
"Carolina Wren",
"Blue-gray Gnatcatcher",
"Blue-winged Warbler",
"Gray Flycatcher"))

ss_track <- NULL
for(j in 1:nrow(species_to_drop)){
  ss <- naturecounts::search_species(as.character(species_to_drop[j,"english_name"]))

  if(nrow(ss) > 1){
    ss_sel <- which(ss$english_name == as.character(species_to_drop[j,"english_name"]))
  }else{
    ss_sel <- 1
  }
  species_to_drop[j,"speciesID"] <- ss[ss_sel,"species_id"]
  species_to_drop[j,"english_name"] <- ss[ss_sel,"english_name"]

  ss_track <- bind_rows(ss_track,ss)
}


species_to_drop

species_groups <- species_groups %>%
  mutate(include = ifelse(speciesID %in% species_non_native$speciesID,
                          "N",include))

species_groups <- species_groups %>%
  filter(popType == 1) %>%
  mutate(subgroupID = ifelse(groupName == "Wetland Birds: All",
                             0,subgroupID))
# species_groups <- read_xlsx("data/SOCB_AnalysisGroups.xlsx") %>%
#   mutate(groupName = gsub("/",x = groupName,
#                            replacement = "_",fixed = TRUE)) #remvoing the special character in "Edge/Early"
#

groups_to_fit <- unique(species_groups$groupName)

# Group loop --------------------------------------------------------------
composite_plots <- vector('list',length(groups_to_fit))
names(composite_plots) <- groups_to_fit
annual_status_combine <- NULL

re_fit <- TRUE


pdf("figures/composite_summary_plots.pdf",
    width = 8.5,
    height = 11)
for(grp in groups_to_fit){

  sub_groups_to_fit <- species_groups %>%
    filter(groupName == grp,
           !speciesID %in% species_to_drop$speciesID) %>%
    select(subgroupID) %>%
    distinct() %>%
    #filter(subgroupID!= "NULL") %>%
    unlist()

  for(sub_grp in sub_groups_to_fit){
    grp_labl <- species_groups %>%
      filter(groupName == grp,
             subgroupID == sub_grp) %>%
      select(groupName,groupNameFr) %>%
      distinct() %>%
      unlist()
    grp_labl <- gsub("[[:punct:]]","",x = grp_labl)
    grp_labl <- gsub("[[:blank:]]","_",x = grp_labl)
    grp_labl <- paste(grp_labl,collapse = "-")



  species_sel <- species_groups %>%
    filter(groupName == grp,
           subgroupID == sub_grp,
           include == "Y",
           !speciesID %in% species_to_drop$speciesID) %>%
    distinct()



  if(nrow(species_sel) == 0){next}
  if(nrow(species_sel) != length(unique(species_sel$speciesID))){
    stop(paste("Problem with species list for",grp,sub_grp))
  }
  n_sp_w_data <- length(which(species_sel$speciesID %in% all_smoothed_indices$speciesID))

  print(paste("There are data for",n_sp_w_data,"of",nrow(species_sel),
              "in the",grp,sub_grp,"group"))

  species_sel <- species_sel %>%
    mutate(sufficient_data = ifelse(speciesID %in% all_smoothed_indices$speciesID,
                                    TRUE,
                                    FALSE)) %>%
    left_join(.,species_names,
              by = "speciesID") %>%
    select(-taxon_group)

if(n_sp_w_data/nrow(species_sel) < 0.2){
  print(paste("Skipping",grp,sub_grp,"because only",round(n_sp_w_data/nrow(species_sel),2)*100,"% of species have data"))
    next}

  ### Drop the base-year values and other species
inds_all <- all_smoothed_indices %>%
  filter(speciesID %in% species_sel$speciesID) %>%
  mutate(species_ind = as.integer(factor(speciesID)))#

years_w_GT_50 <- inds_all %>%
  group_by(year) %>%
  summarise(n_sp = n(),.groups = "drop") %>%
  mutate(p_sp = n_sp/max(n_sp,na.rm = TRUE)) %>%
  filter(p_sp > 0.5)
min_yr <- range(years_w_GT_50$year)[1]
max_yr <- range(years_w_GT_50$year)[2]


base_yr <- max(base_year,min_yr)

inds <- inds_all %>%
  group_by(speciesID,species_ind) %>%
  mutate(yearn2 = year-base_yr) %>%
  filter(year <= max_yr,
         year > base_yr,
         scaled_log_status_sd > 0,
         scaled_status_sd > 0,
         annual_diff_sd > 0) %>%
  arrange(species_ind,yearn2)

## track the start and end years for each species
sp_y <- inds %>%
  group_by(species_ind,speciesID) %>%
  summarise(first_year = min(year),
            last_year = max(year),
            first_yearn2 = min(yearn2),
            last_yearn2 = max(yearn2),
            .groups = "drop")
if(re_fit){

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
  fit2 <- mod2$sample(data = stan_data2,
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
  mx_rhat2 <- max(sum2$rhat,na.rm = TRUE)
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


sigma2 <- sum2 %>%
  filter(grepl("sigma",variable))

annual_status_difference <- annual_status_difference %>%
  mutate(groupName = grp,
         subgroupID = sub_grp)
inds_all <- inds_all %>%
  mutate(groupName = grp,
         subgroupID = sub_grp)
species_sel <- species_sel %>%
  mutate(groupName = grp,
         subgroupID = sub_grp)

  saveRDS(annual_status_difference,paste0("output/composite_fit_",grp_labl,".rds"))
  saveRDS(inds_all,paste0("output/composite_data_",grp_labl,".rds"))
  saveRDS(species_sel,paste0("output/composite_species_list_",grp_labl,".rds"))

    }else{#end if re-fit
      annual_status_difference <- readRDS(paste0("output/composite_fit_",grp_labl,".rds"))
      inds_all <- readRDS(paste0("output/composite_data_",grp_labl,".rds"))
      species_sel <- readRDS(paste0("output/composite_species_list_",grp_labl,".rds"))


    }
# alternate original model ------------------------------------------------
if(fit_alternate_socb_model){



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

  saveRDS(annual_status_standard,paste0("output/composite_fit_standard_",grp_labl,".rds"))

  annual_status_combine <- bind_rows(annual_status_difference,
                                     annual_status_standard)

}else{
  annual_status_combine <- bind_rows(annual_status_combine,annual_status_difference)
}
  brks_pch <- c(-95,-90,-75,-50,-40,-25,0,25,50,100,200,300,500)
  brks_log <- log((brks_pch/100)+1)
  brks_labs <- paste0(brks_pch,"%")

species_miss_en <- species_sel %>%
  filter(sufficient_data == FALSE) %>%
  select(english_name) %>%
  unlist()

species_miss_fr <- species_sel %>%
  filter(sufficient_data == FALSE) %>%
  select(french_name) %>%
  unlist()


capt_en <- paste("No data:",paste(species_miss_en,collapse = ", "))
capt_fr <- paste("Pas de données:",paste(species_miss_fr,collapse = ", "))

capt <- paste(capt_en,"\n",capt_fr)

ylimu <- max(annual_status_difference$q97_5,1.1)
yliml <- min(annual_status_difference$q2_5,-1.1)


inds_label <- inds_all %>%
    inner_join(.,sp_y,
               by = c("species_ind",
                      "speciesID",
                      "year" = "last_year"))

  tst <- ggplot(data = annual_status_difference,
                aes(x = year,y = mean))+
    geom_hline(yintercept = 0)+
    geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
                alpha = 0.25)+
    geom_line()+
    scale_y_continuous(breaks = brks_log,
                       labels = brks_labs,
                       limits = c(yliml,ylimu))+
    labs(title = paste(grp_labl),
         caption = capt)+
    xlab("")+
    theme_bw()+
    theme(plot.caption = element_text(size = 5))

  tst2 <- ggplot(data = annual_status_difference,
                aes(x = year,y = mean))+
    geom_hline(yintercept = 0)+
    geom_line(data = inds_all,
              aes(x = year,y = scaled_status,
                  group = species_ind),
              alpha = 0.2,
              inherit.aes = FALSE)+
    geom_ribbon(aes(ymin = q2_5,ymax = q97_5),
                alpha = 0.5)+
    geom_line()+
    ggrepel::geom_text_repel(data = inds_label,
              aes(x = year,y = scaled_status,
                  label = english_name),
              size = 1.5,
              max.overlaps = 30,
              min.segment.length = 0,
              nudge_x = 5,
              alpha = 0.8,
              box.padding = 0.1,
              segment.alpha = 0.3,
              segment.size = 0.2)+
    scale_y_continuous(breaks = brks_log,
                       labels = brks_labs)+
    scale_x_continuous(limits = c(1970,2062))+
    labs(title = paste(grp_labl))+
    xlab("")+
    theme_bw()

  #tst2

  print(tst/
          tst2)

composite_plots[[grp]] <- tst2


}

}


dev.off()
#

saveRDS(annual_status_combine,"output/annual_status_combine.rds")
