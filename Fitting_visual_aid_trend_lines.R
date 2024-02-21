### fitting the composite group indicator models for the State of Canada's Birds

library(tidyverse)
library(cmdstanr)
library(readxl)
library(naturecounts)
library(patchwork)

# Quantile functions to generate 95% CI instead of stan's default 90
q2_5 <- function(x)c(q2_5 = quantile(x,probs = c(0.025),
                                     names = FALSE,
                                     na.rm = TRUE))
q97_5 <- function(x)c(q97_5 = quantile(x,probs = c(0.975),
                                       names = FALSE,
                                       na.rm = TRUE))




# extract state of canada's birds -----------------------------------------
re_download <- FALSE

if(re_download){
  #Data on species populations#
  nc_query_table(username = "adam.smith",
                 table = "vwResultsSocbSpecies") -> sp_tbl
  saveRDS(sp_tbl,"data/vwResultsSocbSpecies.rds")

  #Information on species' group#
  nc_query_table(username = "adam.smith",
                 table = "vwResultsGroups") -> group_tbl
  saveRDS(group_tbl,"data/vwResultsGroups.rds")

  #Data on all trends#
  nc_query_table(username = "adam.smith",
                 table = "vwResultsTrendsAll") -> trend_tbl
  saveRDS(trend_tbl,"data/vwResultsTrendsAll.rds")

  #Data on species trends used for SoCB rank 1-3#
  nc_query_table(username = "adam.smith",
                 table = "vwResultsSocbTrendRank") -> rank_tbl
  saveRDS(rank_tbl,"data/vwResultsSocbTrendRank.rds")

  #Data on species annual indices of abundance #
  nc_query_table(username = "adam.smith",
                 table = "vwResultsTrendsIndices") -> indices_tbl
  saveRDS(indices_tbl,"data/vwResultsTrendsIndices.rds")
}else{
  sp_tbl <- readRDS("data/vwResultsSocbSpecies.rds")
  group_tbl <- readRDS("data/vwResultsGroups.rds")
  trend_tbl <- readRDS("data/vwResultsTrendsAll.rds")
  rank_tbl <- readRDS("data/vwResultsSocbTrendRank.rds")
  indices_tbl <- readRDS("data/vwResultsTrendsIndices.rds")

}

# Generate species-level smooths ------------------------------------------

sel_surveys <- rank_tbl %>%
  filter(
    #rank == 1,
         period %in% c("long-term","Long-term"),
         areaCode %in% c("Canada","USACAN", "NorthAmerica"),
         resultsCode == "CBC") %>%
  select(speciesID,resultsCode,period,areaCode) %>%
  distinct()

all_inds <- indices_tbl %>%
  inner_join(.,sel_surveys,
             by = c("speciesId" = "speciesID",
                    "resultsCode",
                    "period",
                    "areaCode"))

#all_inds <- readRDS("socb_smoothed_indices.rds")

all_inds <- all_inds %>%
  mutate(ln_index = log(index),
         ln_lci = log(indexLowerCI),
         ln_uci = log(indexUpperCI),
         ln_index_sd = (ln_uci-ln_lci)/3.9) %>%
  group_by(speciesId) %>%
  arrange(speciesId,year)

all_sp <- unique(all_inds$speciesId)

all_smoothed_indices <- NULL

pdf("figures/linear_display_trend_examples_CBC.pdf")
for(species in all_sp){

  inds <- all_inds %>%
    filter(speciesId == species) %>%
    select(speciesId,year,
           ln_index,ln_index_sd,ln_lci,ln_uci) %>%
    distinct()

  trnd <- all_inds %>%
    filter(speciesId == species) %>%
    select(speciesId,
           trnd,lowerCI,upperCI) %>%
    distinct() %>%
    mutate(trnd = log((trnd/100)+1),
           lowerCI = log((lowerCI/100)+1),
           upperCI = log((upperCI/100)+1))


  sp_code <- unique(sp_tbl[which(sp_tbl$speciesID == species),"speciesCode"])

  n_years <- as.integer(length(min(inds$year):max(inds$year)))
  n_indices <- as.integer(nrow(inds))
  mid_year <- floor(mean(inds$year))
  base_year <- min(inds$year)-1
  yearn <- inds$year-base_year
  year_centered <- c(min(inds$year):max(inds$year))-mid_year

  stan_data <- list(
    n_years = n_years,
    n_indices = n_indices,
    year = yearn,
    year_centered = year_centered,
    ln_index = inds$ln_index,
    ln_index_sd = inds$ln_index_sd
  )
  ## fit model with cmdstanr
  file <- "models/Log_linear_trend_smooth_model.stan"
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
    filter(grepl("smooth_inds",variable)) %>%
    mutate(year = as.integer(str_extract_all(variable,
                                              "[[:digit:]]{1,}",
                                              simplify = TRUE)) + base_year,
           speciesId = species,
           smooth_ind = mean,
           smooth_ind_sd = sd,
           smooth_ind_lci = q2_5,
           smooth_ind_uci = q97_5) %>%
    select(speciesId,year,smooth_ind,smooth_ind_sd,
           smooth_ind_lci,smooth_ind_uci)

  tmp <- smooth_inds %>%
    left_join(.,inds,
              by = "year")
  intercept <- tmp %>%
    summarise(smooth_ind = mean(smooth_ind))

  tmp <- tmp%>%
    mutate(trend_line = (trnd$trnd*(year-mid_year))+intercept$smooth_ind,
           trend_line_lci = (trnd$lowerCI*(year-mid_year))+intercept$smooth_ind,
           trend_line_uci = (trnd$upperCI*(year-mid_year))+intercept$smooth_ind)


   preds_i <- NULL
  for(i in 1:500){
    i_sel <- data.frame(y = rnorm(nrow(inds),
                   inds$ln_index,
                   inds$ln_index_sd),
                   year = inds$year)
    m1 <- lm(y~year,data = i_sel)
    #preds_i[,paste0("i_",i)] <- predict(m1,se.fit = FALSE)
    p_tmp <- data.frame(pred = as.numeric(predict(m1,se.fit = FALSE)),
                        year = inds$year,
                        it = i)
    preds_i <- bind_rows(preds_i,p_tmp)

  }

  preds_i_sum <- preds_i %>%
    group_by(year) %>%
    summarise(pred_m = mean(pred),
              pred_lci = quantile(pred,0.025,names = FALSE),
              pred_uci = quantile(pred,0.975,names = FALSE))

  # tst_gg <- ggplot(data = preds_i,
  #                  aes(x = year, y = pred, colour = it, group = it))+
  #   geom_line(alpha = 0.4)+
  #   geom_pointrange(data = inds,
  #                   aes(x = year,ymin = ln_lci,
  #                       ymax = ln_uci,
  #                       y = ln_index),
  #                   inherit.aes = FALSE)+
  #   scale_colour_viridis_c()
  # tst_gg

  upy1 <- min(max(tmp$ln_index+(0.1*max(tmp$ln_index))),
             max(tmp$ln_uci))

  upy <- min(max(tmp$ln_index+(1.1*max(tmp$ln_index))),
             max(tmp$ln_uci))

  tst <- ggplot(data = tmp)+
    geom_ribbon(data = preds_i_sum,
              aes(x = year,ymin = pred_lci,
                  ymax = pred_uci),
              alpha = 0.3,
              fill = "darkorange")+
    geom_line(data = preds_i_sum,
              aes(x = year,y = pred_m),
              colour = "darkorange")+
    geom_ribbon(aes(x = year, ymin = smooth_ind_lci,
                    ymax = smooth_ind_uci),
                alpha = 0.3)+
    geom_pointrange(aes(x = year,ymin = ln_lci,
                        ymax = ln_uci,
                        y = ln_index))+
    geom_line(aes(x = year, y = smooth_ind))+
    # geom_line(data = preds_i,
    #           aes(x = year,y = pred,group = it),
    #           alpha = 0.1,
    #           colour = "blue")+
    geom_line(aes(x = year,y = trend_line),colour = "red",
              linetype = 2)+
    geom_line(aes(x = year,y = trend_line_lci),colour = "red",
              linetype = 3)+
    geom_line(aes(x = year,y = trend_line_uci),colour = "red",
              linetype = 3)+
    labs(title = sp_code)+
    coord_cartesian(ylim = c(NA,upy1))+
    scale_y_continuous(trans = "exp")+
    theme_bw()

  tst2 <- ggplot(data = tmp)+
    geom_ribbon(data = preds_i_sum,
                aes(x = year,ymin = pred_lci,
                    ymax = pred_uci),
                alpha = 0.3,
                fill = "darkorange")+
    geom_line(data = preds_i_sum,
              aes(x = year,y = pred_m),
              colour = "darkorange")+
    geom_ribbon(aes(x = year, ymin = smooth_ind_lci,
                    ymax = smooth_ind_uci),
                alpha = 0.3)+
    geom_pointrange(aes(x = year,ymin = ln_lci,
                        ymax = ln_uci,
                        y = ln_index))+
    geom_line(aes(x = year, y = smooth_ind))+
    geom_line(aes(x = year,y = trend_line),colour = "red",
              linetype = 2)+
    geom_line(aes(x = year,y = trend_line_lci),colour = "red",
              linetype = 3)+
    geom_line(aes(x = year,y = trend_line_uci),colour = "red",
              linetype = 3)+
    coord_cartesian(ylim = c(NA,upy))+
    labs(title = sp_code)+
    theme_bw()

  print(tst / tst2)

  all_smoothed_indices <- bind_rows(all_smoothed_indices,smooth_inds)


  print(paste(species,"complete",round(which(all_sp == species)/length(all_sp),2)))

}
dev.off()

saveRDS(all_smoothed_indices,"cbc_linear_display_indices.rds")

sp_tbl_sel <- sp_tbl %>%
  select("speciesID","speciesCode") %>%
  distinct()


all_smooth_out <- all_inds %>%
  inner_join(.,all_smoothed_indices,
             by = c("year","speciesId")) %>%
  inner_join(.,sp_tbl_sel,
             by = c("speciesId" = "speciesID"))



