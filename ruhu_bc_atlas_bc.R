

library(naturecounts)
library(tidyverse)
# bc_atlas_pc <- nc_data_dl(collections = "BCATLAS1PC",
#                           username = "adam.smith",
#                           info = "Rufous Hummingbird project")

#saveRDS(bc_atlas_pc,"bc_atlas_pc.rds")

bc_atlas_pc <- readRDS("bc_atlas_pc.rds")

#species is for ruhu
ruhu_n <- search_species("Rufous Hummingbird") %>%
  filter(!grepl("hybrid",english_name))

#only pc with ruhu observations
bc_atlas_pres <- bc_atlas_pc %>%
  filter(species_id == ruhu_n$species_id) %>%
  select(SamplingEventIdentifier,ObservationCount,species_id) %>%
  filter(ObservationCount > 0) %>%
  rename(Count_RUFU = ObservationCount)

# list of all unique sampling event locations and times
# join to sampling events with rufu counts
# fill in zero-counts (absences)
pc_atlas_pc_rufu <- bc_atlas_pc %>%
  select(SamplingEventIdentifier,latitude,longitude,bcr,survey_year,survey_month,survey_day,
         TimeCollected,DurationInHours) %>%
  distinct() %>%
  left_join(.,bc_atlas_pres,
            by = "SamplingEventIdentifier") %>%
  mutate(Count_RUFU = ifelse(!is.na(Count_RUFU),Count_RUFU,0)) #zero-fill

write_csv(pc_atlas_pc_rufu,"BC_atlas_rufu_counts_absences.csv")




