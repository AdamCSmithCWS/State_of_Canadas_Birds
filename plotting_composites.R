### combine composite indicators into single graph


library(tidyverse)
library(readxl)
library(patchwork)



base_year <- 1970


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

all_composites <- NULL

for(grp in groups_to_fit[c(8,11,12,13,14)]){

file_grp <- paste0("output/composite_fit_",grp,".rds")
if(file.exists(file_grp)){
annual_status_difference <- readRDS(file_grp) %>%
  mutate(composite_group = grp)

all_composites <- bind_rows(all_composites,annual_status_difference)

}

}

out_composites <- all_composites %>%
  select(year,mean,q2_5,q97_5,percent_diff,percent_diff_lci,percent_diff_uci,
         composite_group)
write_csv(out_composites,
          "output/saved_draft_composite_trajectories.csv")
groups_to_fit

final_years <- all_composites %>%
  group_by(composite_group) %>%
  summarise(last_year = max(year))

names_plot <- all_composites %>%
  inner_join(.,final_years,
             by = c("composite_group",
                    "year" = "last_year")) %>%
  mutate(lbl = paste(composite_group,round(percent_diff),"%"))




brks_pch <- c(-95,-90,-75,-50,-40,-25,0,25,50,100,200,300,500) # values of  percent change I’d like to show on the y-axis
brks_log <- log((brks_pch/100)+1) # above values transformed to original log-scale – used to set the breaks in the log-scale graph below.
brks_labs <- paste0(brks_pch,"%") # text labels to add to the y-axis



overview <- ggplot(data = all_composites,
                   aes(x = year,y = mean, group = composite_group,
                       colour = composite_group))+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_text(data = names_plot,
            aes(label = lbl),nudge_x = 10,
            size = 2)+
  coord_cartesian(xlim = c(1970,2040))+
  # scale_color_brewer(type = "qual",
  #                    palette = "Dark2")+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs)+
  theme(legend.position = "none")

pdf("figures/Initial_composite_trends.pdf",
    width = 11,
    height = 8.5)
print(overview)
dev.off()





# simpler plot ------------------------------------------------------------

all_composites <- read_csv("output/saved_draft_composite_trajectories.csv")




brks_pch <- c(-95,-90,-75,-50,-40,-25,0,25,50,100,200,300,500) # values of  percent change I’d like to show on the y-axis
brks_log <- log((brks_pch/100)+1) # above values transformed to original log-scale – used to set the breaks in the log-scale graph below.
brks_labs <- paste0(brks_pch,"%") # text labels to add to the y-axis



overview <- ggplot(data = all_composites,
                   aes(x = year,y = mean, group = composite_group,
                       colour = composite_group))+
  geom_hline(yintercept = 0)+
  geom_line()+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_y_continuous(breaks = brks_log,
                     labels = brks_labs)+
  theme(legend.position = "none")

pdf("figures/example_composite_trends.pdf",
    width = 11,
    height = 8.5)
print(overview)
dev.off()





