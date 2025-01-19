library(tidyverse)
library(readr)
library(data.table)
library(vegan)
library(plotly)
library(RColorBrewer)
library(emmeans)
library(mgcv)
library(ggplot2)
library(ggpubr)
library(zoo)
library(gridExtra)
library(gratia)
library(corrplot)
library(mgcViz)
library(viridis)

seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered <- read_csv("WP2A_ASV_table_long_filtered_family.csv")

#aggregating replicates
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=mean(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads)

#REMOVE CLEAR FISH THAT ARE LIEKLY CONSUMED BY HUAMNS AND NOT NATIVE
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Gadus_macrocephalus"),]#Pacific cod
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Hippoglossus_stenolepis"),]#Pacific halibut
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Ammodytes_personatus"),]#Pacific Sand Lance 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Clupea_harengus"),]#Atlantic herring 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Clupea_pallasii"),]#Pacific herring 
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Copadichromis_virginalis"),]#Haplochromine cichlid  
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Favonigobius_gymnauchen"),]#Sand goby   
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Pholis_crassispina"),]#Rock gunnel  
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Melanogrammus_aeglefinus"),]#Haddock
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Platichthys_stellatus"),]#Starry flounder
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species=="Scomber_scombrus"),]#Atlantic mackerel - POSITIVE CONTROL
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salmo_obtusirostris','Salmo_trutta',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Anguilla_rostrata','Anguilla_anguilla',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species<-sub('Salvelinus_fontinalis_x_Salvelinus_malma','Salvelinus_malma',seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$species)

#Remove negative controls
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T10"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T17"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T18"),]
seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered[!(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered$SampleSite_time=="ENC_T19"),]

#CONVERTING FILTERED LONG FORMAT INTO WIDE FORMAT
seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered
seqtabNoC_WP2A_12S_tax_wide_filtered$read_filter <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Days <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Date_sampled <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$SampleSite_code <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$WP <- NULL 
seqtabNoC_WP2A_12S_tax_wide_filtered$Seq_length <- NULL 

seqtabNoC_WP2A_12S_tax_wide_filtered<-seqtabNoC_WP2A_12S_tax_wide_filtered %>% 
  group_by(SampleSite_time,species) %>% 
  summarize(reads=sum(reads)) %>%
  rename(SampleSite_time=SampleSite_time,species=species, reads= reads,)

seqtabNoC_WP2A_12S_tax_wide_filtered <- spread(seqtabNoC_WP2A_12S_tax_wide_filtered,SampleSite_time,reads)
seqtabNoC_WP2A_12S_tax_wide_filtered[is.na(seqtabNoC_WP2A_12S_tax_wide_filtered)] <- 0
seqtabNoC_WP2A_12S_tax_wide_filtered <- data.frame(t(seqtabNoC_WP2A_12S_tax_wide_filtered))

names(seqtabNoC_WP2A_12S_tax_wide_filtered) <- as.matrix(seqtabNoC_WP2A_12S_tax_wide_filtered[1, ])
seqtabNoC_WP2A_12S_tax_wide_filtered <- seqtabNoC_WP2A_12S_tax_wide_filtered[-1, ]
seqtabNoC_WP2A_12S_tax_wide_filtered[] <- lapply(seqtabNoC_WP2A_12S_tax_wide_filtered, function(x) type.convert(as.character(x)))

##SPLITTING BY SPECIES
Salmo_salar<-seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered %>% filter(species == "Salmo_salar")

#NORMALISING READS
seqtabNoC_WP2A_12S_tax_col_sum2<-rowSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
WP2A_12S_sps_read_summary<-data.frame(WP2A_12S_sps_read_summary)

seqtabNoC_WP2A_12S_tax_col_sum2<-data.frame(seqtabNoC_WP2A_12S_tax_col_sum2)
seqtabNoC_WP2A_12S_tax_col_sum2$SampleSite_time <- rownames(seqtabNoC_WP2A_12S_tax_col_sum2)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")

#GETTING THE SUM OF EACH SAMPLE FOR NORMALISING READS
seqtabNoC_WP2A_12S_tax_col_sum3<-colSums(seqtabNoC_WP2A_12S_tax_wide_filtered[c(1:16)])
seqtabNoC_WP2A_12S_tax_col_sum3<-data.frame(seqtabNoC_WP2A_12S_tax_col_sum3)

seqtabNoC_WP2A_12S_tax_col_sum3$species <- rownames(seqtabNoC_WP2A_12S_tax_col_sum3)

seqtabNoC_WP2A_12S_tax_normalised<- merge(seqtabNoC_WP2A_12S_tax_col_sum2[, c("SampleSite_time", "seqtabNoC_WP2A_12S_tax_col_sum2")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="SampleSite_time")
seqtabNoC_WP2A_12S_tax_top<- merge(seqtabNoC_WP2A_12S_tax_col_sum3[, c("species", "seqtabNoC_WP2A_12S_tax_col_sum3")],seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered, by="species")

seqtabNoC_WP2A_12S_tax_normalised <- transform(seqtabNoC_WP2A_12S_tax_normalised, normalised_reads = reads / seqtabNoC_WP2A_12S_tax_col_sum2)

ggplot(seqtabNoC_WP2A_12S_tax_normalised , aes(x = SampleSite_time, fill = species, y = normalised_reads)) + 
  geom_bar(stat = "identity", colour = "black")+
  theme(axis.text.x = element_text(angle = 90),legend.position = "none")

lofresh_metadata_WP2_3 <- read_csv("lofresh_metadata_WP2&3.csv")
lofresh_metadata_WP2_3_summary = lofresh_metadata_WP2_3[!duplicated(lofresh_metadata_WP2_3$SampleSite_time),]

e12_values <- lofresh_metadata_WP2_3_summary %>%
  filter(SampleSite_code == "E12") %>%
  select(Date_sampled, daily_flow_by_catchment)

lofresh_metadata_WP2_3_summary <- lofresh_metadata_WP2_3_summary %>%
  left_join(e12_values, by = "Date_sampled", suffix = c("", "_E12")) %>%
  mutate(daily_flow_by_catchment = case_when(
      SampleSite_code %in% c("E13", "E14") ~ daily_flow_by_catchment_E12,
      TRUE ~ daily_flow_by_catchment)) %>%
  select(-daily_flow_by_catchment_E12)

seqtabNoC_WP2A_12S_tax_normalised <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time", "Days", "SampleSite_code")],seqtabNoC_WP2A_12S_tax_normalised, by="SampleSite_time")

rm(seqtabNoC_WP2A_12S_tax_top)
rm(seqtabNoC_WP2A_12S_tax_long_meta_clean_filtered)
rm(seqtabNoC_WP2A_12S_tax_col_sum2)
rm(seqtabNoC_WP2A_12S_tax_col_sum3)

#######SALMON CASE STUDY######
Salmo_salar_norm<-seqtabNoC_WP2A_12S_tax_normalised %>% filter(species == "Salmo_salar")
Salmo_salar_norm <- merge(lofresh_metadata_WP2_3_summary[, c("SampleSite_time","Month.Year", "season","Dist_from_lake","Date_sampled","pH",
                                                             "conductivity_uS_cm","gran_alk_ueq_L","daily_flow_by_catchment","monthly_flow_by_catchment","seasonal_flow_by_catchment",
                                                             "week_before_flow_by_catchment","river_width_m","river_gradient_percent","rainfall_monthly_average",
                                                             "temp_monthly_average")],Salmo_salar_norm, by="SampleSite_time")
Salmo_salar_norm$Days<-as.numeric(Salmo_salar_norm$Days)
Salmo_salar_norm$Month.Year<-as.factor(Salmo_salar_norm$Month.Year)

Salmo_salar_norm$Date_sampled <- as.Date(Salmo_salar_norm$Date_sampled, format = "%d.%m.%y")

ggplot(Salmo_salar_norm, aes(x=Date_sampled, y=normalised_reads)) + 
  geom_point()+
  geom_smooth(method = "gam", formula = y ~ s(x,k=6))+
  scale_x_date(labels = function(d) format(d, "%b %Y"), 
               breaks = "1 month") +theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~SampleSite_code, ncol=2, nrow=4)+labs(x ="Date", y = "Normalised reads")

write.csv(Salmo_salar_norm,"Salmo_salar_norm.csv")

salmon_estimates <- read_csv("salmon_estimates_Dec24.csv")
salmon_estimates$Date_sampled <- as.Date(salmon_estimates$Date_sampled, format = "%d/%m/%Y")

Salmo_salar_norm_catch_compare <- merge(salmon_estimates, Salmo_salar_norm, by = "Date_sampled", all.x = TRUE)

Salmo_salar_norm_catch_compare$Days<-as.numeric(Salmo_salar_norm_catch_compare$Days)
Salmo_salar_norm_catch_compare<-distinct(Salmo_salar_norm_catch_compare)

Salmo_salar_norm_catch_compare$SampleSite_code <- as.factor(Salmo_salar_norm_catch_compare$SampleSite_code)
Salmo_salar_norm_catch_compare <- Salmo_salar_norm_catch_compare %>%
  rename(Site = SampleSite_code)


Salmo_salar_norm_catch_compare$river_section<-Salmo_salar_norm_catch_compare$Site
Salmo_salar_norm_catch_compare$river_section <- gsub("E07|E08|E09", "Upper river salmon eDNA (normalised reads)", Salmo_salar_norm_catch_compare$river_section)
Salmo_salar_norm_catch_compare$river_section <- gsub("E10", "Conwy falls salmon eDNA (normalised reads)", Salmo_salar_norm_catch_compare$river_section)
Salmo_salar_norm_catch_compare$river_section <- gsub("E11|E12|E13|E14", "Lower river salmon eDNA (normalised reads)", Salmo_salar_norm_catch_compare$river_section)
desired_order <- c("Upper river salmon eDNA (normalised reads)", "Conwy falls salmon eDNA (normalised reads)", "Lower river salmon eDNA (normalised reads)")
Salmo_salar_norm_catch_compare$river_section <- factor(Salmo_salar_norm_catch_compare$river_section, levels = desired_order)

Salmo_salar_norm_catch_compare_biomass_plot<-Salmo_salar_norm_catch_compare[!duplicated(Salmo_salar_norm_catch_compare$interp_Adult_eggsBM), ]
Salmo_salar_norm_catch_compare_biomass_plot$river_section <- "Whole river biomass (kg)"
Salmo_salar_norm_catch_compare_biomass_plot$Site<-"Adult and egg biomass (kg)"

Salmo_salar_norm_catch_compare_biomass_plot$normalised_reads <- Salmo_salar_norm_catch_compare_biomass_plot$interp_Adult_eggsBM
desired_order2 <- c("Whole river biomass (kg)","Upper river salmon eDNA (normalised reads)", "Conwy falls salmon eDNA (normalised reads)", "Lower river salmon eDNA (normalised reads)")
desired_order3 <- c("E07","E08", "E09", "E10","E11","E12","E13","E14",
                    "Juvenile biomass","Adult and egg biomass (kg)")

Salmo_salar_norm_catch_compare_biomass_plot<-rbind(Salmo_salar_norm_catch_compare_biomass_plot,Salmo_salar_norm_catch_compare)
Salmo_salar_norm_catch_compare_biomass_plot$river_section <- factor(Salmo_salar_norm_catch_compare_biomass_plot$river_section, levels = desired_order2)
Salmo_salar_norm_catch_compare_biomass_plot$Site <- factor(Salmo_salar_norm_catch_compare_biomass_plot$Site , levels = desired_order3)

juvenile_data <- Salmo_salar_norm_catch_compare_biomass_plot
juvenile_data$normalised_reads <- Salmo_salar_norm_catch_compare_biomass_plot$interp_Juvenile
juvenile_data$river_section <- "Whole river biomass (kg)"
juvenile_data$Site <- "Juvenile biomass (kg)"

combined_data <- Salmo_salar_norm_catch_compare_biomass_plot
combined_data$normalised_reads <- Salmo_salar_norm_catch_compare_biomass_plot$Interp_combined_adult_juvenle_BM
combined_data$river_section <- "Whole river biomass (kg)"
combined_data$Site <- "Adult, egg and juvenile biomass (kg)"

Salmo_salar_norm_catch_compare_biomass_plot <- rbind(
  Salmo_salar_norm_catch_compare_biomass_plot,
  juvenile_data,combined_data)

write.csv(Salmo_salar_norm_catch_compare,"Salmo_salar_norm_catch_compare.csv", row.names=FALSE)

biomass_eDNA_lineplot<-ggplot(Salmo_salar_norm_catch_compare_biomass_plot, aes(x = Date_sampled, 
                              y = normalised_reads, color = Site)) + 
  geom_line(size = 1.1) + 
  scale_color_manual(values = c("#999999", "navy", "#0072B2", "#009E73", "#CC79A7", "plum4", 
                                "purple", "magenta", "#D55E00","#E69F00","red2")) +
  scale_x_date(labels = function(d) format(d, "%b %Y"), breaks = "1 month") +
  scale_y_continuous(name = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "", y = "") +
  facet_wrap(~ river_section, ncol = 1, scales = "free_y")+
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    text = element_text(size = 16),         
    axis.title = element_text(size = 18),     
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16))

ggsave("plot_biomass_eDNA_lineplot.png", plot = biomass_eDNA_lineplot, width = 9, height = 13, dpi = 300)


######GAMS######

gam_salmon_rod_full <- gam(normalised_reads ~ 
                             s(temp_monthly_average,k=6) + 
                             s(Days,k=6)+
                             s(Days, by = Site,k=6)+
                             s(daily_flow_by_catchment,k=6) + 
                             s(daily_flow_by_catchment, by = Site,k=6) + 
                             s(pH,k=6)+
                             s(Dist_from_lake,k=6)+  
                             s(interp_Adult_eggsBM, by = Site,k=6)+
                             s(interp_Juvenile, by = Site,k=6)+
                             s(Interp_combined_adult_juvenle_BM, by = Site,k=6),
                           data = Salmo_salar_norm_catch_compare)
summary(gam_salmon_rod_full)

a <- getViz(gam_salmon_rod_full )
check(a)

plot(gam_salmon_rod_full, rug=TRUE, residuals = TRUE, pch =1,shade=TRUE,
     shade.col="lightblue", ylim = c(-0.5, 0.5))#temp/flow/pH/dist_from_lake

plot(gam_salmon_rod_full, rug=TRUE, residuals = TRUE, pch =1,shade=TRUE,
     shade.col="lightblue")

gam_salmon_rod <- gam(normalised_reads ~ 
                        #s(temp_monthly_average,k=6) + 
                        s(Days,k=6)+
                        s(Days, by = Site,k=6)+
                        s(daily_flow_by_catchment,k=6) + 
                        s(daily_flow_by_catchment, by = Site,k=6) + 
                        s(pH,k=6)+
                        s(Dist_from_lake,k=6)+  
                        s(interp_Adult_eggsBM, by = Site,k=6)+
                        s(interp_Juvenile, by = Site,k=6)+
                        s(Interp_combined_adult_juvenle_BM, by = Site,k=6), 
                      data = Salmo_salar_norm_catch_compare)

gam_summary<-summary(gam_salmon_rod)
gam_summary

b <- getViz(gam_salmon_rod)
check(b)

plot(gam_salmon_rod, rug=TRUE, residuals = TRUE, pch =1,shade=TRUE,
     shade.col="lightblue", ylim = c(-0.5, 0.5))#temp/flow/pH/dist_from_lake

plot(gam_salmon_rod, rug=TRUE, residuals = TRUE, pch =1,shade=TRUE,shade.col="lightblue")

data_subset <- Salmo_salar_norm_catch_compare[, c("pH", "Days", "daily_flow_by_catchment")]
cor_matrix <- cor(data_subset, use = "complete.obs")
corrplot(cor_matrix, method = "color", 
         type = "upper",
         order = "hclust", 
         addCoef.col = "black", 
         tl.col = "black",
         tl.srt = 45, 
         diag = FALSE)

average_reads_date <- Salmo_salar_norm_catch_compare %>%
  group_by(Date_sampled) %>%
  summarise(avg_normalised_reads = mean(normalised_reads, na.rm = TRUE))

ggplot() +
  geom_point(data=average_reads_date,aes(x = Date_sampled,y = avg_normalised_reads*40),colour="red") +
  geom_point(data=Salmo_salar_norm_catch_compare, aes(x = Date_sampled,y = temp_monthly_average)) +
  labs(x = "Date", color = "Legend") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(data=Salmo_salar_norm_catch_compare,aes(x = Days,y = pH, colour=Site)) +
  geom_point() +
  geom_smooth()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggplot(data=Salmo_salar_norm_catch_compare,aes(x = temp_monthly_average,y = normalised_reads)) +
  geom_point() +
  geom_smooth()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

data_subset <- Salmo_salar_norm_catch_compare[, c("pH", "Days", "temp_monthly_average", "daily_flow_by_catchment")]

# Compute the correlation matrix
cor_matrix <- cor(data_subset, use = "complete.obs")

# Plot the correlation matrix
corrplot(cor_matrix, method = "color", 
         type = "upper", 
         order = "hclust", 
         addCoef.col = "black", 
         tl.col = "black", 
         tl.srt = 45, 
         diag = FALSE) 

x_labels <- c("Days", "Days","Days","Days","Days","Days","Days","Days","Days",
              "Flow", "Flow","Flow","Flow","Flow","Flow","Flow","Flow","Flow",
              "pH", "Distance from Lake (m)",#"Adult & egg biomass (kg)",
              "Adult & egg biomass (kg)","Adult & egg biomass (kg)",
              "Adult & egg biomass (kg)" ,"Adult & egg biomass (kg)" ,
              "Adult & egg biomass (kg)" ,"Adult & egg biomass (kg)",
              "Adult & egg biomass (kg)","Adult & egg biomass (kg)",#"Juvenile biomass (kg)",
               "Juvenile biomass (kg)", "Juvenile biomass (kg)", "Juvenile biomass (kg)", 
              "Juvenile biomass (kg)","Juvenile biomass (kg)","Juvenile biomass (kg)",
               "Juvenile biomass (kg)", "Juvenile biomass (kg)",#"Adult, egg, juv. biomass (kg)",
              "Adult, egg, juv. biomass (kg)","Adult, egg, juv. biomass (kg)",
              "Adult, egg, juv. biomass (kg)","Adult, egg, juv. biomass (kg)",
              "Adult, egg, juv. biomass (kg)","Adult, egg, juv. biomass (kg)",
              "Adult, egg, juv. biomass (kg)","Adult, egg, juv. biomass (kg)")

plots <- lapply(1:length(gam_salmon_rod$smooth), function(i) {
  p <- draw(gam_salmon_rod, select = i) 
  p <- p + 
    labs(y = "", x = x_labels[i], title = NULL) + 
    theme(legend.position = "none",
          plot.caption = element_blank(),
          plot.title = element_blank())
  return(p)
})
blank_plot <- ggplot() + 
  theme_void()  
plots <- plots[21:length(plots)]
grid.arrange(grobs = plots, nrow = 8)

nrow <- 8
ncol <- ceiling(length(plots) / nrow)  
plots <- plots[order(rep(1:nrow, length.out = length(plots)))]

p_GAM_eDNA_biomass<-grid.arrange(grobs = plots, nrow = nrow)
ggsave("p_GAM_eDNA_biomass.png", p_GAM_eDNA_biomass, width = 8, height = 12, dpi = 300)

#LINE PLOT WITH EDNA AND BIOMASS BY SITE

F_values <- gam_summary$s.table[, "F"]
F_values_df <- data.frame(Term = rownames(gam_summary$s.table), F_value = F_values)

F_values_df <- F_values_df %>%
  filter(!grepl("Days|pH|Dist_from_lake", Term))
F_values_df$Term <- gsub("s\\(|\\)", "", F_values_df$Term)
F_values_df <- F_values_df %>%
  separate(Term, into = c("Measure", "Site"), sep = ":", fill = "right")
F_values_df$Site <- gsub("Site", "", F_values_df$Site)
F_values_df$Measure <- gsub("interp_Adult_eggsBM", "F_value_Adult", F_values_df$Measure)
F_values_df$Measure <- gsub("interp_Juvenile", "F_value_Juv", F_values_df$Measure)
F_values_df$Measure<- gsub("Interp_combined_adult_juvenle_BM", "F_value_Comb", F_values_df$Measure)
F_values_df <- F_values_df %>%
  pivot_wider(names_from = Measure, values_from = F_value)

Salmo_salar_norm_catch_compare <- merge(Salmo_salar_norm_catch_compare, F_values_df, by = "Site", all.x = TRUE)

biomass_eDNA_sites_lineplot<-
  ggplot(Salmo_salar_norm_catch_compare, aes(x = Date_sampled)) + 
  geom_line(aes(y = interp_Adult_eggsBM/6000, color = "Adult and egg biomass", alpha = F_value_Adult), size = 1.1) + 
  geom_line(aes(y = interp_Juvenile/6000, color = "Juvenile biomass", alpha = F_value_Juv), size = 1.1) + 
  geom_line(aes(y = Interp_combined_adult_juvenle_BM/6000, color = "Adult, eggs and juvenile biomass", alpha = F_value_Comb), size = 1.1) +
  geom_line(aes(y = normalised_reads, color = "Normalised reads"), size = 1.1) +
  scale_color_manual(values = c(
    "Normalised reads" = "#56B4E9",
    "Adult and egg biomass" = "#D55E00",
    "Juvenile biomass" = "#E69F00",
    "Adult, eggs and juvenile biomass" = "red2"
  )) +
  scale_x_date(labels = function(d) format(d, "%b %Y"), 
               breaks = "1 month") +
  scale_y_continuous(name = "", 
                     sec.axis = sec_axis(~ . * 6000, name = "Biomass (kg)")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_wrap(~Site, ncol = 2, nrow = 4) +
  labs(x = "", y = "", color = "Measure",alpha="GAM F value") +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    text = element_text(size = 16),         
    axis.title = element_text(size = 18),     
    axis.text = element_text(size = 14),   
    strip.text = element_text(size = 16)
  )
  
ggsave("plot_biomass_eDNA_sites.png", plot = biomass_eDNA_sites_lineplot, width = 13, height = 13, dpi = 300)

#WORKING OUT R SQUARED FOR SITES
r_squared_values <- numeric()
p_values <- numeric()
sample_sites <- unique(Salmo_salar_norm_catch_compare$Site)

for (site in sample_sites) {
  subset_data <- subset(Salmo_salar_norm_catch_compare, Site == site)
  lm_model <- lm(normalised_reads ~ Interp_combined_adult_juvenle_BM,
                 data = subset_data)
  r_squared <- summary(lm_model)$r.squared
  
  p_value <- summary(lm_model)$coefficients["Interp_combined_adult_juvenle_BM", "Pr(>|t|)"]
  cat("Site:", site, "R-squared:", r_squared, "p-value:", p_value, "\n")
  r_squared_values <- c(r_squared_values, r_squared)
  p_values <- c(p_values, p_value)
}

#heatmap
heatmap_data <- Salmo_salar_norm_catch_compare[, c("Site", "Date_sampled", 
          "normalised_reads","interp_Adult_eggsBM","interp_Juvenile","Interp_combined_adult_juvenle_BM")]

heatmap_data <- heatmap_data %>%
  pivot_longer(
    cols = normalised_reads:ncol(heatmap_data),  
    names_to = "variable",                        
    values_to = "value")

heatmap_data <- heatmap_data %>%
  mutate(Site = as.character(Site))%>%
  mutate(Site = ifelse(variable == "interp_Adult_eggsBM", "Adult and egg biomass", Site))%>%
  mutate(Site = ifelse(variable == "interp_Juvenile", "Juvenile biomass", Site))%>%
  mutate(Site = ifelse(variable == "Interp_combined_adult_juvenle_BM", "Adult, egg, juvenile biomass", Site))

heatmap_data<-unique(heatmap_data)
heatmap_data$Site <- gsub("interp_Adult_eggsBM", "Adult and egg biomass", heatmap_data$Site)

heatmap_data <- heatmap_data %>%
  group_by(Site) %>%
  mutate(value = value / sum(value)) %>%
  ungroup()

heatmap_data <- heatmap_data %>%
  complete(Site, Date_sampled = seq(min(Date_sampled), max(Date_sampled), by = "day"))
heatmap_data<- heatmap_data %>%
  group_by(Site) %>%
  mutate(value = na.approx(value, Date_sampled, na.rm = FALSE))

heatmap_data$Measurement <-heatmap_data$Site
heatmap_data <- heatmap_data %>%
  mutate(Measurement = if_else(startsWith(Measurement, "E"), "eDNA", Measurement))

desired_order4<-c("Adult and egg biomass","Juvenile biomass","Adult, egg, juvenile biomass",
                  "E07","E08", "E09", "E10","E11","E12","E13","E14")
heatmap_data$Site <- factor(heatmap_data $Site , levels = desired_order4)
desired_order5<-c("Adult and egg biomass","Juvenile biomass","Adult, egg, juvenile biomass",
                  "eDNA")
heatmap_data$Measurement <- factor(heatmap_data$Measurement, levels = desired_order5)

plot_abundances_heatmap<-ggplot(heatmap_data, aes(x = Site, y = Date_sampled, fill = value)) +
  geom_tile(width = 1.1, height = 1.1)+  
  labs(x = "Site", y = "Date Sampled") +  
  theme_minimal() +  
  scale_fill_viridis(discrete=FALSE) +
  scale_y_date(date_breaks = "1 month", date_labels = "%b %Y") +
  ylab(NULL) +
  xlab(NULL)+labs(fill = "Column\nnormalised value")+
  theme(plot.background = element_rect(fill = "white", colour = NA),
    text = element_text(size = 16),         
    axis.title = element_text(size = 18),     
    axis.text = element_text(size = 14),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 16),
    axis.ticks.y = element_line())

ggsave("plot_abundances_heatmap.png", plot = plot_abundances_heatmap, width = 10, height = 12, dpi = 300)
