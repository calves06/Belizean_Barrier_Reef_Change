## Belize Benthic Change paper
## R Code for public access
library(readr)
library(dplyr)
library(tidyr)
library(scales)
library(ggsci)
library(ggplot2)
library(blme)
library(sjPlot)


### DATA CLEANING ####
## Read data
  MasterDF <- readr::read_csv("Data/Processed/Long.Master.Species.Groups.csv")
  ProtectionDF <- readr::read_csv("Data/Site/Belize_site_coord_protection.csv")
  
## Check that all benthic images add to 100
 MasterDF %>%
   group_by(Site, 
            Year,
            Transect,
            Image.Code) %>% 
   summarise(photocover= sum(Cover)) %>% 
   filter(photocover != 100)

## Add protection level information to master df
 MasterDF_protection <- dplyr::left_join(MasterDF,
                                         ProtectionDF,
                                         by = "Site")
 
 
## Remove sites with only one survey
  master_df <- MasterDF_protection %>% 
    filter(!(Site %in% c("Chapel", 
                        "Glovers_Control", 
                        "Mitchel_Rock", 
                        "Nicholas_Control", 
                        "Pompian_Control", 
                        "Tobacco_Caye"))
           )
  
## Check coral names 
  unique(master_df$ID)
  
### EXPLORATORY PLOTS ####
  
#recode coral species names to standarize names over the years
  master_df$newID <-
    dplyr::recode_factor(master_df$ID,
                         "Acropora.cervicornis" = "Acropora spp.",
                         "Acropora.palmata" =  "Acropora spp.",
                         "Acropora.species" = "Acropora spp.",
                         "Acropora.prolifera" = "Acropora spp.",
                         "Agaricia.agaricites.complex" = "Agaricia agaricites",
                         "Agaricia.grahamae" = "Agaricia spp.",
                         "Agaricia.humilis" = "Agaricia spp.",
                         "Agaricia.lamarcki" = "Agaricia spp.",
                         "Agaricia.fragilis" = "Agaricia spp.",
                         "Agaricia.tenuifolia" = "Agaricia tenuifolia",
                         "Diploria.labyrinthiformis" = "D. labyrinthiformis",
                         "Diploria.strigosa" = "Pseudodiploria spp.",
                         "Diploria.clivosa" = "Pseudodiploria spp.",
                         "Pseudodiploria.clivosa" = "Pseudodiploria spp.",
                         "Pseudodiploria.strigosa" = "Pseudodiploria spp.",
                         "Orbicella.annularis" = "Orbicella spp.",
                         "Orbicella.faveolata" = "Orbicella spp.",
                         "Orbicella.franksi" = "Orbicella spp.",
                         "Porites.astreoides" = "P. astreoides",
                         "Porites.divaricata" = "Porites spp.",
                         "Porites.furcata" = "Porites spp.",
                         "Porites.porites" = "Porites spp.",
                         "Porites.colonensis" = "Porites spp.",
                         "Siderastrea.species" =  "Siderastrea spp.",
                         "Siderastrea.radians" = "Siderastrea spp.",
                         "Siderastrea.siderea" = "Siderastrea spp.",
                         "Millepora" = "Millepora spp.", 
                         "Millepora.complenata" = "Millepora spp.",
                         "Millepora.alcicornis" = "Millepora spp.",
                         "Mycetophyllia" = "Mycetophyllia spp.",
                         "Mycetophyllia.aliciae" = "Mycetophyllia spp.",
                         "Mycetophyllia.danaana" = "Mycetophyllia spp.",
                         "Mycetophyllia.ferox" = "Mycetophyllia spp.",
                         "Mycetophyllia.lamarckiana" = "Mycetophyllia spp.",
                         "Isophyllia.sinuosa" = "Isophyllia spp.",
                         "Isophyllia.rigida" = "Isophyllia spp.",
                         "Scolymia" = "Scolymia spp.",
                         "Scolymia.cubensis" = "Scolymia spp.",
                         "Madracis.auretenra" = "Madracis spp.",
                         "Madracis.decactis" = "Madracis spp.",
                         "Madracis.mirabilis" = "Madracis spp.",
                         "Colpophyllia.natans" = "Colpophyllia natans",
                         "Dendrogyra.cylindrus" = "Dendrogyra cylindrus",
                         "Dichocoenia.stokesi" = "Dichocoenia stokesi",
                         "Eusmilia.fastigiata" = "Eusmilia fastigiata",
                         "Dichocoenia.stokesii" =  "Dichocoenia stokesii",
                         "Favia.fragum" = "Favia fragum",
                         "Leptoseris.cucullata" = "Leptoseris cucullata",
                         "Manicina.areolata" = "Manicina areolata",
                         "Meandrina.meandrites" = "Meandrina meandrites",
                         "Montastrea.cavernosa" = "Montastrea cavernosa",
                         "Mussa.angulosa" = "Mussa angulosa",
                         "Solenastrea.bournoni" = "Solenastrea bournoni",
                         "Stephanocoenia.intersepta" = "Stephanocoenia intersepta"
                         )
  ## change ID to group coral species with less than 0.5 % cover 
  master_df$newID_grouped <-
    dplyr::recode_factor(master_df$ID,
                         "Acropora.cervicornis" = "Acropora spp.",
                         "Acropora.palmata" =  "Acropora spp.",
                         "Acropora.species" = "Acropora spp.",
                         "Acropora.prolifera" = "Acropora spp.",
                         "Agaricia.agaricites.complex" = "Agaricia agaricites",
                         "Agaricia.grahamae" = "Other coral species",
                         "Agaricia.humilis" = "Other coral species",
                         "Agaricia.lamarcki" = "Other coral species",
                         "Agaricia.fragilis" = "Other coral species",
                         "Agaricia.tenuifolia" = "Agaricia tenuifolia",
                         "Diploria.labyrinthiformis" = "Diploria/Pseudodiploria spp.",
                         "Diploria.strigosa" = "Diploria/Pseudodiploria spp.",
                         "Diploria.clivosa" = "Diploria/Pseudodiploria spp.",
                         "Pseudodiploria.clivosa" = "Diploria/Pseudodiploria spp.",
                         "Pseudodiploria.strigosa" = "Diploria/Pseudodiploria spp.",
                         "Orbicella.annularis" = "Orbicella spp.",
                         "Orbicella.faveolata" = "Orbicella spp.",
                         "Orbicella.franksi" = "Orbicella spp.",
                         "Porites.astreoides" = "Porites astreoides",
                         "Porites.divaricata" = "Porites spp.",
                         "Porites.furcata" = "Porites spp.",
                         "Porites.porites" = "Porites spp.",
                         "Porites.colonensis" = "Porites spp.",
                         "Siderastrea.species" =  "Siderastrea spp.",
                         "Siderastrea.radians" = "Siderastrea spp.",
                         "Siderastrea.siderea" = "Siderastrea spp.",
                         "Millepora" = "Other coral species",
                         "Millepora.complenata" = "Other coral species",
                         "Millepora.alcicornis" = "Other coral species",
                         "Mycetophyllia" = "Other coral species",
                         "Mycetophyllia.aliciae" = "Other coral species",
                         "Mycetophyllia.danaana" = "Other coral species",
                         "Mycetophyllia.ferox" = "Other coral species",
                         "Mycetophyllia.lamarckiana" = "Other coral species",
                         "Isophyllia.sinuosa" = "Other coral species",
                         "Isophyllia.rigida" = "Other coral species",
                         "Scolymia" = "Other coral species",
                         "Scolymia.cubensis" = "Other coral species",
                         "Madracis.auretenra" = "Other coral species",
                         "Madracis.decactis" = "Other coral species",
                         "Madracis.mirabilis" = "Other coral species",
                         "Colpophyllia.natans" = "Colpophyllia natans",
                         "Dendrogyra.cylindrus" = "Dendrogyra cylindrus",
                         "Dichocoenia.stokesi" = "Other coral species",
                         "Eusmilia.fastigiata" = "Other coral species",
                         "Dichocoenia.stokesii" =  "Other coral species",
                         "Favia.fragum" = "Other coral species",
                         "Leptoseris.cucullata" = "Other coral species",
                         "Manicina.areolata" = "Other coral species",
                         "Meandrina.meandrites" = "Other coral species",
                         "Montastrea.cavernosa" = "Montastrea cavernosa",
                         "Mussa.angulosa" = "Other coral species",
                         "Solenastrea.bournoni" = "Other coral species",
                         "Stephanocoenia.intersepta" = "Other coral species",
                         "Scleractinia" = "Other coral species"
    )
  
## Check new coral names again
  unique(master_df$newID_grouped)
  
## Calculate cover for coral general by transect
 coral.genera.df <- 
    master_df %>%
      filter(Specific.Type %in% "Hard.coral") %>%
        group_by(Year, Site, Fishing_level, Transect, Image.Code, newID_grouped) %>%
          summarise(total_cover = sum(Cover, na.rm = TRUE)) %>%
             group_by(Year, Site, Transect, Fishing_level, newID_grouped) %>%
                summarise(cover = mean(total_cover, na.rm = TRUE))
        
#mean per site
  mean.coral.genera.df <- coral.genera.df %>% 
      group_by(Site, Year, Fishing_level, newID_grouped) %>% 
        summarise(MeanCover = mean(cover, na.rm=TRUE))

## order names alphabetically
  mean.coral.genera.df$newID_grouped <-
    factor(mean.coral.genera.df$newID_grouped,
           levels = c("Acropora spp.",
                      "Agaricia agaricites",
                      "Agaricia tenuifolia",
                      "Colpophyllia natans",
                      "Dendrogyra cylindrus",
                      "Diploria/Pseudodiploria spp.",
                      "Montastrea cavernosa",
                      "Orbicella spp.",
                      "Porites astreoides",
                      "Porites spp.",
                      "Siderastrea spp.",
                      "Other coral species")
    )
  
## Y axis break function
  count <- 0
  breaks_fun <- function(y) {
    count <<- count + 1L
    switch(
      count,
      c(1, 3, 5, 7, 9),
      c(45, 55),
      c(0, 50, 100),
      seq(0, 8, 0.2)
    )
  }
  

##plot
  coral.genera.plot_20 <- 
    ggplot(data = mean.coral.genera.df,
           aes(x = Year, 
               y = MeanCover))+
    facet_wrap(~ factor(newID_grouped), 
               scales= "free_y")+
    geom_smooth(aes(col = Fishing_level),
                method = "loess",
                formula = 'y ~ x',
                span = 1,
                alpha = 0.05)+
     geom_point(aes(col = Fishing_level),
               position = position_dodge(width = 0.75),alpha = 0.2) +
    scale_y_continuous(limits = c(0,max(mean.coral.genera.df$MeanCover)), #change scale to overlay plots
                       oob = squish) + #to avoid CI with negative values
    scale_color_nejm()+
    scale_fill_nejm()+
     labs(x=NULL, 
         y="Mean cover (%)")+
    theme_bw() +
    theme(legend.direction="horizontal",
          legend.position="bottom",
          axis.text.x=element_text(angle=0, hjust=0.5),
          panel.grid = element_blank(),
          strip.text = element_text(face = "italic"),
          strip.background = element_rect(fill = "grey95")) +
    labs(col ='Fishing level') 
  coral.genera.plot_20
 
## Save plot
 ggsave("Figures/coral_genera_plot_scale20.png",
        device =  "png",
        width = 7.5,
        height = 7, 
        units ="in", 
        dpi = 600)
        
 
 ## Recode and group specific Type to standarize over years 
 unique(master_df$Specific.Type)
 
  master_df$newType <- recode_factor(
    master_df$Specific.Type,
    "Hard.coral" = "Live coral",
    "Fleshy.macroalgae" = "Macroalgae",
    "Halimeda" = "Macroalgae",
    "Soft.coral" = "Gorgonian",
    "Coraline" = "Macroalgae",
    "Calcareous.algae" = "Macroalgae",
    "Anemone" = "Other invertebrates",
    "Annelid" = "Other invertebrates",
    "Ascidian" = "Other invertebrates",
    "Barnacle" = "Other invertebrates",
    "Bivalve" = "Other invertebrates",
    "Bryozoa" = "Other invertebrates",
    "Corallimorph" = "Other invertebrates",
    "Hydroid" = "Other invertebrates",
    "Nudibranch" = "Other invertebrates",
    "Tunicate" =  "Other invertebrates",
    "Zoanthid" = "Other invertebrates",
    "Sea.cucumber" = "Other invertebrates",
    "CCA" = "CTB", # Combined, CCA, Turf, and Bare substrate
    "Substrate" = "CTB",
    "Dead" = "CTB",
    "Rubble" = "Substrate",
    "Sand.sediment" = "Substrate",
    "Bacterial.mat" = "Cyanobacteria",
    "Fish" = "Other",
    "Unknown" = "Other",
    "Water" = "Other",
    "Equipment" = "Other",
    "N.c" = "Other"
    )
  
  ## Calculate cover for benthic groups  by transect
  benthic.groups.df <- 
    master_df %>%
     group_by(Year, Site, Fishing_level, Transect, Image.Code, newType) %>%
      summarise(total_cover = sum(Cover, na.rm = TRUE)) %>%
        group_by(Year, Site, Transect, Fishing_level, newType) %>%
          summarise(cover = mean(total_cover, na.rm = TRUE))
  
  #mean per site and filter for key categories 
  mean.benthic.groups.df <- benthic.groups.df %>% 
     group_by(Site, Year, Fishing_level, newType) %>% 
        summarise(MeanCover = mean(cover, na.rm=TRUE)) %>%
          filter (newType %in% c("Live coral",
                                 "Macroalgae",
                                 "CTB",
                                 "Gorgonian",
                                 "Sponge"))
  #Sort benthic groups
  mean.benthic.groups.df$newType <- 
    factor (mean.benthic.groups.df$newType,
            levels = c("Live coral", 
                       "Macroalgae",
                       "CTB",
                       "Gorgonian",
                       "Sponge"))
  
 #mean hard coral cover in 1997
  mean.benthic.groups.df %>%
    filter (newType == "Live coral" &
              Year == 1997) %>%
    group_by(Year) %>%
     summarise(mean(MeanCover),
               sd(MeanCover))
  
  #mean hard coral cover in 2016
  mean.benthic.groups.df %>%
    filter (newType == "Live coral" &
              Year == 2016) %>%
    group_by(Year) %>%
     summarise(mean(MeanCover),
              sd(MeanCover))
  
##plot
  benthic.groups.plot <- 
    ggplot(data = mean.benthic.groups.df,
           aes(x = Year, 
               y = MeanCover))+
    facet_wrap(~ factor(newType), 
               scales= "free_y")+
    geom_smooth(aes(col = Fishing_level),
              method = "loess",
              formula = 'y ~ x',
              span = 1,
              alpha = 0.05)+
    geom_point(aes(col = Fishing_level),
               position = position_dodge(width = 0.75),
               alpha = 0.2) +
    scale_y_continuous(limits = c(0, 
                                  max(mean.benthic.groups.df$MeanCover)),
                       oob = squish) + #to avoid CI with negative values
    scale_color_nejm()+
    scale_fill_nejm()+
    labs(x=NULL, 
         y="Mean Cover (%)")+
    theme_bw()+
    theme(legend.direction="vertical", 
          legend.position= c(0.8,0.2),
          legend.title = element_text (size = 10),
          axis.text.x=element_text(angle=0, 
                                   hjust=0.5),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95"))+
    labs(col ='Fishing level') 
  
  benthic.groups.plot
  
  ggsave("Figures/Fig.2_benthic_groups_plot_trends.png",
         device =  "png",
         width = 7,
         height = 5, 
         units ="in", 
         dpi = 600)
  

## Calculate cover for coral general by transect
  
  master_df$new_algae_id <- recode_factor(
    master_df$ID,
    "Amphiroa.species" = "branching calcareous",
    "Branching.coraline" = "branching calcareous",
    "Branching.coraline.antillarum" = "branching calcareous",
    "Calcareous" = "branching calcareous",
    "Dictyota" = "fleshy macroalgae",
    "Encrusting.coraline" = "CT",
    "Erect.rhodophyta" = "corticate",
    "Fleshy_Macroalgae" = "fleshy macroalgae",
    "Galaxaura.species" = "corticate",
    "Halimeda" = "branching calcareous",
    "Lobophora" = "corticate", 
    "Macroalgae" = "fleshy macroalgae", 
    "Padina.species" = "corticate",
    "Turf.algae" = "CT",
    "Ventricaria.species" = "corticate",
    "Wrangelia" = "fleshy macroalgae",
    "Crustose.turf.bare" = "CT",
    "CCA" = "CT"
  )
  
  algae.genera.df <- 
    master_df %>%
    filter(Specific.Type %in% c("Macroalgae", 
                                "Coraline", 
                                "Calcareous.algae", 
                                "Fleshy.macroalgae",
                                "Halimeda",
                                "CT"
                                )) %>%
    group_by(Year, Site, Fishing_level, Transect, Image.Code, new_algae_id) %>%
      summarise(total_cover = sum(Cover, na.rm = TRUE)) %>%
        group_by(Year, Site, Transect, Fishing_level, new_algae_id) %>%
          summarise(cover = mean(total_cover, na.rm = TRUE))
  
  unique(algae.genera.df$new_algae_id)
  
  #mean per site
  mean.algae.genera.df <- algae.genera.df %>% 
    group_by(Site, Year, Fishing_level, new_algae_id) %>% 
      summarise(MeanCover = mean(cover, na.rm=TRUE))
                                     
  ##plot
  algae.genera.plot <- 
    ggplot(data =  mean.algae.genera.df,
           aes(x = Year, 
               y = MeanCover))+
    facet_wrap(~ factor(new_algae_id), 
               scales= "free_y")+
    geom_smooth(aes(fill = Fishing_level,
                   col = Fishing_level),
               method = "glm",
               formula = 'y ~ x',
               span = 1,
               alpha = 0.05)+
    geom_point(aes(fill = Fishing_level,
                   col = Fishing_level),
               position = position_dodge(width = 0.75),
               alpha = 0.2) +
    scale_y_continuous(#limits = c(0, max(mean.benthic.groups.df$MeanCover)),
      oob = squish) + #to avoid CI with negative values
    scale_color_nejm()+
    scale_fill_nejm()+
    labs(x=NULL, 
         y="Mean cover of macroalgae (%)")+
    theme_bw() +
    theme(legend.direction="horizontal", 
          legend.position= "bottom", #(0.8, 0.1),
          axis.text.x=element_text(angle=0, hjust=0.5),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey95")) +
    guides(fill=guide_legend(nrow=2, byrow=TRUE))
  
  algae.genera.plot
  

  ### BLMER MODELS #####
  
  # Model Development 
  ## 1. Read in Master data with covariates frame and join with benthic 
  ### It contains Mean Cover, SE, Temperature, Protection & HII data     
  Cover.Prot.Temp.HiiDF <- readr::read_csv("Data/Processed/Cover.Prot.Temp.HiiDF.csv")
  
  ## summarize covariates by site, year, and HII - 50km
  Cover.Prot.Temp.HiiDF_summary <-
    Cover.Prot.Temp.HiiDF[, c("Year", "Site", 
                               "SSTA_Freq",
                               "SSTA_Freq_hist",
                               "TSA_Freq",
                               "TSA_Freq_hist",
                               "TSA_Freq_btw_surveys",
                               "Longitude",
                               "Latitude",
                               "Buffer",
                               "HII.Value")] %>%
        filter (Buffer == 50) %>%
          group_by (Site, Year) %>%
            summarise(SSTA_Freq = mean (SSTA_Freq),
                      SSTA_Freq_hist = mean (SSTA_Freq_hist),
                      TSA_Freq = mean(TSA_Freq),
                      TSA_Freq_hist = mean(TSA_Freq_hist),
                      TSA_Freq_btw_surveys = mean(TSA_Freq_btw_surveys),
                      Longitude = mean(Longitude),
                      Latitude = mean(Latitude),
                      HII_50km = mean(HII.Value))
  
  dim(Cover.Prot.Temp.HiiDF_summary)
  names(Cover.Prot.Temp.HiiDF_summary)
  
  
  
  ### Join benthic.groups.df (from above) with covariates 
  ## use Site means
  benthic.groups.df.cov <-
    left_join (mean.benthic.groups.df,
               Cover.Prot.Temp.HiiDF_summary,
                by = c("Site", "Year"))
               
  dim(benthic.groups.df)
  dim(benthic.groups.df.cov)
  names(benthic.groups.df.cov)
  
  #Check sites 
  unique(benthic.groups.df.cov$Site)
  
  #Check benthic groups
  unique(benthic.groups.df.cov$newType)
  
 
  ##  2. Logit Transform MeanCover and Weight by SE    
  LogitMasterDF <- benthic.groups.df.cov %>%
    dplyr::mutate(Logit.Cover = car::logit(MeanCover,
                                           percents = TRUE)) 
  str(LogitMasterDF)
  
   ## 3. Rescale  predictor variables. We want to subtract the mean then divide by two standard deviations - I think we can do this with the arm::rescale() function
  
  # first center year to help w/ model convergence
  LogitMasterDF$centYear <- LogitMasterDF$Year-2005
  head(LogitMasterDF)
  
  # Rescale all predictor variables, including the centered year:
  LogitMasterDF$RS.centYear <- arm::rescale(LogitMasterDF$centYear)
  LogitMasterDF$RS.Hii.50 <- arm::rescale(LogitMasterDF$HII_50km)
  LogitMasterDF$RS.SSTA_Freq <- arm::rescale(LogitMasterDF$SSTA_Freq)
  LogitMasterDF$RS.SSTA_Freq_hist <- arm::rescale(LogitMasterDF$SSTA_Freq_hist)
  LogitMasterDF$RS.TSA_Freq <- arm::rescale(LogitMasterDF$TSA_Freq)
  LogitMasterDF$RS.TSA_Freq_hist <- arm::rescale(LogitMasterDF$TSA_Freq_hist)
  LogitMasterDF$RS.TSA_Freq_btw_surveys <- arm::rescale(LogitMasterDF$TSA_Freq_btw_surveys)
                
  dim(LogitMasterDF)
  str(LogitMasterDF)
  
  #Create factor year
  LogitMasterDF$factorYear <- as.factor(LogitMasterDF$Year)
 
  ## 5. Fit Bayesian Generalized Linear Mixed-Effects Models
  
  # Model Logit.Cover as the response, 
  #filtering for newType 
  
  live.coral.hist.randsite <- 
    blme::blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site), 
          REML=TRUE,
          data=LogitMasterDF %>% 
              filter(newType == "Live coral")
          ) 
  
  summary(live.coral.hist.randsite)
  AIC(live.coral.hist.randsite)
  plot_model(live.coral.hist.randsite)
  

  # Let's make a model w/ Hii at 50, but an interaction b/w Year and Protection, 
  # and adding a random effect for year as a factor:
  live.coral.int.hist.randsite.randyear <- 
    blmer(Logit.Cover~
            RS.centYear*Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site)+
            (1|factorYear), 
          REML =TRUE,
          data=LogitMasterDF %>% 
            filter(newType == "Live coral")) 
  
  summary(live.coral.int.hist.randsite.randyear)
  AIC(live.coral.int.hist.randsite.randyear)
  plot_model(live.coral.int.hist.randsite.randyear)
  car::vif(live.coral.int.hist.randsite.randyear)
  
  # now let's try a model for SSTA_Freq_hist w/o int b/w year and prot.code:
  live.coral.hist.randsite.randyear <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site)+
            (1|factorYear), 
          data=LogitMasterDF %>%
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.hist.randsite.randyear) 
  AIC(live.coral.hist.randsite.randyear)
  plot_model(live.coral.hist.randsite.randyear)
  
  
  
  # Let's make the same model as above, but use SSTA_Freq instead of hist:
  live.coral.int.freq.randsite.randyear <- 
    blmer(Logit.Cover~
            RS.centYear*Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site)+
            (1|factorYear), 
          data=LogitMasterDF %>% 
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.int.freq.randsite.randyear) 
  plot_model(live.coral.int.freq.randsite.randyear)
  
  # Abel suggested we try a model w/o interaction b/w year and prot.code - for SSTA_Freq:
  live.coral.freq.randsite.randyear <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site)+
            (1|factorYear), 
          data=LogitMasterDF %>% 
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.freq.randsite.randyear) 
  plot_model(live.coral.freq.randsite.randyear)
  
  # same model as above, but without year as a random effect:
  live.coral.freq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.freq.randsite) 
  plot_model(live.coral.freq.randsite)
  
  # model with interaction b/w year and prot, same as above otherwise:
  live.coral.int.freq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear*Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site), 
          data=LogitMasterDF %>%
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.int.freq.randsite)
  plot_model(live.coral.int.freq.randsite)
  
  # model with TSA_Freq NOT SSTA_Freq
  live.coral.TSAfreq.randsite <- 
    bglmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Live coral"),
          REML=TRUE) 
  
  summary(live.coral.TSAfreq.randsite)
  AIC(live.coral.TSAfreq.randsite)
  plot_model(live.coral.TSAfreq.randsite)
  
  # model with TSA_Freq_hist
  live.coral.TSAhist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>%
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.TSAhist.randsite) 
  AIC(live.coral.TSAhist.randsite)
  plot_model(live.coral.TSAhist.randsite)
  
  # model with TSA_Freq_btw_surveys
  live.coral.TSAfreq.surv.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_btw_surveys +
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Live coral"), 
          REML=TRUE) 
  
  summary(live.coral.TSAfreq.surv.randsite) 
  AIC(live.coral.TSAfreq.surv.randsite)
  plot_model(live.coral.TSAfreq.surv.randsite)
  
  ### Compare models
  anova(live.coral.hist.randsite, 
      live.coral.int.hist.randsite.randyear, 
      live.coral.hist.randsite.randyear, 
      live.coral.int.freq.randsite.randyear, 
      live.coral.freq.randsite.randyear, 
      live.coral.freq.randsite, 
      live.coral.TSAfreq.randsite, 
      live.coral.TSAhist.randsite, 
      live.coral.TSAfreq.surv.randsite)
  
  ### Hard coral models comparison for paper
  hard_coral_mod_anova <-
    anova(live.coral.TSAfreq.randsite, 
        live.coral.TSAhist.randsite, 
        live.coral.TSAfreq.surv.randsite)
  write.csv(hard_coral_mod_anova, 
            "Data/Processed/hard_coral_mod_anova.csv")
  
  # Model with lowest AIC is live.coral.TSAfreq.randsite,
  
  # Test for VIF (variance inflation factor) 
  
  # use function 'vif' to calculate variance inflation factors for each variable in model
  car::vif(live.coral.hist.randsite) #yes collinearity
  car::vif(live.coral.int.hist.randsite.randyear)  
  car::vif(live.coral.hist.randsite.randyear)
  car::vif(live.coral.int.freq.randsite.randyear)
  car::vif(live.coral.freq.randsite.randyear)
  car::vif(live.coral.freq.randsite)
  car::vif(live.coral.TSAfreq.randsite)
  car::vif(live.coral.TSAhist.randsite) # yes collinearity
  car::vif(live.coral.TSAfreq.surv.randsite) 
  
  
  ##  **5.**  Model other Specific.Types  
  ## Now we make models for the other Specific.Types of interest using the convention from Hard.coral.freq.randsite.randyear
  
  # models for Macroalgae:
  # using SSTA_freq
  Macro.freq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Macroalgae"), 
          REML=TRUE) 
  
  summary(Macro.freq.randsite) 
  plot_model(Macro.freq.randsite)
  car::vif(Macro.freq.randsite)
  
  # using SSTA_Freq_hist
  Macro.hist.randsite <- 
    blmer(Logit.Cover ~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Macroalgae"), 
          REML=TRUE) 
  
  summary(Macro.hist.randsite)
  plot_model(Macro.hist.randsite)
  
  # using TSA_Freq
  Macro.TSAfreq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Macroalgae"),
          REML=TRUE) 
  
  summary(Macro.TSAfreq.randsite)
  plot_model(Macro.TSAfreq.randsite)
  car::vif(Macro.TSAfreq.randsite)
  
  # model, like above, instead with TSA_Freq_hist
  Macro.TSAhist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Macroalgae"), 
          REML=TRUE) 
  
  summary(Macro.TSAhist.randsite)
  plot_model(Macro.TSAhist.randsite)
  
  # use TSA_Freq_btw_surveys
  Macro.TSAfreq.surv.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_btw_surveys+
            (1|Site), 
          data=LogitMasterDF %>%
            filter(newType == "Macroalgae"), 
          REML=TRUE) 
  
  summary(Macro.TSAfreq.surv.randsite)
  plot_model(Macro.TSAfreq.surv.randsite)
  
  ### Compare models
  anova(Macro.freq.randsite,
      Macro.hist.randsite,
      Macro.TSAfreq.randsite,
      Macro.TSAhist.randsite,
      Macro.TSAfreq.surv.randsite)
  
  macro_mod_anova <-
    anova(Macro.TSAfreq.randsite, 
          Macro.TSAhist.randsite, 
          Macro.TSAfreq.surv.randsite)
  write.csv(macro_mod_anova, 
            "Data/Processed/macro_mod_anova.csv")
      

  # model for CTB - crustose, turf, bare
  # using SSTA_Freq
  CTB.freq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "CTB"), 
          REML=TRUE)
  
  summary(CTB.freq.randsite)
  plot_model(CTB.freq.randsite)
  
  # using SSTA_Freq_hist
  CTB.hist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "CTB"), 
          REML=TRUE) 
  
  summary(CTB.hist.randsite)
  plot_model(CTB.hist.randsite)
  
  # using TSA_Freq
  CTB.TSAfreq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq+
            (1|Site), 
          data=LogitMasterDF %>%
            filter(newType == "CTB"), 
          REML=TRUE)
  
  summary(CTB.TSAfreq.randsite)
  plot_model(CTB.TSAfreq.randsite)
  
  # model, like above, instead with TSA_Freq_hist
  CTB.TSAhist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>%
            filter(newType == "CTB"), 
          REML=TRUE)
  
  summary(CTB.TSAhist.randsite)
  plot_model(CTB.TSAhist.randsite)
  
  # use TSA_Freq_btw_surveys
  CTB.TSAfreq.surv.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_btw_surveys+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "CTB"), 
          REML=TRUE)
  
  summary(CTB.TSAfreq.surv.randsite)
  plot_model(CTB.TSAfreq.surv.randsite)
  
  ## AIC
  anova(CTB.freq.randsite,
      CTB.hist.randsite,
      CTB.TSAfreq.randsite,
      CTB.TSAhist.randsite,
      CTB.TSAfreq.surv.randsite)
  
  ctb_mod_anova <-
    anova(CTB.TSAfreq.randsite, 
          CTB.TSAhist.randsite, 
          CTB.TSAfreq.surv.randsite)
  write.csv(ctb_mod_anova, 
            "Data/Processed/ctb_mod_anova.csv")
  
  # model for gorgonian (formerly soft coral)
  # using SSTA_freq
  gorg.freq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Gorgonian"), 
          REML=TRUE)
  
  summary(gorg.freq.randsite)
  plot_model(gorg.freq.randsite)
  
  # using SSTA_Freq_hist
  gorg.hist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Gorgonian"), 
          REML=TRUE) 
  
  summary(gorg.hist.randsite)
  
  # using TSA_Freq
  gorg.TSAfreq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Gorgonian"), 
          REML=TRUE)
  
  summary(gorg.TSAfreq.randsite)
  plot_model(gorg.TSAfreq.randsite)
  
  
  # model, like above, instead with TSA_Freq_hist
  gorg.TSAhist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>%
            filter(newType == "Gorgonian"), 
          REML=TRUE)
  
  summary(gorg.TSAhist.randsite)
  plot_model(gorg.TSAhist.randsite)
  
  # use TSA_Freq_btw_surveys
  gorg.TSAfreq.surv.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_btw_surveys+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Gorgonian"), 
          REML=TRUE)
  
  summary(gorg.TSAfreq.surv.randsite)
  plot_model(gorg.TSAfreq.surv.randsite)
  
  ## compare models
  anova(gorg.freq.randsite,
        gorg.hist.randsite,
        gorg.TSAfreq.randsite,
        gorg.TSAhist.randsite,
        gorg.TSAfreq.surv.randsite)
  
  gorg_mod_anova <-
    anova(gorg.TSAfreq.randsite, 
          gorg.TSAhist.randsite, 
          gorg.TSAfreq.surv.randsite)
  write.csv(gorg_mod_anova, 
            "Data/Processed/gorg_mod_anova.csv")
  
  # model for sponge
  # using SSTA_Freq
  Sponge.freq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Sponge"), 
          REML=TRUE)
  
  summary(Sponge.freq.randsite)
  plot_model(Sponge.freq.randsite)
  
  # using SSTA_Freq_hist
  Sponge.hist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.SSTA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Sponge"), 
          REML=TRUE) 
  
  summary(Sponge.hist.randsite)
  plot_model(Sponge.hist.randsite)
  
  
  # using TSA_Freq
  Sponge.TSAfreq.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Sponge"), 
          REML=TRUE)
  
  summary(Sponge.TSAfreq.randsite)
  plot_model(Sponge.TSAfreq.randsite)
  
  # model, like above, instead with TSA_Freq_hist
  Sponge.TSAhist.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_hist+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Sponge"), 
          REML=TRUE)
  
  summary(Sponge.TSAhist.randsite)
  plot_model(Sponge.TSAhist.randsite)
  
  # use TSA_Freq_btw_surveys
  Sponge.TSAfreq.surv.randsite <- 
    blmer(Logit.Cover~
            RS.centYear+
            Fishing_level+
            RS.Hii.50+
            RS.TSA_Freq_btw_surveys+
            (1|Site), 
          data=LogitMasterDF %>% 
            filter(newType == "Sponge"), 
          REML=TRUE)
  
  summary(Sponge.TSAfreq.surv.randsite)
  plot_model(Sponge.TSAfreq.surv.randsite)
  
  ## compare models
  anova(Sponge.freq.randsite,
        Sponge.hist.randsite,
        Sponge.TSAfreq.randsite,
        Sponge.TSAhist.randsite,
        Sponge.TSAfreq.surv.randsite)
  
  sponge_mod_anova <-
    anova(Sponge.TSAfreq.randsite, 
          Sponge.TSAhist.randsite, 
          Sponge.TSAfreq.surv.randsite)
  write.csv(sponge_mod_anova, 
            "Data/Processed/sponge_mod_anova.csv")
  
  ## 5. Compare AIC values  
  
  ### **5a.** Plot response versus predictor(s)  
   # plot live coral versus temp - TSAFreq.randsite
  live.coral.tsa.plot <- 
      LogitMasterDF %>% 
          filter(newType == "Live coral") %>% 
            ggplot(aes(x = TSA_Freq,
                       y = boot::inv.logit( ## back transform to percentage
                              predict (live.coral.TSAfreq.randsite,
                                      type = "response"))*100))+ 
              geom_point(color = "grey20", 
                         size = 3,
                         alpha = 0.5)+
              geom_smooth(method="glm",
                          color = "darkblue",
                          fill = "blue",
                          alpha = 0.1)+
                  labs(title = "HARD CORAL",
                        x = "TSA Frequency", 
                       y = "Predicted cover (%)")+
                scale_y_continuous(limits=c(0,30))+
                theme_bw()+
                theme(panel.grid = element_blank(),
                      axis.text = element_text(size = 14),
                      axis.title = element_text(size = 14),
                      plot.title = element_text(hjust = 0.5))
   live.coral.tsa.plot
  
  # plot gorgonian versus TSA freq
   
   gorg.coral.temp.plot <- 
     LogitMasterDF %>% 
     filter(newType == "Gorgonian") %>% 
     ggplot(aes(x = TSA_Freq,
                y = boot::inv.logit( ## back transform to percentage
                  predict (gorg.TSAfreq.randsite,
                           type = "response"))*100))+ 
     geom_point(color = "grey20", 
                size = 3,
                alpha = 0.5)+
     geom_smooth(method="glm",
                 color = "darkblue",
                 fill = "blue",
                 alpha = 0.1)+
     labs(title = "GORGONIAN",
          x = "TSA Frequency", 
          y = NULL)+
     scale_y_continuous(limits=c(0,30))+
     theme_bw()+
     theme(panel.grid = element_blank(),
           axis.text = element_text(size = 14),
           axis.title = element_text(size = 14),
           plot.title = element_text(hjust = 0.5))
   gorg.coral.temp.plot
   
  coral_gorg_tsa <-
    cowplot::plot_grid(live.coral.tsa.plot,
                      gorg.coral.temp.plot)
   
   cowplot::save_plot("Figures/Fig.4 coral_gorg_tsa.png",
                      coral_gorg_tsa,
                      base_height = 4,
                      base_width = 6.5)
                        
                      
                      
 
  
  ###  **6.**  Test for Homoscedascicity  
  #resid vs fitted
  plot(live.coral.TSAfreq.randsite)
  plot(Macro.TSAfreq.randsite)
  plot(CTB.TSAfreq.randsite)
  plot(gorg.TSAfreq.randsite)
  plot(Sponge.TSAfreq.randsite)
  
  #q-q plot
  qq_plot <- function (.model) {
    qqnorm(residuals(.model), main = " ") 
    qqline(residuals(.model))
  }
  
  par(mfrow = c(2,3))
  #Coral
  qq_plot(live.coral.TSAfreq.randsite)
    title("Coral \n Normal Q-Q plot")
  #Macroalgae
  qq_plot(Macro.TSAfreq.randsite)
    title("Macroalgae \n Normal Q-Q plot")
  #CTB
  qq_plot(CTB.TSAfreq.randsite)
    title("CTB \n Normal Q-Q plot")
  #Gorgonian
  qq_plot(gorg.TSAfreq.randsite)
    title("Gorgonia \n Normal Q-Q plot")
  #Sponge
  qq_plot(Sponge.TSAfreq.randsite)
    title("Sponge \n Normal Q-Q plot")
  

  ##  **7.**  Generate Confusion Matrix 
  
  ## Generate mean values of all points and the sigma value for data for each model
    
  selectedModel <- live.coral.TSAfreq.randsite
  benthic_group <- "Live coral" #
    
  nullmean <- predict(selectedModel)
    head(nullmean) 
  nullsigma <- summary(selectedModel)$sigma
  
  # generate new data set based on null model; add a normal deviate to each data point
  null.sim <- sapply(nullmean,
                     function(x) rnorm(1, x, nullsigma)) # x = mean, sigma = null sigma
  
  # fit the model to these simulated data. 
  sim.fit <- blmer(null.sim ~
                     RS.centYear +
                     Fishing_level+
                     RS.Hii.50+
                     RS.TSA_Freq+
                     (1|Site), 
                   data=LogitMasterDF %>% 
                     filter(newType == benthic_group),
                   REML = TRUE) 
  plot_model(sim.fit)

  
  ### COEFFICIENT PLOTS ####
  ## **10.** Format Effect Sizes Plot   
  
  # hard coral
  hard.coral.plot <- 
    sjPlot::plot_model(live.coral.TSAfreq.randsite,
                show.values = TRUE, 
                show.p=TRUE,
                type ="est",
                axis.labels = c("TSA Freq", 
                                "HII - 50 Km", 
                                "Fishing vs. \nNo Fishing", 
                                "Year"), #labels are in order from bottom to top!!
                axis.title = " ",
                title="HARD CORAL",
                color = "bw",
                vline.color = "grey") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,
                                    size = 14),
          panel.grid = element_blank())
  
  # macroalgae
  macro.plot <- 
    sjPlot::plot_model(Macro.TSAfreq.randsite,
                 show.values = TRUE, 
                 show.p=TRUE,
                 type ="est", 
                 axis.labels = c(" ", " ", " ", " "), 
                 axis.title = " ",
                 title="MACROALGAE",
                 color = "bw",
                 vline.color = "grey") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,
                                    size = 14),
          panel.grid = element_blank())
  
  # CTB 
  ctb.plot <- 
    sjPlot::plot_model(CTB.TSAfreq.randsite,
                  show.values = TRUE, 
                  show.p=TRUE,
                  type ="est",
                  axis.labels = c(" ", " ", " "," "),
                  axis.title = "Coefficient Estimates",
                  title="CTB",
                  color = "bw",
                  vline.color = "grey")+
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,
                                    size = 14),
          panel.grid = element_blank())
  
  # gorgonian
  gorg.plot <- 
    sjPlot::plot_model(gorg.TSAfreq.randsite,
                  show.values = TRUE, 
                  show.p=TRUE,
                  type ="est",
                  axis.labels = c(" "," ", " "," "), #labels are in order from bottom to top!!
                  axis.title = " ",
                  title="GORGONIAN",
                  color = "bw",
                  vline.color = "grey") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5,
                                    size = 14),
          panel.grid = element_blank())
  
  # sponge
  sponge.plot <- 
    sjPlot::plot_model(Sponge.TSAfreq.randsite,
                   show.values = TRUE, 
                   show.p=TRUE,
                   type = "est",
                   p.threshold = 0.05,
                   axis.labels =c(" ", " ", " "," "), #labels are in order from bottom to top!!
                   axis.title = " ",
                   title="SPONGE",
                   color = "bw",
                   vline.color = "grey") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 12),
          plot.title = element_text(hjust = 0.5, 
                                    size = 14),
          panel.grid = element_blank())
  
  
  # can i use cowplot to combine?
  effects.combo.plot <- 
    cowplot::plot_grid(hard.coral.plot, 
                       macro.plot,
                       ctb.plot, 
                       gorg.plot, 
                       sponge.plot, 
                       ncol=5,
                       rel_widths = c(1.5,1,1,1,1))
  
  # view plots combined
  effects.combo.plot
  
  # save as PNG
  cowplot::save_plot("Figures/Fig.3_effects.combo.plot_v3.png", 
                     effects.combo.plot, 
                     base_width = 10, 
                     base_height = 5,
                     )
 
  ### BLMER MODEL RESULTS ####
  ## **11.** Make Model Output Table    
  
  # unformatted table 
  sjPlot::tab_model(live.coral.TSAfreq.randsite) # look to see what the predictors are here so you can label them properly below!
  sjPlot::tab_model(Macro.TSAfreq.randsite) # look to see what the predictors are here so you can label them properly below!
  sjPlot::tab_model(CTB.TSAfreq.randsite) # look to see what the predictors are here so you can label them properly below!
  sjPlot::tab_model(gorg.TSAfreq.randsite) # look to see what the predictors are here so you can label them properly below!
  sjPlot::tab_model(Sponge.TSAfreq.randsite) # look to see what the predictors are here so you can label them properly below!
  
  # make a table with all the models
  sjPlot::tab_model(live.coral.TSAfreq.randsite,
                    Macro.TSAfreq.randsite,
                    CTB.TSAfreq.randsite,
                    gorg.TSAfreq.randsite,
                    Sponge.TSAfreq.randsite,
                    show.re.var=TRUE,
                    pred.labels = c("(Intercept)", 
                                    "Year", 
                                    "Fishing vs. No Fishing", 
                                    "HII - 50 Km", 
                                    "TSA Freq"),
                    dv.labels = c("Live Coral", 
                                  "Macroalgae",
                                  "CTB", 
                                  "Gorgonian",
                                  "Sponge"))
  
  # Make tables into csv files so you can copy+paste into MS word:
  
  hc.mod.data <- sjPlot::get_model_data(live.coral.TSAfreq.randsite, 
                                        type = "est",
                                        show.intercept = TRUE)
     write.csv(hc.mod.data, 
               "Data/Processed/hc.mod.df.csv")
  
  macro.mod.data <- sjPlot::get_model_data(Macro.TSAfreq.randsite, 
                                           type = "est",
                                           show.intercept = TRUE)
     write.csv(macro.mod.data, 
               "Data/Processed/macro.mod.df.csv")
     
  ctb.mod.data <-  sjPlot::get_model_data(CTB.TSAfreq.randsite, 
                                          type = "est",
                                          show.intercept = TRUE)
     write.csv(ctb.mod.data, 
               "Data/Processed/ctb.mod.df.csv")
     
  gorg.mod.data <- sjPlot::get_model_data(gorg.TSAfreq.randsite, 
                                          type = "est",
                                          show.intercept = TRUE)
      write.csv(gorg.mod.data, 
                "Data/Processed/gorg.mod.df.csv")
  
  sponge.mod.data <- sjPlot::get_model_data(Sponge.TSAfreq.randsite, 
                                            type = "est",
                                            show.intercept = TRUE)
   write.csv(sponge.mod.data, 
                  "Data/Processed/sponge.mod.df.csv")
  
  
  ## **12.** Plot Model Estimates with Data   
  ### **12a.** Save effect size estimates into a data frame  
  # using the effects pkg, use effect function. term = the fixed effect you want to get data on, mod = name of model
  
  ## for hard coral
  hard_coral_effects_year <- 
    effects::effect(term = "RS.centYear",
                    mod = live.coral.TSAfreq.randsite)
  
  #summary(hard_coral_effects_year)
  # save the effects values as a df:
  hc_x_year <- as.data.frame(hard_coral_effects_year)
  
  # recode RS.centYear to our actual years:
  hc_x_year <- hc_x_year %>% 
    dplyr::mutate(Year = dplyr::recode(RS.centYear, 
                                             `-0.60` = "1997",
                                             `-0.30` = "1999",
                                             `0.06` = "2005",
                                             `0.40` = "2009",
                                             `0.80` = "2016"))
  ## for CCA
  macro_effects_year <- 
    effects::effect(term = "RS.centYear", 
                    mod = Macro.TSAfreq.randsite)
  #summary(CCA_effects_year)
  # save the effects values as a df:
  macro_x_year <- as.data.frame(macro_effects_year)
  
  # recode RS.centYear to our actual years:
  macro_x_year <- macro_x_year %>% 
    dplyr::mutate(Year = dplyr::recode(RS.centYear, 
                                             `-0.60` = "1997",
                                             `-0.30` = "1999",
                                             `0.06` = "2005",
                                             `0.40` = "2009",
                                             `0.80` = "2016"))
  ## for CTB
  CTB_effects_year <- 
    effects::effect(term = "RS.centYear", 
                    mod = CTB.TSAfreq.randsite)
  #summary(CTB_effects_year)
  # save the effects values as a df:
    ctb_x_year <- as.data.frame(CTB_effects_year)
  
  # recode RS.centYear to our actual years:
  ctb_x_year <- ctb_x_year %>% 
    dplyr::mutate(Year = dplyr::recode(RS.centYear, 
                                             `-0.60` = "1997",
                                             `-0.30` = "1999",
                                             `0.06` = "2005",
                                             `0.40` = "2009",
                                             `0.80` = "2016"))
  ## gorgonia
  gorg_effects_year <- 
    effects::effect(term = "RS.centYear", 
                    mod = gorg.TSAfreq.randsite)
  #summary(Soft_coral_effects_year)
  # save the effects values as a df:
  gorg_x_year <- as.data.frame(gorg_effects_year)
  
  # recode RS.centYear to our actual years:
  gorg_x_year <- gorg_x_year  %>% 
    dplyr::mutate(Year = dplyr::recode(RS.centYear, 
                                             `-0.60` = "1997",
                                             `-0.30` = "1999",
                                             `0.06` = "2005",
                                             `0.40` = "2009",
                                             `0.80` = "2016"))
  ## for Sponge
  Sponge_effects_year <- 
    effects::effect(term = "RS.centYear", 
                    mod = Sponge.TSAfreq.randsite)
  #summary(Sponge_effects_year)
  # save the effects values as a df:
  sponge_x_year <- as.data.frame(Sponge_effects_year)
  
  # recode RS.centYear to our actual years:
  sponge_x_year <- sponge_x_year %>% 
    dplyr::mutate(Year = dplyr::recode(RS.centYear, 
                                             `-0.60` = "1997",
                                             `-0.30` = "1999",
                                             `0.06` = "2005",
                                             `0.40` = "2009",
                                             `0.80` = "2016"))
  
  
  ### **12b.** Use effects value dfs (from above) to plot  estimates  
   # hard coral
  hard.coral.year.plot <-
      ggplot(LogitMasterDF %>% 
               filter(newType =="Live coral"))+
         #geom_point(data = hc_x_year, 
         #           aes(x = as.numeric(Year), 
         #               y = boot::inv.logit(fit)*100), 
         #           color="blue")+
         geom_line(data=hc_x_year, 
                    aes(x = as.numeric(Year), 
                        y = boot::inv.logit(fit)*100), 
                    color = "darkblue")+
        geom_ribbon(data = hc_x_year, 
                      aes(x = as.numeric(Year), 
                          ymin = boot::inv.logit(lower)*100, 
                          ymax = boot::inv.logit(upper)*100), 
                      alpha=0.1, 
                      fill="blue")+
          geom_point(aes(x= Year, 
                          y = MeanCover),
                   color = "grey20",
                   size = 3,
                   alpha = 0.5)+
          labs(title = "HARD CORALS",
               x= NULL, 
               y="Hard Coral Cover (%)")+
    theme_bw() +
    theme (panel.grid = element_blank(),
           plot.title = element_text(hjust = 0.5),
           axis.title = element_text (size = 14),
           axis.text = element_text (size = 14))
    
  hard.coral.year.plot
  
  # Macrolage
  macroalgae.year.plot <-
    ggplot(LogitMasterDF %>% 
             filter(newType =="Macroalgae"))+
    #geom_point(data = hc_x_year, 
    #           aes(x = as.numeric(Year), 
    #               y = boot::inv.logit(fit)*100), 
    #           color="blue")+
    geom_line(data=macro_x_year, 
              aes(x = as.numeric(Year), 
                  y = boot::inv.logit(fit)*100), 
              color = "darkblue")+
    geom_ribbon(data = macro_x_year, 
                aes(x = as.numeric(Year), 
                    ymin = boot::inv.logit(lower)*100, 
                    ymax = boot::inv.logit(upper)*100), 
                alpha=0.1, 
                fill="blue")+
    geom_point(aes(x= Year, 
                   y = MeanCover),
               color = "grey20",
               size = 3,
               alpha = 0.5)+
    labs(title = "MACROALGAE",
         x=NULL, 
         y="Macroalgae Cover (%)")+
    scale_y_continuous(limits = c(0,80))+
    theme_bw() +
    theme (panel.grid = element_blank(),
           plot.title = element_text(hjust = 0.5),
           axis.title = element_text (size = 14),
           axis.text = element_text (size = 14))
  macroalgae.year.plot
  
 #CTB
  ctb.year.plot <-
    ggplot(LogitMasterDF %>% 
             filter(newType =="CTB"))+
    #geom_point(data = ctb_x_year, 
    #           aes(x = as.numeric(Year), 
    #               y = boot::inv.logit(fit)*100), 
    #           color="blue")+
    geom_line(data= ctb_x_year, 
              aes(x = as.numeric(Year), 
                  y = boot::inv.logit(fit)*100), #back transform to %
              color = "darkblue")+
    geom_ribbon(data = ctb_x_year, 
                aes(x = as.numeric(Year), 
                    ymin = boot::inv.logit(lower)*100, #back transform to %
                    ymax = boot::inv.logit(upper)*100), #back transform to %
                alpha=0.1, 
                fill="blue")+
    geom_point(aes(x= Year, 
                   y = MeanCover),
               color = "grey20",
               size = 3,
               alpha = 0.5)+
    scale_y_continuous(limits = c(0,80))+
    labs(title = "CTB",
         x=NULL, 
         y="CTB Cover (%)")+
    theme_bw() +
    theme (panel.grid = element_blank(),
           plot.title = element_text(hjust = 0.5),
           axis.title = element_text (size = 14),
           axis.text = element_text (size = 14))
  ctb.year.plot
  
  #Gorgonian
  gorg.year.plot <-
    ggplot(LogitMasterDF %>% 
             filter(newType =="Gorgonian"))+
    #geom_point(data = ctb_x_year, 
    #           aes(x = as.numeric(Year), 
    #               y = boot::inv.logit(fit)*100), 
    #           color="blue")+
    geom_line(data= gorg_x_year, 
              aes(x = as.numeric(Year), 
                  y = boot::inv.logit(fit)*100), #back transform to %
              color = "darkblue")+
    geom_ribbon(data = gorg_x_year, 
                aes(x = as.numeric(Year), 
                    ymin = boot::inv.logit(lower)*100, #back transform to %
                    ymax = boot::inv.logit(upper)*100), #back transform to %
                alpha=0.1, 
                fill="blue")+
    geom_point(aes(x= Year, 
                   y = MeanCover),
               color = "grey20",
               size = 3,
               alpha = 0.5)+
    scale_y_continuous(limits = c(0,32))+
    labs(title = "GORGONIAN",
         x=NULL, 
         y="Gorgonian Cover (%)")+
    theme_bw() +
    theme (panel.grid = element_blank(),
           plot.title = element_text(hjust = 0.5),
           axis.title = element_text (size = 14),
           axis.text = element_text (size = 14))
  gorg.year.plot
  
  #Sponge
  sponge.year.plot <-
    ggplot(LogitMasterDF %>% 
             filter(newType =="Gorgonian"))+
    #geom_point(data = ctb_x_year, 
    #           aes(x = as.numeric(Year), 
    #               y = boot::inv.logit(fit)*100), 
    #           color="blue")+
    geom_line(data= sponge_x_year, 
              aes(x = as.numeric(Year), 
                  y = boot::inv.logit(fit)*100), #back transform to %
              color = "darkblue")+
    geom_ribbon(data = sponge_x_year, 
                aes(x = as.numeric(Year), 
                    ymin = boot::inv.logit(lower)*100, #back transform to %
                    ymax = boot::inv.logit(upper)*100), #back transform to %
                alpha=0.1, 
                fill="blue")+
    geom_point(aes(x= Year, 
                   y = MeanCover),
               color = "grey20",
               size = 3,
               alpha = 0.5)+
    labs(title = "SPONGE",
         x=NULL, 
         y="SPONGE Cover (%)")+
    scale_y_continuous(limits = c(0,32))+
    theme_bw() +
    theme (panel.grid = element_blank(),
           plot.title = element_text(hjust = 0.5),
           axis.title = element_text (size = 14),
           axis.text = element_text (size = 14))
  sponge.year.plot
  
  
  ## Combine all plots in one
  prediction_combo_plot <-
      cowplot::plot_grid(hard.coral.year.plot,
                     macroalgae.year.plot,
                     ctb.year.plot,
                     gorg.year.plot,
                     sponge.year.plot)
  
  cowplot::save_plot("Figures/Fig S4. Predicted_benthic_groups.png",
            prediction_combo_plot, 
            ncol=3,
            base_height = 7,
            base_width = 3)
  
  
  ## MDS ####
  ## regroup by ID making sure appears in every year
  master_df$newID_grouped <-
    dplyr::recode_factor(master_df$ID,
                         "Acropora.cervicornis" = "Acropora spp.",
                         "Acropora.palmata" =  "Acropora spp.",
                         "Acropora.species" = "Acropora spp.",
                         "Acropora.prolifera" = "Acropora spp.",
                         "Agaricia.agaricites.complex" = "Agaricia agaricites",
                         "Agaricia.grahamae" = "Other coral species",
                         "Agaricia.humilis" = "Other coral species",
                         "Agaricia.lamarcki" = "Other coral species",
                         "Agaricia.fragilis" = "Other coral species",
                         "Agaricia.tenuifolia" = "Agaricia tenuifolia",
                         "Diploria.labyrinthiformis" = "Diploria/Pseudodiploria spp.",
                         "Diploria.strigosa" = "Diploria/Pseudodiploria spp.",
                         "Diploria.clivosa" = "Diploria/Pseudodiploria spp.",
                         "Pseudodiploria.clivosa" = "Diploria/Pseudodiploria spp.",
                         "Pseudodiploria.strigosa" = "Diploria/Pseudodiploria spp.",
                         "Orbicella.annularis" = "Orbicella spp.",
                         "Orbicella.faveolata" = "Orbicella spp.",
                         "Orbicella.franksi" = "Orbicella spp.",
                         "Porites.astreoides" = "Porites astreoides",
                         "Porites.divaricata" = "Porites spp.",
                         "Porites.furcata" = "Porites spp.",
                         "Porites.porites" = "Porites spp.",
                         "Porites.colonensis" = "Porites spp.",
                         "Siderastrea.species" =  "Siderastrea spp.",
                         "Siderastrea.radians" = "Siderastrea spp.",
                         "Siderastrea.siderea" = "Siderastrea spp.",
                         "Millepora" = "Other coral species",
                         "Millepora.complenata" = "Other coral species",
                         "Millepora.alcicornis" = "Other coral species",
                         "Mycetophyllia" = "Other coral species",
                         "Mycetophyllia.aliciae" = "Other coral species",
                         "Mycetophyllia.danaana" = "Other coral species",
                         "Mycetophyllia.ferox" = "Other coral species",
                         "Mycetophyllia.lamarckiana" = "Other coral species",
                         "Isophyllia.sinuosa" = "Other coral species",
                         "Isophyllia.rigida" = "Other coral species",
                         "Scolymia" = "Other coral species",
                         "Scolymia.cubensis" = "Other coral species",
                         "Madracis.auretenra" = "Other coral species",
                         "Madracis.decactis" = "Other coral species",
                         "Madracis.mirabilis" = "Other coral species",
                         "Colpophyllia.natans" = "Colpophyllia natans",
                         "Dendrogyra.cylindrus" = "Dendrogyra cylindrus",
                         "Dichocoenia.stokesi" = "Other coral species",
                         "Eusmilia.fastigiata" = "Other coral species",
                         "Dichocoenia.stokesii" =  "Other coral species",
                         "Favia.fragum" = "Other coral species",
                         "Leptoseris.cucullata" = "Other coral species",
                         "Manicina.areolata" = "Other coral species",
                         "Meandrina.meandrites" = "Other coral species",
                         "Montastrea.cavernosa" = "Montastrea cavernosa",
                         "Mussa.angulosa" = "Other coral species",
                         "Solenastrea.bournoni" = "Other coral species",
                         "Stephanocoenia.intersepta" = "Other coral species",
                         "Scleractinia" = "Other coral species",
                         "Amphiroa.species" = "Calcareous algae",
                         "Branching.coraline" = "Calcareous algae",
                         "Branching.coraline.antillarum" = "Calcareous algae",
                         "Calcareous" = "Calcareous algae",
                         "Dictyota" = "Fleshy macroalgae",
                         "Encrusting.coraline" = "CTB",
                         "Erect.rhodophyta" = "Corticate algae",
                         "Fleshy_Macroalgae" = "Fleshy macroalgae",
                         "Galaxaura.species" = "Corticate algae",
                         "Halimeda" = "Calcareous algae",
                         "Lobophora" = "Corticate algae", 
                         "Macroalgae" = "Fleshy macroalgae", 
                         "Padina.species" = "Corticate algae",
                         "Turf.algae" = "CTB",
                         "Ventricaria.species" = "Corticate algae",
                         "Wrangelia" = "Fleshy macroalgae",
                         "Crustose.turf.bare" = "CTB",
                         "CCA" = "CTB",
                         "Anemone" = "Other invertebrates",
                         "Annelida" = "Other invertebrates",
                         "Ascidians" = "Other invertebrates",
                         "Bacterial.mat" = "Cyanobacteria",
                         "Barnacle" = "Other invertebrates",
                         "Bivalve" = "Other invertebrates",
                         "Black.octocorals" = "Other invertebrates",
                         "Bryozoa" = "Other invertebrates",
                         "Dead" = "CTB",
                         "Corallimorph" = "Other invertebrates",
                         "Fish" = "Other",
                         "Hydroid" = "Other invertebrates",
                         "Mat.tunicate" = "Other invertebrates",
                         "N.c" = "Other",
                         "Nudibranch.species" = "Other invertebrates",
                         "Sand.sediment" = "Other",
                         "Sea.cucumber" = "Other invertebrates",
                         "Soft.coral" = "Gorgonian",
                         "Sponge" = "Sponge",
                         "Substrate" = "CTB",
                         "Tunicate" = "Other invertebrates",
                         "Unknown" = "Other",
                         "Water" = "Other",
                         "Equipment" = "Other",
                         "Zoanthid" = "Other invertebrates",
                         "Rubble" = "Other"
    )
  
  master_df <- master_df %>% 
    filter(!(newID_grouped %in% 
               c("Other",
                 "Cyanobacteria",
                 "Other invertebrates")
    ))
  
  #summary per transect
  master_df_mds  <- master_df %>% 
    group_by(Year, 
             Site, 
             Fishing_level, 
             Transect, 
             Image.Code, 
             newID_grouped) %>%
    summarise(total_cover = sum(Cover, 
                                na.rm = TRUE)) %>%
    group_by(Year, 
             Site, 
             Transect, 
             Fishing_level, 
             newID_grouped) %>%
    summarise(cover = mean(total_cover, 
                           na.rm = TRUE))
  
  #mean per site
  mean.master_df_mds <- master_df_mds %>% 
    group_by(Site, 
             Year, 
             Fishing_level, 
             newID_grouped) %>% 
    summarise(MeanCover = mean(cover, 
                               na.rm=TRUE))
  
  # spread by site
  mean.master_df_mds_wide <- mean.master_df_mds %>%
    spread(newID_grouped, 
           MeanCover, 
           fill = 0)
  
  #meam of Orbicella
  mean.master_df_mds_wide %>%  
    filter (Year == 2016) %>% 
     select (`Orbicella spp.`) %>% 
      group_by (Year) %>%
        summarise(mean(`Orbicella spp.`),
                  sd (`Orbicella spp.`))
  #mean.master_df_mds_wide <- master_df_mds %>%
    #spread(newID_grouped, 
   #        cover, 
    #       fill = 0)
  
  #View(mean.master_df_mds_wide)
  
  # Now we make a data frame with only the abundance values for each sample
  ncol.WMCID <- ncol(mean.master_df_mds_wide)
  FG.community <- mean.master_df_mds_wide[4:ncol.WMCID] # 5 for transect level data
  
  # use metaMDS and bray curtis distance
  FG.community.mds <- vegan::metaMDS(comm = FG.community,
                                     k = 2,
                                     distance = "bray", 
                                     engine = "monoMDS",
                                     try = 100,
                                     trace = 2, 
                                     autotransform = FALSE,
                                     zerodist = "ignore",
                                     expand = TRUE,
                                     plot = TRUE) 
  
  # display the results
  FG.community.mds # stress = 0.0976 
  
  # call the plot function on the x-y coords of the metaMDS output makes an MDS ordination plot:
  plot(FG.community.mds$points) 
  
  # look at the loadings for the MDS1 and MDS2 
  FG.community.mds$species 
  
  # extract the x and y coords from the MDS plot into a new df
  MDS_xy <- data.frame(FG.community.mds$points)
  
  # add Site vs Year factors to those coordinates (the first 2 columns of original data set)
  MDS_xy$Site <- mean.master_df_mds_wide$Site
  MDS_xy$Year <- mean.master_df_mds_wide$Year
  MDS_xy$Fishing_level <- mean.master_df_mds_wide$Fishing_level
  
  # create MDS_species data frame with ID names for overlaying:
  MDS_species <- tibble::rownames_to_column(
    data.frame(
      FG.community.mds$species
    ),
    var="species")
  
  ## Fit covariates to MDS  results #####
  permanova_data <-
    left_join(mean.master_df_mds_wide,
              Cover.Prot.Temp.HiiDF_summary,
              by = c("Site", "Year"))
  
  dim(permanova_data)
  dim(FG.community)
  names(permanova_data)
  
  #Select env data
  env <- permanova_data[, c("Site", "Year", "Fishing_level", "TSA_Freq", "HII_50km")]
  
  ## EnvFit
  mds_env <- vegan::envfit(FG.community.mds ~  
                             scale(Year) + 
                             scale(TSA_Freq) + 
                             scale(HII_50km) + 
                             Fishing_level,
                           strata = env$Site,
                           data = env,
                           permutations = 999, 
                           na.rm = TRUE)
  mds_env
  
  #Plot
  plot(mds_env)
  
  
  ## plot points and color by Year
  mds.year.protection <- 
    ggplot(MDS_xy,
           aes(- MDS1, MDS2,
               color = as.factor(Year),
               fill = as.factor(Year)))+
    geom_hline(yintercept = 0, lty = 3, col = "grey50")+ 
    geom_vline(xintercept = 00, lty = 3, col = "grey50")+
    geom_point(aes(pch = Fishing_level),
               cex = 3,
               alpha = 0.6)+
    geom_segment(x = 0, xend = mds_env$vectors$arrows[1,1]/2, 
                y = 0, yend = mds_env$vectors$arrows[1,2]/2,
                arrow = arrow(length = unit(0.1,"in")),
                col= "grey70") +
    annotate("text",  x = mds_env$vectors$arrows[1,1]/1.7,
                    y = mds_env$vectors$arrows[1,2]/1.7,
             color = "grey20",
             size = 3,
             label = "Year")+
    geom_segment(x = 0, xend = mds_env$vectors$arrows[2,1]/2, 
                 y = 0, yend = mds_env$vectors$arrows[2,2]/2,
                 arrow = arrow(length = unit(0.1,"in")),
                 col= "grey70") +
    annotate("text",  x = mds_env$vectors$arrows[2,1]/1.7,
             y = mds_env$vectors$arrows[2,2]/1.7,
             color = "grey20",
             size = 3,
             label = "TSA Freq")+
    geom_segment(x = 0, xend = mds_env$vectors$arrows[3,1]/2, 
                 y = 0, yend = mds_env$vectors$arrows[3,2]/2,
                 arrow = arrow(length = unit(0.1,"in")),
                 col= "grey70") +
    annotate("text",  x = mds_env$vectors$arrows[3,1]/1.7,
             y = mds_env$vectors$arrows[3,2]/1.7,
             color = "grey20",
             size = 3,
             label = "HII_50km")+
   scale_shape_manual(values = c(21,22))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10))+
    scale_color_nejm()+
    scale_fill_nejm(guide = 'none')+
    stat_ellipse(alpha = 0.3,
                 type = "t",
                 level = 0.95) +
    guides(color = guide_legend(title = "Year"),
           shape = guide_legend(title = "Fishing level")) +
    scale_x_continuous(limits = c(-1.1,1.2)) +
    #scale_y_continuous(limits = c(-1,0) +
    annotate("text", 
             x =1, y = 0.8,
             color = "grey20",
             size = 3,
             label = paste ("stress = ", 
                            round(FG.community.mds$stress,3)))
  
  mds.year.protection
  
  
  mds.species <- ggplot(MDS_species, 
                        aes(- MDS1, MDS2,
                            color = "                "),
                        guide = 'none')+
    geom_point()+
    geom_hline(yintercept = 0, lty = 3, col = "grey50")+ 
    geom_vline(xintercept = 00, lty = 3, col = "grey50")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 10))+
    geom_segment(data = MDS_species , 
                 aes(x = 0, 
                     xend = -MDS1, 
                     y = 0, 
                     yend = MDS2,
                     color = species), 
                 arrow = arrow(length = unit(0.15, "cm")), 
                 lwd = 0.2,
                 colour="grey") +
    geom_text(data = MDS_species, 
              aes(x = -MDS1*1, 
                  y = MDS2*1.2, 
                  label = species), 
              position = "nudge",
              size = 3,
              colour = "grey20",
              check_overlap = FALSE)+ ## this prevent overlaping names
    scale_color_manual(values = "#FFFFFF") +
    scale_x_continuous(limits = c(-1.1,1.2)) +
    scale_y_continuous(limits = c(-0.6,0.6)) 
  
  mds.species
  
  mds_plot <- gridExtra::grid.arrange (mds.year.protection, 
                                       mds.species, 
                                       ncol = 1)
  mds_plot
  
  ggsave("./Figures/Fig.4_mds_plot.png",
         mds_plot,
         device = "png",
         width = 6,
         height = 6.5, 
         units ="in", 
         dpi = 600)
  
  
  ### Calculating different aspects of MDS - helps us interpret data!
  ## A) stress: 
  FG.community.mds$stress 
  
  ## B) look for distance matrix w/in FG.community.mds
  FG.community.mds.dist <- FG.community.mds$dist
  
  #### PERMANOVA ####
   ## Join dataset for mds (mean.master_df_mds_wide) with covariates dataset
  
  #Set strata
  perm <- permute::how(nperm = 999)
  permute::setBlocks(perm) <- permanova_data$Site

  
  # Run adonis: 
  permanova_results <-
      vegan::adonis2(permanova_data[,c(4:21)] ~ 
                     arm::rescale(Year)+
                     arm::rescale(TSA_Freq)+
                    arm::rescale(HII_50km)+
                     Fishing_level, 
                   sqrt.dist = FALSE,
                   method = "bray", 
                   by = "terms",
                   permutations = 10000,
                   data = permanova_data)
  
  #See reults
  print(permanova_results)
  
  #save results as csv file
  write.csv(permanova_results, "Data/Processed/permanova_results.csv")
  
  
  #### BLMER FOR CORAL GENERA  #####
  
  # Join mean coral genera with covariates
  coral.genera.df.cov <-
    left_join (coral.genera.df,
               Cover.Prot.Temp.HiiDF_summary,
               by = c("Site", "Year"))
  
  unique(coral.genera.df.cov$newID_grouped)
  names(coral.genera.df.cov)
  
  ##Logit Transform MeanCover and Weight by SE    
  coral.genera.df.cov <- 
    coral.genera.df.cov %>%
          dplyr::mutate(Logit.Cover = car::logit(cover,
                                           percents = TRUE)) 
  
  ##  Rescale predictor variable and center year
  coral.genera.df.cov$centYear <- coral.genera.df.cov$Year-2005
  
  # Rescale all predictor variables, including the centered year:
  coral.genera.df.cov$RS.centYear <- arm::rescale(coral.genera.df.cov$centYear)
  coral.genera.df.cov$RS.Hii.50 <- arm::rescale(coral.genera.df.cov$HII_50km)
  coral.genera.df.cov$RS.TSA_Freq <- arm::rescale(coral.genera.df.cov$TSA_Freq)
  
  ## Run models for coral taxa and create table with results
  model_results <- tibble (term = "",
                           estimate = "",
                           std.error = "",
                           statistic = "",
                           p.value = "",
                           p.stars = "",
                           group = "")
  
  for (i in unique(coral.genera.df.cov$newID_grouped)) {
    coral.general.model <- 
        blme::blmer(Logit.Cover~
                  RS.centYear+
                  Fishing_level+
                  RS.Hii.50+
                  RS.TSA_Freq+
                  (1|Site), 
                REML = TRUE,
                data = coral.genera.df.cov %>% 
                  filter(newID_grouped == i)
        )
  
      print(summary(coral.general.model))
      print(sjPlot::plot_model(coral.general.model,
                     title = i)+
              theme_bw())
      
      model_result <- sjPlot::get_model_data(coral.general.model, 
                             type = "est",
                             show.r2 = TRUE,
                             show.intercept = TRUE) %>%
                      dplyr::select(term, 
                                    estimate, 
                                    std.error, 
                                    statistic, 
                                    p.value, 
                                    p.stars, 
                                    group)
      
      model_results <- rbind(model_results, 
                             c(i, rep('',6)),  
                             model_result,
                             c("Marginal R2", 
                               performance::r2(coral.general.model)$R2_marginal, 
                               rep('',5)),
                             c("Conditional R2", 
                               performance::r2(coral.general.model)$R2_conditional, 
                               rep('',5))
                             )
  }
  
  write_csv(model_results, "Data/Processed/coral_genera_blmer_results.csv")
  
  
  
  
  

  
 