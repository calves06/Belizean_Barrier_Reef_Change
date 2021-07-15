# Belizean_Barrier_Reef_Change
Pubic repository containing the code and data associated with analysis for Alves et al (2021): Twenty years of change in benthic communities across the Belizean Barrier Reef.  
DOI: https://doi.org/10.1101/2021.03.15.435443  

**Authors:** Catherine Alves, Abel Valdivia, Richard B. Aronson, Nadia Bood, Karl D. Castillo, Courtney Cox, Clare Fieseler, Zachary Locklear, Melanie McField, Laura Mudge, James Umbanhowar, John F. Bruno

**Abstract:** Disease, storms, ocean warming, and pollution have caused the mass mortality of reef-building corals across the Caribbean over the last four decades. Subsequently, stony corals have been replaced by macroalgae, bacterial mats, and other invertebrates including soft corals and sponges, causing changes to the functioning of coral reef ecosystems. Here we describe changes in the absolute cover of benthic reef taxa including corals, sponges, and algae at 15 forereef sites (12-15m depth) across the Belizean Barrier Reef (BBR) from 1997 to 2016. We also tested whether Marine Protected Areas (MPAs) in which fishing was prohibited (but likely still occurred) mitigated these changes. Additionally, we determined whether ocean temperature anomalies (measured via satellite) or local human impacts (estimated using the Human Influence Index (HII), were related to changes in benthic community structure. We observed a reduction in the cover of reef-building corals including of the long-lived, massive coral *Orbicella* spp. (from 13 to 2 %), and an increase in fleshy and corticate macroalgae across most sites. These and other changes to reef benthic communities were unaffected by local protection. The covers of hard-coral taxa, including *Acropora* spp., *M. cavernosa*, *Orbicella* spp., and *Porites* spp., were negatively related to the frequency of ocean temperature anomalies. Only gorgonian cover was related, negatively, to our metric of the magnitude of local impacts, HII. Our results suggest that benthic communities along the BBR have experienced disturbances that are beyond the capacity of the current management structure to mitigate. We recommend that managers devote greater resources and capacity to enforcing and expanding existing marine protected areas and to mitigating local stressors, and most importantly that government, industry, and the public act immediately to reduce global carbon emissions. 

**The Repository contains the following code and data files:**  

[BRC_R_code](BRC_R_code.R) - is the R code for the data analysis and visualization  

[MasterDF](Data/Processed/Long.Master.Species.Groups.csv) - is the master data frame, containing the raw cover data, and the following variables (columns):  
* **Year** - Year of benthic surveys (1997, 1999, 2005, 2009, and 2016)  
* **File.Name** - Name of image file  
* **Site** - Name of study site along Belizean Barrier Reef  
* **Image.Code** - Denotes the *Year_Site_Transect_Image*, where *97_Bacalar_Chico_1_1* represents image 1 from transect 1 from 1997 at Bacalar Chico  
* **ID** - List of benthic species and identifications from image analysis  
* **Cover** - Percent cover of each ID per image, transect, site, and year. These values are whole, two-digit numbers divisible by 10, ranging from 0 - 100  
* **Specific.Type** - Denotes the benthic category of the corresponding ID (e.g. the *Acropora.cervicornis* ID corresponds to the *Hard.coral* Specific.Type)  
* **General.Type** - Denotes a more general category corresponding to the ID (e.g. *Acropora.cervicornis* ID corresponds to the *Coral* General.Type)  

[ProtectionDF](Data/Site/Belize_site_coord_protection.csv) - this data frame contains the following site-level variables:
* **Site** - Name of study site along Belizean Barrier Reef  
* **Site_Label** - Name of study site, with spaces between words  
* **Longitude** - Longitude of study site, in decimal degrees  
* **Latitude** - Latitude of study site, in decimal degrees  
* **Protection** - Protection code of study site; where NP = Not Protected, FP = Fully Protected, and GU = General Use  
* **Fishing_level** - Fishing level of study site, determined based off *Protection*; where NP and GU = Fishing, and FP = No fishing  

[Cover.Prot.Temp.HiiDF](Data/Processed/Cover.Prot.Temp.HiiDF.csv) - this data frame contains the site-level mean cover data and predictor variable data, which was used for the analyses, and contains the following variables:  
* **Year** - Year of benthic surveys (1997, 1999, 2005, 2009, and 2016)  
* **Site** - Name of study site along Belizean Barrier Reef  
* **Specific.Type** - Denotes the benthic category of the corresponding ID (e.g. the *Acropora.cervicornis* ID corresponds to the *Hard.coral* Specific.Type)  
* **General.Type** - Denotes a more general category corresponding to the ID (e.g. *Acropora.cervicornis* ID corresponds to the *Coral* General.Type)  
* **MeanCover** - Mean percent cover at site and year level (taken from raw Cover data)  
* **SE** - Standard error from mean percent cover calculation  
* **Protection.Code** - Protection code of study site; where NP = Not Protected, FP = Fully Protected, and GU = General Use  
* **SSTA_Freq** - Frequency of Sea Surface Temperature Anomalies (SSTA), or the number of instances SSTA was over 1 deg. C over the previous 52 weeks, obtained from the [Coral Reef Temperature Anomaly Database (CoRTAD)](https://www.nodc.noaa.gov/sog/cortad/) 
* **SSTA_Freq_hist** - Accumulative SSTA, represented by the number of times since the beginning of the dataset (1982) to survey year thaat SSTA was over 1 deg. C.
* **TSA_Freq** - Site-specific frequency of Thermal Stress Anomalies (TSA), defined by the number of instances TSA was over 1 deg. C over the previous 52 weeks, obtained from the [Coral Reef Temperature Anomaly Database (CoRTAD)](https://www.nodc.noaa.gov/sog/cortad/)  
* **TSA_Freq_hist** - Accumulative TSA, represented by the number of times since the beginning of the dataset (1982) to survey year thaat TSA was over 1 deg. C.
* **TSA_Freq_btw_surveys** - The frequency of TSA between survey years, defined as the number of instances since the previous survey year that TSA was over 1 deg. C.
* **Longitude** - Longitude of study site, in decimal degrees  
* **Latitude** - Latitude of study site, in decimal degrees  
* **Buffer** - Radius (in km) around study site, within which the Human Influence Index (HII) was calculated.  
* **HII.Value** - Human Influence Index (HII) value for each *Site* and *Buffer*, obtained from [NASA's Socioeconomic Data and Applications Center (SEDAC) database](https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-influence-index-geographic)  
