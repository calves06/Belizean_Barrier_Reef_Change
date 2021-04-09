# Belizean_Barrier_Reef_Change
Pubic repository containing the code and data associated with analysis for Alves et al (2021): Twenty years of change in benthic communities across the Belizean Barrier Reef.  
DOI: https://doi.org/10.1101/2021.03.15.435443

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
