##--- PACKAGES ---#### 
library(tibble)
library(plyr)  
library(dplyr) 
library(tidyr)
library(data.table) 
library(sp)
library(rgbif)
library(reshape2) 
library(raster)   
library(rgeos)   
library(maptools)   
library(sf)      
library(rgdal)   
library(ggplot2)
library(gplots) 
library(ggpubr)
library(ggExtra)
library(gridExtra)
library(trend)
library(stringr)
library(FSA)
library(MASS)
library(lme4)
library(ciTools)
library(fmsb)
library(iNEXT)

##--- 1. DATA RETRIEVAL ---####
##--- 1.1 GBIF data     ---####
# "GBIF Occurrence Download 10.15468/dl.dmdxne accessed via GBIF.org on 2019-11-19"
# GBIF download key:
download_key_national <- occ_download(
  'hasGeospatialIssue = FALSE',
  'hasCoordinate = TRUE',
  'country = NO',
  type = "and"
) %>% 
  occ_download_meta

# Define function to retrieve the requested data
download_GBIF_API <- function(download_key,n_try,Sys.sleep_duration,destfile_name){
  start_time <- Sys.time()
  n_try_count <- 1
  
  download_url <- paste("http://api.gbif.org/v1/occurrence/download/request/",
                        download_key[1],sep="")
  
  try_download <- try(download.file(url=download_url,destfile=destfile_name,
                                    quiet=TRUE, mode="wb"),silent = TRUE)
  
  while (inherits(try_download, "try-error") & n_try_count < n_try) {   
    Sys.sleep(Sys.sleep_duration)
    n_try_count <- n_try_count+1
    try_download <- try(download.file(url=download_url,destfile=destfile_name,
                                      quiet=TRUE),silent = TRUE)
    print(paste("trying... Download link not ready. Time elapsed (min):",
                round(as.numeric(paste(difftime(Sys.time(),start_time, units = "mins"))),2)))
  }
}

# Call function,  create and unzip zip-file for download - IGNORE THE WARNING, IT'S FINE!
download_GBIF_API(download_key=download_key_national,destfile_name="tmp.zip",n_try=50,Sys.sleep_duration=180)   # Try 50 times, wait 3 minutes in between each try
archive_files_national <- unzip("tmp.zip", files = "NULL", list = T) 

# Get the occurrence.txt file in as a dataframe (using import from rio)
occurrence_national <- import(unzip("tmp.zip", files="occurrence.txt"),header=T,sep="\t")

# occurrence_national <- fread("occurrence_national.txt")
#paste("GBIF Occurrence Download", download_key_national[2], "accessed via GBIF.org on", Sys.Date())

# Remove "absent" records:
occurrence_national$occurrenceStatus <- as.factor(occurrence_national$occurrenceStatus)
occurrence_national <- occurrence_national[!occurrence_national$occurrenceStatus=="absent",]

# Remove records with no species information (758,950)
occurrence_national <- occurrence_national[!occurrence_national$species=="",]

# Remove records from 2019 or records with no year, due to time lag in digitising specimens
GBIF_nat_2018 <- occurrence_national[!occurrence_national$year==2019,]

# Remove potential duplicate records --> same species, date, basis of record, recorder, coordinates, and coordinate uncertainty
GBIF_nat_2018 <- GBIF_nat_2018[!duplicated(GBIF_nat_2018[,c("species","date","basisOfRecord","decimalLatitude", "decimalLongitude","coordinateUncertaintyInMeters")]),]

# Add a unique 'ID' to each record and drop unused levels:
GBIF_nat_2018$ID <- c(1:nrow(GBIF_nat_2018))
GBIF_nat_2018 <- droplevels(GBIF_nat_2018)

# Make the data spatial
coordinates(GBIF_nat_2018) <- ~decimalLongitude+decimalLatitude
proj4string(GBIF_nat_2018) <- CRS("+init=epsg:4326")
GBIF_nat_2018 <- spTransform(GBIF_nat_2018, CRS("+init=epsg:5130")) 

# Load/get '"Norway", and make sure it matches:
norway <- getData("GADM", country="NO", level=0)

norway_buff <- gBuffer(norway, width = 1000)
norway <- spTransform(norway, GBIF_nat_2018@proj4string) 
norway_buff <- spTransform(norway, GBIF_nat_2018@proj4string) 

GBIF_nat_2018 <- GBIF_nat_2018[norway_buff,]    
GBIF_nat_2018@data <- droplevels(GBIF_nat_2018@data)

# Convert to 'simple features'-object
GBIF_nat_2018<- st_as_sf(GBIF_nat_2018)

# Limit analysis to Animals, Plants and Fungi
GBIF_df <- GBIF_nat_2018[GBIF_nat_2018$kingdom=="Animalia" |
                           GBIF_nat_2018$kingdom=="Plantae" |
                           GBIF_nat_2018$kingdom=="Fungi",] 

# Retain the "HUMAN_OBSERVATION"s - this is done as we technically don't know
# if the collected specimens were collected by citizens and then handed over to an institution:
GBIF_df$basisOfRecord.f <- as.factor(GBIF_df$basisOfRecord)
GBIF_df <- GBIF_df[GBIF_df$basisOfRecord.f=="HUMAN_OBSERVATION",]

##--- 1.2 Red-list and alien species list ---####
# The Alien Species List were downloaded through the Norwegian Biodiversity Information Centre: https://artsdatabanken.no/
# The Norwegian Red List was provided in Excel format by the Norwegian Biodiversity Information Centre, and is identical to the on presented on https://artsdatabanken.no/
redlist <- fread("redlist_threatened_Norway.csv", header = T)
aliens <- fread("Fremmedartslista2018.csv", header=T)
{
  # Fix spaces in the names, which cause problems:
  names(aliens)<-stringr::str_replace_all(names(aliens), c(" " = "_" , "," = "_" ))
  # Remove the species which are alien to Svalbard (not the mainland),species which are expected to be present within 50 years,
  # are not alien at all, or not present in region, or are registered as "taxonIsEvaluatedInHigherRank"
  aliens <- aliens %>% 
    filter(!stringr::str_detect(Ekspertgruppe, '(Svalbard)'))  %>%
    filter(!stringr::str_detect(Utenfor_definisjon, 'notAlienSpecie')) %>%
    filter(!stringr::str_detect(Utenfor_definisjon, 'canNotEstablishWithin50years')) %>%
    filter(!stringr::str_detect(Utenfor_definisjon, 'establishedBefore1800')) %>%
    filter(!stringr::str_detect(Utenfor_definisjon, 'NotPresentInRegion')) %>%
    filter(!stringr::str_detect(Utenfor_definisjon, 'traditionalProductionSpecie')) %>%
    filter(!stringr::str_detect(Utenfor_definisjon, 'taxonIsEvaluatedInHigherRank'))
}

# Compare names with GBIF backbone taxonomy
nb_red <- list()
for(i in 1:nrow(redlist)){
  nb_red[i] <- list(name_backbone(name = c(redlist[i,"Vitenskapelig.navn"])))
}


nb_alien <- list()
for(i in 1:nrow(aliens)){
  nb_alien[i] <- list(name_backbone(name = c(aliens[i,"Vitenskapelig_navn"])))
} 

# Clean up mis-matches in the columns for better subsetting:
{
  # Redlist
  nb_red.df <- as.data.frame(do.call(rbind,lapply(nb_red, `length<-`, max(sapply(nb_red, length)))))  # OBS! The non-accepted names are not in the right columns!
  redlist_synonyms <- (nb_red.df[!(nb_red.df$status=="ACCEPTED" | nb_red.df$status=="DOUBTFUL"),])
  
  rownames(redlist_synonyms) <- NULL
  redlist_synonyms[c(5,84,148,196),]$kingdomKey <- redlist_synonyms[c(5,84,148,196),]$species
  redlist_synonyms[c(90),]$kingdomKey <- redlist_synonyms[c(90),]$species
  redlist_synonyms <- redlist_synonyms[!redlist_synonyms$scientificName=='NONE',]
  
  redlist_synonyms2 <- data.frame(Vitenskapelig.navn = unlist(redlist_synonyms$rank),
                                  AccSyns = as.character(unlist(redlist_synonyms$kingdomKey)))
  redlist_synonyms2$AccSyns <- as.character(redlist_synonyms2$AccSyns)
  
  # Alien list
  nb_alien.df <- as.data.frame(do.call(rbind,lapply(nb_alien, `length<-`, max(sapply(nb_alien, length)))))  # OBS! The non-accepted names are not in the right columns!
  alien_synonyms <- (nb_alien.df[!(nb_alien.df$status=="ACCEPTED" | nb_alien.df$status=="DOUBTFUL"),])
  
  rownames(alien_synonyms) <- NULL
  alien_synonyms <- alien_synonyms[!alien_synonyms$scientificName=='NONE',]
  
  alien_synonyms2 <- data.frame(Vitenskapelig_navn = unlist(alien_synonyms$rank),
                                AccSyns = as.character(unlist(alien_synonyms$kingdomKey)))
  alien_synonyms2$AccSyns <- as.character(alien_synonyms2$AccSyns)
}

# Add the accepted names/synonyms to the lists
redlist <- left_join(redlist, redlist_synonyms2)
redlist[is.na(redlist$AccSyns),]$AccSyns <- redlist[is.na(redlist$AccSyns),]$Vitenskapelig.navn

aliens <- left_join(aliens, alien_synonyms2)
aliens[is.na(aliens$AccSyns),]$AccSyns <- aliens[is.na(aliens$AccSyns),]$Vitenskapelig_navn

# Remove unneccesary files
rm(nb_red.df)
rm(redlist_synonyms2)
rm(redlist_synonyms)
rm(nb_alien.df)
rm(alien_synonyms2)
rm(alien_synonyms)

# Add to the (cleaned) GBIF data, if the species are either alien or threatened:
GBIF_df$list <- ifelse(GBIF_df$species %in% redlist$AccSyns, paste0("threat"),
                       ifelse(GBIF_df$species %in% aliens$AccSyns, paste0("alien"),
                              NA))

##--- 1.3 Subsetting of datasets ---####
# Get the citation (and thus the needed information on the dataset to be re-found in GBIF) - this information can also be found in the download's zip-folder in the .txt-file "citations"
dtsets <- data.frame(datasetKey = unique(GBIF_df$datasetKey),
                     datasetName = NA,
                     citation = NA,
                     URL = NA,
                     Animalia = NA,
                     Fungi = NA,
                     Plantae = NA)
dtsets$datasetKey <- as.character(dtsets$datasetKey)

for(i in 1:nrow(dtsets)){
  print(i)
  # Dataset name
  dtsets$datasetName[i] <- gbif_citation(x=dtsets$datasetKey[i])$citation$title
  # Full citation of dataset
  dtsets$citation[i] <- gbif_citation(x=dtsets$datasetKey[i])$citation$citation
  # Tabulation of #records in each of the three kingdoms
  tab <- as.data.frame(t(as.matrix(table(GBIF_df[GBIF_df$datasetKey==dtsets$datasetKey[i], ]$kingdom))))  # This extra step is needed as some taxa are missing in some of the datasets
  dtsets$Animalia[i] <- ifelse(is.null(tab$Animalia), 0, tab$Animalia)
  dtsets$Fungi[i] <- ifelse(is.null(tab$Fungi), 0, tab$Fungi)
  dtsets$Plantae[i] <- ifelse(is.null(tab$Plantae), 0, tab$Plantae)
  # Retrieve the URL/DOI for easier access later on
  dtsets$URL[i] <- sub(".*?dataset (.*?) accessed.*", "\\1", dtsets$citation[i])  # This retrieves the string between 'dataset ' and ' accessed' 
}

# Add total number of records 
dtsets$Total <- rowSums(dtsets[,c("Animalia","Plantae","Fungi")])
dtsets_long <- gather(dtsets, taxa, n, Animalia, Plantae, Fungi)   # For ggplotting

# Divided by datasetName
ggplot(dtsets_long[dtsets_long$Total>50000,], aes(x=reorder(datasetName, -n), y=n, fill=taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 55, hjust = 1))

# Add name of the publishing institutions (or equivalent) to the dataframe
dtsets$Publisher_name <-ifelse(dtsets$Publisher.host=="Artsdatabanken", paste0("Norwegian Biodiversity Information Centre"),
                               ifelse(dtsets$Publisher.host=="Cornell Lab of Ornithology", paste0("Cornell Lab of Ornithology"),
                                      ifelse(dtsets$Publisher.host=="BioFokus", paste0("BioFokus"),
                                             ifelse(dtsets$Publisher.host=="Biolog J.B. Jordal AS", paste0("Biologist J.B. Jordal AS"),
                                                    ifelse(dtsets$Publisher.host=="Natural History Museum, University of Oslo ", paste0("Natural History Museum, UiO"),
                                                           ifelse(dtsets$Publisher.host=="NTNU University Museum ", paste0("NTNU University Museum"),
                                                                  ifelse(dtsets$Publisher.host=="Agder Museum of Natural History and Botanical Garden ", paste0("Agder Museum of Natural History (...)"), NA)))))))

# Subset to the top-10 datasets (datasets containing >50,000 records)
dt_key <- as.character(dtsets[dtsets$Total>50000,]$datasetKey)
GBIF_top10 <- GBIF_df[GBIF_df$datasetKey %in% dt_key,]   

# Add the better names for the publishing agent and the datasets
GBIF_top10 <- left_join(GBIF_top10, dtsets[,c("datasetKey","Publisher_name")])

GBIF_top10$dataset_Name <- ifelse(GBIF_top10$datasetKey=="2e4cc37b-302e-4f1b-bbbb-1f674ff90e14", paste0("Biofokus"),
                                  ifelse(GBIF_top10$datasetKey=="4fa7b334-ce0d-4e88-aaae-2e0c138d049e", paste0("eBird"),
                                         ifelse(GBIF_top10$datasetKey=="b124e1e0-4755-430f-9eab-894f25a9b59c", paste0("NBIC_CitizenScience"),
                                                ifelse(GBIF_top10$datasetKey=="09c38deb-8674-446e-8be8-3347f6c094ef", paste0("J.B.Jordal"),
                                                       ifelse(GBIF_top10$datasetKey=="0efbb3b3-6473-414a-af91-cdbcbcd8aba1", paste0("AgderMuseum_VascPlants"),
                                                              ifelse(GBIF_top10$datasetKey=="88ea9eca-96ad-4a5a-8a01-99a4acabe050", paste0("UiO_VascPlantsObs"), 
                                                                     ifelse(GBIF_top10$datasetKey=="d34ed8a4-d3cb-473c-a11c-79c5fec4d649", paste0("UiO_VascPlantsNotes"),
                                                                            ifelse(GBIF_top10$datasetKey=="492d63a8-4978-4bc7-acd8-7d0e3ac0e744", paste0("NBIC_other"),
                                                                                   ifelse(GBIF_top10$datasetKey=="d74d5b23-922d-40ff-9120-eb742d1d66f6", paste0("UiO_Lichen"),
                                                                                          ifelse(GBIF_top10$datasetKey=="f06ed729-104f-4e51-9db3-dfd9b228a1be", paste0("NTNU-INH_VascPlants"), NA))))))))))

##--- 1.4 Land-cover data ---####
# AR50 land-cover maps (downloaded from Kartkatalogen on November 23rd 2019 (https://kartkatalog.geonorge.no/metadata/arealressurskart-ar50---arealtyper/41f6b000-c394-41c5-8ebb-07a0a3ec914f)).
# Due to computational limitations, we used a reduced AR50-dataset with only one attribute related to landcover; information on land-cover type,
# forest types were combined in a single attribute column manually in ArcGIS. Likewise, the dataset were 'dissolved' to reduce the number of polygons.
# The dataset was further cropped to the 'norway_buff' shapefile in ArcGIS
AR50_2016_sf_crop <- readOGR(dsn=path.expand("Data_files"), layer="AR50_clip_NObuff")  
AR50_2016_sf_crop<- st_as_sf(AR50_2016_sf_crop)
AR50_2016_sf_crop$m2 <- st_area(AR50_2016_sf_crop)

# Transform the CRS to match the GBIF-data
AR50_2016_sf_crop <- st_transform(AR50_2016_sf_crop, GBIF_nat_2018@proj4string)

##---------------------------####
##--- 2. TAXONOMIC SKEW ---####
##--- 2.1 Overall       ---####
ggplot(GBIF_top10, aes(x=kingdom, fill=dataset_Name)) + geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme()

chi.tax <- chisq.test(table(GBIF_top10$kingdom))
chi.tax$observed
chi.tax$expected
chi.tax$residuals   

##--- 2.2 Within kingdoms ---####
# Animals  (to save time later, subset the data once: GBIF_animal <- as.data.frame(GBIF_top10[GBIF_top10$kingdom=="Animalia",])   )
ggplot(GBIF_animal, aes(x=class, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chisq.test(table(GBIF_animal$class))$residuals    

# Plants  (to save time later, subset the data once: GBIF_plant <- as.data.frame(GBIF_top10[GBIF_top10$kingdom=="Plantae",])   )
ggplot(GBIF_plant, aes(x=phylum, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chisq.test(table(GBIF_plant$phylum))$residuals  
# Vascular plants:
ggplot(GBIF_plant[GBIF_plant$phylum=="Tracheophyta",], aes(x=order, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chisq.test(table(GBIF_plant[GBIF_plant$phylum=="Tracheophyta",]$order))$residuals    

# Fungi  (to save time later, subset the data once: GBIF_fungi <- as.data.frame(GBIF_top10[GBIF_top10$kingdom=="Fungi",])    )
ggplot(GBIF_fungi, aes(x=phylum, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
chisq.test(table(GBIF_fungi$phylum))$residuals  


##--- 2.3 Between the different datasets (between and within kingdoms) ---####
chisq.test(table(GBIF_top10$dataset_Name))

ggplot(GBIF_top10, aes(x=dataset_Name, fill=kingdom)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("#F8766D","#619CFF","#00BA38")) +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

chi.tax_pub <- chisq.test(table(GBIF_top10$dataset_Name,
                                GBIF_top10$kingdom))
chi.tax_pub$observed
chi.tax_pub$expected
chi.tax_pub$residuals

# Animals
chisq.test(table(GBIF_animal$dataset_Name,
                 GBIF_animal$class))$residuals  

ggplot(GBIF_animal, aes(x=class, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ dataset_Name, scales = "free")

# Plants
chisq.test(table(GBIF_plant$datasetKey,
                 GBIF_plant$phylum))$residuals  

ggplot(GBIF_plant, aes(x=phylum, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Vascular Plants
chisq.test(table(GBIF_plant[GBIF_plant$phylum=="Tracheophyta",]$dataset_Name,
                 GBIF_plant[GBIF_plant$phylum=="Tracheophyta",]$order))$residuals 

ggplot(GBIF_plant[GBIF_plant$phylum=="Tracheophyta" ,], aes(x=order, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##--- 2.4 Between threatened and alien species (between and within kingdoms) ---####
ggplot(GBIF_top10, aes(x=list, fill=kingdom)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("#F8766D","#619CFF","#00BA38")) +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

chi.tax_list <- chisq.test(table(GBIF_top10$list,
                                 GBIF_top10$kingdom))
chi.tax_list$observed
chi.tax_list$expected
chi.tax_list$residuals

chi.tp_list <- chisq.test(table(GBIF_top10$dataset_Name,
                                GBIF_top10$list))

chi.tp_list
chi.tp_list$observed
chi.tp_list$expected
chi.tp_list$residuals

ggplot(GBIF_top10, aes(x=dataset_Name, fill=list)) +
  geom_bar(stat = "count") +
  scale_fill_manual(values = c("gray20","#F8766D","gray90")) +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

chi_list <- chisq.test(table(GBIF_top10$list))
chi_list
chi_list$observed
chi_list$expected
chi_list$residuals

##--------------------------------------------####
##--- 3. TEMPORAL SKEW   ---####
##--- 3.1 Overall        ---####
# All records
ggplot(GBIF_top10, aes(x=year)) +
  geom_bar(stat = "count") +
  #geom_density(alpha = 0.25) +    
  theme(axis.text.x = element_text(angle=45))

# Mann-Kendall test - non-parametric test for a monotonic trend
# in the data:
mk <- mk.test(as.data.frame(table(GBIF_top10$year))$Freq)  

##--- 3.2 Between datasets ---####
# Development of no.records over time between the datasets. Make a better dataframe for plotting and testing:
dat.year <- data.frame(year=c(GBIF_top10[GBIF_top10$dataset_Name=="AgderMuseum_VascPlants",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="Biofokus",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="eBird",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="J.B.Jordal",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="NBIC_CitizenScience",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="NBIC_other",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="NTNU-INH_VascPlants",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="UiO_Lichen",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsNotes",]$year,
                              GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsObs",]$year),
                       dataset=factor(c(rep("AgderMuseum_VascPlants", length(GBIF_top10[GBIF_top10$dataset_Name=="AgderMuseum_VascPlants",]$year)),
                                        rep("Biofokus", length(GBIF_top10[GBIF_top10$dataset_Name=="Biofokus",]$year)),
                                        rep("eBird", length(GBIF_top10[GBIF_top10$dataset_Name=="eBird",]$year)),
                                        rep("J.B.Jordal", length(GBIF_top10[GBIF_top10$dataset_Name=="J.B.Jordal",]$year)),
                                        rep("NBIC_CitizenScience", length(GBIF_top10[GBIF_top10$dataset_Name=="NBIC_CitizenScience",]$year)),
                                        rep("NBIC_other", length(GBIF_top10[GBIF_top10$dataset_Name=="NBIC_other",]$year)),
                                        rep("NTNU-INH_VascPlants", length(GBIF_top10[GBIF_top10$dataset_Name=="NTNU-INH_VascPlants",]$year)),
                                        rep("UiO_Lichen", length(GBIF_top10[GBIF_top10$dataset_Name=="UiO_Lichen",]$year)),
                                        rep("UiO_VascPlantsNotes", length(GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsNotes",]$year)),
                                        rep("UiO_VascPlantsObs", length(GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsObs",]$year))),
                                      levels = c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other",
                                                 "NTNU-INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")))

ggplot(dat.year, aes(x=year, color=dataset, fill=dataset, linetype=dataset)) +
  geom_density(alpha=0.25)
ggplot(dat.year, aes(x=year, color=NULL, fill=dataset)) +
  geom_bar(position = "dodge", alpha=0.75) 

# Means and medians:
dplyr::group_by(dat.year, dataset) %>%  # Package specification needed
  dplyr::summarise(
    count = dplyr::n(),
    mean = mean(year, na.rm=T),
    sd = sd(year, na.rm=T),
    median = median(year, na.rm=T),
    IQR = IQR(year, na.rm=T) )

# Convert dataframe format for testing:
dat.year2 <- (as.data.frame.matrix(table(dat.year$year, dat.year$dataset)))
dat.year2 <- tibble::rownames_to_column(dat.year2, var="year")
dat.year2$year <- as.integer(dat.year2$year)
dat.year2 <- melt(dat.year2, id=c("year"))
colnames(dat.year2) <- c("year","dataset", "nrec")
plot(nrec~year, data=dat.year2, col=dat.year2$dataset)

# The data does not meet the assumptions for an ANOVA, therefore Kruskal Wallis test for unmatched samples (non-parametric due to non-normality)
boxplot(year~dataset, data=dat.year, las=2)
kruskal.test(nrec~dataset, data=dat.year2) 
dunnTest(year ~ dataset, data=dat.year, method="bonferroni")

##--- 3.3 Temporal within-year variation ---####
table(GBIF_top10$month, GBIF_top10$dataset_Name)

ggplot(data=GBIF_top10, aes(as.factor(month), fill=kingdom)) +
  geom_bar(stat="count") +
  facet_wrap(~dataset_Name + kingdom) 
ggplot(data=GBIF_top10, aes(as.factor(month), fill=kingdom)) +
  geom_bar(stat="count") +
  facet_wrap(~dataset_Name)

# Create dataframes for radar-plots. Ppreliminary plots showed a spike in January for animals (in all groups) as records with no information is given 1/1 as
# recording date - a independent plot leaving out these are constructed as well

{radar_animal_excl <- GBIF_animal[!((GBIF_animal$month==1 & GBIF_animal$day==1)) ,]  # Removing january 1st
radar_animal_excl <- as.data.frame(table(as.factor(radar_animal_excl$month), radar_animal_excl$dataset_Name))
colnames(radar_animal_excl) <- c("month","dataset_Name","Freq")
radar_animal_excl <- dcast(radar_animal_excl, month ~ dataset_Name, value.var = "Freq")
rownames(radar_animal_excl) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
radar_animal_excl$Biofokus <- radar_animal_excl$Biofokus / sum(radar_animal_excl$Biofokus)
radar_animal_excl$eBird <- radar_animal_excl$eBird / sum(radar_animal_excl$eBird)
radar_animal_excl$NBIC_CitizenScience <- radar_animal_excl$NBIC_CitizenScience / sum(radar_animal_excl$NBIC_CitizenScience)
radar_animal_excl$J.B.Jordal <- radar_animal_excl$J.B.Jordal / sum(radar_animal_excl$J.B.Jordal)
radar_animal_excl$NBIC_other <- radar_animal_excl$NBIC_other / sum(radar_animal_excl$NBIC_other)
radar_animal_excl[,"month"] <- NULL
radar_animal_excl <- as.data.frame(t(radar_animal_excl))
radar_animal_excl <- rbind(rep(0.5,12), rep(0,12), radar_animal_excl)  # Use 0 and 0.5 as min and max to make the plots easier to interpret
}
{radar_animal <- as.data.frame(table(as.factor(GBIF_animal$month), GBIF_animal$dataset_Name))
colnames(radar_animal) <- c("month","dataset_Name","Freq")
radar_animal <- dcast(radar_animal, month ~ dataset_Name, value.var = "Freq")
rownames(radar_animal) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
radar_animal$Biofokus <- radar_animal$Biofokus / sum(radar_animal$Biofokus)
radar_animal$eBird <- radar_animal$eBird / sum(radar_animal$eBird)
radar_animal$NBIC_CitizenScience <- radar_animal$NBIC_CitizenScience / sum(radar_animal$NBIC_CitizenScience)
radar_animal$J.B.Jordal <- radar_animal$J.B.Jordal / sum(radar_animal$J.B.Jordal)
radar_animal$NBIC_other <- radar_animal$NBIC_other / sum(radar_animal$NBIC_other)
radar_animal[,"month"] <- NULL
radar_animal <- as.data.frame(t(radar_animal))
radar_animal <- rbind(rep(0.5,12), rep(0,12), radar_animal)}  
{radar_plant <- as.data.frame(table(as.factor(GBIF_plant$month), GBIF_plant$dataset_Name))
colnames(radar_plant) <- c("month","dataset_Name","Freq")
radar_plant <- dcast(radar_plant, month ~ dataset_Name, value.var = "Freq")
rownames(radar_plant) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
names(radar_plant) <- c("month","AgderMuseum_VascPlants","Biofokus","J.B.Jordal","NBIC_CitizenScience","NBIC_other","NTNU_INH_VascPlants","UiO_VascPlantsNotes","UiO_VascPlantsObs")
radar_plant$AgderMuseum_VascPlants <- radar_plant$AgderMuseum_VascPlants / sum(radar_plant$AgderMuseum_VascPlants)
radar_plant$Biofokus <- radar_plant$Biofokus / sum(radar_plant$Biofokus)
radar_plant$J.B.Jordal <- radar_plant$J.B.Jordal / sum(radar_plant$J.B.Jordal)
radar_plant$NBIC_CitizenScience <- radar_plant$NBIC_CitizenScience / sum(radar_plant$NBIC_CitizenScience)
radar_plant$NBIC_other <- radar_plant$NBIC_other / sum(radar_plant$NBIC_other)
radar_plant$NTNU_INH_VascPlants <- radar_plant$NTNU_INH_VascPlants / sum(radar_plant$NTNU_INH_VascPlants)
radar_plant$UiO_VascPlantsNotes <- radar_plant$UiO_VascPlantsNotes / sum(radar_plant$UiO_VascPlantsNotes)
radar_plant$UiO_VascPlantsObs <- radar_plant$UiO_VascPlantsObs / sum(radar_plant$UiO_VascPlantsObs)
radar_plant[,"month"] <- NULL
radar_plant <- as.data.frame(t(radar_plant))
radar_plant <- rbind(rep(0.5,12), rep(0,12), radar_plant)}  
{radar_fungi <- as.data.frame(table(as.factor(GBIF_fungi$month), GBIF_fungi$dataset_Name))
colnames(radar_fungi) <- c("month","dataset_Name","Freq")
radar_fungi <- dcast(radar_fungi, month ~ dataset_Name, value.var = "Freq")
rownames(radar_fungi) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
radar_fungi$Biofokus <- radar_fungi$Biofokus / sum(radar_fungi$Biofokus)
radar_fungi$J.B.Jordal <- radar_fungi$J.B.Jordal / sum(radar_fungi$J.B.Jordal)
radar_fungi$NBIC_CitizenScience <- radar_fungi$NBIC_CitizenScience / sum(radar_fungi$NBIC_CitizenScience)
radar_fungi$NBIC_other <- radar_fungi$NBIC_other / sum(radar_fungi$NBIC_other)
radar_fungi$UiO_Lichen <- radar_fungi$UiO_Lichen / sum(radar_fungi$UiO_Lichen)
radar_fungi[,"month"] <- NULL
radar_fungi <- as.data.frame(t(radar_fungi))
radar_fungi <- rbind(rep(0.5,12), rep(0,12), radar_fungi)} 
{radar_all <- as.data.frame(table(as.factor(GBIF_top10$month), GBIF_top10$dataset_Name))
colnames(radar_all) <- c("month","dataset_Name","Freq")
radar_all <- dcast(radar_all, month ~ dataset_Name, value.var = "Freq")
rownames(radar_all) <- c("January","February","March","April","May","June","July","August","September","October","November","December")
names(radar_all)[names(radar_all) == 'NTNU-INH_VascPlants'] <- 'NTNU_INH_VascPlants'   # Fix a name from hyphen to underscore
radar_all$AgderMuseum_VascPlants <- radar_all$AgderMuseum_VascPlants / sum(radar_all$AgderMuseum_VascPlants)
radar_all$Biofokus <- radar_all$Biofokus / sum(radar_all$Biofokus)
radar_all$eBird <- radar_all$eBird / sum(radar_all$eBird)
radar_all$J.B.Jordal <- radar_all$J.B.Jordal / sum(radar_all$J.B.Jordal)
radar_all$NBIC_CitizenScience <- radar_all$NBIC_CitizenScience / sum(radar_all$NBIC_CitizenScience)
radar_all$NBIC_other <- radar_all$NBIC_other / sum(radar_all$NBIC_other)
radar_all$NTNU_INH_VascPlants <- radar_all$NTNU_INH_VascPlants / sum(radar_all$NTNU_INH_VascPlants)
radar_all$UiO_Lichen <- radar_all$UiO_Lichen / sum(radar_all$UiO_Lichen)
radar_all$UiO_VascPlantsNotes <- radar_all$UiO_VascPlantsNotes / sum(radar_all$UiO_VascPlantsNotes)
radar_all$UiO_VascPlantsObs <- radar_all$UiO_VascPlantsObs / sum(radar_all$UiO_VascPlantsObs)
radar_all[,"month"] <- NULL
radar_all <- as.data.frame(t(radar_all))
radar_all <- rbind(rep(0.5,12), rep(0,12), radar_all)}  
layout(matrix(c(1, 2, 5,
                3, 4, 6), nrow=2, byrow=TRUE), widths = c(1,1,1))
par(mar=c(1,1,3,1))
radarchart(radar_animal, pcol = c("#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c"), plwd = 2, cglcol="grey", cglty = 1, seg=5, axislabcol="grey", vlcex = 0.75, title = "Animals")
radarchart(radar_animal_excl, pcol = c("#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c"), plwd = 2, cglcol="grey", cglty = 1, seg=5, axislabcol="grey", vlcex = 0.75, title = "Animals (excl. 1/1)")
radarchart(radar_plant, pcol = c("#a6cee3","#1f78b4","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#cab2d6","#6a3d9a"), plwd = 2, cglcol="grey", cglty = 1, seg=5, axislabcol="grey", vlcex = 0.75, title = "Plants")
radarchart(radar_fungi, pcol = c("#1f78b4","#33a02c","#fb9a99","#e31a1c","#ff7f00"), plwd = 2, cglcol="grey", cglty = 1, seg=5, axislabcol="grey", vlcex = 0.75, title = "Fungi")
radarchart(radar_all, pcol = c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"), plwd = 2, cglcol="grey", cglty = 1, seg=5, axislabcol="grey", vlcex = 0.75, title = "All records")
plot(0,type='n',axes=FALSE,ann=FALSE)
legend("center", legend = c("AgderMuseum_VascPlants","Biofokus","eBird","J.B.Jordal","NBIC_CitizenScience","NBIC_other","NTNU_INH_VascPlants","UiO_Lichen","UiO_VascPlantsNotes","UiO_VascPlantsObs"),
       bty="n", pch=20, col=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"),
       text.col = "black", cex=.75, pt.cex=2, y.intersp = 2)

# Testing within-year variation - Kruskal-Wallis Test. Construct dataframes in format appropriate for plotting and testing
dat.month_animal <- {data.frame(month=c(GBIF_animal[GBIF_animal$dataset_Name=="Biofokus",]$month,
                                        GBIF_animal[GBIF_animal$dataset_Name=="eBird",]$month,
                                        GBIF_animal[GBIF_animal$dataset_Name=="J.B.Jordal",]$month,
                                        GBIF_animal[GBIF_animal$dataset_Name=="NBIC_CitizenScience",]$month,
                                        GBIF_animal[GBIF_animal$dataset_Name=="NBIC_other",]$month),
                                dataset=factor(c(rep("Biofokus", length(GBIF_animal[GBIF_animal$dataset_Name=="Biofokus",]$month)),
                                                 rep("eBird", length(GBIF_animal[GBIF_animal$dataset_Name=="eBird",]$month)),
                                                 rep("J.B.Jordal", length(GBIF_animal[GBIF_animal$dataset_Name=="J.B.Jordal",]$month)),
                                                 rep("NBIC_CitizenScience", length(GBIF_animal[GBIF_animal$dataset_Name=="NBIC_CitizenScience",]$month)),
                                                 rep("NBIC_other", length(GBIF_animal[GBIF_animal$dataset_Name=="NBIC_other",]$month))),
                                               levels = c("Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other")))}
dat.month_plant <- {data.frame(month=c(GBIF_plant[GBIF_plant$dataset_Name=="AgderMuseum_VascPlants",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="Biofokus",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="J.B.Jordal",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="NBIC_CitizenScience",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="NBIC_other",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="NTNU-INH_VascPlants",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="UiO_VascPlantsNotes",]$month,
                                       GBIF_plant[GBIF_plant$dataset_Name=="UiO_VascPlantsObs",]$month),
                               dataset=factor(c(rep("AgderMuseum_VascPlants", length(GBIF_plant[GBIF_plant$dataset_Name=="AgderMuseum_VascPlants",]$month)),
                                                rep("Biofokus", length(GBIF_plant[GBIF_plant$dataset_Name=="Biofokus",]$month)),
                                                rep("J.B.Jordal", length(GBIF_plant[GBIF_plant$dataset_Name=="J.B.Jordal",]$month)),
                                                rep("NBIC_CitizenScience", length(GBIF_plant[GBIF_plant$dataset_Name=="NBIC_CitizenScience",]$month)),
                                                rep("NBIC_other", length(GBIF_plant[GBIF_plant$dataset_Name=="NBIC_other",]$month)),
                                                rep("NTNU-INH_VascPlants", length(GBIF_plant[GBIF_plant$dataset_Name=="NTNU-INH_VascPlants",]$month)),
                                                rep("UiO_VascPlantsNotes", length(GBIF_plant[GBIF_plant$dataset_Name=="UiO_VascPlantsNotes",]$month)),
                                                rep("UiO_VascPlantsObs", length(GBIF_plant[GBIF_plant$dataset_Name=="UiO_VascPlantsObs",]$month))),
                                              levels = c("AgderMuseum_VascPlants", "Biofokus", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other",
                                                         "NTNU-INH_VascPlants", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")))}
dat.month_fungi <- {data.frame(month=c(GBIF_fungi[GBIF_fungi$dataset_Name=="Biofokus",]$month,
                                       GBIF_fungi[GBIF_fungi$dataset_Name=="J.B.Jordal",]$month,
                                       GBIF_fungi[GBIF_fungi$dataset_Name=="NBIC_CitizenScience",]$month,
                                       GBIF_fungi[GBIF_fungi$dataset_Name=="NBIC_other",]$month,
                                       GBIF_fungi[GBIF_fungi$dataset_Name=="UiO_Lichen",]$month),
                               dataset=factor(c(rep("Biofokus", length(GBIF_fungi[GBIF_fungi$dataset_Name=="Biofokus",]$month)),
                                                rep("J.B.Jordal", length(GBIF_fungi[GBIF_fungi$dataset_Name=="J.B.Jordal",]$month)),
                                                rep("NBIC_CitizenScience", length(GBIF_fungi[GBIF_fungi$dataset_Name=="NBIC_CitizenScience",]$month)),
                                                rep("NBIC_other", length(GBIF_fungi[GBIF_fungi$dataset_Name=="NBIC_other",]$month)),
                                                rep("UiO_Lichen", length(GBIF_fungi[GBIF_fungi$dataset_Name=="UiO_Lichen",]$month))),
                                              levels = c("Biofokus", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other",
                                                         "UiO_Lichen")))}
dat.month_all <- {data.frame(month=c(GBIF_top10[GBIF_top10$dataset_Name=="AgderMuseum_VascPlants",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="Biofokus",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="eBird",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="J.B.Jordal",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="NBIC_CitizenScience",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="NBIC_other",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="NTNU-INH_VascPlants",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="UiO_Lichen",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsNotes",]$month,
                                     GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsObs",]$month),
                             dataset=factor(c(rep("AgderMuseum_VascPlants", length(GBIF_top10[GBIF_top10$dataset_Name=="AgderMuseum_VascPlants",]$month)),
                                              rep("Biofokus", length(GBIF_top10[GBIF_top10$dataset_Name=="Biofokus",]$month)),
                                              rep("eBird", length(GBIF_top10[GBIF_top10$dataset_Name=="eBird",]$month)),
                                              rep("J.B.Jordal", length(GBIF_top10[GBIF_top10$dataset_Name=="J.B.Jordal",]$month)),
                                              rep("NBIC_CitizenScience", length(GBIF_top10[GBIF_top10$dataset_Name=="NBIC_CitizenScience",]$month)),
                                              rep("NBIC_other", length(GBIF_top10[GBIF_top10$dataset_Name=="NBIC_other",]$month)),
                                              rep("NTNU-INH_VascPlants", length(GBIF_top10[GBIF_top10$dataset_Name=="NTNU-INH_VascPlants",]$month)),
                                              rep("UiO_Lichen", length(GBIF_top10[GBIF_top10$dataset_Name=="UiO_Lichen",]$month)),
                                              rep("UiO_VascPlantsNotes", length(GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsNotes",]$month)),
                                              rep("UiO_VascPlantsObs", length(GBIF_top10[GBIF_top10$dataset_Name=="UiO_VascPlantsObs",]$month))),
                                            levels = c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other",
                                                       "NTNU-INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")))}
{dat.month_animal2 <- (as.data.frame.matrix(table(dat.month_animal$month, dat.month_animal$dataset)))
  dat.month_animal2 <- tibble::rownames_to_column(dat.month_animal2, var="month")
  dat.month_animal2$month <- as.integer(dat.month_animal2$month)
  dat.month_animal2 <- melt(dat.month_animal2, id=c("month"))
  colnames(dat.month_animal2) <- c("month","dataset", "nrec")
  plot(nrec~month, data=dat.month_animal2, col=dat.year2$dat.month_animal2)}
{dat.month_plant2 <- (as.data.frame.matrix(table(dat.month_plant$month, dat.month_plant$dataset)))
  dat.month_plant2 <- tibble::rownames_to_column(dat.month_plant2, var="month")
  dat.month_plant2$month <- as.integer(dat.month_plant2$month)
  dat.month_plant2 <- melt(dat.month_plant2, id=c("month"))
  colnames(dat.month_plant2) <- c("month","dataset", "nrec")
  plot(nrec~month, data=dat.month_plant2, col=dat.month_plant2$dataset)}
{dat.month_fungi2 <- (as.data.frame.matrix(table(dat.month_fungi$month, dat.month_fungi$dataset)))
  dat.month_fungi2 <- tibble::rownames_to_column(dat.month_fungi2, var="month")
  dat.month_fungi2$month <- as.integer(dat.month_fungi2$month)
  dat.month_fungi2 <- melt(dat.month_fungi2, id=c("month"))
  colnames(dat.month_fungi2) <- c("month","dataset", "nrec")
  plot(nrec~month, data=dat.month_fungi2, col=dat.month_fungi2$dataset)}
{dat.month_all2 <- (as.data.frame.matrix(table(dat.month_all$month, dat.month_all$dataset)))
  dat.month_all2 <- tibble::rownames_to_column(dat.month_all2, var="month")
  dat.month_all2$month <- as.integer(dat.month_all2$month)
  dat.month_all2 <- melt(dat.month_all2, id=c("month"))
  colnames(dat.month_all2) <- c("month","dataset", "nrec")
  plot(nrec~month, data=dat.month_all2, col=dat.month_all2$dataset)}

dplyr::group_by(dat.month_animal, dataset) %>%  # Package specification needed, for some reason
  dplyr::summarise(
    count = dplyr::n(),
    mean = mean(month, na.rm=T),
    sd = sd(month, na.rm=T),
    median = median(month, na.rm=T),
    IQR = IQR(month, na.rm=T) )

dplyr::group_by(dat.month_plant, dataset) %>%  
  dplyr::summarise(
    count = dplyr::n(),
    mean = mean(month, na.rm=T),
    sd = sd(month, na.rm=T),
    median = median(month, na.rm=T),
    IQR = IQR(month, na.rm=T) )

dplyr::group_by(dat.month_fungi, dataset) %>%  
  dplyr::summarise(
    count = dplyr::n(),
    mean = mean(month, na.rm=T),
    sd = sd(month, na.rm=T),
    median = median(month, na.rm=T),
    IQR = IQR(month, na.rm=T) )

dplyr::group_by(dat.month_all, dataset) %>% 
  dplyr::summarise(
    count = dplyr::n(),
    mean = mean(month, na.rm=T),
    sd = sd(month, na.rm=T),
    median = median(month, na.rm=T),
    IQR = IQR(month, na.rm=T) )

# The data does not meet the assumptions for an ANOVA
# Kruskal Wallis test for unmatched samples (non-parametric due to non-normality)
boxplot(month~dataset, data=dat.month_animal, las=2)
boxplot(month~dataset, data=dat.month_plant, las=2)
boxplot(month~dataset, data=dat.month_fungi, las=2)
boxplot(month~dataset, data=dat.month_all, las=2)

kruskal.test(nrec~dataset, data=dat.month_animal2)   
kruskal.test(nrec~dataset, data=dat.month_plant2)    
kruskal.test(nrec~dataset, data=dat.month_fungi2)    
kruskal.test(nrec~dataset, data=dat.month_all2)     

dunnTest(month ~ dataset, data=dat.month_animal, method="bonferroni")
dunnTest(month ~ dataset, data=dat.month_plant, method="bonferroni")
dunnTest(month ~ dataset, data=dat.month_fungi, method="bonferroni")
dunnTest(month ~ dataset, data=dat.month_all, method="bonferroni")




##----------------------------####
##--- 4. SKEW IN LAND-USE  ---####
##--- 4.1 Assign land-cover data to points  ---####
# Assess whether some land cover types have more records than others.
# Assign land-cover data to the datapoints:

GBIF_top10 <- st_join(st_as_sf(GBIF_top10[,c("ID","year","species","list","Publisher_name","dataset_Name","geometry")]),
                                st_as_sf(AR50_2016_sf_crop[,"artype_v2"]), join=st_intersects)
GBIF_top10$cvr <- GBIF_top10$artype_v2  # Rename column

# Add better names to the land-cover:
GBIF_top10$cvr <- {factor(ifelse(GBIF_top10$artype_v2=="10", paste0("Developed"),
                                 ifelse(GBIF_top10$artype_v2=="20", paste0("Agriculture"),
                                        ifelse(GBIF_top10$artype_v2=="24", paste0("Agri.cultivated"),
                                               ifelse(GBIF_top10$artype_v2=="25", paste0("Agri.grazing"),
                                                      ifelse(GBIF_top10$artype_v2=="30", paste0("Forest"),
                                                             ifelse(GBIF_top10$artype_v2=="31", paste0("Forest.conif"),
                                                                    ifelse(GBIF_top10$artype_v2=="32", paste0("Forest.decid"),
                                                                           ifelse(GBIF_top10$artype_v2=="33", paste0("Forest.mix"),
                                                                                  ifelse(GBIF_top10$artype_v2=="50", paste0("Firmground"),
                                                                                         ifelse(GBIF_top10$artype_v2=="60", paste0("Mire"),
                                                                                                ifelse(GBIF_top10$artype_v2=="70", paste0("Snow.ice"),
                                                                                                       ifelse(GBIF_top10$artype_v2=="81", paste0("Freshwater"),
                                                                                                              ifelse(GBIF_top10$artype_v2=="82", paste0("Ocean"),
                                                                                                                     ifelse(GBIF_top10$artype_v2=="99", paste0("Not.mapped"), NA)))))))))))))))}


ggplot(GBIF_top10, aes(x=cvr, fill=dataset_Name)) +
  geom_bar(stat = "count") +
  scale_y_continuous(expand = c(0,0))   +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sum area of each land-cover type
sum_area <- stats::aggregate(AR50_2016_sf_crop$m2, by=list(AR50_2016_sf_crop$artype_v2), FUN=sum)
sum_area$area.prop <- sum_area$area / sum(sum_area$area)  # Calculate proportion of area rather than raw area measure

# Is the number of records from each land-cover type linearly correlated with area.  As the landcover map was latest updated in 2016, the data should be subsetted further.
# The chosen option is to use a 15 year period (2004-2018, both incl.):
GBIF_top10_2004 <- GBIF_top10[GBIF_top10$year>=2004,]

# OBS! Some of the points have 'NA' as cvr - rename these to '99' to be certain that they are classified as unknown rather than discarded
GBIF_top10_2004[is.na(GBIF_top10_2004$cvr),"cvr"] <- 99

sum_area$nrec <- as.data.frame(table(GBIF_top10_2004$cvr))$Freq  # Add numer of records found within each land-cover type

ggplot(data=sum_area, aes(x=cvr, y=area)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(vjust=.75, angle = 90))
ggplot(data=sum_area, aes(x=area, y=nrec, color=cvr)) +
  geom_point(shape=20, size=3) +
  stat_smooth(method = 'lm', color="red", size=0.5, linetype="dashed") +
  scale_color_manual(values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen",
                              "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) 

# Calculate the number of records from each dataset within in each land-cover type
sum_area$nrec_Biofokus <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="Biofokus",]$cvr))$Freq
sum_area$nrec_eBird <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="eBird",]$cvr))$Freq
sum_area$nrec_NBIC_CitizenScience <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="NBIC_CitizenScience",]$cvr))$Freq
sum_area$nrec_J.B.Jordal <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="J.B.Jordal",]$cvr))$Freq
sum_area$nrec_AgderMuseum_VascPlants <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="AgderMuseum_VascPlants",]$cvr))$Freq
sum_area$nrec_UiO_VascPlantsObs <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="UiO_VascPlantsObs",]$cvr))$Freq
sum_area$nrec_UiO_VascPlantsNotes <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="UiO_VascPlantsNotes",]$cvr))$Freq
sum_area$nrec_NBIC_other <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="NBIC_other",]$cvr))$Freq
sum_area$nrec_UiO_Lichen <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="UiO_Lichen",]$cvr))$Freq
sum_area$nrec_NTNU_INH_VascPlants <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$dataset_Name=="NTNU-INH_VascPlants",]$cvr))$Freq

# Calculate number of Red-listed and alien records within each land-cover type
sum_area$nrec_threat <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$list=="threat",]$cvr))$Freq
sum_area$nrec_alien <- as.data.frame(table(st_drop_geometry(GBIF_top10_2004)[st_drop_geometry(GBIF_top10_2004)$list=="alien",]$cvr))$Freq





##--- 4.2 Simulation//random placement of points ---####
# Retrieve numbers needed in the loop:
nrow(GBIF_top10_2004) 
table(GBIF_top10_2004$dataset_Name)
table(GBIF_top10_2004$list)    # Number of NA's: nrow(GBIF_top10_2004) - sum(table(GBIF_top10_2004$list))
table(GBIF_top10_2004$list, GBIF_top10_2004$dataset_Name)  # OBS! The sums do not add up because the 'NA's are not counted - these are the remaining ones
View(as.data.frame(table(GBIF_top10_2004$dataset_Name, GBIF_top10_2004$list, GBIF_top10_2004$cvr)))  # Observed distribution

for(i in 1:100){
  print(i)    # Print progress (simulation number)
  df <- st_sample(norway_buff, nrow(GBIF_top10_2004), type="random", exact=TRUE) # Create random points 
  df <- st_join(st_as_sf(df),       # Overlay the random points with the land-cover map:
                st_as_sf(AR50_2016_sf_crop[,"artype_v2"]),
                join=st_intersects) 
  df$dataset_Name <- NA    # Assign dataset_Name in the same numbers as the observed dataset
  df$dataset_Name[sample(1:nrow(df), nrow(df), FALSE)] <- rep(c("AgderMuseum_VascPlants","Biofokus","eBird","J.B.Jordal","NBIC_CitizenScience","NBIC_other","NTNU_INH_VascPlants","UiO_Lichen","UiO_VascPlantsNotes","UiO_VascPlantsObs"),
                                                              c(13373,282252,77373,111452,4465574,530913,8952,5083,44760,82528))
  df$list <- NA   # Assign list status - this is slightly convoluted, as we have to ensure the same proportions of threatened/alien species within datasets:
  print("AgderMuseum_VascPlants")     # Print which dataset is being simulated
  df[df$dataset_Name=="AgderMuseum_VascPlants",]$list[sample(1:13373, 13373, FALSE)] <- rep(c("threat","alien",NA), c(93,365,12915))  # Shuffle the rows randomly and assign list status  
  print("Biofokus") 
  df[df$dataset_Name=="Biofokus",]$list[sample(1:282252, 282252, FALSE)] <- rep(c("threat","alien",NA), c(42577,9317,230358))  
  print("eBird")  
  df[df$dataset_Name=="eBird",]$list[sample(1:77373, 77373, FALSE)] <- rep(c("threat","alien",NA),c(14086,369,62918)) 
  print("J.B.Jordal")  
  df[df$dataset_Name=="J.B.Jordal",]$list[sample(1:111452, 111452, FALSE)] <- rep(c("threat","alien",NA), c(5368,2463,103621))
  print("NBIC_CitizenScience") 
  df[df$dataset_Name=="NBIC_CitizenScience",]$list[sample(1:4465574, 4465574, FALSE)] <- rep(c("threat","alien",NA), c(419880,190406,3855288))
  print("NBIC_other") 
  df[df$dataset_Name=="NBIC_other",]$list[sample(1:530913, 530913, FALSE)] <- rep(c("threat","alien",NA),c(54952,3146,472815))   
  print("NTNU_INH_VascPlants") 
  df[df$dataset_Name=="NTNU_INH_VascPlants",]$list[sample(1:8952, 8952, FALSE)] <- rep(c("threat","alien",NA), c(103,87,8762)) 
  print("UiO_Lichen")  
  df[df$dataset_Name=="UiO_Lichen",]$list[sample(1:5083, 5083, FALSE)] <- rep(c("threat","alien",NA),c(272,0,4811))   
  print("UiO_VascPlantsNotes")  
  df[df$dataset_Name=="UiO_VascPlantsNotes",]$list[sample(1:44760, 44760, FALSE)] <- rep(c("threat","alien",NA),c(670,1613,42477))   
  print("UiO_VascPlantsObs") 
  df[df$dataset_Name=="UiO_VascPlantsObs",]$list[sample(1:82528, 82528, FALSE)] <- rep(c("threat","alien",NA),c(10063,36244,36221))  
  
  # Summarize in a single table
  sum_area_x <- sum_area[,c("artype","area","cvr","area.prop")]    # Retrieve the constant areas of the map
  sum_area_x$nrec <- as.data.frame(table(df$artype_v2))$Freq   # Calculate no. records in each land cover
  sum_area_x$nrec_AgderMuseum_VascPlants <- as.data.frame(table(df[df$dataset_Name=="AgderMuseum_VascPlants",]$artype_v2))$Freq  
  sum_area_x$nrec_Biofokus <- as.data.frame(table(df[df$dataset_Name=="Biofokus",]$artype_v2))$Freq  
  sum_area_x$nrec_eBird <- as.data.frame(table(df[df$dataset_Name=="eBird",]$artype_v2))$Freq  
  sum_area_x$nrec_J.B.Jordal <- as.data.frame(table(df[df$dataset_Name=="J.B.Jordal",]$artype_v2))$Freq  
  sum_area_x$nrec_NBIC_CitizenScience <- as.data.frame(table(df[df$dataset_Name=="NBIC_CitizenScience",]$artype_v2))$Freq  
  sum_area_x$nrec_NBIC_other <- as.data.frame(table(df[df$dataset_Name=="NBIC_other",]$artype_v2))$Freq  
  sum_area_x$nrec_NTNU_INH_VascPlants <- as.data.frame(table(df[df$dataset_Name=="NTNU_INH_VascPlants",]$artype_v2))$Freq  
  sum_area_x$nrec_UiO_Lichen <- as.data.frame(table(df[df$dataset_Name=="UiO_Lichen",]$artype_v2))$Freq  
  sum_area_x$nrec_UiO_VascPlantsNotes <- as.data.frame(table(df[df$dataset_Name=="UiO_VascPlantsNotes",]$artype_v2))$Freq  
  sum_area_x$nrec_UiO_VascPlantsObs <- as.data.frame(table(df[df$dataset_Name=="UiO_VascPlantsObs",]$artype_v2))$Freq  
  sum_area_x$nrec_threat <- as.data.frame(table(df[df$list=="threat",]$artype_v2))$Freq  # No. redlisted records
  sum_area_x$nrec_alien <- as.data.frame(table(df[df$list=="alien",]$artype_v2))$Freq  # No. alien records
  # Create name vectors
  name_x <- paste('sum_rev/sum_area_',i,'.txt',sep='')  # Create a name for the export
  name_y <- paste('raw_rev/raw_numbers_',i,'.txt',sep='')  # Create a name for the export
  # Export the tables
  write.table(sum_area_x, file = name_x, row.names = TRUE, col.names = TRUE)  # Export as a textfile
  write.table(as.data.frame(table(df$dataset_Name, df$list, df$artype_v2)),
              file = name_y, row.names = TRUE, col.names = TRUE)  # Export raw numbers as a textfile
}

### To import the files exported above:
filelist <- list.files(path = "sum_rev", recursive = FALSE, pattern = "\\.txt$", full.names = TRUE)
datalist <- lapply(filelist, FUN=read.table, header=TRUE)
datafr <- do.call("rbind", datalist)                

##--- 4.3 Model the simulations against data ---####
##--- 4.3.1 Overall ---####
# Preliminary analyses proved Poisson distribution with identity link to give the best fit.
m_sim_poisid <- glm(nrec ~ area.prop, data=datafr, family=poisson(link="identity"))  

# Model validation  
par(mfrow = c(1,1), cex.lab = 1.5)
plot(x = fitted(m_sim_poisid), y = rstandard(m_sim_poisid), xlab = "Fitted values", ylab = "Residuals", main = "Homogeneity?") # Residuals vs. fittes values
abline(h = 0, lty = 2)
hist(rstandard(m_sim_poisid), main = "Normality", breaks=10)  # Normality of residuals
plot(x = datafr$area, y = rstandard(m_sim_poisid))  # Residuals vs. covariate
abline(h = 0, lty = 2)

ggplot(data=datafr, aes(x=area.prop, y=nrec)) +
  # Simulated data
  geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), linetype="dashed", size=.5) +  # Poisson model with identity link
  geom_point(aes(color=cvr), shape=4) +
  # Observed data
  geom_point(data=sum_area, aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
  geom_smooth(data=sum_area, aes(x=area.prop, y=nrec),
              method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity link 
  #coord_cartesian(xlim=c(0.0520, 0.05225), ylim = c(305460, 308405))   +  # For zooming in on a specified section
  #coord_cartesian(xlim=c(0.0521, 0.05215), ylim = c(283500, 284500))   +  # Even more zoomed in, just to see the CI
  theme_minimal(base_size = 10) +
  scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                 "green3","lightgreen","blue","cyan","gray","navy","white")) 

# Make a dataframe of the observed, predicted and CI-values for the model 
X <- model.matrix(~ area.prop, data = datafr[1:14,])    # Matrix of the needed values
rownames(X) <- datafr[1:14,"cvr"]

pred_overall_CI = predict(m_sim_poisid, newdata=as.data.frame(X), type='response', se.fit=TRUE) 
pred_overall <- data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec, predCI_all = pred_overall_CI$fit, lowerCI_all = pred_overall_CI$fit - 1.96 * pred_overall_CI$se.fit, upperCI_all = pred_overall_CI$fit + 1.96 * pred_overall_CI$se.fit)

##--- 4.3.2 Between datasets ---####
# Generalized Linear models of the simulated data - one for each of the datasets. In some cases, we have to supply starting values for the models to converge
# Gather the data, models and predictions in lists
datasets.list <- {list(data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_AgderMuseum_VascPlants),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_Biofokus),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_eBird),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_J.B.Jordal),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_NBIC_CitizenScience),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_NBIC_other),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_NTNU_INH_VascPlants),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_UiO_Lichen),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_UiO_VascPlantsNotes),
                       data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_UiO_VascPlantsObs))}
names(datasets.list) <- {c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", 
                           "NBIC_other", "NTNU_INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")}
sum_area_ls1 <- {list(data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_AgderMuseum_VascPlants),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_Biofokus),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_eBird),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_J.B.Jordal),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_NBIC_CitizenScience),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_NBIC_other),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_NTNU_INH_VascPlants),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_UiO_Lichen),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_UiO_VascPlantsNotes),
                      data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec_UiO_VascPlantsObs))}
names(sum_area_ls1) <- {c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other", "NTNU_INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")}

# Individual models for each dataset:
{
  # AgderMuseum_VascPlants
  {m_sim_1 <- glm(nrec ~ area.prop, data=datasets.list[["AgderMuseum_VascPlants"]], family=poisson(link="identity"), start=c(0.1,1000))
  summary(m_sim_1)
  step(m_sim_1)}
  # Biofokus
  {
    m_sim_2 <- glm(nrec ~ area.prop, data=datasets.list[["Biofokus"]], family=poisson(link="identity"))
    summary(m_sim_2)
    step(m_sim_2)
  }
  # eBird
  {
    m_sim_3 <- glm(nrec ~ area.prop, data=datasets.list[["eBird"]], family=poisson(link="identity"))
    summary(m_sim_3)
    step(m_sim_3)
  }
  # J.B.Jordal
  {
    m_sim_4 <- glm(nrec ~ area.prop, data=datasets.list[["J.B.Jordal"]], family=poisson(link="identity"))
    summary(m_sim_4)
    step(m_sim_4)
  }
  # NBIC_CitizenScience
  {
    m_sim_5 <- glm(nrec ~ area.prop, data=datasets.list[["NBIC_CitizenScience"]], family=poisson(link="identity"))
    summary(m_sim_5)
    step(m_sim_5)
  }
  # nrec_NBIC_other
  {
    m_sim_6 <- glm(nrec ~ area.prop, data=datasets.list[["NBIC_other"]], family=poisson(link="identity"))
    summary(m_sim_6)
    step(m_sim_6)
  }
  # NTNU_INH_VascPlants
  {
    m_sim_7 <- glm(nrec ~ area.prop, data=datasets.list[["NTNU_INH_VascPlants"]], family=poisson(link="identity"), start=c(0.1,1000))
    summary(m_sim_7)
    step(m_sim_7)
  }
  # UiO_Lichen
  {
    m_sim_8 <- glm(nrec ~ area.prop, data=datasets.list[["UiO_Lichen"]], family=poisson(link="identity"), start = c(0.1,1000))
    summary(m_sim_8)
    step(m_sim_8)
  }
  # UiO_VascPlantsNotes
  {
    m_sim_9 <- glm(nrec ~ area.prop, data=datasets.list[["UiO_VascPlantsNotes"]], family=poisson(link="identity"), start = c(0.1,1000))
    summary(m_sim_9)
    step(m_sim_9)
  }
  # UiO_VascPlantsObs
  {
    m_sim_10 <- glm(nrec ~ area.prop, data=datasets.list[["UiO_VascPlantsObs"]], family=poisson(link="identity"))
    summary(m_sim_10)
    step(m_sim_10)
  }
}
models_list <- list(m_sim_1,m_sim_2,m_sim_3,m_sim_4,m_sim_5,m_sim_6,m_sim_7,m_sim_8,m_sim_9,m_sim_10)
names(models_list) <- {c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience","NBIC_other", "NTNU_INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")}

# Plot of the models, simulated data and the original data -  Separate plots using grid.arrange
{
  p1 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_AgderMuseum_VascPlants, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,1000)), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_AgderMuseum_VascPlants, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p2 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_Biofokus, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_Biofokus, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p3 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_eBird, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_eBird, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p4 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_J.B.Jordal, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_J.B.Jordal, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p5 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_NBIC_CitizenScience, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_NBIC_CitizenScience, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p6 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_NBIC_other, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_NBIC_other, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p7 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_NTNU_INH_VascPlants, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,1000)), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_NTNU_INH_VascPlants, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p8 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_UiO_Lichen, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,1000)), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_UiO_Lichen, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p9 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_UiO_VascPlantsNotes, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,1000)), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_UiO_VascPlantsNotes, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p10 <-ggplot(data=datafr, aes(x=area.prop, y=nrec_UiO_VascPlantsObs, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_UiO_VascPlantsObs, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
}
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow=3, ncol=4)

# Predictions for each model/dataset 
X <- model.matrix(~ area.prop, data = datafr[1:14,])
rownames(X) <- datafr[1:14,"cvr"]

preds_dataset <- {list(predict(models_list[["AgderMuseum_VascPlants"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["Biofokus"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["eBird"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["J.B.Jordal"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["NBIC_CitizenScience"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["NBIC_other"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["NTNU_INH_VascPlants"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["UiO_Lichen"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["UiO_VascPlantsNotes"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
                       predict(models_list[["UiO_VascPlantsObs"]], newdata=as.data.frame(X), type='response', se.fit=TRUE))}
names(preds_dataset) <- {c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other", "NTNU_INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")}

model_pred_dataset <- {list(data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["AgderMuseum_VascPlants"]]$area.prop,
                                       nrec = sum_area_ls1[["AgderMuseum_VascPlants"]]$nrec,
                                       predCI = preds_dataset[["AgderMuseum_VascPlants"]]$fit,
                                       lowerCI = preds_dataset[["AgderMuseum_VascPlants"]]$fit - 1.96 * preds_dataset[["AgderMuseum_VascPlants"]]$se.fit,
                                       upperCI = preds_dataset[["AgderMuseum_VascPlants"]]$fit + 1.96 * preds_dataset[["AgderMuseum_VascPlants"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["Biofokus"]]$area.prop,
                                       nrec = sum_area_ls1[["Biofokus"]]$nrec,
                                       predCI = preds_dataset[["Biofokus"]]$fit,
                                       lowerCI = preds_dataset[["Biofokus"]]$fit - 1.96 * preds_dataset[["Biofokus"]]$se.fit,
                                       upperCI = preds_dataset[["Biofokus"]]$fit + 1.96 * preds_dataset[["Biofokus"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["eBird"]]$area.prop,
                                       nrec = sum_area_ls1[["eBird"]]$nrec,
                                       predCI = preds_dataset[["eBird"]]$fit,
                                       lowerCI = preds_dataset[["eBird"]]$fit - 1.96 * preds_dataset[["eBird"]]$se.fit,
                                       upperCI = preds_dataset[["eBird"]]$fit + 1.96 * preds_dataset[["eBird"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["J.B.Jordal"]]$area.prop,
                                       nrec = sum_area_ls1[["J.B.Jordal"]]$nrec,
                                       predCI = preds_dataset[["J.B.Jordal"]]$fit,
                                       lowerCI = preds_dataset[["J.B.Jordal"]]$fit - 1.96 * preds_dataset[["J.B.Jordal"]]$se.fit,
                                       upperCI = preds_dataset[["J.B.Jordal"]]$fit + 1.96 * preds_dataset[["J.B.Jordal"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["NBIC_CitizenScience"]]$area.prop,
                                       nrec = sum_area_ls1[["NBIC_CitizenScience"]]$nrec,
                                       predCI = preds_dataset[["NBIC_CitizenScience"]]$fit,
                                       lowerCI = preds_dataset[["NBIC_CitizenScience"]]$fit - 1.96 * preds_dataset[["NBIC_CitizenScience"]]$se.fit,
                                       upperCI = preds_dataset[["NBIC_CitizenScience"]]$fit + 1.96 * preds_dataset[["NBIC_CitizenScience"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["NBIC_other"]]$area.prop,
                                       nrec = sum_area_ls1[["NBIC_other"]]$nrec,
                                       predCI = preds_dataset[["NBIC_other"]]$fit,
                                       lowerCI = preds_dataset[["NBIC_other"]]$fit - 1.96 * preds_dataset[["NBIC_other"]]$se.fit,
                                       upperCI = preds_dataset[["NBIC_other"]]$fit + 1.96 * preds_dataset[["NBIC_other"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["NTNU_INH_VascPlants"]]$area.prop,
                                       nrec = sum_area_ls1[["NTNU_INH_VascPlants"]]$nrec,
                                       predCI = preds_dataset[["NTNU_INH_VascPlants"]]$fit,
                                       lowerCI = preds_dataset[["NTNU_INH_VascPlants"]]$fit - 1.96 * preds_dataset[["NTNU_INH_VascPlants"]]$se.fit,
                                       upperCI = preds_dataset[["NTNU_INH_VascPlants"]]$fit + 1.96 * preds_dataset[["NTNU_INH_VascPlants"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["UiO_Lichen"]]$area.prop,
                                       nrec = sum_area_ls1[["UiO_Lichen"]]$nrec,
                                       predCI = preds_dataset[["UiO_Lichen"]]$fit,
                                       lowerCI = preds_dataset[["UiO_Lichen"]]$fit - 1.96 * preds_dataset[["UiO_Lichen"]]$se.fit,
                                       upperCI = preds_dataset[["UiO_Lichen"]]$fit + 1.96 * preds_dataset[["UiO_Lichen"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["UiO_VascPlantsNotes"]]$area.prop,
                                       nrec = sum_area_ls1[["UiO_VascPlantsNotes"]]$nrec,
                                       predCI = preds_dataset[["UiO_VascPlantsNotes"]]$fit,
                                       lowerCI = preds_dataset[["UiO_VascPlantsNotes"]]$fit - 1.96 * preds_dataset[["UiO_VascPlantsNotes"]]$se.fit,
                                       upperCI = preds_dataset[["UiO_VascPlantsNotes"]]$fit + 1.96 * preds_dataset[["UiO_VascPlantsNotes"]]$se.fit),
                            data.frame(cvr = sum_area$cvr, 
                                       area.prop = sum_area_ls1[["UiO_VascPlantsObs"]]$area.prop,
                                       nrec = sum_area_ls1[["UiO_VascPlantsObs"]]$nrec,
                                       predCI = preds_dataset[["UiO_VascPlantsObs"]]$fit,
                                       lowerCI = preds_dataset[["UiO_VascPlantsObs"]]$fit - 1.96 * preds_dataset[["UiO_VascPlantsObs"]]$se.fit,
                                       upperCI = preds_dataset[["UiO_VascPlantsObs"]]$fit + 1.96 * preds_dataset[["UiO_VascPlantsObs"]]$se.fit))}
names(model_pred_dataset) <- {c("AgderMuseum_VascPlants", "Biofokus", "eBird", "J.B.Jordal", "NBIC_CitizenScience", "NBIC_other", "NTNU_INH_VascPlants", "UiO_Lichen", "UiO_VascPlantsNotes", "UiO_VascPlantsObs")}

# Plots
{
## AgderMuseum_VascPlants
{ggplot(data=datasets.list[["AgderMuseum_VascPlants"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["AgderMuseum_VascPlants"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["AgderMuseum_VascPlants"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["AgderMuseum_VascPlants"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## Biofokus
{ggplot(data=datasets.list[["Biofokus"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["Biofokus"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["Biofokus"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["Biofokus"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## eBird
{ggplot(data=datasets.list[["eBird"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["eBird"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["eBird"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["eBird"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## J.B.Jordal
{ggplot(data=datasets.list[["J.B.Jordal"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["J.B.Jordal"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["J.B.Jordal"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["J.B.Jordal"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## NBIC_CitizenScience
{ggplot(data=datasets.list[["NBIC_CitizenScience"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["NBIC_CitizenScience"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["NBIC_CitizenScience"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["NBIC_CitizenScience"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## NBIC_other
{ggplot(data=datasets.list[["NBIC_other"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["NBIC_other"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["NBIC_other"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["NBIC_other"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## NTNU_INH_VascPlants
{ggplot(data=datasets.list[["NTNU_INH_VascPlants"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["NTNU_INH_VascPlants"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,10000)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["NTNU_INH_VascPlants"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["NTNU_INH_VascPlants"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## UiO_Lichen
{ggplot(data=datasets.list[["UiO_Lichen"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["UiO_Lichen"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["UiO_Lichen"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["UiO_Lichen"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## UiO_VascPlantsNotes
{ggplot(data=datasets.list[["UiO_VascPlantsNotes"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["UiO_VascPlantsNotes"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["UiO_VascPlantsNotes"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["UiO_VascPlantsNotes"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
## UiO_VascPlantsObs
{ggplot(data=datasets.list[["UiO_VascPlantsObs"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset[["UiO_VascPlantsObs"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls1[["UiO_VascPlantsObs"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls1[["UiO_VascPlantsObs"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen",
                                                   "green3","lightgreen","blue","cyan","gray","navy","white")) }
}

### Comparison of residuals
# Convert the predictions to a long dataframe
model_pred_dataset.2 <- {data.frame(cvr = factor(c(rep(as.character(model_pred_dataset[[1]]$cvr),10)),
                                                 levels=c("Not.mapped", "Agriculture", "Developed", "Agri.grazing", "Snow.ice", "Forest.mix", "Agri.cultivated", "Freshwater", "Forest", "Mire","Forest.decid", "Ocean", "Forest.conif", "Firmground"), ordered = T),   # Ordered by area
                                    area.prop = c(rep(model_pred_dataset[[1]]$area.prop,10)),
                                    dataset = c(rep("AgderMuseum_VascPlants",14), rep("Biofokus",14), rep("eBird",14), rep("J.B.Jordal",14), rep("NBIC_CitizenScience",14), rep("NBIC_other",14),rep("NTNU_INH_VascPlants",14), rep("UiO_Lichen",14), rep("UiO_VascPlantsNotes",14),rep("UiO_VascPlantsObs",14)),
                                    nrec = c(model_pred_dataset[["AgderMuseum_VascPlants"]]$nrec, model_pred_dataset[["Biofokus"]]$nrec, model_pred_dataset[["eBird"]]$nrec,model_pred_dataset[["J.B.Jordal"]]$nrec, model_pred_dataset[["NBIC_CitizenScience"]]$nrec, model_pred_dataset[["NBIC_other"]]$nrec,model_pred_dataset[["NTNU_INH_VascPlants"]]$nrec, model_pred_dataset[["UiO_Lichen"]]$nrec, model_pred_dataset[["UiO_VascPlantsNotes"]]$nrec, model_pred_dataset[["UiO_VascPlantsObs"]]$nrec),
                                    pred = c(model_pred_dataset[["AgderMuseum_VascPlants"]]$predCI, model_pred_dataset[["Biofokus"]]$predCI, model_pred_dataset[["eBird"]]$predCI, model_pred_dataset[["J.B.Jordal"]]$predCI, model_pred_dataset[["NBIC_CitizenScience"]]$predCI, model_pred_dataset[["NBIC_other"]]$predCI, model_pred_dataset[["NTNU_INH_VascPlants"]]$predCI, model_pred_dataset[["UiO_Lichen"]]$predCI, model_pred_dataset[["UiO_VascPlantsNotes"]]$predCI, model_pred_dataset[["UiO_VascPlantsObs"]]$predCI),
                                    lwrCI =c(model_pred_dataset[["AgderMuseum_VascPlants"]]$lowerCI, model_pred_dataset[["Biofokus"]]$lowerCI, model_pred_dataset[["eBird"]]$lowerCI, model_pred_dataset[["J.B.Jordal"]]$lowerCI, model_pred_dataset[["NBIC_CitizenScience"]]$lowerCI, model_pred_dataset[["NBIC_other"]]$lowerCI, model_pred_dataset[["NTNU_INH_VascPlants"]]$lowerCI, model_pred_dataset[["UiO_Lichen"]]$lowerCI, model_pred_dataset[["UiO_VascPlantsNotes"]]$lowerCI, model_pred_dataset[["UiO_VascPlantsObs"]]$lowerCI),
                                    uprCI = c(model_pred_dataset[["AgderMuseum_VascPlants"]]$upperCI, model_pred_dataset[["Biofokus"]]$upperCI, model_pred_dataset[["eBird"]]$upperCI, model_pred_dataset[["J.B.Jordal"]]$upperCI, model_pred_dataset[["NBIC_CitizenScience"]]$upperCI, model_pred_dataset[["NBIC_other"]]$upperCI, model_pred_dataset[["NTNU_INH_VascPlants"]]$upperCI, model_pred_dataset[["UiO_Lichen"]]$upperCI, model_pred_dataset[["UiO_VascPlantsNotes"]]$upperCI, model_pred_dataset[["UiO_VascPlantsObs"]]$upperCI)) }

# Calculate absolute and relative residuals
model_pred_dataset.2$abs.resid <- model_pred_dataset.2$nrec - model_pred_dataset.2$pred
model_pred_dataset.2$rel.resid_final <- NA
for(i in 1:nrow(model_pred_dataset.2)){
  model_pred_dataset.2$rel.resid_final[i] <- model_pred_dataset.2$abs.resid[i] / mean(c(model_pred_dataset.2$nrec[i], model_pred_dataset.2$pred[i])) 
}

# Absolute residuals
ggplot(data=model_pred_dataset.2, aes(x=cvr, y=abs.resid, color=dataset, shape=dataset)) +
  geom_point(alpha=.75, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")) +
  scale_shape_manual(values = c(15,16,17,18,3,4,5,6,7,8)) +
  xlab("Land-cover") +
  ylab("Absolute residual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom", legend.box = "horizontal") 

# Relative residuals
ggplot(data=model_pred_dataset.2, aes(x=cvr, y=rel.resid_final, color=dataset, shape=dataset)) +
  geom_point(alpha=.75, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")) +
  scale_shape_manual(values = c(15,16,17,18,3,4,5,6,7,8)) +
  xlab("Land-cover") +
  ylab("Relative residual (v2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom", legend.box = "horizontal")


##--- 4.3.3 Between lists ---####
# Make linear models of the simulated data - one for each of the list-types
datasets.list2 <- {list(data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_threat), data.frame(cvr = datafr$cvr, area.prop = datafr$area.prop, nrec = datafr$nrec_alien))}
names(datasets.list2) <- {c("threat", "alien")}
sum_area_ls2 <- {list(data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec), data.frame(cvr = sum_area$cvr, area.prop = sum_area$area.prop, nrec = sum_area$nrec))}
names(sum_area_ls2) <- {c("threat", "alien")}

# Individual models for each list:
{
  # Redlisted
  {m_threat <- glm(nrec ~ area.prop, data=datasets.list2[["threat"]], family=poisson(link="identity"))
  summary(m_sim_1)
  step(m_sim_1)}
  # Alien
  {
    m_alien <- glm(nrec ~ area.prop, data=datasets.list2[["alien"]], family=poisson(link="identity"))
    summary(m_sim_2)
    step(m_sim_2)
  }
}
models_list2 <- list(m_threat,m_alien)
names(models_list2) <- {c("threat", "alien")}

# Plot  models, simulated data and the original data -  Separate plots using grid.arrange
{
  p_threat <-ggplot(data=datafr, aes(x=area.prop, y=nrec_threat, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_AgderMuseum_VascPlants, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen", "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
  p_alien <-ggplot(data=datafr, aes(x=area.prop, y=nrec_alien, color=cvr)) + 
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity")), se=TRUE, size=.5, color="blue", linetype="dashed") +   
    geom_point(aes(color=cvr), shape=2, alpha=.5, size=3) +  
    geom_point(data=sum_area, aes(x=area.prop, y=nrec_Biofokus, color=cvr), shape=17, size=3) +
    scale_color_manual(values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen", "green3","lightgreen","blue","cyan","gray","navy","white")) +
    theme(legend.position="none")
}
grid.arrange(p_threat,p_alien, nrow=1, ncol=2)

# Make predictions for each model/dataset:
X <- model.matrix(~ area.prop, data = datafr[1:14,])          # A matrix of the needed values - obs on using the NOT 'raw' dataframe
rownames(X) <- datafr[1:14,"cvr"]

preds_dataset2 <- {list(predict(models_list2[["threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE), predict(models_list2[["alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE))}
names(preds_dataset2) <- {c("threat", "alien")}

model_pred_dataset2 <- {list(data.frame(cvr = sum_area$cvr, area.prop = sum_area_ls2[["threat"]]$area.prop, nrec = sum_area_ls2[["threat"]]$nrec, predCI = preds_dataset2[["threat"]]$fit, lowerCI = preds_dataset2[["threat"]]$fit - 1.96 * preds_dataset2[["threat"]]$se.fit, upperCI = preds_dataset2[["threat"]]$fit + 1.96 * preds_dataset2[["threat"]]$se.fit),
                             data.frame(cvr = sum_area$cvr, area.prop = sum_area_ls2[["alien"]]$area.prop, nrec = sum_area_ls2[["alien"]]$nrec, predCI = preds_dataset2[["alien"]]$fit, lowerCI = preds_dataset2[["alien"]]$fit - 1.96 * preds_dataset2[["alien"]]$se.fit, upperCI = preds_dataset2[["alien"]]$fit + 1.96 * preds_dataset2[["alien"]]$se.fit))}
names(model_pred_dataset2) <- {c("threat", "alien")}

# Plots
## Redlisted
{ggplot(data=datasets.list2[["threat"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset2[["threat"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls2[["threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls2[["threat"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen", "green3","lightgreen","blue","cyan","gray","navy","white")) }
## Alien
{ggplot(data=datasets.list2[["alien"]], aes(x=area.prop, y=nrec)) +
    # CI ribbon
    geom_ribbon(data=model_pred_dataset2[["alien"]], aes(x=area.prop, ymin=lowerCI, ymax=upperCI), color="gray", alpha=0.5) +
    # Simulated data
    geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
    geom_point(aes(color=cvr), shape=4) +
    # Observed data
    geom_point(data=sum_area_ls2[["alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
    geom_smooth(data=sum_area_ls2[["alien"]], aes(x=area.prop, y=nrec),
                method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with identity-link link
    theme_minimal(base_size = 10) +
    #coord_cartesian(xlim=c(0.0517, 0.0525), ylim = c(282500, 285000))   +  # Even more zoomed in, just to see the PI
    scale_color_manual(name="Land-cover", values=c("darkorange","khaki1","orange","hotpink","tan3","darkgreen","forestgreen", "green3","lightgreen","blue","cyan","gray","navy","white")) }

### COMPARE RESIDUALS
# Convert the predictions to a long dataframe
model_pred_dataset.2.2 <- {data.frame(cvr = factor(c(rep(as.character(model_pred_dataset2[[1]]$cvr),2)),
                                                   levels=c("Not.mapped", "Agriculture", "Developed", "Agri.grazing", "Snow.ice", "Forest.mix", "Agri.cultivated", "Freshwater", "Forest", "Mire", "Forest.decid", "Ocean", "Forest.conif", "Firmground"), ordered = T),   # Ordered by area
                                      area.prop = c(rep(model_pred_dataset2[[1]]$area.prop,2)),
                                      dataset = c(rep("threat",14), rep("alien",14)),
                                      nrec = c(model_pred_dataset2[["threat"]]$nrec, model_pred_dataset2[["alien"]]$nrec),
                                      pred = c(model_pred_dataset2[["threat"]]$predCI, model_pred_dataset2[["alien"]]$predCI),
                                      lwrCI =c(model_pred_dataset2[["threat"]]$lowerCI, model_pred_dataset2[["alien"]]$lowerCI),
                                      uprCI = c(model_pred_dataset2[["threat"]]$upperCI, model_pred_dataset2[["alien"]]$upperCI)) }

# Calculate the absolute and the relative residuals
model_pred_dataset.2.2$abs.resid <- model_pred_dataset.2.2$nrec - model_pred_dataset.2.2$pred
model_pred_dataset.2.2$rel.resid_final <- NA
for(i in 1:nrow(model_pred_dataset.2.2)){
  model_pred_dataset.2.2$rel.resid_final[i] <- model_pred_dataset.2.2$abs.resid[i] / mean(c(model_pred_dataset.2.2$nrec[i], model_pred_dataset.2.2$pred[i])) 
}

# Absolute residuals
ggplot(data=model_pred_dataset.2.2, aes(x=cvr, y=abs.resid, color=dataset, shape=dataset)) +
  geom_point(alpha=.75, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("red","gray20")) +
  scale_shape_manual(values = c(15,19)) +
  xlab("Land-cover") +
  ylab("Absolute residual") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom", legend.box = "horizontal") 

# Relative residuals
ggplot(data=model_pred_dataset.2.2, aes(x=cvr, y=rel.resid_final, color=dataset, shape=dataset)) +
  geom_point(alpha=.75, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("red","gray20")) +
  scale_shape_manual(values = c(15,19)) +
  xlab("Land-cover") +
  ylab("Relative residual (v2)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom", legend.box = "horizontal")

##--- 4.4 Datasets AND conservation status simultaneously ---####
# List the files in the folder
filelist_raw <- list.files(path = "/home/ahomez/t/tanjakp/export/GBIF_bias/raw_rev", recursive = FALSE, pattern = "\\.txt$", full.names = TRUE)
datalist_raw <- lapply(filelist_raw, FUN=read.table, header=TRUE)
datafr_raw <- do.call("rbind", datalist_raw)                
colnames(datafr_raw) <- c("dataset","list","artype","nrec")
datafr_raw$artype <- as.factor(datafr_raw$artype)

sum_area2 <- sum_area[,c("artype","cvr","area","area.prop")]
sum_area2 <- full_join(sum_area2,  as.data.frame(table(GBIF_top10_2004$list, GBIF_top10_2004$cvr, GBIF_top10_2004$dataset_Name)), by = c("artype" = "Var2"))
colnames(sum_area2) <- c("artype","cvr","area","area.prop","list","dataset","nrec")
sum_area2$dataset <- revalue(sum_area2$dataset, c("NTNU-INH_VascPlants"="NTNU_INH_VascPlants")) # Correct a mistake in the datasetNames

datafr_raw <- full_join(datafr_raw, sum_area2[,c("artype","cvr","area","area.prop","list","dataset")])

##--- 4.4.1 Models, subdivisions ---####
# Subset into a list of smaller dataframes, divided by  dataset and conservation status
datals <- {list((subset(datafr_raw, dataset=="AgderMuseum_VascPlants" & list=="threat")), (subset(datafr_raw, dataset=="AgderMuseum_VascPlants" & list=="alien")),
                (subset(datafr_raw, dataset=="Biofokus" & list=="threat")), (subset(datafr_raw, dataset=="Biofokus" & list=="alien")),
                (subset(datafr_raw, dataset=="eBird" & list=="threat")), (subset(datafr_raw, dataset=="eBird" & list=="alien")),
                (subset(datafr_raw, dataset=="J.B.Jordal" & list=="threat")), (subset(datafr_raw, dataset=="J.B.Jordal" & list=="alien")),
                (subset(datafr_raw, dataset=="NBIC_CitizenScience" & list=="threat")), (subset(datafr_raw, dataset=="NBIC_CitizenScience" & list=="alien")),
                (subset(datafr_raw, dataset=="NBIC_other" & list=="threat")), (subset(datafr_raw, dataset=="NBIC_other" & list=="alien")),
                (subset(datafr_raw, dataset=="NTNU_INH_VascPlants" & list=="threat")), (subset(datafr_raw, dataset=="NTNU_INH_VascPlants" & list=="alien")),
                (subset(datafr_raw, dataset=="UiO_Lichen" & list=="threat")), (subset(datafr_raw, dataset=="UiO_Lichen" & list=="alien")),
                (subset(datafr_raw, dataset=="UiO_VascPlantsNotes" & list=="threat")), (subset(datafr_raw, dataset=="UiO_VascPlantsNotes" & list=="alien")),
                (subset(datafr_raw, dataset=="UiO_VascPlantsObs" & list=="threat")), (subset(datafr_raw, dataset=="UiO_VascPlantsObs" & list=="alien")))}
names(datals) <- {c("AgderMuseum_VascPlants_threat", "AgderMuseum_VascPlants_alien", "Biofokus_threat", "Biofokus_alien", "eBird_threat", "eBird_alien", "J.B.Jordal_threat", "J.B.Jordal_alien", "NBIC_CitizenScience_threat", "NBIC_CitizenScience_alien","NBIC_other_threat", "NBIC_other_alien", "NTNU_INH_VascPlants_threat", "NTNU_INH_VascPlants_alien", "UiO_Lichen_threat","UiO_Lichen_alien", "UiO_VascPlantsNotes_threat", "UiO_VascPlantsNotes_alien", "UiO_VascPlantsObs_threat", "UiO_VascPlantsObs_alien")}
sum_area_ls <- {list(subset(sum_area2, dataset=="AgderMuseum_VascPlants" & list=="threat"), subset(sum_area2, dataset=="AgderMuseum_VascPlants" & list=="alien"),
                     subset(sum_area2, dataset=="Biofokus" & list=="threat"), subset(sum_area2, dataset=="Biofokus" & list=="alien"),
                     subset(sum_area2, dataset=="eBird" & list=="threat"), subset(sum_area2, dataset=="eBird" & list=="alien"),
                     subset(sum_area2, dataset=="J.B.Jordal" & list=="threat"), subset(sum_area2, dataset=="J.B.Jordal" & list=="alien"),
                     subset(sum_area2, dataset=="NBIC_CitizenScience" & list=="threat"), subset(sum_area2, dataset=="NBIC_CitizenScience" & list=="alien"),
                     subset(sum_area2, dataset=="NBIC_other" & list=="threat"), subset(sum_area2, dataset=="NBIC_other" & list=="alien"),
                     subset(sum_area2, dataset=="NTNU_INH_VascPlants" & list=="threat"), subset(sum_area2, dataset=="NTNU_INH_VascPlants" & list=="alien"),
                     subset(sum_area2, dataset=="UiO_Lichen" & list=="threat"), subset(sum_area2, dataset=="UiO_Lichen" & list=="alien"),
                     subset(sum_area2, dataset=="UiO_VascPlantsNotes" & list=="threat"), subset(sum_area2, dataset=="UiO_VascPlantsNotes" & list=="alien"),
                     subset(sum_area2, dataset=="UiO_VascPlantsObs" & list=="threat"), subset(sum_area2, dataset=="UiO_VascPlantsObs" & list=="alien"))}
names(sum_area_ls) <- {c("AgderMuseum_VascPlants_threat", "AgderMuseum_VascPlants_alien", "Biofokus_threat", "Biofokus_alien", "eBird_threat", "eBird_alien", "J.B.Jordal_threat", "J.B.Jordal_alien", "NBIC_CitizenScience_threat", "NBIC_CitizenScience_alien", "NBIC_other_threat", "NBIC_other_alien", "NTNU_INH_VascPlants_threat", "NTNU_INH_VascPlants_alien", "UiO_Lichen_threat","UiO_Lichen_alien", "UiO_VascPlantsNotes_threat", "UiO_VascPlantsNotes_alien", "UiO_VascPlantsObs_threat", "UiO_VascPlantsObs_alien")}

### Make the models and plot
# AgderMuseum_VascPlants
{
  m1.threat <- glm(nrec ~ area.prop, data=datals[["AgderMuseum_VascPlants_threat"]], family=poisson(link="identity"), start=c(0.0001,100))
  summary(m1.threat)
  {ggplot(data=datals[["AgderMuseum_VascPlants_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.0001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("AgderMuseum_VascPlants_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["AgderMuseum_VascPlants_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["AgderMuseum_VascPlants_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m1.alien <- glm(nrec ~ area.prop, data=datals[["AgderMuseum_VascPlants_alien"]], family=poisson(link="identity"), start=c(0.001,100))  # Make a Poisson model with an identity link
  summary(m1.alien)
  {ggplot(data=datals[["AgderMuseum_VascPlants_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("AgderMuseum_VascPlants_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["AgderMuseum_VascPlants_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["AgderMuseum_VascPlants_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# Biofokus
{
  m2.threat <- glm(nrec ~ area.prop, data=datals[["Biofokus_threat"]], family=poisson(link="identity"), start=c(0.1,10000))  # Make a Poisson model with an identity link
  summary(m2.threat)
  {ggplot(data=datals[["Biofokus_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,10000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("Biofokus_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["Biofokus_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["Biofokus_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m2.alien <- glm(nrec ~ area.prop, data=datals[["Biofokus_alien"]], family=poisson(link="identity"), start=c(0.01,100))  # Make a Poisson model with an identity link
  summary(m1.alien)
  {ggplot(data=datals[["Biofokus_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("Biofokus_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["Biofokus_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["Biofokus_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# eBird
{
  m3.threat <- glm(nrec ~ area.prop, data=datals[["eBird_threat"]], family=poisson(link="identity"), start=c(0.001,10000))  # Make a Poisson model with an identity link
  summary(m3.threat)
  {ggplot(data=datals[["eBird_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,10000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("eBird_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["eBird_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["eBird_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m3.alien <- glm(nrec ~ area.prop, data=datals[["eBird_alien"]], family=poisson(link="identity"), start=c(0.01,300))  # Make a Poisson model with an identity link
  summary(m3.alien)
  {ggplot(data=datals[["eBird_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,300)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("eBird_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["eBird_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["eBird_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# J.B.Jordal_threat
{
  m4.threat <- glm(nrec ~ area.prop, data=datals[["J.B.Jordal_threat"]], family=poisson(link="identity"), start=c(0.01,5000))  # Make a Poisson model with an identity link
  summary(m4.threat)
  {ggplot(data=datals[["J.B.Jordal_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,5000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("J.B.Jordal_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["J.B.Jordal_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["J.B.Jordal_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m4.alien <- glm(nrec ~ area.prop, data=datals[["J.B.Jordal_alien"]], family=poisson(link="identity"), start=c(0.001,300))  # Make a Poisson model with an identity link
  summary(m1.alien)
  {ggplot(data=datals[["J.B.Jordal_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,300)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("J.B.Jordal_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["J.B.Jordal_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["J.B.Jordal_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# NBIC_CitizenScience
{
  m5.threat <- glm(nrec ~ area.prop, data=datals[["NBIC_CitizenScience_threat"]], family=poisson(link="identity"), start=c(0.1,100000))  # Make a Poisson model with an identity link
  summary(m5.threat)
  {ggplot(data=datals[["NBIC_CitizenScience_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,100000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("NBIC_CitizenScience_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["NBIC_CitizenScience_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["NBIC_CitizenScience_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m5.alien <- glm(nrec ~ area.prop, data=datals[["NBIC_CitizenScience_alien"]], family=poisson(link="identity"), start=c(0.0001,100000))  # Make a Poisson model with an identity link
  summary(m5.alien)
  {ggplot(data=datals[["NBIC_CitizenScience_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.0001,100000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("NBIC_CitizenScience_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["NBIC_CitizenScience_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["NBIC_CitizenScience_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# NBIC_other
{
  m6.threat <- glm(nrec ~ area.prop, data=datals[["NBIC_other_threat"]], family=poisson(link="identity"), start=c(0.001,10000))  # Make a Poisson model with an identity link
  summary(m6.threat)
  {ggplot(data=datals[["NBIC_other_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,10000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("NBIC_other_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["NBIC_other_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["NBIC_other_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m6.alien <- glm(nrec ~ area.prop, data=datals[["NBIC_other_alien"]], family=poisson(link="identity"), start=c(0.01,1000))  # Make a Poisson model with an identity link
  summary(m6.alien)
  {ggplot(data=datals[["NBIC_other_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,1000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("NBIC_other_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["NBIC_other_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["NBIC_other_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# NTNU_INH_VascPlants
{
  m7.threat <- glm(nrec ~ area.prop, data=datals[["NTNU_INH_VascPlants_threat"]], family=poisson(link="identity"), start=c(-0.0001,100))  # Make a Poisson model with an identity link
  summary(m7.threat)
  {ggplot(data=datals[["NTNU_INH_VascPlants_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(-0.0001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("NTNU_INH_VascPlants_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["NTNU_INH_VascPlants_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["NTNU_INH_VascPlants_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m7.alien <- glm(nrec ~ area.prop, data=datals[["NTNU_INH_VascPlants_alien"]], family=poisson(link="identity"), start=c(0.0001,100))  # Make a Poisson model with an identity link
  summary(m7.alien)
  {ggplot(data=datals[["NTNU_INH_VascPlants_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.0001,100)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("NTNU_INH_VascPlants_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["NTNU_INH_VascPlants_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["NTNU_INH_VascPlants_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# UiO_Lichen
{
  m8.threat <- glm(nrec ~ area.prop, data=datals[["UiO_Lichen_threat"]], family=poisson(link="identity"), start=c(0.001,200))  # Make a Poisson model with an identity link
  summary(m8.threat)
  {ggplot(data=datals[["UiO_Lichen_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,200)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("UiO_Lichen_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["UiO_Lichen_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["UiO_Lichen_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  # m8.alien - as there are no observation of alien lichens, this i not relevant
}
# UiO_VascPlantsNotes
{
  m9.threat <- glm(nrec ~ area.prop, data=datals[["UiO_VascPlantsNotes_threat"]], family=poisson(link="identity"), start=c(0.001,1000))  # Make a Poisson model with an identity link
  summary(m9.threat)
  {ggplot(data=datals[["UiO_VascPlantsNotes_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,1000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("UiO_VascPlantsNotes_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["UiO_VascPlantsNotes_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["UiO_VascPlantsNotes_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m9.alien <- glm(nrec ~ area.prop, data=datals[["UiO_VascPlantsNotes_alien"]], family=poisson(link="identity"), start=c(0.001,1000))  # Make a Poisson model with an identity link
  summary(m9.alien)
  {ggplot(data=datals[["UiO_VascPlantsNotes_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.001,1000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("UiO_VascPlantsNotes_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["UiO_VascPlantsNotes_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["UiO_VascPlantsNotes_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
# UiO_VascPlantsObs
{
  m10.threat <- glm(nrec ~ area.prop, data=datals[["UiO_VascPlantsObs_threat"]], family=poisson(link="identity"), start=c(0.01,10000))  # Make a Poisson model with an identity link
  summary(m10.threat)
  {ggplot(data=datals[["UiO_VascPlantsObs_threat"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,10000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("UiO_VascPlantsObs_threat") +
      # Observed data
      geom_point(data=sum_area_ls[["UiO_VascPlantsObs_threat"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["UiO_VascPlantsObs_threat"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(0.1,100)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
  
  m10.alien <- glm(nrec ~ area.prop, data=datals[["UiO_VascPlantsObs_alien"]], family=poisson(link="identity"), start=c(0.01,10000))  # Make a Poisson model with an identity link
  summary(m10.alien)
  {ggplot(data=datals[["UiO_VascPlantsObs_alien"]], aes(x=area.prop, y=nrec)) +
      # Simulated data
      geom_smooth(method="glm", method.args = list(family = poisson(link="identity"), start=c(0.01,10000)), linetype="dashed", size=.5) +  # Poisson model with identity link
      geom_point(aes(color=cvr), shape=4) +
      ylab("UiO_VascPlantsObs_alien") +
      # Observed data
      geom_point(data=sum_area_ls[["UiO_VascPlantsObs_alien"]], aes(x=area.prop, y=nrec, color=cvr), shape=19, size=2) +
      geom_smooth(data=sum_area_ls[["UiO_VascPlantsObs_alien"]], aes(x=area.prop, y=nrec),
                  method="glm", method.args = list(family = poisson(link="identity"), start=c(1,1)), linetype="dotted", size=.5, se=T) +  # Poisson model with log-link link
      scale_color_manual(name="Land-cover", values=c("hotpink","orange","darkorange","khaki1","darkgreen","forestgreen", "green3","lightgreen","tan3","cyan","white","blue", "navy", "gray")) }
}
models <- list(m1.threat,m1.alien,m2.threat,m2.alien,m3.threat,m3.alien,m4.threat,m4.alien,m5.threat,m5.alien,
               m6.threat,m6.alien,m7.threat,m7.alien,m8.threat,m9.threat,m9.alien,m10.threat,m10.alien)
names(models) <- c("AgderMuseum_VascPlants_threat", "AgderMuseum_VascPlants_alien", "Biofokus_threat", "Biofokus_alien", "eBird_threat", "eBird_alien", "J.B.Jordal_threat", "J.B.Jordal_alien", "NBIC_CitizenScience_threat", "NBIC_CitizenScience_alien", "NBIC_other_threat", "NBIC_other_alien", "NTNU_INH_VascPlants_threat", "NTNU_INH_VascPlants_alien", "UiO_Lichen_threat", "UiO_VascPlantsNotes_threat", "UiO_VascPlantsNotes_alien", "UiO_VascPlantsObs_threat", "UiO_VascPlantsObs_alien")

# Model validation - perform this chunk for all models 
{
  par(mfrow = c(2,2), mar=c(5,4,4,2), cex.lab = 1.5)
  plot(x = fitted(models[["UiO_VascPlantsObs_alien"]]), y = rstandard(models[["UiO_VascPlantsObs_alien"]]),   # Residuals vs. fittes values
       xlab = "Fitted values", ylab = "Residuals", main = "Homogeneity?")
  abline(h = 0, lty = 2)
  hist(rstandard(models[["UiO_VascPlantsObs_alien"]]), main = "Normality", breaks=10)  # Normality of residuals
  plot(x = datals[["UiO_VascPlantsObs_alien"]]$area.prop, y = rstandard(models[["UiO_VascPlantsObs_alien"]]),
       xlab = "Covariate", ylab = "Residuals", main = "Residuals vs. covariate")  # Residuals vs. covariate
  abline(h = 0, lty = 2) 
  plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes=FALSE)
  text(x=5,y=5,labels="UiO_VascPlantsObs_alien", cex=1.5)
}

# Get predicted values for all models:
X <- model.matrix(~ area.prop, data = datafr[1:14,]) 
rownames(X) <- datafr[1:14,"cvr"]

preds <- {list(predict(models[["AgderMuseum_VascPlants_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["AgderMuseum_VascPlants_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["Biofokus_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["Biofokus_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["eBird_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["eBird_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["J.B.Jordal_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["J.B.Jordal_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["NBIC_CitizenScience_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["NBIC_CitizenScience_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["NBIC_other_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["NBIC_other_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["NTNU_INH_VascPlants_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["NTNU_INH_VascPlants_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["UiO_Lichen_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["UiO_VascPlantsNotes_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["UiO_VascPlantsNotes_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["UiO_VascPlantsObs_threat"]], newdata=as.data.frame(X), type='response', se.fit=TRUE),
               predict(models[["UiO_VascPlantsObs_alien"]], newdata=as.data.frame(X), type='response', se.fit=TRUE))}
names(preds) <- {c("AgderMuseum_VascPlants_threat", "AgderMuseum_VascPlants_alien", "Biofokus_threat", "Biofokus_alien","eBird_threat", "eBird_alien", "J.B.Jordal_threat", "J.B.Jordal_alien", "NBIC_CitizenScience_threat", "NBIC_CitizenScience_alien","NBIC_other_threat", "NBIC_other_alien", "NTNU_INH_VascPlants_threat", "NTNU_INH_VascPlants_alien", "UiO_Lichen_threat", "UiO_VascPlantsNotes_threat", "UiO_VascPlantsNotes_alien", "UiO_VascPlantsObs_threat", "UiO_VascPlantsObs_alien")}

model_pred <- {list(data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["AgderMuseum_VascPlants_threat"]]$area.prop,
                               nrec = sum_area_ls[["AgderMuseum_VascPlants_threat"]]$nrec,
                               predCI = preds[["AgderMuseum_VascPlants_threat"]]$fit,
                               lowerCI = preds[["AgderMuseum_VascPlants_threat"]]$fit - 1.96 * preds[["AgderMuseum_VascPlants_threat"]]$se.fit,
                               upperCI = preds[["AgderMuseum_VascPlants_threat"]]$fit + 1.96 * preds[["AgderMuseum_VascPlants_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["AgderMuseum_VascPlants_alien"]]$area.prop,
                               nrec = sum_area_ls[["AgderMuseum_VascPlants_alien"]]$nrec,
                               predCI = preds[["AgderMuseum_VascPlants_alien"]]$fit,
                               lowerCI = preds[["AgderMuseum_VascPlants_alien"]]$fit - 1.96 * preds[["AgderMuseum_VascPlants_alien"]]$se.fit,
                               upperCI = preds[["AgderMuseum_VascPlants_alien"]]$fit + 1.96 * preds[["AgderMuseum_VascPlants_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["Biofokus_threat"]]$area.prop,
                               nrec = sum_area_ls[["Biofokus_threat"]]$nrec,
                               predCI = preds[["Biofokus_threat"]]$fit,
                               lowerCI = preds[["Biofokus_threat"]]$fit - 1.96 * preds[["Biofokus_threat"]]$se.fit,
                               upperCI = preds[["Biofokus_threat"]]$fit + 1.96 * preds[["Biofokus_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["Biofokus_alien"]]$area.prop,
                               nrec = sum_area_ls[["Biofokus_alien"]]$nrec,
                               predCI = preds[["Biofokus_alien"]]$fit,
                               lowerCI = preds[["Biofokus_alien"]]$fit - 1.96 * preds[["Biofokus_alien"]]$se.fit,
                               upperCI = preds[["Biofokus_alien"]]$fit + 1.96 * preds[["Biofokus_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["eBird_threat"]]$area.prop,
                               nrec = sum_area_ls[["eBird_threat"]]$nrec,
                               predCI = preds[["eBird_threat"]]$fit,
                               lowerCI = preds[["eBird_threat"]]$fit - 1.96 * preds[["eBird_threat"]]$se.fit,
                               upperCI = preds[["eBird_threat"]]$fit + 1.96 * preds[["eBird_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["eBird_alien"]]$area.prop,
                               nrec = sum_area_ls[["eBird_alien"]]$nrec,
                               predCI = preds[["eBird_alien"]]$fit,
                               lowerCI = preds[["eBird_alien"]]$fit - 1.96 * preds[["eBird_alien"]]$se.fit,
                               upperCI = preds[["eBird_alien"]]$fit + 1.96 * preds[["eBird_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["J.B.Jordal_threat"]]$area.prop,
                               nrec = sum_area_ls[["J.B.Jordal_threat"]]$nrec,
                               predCI = preds[["J.B.Jordal_threat"]]$fit,
                               lowerCI = preds[["J.B.Jordal_threat"]]$fit - 1.96 * preds[["J.B.Jordal_threat"]]$se.fit,
                               upperCI = preds[["J.B.Jordal_threat"]]$fit + 1.96 * preds[["J.B.Jordal_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["J.B.Jordal_alien"]]$area.prop,
                               nrec = sum_area_ls[["J.B.Jordal_alien"]]$nrec,
                               predCI = preds[["J.B.Jordal_alien"]]$fit,
                               lowerCI = preds[["J.B.Jordal_alien"]]$fit - 1.96 * preds[["J.B.Jordal_alien"]]$se.fit,
                               upperCI = preds[["J.B.Jordal_alien"]]$fit + 1.96 * preds[["J.B.Jordal_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["NBIC_CitizenScience_threat"]]$area.prop,
                               nrec = sum_area_ls[["NBIC_CitizenScience_threat"]]$nrec,
                               predCI = preds[["NBIC_CitizenScience_threat"]]$fit,
                               lowerCI = preds[["NBIC_CitizenScience_threat"]]$fit - 1.96 * preds[["NBIC_CitizenScience_threat"]]$se.fit,
                               upperCI = preds[["NBIC_CitizenScience_threat"]]$fit + 1.96 * preds[["NBIC_CitizenScience_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["NBIC_CitizenScience_alien"]]$area.prop,
                               nrec = sum_area_ls[["NBIC_CitizenScience_alien"]]$nrec,
                               predCI = preds[["NBIC_CitizenScience_alien"]]$fit,
                               lowerCI = preds[["NBIC_CitizenScience_alien"]]$fit - 1.96 * preds[["NBIC_CitizenScience_alien"]]$se.fit,
                               upperCI = preds[["NBIC_CitizenScience_alien"]]$fit + 1.96 * preds[["NBIC_CitizenScience_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["NBIC_other_threat"]]$area.prop,
                               nrec = sum_area_ls[["NBIC_other_threat"]]$nrec,
                               predCI = preds[["NBIC_other_threat"]]$fit,
                               lowerCI = preds[["NBIC_other_threat"]]$fit - 1.96 * preds[["NBIC_other_threat"]]$se.fit,
                               upperCI = preds[["NBIC_other_threat"]]$fit + 1.96 * preds[["NBIC_other_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["NBIC_other_alien"]]$area.prop,
                               nrec = sum_area_ls[["NBIC_other_alien"]]$nrec,
                               predCI = preds[["NBIC_other_alien"]]$fit,
                               lowerCI = preds[["NBIC_other_alien"]]$fit - 1.96 * preds[["NBIC_other_alien"]]$se.fit,
                               upperCI = preds[["NBIC_other_alien"]]$fit + 1.96 * preds[["NBIC_other_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["NTNU_INH_VascPlants_threat"]]$area.prop,
                               nrec = sum_area_ls[["NTNU_INH_VascPlants_threat"]]$nrec,
                               predCI = preds[["NTNU_INH_VascPlants_threat"]]$fit,
                               lowerCI = preds[["NTNU_INH_VascPlants_threat"]]$fit - 1.96 * preds[["NTNU_INH_VascPlants_threat"]]$se.fit,
                               upperCI = preds[["NTNU_INH_VascPlants_threat"]]$fit + 1.96 * preds[["NTNU_INH_VascPlants_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["NTNU_INH_VascPlants_alien"]]$area.prop,
                               nrec = sum_area_ls[["NTNU_INH_VascPlants_alien"]]$nrec,
                               predCI = preds[["NTNU_INH_VascPlants_alien"]]$fit,
                               lowerCI = preds[["NTNU_INH_VascPlants_alien"]]$fit - 1.96 * preds[["NTNU_INH_VascPlants_alien"]]$se.fit,
                               upperCI = preds[["NTNU_INH_VascPlants_alien"]]$fit + 1.96 * preds[["NTNU_INH_VascPlants_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["UiO_Lichen_threat"]]$area.prop,
                               nrec = sum_area_ls[["UiO_Lichen_threat"]]$nrec,
                               predCI = preds[["UiO_Lichen_threat"]]$fit,
                               lowerCI = preds[["UiO_Lichen_threat"]]$fit - 1.96 * preds[["UiO_Lichen_threat"]]$se.fit,
                               upperCI = preds[["UiO_Lichen_threat"]]$fit + 1.96 * preds[["UiO_Lichen_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["UiO_VascPlantsNotes_threat"]]$area.prop,
                               nrec = sum_area_ls[["UiO_VascPlantsNotes_threat"]]$nrec,
                               predCI = preds[["UiO_VascPlantsNotes_threat"]]$fit,
                               lowerCI = preds[["UiO_VascPlantsNotes_threat"]]$fit - 1.96 * preds[["UiO_VascPlantsNotes_threat"]]$se.fit,
                               upperCI = preds[["UiO_VascPlantsNotes_threat"]]$fit + 1.96 * preds[["UiO_VascPlantsNotes_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["UiO_VascPlantsNotes_alien"]]$area.prop,
                               nrec = sum_area_ls[["UiO_VascPlantsNotes_alien"]]$nrec,
                               predCI = preds[["UiO_VascPlantsNotes_alien"]]$fit,
                               lowerCI = preds[["UiO_VascPlantsNotes_alien"]]$fit - 1.96 * preds[["UiO_VascPlantsNotes_alien"]]$se.fit,
                               upperCI = preds[["UiO_VascPlantsNotes_alien"]]$fit + 1.96 * preds[["UiO_VascPlantsNotes_alien"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["UiO_VascPlantsObs_threat"]]$area.prop,
                               nrec = sum_area_ls[["UiO_VascPlantsObs_threat"]]$nrec,
                               predCI = preds[["UiO_VascPlantsObs_threat"]]$fit,
                               lowerCI = preds[["UiO_VascPlantsObs_threat"]]$fit - 1.96 * preds[["UiO_VascPlantsObs_threat"]]$se.fit,
                               upperCI = preds[["UiO_VascPlantsObs_threat"]]$fit + 1.96 * preds[["UiO_VascPlantsObs_threat"]]$se.fit),
                    data.frame(cvr = sum_area$cvr, 
                               area.prop = sum_area_ls[["UiO_VascPlantsObs_alien"]]$area.prop,
                               nrec = sum_area_ls[["UiO_VascPlantsObs_alien"]]$nrec,
                               predCI = preds[["UiO_VascPlantsObs_alien"]]$fit,
                               lowerCI = preds[["UiO_VascPlantsObs_alien"]]$fit - 1.96 * preds[["UiO_VascPlantsObs_alien"]]$se.fit,
                               upperCI = preds[["UiO_VascPlantsObs_alien"]]$fit + 1.96 * preds[["UiO_VascPlantsObs_alien"]]$se.fit))}
names(model_pred) <- {c("AgderMuseum_VascPlants_threat", "AgderMuseum_VascPlants_alien", "Biofokus_threat", "Biofokus_alien","eBird_threat", "eBird_alien", "J.B.Jordal_threat", "J.B.Jordal_alien", "NBIC_CitizenScience_threat", "NBIC_CitizenScience_alien","NBIC_other_threat", "NBIC_other_alien", "NTNU_INH_VascPlants_threat", "NTNU_INH_VascPlants_alien", "UiO_Lichen_threat","UiO_VascPlantsNotes_threat", "UiO_VascPlantsNotes_alien", "UiO_VascPlantsObs_threat", "UiO_VascPlantsObs_alien")}

# Make it a long dataframe for further testing and plotting 
model_pred.2 <- {data.frame(cvr = factor(c(rep(as.character(model_pred[[1]]$cvr),19)),
                                         levels=c("Not.mapped", "Agriculture", "Developed", "Agri.grazing", "Snow.ice", "Forest.mix", "Agri.cultivated", "Freshwater", "Forest", "Mire", "Forest.decid", "Ocean", "Forest.conif", "Firmground"), ordered = T),   # Ordered by area
                            area.prop = c(rep(model_pred[[1]]$area.prop,19)),
                            dataset = c(rep("AgderMuseum_VascPlants",28), rep("Biofokus",28), rep("eBird",28), rep("J.B.Jordal",28), rep("NBIC_CitizenScience",28), rep("NBIC_other",28), rep("NTNU_INH_VascPlants",28), rep("UiO_Lichen",14), rep("UiO_VascPlantsNotes",28),rep("UiO_VascPlantsObs",28)),
                            list = c(rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("threat",14), rep("alien",14), rep("threat",14), rep("alien",14)),
                            nrec = c(model_pred[["AgderMuseum_VascPlants_threat"]]$nrec, model_pred[["AgderMuseum_VascPlants_alien"]]$nrec, model_pred[["Biofokus_threat"]]$nrec, model_pred[["Biofokus_alien"]]$nrec, model_pred[["eBird_threat"]]$nrec, model_pred[["eBird_alien"]]$nrec, model_pred[["J.B.Jordal_threat"]]$nrec, model_pred[["J.B.Jordal_alien"]]$nrec, model_pred[["NBIC_CitizenScience_threat"]]$nrec, model_pred[["NBIC_CitizenScience_alien"]]$nrec, model_pred[["NBIC_other_threat"]]$nrec, model_pred[["NBIC_other_alien"]]$nrec, model_pred[["NTNU_INH_VascPlants_threat"]]$nrec, model_pred[["NTNU_INH_VascPlants_alien"]]$nrec, model_pred[["UiO_Lichen_threat"]]$nrec, model_pred[["UiO_VascPlantsNotes_threat"]]$nrec, model_pred[["UiO_VascPlantsNotes_alien"]]$nrec, model_pred[["UiO_VascPlantsObs_threat"]]$nrec, model_pred[["UiO_VascPlantsObs_alien"]]$nrec),
                            pred = c(model_pred[["AgderMuseum_VascPlants_threat"]]$predCI, model_pred[["AgderMuseum_VascPlants_alien"]]$predCI, model_pred[["Biofokus_threat"]]$predCI, model_pred[["Biofokus_alien"]]$predCI, model_pred[["eBird_threat"]]$predCI, model_pred[["eBird_alien"]]$predCI, model_pred[["J.B.Jordal_threat"]]$predCI, model_pred[["J.B.Jordal_alien"]]$predCI, model_pred[["NBIC_CitizenScience_threat"]]$predCI, model_pred[["NBIC_CitizenScience_alien"]]$predCI, model_pred[["NBIC_other_threat"]]$predCI, model_pred[["NBIC_other_alien"]]$predCI, model_pred[["NTNU_INH_VascPlants_threat"]]$predCI, model_pred[["NTNU_INH_VascPlants_alien"]]$predCI, model_pred[["UiO_Lichen_threat"]]$predCI, model_pred[["UiO_VascPlantsNotes_threat"]]$predCI, model_pred[["UiO_VascPlantsNotes_alien"]]$predCI, model_pred[["UiO_VascPlantsObs_threat"]]$predCI, model_pred[["UiO_VascPlantsObs_alien"]]$predCI),
                            lwrCI =c(model_pred[["AgderMuseum_VascPlants_threat"]]$lowerCI, model_pred[["AgderMuseum_VascPlants_alien"]]$lowerCI, model_pred[["Biofokus_threat"]]$lowerCI, model_pred[["Biofokus_alien"]]$lowerCI, model_pred[["eBird_threat"]]$lowerCI, model_pred[["eBird_alien"]]$lowerCI, model_pred[["J.B.Jordal_threat"]]$lowerCI, model_pred[["J.B.Jordal_alien"]]$lowerCI, model_pred[["NBIC_CitizenScience_threat"]]$lowerCI, model_pred[["NBIC_CitizenScience_alien"]]$lowerCI, model_pred[["NBIC_other_threat"]]$lowerCI, model_pred[["NBIC_other_alien"]]$lowerCI, model_pred[["NTNU_INH_VascPlants_threat"]]$lowerCI, model_pred[["NTNU_INH_VascPlants_alien"]]$lowerCI, model_pred[["UiO_Lichen_threat"]]$lowerCI, model_pred[["UiO_VascPlantsNotes_threat"]]$lowerCI, model_pred[["UiO_VascPlantsNotes_alien"]]$lowerCI, model_pred[["UiO_VascPlantsObs_threat"]]$lowerCI, model_pred[["UiO_VascPlantsObs_alien"]]$lowerCI),
                            uprCI = c(model_pred[["AgderMuseum_VascPlants_threat"]]$upperCI, model_pred[["AgderMuseum_VascPlants_alien"]]$upperCI, model_pred[["Biofokus_threat"]]$upperCI, model_pred[["Biofokus_alien"]]$upperCI, model_pred[["eBird_threat"]]$upperCI, model_pred[["eBird_alien"]]$upperCI, model_pred[["J.B.Jordal_threat"]]$upperCI, model_pred[["J.B.Jordal_alien"]]$upperCI, model_pred[["NBIC_CitizenScience_threat"]]$upperCI, model_pred[["NBIC_CitizenScience_alien"]]$upperCI, model_pred[["NBIC_other_threat"]]$upperCI, model_pred[["NBIC_other_alien"]]$upperCI, model_pred[["NTNU_INH_VascPlants_threat"]]$upperCI, model_pred[["NTNU_INH_VascPlants_alien"]]$upperCI, model_pred[["UiO_Lichen_threat"]]$upperCI, model_pred[["UiO_VascPlantsNotes_threat"]]$upperCI, model_pred[["UiO_VascPlantsNotes_alien"]]$upperCI, model_pred[["UiO_VascPlantsObs_threat"]]$upperCI, model_pred[["UiO_VascPlantsObs_alien"]]$upperCI))}

##--- 4.5 Comparisons of residuals ---####
# Calculate the absolute and the relative residuals
model_pred.2$abs.resid <- model_pred.2$nrec - model_pred.2$pred
model_pred.2$rel.resid_final <- NA
for(i in 1:nrow(model_pred.2)){
  model_pred.2$rel.resid_final[i] <- model_pred.2$abs.resid[i] / mean(c(model_pred.2$nrec[i], model_pred.2$pred[i])) 
}

# Plot the residuals
p_abs <- ggplot(data=model_pred.2, aes(x=cvr, y=abs.resid, color=list, shape=dataset)) +
  geom_point(alpha=.5, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("black","red")) +
  scale_shape_manual(values = c(15,16,17,18,3,4,5,6,7,8)) +
  xlab("Land-cover") +
  ylab("Absolute residual") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  #theme(legend.position="none") +
  facet_grid( list ~ .)

# Relative residuals 
p_rel2 <- ggplot(data=model_pred.2, aes(x=cvr, y=rel.resid_final, color=list, shape=dataset)) +
  geom_point(alpha=.5, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("black","red")) +
  scale_shape_manual(values = c(15,16,17,18,3,4,5,6,7,8)) +
  xlab("Land-cover") +
  ylab("Relative residual (v2)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  facet_grid( list ~ .)

grid.arrange(p_abs, p_rel2, nrow=2, ncol=1)

# Make the plots coloured by dominance of citizen science, museum records or otherwise profesional records:
model_pred.2$inst <- c(rep("mus",28),rep("prof",28),rep("cit",28),rep("prof",28),rep("cit",28),rep("prof",28),rep("mus",98))

ggplot(data=model_pred.2, aes(x=cvr, y=rel.resid_final, color=inst, shape=dataset)) +
  geom_point(alpha=.5, size=5) +
  geom_hline(yintercept = 0, linetype="dashed") +
  scale_colour_manual(values=c("darkorange","purple1","orchid1")) +
  scale_shape_manual(values = c(15,16,17,18,3,4,5,6,7,8)) +
  xlab("Land-cover") +
  ylab("Relative residual (v2)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  facet_grid( list ~ .)

##---------------------------------####
##--- 5. RAREFACTION CURVES ---####
# Prepare dataframes for rarefaction curves
raref.data <- {list("AgderMuseum_VascPlants, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="AgderMuseum_VascPlants", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "AgderMuseum_VascPlants, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="AgderMuseum_VascPlants", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "Biofokus, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="Biofokus", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "Biofokus, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="Biofokus", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "eBird, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="eBird", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "eBird, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="eBird", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "J.B.Jordal, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="J.B.Jordal", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "J.B.Jordal, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="J.B.Jordal", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "NBIC_CitizenScience, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="NBIC_CitizenScience", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "NBIC_CitizenScience, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="NBIC_CitizenScience", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "NBIC_other, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="NBIC_other", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "NBIC_other, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="NBIC_other", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "NTNU-INH_VascPlants, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="NTNU-INH_VascPlants", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "NTNU-INH_VascPlants, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="NTNU-INH_VascPlants", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "UiO_Lichen, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="UiO_Lichen", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "UiO_VascPlantsNotes, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="UiO_VascPlantsNotes", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "UiO_VascPlantsNotes, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="UiO_VascPlantsNotes", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "UiO_VascPlantsObs, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="UiO_VascPlantsObs", list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "UiO_VascPlantsObs, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(dataset_Name=="UiO_VascPlantsObs", list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "Combined, threatened" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(list=="threat") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec),
                    "Combined, alien" = st_drop_geometry(GBIF_top10) %>%
                      dplyr::select(species, dataset_Name, list) %>%
                      dplyr::filter(list=="alien") %>%
                      dplyr::group_by(species) %>%
                      dplyr::summarise(nrec = n()) %>%
                      dplyr::arrange(desc(nrec)) %>%
                      dplyr::pull(nrec))}
raref <- iNEXT(raref.data, datatype="abundance", endpoint = 2500000)  

# Plot the rarefaction curves
rarefcurve <- ggiNEXT(raref, type=1, color.var = "site", se=TRUE) +
  labs(x="Number of records") +
  scale_fill_manual(values=c("#a6cee3","#a6cee3","#1f78b4","#1f78b4","black","black","#b2df8a","#b2df8a","#33a02c","#33a02c","#fb9a99","#fb9a99","#e31a1c","#e31a1c","#fdbf6f","#fdbf6f","#ff7f00","#cab2d6","#cab2d6","#6a3d9a","#6a3d9a")) +
  scale_color_manual(values=c("gray20","#F8766D","gray20","#F8766D","gray20","#F8766D","gray20","#F8766D","gray20","#F8766D","gray20","#F8766D","gray20","#F8766D","gray20","#F8766D","#F8766D","gray20","#F8766D","gray20","#F8766D")) +
  scale_shape_manual(values=c(rep(19,21))) +
  theme_gray(base_size = 10) 

# Change  linewidth etc.
rarefcurve <- ggplot_build(rarefcurve)
rarefcurve$data[[2]]$size <- 0.75         # Lines are drawn in layer 2. Default is 1.5
rarefcurve$data[[1]]$size <- 1.5         # Points are drawn in layer 1. Default size is 5
rarefcurve$data[[1]]$shape <- 19       # Change the point shape to be the same for all
rarefcurve <- ggplot_gtable(rarefcurve)
grid.arrange(rarefcurve)