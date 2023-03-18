library(rgbif)
library(spocc)
library(sp)
library(rgdal)
library(readxl)
library(rndop)#https://github.com/kalab-oto/rndop

species <- "Gryllotalpa gryllotalpa"

occ_dir <- file.path("..", "..", "data", "occ")
dir.create(occ_dir)

# NDOP
ndop_occs <- ndop_download(species)

ndop_path <- file.path("data", "occ", "gryll_ndop_tab.csv")
# ndop_occs <- read.csv(ndop_path, dec = ", ")

write.csv(ndop_occs, ndop_path, row.names = F)

# GBIF
dat_gbif <- occ_download_get("0209862-200613084148143", overwrite = T)
occs_gbif <- occ_download_import(dat_gbif)

write.csv(occs_gbif, file.path(occ_dir, "gryll_gbif.csv"), row.names = FALSE)

# iNat
sp_inat <- occ(query = species, from = 'inat', limit = 10000)
gdata <- sp_inat$inat$data$Gryllotalpa_gryllotalpa

gdata <- gdata[gdata$name == "Gryllotalpa gryllotalpa", ]

occs_inat <- cbind(as.data.frame(occs_inat), gdata$quality_grade, gdata$positional_accuracy, gdata$uri, gdata$id, gdata$observed_on, gdata$user.name)
names(occs_inat) <- c('lon', 'lat', "quality_grade", "positional_accuracy", "uri", "id", "observed_on", "user.name" )

write.csv(occs_inat, file.path(occ_dir, "gryll_inat.csv"), row.names = FALSE)

# Biomonitoring.sk 
# Due to unclear permissions for sharing the data we provide the GeoPackage `gryll_sk.gpkg` only with point geometry and url. The geometry was created manually all metadata and GPS positions can be checked at the biomonitoring.sk. IDs of used records are stored in `occs_meta.csv` file. The url are stored in GeoPackage or can be constructed with ID from csv file:

# http://www.biomonitoring.sk/OccurenceData/ZoologicalOccurenceRecord/Detail/<ID_OF_THE_RECORDS>?ReturnPage = oc_ZooGallery

# for example:
# http://www.biomonitoring.sk/OccurenceData/ZoologicalOccurenceRecord/Detail/1493970?ReturnPage = oc_ZooGallery


#' # Iorgu et al 2016 (Zookeys)
download.file(url = "http://zookeys.pensoft.net//lib/ajax_srv/article_elements_srv.php?action = download_suppl_file&instance_id = 31&article_id = 8804", destfile = "iorgu.xlsx", mode = "wb")
my_data <- as.data.frame(read_excel("iorgu.xlsx", skip = 1))
iorgu <- my_data[my_data$Species == species , ]

iorgu$Latitude <- gsub("°N", "", iorgu$Latitude)
iorgu$Longitude <- gsub("°E", "", iorgu$Longitude)
iorgu$Data <- substring(iorgu$Data, 1, 4)
names(iorgu)[8] <- "year"

write.csv(iorgu, file.path(occ_dir, "zookeys-605-073-s001.csv"))

unlink("iorgu.xlsx")

