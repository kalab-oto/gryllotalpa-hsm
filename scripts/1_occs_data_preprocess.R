library(rgdal)
library(spThin)
library(raster)

# Data
crs_5514 <- "+proj = krovak +lat_0 = 49.5 +lon_0 = 24.83333333333333 +alpha = 30.28813972222222 +k = 0.9999 +x_0 = 0 +y_0 = 0 +ellps = bessel +towgs84 = 589, 76, 480, 0, 0, 0, 0 +units = m +no_defs" 
crs_4326 <- "+proj = longlat +datum = WGS84 +no_defs"

occs_meta_path <- file.path("..", "data", "occ", "occs_meta.csv")
occs_meta <- read.csv(occs_meta_path)

inat_path <- file.path("..", "data", "occ", "gryll_inat.csv")
gbif_path <- file.path("..", "data", "occ", "gryll_gbif.csv")
iorgu_path <- file.path("..", "data", "occ", "zookeys-605-073-s001.csv")
biom_sk_path <- file.path("..", "data", "occ", "gryll_sk.gpkg")
cz_path <- file.path("..", "data", "occ", "gryll_cz.gpkg")
gbif_occs <- read.csv(gbif_path)
inat_occs <- read.csv(inat_path)
iorgu_occs <- read.csv(iorgu_path)

occs_names <- c("lat", "lon", "year", "accuracy", "source", "id")
occs <- setNames(data.frame(matrix(ncol = length(occs_names), nrow = 0)), occs_names)

## GBIF
gbif_occs$source <- "gbif"
gbif_occs <- gbif_occs[, c("decimalLatitude", "decimalLongitude", "year", "coordinateUncertaintyInMeters", "source", "occurrenceID")]
names(gbif_occs) <- occs_names

## iNat
inat_occs <- inat_occs[inat_occs$quality_grade == "research", ]
inat_occs$year <- as.numeric(substring(inat_occs$observed, 1, 4))
inat_occs$source <- "inat"
inat_occs <- inat_occs[, c("lat", "lon", "year", "positional_accuracy", "source", "id")]
inat_occs <- inat_occs[inat_occs$id %in% occs_meta[occs_meta$source == "inat", 2], ]
names(inat_occs) <- occs_names

## iorgu
iorgu_occs$source <-"iorgu2016"
iorgu_occs$accuracy <- 1000 #all GPS postions are at least for two digits, we assume that that its precise about 1 km
iorgu_occs$id <- NA
iorgu$year <- as.numeric(iorgu$year)
iorgu_occs <- iorgu_occs[, c("Latitude", "Longitude", "year", "accuracy", "source", "id")]
names(iorgu_occs) <- occs_names

## biomonitoring.sk
biom_sk_occs <- readOGR(biom_sk_path)
biom_sk_occs$lon <- biom_sk_occs@coords[, 1]
biom_sk_occs$lat <- biom_sk_occs@coords[, 2]
biom_sk_occs$source <- "biom_sk"

# due to insufficient or unclear permission conditions, accuracy and year are defined only for to meet the later condition !!!
biom_sk_occs$accuracy <- 1000
biom_sk_occs$year <- 2000

## Data CZ
cz_occs <- readOGR(cz_path)
cz_occs$source <- "cz"
cz_occs$id <- NA
cz_occs$lon <- cz_occs@coords[, 1]
cz_occs$lat <- cz_occs@coords[, 2]
# due to insufficient or unclear permission conditions, accuracy and year are defined only for to meet the later condition !!!
cz_occs$accuracy <- 1000
cz_occs$year <- 2000

occs_4326 <- rbind(inat_occs, gbif_occs, iorgu_occs)

occs_4326 <- occs_4326[!is.na(occs_4326$lat) & !is.na(occs_4326$lon), ]
occs_4326 <- SpatialPointsDataFrame(
                        coords = occs_4326[, c("lon", "lat")], 
                        data = occs_4326, 
                        proj4string = CRS(crs_4326)
                    )
                    
occs_4326 <- rbind(occs_4326, cz_occs, biom_sk_occs)

## NDOP
ndop_path <- file.path("..", "data", "occ", "gryll_ndop_tab.csv")
ndop_occs <- read.csv(ndop_path)
ndop_occs$source <- "ndop"
ndop_occs$year <- as.numeric(substring(ndop_occs$DATUM_OD, 1, 4))
ndop_occs <- ndop_occs[, c("X", "Y", "year", "CXPRESNOST", "source", "ID_NALEZ")]
ndop_occs <- ndop_occs[ndop_occs$ID_NALEZ %in% occs_meta[occs_meta$source == "ndop", 2], ]

ndop_occs <- SpatialPointsDataFrame(
                        coords = ndop_occs[c("X", "Y")], 
                        data = ndop_occs, 
                        proj4string = CRS(crs_5514)
                    )
ndop_occs <- spTransform(ndop_occs, CRS(crs_4326))
names(ndop_occs) <- occs_names
ndop_occs$lon <- ndop_occs@coords[, 1]
ndop_occs$lat <- ndop_occs@coords[, 2]

#############################################################
occs <- rbind(ndop_occs, occs_4326)

occs <- occs[!is.na(occs$year), ]
occs <- occs[!is.na(occs$accuracy), ]
occs <- occs[!is.na(occs$lat), ]

occs <- occs[occs$accuracy < 1001 & occs$year > 1978 & occs$year < 2014, ]

writeOGR(occs, file.path("..", "processed_data", paste0("occs.gpkg")), layer = "occs", driver = "GPKG", overwrite = T)

# Spatial filtering
# To reduce spatial bias in occurence data, we perform spatial thining, method, but at first we removed occurence points that coresponds to NA value at least in one of the environmental raster layers

na_mask <- raster(list.files(file.path("..", "processed_data", "env"), full.names = T)[1])
occs <- occs[which(extract(!is.na(na_mask), occs) == T), ]

occs$species <- "gryllotalpa"
rownames(occs@data) <- 1:nrow(occs)
occs$rid <- rownames(occs@data)

thinned_sp <- spThin::thin(loc.data =  as.data.frame(occs), 
            lat.col = "X", 
            long.col = "Y", 
            spec.col = "species", 
            thin.par = 10, #km
            reps = 1, 
            out.dir = file.path("..", "processed_data", "tmp"), 
            locs.thinned.list.return = TRUE, 
            write.files = FALSE
        )
thinned_sp <- as.data.frame(thinned_sp)
thinned_sp$rid <- rownames(thinned_sp)
thinned_sp <- merge(thinned_sp, occs, by = "rid")
thinned_sp <- SpatialPointsDataFrame(
                        coords = thinned_sp[c("X", "Y")], 
                        data = data.frame(thinned_sp[c("source", "species")]), 
                        proj4string = CRS(crs_4326)
                    )

writeOGR(thinned_sp, file.path("..", "processed_data", paste0("thinned_sp_10.gpkg")), layer = "thinned_sp", driver = "GPKG", overwrite = T)
