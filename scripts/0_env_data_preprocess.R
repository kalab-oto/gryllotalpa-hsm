library(raster)
library(rgdal)
library(envirem)

extent_eu <- extent(-10, 30, 45, 60)
processed_data <- file.path("..", "processed_data", "env")
dir.create(processed_data, showWarnings = F, recursive = T)

# CHELSA
chelsa_list <- list.files(file.path("..", "data", "env", "climatologies"), full.names = T, pattern = "CHELSA")

for (f in chelsa_list){
    chelsa <- raster(f)
    ch <- crop(chelsa, extent_eu)
    writeRaster(ch, file.path(processed_data, names(ch)), format = "GTiff", overwrite = T)
    print(paste0(names(ch), "-", "done"))
    }

ref_raster <- ch

## CHELSA rh_mean
rh_list <- list.files(processed_data, full.names = T, pattern = "_rh_")[1:12]
rh_mean <- mean(stack(rh_list))
names(rh_mean) <- "CHELSA_rh_mean"
writeRaster(rh_mean, file.path(processed_data, names(rh_mean)), format = "GTiff", overwrite = T)
unlink(rh_list)

## CHELSA ENVIREM
clim_path <- list.files(processed_data, full.names = T, pattern = "tmax|tmin|prec")

assignNames(
            tmax = "CHELSA_tmax10_##_1979.2013_V1",
            tmin = "CHELSA_tmin10_##_1979.2013_V1",
            precip = "CHELSA_prec_##_V1",
            tmean = "CHELSA_tmean10_##_1979.2013_V1"
            )

tmax <- stack(grep(clim_path, pattern = 'tmax', value = TRUE)[1:12]) / 10
tmin <- stack(grep(clim_path, pattern = 'tmin', value = TRUE)[1:12]) / 10

tmean <- stack(tmin)
for (i in 1:12){
    tmean[[i]] <- mean(tmax[[i]], tmin[[i]])
    }

names(tmean) <- names(tmin)
names(tmean) <- gsub("tmin", "tmean", names(tmean))

#srad 
prec <- stack(grep(clim_path, pattern =  'prec', value = TRUE)[-13])
srad_envirem <- ETsolradRasters(prec[[1]], 50)

#pet
pet <- monthlyPET(tmean, srad_envirem, tmax-tmin)


writeRaster(pet, filename = file.path(processed_data, paste0("CHELSA_ENVIREM_", names(pet))), bylayer = TRUE, format = "GTiff")

#aridity

arid <- aridityIndexThornthwaite(prec, pet)

writeRaster(arid, filename = file.path(processed_data, paste0("CHELSA_ENVIREM_", names(arid))), bylayer = TRUE, format = "GTiff", overwrite = T)

#moist

moist <- climaticMoistureIndex(sum(prec), sum(pet))

writeRaster(moist, filename = file.path(processed_data, paste0("CHELSA_ENVIREM_", names(moist))), bylayer = TRUE, format = "GTiff", overwrite = T)

unlink(clim_path)

# Soilgrids & Geomorpho90 (.vrts)

vrts <- list.files(file.path("..", "data", "env"), pattern = "*.vrt", full.names = T)

for (i in vrts){
    r_i <- crop(raster(i), ref_raster)
    r_i <- aggregate(r_i, 10)
    r_i <- resample(r_i, ref_raster)
    writeRaster(r_i, file.path(processed_data, names(r_i)), format = "GTiff", overwrite = T)
    print(paste0(names(r_i), "-", "done"))
}

# SoilTemp
SoilTemp_list <- list.files(file.path("..", "data", "env"), full.names = T, pattern = "SBIO")

r_SoilTemp <- resample(stack(SoilTemp_list), ref_raster)

writeRaster(r_SoilTemp, file.path(processed_data, names(r_SoilTemp)), format = "GTiff", overwrite = T, bylayer = T)

# WTD
wtd <- raster(file.path("..", "data", "env", "EURASIA_WTD_annualmean.nc"), varname = "WTD")

wtd_cl <- crop(wtd, ref_raster)
writeRaster(wtd_cl, file.path(processed_data, "EURASIA_WTD_annualmean.tif"))

# Mask
env_source <- list.files(file.path("..", "processed_data", "env"), full.names = T)

env_data <- stack(env_source)
na_mask <- sum(is.na(env_data)) > 0
na_mask[na_mask == 1] <- NA
env_data <- mask(env_data, na_mask)

writeRaster(env_data, file.path(processed_data, names(env_data)), format = "GTiff", overwrite = T, bylayer = T)
