library(terra)
library(raster)
env_path <- file.path("..", "..", "data", "env")
dir.create(env_path, showWarnings = F)

if_download.file <- function (url, destfile, method){
    if (!file.exists(destfile)) {
    download.file(url, destfile, method = method) }
}

# CHELSA

##  monthly climatologies
dir.create(file.path(env_path, "climatologies"), showWarnings = F)
month <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
url_base <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/climatologies/"
url_base_ex <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/exchelsa/"

for (m in month){
        if_download.file(paste0(url_base, "tmin/CHELSA_tmin10_", m, "_1979-2013_V1.2_land.tif"), file.path(env_path, "climatologies", paste0("CHELSA_tmin10_", m, "_1979-2013_V1.2_land.tif")), "wget")
        if_download.file(paste0(url_base, "tmax/CHELSA_tmax10_", m, "_1979-2013_V1.2_land.tif"), file.path(env_path, "climatologies", paste0("CHELSA_tmax10_", m, "_1979-2013_V1.2_land.tif")), "wget")
        if_download.file(paste0(url_base, "prec/CHELSA_prec_", m, "_V1.2_land.tif"), file.path(env_path, "climatologies", paste0("CHELSA_prec_", m, "_V1.2_land.tif")), "wget")
}

for (m in 1:12){
        if_download.file(paste0(url_base_ex, "rh/CHELSA_rh_", m, "_1979-2013_V1.0.tif"), file.path(env_path, "climatologies", paste0("CHELSA_rh_", m, "_1979-2013_V1.0.tif")), "wget")
}

# Geomorpho90

geom_tiles <- c("n30e000", "n30e030", "n30w030", "n60e000", "n60e030", "n60w030")
gv <- c("cti")


dir.create(file.path(env_path, gv), showWarnings = F)
for (tl in geom_tiles){
    if_download.file(paste0("https://opentopography.s3.sdsc.edu/dataspace/OTDS.012020.4326.1/raster/", gv, "/", gv, "_90M_", tl, ".tar.gz"), file.path(env_path, gv, paste0(gv, "_90M_", tl, ".tar.gz")), "wget")
}

# ## Build `.vrt`
geom_vrt <- function (var){
    tif_l <- NULL
    l <- normalizePath(list.files(file.path(env_path, var), full.names = T))
    for (tr in l){
        tif_l <- c(tif_l, paste0("/vsitar/", file.path(tr, untar(tr, list = T))))
    }
    print(tif_l)
    terra::vrt(tif_l, file.path(env_path, paste0(var, ".vrt")), overwrite = T)
}

for (i in c("cti")){
    geom_vrt(i)
}

# wtd
if_download.file("http://thredds-gfnl.usc.es/thredds/fileServer/GLOBALWTDFTP/annualmeans/EURASIA_WTD_annualmean.nc", file.path(env_path, "EURASIA_WTD_annualmean.nc"))

# soilgrids
eur_grids <- c(66:70, 105:113, 149:156, 193:199, 235:242)
soils <- c("Gleysols", "Fluvisols")

for (l in soils){
    dir.create(file.path(env_path, l), showWarnings = F)
    for (i in eur_grids){
        url <- paste0("https://files.isric.org/soilgrids/latest/data/wrb/", l, "/", i, ".geotiff")
        try(if_download.file(url, file.path(env_path, l, paste0(i, ".geotiff"))))
    }
    terra::vrt(list.files(file.path(env_path, l), pattern = "*.geotiff", full.names = T), file.path(env_path, paste0(l, ".vrt")))
}

# SoilTemp
if_download.file("https://zenodo.org/record/4558732/files/SBIO6_Min_Temperature_of_Coldest_Month_5_15cm.tif?download = 1", file.path(env_path, "SBIO6_5_15.tif"))
if_download.file("https://zenodo.org/record/4558732/files/SBIO10_Mean_Temperature_of_Warmest_Quarter_5_15cm.tif?download = 1", file.path(env_path, "SBIO10_5_15.tif"))

# Admin boundaries
admin_path <- file.path("..", "..", "data", "admin")
dir.create(admin_path, showWarnings = F)
getData("GADM", country = "CZE", level = 0, path = admin_path)
