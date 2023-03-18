`data_download/occurence_download.R`

 - download occurence data

`data_download/env_download.R`

 - download environmental raster data

`0_env_data_preprocess.R`

 - crop and resample downloaded data to study area extent and CHELSA dataset resolution

`1_occs_data_preprocess.R`

 - crop downloaded occurence data to study area extent, and filter by temporal and spatial resolution of CHELSA, and perfrom spatial thinning (`spThin`)

`2_models.R`

 - perform analysis (VIF, ExDet, ENMeval), run and evaluate model, export resutls

`3_plots.R`

 - make maps and plots



