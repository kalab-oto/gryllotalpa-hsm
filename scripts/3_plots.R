library(raster)
library(RColorBrewer)
library(rgdal)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(sdm)
library(sf)
library(ggExtra)

export_dpi <- 100

cz_bound <- readRDS(file.path("..", "data", "admin", "gadm36_CZE_0_sp.rds"))
cz_bound <- spTransform(cz_bound, crs(st_crs(32633)$"proj4string"))

unc <- raster(file.path("..", "results", "unc_sd.grd"))
pred <- raster(file.path("..", "results", "ensem.grd"))
pred_tr <- raster(file.path("..", "results", "ensem_th.grd"))

pred_stack <- mask(projectRaster(stack(pred, pred_tr, pred_tr == 1, unc), crs = st_crs(32633)$"proj4string", method = "ngb"), cz_bound)
names(pred_stack) <- c("pred", "pred_tr", "pred_ens", "unc")

m <- readRDS(file.path("..", "results", "model.RDS"))

pdata <- merge(m@data@info@coords, m@data@features, "rID")
dim(pdata)
pdata[pdata[, 1] %in% m@data@species$sp_name@presence & pdata[, 1] %in% m@data@groups$training@indices$train, "bg_occ"] <- "occs_train"
pdata[pdata[, 1] %in% m@data@species$sp_name@presence & pdata[, 1] %in% m@data@groups$training@indices$test, "bg_occ"] <- "occs_test"
pdata[pdata[, 1] %in% m@data@species$sp_name@background & pdata[, 1] %in% m@data@groups$training@indices$test, "bg_occ"] <- "bg_test"
pdata[pdata[, 1] %in% m@data@species$sp_name@background & pdata[, 1] %in% m@data@groups$training@indices$train, "bg_occ"] <- "bg_train"
sf_pdata <- st_as_sf(pdata, coords = c("x", "y"), crs = st_crs(4326))

pred_df <- na.omit(as.data.frame(pred_stack, xy = T))

# Maps
labs <- c("predicted area", "testing occurence")

p_pr <- ggplot(pred_df) +
    geom_raster(aes(x = x, y = y, fill = pred)) +
    scale_fill_gradientn(colours = c("#f7fcf5", "#74c476", "#00441b"), breaks = c(min(pred_df$pred), 0.25, 0.5, 0.75, max(pred_df$pred)), 
    labels = c("0.00", 0.25, 0.5, 0.75, "1.00"), guide = guide_colourbar(direction = "horizontal", barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5, title.vjust = 2), name = "suitability") +
    geom_sf(data = st_transform(st_as_sf(cz_bound), st_crs(32633)), color = "grey30", fill = NA, show.legend = FALSE, lwd = 0.2) +
    theme_void() + theme(legend.position = c(0.8, 0.9))

p_prtr <- ggplot(pred_df) +
    geom_raster(aes(x = x, y = y, fill = pred_tr)) +
  scale_fill_gradientn(colours = c("#fee6ce", "#fdae6b", "#e6550d"), guide = guide_colourbar(direction = "horizontal", barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5, title.vjust = 2), name = "mean of binary predictions") +
  geom_sf(data = st_transform(st_as_sf(cz_bound), st_crs(32633)), color = "grey30", fill = NA, show.legend = FALSE, lwd = 0.2) +
    theme_void() + theme(legend.position = c(0.8, 0.9))

p_unc <- ggplot(pred_df) +
  geom_raster(aes(x = x, y = y, fill = unc)) +
  scale_fill_gradientn(
    colours = c("#fcfbfd", "#9e9ac8", "#3f007d"),
    breaks = c(min(pred_df$unc), 0.1, max(pred_df$unc)),
    labels = c(0, 0.1, round(max(na.omit(getValues(unc))), 2)),
    guide = guide_colourbar(direction = "horizontal", barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5, title.vjust = 2), 
    name = "standard deviation") +
  geom_sf(data = st_transform(st_as_sf(cz_bound), st_crs(32633)), color = "grey30", fill = NA, show.legend = FALSE, lwd = 0.2) +theme_void() + theme(legend.position = c(0.8, 0.9))

p_ens <- ggplot(pred_df[pred_df$pred_ens == 1, ]) +
    geom_tile(aes(x = x, y = y, fill = as.factor(pred_ens))) +
    geom_sf(data = st_transform(sf_pdata[sf_pdata$bg_occ == "occs_test", ], st_crs(32633)), aes(color = bg_occ)) +
    scale_color_manual(labels = labs[2], values = c("darkred")) +
    scale_fill_manual(labels = labs[1], values = c("grey")) +
    geom_sf(data = st_transform(st_as_sf(cz_bound), st_crs(32633)), color = "grey30", fill = NA, show.legend = FALSE, lwd = 0.2) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 2)) +
    theme_void() + theme(legend.position = c(0.8, 0.9), legend.spacing.y = unit(-0.2, 'cm')) + labs(fill = "", color = "")

p <- p_pr + p_prtr + p_unc + p_ens + plot_annotation(tag_levels = "A")

ggsave(
  file.path("..", "results", "plots", paste0("preds_", export_dpi, ".jpeg")), 
  plot = p,
  width = 16.1,
  height = 8.7,
  dpi = export_dpi,
)

# Bivariate map/plot
# R. Hijmans - https://stackoverflow.com/a/56269165/3984070
makeCM <- function(breaks = 10, upperleft, upperright, lowerleft, lowerright) {
    m <- matrix(ncol = breaks, nrow = breaks)
    b <- breaks - 1
    b <- (0:b) / b
    col1 <- rgb(colorRamp(c(upperleft, lowerleft))(b), max = 255)
    col2 <- rgb(colorRamp(c(upperright, lowerright))(b), max = 255)
    cm <- apply(cbind(col1, col2), 1, function(i) rgb(colorRamp(i)(b), max = 255))
    cm[, ncol(cm):1 ]
}

plotCM <- function(x, y, cm, xlab = "", ylab = "", main = "", axes = TRUE) {
    n <- cm
    n <- matrix(1:length(cm), nrow = nrow(cm), byrow = TRUE)
    r <- raster(n)
    cm <- cm[, ncol(cm):1]
    extent(r)[4] <- max(na.omit(getValues(y)))
    extent(r)[2] <- max(na.omit(getValues(x)))
    image(r, col = cm, axes = axes, xlab = xlab, ylab = ylab, main = main)
    return(r)
}

breaks <- 3
col_lu <- rgb(118, 42, 131, max = 255)
col_rd <- rgb(27, 120, 55, max = 255)
col_ru <- rgb(27, 42, 55, max = 255)
col_ld <- rgb(247, 247, 247, max = 255)

cmat <- makeCM(breaks, col_lu, col_ru, col_ld, col_rd)
r <- plotCM(pred_stack[["pred"]], pred_stack[["unc"]], cmat)

xy_df <- pred_df[, c("x", "y", "pred", "unc")]
xy_df[5]<-extract(r, xy_df[, 3:4])

t_cmat <- cmat[, ncol(cmat):1 ]
n_cmat <- data.frame(n = as.vector(matrix(1:length(t_cmat), nrow = nrow(t_cmat))), cmat = as.vector(t_cmat))
m_test <- merge(xy_df, n_cmat, by.x = "V5", by.y = "n", all.x = T)

pr_unc <- ggplot(m_test) +
  geom_raster(aes(x = x, y = y, fill = as.factor(V5)), show.legend = F) +
  scale_fill_manual(values = as.character(unique(m_test$cmat))) +
  geom_sf(data = st_transform(st_as_sf(cz_bound), st_crs(32633)), color = "grey30", fill = NA, show.legend = FALSE, lwd = 0.5) +
  theme_void()

unc_leg <- ggplot(as.data.frame(r, xy = T)) +
  geom_raster(aes(x = x, y = y, fill = as.factor(layer)), show.legend = F) +
  scale_fill_manual(values = as.character(unique(n_cmat$cmat))) +
  geom_point(data = xy_df, aes(x = pred, y = unc), size = 0.01, alpha = 0.05, color = "grey90") + labs(x = 'suitability', y = 'standard deviation') +
  geom_rangeframe(data = data.frame(x = c(0, 1), y = c(0, round(max(na.omit(getValues(unc))), 2))), aes(x, y)) +
  scale_x_continuous(breaks = c(0, 1)) +
  scale_y_continuous(breaks = c(0, round(max(na.omit(getValues(unc))), 2))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())

unc_leg <- ggMarginal(unc_leg, type = "density", size = 5, color = "grey20", lwd = 0.1)

pr_unc_leg <- pr_unc + inset_element(unc_leg, left = 0.8, bottom = 0.7, right = 1, top = 1)

ggsave(
  file.path("..", "results", "plots", paste0("unc_bivar_", export_dpi, ".jpeg")),
  plot = pr_unc_leg,
  width = 16.1,
  height = 8.7,
  dpi = export_dpi,
)

# Env and geo plot
var_imp <- getVarImp(m, wtest = "train")@varImportanceMean$corTest[1:2]
var_imp <- var_imp[order(var_imp[2], decreasing = T), ]

rownames(var_imp)[1] -> var1
rownames(var_imp)[2] -> var2
gg_pdata <- pdata[, c("x", "y", var1, var2, "bg_occ")]
sf_pdata$bg_occ <- factor(sf_pdata$bg_occ, 
                          levels =  c("bg_train", "occs_train", "bg_test", "occs_test"),
                          labels = c("background training points", "training occurrences", "background testing points", "testing occurrences")
                          )

env_bg_train <- ggplot(gg_pdata[gg_pdata$bg_occ == "bg_train", ]) +
    geom_point(aes_string(x = var1, y = var2), color = "grey", size = 0.04) +
    geom_point(gg_pdata[gg_pdata$bg_occ == "occs_train", ], mapping = aes_string(x = var1, y = var2), color = "darkblue", size = 0.5) +
     xlim(min(gg_pdata[var1]), max(gg_pdata[var1])) +
     ylim(min(gg_pdata[var2]), max(gg_pdata[var2])) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    ylab("Minimum soil temperature of the coldest month") +
    xlab(" ")
    
env_bg_test <- ggplot(gg_pdata[gg_pdata$bg_occ == "bg_test", ]) +
    geom_point(aes_string(x = var1, y = var2), color = "grey40", size = 0.04) +
    geom_point(gg_pdata[gg_pdata$bg_occ == "occs_test", ], mapping = aes_string(x = var1, y = var2), color = "darkred", size = 0.5) +
     xlim(min(gg_pdata[var1]), max(gg_pdata[var1])) +
     ylim(min(gg_pdata[var2]), max(gg_pdata[var2])) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab("Mean relative humidity") +
    ylab(" ")


env_bg <- ggplot() +
    geom_point(gg_pdata[gg_pdata$bg_occ == "bg_train", ], mapping = aes_string(x = var1, y = var2), color = "grey", size = 0.04) +
    geom_point(gg_pdata[gg_pdata$bg_occ == "bg_test", ], mapping = aes_string(x = var1, y = var2), color = "grey40", size = 0.04) +
    geom_point(gg_pdata[gg_pdata$bg_occ == "occs_train", ], mapping = aes_string(x = var1, y = var2), color = "darkblue", size = 0.5) +
    geom_point(gg_pdata[gg_pdata$bg_occ == "occs_test", ], mapping = aes_string(x = var1, y = var2), color = "darkred", size = 0.5) +
     xlim(min(gg_pdata[var1]), max(gg_pdata[var1])) +
     ylim(min(gg_pdata[var2]), max(gg_pdata[var2])) +
    theme_bw() +
    theme(aspect.ratio = 1) +
    xlab(" ") +
    ylab(" ")

geo_bg_train <- ggplot(sf_pdata[sf_pdata$bg_occ %in% c("background training points", "training occurrences"), ]) +
    geom_sf(aes(color = bg_occ, size = bg_occ), show.legend = F) +
    scale_size_manual(values = c(0.04, 0.4)) +
    scale_color_manual(values = c("grey", "darkblue")) +
    coord_sf(xlim = 
        c(min(st_coordinates(sf_pdata)[, 1]),
        max(st_coordinates(sf_pdata)[, 1])),
        ylim = c(min(st_coordinates(sf_pdata)[, 2]),
        max(st_coordinates(sf_pdata)[, 2]))
        )   +
    theme_bw() +
    xlab("") +
    ylab("Latitude")

geo_bg_test <- ggplot(sf_pdata[sf_pdata$bg_occ %in% c("background testing points", "testing occurrences"), ]) +
    geom_sf(aes(color = bg_occ, size = bg_occ), show.legend = F) +
    scale_size_manual(values = c(0.04, 0.4)) +
    scale_color_manual(values = c("grey40", "darkred")) +
    coord_sf(xlim = 
        c(min(st_coordinates(sf_pdata)[, 1]),
        max(st_coordinates(sf_pdata)[, 1])),
        ylim = c(min(st_coordinates(sf_pdata)[, 2]),
        max(st_coordinates(sf_pdata)[, 2]))
        ) +
    theme_bw() +
    xlab("Longitude") +
    ylab("")


geo_bg <- ggplot(sf_pdata) +
    geom_sf(aes(color = bg_occ, size = bg_occ), show.legend = T) +
    scale_size_manual(name = "", values = c(0.04, 0.04, 0.4, 0.4)) +
    scale_color_manual(name = "", values = c("grey", "darkblue", "grey40", "darkred")) +
    theme_bw() + 
    guides(colour = guide_legend(override.aes = list(size = 5)))


p_env_geo <- wrap_plots(geo_bg_train,
                        geo_bg_test,
                        geo_bg,
                        ggMarginal(env_bg_train, type = "violin", size = 30, fill = "grey", lwd = 0.1),
                        ggMarginal(env_bg_test, type = "violin", size = 30, fill = "grey40", lwd = 0.1),
                        ggMarginal(env_bg, type = "violin", size = 30, fill = "white", color = "white", lwd = 0.1)) +
                        plot_annotation(tag_levels = "A") +
                        plot_layout(guides = 'collect') & theme(legend.position = 'bottom')


ggsave(
  file.path("..", "results", "plots", paste0("env_geo_", export_dpi, ".jpeg")),
  plot = p_env_geo,
  width = 16.1,
  height = 9.7,
  dpi = export_dpi,
)

#########################################
# Response curves

library(reshape)
rcs <- getResponseCurve(m)
test_data <- data.frame(var_name = character(), var = double(), m_mean = double(), m_ci_m = double(), m_ci_l = double())
for (plot_var in rcs@variables){
    var_data <- data.frame(var = rcs@response[[plot_var]][, 1])
    var_data$var_name <- plot_var 
    var_data$m_mean <- rowMeans(rcs@response[[plot_var]][-1], na.rm = FALSE)
    m_ci <- 1.96*apply(rcs@response[[plot_var]][-1], 1, sd, na.rm = TRUE)/sqrt(100)
    var_data$m_ci_m <- var_data$m_mean + m_ci
    var_data$m_ci_l <- var_data$m_mean - m_ci
    test_data <- rbind(var_data, test_data)
    }

pdata_train <- pdata[pdata$bg_occ %in% c("bg_train", "occs_train"), -c(1:3)]
pdata_train <- melt(pdata_train)
names(pdata_train)[2] <- "var_name"

ggcti_rc <- ggplot(test_data, aes(x = var, y = m_mean)) +
    geom_ribbon(aes(ymin = m_ci_l, ymax = m_ci_m), fill = "grey70") + geom_line() +
    geom_rug(data = pdata_train[pdata_train$bg_occ == "occs_train", ], mapping = aes(x = value), color = "darkblue", inherit.aes = F, sides = "t", alpha = 0.1) +
    geom_rug(data = pdata_train[pdata_train$bg_occ == "bg_train", ], mapping = aes(x = value), color = "grey50", inherit.aes = F, sides = "b", alpha = 0.1)

resp <- ggcti_rc  + facet_wrap(~var_name, scales = "free_x", nrow = 2) +labs(caption = "Response curves of variables. Grey rugs - train background values density; Blue rugs - train precence values density")

ggsave(
  file.path("..", "results", "plots", paste0("response_", export_dpi, ".jpeg")),
  plot = resp,
  width = 16.1,
  height = 8.7,
  dpi = export_dpi,
)
