library(raster)
library(dismo)
library(rgdal)
library(usdm)
library(sdm)
library(ecospat)
library(ENMeval)
require(devtools)
    source_url("https://raw.githubusercontent.com/RWunderlich/SEDI/master/R/sedi.R")

cz_bound <- readRDS(file.path("..", "data", "admin", "gadm36_CZE_0_sp.rds"))

# Env
env_source <- list.files(file.path("..", "processed_data", "env"), full.names = T)
env_data <- stack(env_source)

cz <- crop(env_data, cz_bound)
env_data_bg <- mask(env_data, cz_bound, inverse = T)

# VIF
env_vif <- vifcor(env_data_bg, th = .7)

env_data_bg <- exclude(env_data_bg, env_vif)
env_data <- exclude(env_data, env_vif)
cz <- exclude(cz, env_vif)

# Occ
spec_data <- readOGR(file.path("..", "processed_data", "thinned_sp_10.gpkg"))
colnames(spec_data@coords) <- c("x", "y")
spec_data$sp_name <- 1
sp_train <- spec_data[, "sp_name"]

sp_test_cz <- intersect(sp_train, cz_bound)[, "sp_name"]
sp_train_eur <- sp_train[!(sp_train@coords[, 1] %in% sp_test_cz@coords[, 1])&!(sp_train@coords[, 2] %in% sp_test_cz@coords[, 2]), ][, "sp_name"]

# Bg points
bg_sp <- randomPoints(env_data_bg, 30000)
bg <- extract(env_data_bg, bg_sp)
bg <- cbind(bg_sp, bg)
bg <- as.data.frame(bg)

bg_testing <- 10000
test_bg <- randomPoints(mask(crop(cz, cz_bound), cz_bound), bg_testing)
test_bg <- SpatialPointsDataFrame(test_bg, data = data.frame("sp_name" = replicate(bg_testing, 0)))
crs(test_bg) <- crs(sp_test_cz)
sp_test_cz <- rbind(sp_test_cz, test_bg)

# ExDet
test_bg_sp <- mask(cz, cz_bound)
climan <- ecospat.climan(bg[-c(1:2)], getValues(test_bg_sp))

climan_r <- test_bg_sp[[1]]*0
names(climan_r) <- "climan"
climan_r@data@values <- climan 
summary(climan)

# ENMeval
rms <- seq(0.5, 10, 0.5)
# for better reproducibility we rewrote the code for ENMeval version 2, however we need to remove "LQTP" feature class, because it is not accepted for some reason, this not affect the results since the final model use "LQH" features

fcs <- c("LQ", "LQP", "LQT", "LQH", "LQHT", "LQHP", "LQHPT")
kfold_part <- get.randomkfold(sp_train@coords, bg_sp, 5)

enm_eval <- ENMevaluate(occ = sp_train@coords,
                        env = env_data,
                        bg.coords = bg_sp,
                        method = 'randomkfold',
                        partitions = kfold_part,
                        tune.args = list(fc = fcs, rm = rms),
                        algorithm = 'maxent.jar',
                        parallel = F)


enm_eval_results <- eval.results(enm_eval)

# Code used for model in published paper, needs ENMeval version <2

# fcs <- c("LQ", "LQP", "LQT", "LQH", "LQHT", "LQTP", "LQHP", "LQHPT")
# enm_eval_results <- data.frame()
# for (i in rms){
#     enm_eval <- ENMevaluate(occ = sp_train@coords, env = env_data, bg.coords = bg_sp, method = 'randomkfold', kfolds = 5, fc = fcs, RMvalues = i, algorithm = 'maxent.jar', parallel = F)
#     enm_eval_results <- rbind(enm_eval@results, enm_eval_results)
#     enm_eval_results$delta.AICc <- enm_eval_results$AICc - min(enm_eval_results$AICc)
#     print(enm_eval_results[c(1, 2, 13, 14)])
#     write.csv(enm_eval_results, file.path("..", "results", "ENMeval.csv")) 
#     }

write.csv(enm_eval_results, file.path("..", "results", "ENMeval.csv"))

comb <- enm_eval_results$features[enm_eval_results$delta.AICc == 0]
b_beta <- enm_eval_results$rm[enm_eval_results$delta.AICc == 0]

# tune arguments used in fnial model:
# comb <- "LQH"
# b_beta <- 1.5

b_args <- c()
if (!grepl("L", comb)){
    b_args <- c(b_args, "nolinear")
}
if (!grepl("H", comb)){
    b_args <- c(b_args, "nohinge")
}
if (!grepl("Q", comb)){
    b_args <- c(b_args, "noquadratic")
}
if (!grepl("P", comb)){
    b_args <- c(b_args, "noproduct")
}
if (grepl("T", comb)){
    b_args <- c(b_args, "threshold")
}

# Model

# For load published model just run: 
# m <- readRDS(file.path("..", "results", "model.RDS"))

d <- sdmData(train = sp_train_eur, test = sp_test_cz, predictors = env_data, bg = bg)

m <- sdm(sp_name~., data = d,
        methods = c('maxent'),
        replication = c('cv'),
        cv.folds = 5, 
        modelSettings = list(maxent = list(beta = b_beta, args = c(b_args, "noremoveDuplicates", "noautofeature"))), 
        n = 10
        )

saveRDS(m, file.path("..", "results", "model.RDS")) 
m

## Evaluation
mean_ev <- function (mdl = m, testing_set = "test.dep", measures = measur){
    df <- as.data.frame(apply(getEvaluation(mdl, w = (getModelId(mdl)), stat = measures, opt = 2, wtest = testing_set)[-1], 2, mean))
    df <- cbind(df, as.data.frame(apply(getEvaluation(mdl, w = (getModelId(mdl)), stat = measures, opt = 2, wtest = testing_set)[-1], 2, sd)))
    df <- cbind(df, as.data.frame(apply(getEvaluation(mdl, w = (getModelId(mdl)), stat = measures, opt = 2, wtest = testing_set)[-1], 2, max)))
    df <- cbind(df, as.data.frame(apply(getEvaluation(mdl, w = (getModelId(mdl)), stat = measures, opt = 2, wtest = testing_set)[-1], 2, min)))
    
    colnames(df) <- c("mean", "sd", "max", "min")
    if ("boyce" %in% measures){
        print("calculating boyce")
        boyce_li <- c()
        for (repl in getModelId(mdl)){
            cv_bg <- c(mdl@models$sp_name$maxent[[repl]]@evaluation$train@predicted[mdl@models$sp_name$maxent[[repl]]@evaluation$train@observed == 0], mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@predicted[mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@observed == 0])
            cv_pres <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@predicted[mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@observed == 1]
            b <- ecospat.boyce(cv_bg, cv_pres, PEplot = F)
            boyce_li <- c(boyce_li, b$cor)
        }
        df <- rbind(df, c(mean(boyce_li), sd(boyce_li), max(boyce_li), min(boyce_li)))
        rownames(df)[nrow(df)] <- "boyce" 
    }
    if ("SEDI" %in% measures){
        print("calculating SEDI")
        sedi_li <- c()

        for (repl in getModelId(mdl)){
            th_val <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@threshold_based[2, 2]
            o <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@observed
            p <- mdl@models$sp_name$maxent[[repl]]@evaluation[[testing_set]]@predicted
            cmx <- sdm:::.cmx(o>th_val, p>th_val) 
            sedi_li <- c(sedi_li, sedi(cmx[1, 1], cmx[1, 2], cmx[2, 1], cmx[2, 2])[[2]])
        }
        df <- rbind(df, c(mean(sedi_li), sd(sedi_li), max(sedi_li), min(sedi_li)))
        rownames(df)[nrow(df)] <- "sedi_li" 
    }

    return(df)
}

measur <- c('boyce', 'SEDI', 'sensitivity', 'threshold')

eval_dep <- mean_ev(m, "test.dep")
eval_dep
eval_ind <- mean_ev(m, "test.indep")
eval_ind

write.csv(eval_dep, file.path("..", "results", "eval_dep.csv"))
write.csv(eval_ind, file.path("..", "results", "eval_indep.csv"))


## Variable imprtance check
getVarImp(m, wtest = "test.dep")
## Response curves check
rcurve(m, id = getModelId(m))

## Predict - ensemble

e <- ensemble(m, cz, setting = list(method = 'weighted', stat = 'AUC', opt = 2))
plot(e)
writeRaster(e_tr, file.path("..", "results", "ensem_th.grd"), overwrite = T)

e_tr <- ensemble(m, cz, setting = list(method = 'pa', opt = 2))
plot(e_tr)
writeRaster(e, file.path("..", "results", "ensem.grd"), overwrite = T)

pr_st <- stack()
for (i in getModelId(m)){
    pr <- predict(m, cz, w = i)
    pr_st <- stack(pr_st, pr)
    }

pr_sd <- calc(pr_st, fun = sd)

writeRaster(pr_sd, file.path("..", "results", "unc_sd.grd"), overwrite = T)

pr_sd_stats <- data.frame(stat = c("min", "max", "mean", "sd"), 
                        val = c(min(getValues(pr_sd), na.rm = T), 
                              max(getValues(pr_sd), na.rm = T), 
                              mean(getValues(pr_sd), na.rm = T), 
                              sd(getValues(pr_sd), na.rm = T)))

write.csv(pr_sd_stats, file.path("..", "results", "pr_sd_stats.csv"))

### Testing ensembles
# binomial test
pred <- e_tr == 1
pred <- mask(pred, cz_bound)
plot(pred, col = c("white", "grey"));plot(cz_bound, add = T, border = "grey40");plot(sp_test_cz[sp_test_cz$sp_name == 1, ], add = T, pch = 20, col = "darkred")

prop <- table(pred@data@values)[2]/length(na.omit(pred@data@values))
rate <- table(extract(pred, sp_test_cz[sp_test_cz$sp_name == 1, ]))

b_t <- binom.test(rate[[2]], rate[[2]]+rate[[1]], prop)

# boyce
b <- ecospat.boyce(e, sp_test_cz[sp_test_cz$sp_name == 1, ]@coords, PEplot = F)

# SEDI
o <- sp_test_cz$sp_name
p <- extract(pred, sp_test_cz)
cmx <- sdm:::.cmx(o, p) 
SEDI <- sedi(cmx[1, 1], cmx[1, 2], cmx[2, 1], cmx[2, 2])[[2]]

ens_tests <- c(b_t$p.value, b_t$null.value[[1]], b_t$estimate[[1]], b$cor, SEDI)
eval_nam <- c("binom_p", "binom_prob", "sensitivity", "boyce", "SEDI")
eval_ens <- data.frame(round(ens_tests, 6), row.names = eval_nam)

eval_ens
write.csv(eval_ens, file.path("..", "results", "eval_ens.csv"))
