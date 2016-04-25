# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
library(earth)
library(e1071)
library(randomForest)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE) 
# 
# 
# ##----------- Parameters for filtering data --------
# gene_expr_thresh <- FALSE
# gene_outl_thresh <- TRUE
# gene_log2_transf <- TRUE
# max_outl <- 600
# 
# 
# 
# ##----------------DF OLD mice results --------------
# load("../files/corr_df_old_SunFeb211527.RData")
# df_old_basis_prof <- basis_prof
# df_old_basis_mean <- basis_mean
# df_old_HTS_data <- HTS_data
# df_old_out_mean <- out_mean
# df_old_out_prof <- out_prof
# 
# proc_data <- preprocess_data(HTS_data = df_old_HTS_data,
#                              max_outl = max_outl,
#                              gene_expr_thresh = gene_expr_thresh,
#                              gene_outl_thresh = gene_outl_thresh,
#                              gene_log2_transf = gene_log2_transf)
# df_old_obs <- proc_data$obs
# df_old_Y <- proc_data$Y
# 
# 
# 
# ##-----------------DF YOUNG mice results ------------
# load("../files/corr_df_young_SunFeb211525.RData")
# df_young_basis_prof <- basis_prof
# df_young_basis_mean <- basis_mean
# df_young_HTS_data <- HTS_data
# df_young_out_mean <- out_mean
# df_young_out_prof <- out_prof
# 
# proc_data <- preprocess_data(HTS_data = df_young_HTS_data,
#                              max_outl = max_outl,
#                              gene_expr_thresh = gene_expr_thresh,
#                              gene_outl_thresh = gene_outl_thresh,
#                              gene_log2_transf = gene_log2_transf)
# df_young_obs <- proc_data$obs
# df_young_Y <- proc_data$Y
# 
# 
# 
# ##------------- Normal OLD mice results ------------
# load("../files/corr_N_old_MonFeb221046.RData")
# N_old_basis_prof <- basis_prof
# N_old_basis_mean <- basis_mean
# N_old_HTS_data <- HTS_data
# N_old_out_mean <- out_mean
# N_old_out_prof <- out_prof
# 
# proc_data <- preprocess_data(HTS_data = N_old_HTS_data,
#                              max_outl = max_outl,
#                              gene_expr_thresh = gene_expr_thresh,
#                              gene_outl_thresh = gene_outl_thresh,
#                              gene_log2_transf = gene_log2_transf)
# N_old_obs <- proc_data$obs
# N_old_Y <- proc_data$Y
# 
# 
# 
# ##------------ Normal YOUNG mice results ----------
# load("../files/corr_N_young_MonFeb221047.RData")
# N_young_basis_prof <- basis_prof
# N_young_basis_mean <- basis_mean
# N_young_HTS_data <- HTS_data
# N_young_out_mean <- out_mean
# N_young_out_prof <- out_prof
# 
# proc_data <- preprocess_data(HTS_data = N_young_HTS_data,
#                              max_outl = max_outl,
#                              gene_expr_thresh = gene_expr_thresh,
#                              gene_outl_thresh = gene_outl_thresh,
#                              gene_log2_transf = gene_log2_transf)
# N_young_obs <- proc_data$obs
# N_young_Y <- proc_data$Y


# DF Old parameters
df_old_W <- data.frame(x = df_old_out_prof$W_opt,
                       y = df_old_Y)

# DF Young parameters
df_young_W <- data.frame(x = df_young_out_prof$W_opt,
                         y = df_young_Y)

# N Old parameters
N_old_W <- data.frame(x = N_old_out_prof$W_opt,
                      y = N_old_Y)

# N Young parameters
N_young_W <- data.frame(x = N_young_out_prof$W_opt,
                        y = N_young_Y)


# From DF OLD
message("Predicting from DF OLD")
DO_predict_DY <- predict_model_gex(model      = df_old_out_prof$gex_model,
                                   test       = df_young_W,
                                   is_summary = FALSE)

DO_predict_NO <- predict_model_gex(model      = df_old_out_prof$gex_model,
                                   test       = N_old_W,
                                   is_summary = FALSE)

DO_predict_NY <- predict_model_gex(model      = df_old_out_prof$gex_model,
                                   test       = N_young_W,
                                   is_summary = FALSE)


# From DF Young
message("Predicting from DF YOUNG")
DY_predict_DO <- predict_model_gex(model      = df_young_out_prof$gex_model,
                                   test       = df_old_W,
                                   is_summary = FALSE)

DY_predict_NO <- predict_model_gex(model      = df_young_out_prof$gex_model,
                                   test       = N_old_W,
                                   is_summary = FALSE)

DY_predict_NY <- predict_model_gex(model      = df_young_out_prof$gex_model,
                                   test       = N_young_W,
                                   is_summary = FALSE)



# From N OLD
message("Predicting from N OLD")
NO_predict_NY <- predict_model_gex(model      = N_old_out_prof$gex_model,
                                   test       = N_young_W,
                                   is_summary = FALSE)

NO_predict_DO <- predict_model_gex(model      = N_old_out_prof$gex_model,
                                   test       = df_old_W,
                                   is_summary = FALSE)

NO_predict_DY <- predict_model_gex(model      = N_old_out_prof$gex_model,
                                   test       = df_young_W,
                                   is_summary = FALSE)



# From N Young
message("Predicting from N YOUNG")
NY_predict_NO <- predict_model_gex(model      = N_young_out_prof$gex_model,
                                   test       = N_old_W,
                                   is_summary = FALSE)

NY_predict_DO <- predict_model_gex(model      = N_young_out_prof$gex_model,
                                   test       = df_old_W,
                                   is_summary = FALSE)

NY_predict_DY <- predict_model_gex(model      = N_young_out_prof$gex_model,
                                   test       = df_young_W,
                                   is_summary = FALSE)



#--------------- Create final data for plotting ---------------------
# DO to all other mouse models
out_DO_to_DY <- list(test_pred = DO_predict_DY$test_pred,
                     test = list(y = df_young_Y))

out_DO_to_NO <- list(test_pred = DO_predict_NO$test_pred,
                     test = list(y = N_old_Y))

out_DO_to_NY <- list(test_pred = DO_predict_NY$test_pred,
                     test = list(y = N_young_Y))


# DY to all other mouse models
out_DY_to_DO <- list(test_pred = DY_predict_DO$test_pred,
                     test = list(y = df_old_Y))

out_DY_to_NO <- list(test_pred = DY_predict_NO$test_pred,
                     test = list(y = N_old_Y))

out_DY_to_NY <- list(test_pred = DY_predict_NY$test_pred,
                     test = list(y = N_young_Y))


# NO to all other mouse models
out_NO_to_DY <- list(test_pred = NO_predict_DY$test_pred,
                     test = list(y = df_young_Y))

out_NO_to_DO <- list(test_pred = NO_predict_DO$test_pred,
                     test = list(y = df_old_Y))

out_NO_to_NY <- list(test_pred = NO_predict_NY$test_pred,
                     test = list(y = N_young_Y))


# NY to all other mouse models
out_NY_to_DY <- list(test_pred = NY_predict_DY$test_pred,
                     test = list(y = df_young_Y))

out_NY_to_DO <- list(test_pred = NY_predict_DO$test_pred,
                     test = list(y = df_old_Y))

out_NY_to_NO <- list(test_pred = NY_predict_NO$test_pred,
                     test = list(y = N_old_Y))



round(cor(out_DO_to_DY$test_pred, out_DO_to_DY$test$y),3)
round(cor(out_DO_to_NO$test_pred, out_DO_to_NO$test$y),3)
round(cor(out_DO_to_NY$test_pred, out_DO_to_NY$test$y),3)

round(cor(out_DY_to_DO$test_pred, out_DY_to_DO$test$y),3)
round(cor(out_DY_to_NO$test_pred, out_DY_to_NO$test$y),3)
round(cor(out_DY_to_NY$test_pred, out_DY_to_NY$test$y),3)

round(cor(out_NO_to_DY$test_pred, out_NO_to_DY$test$y),3)
round(cor(out_NO_to_DO$test_pred, out_NO_to_DO$test$y),3)
round(cor(out_NO_to_NY$test_pred, out_NO_to_NY$test$y),3)

round(cor(out_NY_to_DY$test_pred, out_NY_to_DY$test$y),3)
round(cor(out_NY_to_DO$test_pred, out_NY_to_DO$test$y),3)
round(cor(out_NY_to_NO$test_pred, out_NY_to_NO$test$y),3)