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

# DF Old parameters
df_old_W <- data.frame(x = df_old_out_mean$W_opt,
                       y = df_old_Y)
# DF Young parameters
df_young_W <- data.frame(x = df_young_out_mean$W_opt,
                         y = df_young_Y)
# N Old parameters
N_old_W <- data.frame(x = N_old_out_mean$W_opt,
                      y = N_old_Y)
# N Young parameters
N_young_W <- data.frame(x = N_young_out_mean$W_opt,
                        y = N_young_Y)


# From DF OLD
message("Predicting from DF OLD")
DO_predict_DY <- predict_model_gex(model      = df_old_out_mean$gex_model,
                                   test       = df_young_W,
                                   is_summary = FALSE)
DO_predict_NO <- predict_model_gex(model      = df_old_out_mean$gex_model,
                                   test       = N_old_W,
                                   is_summary = FALSE)
DO_predict_NY <- predict_model_gex(model      = df_old_out_mean$gex_model,
                                   test       = N_young_W,
                                   is_summary = FALSE)


# From DF Young
message("Predicting from DF YOUNG")
DY_predict_DO <- predict_model_gex(model      = df_young_out_mean$gex_model,
                                   test       = df_old_W,
                                   is_summary = FALSE)
DY_predict_NO <- predict_model_gex(model      = df_young_out_mean$gex_model,
                                   test       = N_old_W,
                                   is_summary = FALSE)
DY_predict_NY <- predict_model_gex(model      = df_young_out_mean$gex_model,
                                   test       = N_young_W,
                                   is_summary = FALSE)


# From N OLD
message("Predicting from N OLD")
NO_predict_NY <- predict_model_gex(model      = N_old_out_mean$gex_model,
                                   test       = N_young_W,
                                   is_summary = FALSE)
NO_predict_DO <- predict_model_gex(model      = N_old_out_mean$gex_model,
                                   test       = df_old_W,
                                   is_summary = FALSE)
NO_predict_DY <- predict_model_gex(model      = N_old_out_mean$gex_model,
                                   test       = df_young_W,
                                   is_summary = FALSE)


# From N Young
message("Predicting from N YOUNG")
NY_predict_NO <- predict_model_gex(model      = N_young_out_mean$gex_model,
                                   test       = N_old_W,
                                   is_summary = FALSE)
NY_predict_DO <- predict_model_gex(model      = N_young_out_mean$gex_model,
                                   test       = df_old_W,
                                   is_summary = FALSE)
NY_predict_DY <- predict_model_gex(model      = N_young_out_mean$gex_model,
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


# ---------- Make plots -------------------------------------
out_DO_to_DY$test_errors <- DO_predict_DY$test_errors
ggplot_scatt_across_cell_line(output = out_DO_to_DY, 
                              main_lab = expression(Mean~DO %->% DY), 
                              is_margins = TRUE)

confus_mean <- plot_ames_conf_corr_matrix(df_old_out_mean, out_DO_to_DY, out_DO_to_NO, out_DO_to_NY,  
                                          df_young_out_mean, out_DY_to_DO, out_DY_to_NO, out_DY_to_NY,
                                          N_old_out_mean, out_NO_to_DY, out_NO_to_DO, out_NO_to_NY, 
                                          N_young_out_mean, out_NY_to_DY, out_NY_to_DO, out_NY_to_NO, 
                                          title_lab = "Mean Methylation Correlation")