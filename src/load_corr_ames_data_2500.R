# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE)


##----------------DF OLD mice results --------------
load("../files/corr_df_old_2500_WedMay111737.RData")
df_old_HTS_data <- HTS_data
df_old_out_mean <- out_mean
df_old_out_prof <- out_prof
df_old_X        <- proc_data$X
df_old_Y        <- proc_data$Y
df_old_genes    <- proc_data$genes

##-----------------DF YOUNG mice results ------------
load("../files/corr_df_young_2500_WedMay111737.RData")
df_young_HTS_data <- HTS_data
df_young_out_mean <- out_mean
df_young_out_prof <- out_prof
df_young_X        <- proc_data$X
df_young_Y        <- proc_data$Y
df_young_genes    <- proc_data$genes


##------------- Normal OLD mice results ------------
load("../files/corr_N_old_2500_WedMay111736.RData")
N_old_HTS_data <- HTS_data
N_old_out_mean <- out_mean
N_old_out_prof <- out_prof
N_old_X        <- proc_data$X
N_old_Y        <- proc_data$Y
N_old_genes    <- proc_data$genes


##------------ Normal YOUNG mice results ----------
load("../files/corr_N_young_2500_WedMay111735.RData")
N_young_HTS_data <- HTS_data
N_young_out_mean <- out_mean
N_young_out_prof <- out_prof
N_young_X        <- proc_data$X
N_young_Y        <- proc_data$Y
N_young_genes    <- proc_data$genes
