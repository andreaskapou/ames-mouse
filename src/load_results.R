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


##----------- Parameters for filtering data --------
gene_expr_thresh <- FALSE
gene_outl_thresh <- TRUE
gene_log2_transf <- TRUE
max_outl <- 600



##----------------DF OLD mice results --------------
load("../files/corr_df_old_SunFeb282233.RData")
df_old_basis_prof <- basis_prof
df_old_basis_mean <- basis_mean
df_old_HTS_data <- HTS_data
df_old_out_mean <- out_mean
df_old_out_prof <- out_prof

proc_data <- preprocess_data(HTS_data = df_old_HTS_data,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
df_old_obs <- proc_data$obs
df_old_Y <- proc_data$Y



##-----------------DF YOUNG mice results ------------
load("../files/corr_df_young_SunFeb282215.RData")
df_young_basis_prof <- basis_prof
df_young_basis_mean <- basis_mean
df_young_HTS_data <- HTS_data
df_young_out_mean <- out_mean
df_young_out_prof <- out_prof

proc_data <- preprocess_data(HTS_data = df_young_HTS_data,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
df_young_obs <- proc_data$obs
df_young_Y <- proc_data$Y



##------------- Normal OLD mice results ------------
load("../files/corr_N_old_MonFeb290004.RData")
N_old_basis_prof <- basis_prof
N_old_basis_mean <- basis_mean
N_old_HTS_data <- HTS_data
N_old_out_mean <- out_mean
N_old_out_prof <- out_prof

proc_data <- preprocess_data(HTS_data = N_old_HTS_data,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
N_old_obs <- proc_data$obs
N_old_Y <- proc_data$Y



##------------ Normal YOUNG mice results ----------
load("../files/corr_N_young_MonFeb290041.RData")
N_young_basis_prof <- basis_prof
N_young_basis_mean <- basis_mean
N_young_HTS_data <- HTS_data
N_young_out_mean <- out_mean
N_young_out_prof <- out_prof

proc_data <- preprocess_data(HTS_data = N_young_HTS_data,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
N_young_obs <- proc_data$obs
N_young_Y <- proc_data$Y