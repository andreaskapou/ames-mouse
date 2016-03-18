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


# ------------------------------------------
# Initialize parameters
# ------------------------------------------
source("init_regr_parameters.R")

bs_files  <- c("../datasets/BEATSON/BS-Seq/1L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/2L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/3L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/4L_df_old.CpG_context.aggregate")
rna_files <- "../datasets/BEATSON/RNA-Seq/processed/RNA_Seq_FPKM_X1L_df_old.bed"


# -----------------------------------------
# Apply the BPR model for regression
# -----------------------------------------
source("bpr_regression.R")


# -----------------------------------------
# Store the results
# -----------------------------------------
filename <- paste0("../files/corr_df_old_", 
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, out_prof, out_mean, file = filename)
