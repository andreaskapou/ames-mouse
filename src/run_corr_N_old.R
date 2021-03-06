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

bs_files  <- c("../datasets/BEATSON/BS-Seq/5L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/6L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/7L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/8L_N_old.CpG_context.aggregate")
rna_files <- "../datasets/BEATSON/RNA-Seq/processed/RNA_Seq_FPKM_X5L_N_old.bed"


# -----------------------------------------
# Apply the BPR model for regression
# -----------------------------------------
source("bpr_regression.R")


# --------------------------------------
# Store the results
# --------------------------------------
filename <- paste0("../files/corr_N_old_", downstream, "_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, out_prof, out_mean, file = filename)
