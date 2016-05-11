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
source("init_cluster_parameters.R")

bs_files  <- c("../datasets/BEATSON/BS-Seq/9L_df_young.CpG_context.aggregate",
                 "../datasets/BEATSON/BS-Seq/10L_df_young.CpG_context.aggregate",
                 "../datasets/BEATSON/BS-Seq/11L_df_young.CpG_context.aggregate",
                 "../datasets/BEATSON/BS-Seq/12L_df_young.CpG_context.aggregate")
rna_files <- "../datasets/BEATSON/RNA-Seq/processed/RNA_Seq_FPKM_X9L_df_young.bed"


# -----------------------------------------
# Apply the BPR model for clustering
# -----------------------------------------
source("bpr_cluster.R")


# -----------------------------------------
# Store the results
# -----------------------------------------
filename <- paste0("../files/cluster_df_young_", downstream, "_",
                   format(Sys.time(), "%a%b%d%H%M"),
                   ".RData")
save(HTS_data, proc_data, bpr_model, file = filename)
