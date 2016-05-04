# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(processHTS)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE)

# ------------------------------------------
# Initialize parameters
# ------------------------------------------
bs_files  <- c("../datasets/BEATSON/BS-Seq/1L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/2L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/3L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/4L_df_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/5L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/6L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/7L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/8L_N_old.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/9L_df_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/10L_df_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/11L_df_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/12L_df_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/13L_N_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/14L_N_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/15L_N_young.CpG_context.aggregate",
               "../datasets/BEATSON/BS-Seq/16L_N_young.CpG_context.aggregate"
)

for (i in seq(1, length(bs_files), 2)){
  bs_data <- pool_bs_bismark_cov_rep(bs_files[i:(i+1)])
  
  filename = paste0(bs_files[i], ".pool.bed")
  final_data = cbind(chr         = as.character(seqnames(bs_data)), 
                     chromStart  = start(bs_data), 
                     chromEnd    = start(bs_data) + 1, 
                     name        = rep(1, length(start(bs_data))), 
                     score       = bs_data$meth_reads + bs_data$unmeth_reads,
                     strand      = rep("*", length(start(bs_data))),
                     thickStart  = rep(1, length(start(bs_data))),
                     thickEnd    = rep(1, length(start(bs_data))), 
                     itemRgb     = rep(1, length(start(bs_data))), 
                     readCount   = bs_data$meth_reads + bs_data$unmeth_reads,
                     percentMeth = ceiling(100 * (bs_data$meth_reads / (bs_data$meth_reads + bs_data$unmeth_reads)))
  )
  write.table(x = final_data, 
              file = filename, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}
rm(final_data)