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
               "../datasets/BEATSON/BS-Seq/4L_df_old.CpG_context.aggregate")

bs_data <- read_bs_bismark_cov(bs_files[1], is_GRanges = FALSE)

filename = paste0(bs_files[1], ".bed")
data = cbind(chr         = bs_data$chr, 
             chromStart  = bs_data$start, 
             chromEnd    = bs_data$start + 1, 
             name        = rep("N", length(bs_data$chr)), 
             score       = bs_data$meth_reads + bs_data$unmeth_reads,
             strand      = rep(".", length(bs_data$chr)),
             thickStart  = rep(".", length(bs_data$chr)),
             thickEnd    = rep(".", length(bs_data$chr)), 
             itemRgb     = rep(".", length(bs_data$chr)), 
             readCount   = bs_data$meth_reads + bs_data$unmeth_reads,
             percentMeth = ceiling(100 * (bs_data$meth_reads / (bs_data$meth_reads + bs_data$unmeth_reads)))
)
write.table(x = data, 
            file = filename, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
