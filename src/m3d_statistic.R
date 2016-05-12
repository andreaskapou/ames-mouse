# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(processHTS)
library(M3Ddevel)
library(BiSeq)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE)


# bs_files  <- c("../datasets/BEATSON/BS-Seq/1L_df_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/2L_df_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/3L_df_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/4L_df_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/5L_N_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/6L_N_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/7L_N_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/8L_N_old.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/9L_df_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/10L_df_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/11L_df_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/12L_df_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/13L_N_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/14L_N_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/15L_N_young.CpG_context.aggregate",
#                "../datasets/BEATSON/BS-Seq/16L_N_young.CpG_context.aggregate")

# bs_files  <- c("../datasets/BEATSON/BS-Seq/processed/1L_df_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/2L_df_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/3L_df_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/4L_df_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/5L_N_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/6L_N_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/7L_N_old.CpG_context.aggregate.bed",
#                "../datasets/BEATSON/BS-Seq/processed/8L_N_old.CpG_context.aggregate.bed")

# Files to read
bs_files  <- c("../datasets/BEATSON/BS-Seq/processed/1L_df_old.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/2L_df_old.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/3L_N_old.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/4L_N_old.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/5L_df_young.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/6L_df_young.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/7L_N_young.pool.bed",
               "../datasets/BEATSON/BS-Seq/processed/8L_N_young.pool.bed")
rna_files  <- "../datasets/BEATSON/RNA-Seq/processed/RNA_Seq_FPKM_X1L_df_old.bed"
#annot_data <- read_rna_beatson("../datasets/BEATSON/RNA-Seq/processed/RNA_Seq_FPKM_X1L_df_old.bed")

# Create annotation data for experiment
group <- factor(c("DF_Old", "DF_Old", 
                  "N_Old", "N_Old"))
samples <- c("DF_Old1", "DF_Old2",
             "N_Old1", "N_Old2")
colData <- DataFrame(group, row.names= samples)

# Read ENCODE data and store them in the required format
wgbs <- readENCODEdata(bs_files[1:4], colData)

# Read BS data
bs_data  <- preprocess_bs_encode_haib(files = bs_files[1])
# Read RNA data
rna_data <- read_rna_beatson(file = rna_files, is_GRanges = TRUE)
# Create promooter regions
downstream <- 5000
prom_regions <- create_prom_region(annot_data = rna_data, 
                                   upstream   = -downstream, 
                                   downstream = downstream)
# Create methylation regions
methyl_region <- create_methyl_region(bs_data       = bs_data,
                                      prom_region   = prom_regions,
                                      cpg_density   = 20)

# Keep only the corresponding gene annotation data
prom_regions <- prom_regions[methyl_region$prom_ind]

# Find overlaps between promoter regions and WGBS data
overlaps <- findOverlaps(prom_regions, rowRanges(wgbs))

# Perform 
MMD_data <- M3D_Wrapper_lite(wgbs, overlaps)

save(bs_files,
     colData,
     wgbs,
     prom_regions,
     overlaps,
     MMD_data, file = paste0("../datasets/BEATSON/M3D_", downstream, ".RData"))

