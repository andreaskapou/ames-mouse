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

##----------------- Initialize -----------------------
downstream <- 7000
# Read BS data
bs_data  <- preprocess_bs_seq(files = bs_files[1],
                              file_format = "encode_rrbs")
# Read RNA data
rna_data <- read_rna_beatson(file = rna_files, is_GRanges = TRUE)
# Create promoter regions
prom_regions <- create_prom_region(annot_data = rna_data, 
                                   upstream   = -downstream, 
                                   downstream = downstream)
# Create methylation regions
methyl_region <- create_methyl_region(bs_data       = bs_data,
                                      prom_region   = prom_regions,
                                      cpg_density   = 20)
# Keep only the corresponding gene annotation data
prom_regions <- prom_regions[methyl_region$prom_ind]


# ##------------------------------------------------------------------------------
# # Create annotation data for experiment DO-NO
# ##------------------------------------------------------------------------------
# group_do_no <- factor(c("DF_Old", "DF_Old", 
#                         "N_Old", "N_Old"))
# samples_do_no <- c("DF_Old1", "DF_Old2",
#                    "N_Old1", "N_Old2")
# colData_do_no <- DataFrame(group_do_no, row.names= samples_do_no)
# 
# # Read ENCODE data and store them in the required format
# wgbs_do_no <- readENCODEdata(bs_files[1:4], colData_do_no)
# 
# # Find overlaps between promoter regions and WGBS data
# overlaps_do_no <- findOverlaps(prom_regions, rowRanges(wgbs_do_no))
# 
# # Perform M3D statistic
# MMD_data_do_no <- M3D_Wrapper_lite(wgbs_do_no, overlaps_do_no)



##------------------------------------------------------------------------------
# Create annotation data for experiment DY-NY
##------------------------------------------------------------------------------
group_dy_ny <- factor(c("DF_Young", "DF_Young", 
                        "N_Young", "N_Young"))
samples_dy_ny <- c("DF_Young1", "DF_Young2",
                   "N_Young1", "N_Young2")
colData_dy_ny <- DataFrame(group_dy_ny, row.names= samples_dy_ny)

# Read ENCODE data and store them in the required format
wgbs_dy_ny <- readENCODEdata(bs_files[5:8], colData_dy_ny)

# Find overlaps between promoter regions and WGBS data
overlaps_dy_ny <- findOverlaps(prom_regions, rowRanges(wgbs_dy_ny))

# Perform M3D statistic
MMD_data_dy_ny <- M3D_Wrapper_lite(wgbs_dy_ny, overlaps_dy_ny)



# ##------------------------------------------------------------------------------
# # Create annotation data for experiment DY-DO
# ##------------------------------------------------------------------------------
# group_dy_do <- factor(c("DF_Young", "DF_Young", 
#                         "DF_Old", "DF_Old"))
# samples_dy_do <- c("DF_Young1", "DF_Young2",
#                    "DF_Old1", "DF_Old2")
# colData_dy_do <- DataFrame(group_dy_do, row.names= samples_dy_do)
# 
# # Read ENCODE data and store them in the required format
# wgbs_dy_do <- readENCODEdata(c(bs_files[5:6], bs_files[1:2]), colData_dy_do)
# 
# # Find overlaps between promoter regions and WGBS data
# overlaps_dy_do <- findOverlaps(prom_regions, rowRanges(wgbs_dy_do))
# 
# # Perform M3D statistic
# MMD_data_dy_do <- M3D_Wrapper_lite(wgbs_dy_do, overlaps_dy_do)
# 
# 
# 
# ##------------------------------------------------------------------------------
# # Create annotation data for experiment NO-NY
# ##------------------------------------------------------------------------------
# group_no_ny <- factor(c("N_Old", "N_Old", 
#                         "N_Young", "N_Young"))
# samples_no_ny <- c("N_Old1", "N_Old2",
#                    "N_Young1", "N_Young2")
# colData_no_ny <- DataFrame(group_no_ny, row.names= samples_no_ny)
# 
# # Read ENCODE data and store them in the required format
# wgbs_no_ny <- readENCODEdata(c(bs_files[3:4], bs_files[7:8]), colData_no_ny)
# 
# # Find overlaps between promoter regions and WGBS data
# overlaps_no_ny <- findOverlaps(prom_regions, rowRanges(wgbs_no_ny))
# 
# # Perform M3D statistic
# MMD_data_no_ny <- M3D_Wrapper_lite(wgbs_no_ny, overlaps_no_ny)



save(bs_files,
     prom_regions,
#      colData_do_no,
#      wgbs_do_no,
#      overlaps_do_no,
#      MMD_data_do_no, 
     colData_dy_ny,
     wgbs_dy_ny,
     overlaps_dy_ny,
     MMD_data_dy_ny,
#      colData_dy_do,
#      wgbs_dy_do,
#      overlaps_dy_do,
#      MMD_data_dy_do,
#      colData_no_ny,
#      wgbs_no_ny,
#      overlaps_no_ny,
#      MMD_data_no_ny,
     file = paste0("../datasets/BEATSON/M3D_", downstream, "_DY_NY.RData"))

