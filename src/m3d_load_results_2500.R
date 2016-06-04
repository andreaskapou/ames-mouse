# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(M3Ddevel)
library(BiSeq)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE)


## ---------------------------------------------------------------

load("../datasets/BEATSON/Diff_Meth/M3D_2500_DO_NO.RData")
load("../datasets/BEATSON/Diff_Meth/genes_2500.RData")

group1 <- unique(colData(wgbs_do_no)$group)[1]
group2 <-unique(colData(wgbs_do_no)$group)[2]
p_values <- pvals_lite(rrbs = wgbs_do_no, 
                       CpGs = prom_regions, 
                       M3D_stat_lite = MMD_data_do_no,
                       group1 = group1,
                       group2 = group2)

called <- which(p_values$FDRmean <= 0.05)

diff_genes_do_no <- gene_data[called]

plotMethProfile(wgbs_do_no, prom_regions, 'DF_Old', 'N_Old', called[1])
rm(wgbs_do_no)
rm(overlaps_do_no)
rm(MMD_data_do_no)
rm(p_values)


## ---------------------------------------------------------------

load("../datasets/BEATSON/Diff_Meth/M3D_2500_DY_DO.RData")
load("../datasets/BEATSON/Diff_Meth/genes_2500.RData")

group1 <- unique(colData(wgbs_dy_do)$group)[1]
group2 <-unique(colData(wgbs_dy_do)$group)[2]
p_values <- pvals_lite(rrbs = wgbs_dy_do, 
                       CpGs = prom_regions, 
                       M3D_stat_lite = MMD_data_dy_do,
                       group1 = group1,
                       group2 = group2)

called <- which(p_values$FDRmean <= 0.05)
diff_genes_dy_do <- gene_data[called]

plotMethProfile(wgbs_dy_do, prom_regions, 'DF_Young', 'DF_Old', called[1])
rm(wgbs_dy_do)
rm(overlaps_dy_do)
rm(MMD_data_dy_do)
rm(p_values)


## ---------------------------------------------------------------

load("../datasets/BEATSON/Diff_Meth/M3D_2500_DY_NY.RData")
load("../datasets/BEATSON/Diff_Meth/genes_2500.RData")

group1 <- unique(colData(wgbs_dy_ny)$group)[1]
group2 <-unique(colData(wgbs_dy_ny)$group)[2]
p_values <- pvals_lite(rrbs = wgbs_dy_ny, 
                       CpGs = prom_regions, 
                       M3D_stat_lite = MMD_data_dy_ny,
                       group1 = group1,
                       group2 = group2)

called <- which(p_values$FDRmean <= 0.05)
diff_genes_dy_ny <- gene_data[called]

plotMethProfile(wgbs_dy_ny, prom_regions, 'DF_Young', 'N_Young', called[1])
rm(wgbs_dy_ny)
rm(overlaps_dy_ny)
rm(MMD_data_dy_ny)
rm(p_values)


## ---------------------------------------------------------------

load("../datasets/BEATSON/Diff_Meth/M3D_2500_NO_NY.RData")
load("../datasets/BEATSON/Diff_Meth/genes_2500.RData")

group1 <- unique(colData(wgbs_no_ny)$group)[1]
group2 <-unique(colData(wgbs_no_ny)$group)[2]
p_values <- pvals_lite(rrbs = wgbs_no_ny, 
                       CpGs = prom_regions, 
                       M3D_stat_lite = MMD_data_no_ny,
                       group1 = group1,
                       group2 = group2)

called <- which(p_values$FDRmean <= 0.05)
diff_genes_no_ny <- gene_data[called]

#plotMethProfile(wgbs_no_ny, prom_regions, 'N_Old', 'N_Young', called[1])
rm(wgbs_no_ny)
rm(overlaps_no_ny)
rm(MMD_data_no_ny)
rm(p_values)

save(diff_genes_do_no, diff_genes_dy_do, diff_genes_dy_ny, diff_genes_no_ny,
     file = "../files/diff_genes_2500.RData")