# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(M3Ddevel)
library(BiSeq)
library(mpgex)
library(processHTS)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE)

load("../files/diff_genes_7000.RData")

if (length(diff_genes_do_no) > 0){
  # DO_NO
  t1_do_no = which(df_old_genes$gene_name == diff_genes_do_no$gene_name[1])
  t2_do_no = which(N_old_genes$gene_name == diff_genes_do_no$gene_name[1])
  gene_name_do_no <- df_old_genes$gene_name[t1_do_no]
  
  # Plot DMRs
  dmr_do_no <- plot_dmr_profiles(region1 = t1_do_no, region2 = t2_do_no, 
                                 X1 = df_old_X, X2 = N_old_X,
                                 prof1 = df_old_out_prof,
                                 prof2 = N_old_out_prof,
                                 title = paste0("DMR Gene ", gene_name_do_no),
                                 leg = c("DO", "NO"),
                                 label = c("-7kb", "TSS", "+7kb"))
}

if (length(diff_genes_dy_do) > 0){
  # DY_DO
  t1_dy_do = which(df_young_genes$gene_name == diff_genes_dy_do$gene_name[1])
  t2_dy_do = which(df_old_genes$gene_name == diff_genes_dy_do$gene_name[1])
  gene_name_dy_do <- df_young_genes$gene_name[t1_dy_do]
  
  # Plot DMRs
  dmr_dy_do <- plot_dmr_profiles(region1 = t1_dy_do, region2 = t2_dy_do, 
                                 X1 = df_young_X, X2 = df_old_X,
                                 prof1 = df_young_out_prof,
                                 prof2 = df_old_out_prof,
                                 title = paste0("DMR Gene ", gene_name_dy_do),
                                 leg = c("DY", "DO"),
                                 label = c("-7kb", "TSS", "+7kb"))
}


if (length(diff_genes_dy_ny) > 0){
  # DY_NY
  t1_dy_ny = which(df_young_genes$gene_name == diff_genes_dy_ny$gene_name[1])
  t2_dy_ny = which(N_young_genes$gene_name == diff_genes_dy_ny$gene_name[1])
  gene_name_dy_ny <- df_young_genes$gene_name[t1_dy_ny]
  
  # Plot DMRs
  dmr_dy_ny <- plot_dmr_profiles(region1 = t1_dy_ny, region2 = t2_dy_ny, 
                                 X1 = df_young_X, X2 = N_young_X,
                                 prof1 = df_young_out_prof,
                                 prof2 = N_young_out_prof,
                                 title = paste0("DMR Gene ", gene_name_dy_ny),
                                 leg = c("DY", "NY"),
                                 label = c("-7kb", "TSS", "+7kb"))
}


if (length(diff_genes_no_ny) > 0){
  # NO_NY
  t1_no_ny = which(N_old_genes$gene_name == diff_genes_no_ny$gene_name[1])
  if (length(t1_no_ny) != 0){
    t2_no_ny = which(N_young_genes$gene_name == diff_genes_no_ny$gene_name[1])
    gene_name_no_ny <- N_old_genes$gene_name[t1_no_ny]
    
    # Plot DMRs
    dmr_no_ny <- plot_dmr_profiles(region1 = t1_no_ny, region2 = t2_no_ny, 
                                   X1 = N_old_X, X2 = N_young_X,
                                   prof1 = N_old_out_prof,
                                   prof2 = N_young_out_prof,
                                   title = paste0("DMR Gene ", gene_name_no_ny),
                                   leg = c("NO", "NY"),
                                   label = c("-7kb", "TSS", "+7kb"))
  }
}
