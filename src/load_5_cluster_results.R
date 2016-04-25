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
load("../files/cluster_15_df_old_MonFeb221823.RData")
df_old_model <- bpr_model
df_old_basis <- basis
df_old_HTS_clust <- HTS_data
K <- df_old_model$K

proc_data <- preprocess_data(HTS_data = df_old_HTS_clust,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
df_old_obs <- proc_data$obs
df_old_Y <- proc_data$Y
df_old_genes <- proc_data$genes

# Cluster labels for DF Old mice
df_old_labels <- list()
df_old_expr <- list()
df_old_gene_ids <- list()
for (i in 1:K){
  df_old_labels[[i]] <- which(df_old_model$labels == i)
  df_old_expr[[i]] <- df_old_Y[df_old_labels[[i]]]
  df_old_gene_ids[[i]] <- df_old_genes$ensembl_id[df_old_labels[[i]]]
  #write(df_old_gene_ids[[i]], paste0("../results/df_old_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}


##-----------------DF YOUNG mice results ------------
load("../files/cluster_15_df_young_MonFeb221817.RData")
df_young_model <- bpr_model
df_young_basis <- basis
df_young_HTS_clust <- HTS_data

proc_data <- preprocess_data(HTS_data = df_young_HTS_clust,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
df_young_obs <- proc_data$obs
df_young_Y <- proc_data$Y
df_young_genes <- proc_data$genes

# Cluster labels for DF Young mice
df_young_labels <- list()
df_young_expr <- list()
df_young_gene_ids <- list()
for (i in 1:K){
  df_young_labels[[i]] <- which(df_young_model$labels == i)
  df_young_expr[[i]] <- df_young_Y[df_young_labels[[i]]]
  df_young_gene_ids[[i]] <- df_young_genes$ensembl_id[df_young_labels[[i]]]
  #write(df_young_gene_ids[[i]], paste0("../results/df_young_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}




##------------- Normal OLD mice results ------------
load("../files/cluster_15_N_old_MonFeb221829.RData")
N_old_model <- bpr_model
N_old_basis <- basis
N_old_HTS_clust <- HTS_data

proc_data <- preprocess_data(HTS_data = N_old_HTS_clust,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
N_old_obs <- proc_data$obs
N_old_Y <- proc_data$Y
N_old_genes <- proc_data$genes

# Cluster labels for DF Old mice
N_old_labels <- list()
N_old_expr <- list()
N_old_gene_ids <- list()
for (i in 1:K){
  N_old_labels[[i]] <- which(N_old_model$labels == i)
  N_old_expr[[i]] <- N_old_Y[N_old_labels[[i]]]
  N_old_gene_ids[[i]] <- N_old_genes$ensembl_id[N_old_labels[[i]]]
  #write(N_old_gene_ids[[i]], paste0("../results/N_old_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}




##------------ Normal YOUNG mice results ----------
load("../files/cluster_15_N_young_MonFeb221826.RData")
N_young_model <- bpr_model
N_young_basis <- basis
N_young_HTS_clust <- HTS_data

proc_data <- preprocess_data(HTS_data = N_young_HTS_clust,
                             max_outl = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)
N_young_obs <- proc_data$obs
N_young_Y <- proc_data$Y
N_young_genes <- proc_data$genes

# Cluster labels for N Young mice
N_young_labels <- list()
N_young_expr <- list()
N_young_gene_ids <- list()
for (i in 1:K){
  N_young_labels[[i]] <- which(N_young_model$labels == i)
  N_young_expr[[i]] <- N_young_Y[N_young_labels[[i]]]
  N_young_gene_ids[[i]] <- N_young_genes$ensembl_id[N_young_labels[[i]]]
  #write(N_young_gene_ids[[i]], paste0("../results/N_young_", K, "_clust_", i, "_", format(Sys.time(), "%a%b%d%H%M"), ".txt"))
}
