# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
library(ggplot2)
library(cowplot)
R.utils::sourceDirectory("lib", modifiedOnly=FALSE)


# -----------------------------------------
# Initialize variables
# -----------------------------------------
do_file <- "../files/cluster_df_old_5000_ThuMay121242.RData"
dy_file <- "../files/cluster_df_young_5000_ThuMay121237.RData"
no_file <- "../files/cluster_N_old_5000_ThuMay121247.RData"
ny_file <- "../files/cluster_N_young_5000_ThuMay121244.RData"

# -----------------------------------------
# Load saved data for DO
# -----------------------------------------
load(do_file)
do_HTS_data <- HTS_data
do_proc_data <- proc_data
do_mix_model <- bpr_model

K <- do_mix_model$K

# Cluster labels for DO cell line
do_labels <- list()
do_expr <- list()
do_gene_ids <- list()

do_clust_ids <- c(3, 5, 4, 1, 2)
for (i in 1:K){
  do_labels[[i]] <- which(do_mix_model$labels == do_clust_ids[i])
  #print(length(do_labels[[i]]))
  do_expr[[i]] <- do_proc_data$Y[do_labels[[i]]]
  do_gene_ids[[i]] <- do_proc_data$genes$ensembl_id[do_labels[[i]]]
  #write(do_gene_ids[[i]], paste0("../results/final_do_5000_", K, "_4_clust_", i, ".txt"))
}



# -----------------------------------------
# Load saved data for DY
# -----------------------------------------
load(dy_file)
dy_HTS_data <- HTS_data
dy_proc_data <- proc_data
dy_mix_model <- bpr_model

# Cluster labels for DY cell line
dy_labels <- list()
dy_expr <- list()
dy_gene_ids <- list()

dy_clust_ids <- c(5, 1, 3, 4, 2)

for (i in 1:K){
  dy_labels[[i]] <- which(dy_mix_model$labels == dy_clust_ids[i])
  #print(length(dy_labels[[i]]))
  dy_expr[[i]] <- dy_proc_data$Y[dy_labels[[i]]]
  dy_gene_ids[[i]] <- dy_proc_data$genes$ensembl_id[dy_labels[[i]]]
  #write(dy_gene_ids[[i]], paste0("../results/final_dy_5000_", K, "_4_clust_", i, ".txt"))
}


# -----------------------------------------
# Load saved data for NO
# -----------------------------------------
load(no_file)
no_HTS_data <- HTS_data
no_proc_data <- proc_data
no_mix_model <- bpr_model

# Cluster labels for NO cell line
no_labels <- list()
no_expr <- list()
no_gene_ids <- list()

no_clust_ids <- c(4, 5, 2, 1, 3)

for (i in 1:K){
  no_labels[[i]] <- which(no_mix_model$labels == no_clust_ids[i])
  #print(length(no_labels[[i]]))
  no_expr[[i]] <- no_proc_data$Y[no_labels[[i]]]
  no_gene_ids[[i]] <- no_proc_data$genes$ensembl_id[no_labels[[i]]]
  #write(no_gene_ids[[i]], paste0("../results/final_no_5000_", K, "_4_clust_", i, ".txt"))
}


# -----------------------------------------
# Load saved data for NY
# -----------------------------------------
load(ny_file)
ny_HTS_data <- HTS_data
ny_proc_data <- proc_data
ny_mix_model <- bpr_model

# Cluster labels for NY cell line
ny_labels <- list()
ny_expr <- list()
ny_gene_ids <- list()

ny_clust_ids <- c(3, 4, 1, 5, 2)

for (i in 1:K){
  ny_labels[[i]] <- which(ny_mix_model$labels == ny_clust_ids[i])
  #print(length(ny_labels[[i]]))
  ny_expr[[i]] <- ny_proc_data$Y[ny_labels[[i]]]
  ny_gene_ids[[i]] <- ny_proc_data$genes$ensembl_id[ny_labels[[i]]]
  #write(ny_gene_ids[[i]], paste0("../results/final_ny_5000_", K, "_4_clust_", i, ".txt"))
}


total_do <- vector(mode = "numeric")
total_dy <- vector(mode = "numeric")
total_no <- vector(mode = "numeric")
total_ny <- vector(mode = "numeric")
for (k in 1:K){
  total_do <- c(total_do, do_gene_ids[[k]])
  total_dy <- c(total_dy, dy_gene_ids[[k]])
  total_no <- c(total_no, no_gene_ids[[k]])
  total_ny <- c(total_ny, ny_gene_ids[[k]])
}

inters_cell <- Reduce(intersect, list(total_do, total_dy, total_no, total_ny))

comm_do_ids <- list()
comm_dy_ids <- list()
comm_no_ids <- list()
comm_ny_ids <- list()
for (k in 1:K){
  comm_do_ids[[k]] <- Reduce(intersect, list(inters_cell, do_gene_ids[[k]]))
  #write(comm_do_ids[[k]], paste0("../results/final_comm_do_5000_", K, "_4_clust_", k, ".txt"))
  comm_dy_ids[[k]] <- Reduce(intersect, list(inters_cell, dy_gene_ids[[k]]))
  #write(comm_dy_ids[[k]], paste0("../results/final_comm_dy_5000_", K, "_4_clust_", k, ".txt"))
  comm_no_ids[[k]] <- Reduce(intersect, list(inters_cell, no_gene_ids[[k]]))
  #write(comm_no_ids[[k]], paste0("../results/final_comm_no_5000_", K, "_4_clust_", k, ".txt"))
  comm_ny_ids[[k]] <- Reduce(intersect, list(inters_cell, ny_gene_ids[[k]]))
  #write(comm_ny_ids[[k]], paste0("../results/final_comm_ny_5000_", K, "_4_clust_", k, ".txt"))
}



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


# -------------------------------------
# Preprocess data for plotting
# -------------------------------------
merged_meth <- list(do_mix_model, dy_mix_model, no_mix_model, ny_mix_model)
merged_expr <- list(do_expr, dy_expr, no_expr, ny_expr)
cell_lines <- c("DO", "DY", "NO", "NY")
x_len   <- 2000

mat_clust_ids <- matrix(cbind(do_clust_ids, dy_clust_ids, no_clust_ids, ny_clust_ids), ncol = 4)

# ----------------------------------------------
# Create data frame containing all experimental 
# output for clustered methylation profiles
# ----------------------------------------------
df_meth <- data.frame(xs = numeric(),
                      y = numeric(), 
                      cluster = character(), 
                      cell_line = character(),
                      stringsAsFactors=TRUE)
for (i in 1:length(merged_meth)){
  for (k in 1:K){
    xs <- seq(-1, 1,len = x_len) # create some values
    y <- as.vector(eval_probit_function(merged_meth[[i]]$basis, 
                                        xs, 
                                        merged_meth[[i]]$w[, mat_clust_ids[k,i]]))
    cluster <- paste("Cluster", k)
    cell_line <- cell_lines[i]
    dd <- data.frame(xs, y, cluster, cell_line, stringsAsFactors = TRUE)
    df_meth <- data.frame(rbind(df_meth, dd))
  }
}


# ----------------------------------------------
# Create data frame containing all experimental 
# output for  clustered gene expression levels
# ----------------------------------------------
df_expr <- data.frame(expr = numeric(), 
                      cluster = character(), 
                      cell_line = character(),
                      stringsAsFactors=TRUE)

for (i in 1:length(merged_expr)){
  for (k in 1:K){
    expr <- merged_expr[[i]][[k]]
    cluster <- paste("Cluster", k)
    cell_line <- cell_lines[i]
    dd <- data.frame(expr, cluster, cell_line, stringsAsFactors = TRUE)
    df_expr <- data.frame(rbind(df_expr, dd))
  }
}


# ----------------------------------------------
# Create plots
# ----------------------------------------------
gg_prof <- ggplot_cluster_prof(df = df_meth, 
                               main_lab = "Clustered methylation profiles",
                               label = c("-5kb", "TSS", "5kb"))
gg_expr <- ggplot_cluster_expr2(df = df_expr, main_lab = "Clustered expression levels")

cluster_plot <- plot_grid(gg_prof, gg_expr, labels = c("A", "B"), 
                          label_size = 20, ncol = 1, nrow = 2)
