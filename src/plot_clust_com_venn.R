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
library(VennDiagram)
R.utils::sourceDirectory("lib", modifiedOnly=FALSE)
#----------------
# Create Venn diagrams for K = 5, RBF = 4
#----------------

# Cluster 1 --------------
do_1 <- "../results/final_comm_do_7000_5_4_clust_1.txt"
dy_1 <- "../results/final_comm_dy_7000_5_4_clust_1.txt"
no_1 <- "../results/final_comm_no_7000_5_4_clust_1.txt"
ny_1 <- "../results/final_comm_ny_7000_5_4_clust_1.txt"

create_quad_venn(do_1, dy_1, no_1, ny_1,
                   filename = "../figures/final_comm_cluster_1_7000_5_4.png")


# Cluster 2 --------------------
do_2 <- "../results/final_comm_do_7000_5_4_clust_2.txt"
dy_2 <- "../results/final_comm_dy_7000_5_4_clust_2.txt"
no_2 <- "../results/final_comm_no_7000_5_4_clust_2.txt"
ny_2 <- "../results/final_comm_ny_7000_5_4_clust_2.txt"

create_quad_venn(do_2, dy_2, no_2, ny_2,
                 filename = "../figures/final_comm_cluster_2_7000_5_4.png")



# Cluster 3 -------------------
do_3 <- "../results/final_comm_do_7000_5_4_clust_3.txt"
dy_3 <- "../results/final_comm_dy_7000_5_4_clust_3.txt"
no_3 <- "../results/final_comm_no_7000_5_4_clust_3.txt"
ny_3 <- "../results/final_comm_ny_7000_5_4_clust_3.txt"

create_quad_venn(do_3, dy_3, no_3, ny_3,
                 filename = "../figures/final_comm_cluster_3_7000_5_4.png")


# Cluster 4 -------------------
do_4 <- "../results/final_comm_do_7000_5_4_clust_4.txt"
dy_4 <- "../results/final_comm_dy_7000_5_4_clust_4.txt"
no_4 <- "../results/final_comm_no_7000_5_4_clust_4.txt"
ny_4 <- "../results/final_comm_ny_7000_5_4_clust_4.txt"

create_quad_venn(do_4, dy_4, no_4, ny_4,
                 filename = "../figures/final_comm_cluster_4_7000_5_4.png")


# Cluster 5 -------------------
do_5 <- "../results/final_comm_do_7000_5_4_clust_5.txt"
dy_5 <- "../results/final_comm_dy_7000_5_4_clust_5.txt"
no_5 <- "../results/final_comm_no_7000_5_4_clust_5.txt"
ny_5 <- "../results/final_comm_ny_7000_5_4_clust_5.txt"

create_quad_venn(do_5, dy_5, no_5, ny_5,
                 filename = "../figures/final_comm_cluster_5_7000_5_4.png")
