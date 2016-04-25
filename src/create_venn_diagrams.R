# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(mpgex)
library(processHTS)
library(VennDiagram)
R.utils::sourceDirectory("lib", modifiedOnly = FALSE) 

#----------------
# Create Venn diagrams for K = 5, RBF = 3
#----------------

# Cluster 1 --------------
df_old_1   <- "../results/df_old_3_rbf5_clust_1_TueFeb230953.txt"
df_young_1 <- "../results/df_young_3_rbf5_clust_1_TueFeb230953.txt"
N_old_1    <- "../results/N_old_3_rbf5_clust_2_TueFeb230953.txt"
N_young_1  <- "../results/N_young_3_rbf5_clust_2_TueFeb230953.txt"

create_quad_venn(df_old_1, df_young_1, N_old_1, N_young_1,
                 filename = "../figures/Cluster_1_k_5_rbf_3.png")


# Cluster 2 --------------------
df_old_2   <- "../results/df_old_3_rbf5_clust_2_TueFeb230953.txt"
df_young_2 <- "../results/df_young_3_rbf5_clust_2_TueFeb230953.txt"
N_old_2    <- "../results/N_old_3_rbf5_clust_4_TueFeb230953.txt"
N_young_2  <- "../results/N_young_3_rbf5_clust_1_TueFeb230953.txt"

create_quad_venn(df_old_2, df_young_2, N_old_2, N_young_2,
                 filename = "../figures/Cluster_2_k_5_rbf_3.png")



# Cluster 3 -------------------
df_old_3   <- "../results/df_old_3_rbf5_clust_3_TueFeb230953.txt"
df_young_3 <- "../results/df_young_3_rbf5_clust_5_TueFeb230953.txt"
N_old_3    <- "../results/N_old_3_rbf5_clust_1_TueFeb230953.txt"
N_young_3  <- "../results/N_young_3_rbf5_clust_3_TueFeb230953.txt"

create_quad_venn(df_old_3, df_young_3, N_old_3, N_young_3,
                 filename = "../figures/Cluster_3_k_5_rbf_3.png")


# Cluster 4 -------------------
df_old_4   <- "../results/df_old_3_rbf5_clust_4_TueFeb230953.txt"
df_young_4 <- "../results/df_young_3_rbf5_clust_4_TueFeb230953.txt"
N_old_4    <- "../results/N_old_3_rbf5_clust_3_TueFeb230953.txt"
N_young_4  <- "../results/N_young_3_rbf5_clust_4_TueFeb230953.txt"

create_quad_venn(df_old_4, df_young_4, N_old_4, N_young_4,
                 filename = "../figures/Cluster_4_k_5_rbf_3.png")


# Cluster 5 -------------------
df_old_5   <- "../results/df_old_3_rbf5_clust_5_TueFeb230953.txt"
df_young_5 <- "../results/df_young_3_rbf5_clust_3_TueFeb230953.txt"
N_old_5    <- "../results/N_old_3_rbf5_clust_5_TueFeb230953.txt"
N_young_5  <- "../results/N_young_3_rbf5_clust_5_TueFeb230953.txt"

create_quad_venn(df_old_5, df_young_5, N_old_5, N_young_5,
                 filename = "../figures/Cluster_5_k_5_rbf_3.png")




#----------------
# Create Venn diagrams for K = 7, RBF = 5
#----------------

# Cluster 1 --------------
df_old_1   <- "../results/df_old_7_clust_1_WedFeb240033.txt"
df_young_1 <- "../results/df_young_7_clust_4_WedFeb240033.txt"
N_old_1    <- "../results/N_old_7_clust_1_WedFeb240033.txt"
N_young_1  <- "../results/N_young_7_clust_3_WedFeb240033.txt"

create_quad_venn(df_old_1, df_young_1, N_old_1, N_young_1,
                 filename = "../figures/Cluster_1_k_7_rbf_5.png")

# Cluster 2 --------------
df_old_2   <- "../results/df_old_7_clust_2_WedFeb240033.txt"
df_young_2 <- "../results/df_young_7_clust_6_WedFeb240033.txt"
N_old_2    <- "../results/N_old_7_clust_6_WedFeb240033.txt"
N_young_2  <- "../results/N_young_7_clust_5_WedFeb240033.txt"

create_quad_venn(df_old_2, df_young_2, N_old_2, N_young_2,
                 filename = "../figures/Cluster_2_k_7_rbf_5.png")


# Cluster 3 --------------
df_old_3   <- "../results/df_old_7_clust_3_WedFeb240033.txt"
df_young_3 <- "../results/df_young_7_clust_1_WedFeb240033.txt"
N_old_3    <- "../results/N_old_7_clust_3_WedFeb240033.txt"
N_young_3  <- "../results/N_young_7_clust_1_WedFeb240033.txt"

create_quad_venn(df_old_3, df_young_3, N_old_3, N_young_3,
                 filename = "../figures/Cluster_3_k_7_rbf_5.png")


# Cluster 4 --------------
df_old_4   <- "../results/df_old_7_clust_4_WedFeb240033.txt"
df_young_4 <- "../results/df_young_7_clust_7_WedFeb240033.txt"
N_old_4    <- "../results/N_old_7_clust_4_WedFeb240033.txt"
N_young_4  <- "../results/N_young_7_clust_7_WedFeb240033.txt"

create_quad_venn(df_old_4, df_young_4, N_old_4, N_young_4,
                 filename = "../figures/Cluster_4_k_7_rbf_5.png")


# Cluster 5 --------------
df_old_5   <- "../results/df_old_7_clust_5_WedFeb240033.txt"
df_young_5 <- "../results/df_young_7_clust_2_WedFeb240033.txt"
N_old_5    <- "../results/N_old_7_clust_5_WedFeb240033.txt"
N_young_5  <- "../results/N_young_7_clust_2_WedFeb240033.txt"

create_quad_venn(df_old_5, df_young_5, N_old_5, N_young_5,
                 filename = "../figures/Cluster_5_k_7_rbf_5.png")


# Cluster 6 --------------
df_old_6   <- "../results/df_old_7_clust_6_WedFeb240033.txt"
df_young_6 <- "../results/df_young_7_clust_5_WedFeb240033.txt"
N_old_6    <- "../results/N_old_7_clust_7_WedFeb240033.txt"
N_young_6  <- "../results/N_young_7_clust_4_WedFeb240033.txt"

create_quad_venn(df_old_6, df_young_6, N_old_6, N_young_6,
                 filename = "../figures/Cluster_6_k_7_rbf_5.png")


# Cluster 7 --------------
df_old_7   <- "../results/df_old_7_clust_7_WedFeb240033.txt"
df_young_7 <- "../results/df_young_7_clust_3_WedFeb240033.txt"
N_old_7    <- "../results/N_old_7_clust_2_WedFeb240033.txt"
N_young_7  <- "../results/N_young_7_clust_6_WedFeb240033.txt"

create_quad_venn(df_old_7, df_young_7, N_old_7, N_young_7,
                 filename = "../figures/Cluster_7_k_7_rbf_5.png")
