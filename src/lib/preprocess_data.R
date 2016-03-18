preprocess_data <-function(HTS_data, max_outl = 600, gene_expr_thresh = FALSE,
                           gene_outl_thresh = FALSE, gene_log2_transf = FALSE){
  X <- HTS_data$methyl_region
  Y   <- HTS_data$rna_data$gene_fpkm
  genes <- HTS_data$rna_data
  
  # Delete all unexpressed genes
  if (gene_expr_thresh){
    # Remove the zero expression genes
    ind <- which(Y == 0)
    Y <- Y[-ind]
    X <- X[-ind]
    genes <- genes[-ind]
  }
  
  # Delete possible outliers / noisy data
  if (gene_outl_thresh){
    ind <- which(Y > max_outl)
    Y <- Y[-ind]
    X <- X[-ind]
    genes <- genes[-ind]
  }
  
  # Log transform the gene expression data
  if (gene_log2_transf){
    Y <- Y + 0.1
    Y <- log2(Y)
  }
  return(list(Y = Y, X = X, genes = genes))
}
