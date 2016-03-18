preprocess_diff_data <-function(HTS_data, max_outl = 600, gene_expr_thresh = FALSE,
                           gene_outl_thresh = FALSE, gene_log2_transf = FALSE){
  
  X <- list(control   = HTS_data$contr$methyl_region,
              treatment = HTS_data$treat$methyl_region)
  Y   <- list(control   = HTS_data$contr$rna_data$gene_fpkm,
              treatment = HTS_data$treat$rna_data$gene_fpkm)
  genes <- list(control   = HTS_data$contr$rna_data,
                treatment = HTS_data$treat$rna_data)

  # Delete all unexpressed genes
  if (gene_expr_thresh){
    # Remove the zero expression genes from control
    ind <- which(Y$control == 0)
    
    Y$control <- Y$control[-ind]
    Y$treatment <- Y$treatment[-ind]
    
    X$control <- X$control[-ind]
    X$treatment <- X$treatment[-ind]
    
    genes$control <- genes$control[-ind]
    genes$treatment <- genes$treatment[-ind]
    
    
    # Remove the zero expression genes from treatment
    ind <- which(Y$treatment == 0)
    
    Y$control <- Y$control[-ind]
    Y$treatment <- Y$treatment[-ind]
    
    X$control <- X$control[-ind]
    X$treatment <- X$treatment[-ind]
    
    genes$control <- genes$control[-ind]
    genes$treatment <- genes$treatment[-ind]
  }

  # Delete possible outliers / noisy data
  if (gene_outl_thresh){
    # Outliers from control
    ind <- which(Y$control > max_outl)
    
    Y$control <- Y$control[-ind]
    Y$treatment <- Y$treatment[-ind]
    
    X$control <- X$control[-ind]
    X$treatment <- X$treatment[-ind]
    
    genes$control <- genes$control[-ind]
    genes$treatment <- genes$treatment[-ind]
    
    
    # Outliers from treatment
    ind <- which(Y$treatment > max_outl)
    
    Y$control <- Y$control[-ind]
    Y$treatment <- Y$treatment[-ind]
    
    X$control <- X$control[-ind]
    X$treatment <- X$treatment[-ind]
    
    genes$control <- genes$control[-ind]
    genes$treatment <- genes$treatment[-ind]
  }

  # Log transform the gene expression data
  if (gene_log2_transf){
    Y$control <- Y$control + 0.1
    Y$treatment <- Y$treatment + 0.1
    
    Y$control <- log2(Y$control)
    Y$treatment <- log2(Y$treatment)
  }
  return(list(Y = Y, X = X, genes = genes))
}
