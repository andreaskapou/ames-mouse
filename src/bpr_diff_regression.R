# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
HTS_data <- process_diff_beatson(bs_contr_files  = bs_contr_files,
                                 bs_treat_files  = bs_treat_files,
                                 rna_contr_files = rna_contr_files,
                                 rna_treat_files = rna_treat_files,
                                 chr_discarded   = chr_discarded,
                                 upstream        = upstream,
                                 downstream      = downstream,
                                 cpg_density     = cpg_density,
                                 sd_thresh       = sd_thresh,
                                 min_bs_cov      = min_bs_cov,
                                 ignore_strand   = ignore_strand)

proc_data <- preprocess_diff_data(HTS_data         = HTS_data,
                                  max_outl         = max_outl,
                                  gene_expr_thresh = gene_expr_thresh,
                                  gene_outl_thresh = gene_outl_thresh,
                                  gene_log2_transf = gene_log2_transf)


# ------------------------------------------
# Apply methylation profiles
# ------------------------------------------
set.seed(seed)
out_prof <- mpgex_differential_regr(formula       = formula,
                                    x             = proc_data$X,
                                    y             = proc_data$Y,
                                    model_name    = model_name,
                                    basis         = basis_prof,
                                    train_perc    = train_perc,
                                    opt_method    = opt_method,
                                    opt_itnmax    = opt_itnmax,
                                    is_parallel   = is_parallel,
                                    no_cores      = no_cores,
                                    is_summary    = is_summary)


# ------------------------------------------
# Apply mean methylation
# ------------------------------------------
set.seed(seed)
out_mean <- mpgex_differential_regr(formula       = formula,
                                    x             = proc_data$X,
                                    y             = proc_data$Y,
                                    model_name    = model_name,
                                    basis         = basis_mean,
                                    train_perc    = train_perc,
                                    opt_method    = opt_method,
                                    opt_itnmax    = opt_itnmax,
                                    is_parallel   = is_parallel,
                                    no_cores      = no_cores,
                                    is_summary    = is_summary)
