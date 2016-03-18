# ------------------------------------------
# Read and preprocess HTS files
# ------------------------------------------
HTS_data <- process_beatson(bs_files        = bs_files,
                            rna_files       = rna_files,
                            chr_discarded   = chr_discarded,
                            upstream        = upstream,
                            downstream      = downstream,
                            cpg_density     = cpg_density,
                            sd_thresh       = sd_thresh,
                            min_bs_cov      = min_bs_cov,
                            ignore_strand   = ignore_strand)

proc_data <- preprocess_data(HTS_data         = HTS_data,
                             max_outl         = max_outl,
                             gene_expr_thresh = gene_expr_thresh,
                             gene_outl_thresh = gene_outl_thresh,
                             gene_log2_transf = gene_log2_transf)


# ------------------------------------------
# Apply EM algorithm on BPR model
# ------------------------------------------
set.seed(seed)
bpr_model <- mpgex_cluster(x               = proc_data$X,
                           K               = K,
                           pi_k            = pi_k,
                           w               = w,
                           basis           = basis,
                           em_max_iter     = em_max_iter,
                           epsilon_conv    = epsilon_conv,
                           opt_method      = opt_method,
                           opt_itnmax      = opt_itnmax,
                           init_opt_itnmax = init_opt_itnmax,
                           is_parallel     = is_parallel,
                           no_cores        = no_cores,
                           is_verbose      = is_verbose)
