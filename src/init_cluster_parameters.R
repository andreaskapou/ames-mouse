# -------------------------------------------------
# Script for initializing experimental parameters
# for the ames-mouse dataset when performing regression.
# -------------------------------------------------


# -------------------------------------------------
# Initialize parameters for processing HTS data
# -------------------------------------------------
upstream    <- -7000
downstream  <- 7000
cpg_density <- 15
sd_thresh   <- 10e-02
min_bs_cov  <- 4
chr_discarded <- c("chrLambda", "chrM")
ignore_strand <- TRUE

gene_expr_thresh <- FALSE
gene_outl_thresh <- TRUE
gene_log2_transf <- TRUE
max_outl <- 600


# -------------------------------------------------
# Initialize parameters for EM BPR model
# -------------------------------------------------
seed            <- 1234
K               <- 5
w               <- NULL
pi_k            <- NULL
em_max_iter     <- 30
epsilon_conv    <- 1e-4
opt_method      <- "CG"
opt_itnmax      <- 50
init_opt_itnmax <- 70
is_parallel     <- TRUE
no_cores        <- 5
is_verbose      <- TRUE


basis <- rbf.object(M = 4)
