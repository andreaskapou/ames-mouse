# -------------------------------------------------
# Script for initializing experimental parameters
# for the ames-mouse dataset when performing regression.
# -------------------------------------------------


# -------------------------------------------------
# Initialize parameters for processing HTS data
# -------------------------------------------------
upstream    <- -5000
downstream  <- 5000
cpg_density <- 12
sd_thresh   <- 7e-02
min_bs_cov  <- 6
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
em_max_iter     <- 20
epsilon_conv    <- 1e-4
opt_method      <- "CG"
opt_itnmax      <- 50
init_opt_itnmax <- 100
is_parallel     <- TRUE
no_cores        <- 3
is_verbose      <- TRUE


basis <- rbf.object(M = 3)