# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(processHTS)
library(data.table)

# Enhancer file
enhancer_file <- "../datasets/BEATSON/super_enhancers.bed"
upstream <- 15000
downstream <- 15000

# # Column names of Enhancer file
# col_names <- c("chr", "start", "end", "strand", "gene_name", "ensembl_id")
# col_classes <- c("character", numeric(), numeric(), "character", "character", "character")
# 
# # Read RNA-Seq data and store them in a data.frame object
# enhancer_data <- read.table(file = enhancer_file, 
#                             header = FALSE, 
#                             sep = "\t", 
#                             col.names = col_names,
#                             colClasses = col_classes,
#                             strip.white = TRUE,
#                             comment.char = "",
#                             stringsAsFactors = FALSE)

enhancer_data <- data.table::fread(input = enhancer_file,
                              sep = "\t",
                              header = TRUE,
                              col.names = c("chr", "start", "end", "strand",
                                            "gene_name", "ensembl_id"))

hist(as.numeric(enhancer_data$end) - as.numeric(enhancer_data$start))

dupl <- which(duplicated((enhancer_data$gene_name)))
enhancers <- enhancer_data[-dupl,]
dupl <- which(duplicated(enhancers$start))
enhancers <- enhancers[-dupl,]
rownames(enhancers) <- NULL


annot_file <- "../datasets/BEATSON/Ensembl_Known_Coding_67_Promoters_Annotated.bed"
# Read gene annotation file
annot_data <- read_annot_beatson(file = annot_file, 
                                 chr_discarded = NULL, 
                                 is_GRanges = FALSE)
# Remove duplicate gene names
annot_data <- annot_data[-which(duplicated(annot_data$gene_name)),]

# Keep only genes that overlap with enhancers
gene_annot <- annot_data[annot_data$gene_name %in% enhancers$gene_name]
rownames(gene_annot) <- NULL

# Get strand information in order to create the appropriate methylation profiles
enhancers$strand <- gene_annot$strand

# Create GRanges object
message("Creating GRanges object for enhancer data ...")
enh_data <- GenomicRanges::GRanges(seqnames = enhancers$chr,
                           ranges   = IRanges::IRanges(start = enhancers$start,
                                                       end   = enhancers$end),
                           strand   = enhancers$strand,
                           ensembl_id = enhancers$ensembl_id,
                           gene_name = enhancers$gene_name)
