# ------------------------------------------
# Set working directory and load libraries
# ------------------------------------------
if (interactive()){
  cur.dir <- dirname(parent.frame(2)$ofile)
  setwd(cur.dir)
}
library(processHTS)

# RNA-Seq and gene annotation files
rna_file <- "../datasets/BEATSON/RNA-Seq/raw/FPKM_Matrix_genes.csv"
annot_file <- "../datasets/BEATSON/Ensembl_Known_Coding_67_Promoters_Annotated.bed"

# Read gene annotation file
annot_data <- read_annot_beatson(file = annot_file, 
                                 chr_discarded = NULL, 
                                 is_GRanges = FALSE)

# Column names of RNA-Seq file
col_names <- c("ensembl_id", "gene_name", "9L_df_young", "10L_df_young", 
               "34L_df_young", "35L_df_young", "1L_df_old", "2L_df_old", 
               "3L_df_old", "4L_df_old", "13L_N_young", "14L_N_young", 
               "15L_N_young", "16L_N_young", "5L_N_old", "6L_N_old", 
               "7L_N_old", "8L_N_old")

col_classes <- c("character", "character", numeric(), numeric(), numeric(), 
                 numeric(), numeric(), numeric(), numeric(), numeric(), 
                 numeric(), numeric(), numeric(), numeric(), numeric(), 
                 numeric(), numeric(), numeric())

# Read RNA-Seq data and store them in a data.frame object
rna_data <- read.table(file = rna_file, 
                       header = FALSE, 
                       sep = "\t", 
                       col.names = col_names,
                       colClasses = col_classes,
                       strip.white = TRUE,
                       comment.char = "",
                       stringsAsFactors = FALSE)

# Discard the gene name column
rna_data$gene_name = NULL

# Merge gene annotation data and RNA-Seq data by Ensembl ID
merged_data = merge(as.data.frame(annot_data), rna_data, by = "ensembl_id")

# Sorting data -----------------------------------------------
# With order priority: 1. chr, 2. start, 3. strand
message("Sorting RNA-Seq data ...")
merged_data <- merged_data[with(merged_data, order(merged_data$chr,
                                                   merged_data$start,
                                                   merged_data$strand)), ]


for (i in 1:16){
  filename = paste0("../datasets/BEATSON/RNA-Seq/processed/RNA_Seq_FPKM_", 
                    colnames(merged_data)[i + 6], ".bed")
  data = cbind(chr=merged_data[2], 
               start=merged_data[5], 
               end=merged_data[6], 
               strand=merged_data[4], 
               gene_name=merged_data[3],
               ensembl_id=merged_data[1], 
               score=merged_data[,i+6]
               )
  write.table(x = data, 
              file = filename, 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)
}