#!/usr/bin/Rscript

#Load required libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(apeglm))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(optparse))

#parse arguments
 
option_list = list(
    make_option(c("-c", "--counts"), type="character", default=NULL, 
              help="counts table input"),
    make_option(c("-n", "--num_reps"), type="integer", default=NULL, 
              help="number of replicates"),
    make_option(c("-p", "--cores"), type="integer", default=NULL, 
              help="number of cores to use"),
    make_option(c("-q", "--FDR"), type="numeric", default=NULL, 
              help="false-discovery rate"),
    make_option(c("-f", "--count_filter"), type="numeric", default=1, 
              help="number of counts to filter against"),
    make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="out directory")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#assign options
input_file <- opt$counts
rep_num <- opt$num_reps
cores <- opt$cores
FDR <- opt$FDR
cf <- opt$count_filter
out_dir <- opt$out_dir

#set up parallelization. 
register(MulticoreParam(cores))

#Assign sample info and set col-names of matrix based on replicate number argument. 
if (rep_num == 2) {
    cts_col_names <- c("Bin_ID", "Chr", "Start", "End", "Strand", "Length", "DNA1", "DNA2", "RNA1", "RNA2")
    condition <- c("DNA", "DNA", "RNA", "RNA")
    RNames <- c("DNA1", "DNA2", "RNA1", "RNA2")
}

if (rep_num == 3) {
    cts_col_names <- c("Bin_ID", "Chr", "Start", "End", "Strand", "Length", "DNA1", "DNA2", "DNA3", "RNA1", "RNA2", "RNA3")
    condition <- c("DNA", "DNA", "DNA", "RNA", "RNA", "RNA")
    RNames <- c("DNA1", "DNA2", "DNA3", "RNA1", "RNA2", "RNA3")
}

if (rep_num == 4) {
    cts_col_names <- c("Bin_ID", "Chr", "Start", "End", "Strand", "Length", "DNA1", "DNA2", "DNA3", "DNA4", "RNA1", "RNA2", "RNA3", "RNA4")
    condition <- c("DNA", "DNA", "DNA", "DNA", "RNA", "RNA", "RNA", "RNA")
    RNames <- c("DNA1", "DNA2", "DNA3", "DNA4", "RNA1", "RNA2", "RNA3", "RNA4")
}

if (rep_num == 5) {
    cts_col_names <- c("Bin_ID", "Chr", "Start", "End", "Strand", "Length", "DNA1", "DNA2", "DNA3", "DNA4", "DNA5", "RNA1", "RNA2", "RNA3", "RNA4", "RNA5")
    condition <- c("DNA", "DNA", "DNA", "DNA", "DNA", "RNA", "RNA", "RNA", "RNA", "RNA")
    RNames <- c( "DNA1", "DNA2", "DNA3", "DNA4", "DNA5", "RNA1", "RNA2", "RNA3", "RNA4", "RNA5")
}

#Read in counts matrix. Skip the first two rows, which is metadata. 
cts_df <- read_tsv(input_file, col_names = cts_col_names, skip = 2) 

#Convert to a matrix by selecting read count columns and assigning bin ID as rownames. In addition, apply count fiter. 
cts_matrix <- dplyr::select(cts_df, -Chr, -Start, -End, -Strand, -Length) %>% 
    column_to_rownames(var = "Bin_ID") %>% filter_if(is.numeric, all_vars((.) >= cf)) %>% 
    as.matrix()

#Prepare dataframe of sample info using the vectors defined in the if statements above. 
coldata <- data.frame(row.names = RNames, condition)

#Check that coldata matches cts. If false, stop. 
if (all(rownames(coldata) == colnames(cts_matrix)) == FALSE) {
    stop("Error: coldata does not match cts matrix")
}

#Set up experimental design for DEseq2
dds <- DESeqDataSetFromMatrix(countData = cts_matrix, colData = coldata, design = ~ condition)

#Differential expression analysis. Use local. I tried all three, mean has a poor fit and parametric defaults to local anyway. 
#This takes awhile so run across multiple cores. 
dds <- DESeq(dds, fitType="local", parallel=TRUE, BPPARAM=MulticoreParam(cores))

#Extract results but shrink L2FC using the apeglm model to remove low count confounders from analysis. 
#This also takes awhile so run across multiple cores. 
res <- results(dds, name="condition_RNA_vs_DNA", parallel=TRUE, BPPARAM=MulticoreParam(cores)) 

#Convert to df and add Bin_ID column from rownames. 
res_df <- as.data.frame(res) %>% rownames_to_column(var = "Bin_ID")

#Obtain active peaks using the supplied padj filter and l2fc > 0.
#In the same pipe, add positional information back by joining with cts_df. Then format to bed (log2Fc as score) and sort by position. Convert to 0-based format since the SAF input is 1-based. 
active <- filter(res_df, padj <= FDR) %>% filter(log2FoldChange > 0) %>% 
    left_join(cts_df, by = c("Bin_ID" = "Bin_ID"), keep = FALSE) %>% 
    dplyr::mutate(Start = Start - 1) %>%
    dplyr::select(Chr, Start, End, Bin_ID, log2FoldChange, Strand) %>% 
    dplyr::arrange(Chr, Start, End)

#repeat for silent, but only filter for less than 0 L2FC.
silent <- filter(res_df, padj <= FDR) %>% filter(log2FoldChange < 0) %>% 
    left_join(cts_df, by = c("Bin_ID" = "Bin_ID"), keep = FALSE) %>% 
    dplyr::mutate(Start = Start - 1) %>%
    dplyr::select(Chr, Start, End, Bin_ID, log2FoldChange, Strand) %>% 
    dplyr::arrange(Chr, Start, End)

#Save active and silent bed bin files. Also save r objects.
write_tsv(active, paste0(out_dir, "/active_bins.bed"), col_names = FALSE)
write_tsv(silent, paste0(out_dir, "/silent_bins.bed"), col_names = FALSE)
save(dds, file = paste0(out_dir, "/dds.Rdata"))
save(res, file = paste0(out_dir, "/res.Rdata"))
