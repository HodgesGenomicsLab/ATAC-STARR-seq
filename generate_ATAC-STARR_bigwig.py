"""
name    | Tyler Hansen
created | 2021.12.15
updated |

Description: Generates an ATAC-STARR activity bigWig file. This script consists of two major parts: 
    1) generate cpm normalized read count bigwig files for the RNA and DNA bams individually
    2) compare the RNA and DNA bigwig to create either a subtracted or log2 ratio comparison bigwig file. 

Input:
   DNA bam file with duplicates for the given replicate/condition
   RNA bam file with duplicates for the given replicate/condition

Output:
   4 bigWig files - DNA only, RNA only, and RNA/DNA comparisons either log2 or subtracted

Parameters:
-d (--DNA_bam): reisolated plasmid DNA bam file 
-r (--RNA_bam): reporter RNA bam file
-n (--name): name of the file prefix to be added to the bigwig files 
-o (--out_dir): output directory for all of the results 
-bs (--bin_size): binsize for the bigwig generation
-p (--num_threads): number of threads
-h (--help): display this help

Required Software downloaded via conda:
    Python v3.6.13
    deeptools v3.5.1  
"""

###
# software
import os
import pandas as pd
import sys
import subprocess
import argparse

###
#   arguments. Make all arguments required. 
#make arg_parser object
arg_parser = argparse.ArgumentParser(description="Generate ATAC-STARR bigWig")

#make a required arguments group
required = arg_parser.add_argument_group('required arguments')

#Add arguments 
required.add_argument("-d", "--DNA_bam", type=str,
                        help="reisolated plasmid DNA bam file")

required.add_argument("-r", "--RNA_bam", type=str,
                        help="reporter RNA bam file")

required.add_argument("-n", "--name", type=str,
                        help="name of file prefix for output files")

required.add_argument("-o", "--out_dir", type=str,
                        help='output directory')

arg_parser.add_argument("-bs", "--bin_size", type=int, default=50,
                        help='binsize parameter for the bigwig generation (default: 50bp)')

arg_parser.add_argument("-p", "--threads", type=int, default=1,
                        help='number of threads (default: 1)')

#order argument groups for help display 
groups_order = {
    'positional arguments': 0,
    'required arguments': 1,
    'optional arguments': 2
}

arg_parser._action_groups.sort(key=lambda g: groups_order[g.title])

#Assign arguments to arg object
args = arg_parser.parse_args()

# save parameters to variables
DNA_FILE = args.DNA_bam
RNA_FILE  = args.RNA_bam
NAME = args.name
OUT_DIR = args.out_dir
BIN_SIZE = args.bin_size
NUM_THREADS = args.threads

# print parameters for troubleshooting purposes
print("PARAMETERS:")
print("DNA File:", DNA_FILE)
print("RNA File:", RNA_FILE)
print("Output Prefix:", NAME)
print("Output Directory:", OUT_DIR)
print("Bin Size:", BIN_SIZE)
print("Using", NUM_THREADS, "cores for processing")
###

###
#generate DNA bigwig
cmd = f"bamCoverage -b {DNA_FILE} -of bigwig -bs {BIN_SIZE} -p {NUM_THREADS} --normalizeUsing CPM -o {OUT_DIR}/{NAME}_DNA-coverage.bw"
subprocess.call(cmd, shell = True)

#generate RNA bigwig
cmd = f"bamCoverage -b {RNA_FILE} -of bigwig -bs {BIN_SIZE} -p {NUM_THREADS} --normalizeUsing CPM -o {OUT_DIR}/{NAME}_RNA-coverage.bw"
subprocess.call(cmd, shell = True)

#compare bigwigs using log2 ratio operation. Apply pseudocount of 1 to avoid x/0 regions. 
cmd = f"bigwigCompare -b1 {OUT_DIR}/{NAME}_RNA-coverage.bw -b2 {OUT_DIR}/{NAME}_DNA-coverage.bw \
    --pseudocount 1 --operation log2 --skipZeroOverZero -p {NUM_THREADS} -of bigwig -o {OUT_DIR}/{NAME}_RNA-to-DNA_log2.bw"
subprocess.call(cmd, shell = True)