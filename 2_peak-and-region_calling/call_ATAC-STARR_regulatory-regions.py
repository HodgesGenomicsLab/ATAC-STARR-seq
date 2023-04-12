"""
name    | Tyler Hansen
created | 2021.09.28
updated | 2023.04.12

Description: Active and silent region calling from ATAC-STARR data. This script consists of four major parts: 
    1) create sliding windows from a region file. 
    2) counting reads for all sliding windows 
    3) determine if sliding window bins are active or silent using RNA-DNA differential analysis with DESeq2. An R script I wrote will be executed.
    4) collapse overlapping bins into a single region and average regulatory score for merged bins. 

Input:
   region file to be converted into windows
   all of the bam files for all samples to have reads mapped to the windows

Output:
   5 BED files: bins analyzed, active bins, silent bins, active regions, silent regions
   The cts matrix analyzed
   R objects for plotting/troubleshooting: deseqdataset(dds) and lfc-shrunk results (res)

Parameters:
-i (--ChrAcc_peaks): chromatin accessibility peaks to analyze (in narrow-peak format)
-d (--DNA_bams): reisolated plasmid DNA bam files. This script will detect and print the number of replicates provided in the analysis. 
-r (--RNA_bams): reporter RNA bam files
-o (--out_dir): output directory for all of the results 
-bs (--bin_size): bin size of windows made
-ss (--step_size): step size of sliding windows
-q (--fdr): False Discovery Rate threshold to apply
-cf (--count_filter): number of raw read counts per bin to filter against
-h (--help): display this help

Required Software downloaded via conda:
    Python v3.6.13
    pybedtools v0.8.2
    subread v2.0.1
    r-essentials v4.1
    bioconductor-deseq2 1.32.0
    bioconductor-apeglm 1.14.0
    bioconductor-biocparallel 1.26.0    
"""

###
# software
import os
import pandas as pd
import sys
import pybedtools
import subprocess
import argparse

#Check that the R script is located in the same directory as this one. If not, error out. 
if os.path.isfile('RNA-to-DNA_differential-analysis.r') == False:
    sys.exit("ERROR: RNA-to-DNA_differential-analysis.r script cannot be located. Make sure it is contained in the same directory as this python script!")
###

###
#   arguments. Make all arguments required. 
#make arg_parser object
arg_parser = argparse.ArgumentParser(description="Call ATAC-STARR Regulatory Regions")

#make a required arguments group
required = arg_parser.add_argument_group('required arguments')

#Add arguments 
required.add_argument("-i", "--ChrAcc_peaks", type=str, default=False,
                        help="accessibility peaks (narrow-peak format)")

required.add_argument("-d", "--DNA_bams", type=str, nargs='+',
                        help="reisolated plasmid DNA bam files")

required.add_argument("-r", "--RNA_bams", type=str, nargs='+',
                        help="reporter RNA bam files")

required.add_argument("-o", "--out_dir", type=str,
                        help='output directory')

arg_parser.add_argument("-bs", "--bin_size", type=int, default=50,
                        help='bin size of windows made (default: 50)')

arg_parser.add_argument("-ss", "--step_size", type=int, default=10,
                        help='step size of sliding windows (default: 10) (turn off tiling by setting -ss = -bs)')

arg_parser.add_argument("-q", "--fdr", type=float, default=0.05,
                        help='false discovery rate for differential analysis (default: 0.05)')

arg_parser.add_argument("-cf", "--count_filter", type=int, default=0,
                        help='number of raw read counts per bin to filter against (default: 0)')

arg_parser.add_argument("-n", "--threads", type=int, default=1,
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
ACC_PEAKS = args.ChrAcc_peaks
DNA_FILES = args.DNA_bams
RNA_FILES = args.RNA_bams
BIN_SIZE = args.bin_size
STEP_SIZE = args.step_size
FDR = args.fdr
CF = args.count_filter
OUT_DIR = args.out_dir
NUM_THREADS = args.threads

# print parameters for troubleshooting purposes
print("PARAMETERS:")
print("ChrAcc Peaks File:", ACC_PEAKS)
print("DNA Files:", *DNA_FILES)
print("RNA Files:", *RNA_FILES)
print("Output Directory:", OUT_DIR)
print("Window Bin Size:", BIN_SIZE)
print("Sliding Window Step Size:", STEP_SIZE)
print("False-Discovery Rate:", FDR)
print("Read Count Filter:", CF)
print("Using", NUM_THREADS, "cores for processing")
###

###
#   determine internal parameters

# calculate the number of replicates and give errors if replicate number is differnt for RNA and DNA or if they are > 5 or < 2.
if len(DNA_FILES) != len(RNA_FILES):
    sys.exit("ERROR: The number of DNA replicates provided does not match the number of RNA replicates provided")

if len(DNA_FILES) == 1:
    sys.exit("ERROR: Only one replicate detected. At least two replicates are required.")

if len(DNA_FILES) > 5:
    sys.exit("ERROR: More than five replicates detected. This script only allows up to five replicates in its current form.")

if len(DNA_FILES) == len(RNA_FILES):
    REPS = len(DNA_FILES)
    print(REPS, "replicates provided")
    #add extra line so the parameters stand out. 
    print("")
###

###
# Disable the following warning:
# A value is trying to be set on a copy of a slice from a DataFrame.
# Try using .loc[row_indexer,col_indexer] = value instead
pd.options.mode.chained_assignment = None
###

###
#   Functions

# # Function: generate sliding windows for downstream analysis 
def make_windows(peaks_file, bins_file, bin_size, step_size):

    ## Prepare table
    #read in narrow peaks file to analyze for activity
    peaks_np = pd.read_table(peaks_file, names=['Chr', 'Start', 'End', 'PeakID', 'Score', 'Strand', 'signalValue', 'pValue', 'qValue', 'peak'])
    #convert to 6 column bedfile. 
    peaks_bed = peaks_np.filter(items=['Chr', 'Start', 'End', 'PeakID', 'Score', 'Strand'])
    #convert to bedtool
    peaks_bedTool = pybedtools.BedTool.from_dataframe(peaks_bed)
    
    ## Create and filter sliding windows
    #make sliding windows with bedtools 
    bins = peaks_bedTool.window_maker(b=peaks_bedTool, w=bin_size, s=step_size)
    #covert to pandas dataframe
    bins_df = pd.read_table(bins.fn, names=['Chr', 'Start', 'End'])
    #Calculate bin size.
    bins_df['size'] = bins_df.End-bins_df.Start
    #keep bins that are the bin size.
    bins_filtered = bins_df[bins_df['size'] == bin_size]
    
    ## Format to SAF
    #adjust start to 1-based coordinate system. FeatureCounts uses 1-bases coordinates. 
    bins_filtered['Start'] =  bins_filtered['Start'] + 1
    #Add strand column. 
    bins_filtered['Strand'] = "."
    #Add the index plus 1 (so it is starting at 1 vs 0) to a bin_ string. 
    bins_filtered['GeneID'] = "bin_" + (bins_filtered.index + 1).astype(str)
    #Select the saf columns
    bins_saf = bins_filtered.filter(items=['GeneID', 'Chr', 'Start', 'End', 'Strand'])
    
    # Write saf to provided bins_file. 
    bins_saf.to_csv(bins_file, sep='\t', header=True, index = False)

# # Function: create a counts table using 2 or 3 ATAC-STARR replicates
def make_cts_table(bins_file, DNA_bams, RNA_bams, cts_file, replicates, threads):
    #define bam files by specifying the list order from the argparse generated list. For reps 3 through 5, assign in the if statement. 
    DNA1 = DNA_bams[0]
    DNA2 = DNA_bams[1]
    RNA1 = RNA_bams[0]
    RNA2 = RNA_bams[1]
    
    #assign feature counts command based on replicate number
    if replicates == 2:
        cmd = f"featureCounts -p -O -B --minOverlap 1 -T {threads} -F SAF -a {bins_file} -o {cts_file} {DNA1} {DNA2} {RNA1} {RNA2}"
    
    if replicates == 3:
        DNA3 = DNA_bams[2]
        RNA3 = RNA_bams[2]
        cmd = f"featureCounts -p -O -B --minOverlap 1 -T {threads} -F SAF -a {bins_file} -o {cts_file} {DNA1} {DNA2} {DNA3} {RNA1} {RNA2} {RNA3}"
    
    if replicates == 4:
        DNA3 = DNA_bams[2]
        RNA3 = RNA_bams[2]
        DNA4 = DNA_bams[3]
        RNA4 = RNA_bams[3]
        cmd = f"featureCounts -p -O -B --minOverlap 1 -T {threads} -F SAF -a {bins_file} -o {cts_file} {DNA1} {DNA2} {DNA3} {DNA4} {RNA1} {RNA2} {RNA3} {RNA4}"

    if replicates == 5:
        DNA3 = DNA_bams[2]
        RNA3 = RNA_bams[2]
        DNA4 = DNA_bams[3]
        RNA4 = RNA_bams[3]
        DNA5 = DNA_bams[4]
        RNA5 = RNA_bams[4]
        cmd = f"featureCounts -p -O -B --minOverlap 1 -T {threads} -F SAF -a {bins_file} -o {cts_file} {DNA1} {DNA2} {DNA3} {DNA4} {DNA5} {RNA1} {RNA2} {RNA3} {RNA4} {RNA5}"

    #run the command
    subprocess.call(cmd, shell = True)

# # Function: call active and silent bins using DESeq2. Here I call an R script that I wrote. 
def call_regulatory_bins(cts_file, out_dir, replicates, FDR, threads, cf):
    
    #assign Rscript as cmd. R script must be in the same directory as this python script
    cmd = f"Rscript RNA-to-DNA_differential-analysis.r --counts {cts_file} --num_reps {replicates} --cores {threads} --count_filter {cf} --FDR {FDR} --out_dir {out_dir}"
    #run the command
    subprocess.call(cmd, shell = True)

# # Function: collapse active and silent bins into regions using bedtools merge.
def merge_bins(out_dir):
    
    #make bedtool objects for each. 
    active = pybedtools.BedTool(os.path.join(out_dir,"active_bins.bed"))
    silent = pybedtools.BedTool(os.path.join(out_dir,"silent_bins.bed"))

    #merge each and isolate the average log2FC score. Sort by position.  
    active_regions = active.merge(c=5, o="mean").sort()
    silent_regions = silent.merge(c=5, o="mean").sort()

    #convert to pandas df. 
    active_regions_df = pd.read_table(active_regions.fn, names=['Chr', 'Start', 'End', 'Score'])
    silent_regions_df = pd.read_table(silent_regions.fn, names=['Chr', 'Start', 'End', 'Score'])

    #add a name column
    #Add the index plus 1 (so it is starting at 1 vs 0) to a region_ string. 
    active_regions_df['Name'] = "region_" + (active_regions_df.index + 1).astype(str)
    silent_regions_df['Name'] = "region_" + (silent_regions_df.index + 1).astype(str)

    #add a strand column 
    active_regions_df['Strand'] = "."
    silent_regions_df['Strand'] = "."

    #rearrange column order
    active_regions_df = active_regions_df[['Chr', 'Start', 'End', 'Name', 'Score', 'Strand']]
    silent_regions_df = silent_regions_df[['Chr', 'Start', 'End', 'Name', 'Score', 'Strand']]

    # Write to out_dir. 
    active_regions_df.to_csv(os.path.join(out_dir, "active_regions.bed"), sep='\t', header=False, index = False)
    silent_regions_df.to_csv(os.path.join(out_dir, "silent_regions.bed"), sep='\t', header=False, index = False)

# # Wrapper Function:  
def call_AS_regions_wrapper(peaks, out_dir, DNA, RNA, reps, num_threads, FDR, bin_size, step_size, cf):
    #define path variables:
    bins = os.path.join(out_dir, "bins.saf")
    cts = os.path.join(out_dir, "cts.tsv")

    #make windows: part 1
    print("***Step 1: Generating sliding windows from the provided ChrAcc peaks file")
    make_windows(peaks_file = peaks, bins_file = bins, bin_size=bin_size, step_size=step_size)
    
    #make cts table: part 2
    print("")
    print("***Step 2: Assigning ATAC-STARR reads to sliding window bins")
    make_cts_table(bins_file = bins, DNA_bams = DNA, RNA_bams = RNA, cts_file = cts, replicates = reps, threads = num_threads)
    
    #call an associated R script to perform differential analysis on counted bins. Part 3. 
    print("")
    print("***Step 3: Performing differential analysis on counted bins")
    call_regulatory_bins(cts_file = cts, out_dir = out_dir, replicates = reps, FDR=FDR, threads = num_threads, cf = cf)
    
    #Merge active and silent bins: part 4
    print("")
    print("***Step 4: Merging active and silent bins into active/silent regulatory regions")
    merge_bins(out_dir = out_dir)
    print("Done!")
###

#Execute wrapper function
call_AS_regions_wrapper(peaks = ACC_PEAKS, out_dir = OUT_DIR, 
                        DNA = DNA_FILES, RNA = RNA_FILES, reps = REPS, 
                        num_threads = NUM_THREADS, FDR=FDR, bin_size = BIN_SIZE, 
                        step_size = STEP_SIZE, cf = CF)