"""
name    | Tyler Hansen
created | 2021.10.15
updated |

This is a python script for processing ATAC-STARR fastq files. 

There are seven functions that each carry out discrete steps of this processing and these are wrapped into one function called process_fastq.

Input:
    fastq read 1 file
    fastq read 2 file

Output:
    trimmed reads with user-defined basename deposited in a new "trimmed_reads" directory
    fastqc analysis results and trimming reports deposited in a user-specified results directory
    bowtie2-mapped reads deposited in a user-defined directory (mapped to the user-specified genome) (.sam)
    filtered bam files (specified file-ending):
        high-quality reads (MAPQ > 30) (.bam)
        remove ChrM reads (.no_chrM.bam)
        remove reads mapping to ENCODE blacklist regions (.filtered.bam)
        remove duplicates (.unique.bam)
    txt files of read counts for each bam file generated above deposited in a user-specified results directory

Parameters:
-a (--input_R1): raw fastq paired-end read 1 file
-b (--input_R2): raw fastq paired-end read 2 file
-n (--name): basename to assign to sample
-e (--blacklist_bed): bedfile of encode blacklist regions (not required) (ex:  ~/hg38_encode_blacklist_ENCFF356LFX.bed)
-t (--output_trimmed): directory to deposit trimmed reads
-o (--output_bam): directory to deposit mapped reads
-r (--output_qc): directory to deposit QC and readcount files 
-g (--genome): bowtie index of genome to map to (ex: /data/hodges_lab/hg38_genome/bt2/hg38)
-h (--help) - display this help

Required Software downloaded via conda:
    Python v3.6.13
    trim-galore v0.6.7
    bowtie2 v2.3.5.1 (also need to downgrade tbb to v2020.2: conda install tbb=2020.2)
    samtools v1.13
    picard v2.26.3 (jar file location:/path/to/venv/share/picard-2.26.3-0/picard.jar) (works for most version numbers using this path: /path/to/venv/share/picard*/picard.jar)
"""

###
# software
import os
import subprocess
import argparse

###
#   Arguments
#make arg_parser object
arg_parser = argparse.ArgumentParser(description="process ATAC-STARR fastq files")

#make a required arguments group
required = arg_parser.add_argument_group('required arguments')

#Add arguments 
required.add_argument("-a", "--input_R1", type=str, default=False,
                        help="raw fastq paired-end read 1 file")

required.add_argument("-b", "--input_R2", type=str, default=False,
                        help="raw fastq paired-end read 2 file")

required.add_argument("-n", "--name", type=str, default=False,
                        help="basename to assign to sample")

arg_parser.add_argument("-e", "--blacklist_bed", type=str, default='none',
                        help="bedfile of ENCODE blacklist regions (default is none and skips this step)")

required.add_argument("-t", "--output_trimmed", type=str, default=False,
                        help='directory to deposit trimmed reads')

required.add_argument("-o", "--output_bam", type=str, default=False,
                        help='directory to deposit mapped reads')

required.add_argument("-r", "--output_qc", type=str, default=False,
                        help='directory to deposit QC and readcount files')

required.add_argument("-g", "--genome", type=str, default=False,
                        help='bowtie index of genome to map to')

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
R1 = args.input_R1
R2 = args.input_R2
NAME = args.name
BLACKLIST_BED = args.blacklist_bed
OUT_TRIM = args.output_trimmed
OUT_BAM = args.output_bam
OUT_QC = args.output_qc
GENOME = args.genome

# print parameters for troubleshooting purposes
print("PARAMETERS:")
print("Fastq R1:", R1)
print("Fastq R2:", R2)
print("Sample name:", NAME)
print("Blacklist Bedfile", BLACKLIST_BED)
print("Trimmed reads directory:", OUT_TRIM)
print("Mapped reads directory:", OUT_BAM)
print("QC results directory:", OUT_QC)
print("Genome:", GENOME)
###

###
#   determine internal parameters

# calculate the number of threads
NUM_THREADS = int(os.getenv('SLURM_CPUS_PER_TASK', 1))
print("Detected and using", NUM_THREADS, "cores for processing")

###
#   Functions

# # Function: trim reads 
def trim(R1, R2, name, output_qc, output_trimmed):

    #Trim reads and assess quality with Trim Galore! The input files can be gzipped.
    #define command
    cmd = f'trim_galore --fastqc --fastqc_args "--outdir {output_qc}" --paired --dont_gzip --basename {name} --output_dir {output_trimmed} {R1} {R2}'
    #execute command
    subprocess.call(cmd, shell = True)

    #Move trimming reports to their directory                                                                                                                              
    cmd = f'mv {output_trimmed}/*trimming_report* {output_qc}'
    subprocess.call(cmd, shell = True)

# # Function: map reads 
def map(trimmed_dir, name, bowtie_index, output_bam, threads):

    #map reads and assess quality with Trim Galore! The input files can be gzipped.
    #define command
    cmd = f'bowtie2 -p {threads} -X 500 --sensitive --no-discordant --no-mixed -x {bowtie_index} -1 {trimmed_dir}/{name}_val_1.fq -2 {trimmed_dir}/{name}_val_2.fq -S {output_bam}/{name}.sam'
    #execute command
    subprocess.call(cmd, shell = True)

# # Function: filter mapped reads
def filter_bam(bam_dir, name, blacklist_bed, threads, output_qc):

    #Convert Sam to Bam and filter for MAPQ > 30, then sort by chromosome. -q filters for mapping quality scores > 30.
    cmd = f'samtools view -@ {threads} -S -b -q 30 {bam_dir}/{name}.sam | samtools sort -@ {threads} - -o {bam_dir}/{name}.bam'
    subprocess.call(cmd, shell = True)

    #Index sorted bam files
    cmd = f'samtools index -@ {threads} -b {bam_dir}/{name}.bam {bam_dir}/{name}.bam.bai'
    subprocess.call(cmd, shell = True)

    #Remove mtDNA reads
    cmd = f"""samtools view -@ {threads}  -h {bam_dir}/{name}.bam | awk '$3 != "chrM"' - | samtools view -S -b - > {bam_dir}/{name}.no_chrM.unsorted.bam"""
    subprocess.call(cmd, shell = True)

    #Sort no chrM bam files
    cmd = f' samtools sort -@ {threads} {bam_dir}/{name}.no_chrM.unsorted.bam -o {bam_dir}/{name}.no_chrM.bam'
    subprocess.call(cmd, shell = True)
    
    #Index no chrM bam files
    cmd = f'samtools index -@ {threads} -b {bam_dir}/{name}.no_chrM.bam {bam_dir}/{name}.no_chrM.bam.bai'
    subprocess.call(cmd, shell = True)

    #filter against blacklisted regions. Use if statement to parse behavior whether blacklist file was included. Default behavior is to assign blacklist_bed as 'none'. 
    if blacklist_bed != 'none':
        cmd = f'samtools view -@ {threads} -b -L {blacklist_bed} -U {bam_dir}/{name}.filtered.unsorted.bam {bam_dir}/{name}.no_chrM.bam > {bam_dir}/{name}.blacklisted.bam'
        subprocess.call(cmd, shell = True)
    if blacklist_bed == 'none':
        cmd = f'samtools view -@ {threads} -b {bam_dir}/{name}.no_chrM.bam > {bam_dir}/{name}.filtered.unsorted.bam'
        subprocess.call(cmd, shell = True)

    #sort filtered bam
    cmd = f'samtools sort -@ {threads} {bam_dir}/{name}.filtered.unsorted.bam > {bam_dir}/{name}.filtered.pos-sorted.bam'
    subprocess.call(cmd, shell = True)

    #index filtered bam
    cmd = f'samtools index -@ {threads} -b {bam_dir}/{name}.filtered.pos-sorted.bam {bam_dir}/{name}.filtered.pos-sorted.bam.bai'
    subprocess.call(cmd, shell = True)

    #name sort filtered bam
    cmd = f'samtools sort -n -@ {threads} {bam_dir}/{name}.filtered.pos-sorted.bam > {bam_dir}/{name}.filtered.n-sorted.bam'
    subprocess.call(cmd, shell = True)

    #remove duplicates
    cmd = f'java -jar /home/hansetj1/.conda/envs/jupyter/share/picard*/picard.jar MarkDuplicates I={bam_dir}/{name}.filtered.pos-sorted.bam \
        O={bam_dir}/{name}.unique.unsorted.bam M={output_qc}/{name}_marked_dup_metrics-all.txt REMOVE_DUPLICATES=TRUE'
    subprocess.call(cmd, shell = True)
    
    #sort unique bam
    cmd = f'samtools sort -@ {threads} {bam_dir}/{name}.unique.unsorted.bam > {bam_dir}/{name}.unique.pos-sorted.bam'
    subprocess.call(cmd, shell = True)

    #index unique bam
    cmd = f'samtools index -@ {threads} -b {bam_dir}/{name}.unique.pos-sorted.bam {bam_dir}/{name}.unique.pos-sorted.bam.bai'
    subprocess.call(cmd, shell = True)

    #name sort unique bam
    cmd = f'samtools sort -n -@ {threads} {bam_dir}/{name}.unique.pos-sorted.bam > {bam_dir}/{name}.unique.n-sorted.bam'
    subprocess.call(cmd, shell = True)


# # Function: count reads with samtools flagstat
def count_reads(bam_dir, name, threads, qc_dir):

    #run flagstat on all bam/sam files to get read counts and save to a txt file in the results directory.  
    #sam
    cmd = f'samtools flagstat -@ {threads} {bam_dir}/{name}.sam > {qc_dir}/{name}.sam_flagstat-results.txt'
    subprocess.call(cmd, shell = True)

    #bam
    cmd = f'samtools flagstat -@ {threads} {bam_dir}/{name}.bam > {qc_dir}/{name}.bam_flagstat-results.txt'
    subprocess.call(cmd, shell = True)

    #no chrM
    cmd = f'samtools flagstat -@ {threads} {bam_dir}/{name}.no_chrM.bam > {qc_dir}/{name}.no_chrM.bam_flagstat-results.txt'
    subprocess.call(cmd, shell = True)

    #filtered
    cmd = f'samtools flagstat -@ {threads} {bam_dir}/{name}.filtered.pos-sorted.bam > {qc_dir}/{name}.filtered.pos-sorted.bam_flagstat-results.txt'
    subprocess.call(cmd, shell = True)

    #unique
    cmd = f'samtools flagstat -@ {threads} {bam_dir}/{name}.unique.pos-sorted.bam > {qc_dir}/{name}.unique.pos-sorted.bam_flagstat-results.txt'
    subprocess.call(cmd, shell = True)


#Master function
def process_fastq(R1, R2, name, output_qc, output_trimmed, output_bam, genome, blacklist_bed, threads) :
    
    #Trim
    print("")
    print("***Trimming raw reads")
    trim(R1=R1, R2=R2, name=name, output_qc=output_qc, output_trimmed=output_trimmed)
    
    #Map
    print("")
    print("***Mapping trimmed reads")
    map(trimmed_dir=output_trimmed, name=name, bowtie_index=genome, output_bam=output_bam, threads=threads)

    #filter
    print("")
    print("***Filtering mapped reads")
    filter_bam(bam_dir=output_bam, name=name, blacklist_bed=blacklist_bed, threads=threads, output_qc=output_qc)
 
    #count reads
    print("")
    print("***Counting mapped reads")
    count_reads(bam_dir=output_bam, name=name, threads=threads, qc_dir=output_qc)
    print("Done!")


### Execute Functions ###
process_fastq(R1=R1, R2=R2, name=NAME, output_qc=OUT_QC, output_trimmed=OUT_TRIM, output_bam=OUT_BAM, genome=GENOME, blacklist_bed=BLACKLIST_BED, threads=NUM_THREADS)
