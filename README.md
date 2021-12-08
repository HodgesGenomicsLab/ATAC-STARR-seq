# ATAC-STARR Analysis 
ATAC-STARR-seq is a massively parallel reporter assay (MPRA) that quantifies regulatory activity within accessible chromatin. In addition to regulatory activty, ATAC-STARR quantifies chromatin accessibility and transcription factor binding. In Hansen, T.J. & Hodges, E. 2021, we present the ATAC-STARR method. This repository serves as a complement to that manuscript. 

In this repository, we share our code for all steps of ATAC-STARR-seq data analysis. The code used at each step is detailed as a Jupyter notebook in the sections below. We also supply some python scripts. 

## Citation
__Please cite our ATAC-STARR-seq method article:__  
TBD - manuscript in revision  
__This repository can also be cited:__    
[![DOI](https://zenodo.org/badge/377602391.svg)](https://zenodo.org/badge/latestdoi/377602391)  
Tyler J Hansen, & Emily Hodges. (2021). HodgesGenomicsLab/ATAC-STARR-seq: Simultaneous profiling of regulatory activity, chromatin accessibility, and transcription factor occupancy with ATAC-STARR-seq (V1.0.0). Zenodo. https://doi.org/10.5281/zenodo.5764666
## Software
Below is a summary of most of the software packages used in ATAC-STARR-seq: 
```
Command-line:       R-based:
  Trim Galore!        R
  bowtie2             tidyverse
  SAMtools            Sushi  
  Picard              DESeq2 
  Preseq              ChIPSeeker
  Genrich             pheatmap
  BEDtools            ggsci
  deepTools           ReactomePA
  HOMER         
  TOBIAS        
  bedops        
  Subread       
```
### Conda 
We use conda to install our packages. We provide a text file of the software packages and their version number used in our conda environment as well as a yml file to clone our environment.
```
pkgs.txt
conda_env.yml
```
## Read Processing and Quality Control (1_read-processing)
We provide a jupyter notebook describing this portion of the analysis. In summary, raw reads are trimmed, assessed for quality, mapped to hg38, and filtered to remove regions mapping to ChrM and ENCODE blacklist regions. Reads with MAPQ scores less than 30 are also removed. For the accessibility analysis, we also generated deduplicated mapped read files. We also estimate library complexity using the preseq package and determine the distribution of insert size using Picard: (Supplementary Figure 1A, Supplementary Figure 2B,C). 
```
fastq-processing-and-QC.ipynb
```
### FASTQ Processing Python Script
To automate this step, I developed and utilitezed the following python script. 
```
fastq_processing.py
```
#### Syntax and Arguments
```
$ python3 fastq_processing.py --help
usage: fastq_processing.py [-h] [-a INPUT_R1] [-b INPUT_R2] [-n NAME]
                           [-e BLACKLIST_BED] [-t OUTPUT_TRIMMED]
                           [-o OUTPUT_BAM] [-r OUTPUT_QC] [-g GENOME]

process ATAC-STARR fastq files

required arguments:
  -a INPUT_R1, --input_R1 INPUT_R1
                        raw fastq paired-end read 1 file
  -b INPUT_R2, --input_R2 INPUT_R2
                        raw fastq paired-end read 2 file
  -n NAME, --name NAME  basename to assign to sample
  -t OUTPUT_TRIMMED, --output_trimmed OUTPUT_TRIMMED
                        directory to deposit trimmed reads
  -o OUTPUT_BAM, --output_bam OUTPUT_BAM
                        directory to deposit mapped reads
  -r OUTPUT_QC, --output_qc OUTPUT_QC
                        directory to deposit QC and readcount files
  -g GENOME, --genome GENOME
                        bowtie index of genome to map to

optional arguments:
  -h, --help            show this help message and exit
  -e BLACKLIST_BED, --blacklist_bed BLACKLIST_BED
                        bedfile of ENCODE blacklist regions (default is none
                        and skips this step)
```
## Accessibility Peak and Regulatory Region Calling (2_peak-and-region_calling)
We provide a jupyter notepook detailing the entire peak and region calling process. The scheme for this is shown in Figure 3A. In brief, de-duplicated DNA bam files analyzed with Genrich to call accessible chromatin (ChrAcc) peaks. These peaks are converted into 50bp sliding window bins (10bp step size) and RNA and DNA read counts for each bin are determined using the duplicated bam files (featureCounts). Differential analysis is performed with DESeq2 to identify active and silent bins. Overlapping bins are then merged into regions. This notebook also generates the plots in Figure 3 and Supp Figure 3. It also performs the de-duplicated and fragment group (Supp Fig 4/5) analyses.  
```
ChrAcc_and_Activity_Region-Calling.ipynb
```
### Region Calling Python Script
As mentioned in the notebook, we developed a python script to call regulatory regions from ATAC-STARR-seq data. 
```
call_ATAC-STARR_regulatory-regions.py
RNA-to-DNA_differential-analysis.r
```
Note: This script calls an R script to do the differential analysis and requires this R script to be in the same directory that the python script is called from. Also, the python script only works for replicate counts 2-5. Anything greater than 5 will require minor revisions to the script.  
#### Syntax and Arguments
```
$ python3 call_ATAC-STARR_regulatory-regions.py --help
usage: call_ATAC-STARR_regulatory-regions.py [-h] [-i CHRACC_PEAKS]
                                             [-d DNA_BAMS [DNA_BAMS ...]]
                                             [-r RNA_BAMS [RNA_BAMS ...]]
                                             [-o OUT_DIR] [-q FDR]
                                             [-n THREADS]

Call ATAC-STARR Regulatory Regions

required arguments:
  -i CHRACC_PEAKS, --ChrAcc_peaks CHRACC_PEAKS
                        accessibility peaks (narrow-peak format)
  -d DNA_BAMS [DNA_BAMS ...], --DNA_bams DNA_BAMS [DNA_BAMS ...]
                        reisolated plasmid DNA bam files
  -r RNA_BAMS [RNA_BAMS ...], --RNA_bams RNA_BAMS [RNA_BAMS ...]
                        reporter RNA bam files
  -o OUT_DIR, --out_dir OUT_DIR
                        output directory

optional arguments:
  -h, --help            show this help message and exit
  -q FDR, --fdr FDR     false discovery rate for differential analysis
                        (default: 0.05)
  -n THREADS, --threads THREADS
                        number of threads (default: 1)
```
## Accessibility Peak and Regulatory Region Characterization (3_peak-and-region_characterization)
We compared our ChrAcc peaks to data from Buenrostro et al. 2013 (Figrue 2), which is detailed in the following notebook:
```
Compare_Buenrostro_and_ATAC-STARR.ipynb
```
We also characterized our active and silent regions in the following notebook (Figure 3):
```
Characterize_Regulatory_Regions.ipynb
```
## TF Footprinting (4_TF_footprinting) 
TF footprinting is an established method to identify TF footprints from an ATAC-seq or DNase-seq dataset. We performed footprinting using the TOBIAS software package (Figure 5). This and the plotting schemes are displayed in the following notebook:
```
TF_footprinting.ipynb
```
## Multiomic Analysis (5_multiomic_analysis)
We created a genome browser view using Sushi  (Figure 6). We also clustered active and silent regions based on the presence and absence of TF footprints and then performed reactome pathway analysis of each cluster (Figure 7). This analysis is displayed in the following notebook:
```
multiomic_analysis.ipynb
```
## Supplementary Analyses (6_supplementary_analysis)
We looked at orientation bias, HiDRA data, and performed a psuedo-replicate analysis. I wrapped these supplementary analyses up into one notebook:
```
Supplementary_Analysis.ipynb
```
