# ATAC-STARR Analysis 
ATAC-STARR-seq is a massively parallel reporter assay (MPRA) that quantifies regulatory activity within accessible chromatin. In addition to regulatory activty, ATAC-STARR quantifies chromatin accessibility and transcription factor binding. In Hansen, T.J. & Hodges, E. Nucleic Acids Research. 2021, we present the ATAC-STARR method. This repository serves as a complement to that publication. 

In this repository, we share our code for all steps of ATAC-STARR-seq data analysis. A flowchart of the analysis pipeline is below:


![github_figure_flowchart](https://user-images.githubusercontent.com/61889919/122470432-ac1a4400-cf83-11eb-8207-a6e3c1428e56.png)


The scripts and their purpose within each step are detailed in the sections below. Every script in this repository is either a slurm script or an R markdown. Slurm scripts are speciallized bash scripts used to submit jobs to the Vanderbilt computer cluster, ACCRE, which is uses the CentOS operating system and the SLURM scheduler. I realize many users may not need this detail, but I left it in because it helps describe the resoruces required to run each script. To convert the slurm scripts into bash scripts, just remove the ##SBATCH lines. 

## Software
ATAC-STARR-seq uses the following publicly available software packages: 

Command line:
```
Trim Galore!  0.6.6 
bowtie2       2.3.4.1
SAMtools      1.9
Picard        2.18.27
Preseq        2.0.0
Genrich       0.5
BEDtools      2.28.0
deepTools     3.3.0
HOMER         4.10
TOBIAS        0.12.7
MACS2         2.1.2
CRADLE        
bedops        2.4.28
Subread       2.0.0
```
RStudio:
```
R             3.6.1
tidyverse     1.3.0
Sushi         1.22.0
DESeq2        1.24.0
apeglm        1.6.0
ChIPSeeker    1.20.0
pheatmap      1.0.12 
```
## Read Processing and Quality Control
For all of the ATAC-STARR analysis, raw reads are trimmed, assessed for quality, mapped to hg38, and filtered to remove regions mapping to ChrM and ENCODE blacklist regions. Reads with MAPQ scores less than 30 are also removed. For the accessibility analysis, we also generated deduplicated mapped read files. 
```
fastq_processing.slrm
```
We estimate library complexity using the preseq package and determine the distribution of insert size using Picard. These two processes are performed in the followiong set of scripts: 
```
complexity+InsertSize.slrm
complexity+InsertSize.Rmd
```
## Accessibility Analysis
The inserts of the ATAC-STARR plasmid library are accessible chromatin. In ATAC-seq, these inserts are sent for sequencing rather than massively parallel cloning, which is how they were used in ATAC-STARR. Because we sequence the plasmid DNA sample in ATAC-STARR, we wondered if we could use this sequence file as a proxy for ATAC-seq. We turned to the buenrostro et al. 2013 dataset for benchmarking. Importantly, we downloaded their raw read files and processed using the same read processing methods above. We then called peaks and compared the peaksets using a variety of analyses. In brief, ATAC-STARR can be used to measure chromatin accessibility. 

Peak Calling: While many peak calling methods exist for ATAC-seq, we prefer Genrich because it includes an ATAC-seq mode and handles biological replicates in a more streamlined fashion. Moreover, Generich is the recommended ATAC-seq peak caller by Harvard FAS Informatics (https://informatics.fas.harvard.edu/atac-seq-guidelines.html). Mapped read files without duplicates are used.  
```
genrich.slrm
```
Calculate peak set overlap:
```
jaccard.slrm
euler-plot.Rmd
```
Correlation between experiments:
```
correlation.slrm
correlation_plotting.Rmd
```
Motif enrichment analysis:
```
motif-enrichment.slrm
motif-enrichment_plotting.Rmd
```
Signal track visualization:
```
generate-bigwigs.slrm
visualize_bigwigs_Sushi.Rmd
```
Calculate FRiP scores:
```
calculate_FRiP_score.slrm
calculate_FRiP_score_plotting.Rmd
```
## TF Footprinting
TF footprinting is an established computational method to identify TF footprints from ATAC-seq or DNase-seq data. We performed footprinting on the same ATAC-STARR plasmid DNA sample used in the accessibility analysis above using the TOBIAS software package. 
```
footprinting.slrm
footprinting_heatmap.slrm
```
## Activity Peak Calling (Sliding Window)
We called active and silent peaks in two ways. The results of these methods are highlighted in the paper. We suggest using the "sliding Windows" approach for ATAC-STARR analysis; however, scripts for both methods are listed below. The two methods are similar because they use DESeq2 to call peaks. Briefly, they work in the following manner: 1) a set of genomic regions is defined, 2) overlapping RNA and DNA reads are assigned to those regions and counted, and finally 3) DESeq2 is used to identify regions where the RNA count is statistically different from the DNA count at a false-discovery rate of ≤ 0.1. 

<img src="https://user-images.githubusercontent.com/61889919/123313257-f6a14080-d4ee-11eb-9815-e6db80ff06d9.png" height="350" width="200px">

The two methods deviate at step 1, defining regions. The "sliding window" approach creates 50bp bins at 10bp step size within each accessiblilty peak. See the figure above. 

The "fragment groups" approach follows the strategy from Wang et al. 2018. Simply, fragments with 75% overlap are combined to create regions representative of overlapping fragments. The scripts for this method are listed below:
```


```
## Active and Silent Region Characterization

## Multiomic Analysis
