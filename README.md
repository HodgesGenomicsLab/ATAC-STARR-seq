# NAR_ATAC-STARR
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

## Accessibility Analysis

## TF Footprinting

## Activity Peak Calling (Sliding Window)

## Activity Peak Calling (Other)

## Active and Silent Region Characterization

## Multiomic Analysis
