# ATAC-STARR Analysis 
ATAC-STARR-seq is a massively parallel reporter assay (MPRA) that quantifies regulatory activity within accessible chromatin. In addition to regulatory activty, ATAC-STARR quantifies chromatin accessibility and transcription factor binding. In Hansen, T.J. & Hodges, E. Nucleic Acids Research. 2021, we present the ATAC-STARR method. This repository serves as a complement to that publication. 

In this repository, we share our code for all steps of ATAC-STARR-seq data analysis. A flowchart of the analysis pipeline is below:


![github_figure_flowchart](https://user-images.githubusercontent.com/61889919/122470432-ac1a4400-cf83-11eb-8207-a6e3c1428e56.png)


The scripts used at each step are detailed in the sections below. Every script in this repository is either a slurm script or an R markdown. Slurm scripts are speciallized bash scripts used to submit jobs to the Vanderbilt computer cluster, ACCRE, which is uses the CentOS operating system and the SLURM scheduler. I realize many users may not need this detail, but I left it in because it helps describe the resources required to run each script. To convert the slurm scripts into bash scripts, just remove the ##SBATCH lines. 

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
We estimate library complexity using the preseq package and determine the distribution of insert size using Picard.
```
complexity+InsertSize.slrm
complexity+InsertSize.Rmd
```
## Accessibility Analysis
The inserts of the ATAC-STARR plasmid library are accessible chromatin. In ATAC-seq, these inserts are sent for sequencing rather than massively parallel cloning. Because we sequence the plasmid DNA sample in ATAC-STARR, we wondered if we could use this sequence file as a proxy for ATAC-seq. We turned to the Buenrostro et al. 2013 dataset for benchmarking. Importantly, we downloaded their raw read files and processed using the same read processing methods above. We then called peaks and compared the peaksets using a variety of analyses. In brief, ATAC-STARR can be used to measure chromatin accessibility. 

#### Peak Calling: 
While many peak calling methods exist for ATAC-seq, we prefer Genrich because it includes an ATAC-seq mode and handles biological replicates in a more streamlined fashion. Moreover, Generich is the recommended ATAC-seq peak caller by Harvard FAS Informatics (https://informatics.fas.harvard.edu/atac-seq-guidelines.html). Mapped read files without duplicates are used.  
```
genrich.slrm
```
#### Calculate peak set overlap:
```
jaccard.slrm
euler-plot.Rmd
```
#### Correlation between experiments:
```
correlation.slrm
correlation_plotting.Rmd
```
#### Motif enrichment analysis:
```
motif-enrichment.slrm
motif-enrichment_plotting.Rmd
```
#### Signal track visualization:
```
generate-bigwigs.slrm
visualize_bigwigs_Sushi.Rmd
```
#### Calculate FRiP scores:
```
calculate_FRiP_score.slrm
calculate_FRiP_score_plotting.Rmd
```
## TF Footprinting
TF footprinting is an established method to identify TF footprints from an ATAC-seq or DNase-seq dataset. We performed footprinting on the same ATAC-STARR plasmid DNA sample used in the accessibility analysis above using the TOBIAS software package. 
```
footprinting.slrm
footprinting_heatmap.slrm
```
## Activity Peak Calling
We called active and silent peaks in two ways. We suggest using the "sliding Windows" approach for ATAC-STARR analysis; however, scripts for both methods are listed below. The "sliding window" and "fragment groups" methods both use DESeq2 to call peaks. This works in the following manner: 1) a set of genomic regions is defined, 2) overlapping RNA and DNA reads are assigned to those regions, and 3) DESeq2 identifies regions where the RNA count is statistically different from the DNA count at a false-discovery rate of ≤ 0.1. 
The two methods are different in how they define the input regions at step 1. The "sliding window" approach creates 50bp bins at 10bp step size within each accessiblilty peak. Whereas, the "fragment groups" approach creates regions by combining paired-end sequence fragments that overlap by 75%. The fragment groups approach follows the strategy used in HiDRA (Wang et al. 2018). Scripts for both approaches are provided below. 
#### Sliding Windows
```
make-bins+counts-matrix.slrm
call-peaks_DESeq2.Rmd
merge_active+silent_bins.slrm
```
#### Fragment Groups
```
merge_bams.slrm
create_fragment_groups.slrm
featureCounts_fragmentGroups.slrm
call-peaks_DESeq2.Rmd
remove_redundancy.slrm
```
We compared the two methods by calculating FRiP and plotting region count with the following scripts:
```
compare_peak-callers_FRiP_scores.slrm
compare_peak-callers_plotting.Rmd
```
## Active and Silent Region Characterization
We performed a set of analyses to characterize our active and silent peak sets. 

#### Peak size distribution:
```
size_distribution.Rmd
```
#### Annotation and ChromHMM assignment:
```
annotation.Rmd
ChromHMM_assignment.slrm
ChromHMM_assignment_plotting.Rmd
```
#### Heatmaps of activity and ENCODE signal:
```
generate_ATAC-STARR-activity_bigwig.slrm
generate_signal_heatmap.slrm
```
#### Motif enrichment analysis:
```
motif-enrichment.slrm
motif-enrichment_plotting.Rmd
```
## Multiomic Analysis
We integrated accessibility, TF footprinting, and regulatory activity obtained from ATAC-STARR.

To visualize signal and regions at the GPR132 locus, we used the Sushi package in R. 
```
visualize_multiomics_GPR132.Rmd
```
We clustered active regions based on the presence and absence of TF footprints. We used the following scripts:
```
clustering_part1.slrm
clustering_part2.Rmd
```
