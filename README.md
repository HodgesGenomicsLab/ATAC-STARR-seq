[![DOI](https://zenodo.org/badge/377602391.svg)](https://zenodo.org/badge/latestdoi/377602391)
# ATAC-STARR Analysis 
ATAC-STARR-seq is a massively parallel reporter assay (MPRA) that quantifies regulatory activity within accessible chromatin. In addition to regulatory activty, ATAC-STARR quantifies chromatin accessibility and transcription factor binding. In Hansen, T.J. & Hodges, E. 2021, we present the ATAC-STARR method. This repository serves as a complement to that manuscript. 

In this repository, we share our code for all steps of ATAC-STARR-seq data analysis. A flowchart of the analysis pipeline is below:

![github_figure_flowchart](https://user-images.githubusercontent.com/61889919/124975325-879a1080-dff3-11eb-9c8f-fcd3b397dd45.png)

The code used at each step is detailed as a Jupyter notebook in the sections below. We also supply some python scripts. 

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
We provide a jupyter notebook describing this portion of the analysis. 
```
fastq-processing-and-QC.ipynb
```
In summary, raw reads are trimmed, assessed for quality, mapped to hg38, and filtered to remove regions mapping to ChrM and ENCODE blacklist regions. Reads with MAPQ scores less than 30 are also removed. For the accessibility analysis, we also generated deduplicated mapped read files. 

To automate this step, I developed and utilitezed the following python script. 
```
fastq_processing.py
```
We also estimate library complexity using the preseq package and determine the distribution of insert size using Picard: (Supplementary Figure 1A, Supplementary Figure 2B,C). 

## Accessibility Analysis
The inserts of the ATAC-STARR plasmid library are accessible chromatin. In ATAC-seq, these inserts are sent for sequencing rather than massively parallel cloning. Because we sequence the plasmid DNA sample in ATAC-STARR, we wondered if we could use this sequence file as a proxy for ATAC-seq. We turned to the Buenrostro et al. 2013 dataset for benchmarking. Importantly, we downloaded their raw read files and processed using the same read processing methods above. We then called peaks and compared the peaksets using a variety of analyses. In brief, ATAC-STARR can be used to measure chromatin accessibility. 

#### Peak Calling: 
While many peak calling methods exist for ATAC-seq, we prefer Genrich because it includes an ATAC-seq mode and handles biological replicates in a more streamlined fashion. Moreover, Generich is the recommended ATAC-seq peak caller by Harvard FAS Informatics (https://informatics.fas.harvard.edu/atac-seq-guidelines.html). Mapped read files without duplicates are used.  
```
call-peaks_genrich.slrm
```
#### Calculate peak set overlap: (Figure 2A)
```
euler-plot.Rmd
```
#### Correlation between experiments: (Figure 2B)
```
correlation.slrm
correlation_plotting.Rmd
```
#### Motif enrichment analysis: (Figure 2C-D)
```
motif-enrichment.slrm
motif-enrichment_plotting.Rmd
```
#### Signal track visualization: (Figure 2E)
```
generate-bigwigs.slrm
visualize_bigwigs_Sushi.Rmd
```
#### Calculate FRiP scores: (Figure 2F)
```
calculate_FRiP_score.slrm
calculate_FRiP_score_plotting.Rmd
```
## Active and Silent Region Calling
We called active and silent regions using two approaches: the "sliding windows" method and the "fragment groups" method. We suggest using the "sliding windows" method for ATAC-STARR analysis; however, scripts for both approaches are listed below. Both methods use DESeq2 to call regions but differ in how they define input regions. In the the "sliding window" approach, 50bp sliding bins are created at a 10bp step size within each accessiblilty peak. On the other hand, the "fragment groups" method creates regions by combining paired-end sequence fragments that overlap by 75%. The fragment groups approach follows the strategy used in HiDRA (Wang et al. 2018). Scripts for both approaches are provided below. 

For both methods, active and silent regions are called in the following manner: 1) a set of genomic regions is defined, 2) overlapping RNA and DNA reads are assigned to those regions, and 3) DESeq2 identifies regions where the RNA count is statistically different from the DNA count at a false-discovery rate of ≤ 0.1. Importantly, the mapped read files used in active and silent region calling contain duplicates. This is contrary to what was done in the accessiblity peak calling step above. 

#### Sliding Window: (Figure 3B-D, Supplementary Figure 3C)
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
We compared the two methods by plotting region count and a venn diagram of overlap: (Supplementary Figure 3A-B,D)
```
compare_peak-callers.Rmd
```
## Active and Silent Region Characterization
We performed a set of analyses to characterize our active and silent region sets. 

#### Region size distribution: (Figure 3F)
```
size_distribution.Rmd
```
#### Annotation and ChromHMM assignment: (Figure 4A,B)
```
annotation.Rmd
ChromHMM_assignment.slrm
ChromHMM_assignment_plotting.Rmd
```
#### Heatmaps of activity and ENCODE signal: (Figure 4C,D)
```
generate_ATAC-STARR-activity_bigwig.slrm
generate_signal_heatmap.slrm
```
#### Motif enrichment analysis: (Figure 4E,F)
```
motif-enrichment.slrm
motif-enrichment_plotting.Rmd
```
## TF Footprinting 
TF footprinting is an established method to identify TF footprints from an ATAC-seq or DNase-seq dataset. We performed footprinting on the same ATAC-STARR plasmid DNA sample used in the accessibility analysis above using the TOBIAS software package: (Figure 5A-C)
```
footprinting.slrm
footprinting_heatmap.slrm
```
## Multiomic Analysis
We made a pie chart of the accessible regions based on overlap of an active, silent, or both regions: (Figure 3E)
```
peaks_within_accessible-peaks.txt
peaks_within_accessible-peaks.Rmd
```
To visualize signal and regions at the GPR132 locus, we used the Sushi package in R: (Figure 6)
```
visualize_multiomics_GPR132.Rmd
```
We clustered active and silent regions based on the presence and absence of TF footprints and then performed reactome pathway analysis of each cluster: (Figure 7A,B) 
```
clustering_part1.slrm
clustering_part2.Rmd
```
