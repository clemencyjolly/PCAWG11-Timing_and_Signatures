# PCAWG11-Timing_and_Signatures
Code for mutational timing of gains and signatures analysis

This repository contains code used for the mutational timing of gains, and the mutational signatures analysis, as reported in _The evolutionary history of 2,658 cancers_, bioRxiv. 
(https://www.biorxiv.org/content/early/2017/07/11/161562)  

## Overview

There are two separate R scripts for the mutational timing of gains, and the timing of signatures. Each script performs the corresponding analysis on one PCAWG sample, and can be run in parallel across the entire dataset.  

The mutational timing of gains proceeds as follows:
* SNVs are assigned as clonal or subclonal
* Clonal SNVs are overlapped with clonal gained segments from the PCAWG-11 consensus copy number (segments level a-d, major allele > 1)
* Timing of segments obtained with "eventTimingOverList" from R package "cancerTiming"
NB - as a prerequisite, this analysis requires the following per sample: clustered SNVs, a consensus copy number profile (with clonal and subclonal segments), estimates for tumour purity and ploidy. 
Per sample, this should take approximately 5-10 minutes.

The timing of signatures proceeds as follows:
* Timed multinomials of SNVs, MNVs, and 4 major classes of indels are extracted from VCFs
* Timed signature weights for SNVs and MNVs are derived using NNLS
* Confidence intervals for signature changes obtained by bootstrapping
NB - as prerequisite, this analysis requires mutational signature compositions, and their activity across samples. For each sample, qualitative time estimates (i.e. early clonal, late clonal, clonal NA and subclonal) for SNVs, MNVs, and indels are required.
Per sample this should take approximately 10-25 minutes.

It should be noted that certain input files are protected data and are not available in this repository. 

## Dependencies

Software packages that are used as part of the analyses. Installing and loading these libraries should take a few minutes.

```
R (version 3.4.1)
```

R libraries for the timing of gains

```
VariantAnnotation
cancerTiming
```

R libraries for the timing of signatures

```
VariantAnnotation
BSgenome.Hsapiens.UCSC.hg19
nnls
reshape2
plyr
```

## Running the R scripts

With all of the required input files available, the two R scripts can be run as follows. 

For the mutational timing of gains:

```
Rscript --no-restore PCAWG_timing.R samplename
```

And similarly, for the timing of signatures:

```
Rscript --no-restore PCAWG_signatures.R samplename
```

## Output

The output of PCAWG_timing.R is a file containing time estimates for all gained segments in a given sample. 
The columns are the default as provided by "cancerTiming" and are as follows:
* "sample" = PCAWG tumour wgs aliquot id
* "pi0" = the time estimate, a value between 0 and 1
* "lCI" = lower 95% confidence interval for time estimate
* "uCI" = upper 95% confidence interval for time estimate
* "N" = no. mutations in segment
* "type" = type of segment
* "segID" = unique segment ID
* "rankInSample" = rank of event based on timing

The output of PCAWG_signature.R are two files, one containing the timed signature weights and their changes (ending in sig_weights.txt), the other containing the bootstrap replicates of the signature changes (ending in sig_change_bootstraps.txt.gz).
The columns in the file containing the changes are as follows:
* "sample" = PCAWG tumour wgs aliquot id
* "signature" = mutational signature
* "clonal [early]", "clonal [late]", "clonal [NA]", "clonal" and "subclonal" = proportional signature weights per epoch
* "early_late" and "clonal_subclonal" = corresponding signature changes
* "el_lCI", "el_uCI", "cs_lCI" and "cs_uCI" = 95% confidence intervals for changes (early-late, clonal-subclonal, respectively)
* "n_early", "n_late", "n_clonalNA", "n_clonal", "n_subclonal" = no. of mutations assigned to each signature per epoch

In the file containing the bootstrap replicates of signature changes, there is one column per signature, and one row per replicate. 
