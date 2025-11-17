---
title: Setup
---

This workshop provides a practical introduction to functional enrichment analysis following differential expression in RNA-seq studies. We will compare two major enrichment strategies, **over-representation analysis (ORA)** and **functional class scoring (FCS)**, and discuss when each approach is most appropriate. Participants will learn how to implement these methods in R using packages including `clusterProfiler`, `fgsea`, `Reg-Enrich` and `STRINGdb`, drawing on pathway and gene-set resources such as **Gene Ontology**, **KEGG Pathway Database** and **Molecular Signatures Database**. By the end of this workshop, you will have a clear understanding of how to interpret enriched pathways in RNA-seq data. 

:::: prereq
•	Installed R and RStudio

•	Have basic R knowledge

•	Completed ‘Intro to R for Biologists’ and ‘RNA-seq: From reads to counts to genes’, or equivalent
::::

## R Packages & Datasets

In this workshop, we will learn how to use `clusterProfiler`, `fgsea`, `RegEnrich` and `STRINGdb` tools,
along with related dependencies `org.Mm.eg.db`, `impute` and `preprocessCore`. 
Please install the following packages:

```r

#Run this code to install BiocManager if you don't have it
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

#Use BiocManager to install relevant packages
BiocManager::install(c("clusterProfiler", "fgsea", "RegEnrich","STRINGdb", "org.Mm.eg.db"))

```

We will use RNA-seq data from Fu et al., 2015 (<https://www.ncbi.nlm.nih.gov/pubmed/25730472>). 
Please download the following datasets from Zenodo:

```r

#To download files from Zenodo

dataurl <- "https://zenodo.org/record/2596382/files/"

debasal <- read.csv(paste0(dataurl,"limma-voom_basalpregnant-basallactate"), header = TRUE, sep = "\t")
deluminal <- read.csv(paste0(dataurl,"limma-voom_luminalpregnant-luminallactate"), header = TRUE, sep = "\t")
seqdata <- read.csv(paste0(dataurl,"seqdata"), header = TRUE, sep = "\t")
load(paste0(dataurl,"mouse_hallmark_sets")) #loads as Mm.H
factordata <- read.table(paste0(dataurl,"factordata"), header = TRUE, sep = "\t")
filteredcounts <- read.csv(paste0(dataurl,"limma-voom_filtered_counts"), header = TRUE, sep = "\t")

```

We have also provided these datasets within the GitHub repo. To access them:

```r

#To retrieve files from GitHub

debasal <- read.csv("data/limma-voom_basalpregnant-basallactate", header = TRUE, sep = "\t")
deluminal <- read.csv("data/limma-voom_luminalpregnant-luminallactate", header = TRUE, sep = "\t")
seqdata <- read.csv("data/seqdata", header = TRUE, sep = "\t")
load("data/mouse_hallmark_sets.RData") #loads as Mm.H
factordata <- read.table("data/factordata", header = TRUE, sep = "\t")
filteredcounts <- read.csv("data/limma-voom_filtered_counts", header = TRUE, sep = "\t")

```

## Summary

:::: checklist

Attendees are required to bring their own laptop computers. **Please ensure you have installed:**

- [Chrome](https://www.google.com/chrome/) or [FireFox](https://www.mozilla.org/en-US/)
- [R](https://cran.ms.unimelb.edu.au/) (Download and install the latest version of R using the UniMelb mirror)
- [RStudio](https://posit.co/download/rstudio-desktop/#download)
- R packages required for this workshop (see above)
- Datasets required for this workshop (see above)

::::


