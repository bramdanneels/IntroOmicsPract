---
title: "Practical 7 - Differential Gene Expression"
authors: "Hakon Tjeldnes, Bram Danneels"
output: html_document
---

```{r include = FALSE}
  knitr::opts_chunk$set(eval = FALSE)
```

# Introduction

Our goal for this experiment is to determine which *Arabidopsis thaliana* (a small flower) genes respond to nitrate. The dataset is a simple experiment where RNA is extracted from roots of independent plants and then sequenced. Two plants were treated with the control (KCl) and two samples were treated with Nitrate (KNO3).

# Input data

The sequences were processed to remove all low quality sequences, trim all low quality nucleotides, and finally aligned against the *Arabidopsis* reference genome using TopHat. Our analysis starts from the .bam files created by the alignment/mapping program.

We have been provided the following files:

```         
4 Bam files: Alignment files, one for each sample

Arabidopsis.gtf file: Arabidopsis thaliana reference genome annotation. This contains information about the genes in Arabidopsis and where they are located in the genome.

Experimental design: A comma separated file containing meta data.

Gene description: Description of the function of the genes in Arabidopsis.
```

# Installing required packages

```{r install}
# If we don't have BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}

# Installs all packages you don't have in this list
package.list <- c("Rsamtools", "GenomicFeatures","GenomicAlignments","DESeq2","GOstats","GO.db","Category","org.At.tair.db")
for (package in package.list) { 
  if (!requireNamespace(package, quietly = TRUE))
    BiocManager::install(package, ask = FALSE)
}
  
```

# Generate counts from bam files

First we want to load all the bam files. The normal way to do this is to create a BamFileList object, that contains the location of the bam files, but does not load them into R. This is because bam files can be huge (often more than 5GB) , and can make your computer run out of memory if you load them all at the same time.

Rember to change the directory bellow to the location of the bam files you downloaded together with the assignment.

```{r}
library(Rsamtools, quietly = TRUE)
# Use your correct directory here
data_folder <- "~/Desktop/Bio_data/SCRIPTS/BIN201_Exercise_7_differential_expression/data/"
bamfiles <- BamFileList(dir(path = data_folder, pattern="sorted.bam", full.names = TRUE))
print(bamfiles)
```

1.   Why do we use pattern = "sorted.bam" above, what if you set pattern = "" ?

We now have a list of bam files, containing the information of where the reads map to in the genome.

# Loading the annotation file (.gtf)

The purpose of a .gtf file is to store the location of genome features, such as genes, exons, and CDS. It is a tab delimited flat file which contains information about all annotated features of the genome (location, length, ...)

2.  What are the different columns in the gtf file and what is their purpose? Use this link to read about the columns: gtf format [https://www.ensembl.org/info/website/upload/gff.html](%22%22).

Now we will load the gtf as a database object (called txdb). Again, remember to link to the correct location of the .gtf file on your disk!

```{r}
library(GenomicFeatures, quietly = TRUE)
# Use your correct directory here
data_folder <- "~/Desktop/Bio_data/SCRIPTS/BIN201_Exercise_5_DEseq/data/"
txdb = makeTxDbFromGFF(paste0(data_folder, "Arabidopsis.gtf"), format="gtf")
```

Now, we need to determine which annotated regions we will include for our counting. Since we are dealing with RNA-seq data, we want to count whole transcripts (5' UTR + CDS + 3' UTR), both for protein-coding as non-protein-coding (e.g. ncRNAs) genes. Additionally, we want to use the canonical gene names, and not the transcript names.

```{r}
ebg = exonsBy(txdb, by="gene")
ebg
```

3.  How many genes are present in the .gtf file? (hint, each gene is represented as a GRangesList)
4.  If would only want to count reads mapping to exons, how should we modify the code? (hint: run ?exonsBy to get information about the function)
5.  The GRangesList object "ebg" has many accessory functions that we can use. What do the start(ebg), end(ebg), seqnames(ebg), and strand(ebg) functions do?

# Counting the reads

Now we will use the GenomicAlignments package to determine where the reads were mapped.

The function summarizeOverlaps uses the bam files and the genome annotation to count the number of reads that map to a region that is annotated as a gene. We use the "union" mode, which means that we will consider all exons of a gene as part of the transcript. We also state that the reads that were mapped to the genome were single-end reads, and not paired-end reads.

```{r}
library("GenomicAlignments", quietly = TRUE)
se = summarizeOverlaps(features=ebg, 
    reads=bamfiles, mode="Union", singleEnd=TRUE, 
    ignore.strand=TRUE ) 
```

6.  What is the class of the "se" object?
7.  Look at the dimensions of the "se" object (run dim(se) in the terminal). What do these numbers correspond to?

Now that we know what reads map to which genes, we can create a count table, showing how many reads are mapped to each gene in every condition.

```{r}
counts = assay(se) 
# Remove version numbers on genes (usually a .1)
countgenenames = gsub("[.][1234567890]", "", row.names(counts)) 

rownames(counts)=countgenenames 
```

Some genes have very few reads. This makes it more difficult to do statistical analyses on them, so we generally filter them out.

```{r}
row_means <- rowMeans(counts)
counts_filtered = counts[row_means >= 10,]
dim(counts); dim(counts_filtered)
```

8.  How many genes were filtered out?

# Implementing the experimental design

In order for the differential gene expression to work, we need to tell the program how our experiment is designed. This means we need to label the columns (= experiments), and tell it which experiments to compare. Here, we use a simple comma separated filecreated using a text editor. In that file, the first column corresponds to the bam filename, and the second column to the treatment category.

We will load in the experimental design using read.csv (read.table, read.csv, and read.delim are all related functions but with different default parameters). This will create a data.frame object we can use later on. Additionally it will convert all text (i.e. non-numbers) to factors, which is useful for our downstream analysis. A factor is basically a label that tells us which files/objects/rows belong together (e.g. underwent the same treatment)

```{r}
# Use your correct directory here
data_folder <- "~/Desktop/Bio_data/SCRIPTS/BIN201_Exercise_5_DEseq/data/"
expdesign = read.csv(paste0(data_folder, "expdesign.txt"), row.names=1, sep=",")
expdesign
```

# Calculate deferentially expressed Genes

Using the filtered count table, and the experimental design, we can now perform differential gene expression analysis. We will use DESeq2 for this. DESeq2 performs a pairwise differential expression test using the negative binomial model.

We start by creating an object that serves as input for DESeq2 using the function newCountDataSet. In order to create this dataset, we need the filtered data frame of read counts and the factor that corresponds to the treatment groups (in our case "condition")

```{r}
library(DESeq2, quietly = TRUE)

ds <- DESeqDataSetFromMatrix(countData=counts_filtered,
colData=expdesign,
design= ~ condition)

# if you would like to try to run without the filtering
# simply commend the above lines and uncomment below.

#cds = DESeqDataSetFromMatrix(countData=counts,
# colData=expdesign,
# design= ~ condition)
```

```{r}
analysis <- DESeq(ds)
```

This will create an analysis object, containing all information about the analysis. We can get the actual results using the following command:

```{r}
res <- results(analysis)
head(res)
```

9.   Describe in 2-3 sentences the purpose of the DESeq2 package
10.  What is a p-value? What does it tell us?
11.  What is the log2FoldChange? Why are some values positive, and others negative?
12.  What is the difference between the p-value and the adjusted p-value (padj)?

Now the analysis is done, we can count the number of significantly differentialy expressed genes between the two treatments.

```{r}
sum(res$padj < 0.05, na.rm=T)
```

The sum command is normally used to calculate the sum of numeric values. However when we use it on logic vectors, TRUE=1 and FALSE=0. As such, we can count the number of TRUE values, i.e. the amount of genes with an adjusted p-value below 0.05.

# Plotting differential expression

We can graphically summarize the results in an MA-plot.

```{r}
plotMA(res, ylim=c(-5,5))
```

Now using the data you created during this tutorial, try answering the following questions:

12. Which gene is the most upregulated (highest Log2FC) when subjected to Nitrate?
13. Which gene is the downregulated (lowest Log2FC) when subjected to Nitrate?
14. How many genes are significantly differentially regulated (padj \< 0.05), and have an absolute fold change of at least 2?

# References:

-   [DESeq2 tutorial](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

-   [Inspiration for this tutorial](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/deseq-2/).
