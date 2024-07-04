# BINF 201 – Practical 2 – Quality Control

In this tutorial we will have a look at quality control of sequencing data. We will have a look at some commonly used QC tools, and what to look out for during QC.

For this practical, you will need to have installed and configured Conda. If this is not the case, see [here](https://docs.conda.io/en/latest/miniconda.html) for installation instructions for your system.

> For Windows users, it is preferable to use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install), and install conda on there using the Linux install instructions.

## Software installation & data retrieval

Open a new terminal (Mac/Linux) or open WSL (Windows)

We will install both Cutadapt (a program to remove adapters from sequencing reads) and FastQC (a QC tool for sequencing data) using conda, in an environment called “QC”

    conda create -n QC -c bioconda -c conda-forge cutadapt fastqc

Once everything is installed, you can activate the environment:

	conda activate QC

Then, download the data for this practical from MittUib (Files->Exercises->Data->Practical2_QC.zip), and unzip the archive. Go into the QC_data directory, and continue with the tutorial.


## FastQC

We will start by running FastQC on some non-processed sequencing reads 

In your terminal run:

    fastqc 1_explore_fastqc_dataset.fastq

This will create 2 files: a `.fastq.html` file containing the fastqc report, and a `.zip` file containing the data used to create the plots in the html file. You can ignore the `.zip` file for this tutorial.

Open the report (`.html`) file and answer following questions:

1.  How many reads are there in the dataset?

2.  What is the read length for this sequencing run?

3.  What is the average %GC (= GC content) of this dataset?

4.  If we would sequence a library containing only “A” bases (e.g. Poly-A tails), what would the %GC be?

Consider the following dummy fastq sequence:

    @ERR899422.1 HS002:1:C600UACXX:2:1101:2590:1987/1
    GCTANTTA
    +
    BBBF#0BF

5.  What do each 4 lines represent?

6. Why must the second and fourth line have the same number of characters?

Go back to the FastQC report

7.  What is the average quality score of the majority of reads?

8.  What does this score mean in terms of sequencing errors? 
(hint: remember what the phred-score stands for)

9.  The “Per Base sequence content” window shows us there seems to be something wrong (red X). What do you think is the problem?

> The observed pattern is very often observed in Illumina data, and nothing to worry about. It is the result of unbiased fragmentation of the DNA (i.e. when fragmenting the DNA, the cuts are not randomly but are more likely to happen after certain bases rather than others)

10.  Look at the “Per Base N content” window. What information can you gather from this plot?

11. What does an “N” signify in a sequence?

12.  In the “Adapter content” window, the adapter content is very slightly going up near the end of the sequence. Is this expected? Why (not)?

## Adapter Trimming

Since we noticed some adapters might still be present in the reads, we need to trim our reads to remove these adapters. Furthermore, as quality drops generally near the end of the sequence, we can also trim low-quality bases while we are at it. Here we will use a tool called Cutadapt that does both at the same time.

Start by running FastQC again, but this time on the `2_adapter_dataset.fastq.gz` file

> FastQC can handle gzipped files, so there is no need to decompress first

1.  Check the adapter content section. Which adapter is present in this dataset?

2.  Is the adapter sequence situated at the 5’ or 3’ end?

As adapter sequences are not informative, and can yield problems in further analyses, we will remove them using Cutadapt. However, to remove the adapters, we first have to know which adapter sequence was used. You can ask your sequencing provider, or most adapter sequences can also be found on Illumina’s website (e.g. [here](https://support-docs.illumina.com/SHARE/AdapterSeq/Content/SHARE/AdapterSeq/AdapterSequencesIntro.htm)). For this practical, an adapter file with some common Nextera adapters has been provided (`adapters.fasta`). 

Now run Cutadapt as follows:

    cutadapt -a file:adapters.fasta -o 2_adapters_trimmed.fastq.gz 2_adapter_dataset.fastq.gz

> Similar to FastQC (and many other bioinformatics software), cutadapt can directly handle gzipped files without you needing to decompress them first.

Once this has completed, run FastQC on the newly created `2_adapters_trimmed.fastq.gz` file

3.  Have the adapters been successfully removed?

4.  FastQC now reports a warning with the Sequence Length Distribution. What is the problem, and why do you think this problem occurred only after running Cutadapt?

Open at the [Cutadapt manual](https://cutadapt.readthedocs.io/en/stable/guide.html) page.

5.  If we had adapter sequences only at the start of our sequence, what would we have to change in the Cutadapt command to make Cutadapt remove them?

## Quality Trimming

In the sequencing lecture we discussed the basics of Illumina sequencing. Briefly, the DNA is fragmented and each fragment is ligated to the flow cell. These fragments are then locally amplified using PCR to produce identical clusters. During sequencing, the same base will be added in all copies of a cluster, yielding a stronger fluorescent signal. This makes the clusters easier to “read” and easier to perform base calling. However, over time some copies get out of sync as they might incorporate two instead of one base, or fail to incorporate a base. As more and more of copies in a cluster get out of sync, the fluorescent signal will become weaker, as a mix of colors will be present instead of just one. This makes it also more difficult to accurately call the bases, leading to lower quality scores in the bases near the end of the reads. To improve the quality of later analyses, it is better to trim bases with low quality scores from the sequencing reads.

We will look at a dataset containing many of these low-quality bases. 
Start by running FastQC on the `3_basequality_dataset.fastq` file.

1.  What is the lowest observed quality score in the dataset?

2.  At which position in the read is this score found?

We will use Cutadapt again, but this time to remove low-quality bases instead of adapters. Run Cutadapt as follows:

    cutadapt -q 30 -o 3_trimmed_output.fastq 3_basequality_dataset.fastq

Now run FastQC on the newly created `3_trimmed_output.fastq` file.

3.  What is now the lowest observed quality score? Is this what you expected based on the Cutadapt command?

4.  Did something change in the read lengths and their distribution? Why would this be the case?

## Nanopore data

In this last part, we will have a look at some RNA transcripts from Illumina and Oxford Nanopore.
Third generation sequencing allows us to generate longer reads than with for example Illumina. Try running FastQC on the `nanopore.bam` file, and open the report.

1.  In the Per Base sequence quality window, there are no boxplots for reads longer than +- 2500 bp. Why would that be the case?

Briefly look at the rest of the report, and notice how it differs from a the Illumina runs you saw before.

> You might have noticed that these files are `.bam` files instead of `.fastq` files. Bam files contain information about read mapping, and the provided bam files (`illumina.bam` and `nanopore.bam`) contain sequenced transcripts (RNA) that were mapped to the human genome. We'll see more about .bam files and mapping reads in the mapping lecture and practical. In any case, since `.bam` files also contain information about the reads and their quality, we can run FastQC on them.

We will visualize the provided `.bam` files using IGV. Launch IGV and load the GRCh38/hg38 version of the human genome. After, load in the `nanopore.bam` and `illumina.bam` files. Once loaded, search for the “TP53” region.

> If you decide to move the `.bam` files before opening them in IGV, make sure to also move the `.bam.bai` file to the same location. If not, you won't be able to load the file in IGV

2. Do you see any mapped reads for Illumina and/or Nanopore to this gene?

3.  You will likely notice that not the whole genomic region is covered. The reads only map to certain regions of the gene. What are these regions, and why do reads only map there?
(hint: what kind of data are we looking at?)

Expand the Refseq Genes section by right-clicking on the row and selecting “expanded”.

4.  Looking at all possible transcripts, you could define 3 different types of transcripts. In what way do these three groups differ from each other?

5.  Are all three of these groups expressed, based on the mapping data?

6.  Which of the two dataset is more useful to assess the previous question, and why?

Zoom in on one of the regions with high coverage (many mapped reads), so you can better see the locations of indels (insertions/deletions) and mismatches.

7.  Which of both datasets seems to be more accurate in sequencing?
