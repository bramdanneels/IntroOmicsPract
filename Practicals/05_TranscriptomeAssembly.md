# BINF201 – Practical 4 – Transcriptome Assembly

In this practical we will perform a transcriptome assembly on a toy dataset containing RNA-seq reads derived from a small region (~100 genes) of the human X chromosome. 

## Software installation and data retrieval

In this tutorial, we will have a look at the following software:

- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) - A _de novo_ RNA-seq transcript assembler
- [Stringtie2](https://github.com/skovaka/stringtie2) - A reference-based RNA-seq transcript assembler
- [Hisat2](https://daehwankimlab.github.io/hisat2/) - A fast and sensitive alignment
- [Megahit](https://github.com/voutcn/megahit) - An ultrafast, memory-efficient De Bruijn Graph assembler for short reads
- [Canu](https://github.com/marbl/canu): A long-read assembler for both Nanopore and early (noisy) PacBio reads
- [Flye](https://github.com/mikolmogorov/Flye): Another long-read assembler for all kinds of long reads (high or low error Nanopore & PacBio)
- [HifiAsm](https://github.com/chhylp123/hifiasm): An assembler specific for PacBio HiFi reads
- [Quast](https://github.com/ablab/quast) - An assembly QC tool that generates assembly statistics
- [BUSCO](https://busco.ezlab.org/) - A tool to assess assembly completeness

### For students using NREC

The data and software has been setup on the NREC server. Before starting the practical, make sure to activate the correct environment before each part of the tutorial 
(`QC` for the QC part, `Assembly` for the assembly part)

You can make a copy of the data you will be working on by running this command from your home directory:
> If your not sure if you are in the home folder, run `cd` or `cd ~` to go to your home directory.

```
mkdir -p Practical3
ln -s /home/brdan9637/Pract3_Assembly/Data/* Practical3/
```
> `mkdir -p` creates a folder called Practical2. The "-p" options tells mkdir to create subdirectories if necessary, and to not give an error if the folder(s) already exist
> `ln -s` creates what we call a "symbolic link". This creates a small file that just says "Instead of this file, use the file that I'm liking to". This allows you to "copy" files without actually having to make a physical copy.

Now, go to the newly created directory (by running `cd Practical3`), and you are ready to start!

### For students running on their own pc

You will first have to setup the correct environment with the necessary tools. See [the intro practical](00_IntroSetup.md) on how to install mamba, and how to create an enviroment for downloading the necessary data.

To create a new environment including the necessary tools for this practical, run the following command:

```
mamba create -n Assembly spades abyss megahit quast busco canu flye hifiasm wget gzip
```
> This will create a new environment called "Assembly" with the necessary tools.
> Apart from the assembly tools (see below), we'll also installing `wget` to retrieve data, and `gzip` (de)compress data.

Create a new folder (e.g. `Practical3`), go into it (`cd Practical3`), and then download the necessary data, by either:

- Download directly from the [Zenodo repository]():

```
wget link_to_zenodo_file(s)
gunzip *gz
```

- Download manually using the commands below:

> Remember to activate the download environment (see [here](00_IntroSetup.md))
> Note: Total download size after decompression is +- 1,25 Gb

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_1.fastq.gz
mv ERR486840_1.fastq.gz MycGen_1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_2.fastq.gz
mv ERR486840_2.fastq.gz MycGen_2.fastq.gz
gunzip *gz

prefetch SRR28208385
fasterq-dump -p --outdir ./ --split-files SRR28208385/SRR28208385.sra
mv SRR28208385.fastq MycOvi_Nano.fastq
rm -r SRR28208385/

prefetch SRR24462972
fasterq-dump -p --outdir ./ --split-files SRR24462972/SRR24462972.sra
seqtk sample -s666 SRR24462972.fastq 50000 > MycOvi_HiFi.fastq
rm -r SRR24462972*
```

## The Data

The data for short-read assembly composed of paired-end Illumina sequencing reads from _Mycoplasmoides genitalium_, a pathogenic bacteria that causes infection of the urinary and genital tracts in humans. 
The data for long-read assembly is composed of Nanopore and PacBio HiFi reads from _Mycoplasma ovipneumoniae_, a pathogenic bacterium causing pneumonia in sheep and goats.
For an overview of where the data comes from, see the [information sheet on Zenodo]().

## Tool installation & data retrieval

We start by creating a new conda environment for the transcriptome assembly. You can install the necessary tools using the following command:

	conda create -n transcriptome -c conda-forge -c bioconda trinity samtools minimap2 gffcompare hisat2 stringtie boost=1.60.0
> In case this installation seems to get stuck on solving the environment, try using Mamba instead of Conda (see the end of this document)

Then we will download the necessary data. Create a new folder for this practical (e.g. Transcriptome), and go into it.
Download the data from MittUiB: Excercises->Data->Practical4_Transcriptome.zip, put it in the newly created folder, and unzip it.

Open the `minigenome.fa` file in a text editor of your choice.

1.	 Some of the nucleotides in this file are written in small, others in capital. Why?

2.	There are also some N’s in the sequence. What do they represent?

Then, open the `minigenome.gtf` file as well in a text-editor

3.	What is GTF file, and what do the different columns represent?

## Quality Control

As usual, we will start by looking at the sequencing reads to see if they look OK. Activate the QC environment (See practical 2), and run FastQC on the read files in the data folder:

    fastqc reads_1.fq reads_2.fq

Then open the report .html files and inspect them

1.	How many reads do we have? And what is the read length?

2.	Based on the report, do we need to do some trimming?

## Reference-assembly with HISAT2 and StringTie

First, we will do a genome-guided assembly of the reads. This is of course only possible if there is already a genome available for the species you are investigating. 
First we will use [Hisat2](http://daehwankimlab.github.io/hisat2/), an RNA-aligner. 
>We will learn more about RNA-aligners in the following lectures and practicals. What you need to remember for now is that an RNA-aligner allows reads to be “split-up” so they can align to multiple exons. 

First, we will use the genome annotation (`.gtf` file) to identify the splice sites, so we know where the exons and where the introns are (remember to activate the transcriptome environment first):

	hisat2_extract_splice_sites.py minigenome.gtf > minigenome.gtf.ss

Open the first few lines of this files as follows:
	
	head minigenome.gtf.ss

1.	What do you think each of the columns in the output represent? 

We can also extract the location of the exons:
	
	hisat2_extract_exons.py minigenome.gtf > minigenome.gtf.exons

Open this file as well using the head command

2.	What do you think each of the columns in the output represent? 

3.	Are the differences between column 2 and 3 larger or smaller than in the splice-site file? Is this expected?

Using the splice site and exon information, we can build a genome index, which will allow easier mapping of reads:

	hisat2-build --exon minigenome.gtf.exons --ss minigenome.gtf.ss minigenome.fa minigenome_idx

Now we have everything ready to start mapping the reads we will:
-	Align the paired-end reads to the genome
-	Pipe the resulting mapping information (in SAM format) to samtools view, which will convert it to a more compact, binary file (BAM format)
-	Pipe the resulting BAM file to samtools sort to sort the mapped locations, making it easier for other tools to use the file

We can do this all in one go with this command:

	hisat2 --dta -x minigenome_idx --max-intronlen 5000 \
	-1 reads_1.fq -2 reads_2.fq \
	| samtools view -Sb \
	| samtools sort -o alignments.hisat2.bam

> In case you are not familiar with linux commands, “|” is what we call a pipe, which tells linux to directly transfer the output from one program into another one, without the need to write it to the disk. The “\” is use to split up the command over multiple lines. 

> In the Hisat command, we specify the `--dta` option to tell the software to produce an output for tailored for transcriptome assembly. If you want to do just mapping, you won't use this option.

> The `-Sb` option of `samtools view` tells samtools that the input is in SAM format (`S`), and that we wan tthe output in BAM format (`b`). However, in the more recent version of samtools, it can automatically detect the input file, so the `S` option is technically obsolete.

4.	Look at what was printed to the screen. Did the mapping go well?

Now we will also create an index for this newly created `.bam` file, for easier access to the file by other software:

	samtools index alignments.hisat2.bam

>The SAM and BAM format are used a lot in genomics and transcriptomics. If you want to know more about these formats, you can read up on them here: https://samtools.github.io/hts-specs/SAMv1.pdf

So now we mapped the reads to the genome. But that doesn’t give us the transcriptome. However, we can use the mapping information to piece the reads together and infer the transcript structure. These assembled transcripts can then tell us what genes/transcripts were expressed in the tissue/condition from which they were derived. To piece the mapped reads together, we will use [StringTie](http://ccb.jhu.edu/software/stringtie/):

	stringtie -o stringtie.gtf alignments.hisat2.bam

This outputs a new `.gtf` file, similar to the one of the reference. This `.gtf` contains the predicted transcripts (including exons and splice-sites) based on the mapped reads.
Now we can compare our transcripts with the reference, and with the mapped reads. 

Open IGV and load the genome (`minigenome.fa`), the reference annotation (`minigenome.gtf`), the new annotation (`stringtie.gtf`), and the mapping information (`alignments.hisat2.bam`).
>`.bam` files need a corresponding index file (`.bam.bai`) to be visualized in a genome browser. So if you copy the .bam file somewhere else to load it in IGV, make sure to also copy the `.bam.bai` file!

Once everything is loaded, try zooming in until you can see the mapping information.

5.	Visually compare the reference annotation with our stringtie annotation. Would you consider the stringtie annotation as good or not?

6.	Scroll a bit to the left or right until you find a gene where the reference transcript has more exons than our stringtie annotation. What might be the cause of these missing exons in our annotation? 

You can leave IGV open for now, as we’ll add some data to it later. 

Another way of comparing annotations is using [gffcompare](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml). As the name suggests, this tool compares two gff files, and calculate some statistics that might help you assess the quality of the annotation. 
> `.gtf` files are a speciel type of gff files, so they can be used by most software accepting `.gff` files

We can run gffcompare as follows:

	gffcompare -r minigenome.gtf -o stringtie_gffcomp stringtie.gtf

This will create some different files, but the one we are interested in is the `stringtie_gffcompare.stats` file. Open it in a text editor (or on the command line using cat, less, or any CL text editor of your choice). 

7.	The document gives sensitivity and precision statistics about our annotation. What do these two terms mean, for example when looking at exons?

8.	Do these statistics make sense given the number of missed & novel exons?

## De novo transcriptome assembly with Trinity

In many cases no reference transcriptome is available, or the reference is bad quality. In that case, we can try to build the transcriptome from scratch (*de novo*) using the reads. For this we will use a software called [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki). Normally, Trinity requires a lot of data, a lot of computing power, a lot of disk space, and a lot of running time. Sadly, in this practical we have none of the above. We will do our best with what we have, but keep in mind that in a real-life setting, you’ll very likely be working with a lot more data, and require a computing cluster to get decent results.

In any case, we’ll try and assemble our mini-dataset using Trinity. This assembly can take some time, so if you have multiple cores available on your PC, and know roughly how much memory you have, you can adapt the `--CPU` and `--max_memory` parameters accordingly. 
We will run Trinity as follows (assuming 2 cores and 10 GB of free RAM):

	Trinity --left reads_1.fq --right reads_2.fq --seqType fq --CPU 2 --max_memory 10G --output trinity

Once this is done, we have our assembled transcripts in the trinity folder (`Trinity.fasta`). However, they are just transcripts: assembled reads into larger contigs. We need to map the transcripts to the genome to know where they come from. That will allow us as well to visualize our transcripts in a genome browser. We’ll map the reads using [minimap2](https://github.com/lh3/minimap2):

	minimap2 -ax asm5 --cs minigenome.fa trinity/Trinity.fasta \
	| samtools view -Sb \
	| samtools sort -o trinity/trinity_denovo.bam
	samtools index trinity/trinity_denovo.bam
	
Now we can load the .bam file in IGV and compare it with our other data. 

>`.bam` files need a corresponding index file (`.bam.bai`) to be visualized in a genome browser. So if you copy the .bam file somewhere else to load it in IGV, make sure to also copy the `.bam.bai` file!

1.	Compare the transcripts from Trinity to the reference genes. Did Trinity do a good job according to you?

Remember, we’re only working on a small dataset, with only a couple of thousand reads. Normally these types of analyses are run with millions of reads.
Lastly, we’ll have a look at the structure of our transcriptome assembly. For this we’ll use a software called [Bandage](https://rrwick.github.io/Bandage/). Start by downloading the program from their website.
Open the software, and load the `Trinity.fasta` file we created earlier.

Most transcripts will just be linear, but some will have special structures.

2. How long is the longest transcript?

3.	Can you find examples of alternative splicing in a transcript?

## Using Mamba instead of Conda

Sometimes conda gets stuck wen installing certain packages, especially when a lot of dependencies need to be resolved. That’s why [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) has been developed. Mamba functions the same as Conda, but handles package installation a lot better. To use mamba instead of conda, just install mamba using conda first:

	conda install mamba
Once that is completed, you can now use mamba install or mamba create to install packages and environments instead of conda. However, activating and deactivating environments still uses conda! 

For example for installing the busco environment:
	
	mamba create -n busco -c bioconda -c conda-forge busco=5.4.7
	conda activate busco
	busco -h

