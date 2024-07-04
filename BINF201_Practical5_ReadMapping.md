
# BINF201 – Practical 5 – Read Mapping

For this practical, you will need to have installed and configured Conda. If this is not the case, see [here](https://docs.conda.io/en/latest/miniconda.html) for installation instructions for your system.

> For Windows users, it is preferable to use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install), and install conda on there using the Linux install instructions.

In this practical we will map different kind of sequencing reads from different types of experiments to the *Sacchoromyces cerevisiae* (budding yeast) genome using different mappers.

## Tool installation & data retrieval

As usual, we will create a new conda environment containing the necessary tools:

	conda create -n mapping -c bioconda -c conda-forge hisat2 sra-tools fastp samtools smalt bwa-mem2 star minimap2 pigz

We will be analyzing multiple datasets in this practical.:

|Data type|	Accession|
|-|-|
|WGS - Illumina SE|	ERR560561|
|WGS - Illumina PE|	SRR25410876|
|WGS - IonTorrent|	SRR5830892|
|WGS - Nanopore|	SRR22575474|
|WGS - PacBio|	ERR10357081|
|RNA - Illumina PE|	SRR12006022|
|RNA - Nanopore|	SRR21958970|
|ChipSeq Illumina PE|	SRR1951324|
|*Candida albicans* PE|	SRR22538153|

If you want to know more about these datasets, you can look them up on [SRA](https://www.ncbi.nlm.nih.gov/sra). As most of these datasets are pretty big, you will not download them yourselves. Instead, a subset of reads (~ 2x1M reads for the largest dataset) is provided on MittUiB. 
You can download these files in Exercises -> Data -> Practical5_Mapping.zip. 
This dataset also includes the *S. cerevisiae* reference sequence, and its annotation.

>If you would want to download an SRA dataset yourself using the command line, you can have a look at SRA-tools (or https://github.com/ncbi/sra-tools/wiki/). If you would want to subset a set of reads yourself, you can take a look at [SeqTK](https://github.com/lh3/seqtk).

Download the data from MittUib, and place the `Practical5_Mapping.zip` file in a new folder. 
>While the data is only a subset of the reads, some files still take up a bit of space. Uncompressed, all data files take up almost 8Gb. If you do not have enough disk space available, consider only uncompressing the necessary files for every tutorial, and removing or compressing them again after.

Uncompress all files using gzip or pigz before continuing:

	pigz -d *gz

or

	gunzip *gz

##	Quality Control

We have many sets of reads, so doing QC on all of them would take some time. To make things simpler, semi-automatic QC toosl exist for sequencing reads, for example [Trim-Galore](https://github.com/FelixKrueger/TrimGalore), which combines Cutadapt with FastQC. 
However, this time we will be using [fastp](https://github.com/OpenGene/fastp), a fast and versatile all-in-one FastQ processing tool. For datasets with only one set of reads, run fastp as follows:

	fastp  -i reads.fastq -o reads.trimmed.fastq -h reads.report.html

For paired-end reads run:

	fastp -i reads_1.fastq -I reads_2.fastq -o reads_1.trimmed.fastq -O reads_2.trimmed.fastq -h reads.report.html

>Make sure to use a unique name for every report, to avoid overwriting old reports

Once everything is done, you can take a quick look at the reports to check if everything looks fine.

1.	How does the IonTorrent quality compare to Illumina (PE/SE) reads?

2.	NanoPore can create very long reads, but do they also tend to have high quality?

3.	When looking at the PacBio reads, what stands out when looking at the quality scores?

## Short-read genomic mappers

In this part we will compare different genomic mappers, both for short and long reads. We will start by assessing the performance of two short-read mappers: [BWA-mem2](https://github.com/bwa-mem2/bwa-mem2) and [Smalt](https://www.sanger.ac.uk/tool/smalt/). We will test the mapping success, and their running time on the whole-genome DNA paired-end reads.

Before mapping, we need to index our reference:

	bwa-mem2 index -p SC_BWA_index sacCer3.fa
First, we’ll map the reads using BWA-mem2. We will also time how long the mapping takes, by placing the `time` command before the bwa-mem command. We will immediately convert the SAM output to BAM, and sort them:

	time bwa-mem2 mem -t 2 SC_BWA_index \
	WGS_PE_1.trimmed.fastq WGS_PE_2.trimmed.fastq \
	| samtools view -bS | samtools sort  > WGS_PE_BWA.sorted.bam
	samtools index WGS_PE_BWA.sorted.bam

Then, we’ll do the same for Smalt. Since Smalt is a different mapper than bwa-mem, we need to create another index of the reference first.

	smalt index -k 15 -s 13 SC_smalt_index sacCer3.fa
	time smalt map -n 2 SC_smalt_index \
	WGS_PE_1.trimmed.fastq WGS_PE_2.trimmed.fastq \
	| samtools view -bS | samtools sort > WGS_PE_smalt.sorted.bam
	samtools index WGS_PE_smalt.sorted.bam
	
> The `-k` option of `smalt index` is the length of the hashed words used in the index, and the `-s` option defines the sampling step size, the distance between successive words that are hashed along the genomic  reference. Both these parameters regulate the trade-off between sensitivity and accuracy on one side, and speed and memory efficiency on the other side. Higher *k* values lead to increased mapping accuracy, but can only be used with longer reads (>100bp).  A low stepsize value increases accuracy, as more parts of the genome are indexed, but leads to an increased memory and runtime. For more detailed information and suggestions on what to put for these parameters: https://software.cqls.oregonstate.edu/updates/docs/smalt/smalt_manual.pdf


1.	Which of both mappers took the least time?

Now we will gather some mapping statistics using samtools and awk. First, we will estimate the average coverage of our reference:

	samtools depth -a WGS_PE_BWA.sorted.bam | awk '{c++;s+=$3}END{print s/c}'

> Samtools depth will show us for every base how much reads cover that base. With the `awk` command, we will go over every line and at every line increase c with 1 (`c++`), and increase s with the value of the third column (`s+=$3`), with the third column containing the coverage for that base. At the end, we will calculate the average coverage as the total coverage ( s ), divided by the total number of bases ( c ): `print s/c`

Then we will calculate which proportion of our reference has at least one read covering it:

	samtools depth -a WGS_PE_BWA.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	
> Similar to previous command, we will again count all the bases with `c++`, and each time we have a base which has a coverage (stored in the 3rd column) higher than zero, we will add 1 to the total (`if ($3>0) total+=1`). In the end we will calculate how much of the genome is covered by dividing the number of bases with coverage > 0 (total) by the total number of bases ( c ): `print (total/c)*100`. We multiply by 100 to get the result in percentage instead of proportion.

Lastly, we can calculate the proportion of the reads in our dataset were mapped:
	
	samtools flagstat WGS_PE_BWA.sorted.bam | grep "mapped ("

> If you want to know more about what these commands do exactly, you can have a look [here](https://sarahpenir.github.io/bioinformatics/awk/calculating-mapping-stats-from-a-bam-file-using-samtools-and-awk/).

2. 	Are there large differences between both mappers? Which one do you prefer?

Choose one of both mappers that you like best, and use that mapper to map the WGS_SE and IonTorrent reads to the reference (you don’t need to re-index the genome). Since we only have one read file, you only need to pass one argument to the mapping command instead of two. 
You don’t need to calculate the mapping statistics for these datasets. 

Once you have the `.bam` files and their indexes, open IGV and load the build in *S. cerevisiae* reference genome (sacCer 3). Then, load the mapping files (`.bam` and `.bam.bai`) of the WGS_PE (only load one of the two you created), WGS_SE, and IonTorrent.
 
Select a chromosome in the dropdown menu, and zoom in until you see the mapping information for all three mapping sets, and look around by scrolling left and right.

3.	Compare the three read datsets. What type of error is more prevalent in IonTorrent, compared to Illumina reads?

4.	Zoom in enough so you can see the bases in the reference sequence. Do you observe a pattern in where the IonTorrent mistakes happen mostly?

## Long Read Mapping

In this part, we will map both NanoPore and PacBio reads to the genome. Long read mappers use different approaches for mapping, as in general we are dealing with fewer, but longer reads. Many long-read mappers use approaches similar to those used for genome-to-genome alignment. One such mapper is [minimap2](https://github.com/lh3/minimap2), which can be used for long-read mapping, but can also be used to align transcriptomes or whole genomes to a reference. 

In contrast to the short-read mappers we used, minimap2 can create the index as part of the mapping pipeline. We can thus go directly to mapping the reads:

	minimap2 -ax map-hifi sacCer3.fa WGS_PacBio.trimmed.fastq | samtools view -bS | samtools sort > WGS_PacBio.sorted.bam
	samtools index WGS_PacBio.sorted.bam
	minimap2 -ax map-ont sacCer3.fa WGS_NanoPore.trimmed.fastq | samtools view -bS | samtools sort > WGS_NanoPore.sorted.bam
	samtools index WGS_NanoPore.sorted.bam

> We are telling minimap2 to do both indexing (`-x`) and aligning (`-a`), and using settings adapted to PacBio HiFi (`map-hifi`) or Oxford Nanopore (`map-ont`) data

We will also calculate the mapping statistics again, but only for PacBio:
	
	samtools depth -a WGS_PacBio.sorted.bam | awk '{c++;s+=$3}END{print s/c}'
	samtools depth -a WGS_PacBio.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	samtools flagstat WGS_PacBio.sorted.bam | grep "mapped ("

1.	Compare the mapping statistics to these of the short read mapping. Are there big differences?

2.	Compare the number reads between the WGS_PacBio and WGS_PE sets. Is there a correlation between the number of reads, and the average base coverage of the genome? Can you explain this correlation, or lack thereof?

Now load both mappings in IGV and compare the NanoPore and PacBio reads.

3.	Based on the error profile of both reads sets, which dataset would you prefer?

4.	Go to Chromosome 1 (`chrI`), and look in the region between 25000 bp and 27500 bp. What is happening here?

5.	This weird region is part of the FLO9 gene. Look up the gene in [UniProt](https://www.uniprot.org/). Look at the “Family & Domains” section of the FLO9 entry. Can you identify a possible reason for the observed pattern in the data? 

## RNA-seq Mapping

Now that we have mapped both short an long DNA reads, we’ll have a look at mapping RNA-derived reads. In the trascriptome practical, we already had a look at HISAT2, and here we’ll a look at another popular RNA-seq aligner: [STAR](https://github.com/alexdobin/STAR). 
As you probably know by now, these aligners are “splice-aware”, and can split up the reads to allow mapping to different regions of the genome (the exons). STAR and HISAT2 are both short-read aligners, so won’t work very well with PacBio or NanoPore reads. 

Before we can map using STAR, we need to put the genome in a separate folder. Just create a folder called “genome”, and make a copy of the genome in that folder:

	mkdir genome && cp sacCer3.fa genome/

Then, we will create a genome index:

	STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ./genome \
	--genomeFastaFiles ./genome/sacCer3.fa --sjdbGTFfile sacCer3_refseq.gtf \
	--sjdbOverhang 100 --genomeSAindexNbases 10

> The sjdb related options refer to splice junction database information

Then, we can map the reads using the following command:

	STAR --runThreadN 2 --genomeDir ./genome --sjdbGTFfile  sacCer3_refseq.gtf --sjdbOverhang 100 --readFilesIn RNA_PE_1.trimmed.fastq RNA_PE_2.trimmed.fastq 
This will create a SAM file named `Aligned.out.sam`, which we will convert to BAM, sort, and index as usual:
	
	samtools view -bS Aligned.out.sam | samtools sort > RNA_PE.sorted.bam
	samtools index RNA_PE.sorted.bam

Then we will delete some intermediate files that STAR created:
	
	rm -r Aligned.out.sam Log.* SJ.out.tab _STARgenome

Then let’s assess the mapping statistics:
	
	samtools depth -a RNA_PE.sorted.bam | awk '{c++;s+=$3}END{print s/c}'
	samtools depth -a RNA_PE.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	samtools flagstat RNA_PE.sorted.bam | grep "mapped ("

1.	Compare the mapping statistics to those of the genomic, paired-end reads. What are the main differences, and what is the reason behind these differences?

Now we will also map the RNA nanopore reads. We will use mimimap2 again for this. We will first identify the splice junctions, which will lead to improved mapping results:

	paftools.js gff2bed sacCer3_refseq.gtf > sacCer3_refseq.bed
	minimap2 -ax splice --junc-bed sacCer3_refseq.bed -uf -k14 sacCer3.fa RNA_NanoPore.trimmed.fastq | samtools -view -bS | samtools sort > RNA_NanoPore.sorted.bam
	samtools index RNA_NanoPore.sorted.bam

> The `-uf` option forces minimap2 to only use the forward strand, as we're dealing with RNA transcripts. The `-k` parameters sets the word length for the genome index. The `splice` option, as you might have guessed, tells minimap2 to do spliced alignment.

Now, open IGV again, and load both RNA-datasets on the reference genome.

2.	 Can you visually see the benefit of using long reads compared to short reads for transcript sequencing?

3.	Look at the gene models of our reference genome. Is this species a good candidate for testing spliced-aware aligners?

4.	Can you identify a negative effect this has on mapping of the short reads?

## ChIP-seq mapping

Whole genome and transcriptome sequencing are probably the most common types of sequencing that is being performed. However, sequencing is used in a lot of other applications. We have discussed Hi-C sequencing during the assembly lecture, where we use DNA-DNA ligation and sequencing to determine which parts of the DNA are close to each other. 

Sequencing is also often used to determine the binding of certain biomolecules to the DNA or RNA. For example, Ribo-Seq is used to determine where on the mRNA molecules the ribosomes tend to bind, which can be used to determine the start and end of translation. Here, we will have a look at a ChIP-seq (ChIP = chromatin immunoprecipitation). With ChIP seq, we investigate the binding of certain protein (e.g. a transcription factor) to the DNA. In short, the proteins of interest are added to the DNA, and “glued” to their binding sites. We can then capture the proteins-DNA complexes, and extract the pieces of DNA that are stuck to the protein. These pieces are then sequenced using general sequencing methods. More info on ChIP-seq: https://en.wikipedia.org/wiki/ChIP_sequencing. Nowadays, a better version, called ChIP-exo has been developed, which gives better results. Since we are only sequencing the binding sites, the fragments are very short. This means we can use very short read sequencing, and only need to do forward reads.

Go to [SRA](https://www.ncbi.nlm.nih.gov/sra), and look up information on the ChIP-seq dataset we will be using (SRR1951324). 
 
1.	Which protein/transcription factor is investigated in this experiment?

Open fastp QC report  of the ChipSeq dataset created at the beginning of the pratical, and look at the duplication rate.

2.	Why is the duplication rate so high you think?

As we are dealing with genomic DNA, we will be mapping the reads using SMALT. We already have the genome index, so we can go straight to the mapping step:

	smalt map -n 2 SC_smalt_index ChipSeq.trimmed.fastq \
	| samtools view -bS | samtools sort > ChipSeq.sorted.bam
	samtools index ChipSeq.sorted.bam

Again, we will calculate the mapping statistics:
	
	samtools depth -a ChipSeq.sorted.bam | awk '{c++;s+=$3}END{print s/c}'
	samtools depth -a ChipSeq.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	samtools flagstat ChipSeq.sorted.bam | grep "mapped ("

3.	Why is only a small portion of the genome covered? 

Now load the bam file into IGV, and zoom in until you see the read mappings.

4.	What type of pattern do you observe in the read mapping.

5.	How would you determine which gene is regulated by TAF1? 
(hint: compare the SHP1 gene (`ChrII`, around 112 000bp) with the PXP1 gene 
(`ChrV`, around 119 000 bp)

## Using the wrong tool

Picking the right tool for the job is important in bioinformatics. In this part we will try using some mappers that are not designed for the job, and look at the results.

First, we will try mapping genomic reads using an RNA-aligner. We will use the WGS_PE set, and map it to the genome using STAR. We still have the genome index from last time, so we can start directly with the mapping:

	STAR --runThreadN 2 --genomeDir ./genome --sjdbGTFfile  sacCer3_refseq.gtf --sjdbOverhang 100 --readFilesIn WGS_PE_1.trimmed.fastq WGS_PE_2.trimmed.fastq 
	samtools view -bS Aligned.out.sam | samtools sort > DNA_STAR.sorted.bam
	samtools index DNA_STAR.sorted.bam
	rm -r Aligned.out.sam Log.* SJ.out.tab _STARgenome
	samtools depth -a DNA_STAR.sorted.bam | awk '{c++;s+=$3}END{print s/c}'
	samtools depth -a DNA_STAR.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	samtools flagstat DNA_STAR.sorted.bam | grep "mapped ("

1.	Compare the mapping statistics with the SMALT mapping. Are there any big differences?

Load the mapping in IGV, together with the original WGS_PE_smalt mapping.

2.	You will see a big difference in the mappings, what is this difference, and why is it there?

Now we will try the opposite. We will try to map the paired RNA-seq reads using a genomic mapper. 

3.	Given the genome’s gene structure, do you think there will be a big problems with the mapping?

Again we already have the index ready, so we can just do the mapping:

	smalt map -n 2 SC_smalt_index RNA_PE_1.trimmed.fastq RNA_PE_2.trimmed.fastq \
	| samtools view -bS | samtools sort > RNA_smalt.sorted.bam
	samtools index RNA_smalt.sorted.bam
	samtools depth -a RNA_smalt.sorted.bam | awk '{c++;s+=$3}END{print s/c}'
	samtools depth -a RNA_smalt.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	samtools flagstat RNA_smalt.sorted.bam | grep "mapped ("
	
As you will see, the mapping statistics are very similar to when we used the RNA-mapper. We will compare the mappings in IGV, so load in both the RNA_smalt and RNA_PE mappings in IGV.

4.	Was your prediction in question 3 accurate? Why (not)?

5.	What would happen if we would do this experiment in the human genome? 

We have also investigated long reads in this practical. Let’s find out what happens when we try mapping long reads to the genome using bwa-mem2:

	bwa-mem2 mem -t 2 SC_BWA_index WGS_PacBio.trimmed.fastq > PacBio_BWA.sam

Sadly, it looks like the software quickly crashes. You could try Smalt instead, but you would quickly realize that the mapping would take forever. Running this mapping on a cluster using 32 cores, will make the software crash after ~20 minutes. So, it is obvious that short read mappers can not handle long reads very well…

Lastly, let’s try the opposite. We will try mapping the WGS_SE reads using minimap2:

> Minimap2 actually has a short-read mapping mode, but for this experiment we will try mapping it in long read mode

	minimap2 -ax map-ont sacCer3.fa WGS_SE.trimmed.fastq | samtools view -bS | samtools sort > WGS_SE_minimap.sorted.bam
	samtools index WGS_SE_minimap.sorted.bam 
	samtools depth -a WGS_SE_minimap.sorted.bam | awk '{c++;s+=$3}END{print s/c}'
	samtools depth -a WGS_SE_minimap.sorted.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
	samtools flagstat WGS_SE_minimap.sorted.bam | grep "mapped ("

The statistics are still looking good, so let’s have a look in IGV, and compare to the Smalt mapping.

7.	How do the mappings compare? 

