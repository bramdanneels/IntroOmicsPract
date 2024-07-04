# BINF201 – Practical 3 – Genome Assembly

For this practical, you will need to have installed and configured Conda. If this is not the case, see [here](https://docs.conda.io/en/latest/miniconda.html) for installation instructions for your system.

> For Windows users, it is preferable to use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install), and install conda on there using the Linux install instructions.

In this practical, you will perform the genome assembly of *Mycoplasma genitalium*, a pathogenic bacterium. It has a very small genome, even for bacteria (+- 580 kbp). 
We will use a dataset of paired-end Illumina sequencing reads (2x150bp). The reads are publicly available in [SRA](https://www.ncbi.nlm.nih.gov/sra), with run accession ERR486840.

## Tool installation and data retrieval

First, we need to install the necessary tools. Similar to last practical, we will create a new conda environment containing the necessary tools.

    conda create -n assembly -c bioconda -c conda-forge spades abyss megahit quast wget gzip unzip curl
Then, activate the correct environment using:

    conda activate assembly

Next, create a new directory called “assembly” (`mkdir assembly`), enter it (`cd assembly`), and download and decompress  the data as follows:

	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_1.fastq.gz
	wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR486/ERR486840/ERR486840_2.fastq.gz
	gunzip *gz

## Quality Control

Before assembling the genome, it is important to take a look at our data. Go back to the QC environment created last practical (use `conda deactivate` first, then activate the QC environment: `conda activate QC`). 
Then, run FastQC on both files, and look at the reports.

1.	How many reads are there in the dataset?

2.	Based on the FastQC reports, do we need to do any trimming (adapter and/or quality)?

3.	How much base coverage do you expect, on average, for this dataset? (hint: you will need the estimated genome size (580 kbp), the number of reads, and the average read length  to calculate this)

## Genome assembly using Spades

We will assemble the genome from the reads using three different assemblers. The first one, [SPAdes](https://github.com/ablab/spades), is one of the most popular genome assemblers for small genomes (viral, bacterial, yeast). It is also very popular for doing metagenome assembly. 

We will run SPAdes in “careful” mode as follows (remember to switch back to the assembly environment first!):

	spades.py -o spades_assembly -1 ERR486840_1.fastq -2 ERR486840_2.fastq --careful

1. What does the “careful” setting do?

The assembly can take some minutes to complete. In the meantime, you can try answering following questions regarding genome sequencing and assembly:

2.	What is the difference between contigs and scaffolds?

3.	What is the difference between paired-end sequencing and mate-pair sequencing?

Once SPAdes is done running, it will have created a new directory (`spades_assembly`) containing the intermediate and final output files. The files we are interested in are the `contigs.fasta` and `scaffolds.fasta` files. 
We can try counting how many sequences are in each file using grep:
	
	grep -c ">" contigs.fasta
	grep -c ">" scaffolds.fasta 

> Here we tell `grep` to count (`-c`) all lines which contain a ">" character

4.	Are there differences in the number of contigs between both files? Is this expected based on the type of data we have available?

## Genome assembly using Abyss

[Abyss](https://github.com/bcgsc/abyss) is a genome assembler that can be used for assemblies of all sizes using short reads. We will also run it using the default settings (remember to go back to the folder containing the reads first: `cd ..`). 
As Abyss doesn’t automatically create the output folder, we have to create it manually first:
	
	mkdir abyss_assembly

Another issue, is that Abyss needs absolute file paths to be specified. We can save the absolute path of our current (working) directory in a variable, and use that in the Abyss call:

	export workingdir=$(pwd)
	abyss-pe k=31 name=abyss_k31 in="${workingdir}/ERR486840_1.fastq ${workingdir}/ERR486840_2.fastq" -c abyss_assembly_k31

> This assembly might take some time, so don’t worry if it seems like nothing is happening for a while. 

Since Abyss works with a fixed *k*-value, we will run a second Abyss assembly, but this time using a higher value of k:
	
	abyss-pe k=77 name=abyss_k77 in="${workingdir}/ERR486840_1.fastq ${workingdir}/ERR486840_2.fastq" -c abyss_assembly_k77

Once both assemblies are finished, the final contig files will be stored in the output directory, as `abyss_kXX-contigs.fa`. 
Use `grep` again to find the number of contigs in the two Abyss assemblies.

1.	Which of both assemblies has the most contigs?

2.	What could be the reason of this difference?

## Genome assembly using Megahit

Lastly, we will use [Megahit](https://github.com/voutcn/megahit) to assemble the reads. Megahit is one of the faster assemblers, and was originally mainly designed for metagenomes. However, it can also be used for single genomes or single-cell assemblies. Similar to the other assemblers, we will run Megahit using the standard settings (think of going back to the assembly folder containing the reads):

	megahit -1 ERR486840_1.fastq -2 ERR486840_2.fastq -o megahit_assembly

The output we care about for this practice is `final.contigs.fa` in the `megahit_assembly` folder.

## Assembly QC

In this part we will use [Quast](https://github.com/ablab/quast) to assess and compare assembly statistics between assemblies, and compare to the published genome assembly of *Mycoplasma genitalum*. 

First, we’ll create a new folder for the QC:

	mkdir assemblyQC
	cd assemblyQC

Then, we’ll download the reference genome from NCBI:

	curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000027325.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000027325.1.zip" -H "Accept: application/zip" > NCBI.zip
	unzip NCBI.zip
	mv ncbi_dataset/data/GCF_000027325.1/GCF_000027325.1_ASM2732v1_genomic.fna ./reference.fasta
	rm -r ncbi_dataset NCBI.zip

Then, we’ll make a copy of all our created assemblies in our new directory:

	cp ../spades_assembly/contigs.fasta ./spades.fasta
	cp ../abyss_assembly_k31/abyss_k31-contigs.fa ./abyss_k31.fasta
	cp ../abyss_assembly_k77/abyss_k77-contigs.fa ./abyss_k77.fasta
	cp ../megahit_assembly/final.contigs.fa ./megahit.fasta

Now we can run Quast to compare to assemblies to each other and to the reference:

	quast -o quast_out -r reference.fasta spades.fasta abyss_k31.fasta abyss_k77.fasta megahit.fasta
Now, open the `report.html` file in the `quast_out` folder, and look at the results.

1.	Which of the assemblers gave the largest contig?

2.	Which of the assemblers yielded the highest number of contigs? Is this good or bad?

3.	Which assembly has the best N50? Is this the highest or lowest value, and why?

4.	Which assembly covers the largest fraction of the reference genome?

5.	Which assembly shows the most mismatches and indels compared to the reference?

6.	Considering the read QC we did, what could be the reason of that many mismatches and indels in that assembly?

7.	Which genome assembler would you prefer, and why?

## BUSCO analysis

In this case we compared our assemblies to the reference genome to assess if we had a good assembly or not. Of course, if you sequence a new species, there will be no reference genome available to compare with, so it gets more difficult to properly assess assembly quality. In that case, tools like [BUSCO](https://busco.ezlab.org/) come in handy. BUSCO checks the presence of lineage specific genes in your assembly. 

We will run BUSCO on the SPAdes assembly, using the general Bacteria dataset containing bacterial core genes. Since BUSCO doesn’t always play nice with other software tools, we will first install it in a separate environment:

	conda create -n busco -c bioconda -c conda-forge busco=5.4.7
	conda activate busco
>Sometimes conda can stall on installing busco, as it has a lot of dependencies to resolve. If conda seems to be stuck on the ‘resolving environment’ step, consider cancelling the installation (`Ctrl+C`), and using the workaround described at the end of the document (“Using Mamba instead of Conda)

Then we can run BUSCO on the SPAdes assembly:

	busco -i spades.fasta -o spades_busco -m genome -l bacteria

> We tell Busco to run in genome mode (`-m genome`), and use the bacterial lineage genes (`-l bacteria`)

Now look at the output (on the terminal, or in the `spades_busco` file

1.	How many Busco genes were checked in total?

2.	What proportion of the conserved bacteria genes is present in the assembly? How many are totally missing?

Look up what kind of bacterium *Mycoplasma genitalium* (syn. *Mycoplasmoides genitalium*) is

3. Can you explain the BUSCO results based on the lifestyle of this bacterium?

One of the strengths of BUSCO is that it has many lineage-specific datasets. Using the general Bacteria dataset is thus not always advised. A list of lineages is available here: https://busco-data.ezlab.org/v5/data/lineages/

*M. genitalium* is part of the Mycoplasmatales, for which there is a dataset available. We’ll try running BUSCO with this set instead:

	busco -i spades.fasta -o spades_busco_myco -m genome -l mycoplasmatales

Look at the results of this new analysis

4.	How did the results change compared to the Bacteria analysis? Are there results better or worse?

Figuring out what lineage your organism belongs to might be tricky sometimes. Luckily, BUSCO has an option to try to automatically try to detect what lineage fits best for your data:

	busco -i spades.fasta -o spades_busco_auto -m genome --auto-lineage-prok

> While the auto-lineage feature is really handy, it sometimes fails to resolve leading to a very long running time

If all went well, BUSCO should have identified the correct lineage for this species (Mycoplasmatales). However, always use this option with caution, and verify the results. If you are investigating a species with no close lineage in the BUSCO dataset, the results might end up somewhat weird!

## Using Mamba instead of Conda

Sometimes conda gets stuck wen installing certain packages, especially when a lot of dependencies need to be resolved. That’s why [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html) has been developed. Mamba functions the same as Conda, but handles package installation a lot better. To use mamba instead of conda, just install mamba using conda first:

	conda install mamba
Once that is completed, you can now use mamba install or mamba create to install packages and environments instead of conda. However, activating and deactivating environments still uses conda! 

For example for installing the busco environment:
	
	mamba create -n busco -c bioconda -c conda-forge busco=5.4.7
	conda activate busco
	busco -h
