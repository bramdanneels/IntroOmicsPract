# BINF201 – Practical 3 – Genome Assembly

In this practical we assemble a bacterial genome using different assembly software and methods, and compare the results.

## Software installation and data retrieval

In this tutorial, we will have a look at the following assembly (or related) software:

- [SPAdes](https://github.com/ablab/spades) - Likely the most popular De Bruijn Graph assembler, used very commonly for bacterial or small eukaryotic genome assembly wih short reads
- [ABySS](https://github.com/bcgsc/abyss) - Another short-read assemblyer using De Bruijn Graphs
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

## Quality Control

Before assembling the genome, it is important to take a look at our data. Go back to the QC environment created in practical 2: `mamba activate QC`.
Then, run FastQC on both files, and look at the reports.

<details>
<summary>How many reads are in this dataset, and how long are the reads?</summary>
_There are 387 568 reads in the dataset, all of them 150bp long._
</details>

<details>
<summary>Based on the report, is adapter and/or quality trimming necessary?</summary>
_No: The quality scores are all above Q30, and there are no adapters detected in both read files_
</details>

<details>
<summary>Given that the genome size of _M. genitalium_ is +- 580 kbp, what coverage do we expect to have with our read set?</summary>
_200x coverage_

_To caculate this, we need to first calculate the total number of bases in our read set. We have two read sets with 387568 reads of 150bp. Thus, the total number of bases is $2*387568*150$ = 116 270 400bp._
_Since we the genome is 580kbp (580 000bp) long, that means we should have an average coverage of $116270400/580000$ or 200x coverage._
</details>

## Genome assembly using Spades

We will assemble the genome from the reads using three different assemblers. The first one, [SPAdes](https://github.com/ablab/spades), is one of the most popular genome assemblers for small genomes (viral, bacterial, yeast). It is also very popular for doing metagenome assembly. 

We will run SPAdes in the "isolate" and the "careful" mode as follows (remember to switch back to the assembly environment first: `mamba deactivate QC && mamba activate Assembly`):
> The `&&` in the above command tells the command line "Do the first command (deactivate), and if that succeeds, run the second command (activate).

```
spades.py -o spades_isolate -1 MycGen_1.fastq -2 MycGen_2.fastq --isolate
spades.py -o spades_careful -1 MycGen_1.fastq -2 MycGen_2.fastq --careful
```

The assembly can take some minutes to complete. In the meantime, you can try answering following questions regarding genome sequencing and assembly:
> You can find the necessary information in the [SPAdes manual](https://ablab.github.io/spades/) and/or by searching on the internet.

<details>
<summary>What do the "isolate" and "careful" options do?</summary>
According to the SPAdes manual:
- _The `--isolate` option: This flag is highly recommended for high-coverage isolate and multi-cell Illumina data; improves the assembly quality and running time. We also recommend trimming your reads prior to the assembly._
- _The `careful` option : Tries to reduce the number of mismatches and short indels. Also runs MismatchCorrector - a post processing tool, which uses BWA tool (comes with SPAdes). This option is recommended only for assembly of small genomes. We strongly recommend not to use it for large and medium-size eukaryotic genomes._
</details>

<details>
<summary>SPAdes will create both contigs and scaffolds. What is the difference between both?</summary>

- _Contigs are contigues sequences that are build from overlapping reads_
- _Scaffolds are a combination of contigs in a certain order and orientation_
</details>

<details>
<summary>What is mate pair sequencing, and how is it similar and/or different from classic paired end sequencing?</summary>
_Mate Pair sequencing is a type of paired end sequnecing where the size of the DNA fragments is significantly larger. Mate pair sequencing is often used in combination with normal paired end sequencing to improve difficult genome assemblies, and to generate longer scaffolds from the contigs._
</details>

Once SPAdes is done running, it will have created two new directory (`spades_isolate` and `spades_careful`) containing the intermediate and final output files. The files we are interested in are the `contigs.fasta` and `scaffolds.fasta` files. 
We will count the number of total sequences in each file using the following commands:

```
grep -c ">" spades_*/contigs.fasta spades_*/scaffolds.fasta
```
> `grep` is a program that performs text searches. Here we tell `grep` to search for the `>` sign in the files we give it (all `contigs.fasta` and `scaffolds.fasta` files in folders starting with `spades_`). 
> The `-c` option tells grep we just want to count the number of times the `>` symbol appears. This command works because every sequence in a `.fasta` file starts with a single `>` symbol, and there can't be `>` symbol in the actual sequence.

<details>
<summary>Are there any differences in the number of contigs or scaffolds between both runs?</summary>
_Yes: the `careful` run has slightly more contigs and scaffolds (57 vs 54/53). There is no difference between contigs and scaffolds in the `careful` run, but there is one fewer scaffold than contigs in the `isolate` run_
> Note: the exact numbers you have may vary from the results given here.
</details>

We will have a more detailed look at these assemblies later in the practical

## Genome assembly using Abyss

[Abyss](https://github.com/bcgsc/abyss) is a genome assembler that can be used for assemblies of all sizes using short reads. We will run it here using the default settings.
As you might have noticed, SPAdes tried assembling the reads using different values for _k_, and then picks the best one. 
ABySS does does not have this capacity, so we will have to rune ABySS by specifying our own _k_-values. We will run ABySS here using a _k_-value of 31 and 75:

```
mkdir -p Abyss_k31 Abyss_k75
abyss-pe k=31 name=abyss_k31  B=1G in="MycGen_1.fastq MycGen_2.fastq"
abyss-pe k=75 name=abyss_k75  B=1G in="MycGen_1.fastq MycGen_2.fastq"
mv abyss_k31* Abyss_k31
mv abyss_k75* Abyss_k75
```
> The `k` option sets the _k_-mer length, the `B` option sets the size of the Bloom filter (a specific datastructure ABySS uses to store the De Bruijn Graph).

Once both assemblies are finished, the final contig files will be stored in the output directory, as `abyss_kXX-contigs.fa`. 
Use `grep` again to find the number of contigs and scaffolds in the two Abyss assemblies.

<details>
<summary>Which of the two ABySS assemblies gives the most contigs? Do you think this is better or worse?</summary>
_The assembly using a k-value of 31 gave a lot more contigs: 570 vs 49_

_Assuming the total assembly length is the same, having more contigs is worse than having fewer contigs, as more contigs means the assembly is more fragmented (more, but smaller contigs)_
</details>

<details>
<summary>Are there differences between scaffolds and contigs in the assemblies?</summary>
_Yes: In general there are fewer scaffolds than contigs (549 vs 570 in k31; 20 vs 49 in k75)_
</details>

<details>
<summary>Why do we observe such a large different in the number of contigs between k31 and k75?</summary>
_Using shorter k-mers (lower k-value) means using shorter sections of our reads to find overlaps. This means it is more difficult to solve repeats using shorter k-mers, leading to higher fragmentation._
</details>

## Genome assembly using Megahit

Lastly, we will use [Megahit](https://github.com/voutcn/megahit) to assemble the reads. Megahit is one of the faster assemblers, and was originally mainly designed for metagenomes. However, it can also be used for single genomes or single-cell assemblies. Similar to the other assemblers, we will run Megahit using the standard settings:
Similar to SPAdes, Megahit will assemble using multiple _k_-values, and pick the best assembly from all of them.

```
megahit -1 MycGen_1.fastq -2 MycGen_2.fastq -o megahit_assembly
```

Megahit does not create scaffolds, and will output the contigs in a file called `final.contigs.fasta`

<details>
<summary>How many contigs are in the megahit assembly?</summary>
_22_
</details>

## Assembly QC

We have now created multiple assemblies using different assemblers and methods. 
Now, we will use [Quast](https://github.com/ablab/quast) to assess and compare assembly statistics between assemblies.
Since there exists a public genome sequence for _M. genitalium_, we will also compare our assemblies to that one.

First, we need to download the _M. genitalium_ reference genome:

```
wget -O ref_genome.zip "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000027325.1/download?include_annotation_type=GENOME_FASTA"
unzip ref_genome.zip
mv ncbi_dataset/data/GCF_000027325.1/GCF_000027325.1_ASM2732v1_genomic.fna ./reference.fasta
rm -r ncbi_dataset/ README.md ref_genome.zip
```

Then we can run Quast on all our assemblies (using the scaffolds file if available) using the following command:

```
quast -o assembly_QC -r ./reference.fasta spades_*/scaffolds.fasta Abyss_k*/*scaffolds.fa megahit_assembly/final.contigs.fa
```

This will create a folder called "assembly_QC" containing the report on the statistics (`report.html`). Go through the report and answer the questions below:
> You can click "extend report" to get more statistics

<details>
<summary>Which of the assemblies gave the largest contig? Which one gave the largest scaffold?</summary>
_The SPAdes assemblies: 417 388bp (you can find this in "Largest contig" in the "Statistics without reference" section._
</details>

<details>
<summary>Which assembler had the most contigs? Which one the fewest?</summary>
_The ABySS assembly using k31 has the most contigs (549), the ABySS assembly using k75 has the fewest contigs (20)._
_You can find this under "# contigs (>= 0bp)" in the "Statistics without reference" section._
</details>

<details>
<summary>Which assembly has the best N50? Is this the highest or lowest value?</summary>
_The SPAdes assemblies have the highest N50: 417 388_

_This is the highest value, as a higher N50 is better. This is because a higher N50 means that the largest contigs (making up at least 50% of the assembly) are larger compared to assemblies with lower N50._

_The N50 value can be found in the "Statistics without reference" section.
</details>

<details>
<summary>Compare the size of the largest contig with the N50 of the SPAdes assembly. Is there anything that stands out? How would you explain this?</summary>
_Yes: they are the same. This is because the largest contig (417kbp) takes up over half of the assembly (577kbp). Thus, half of the assembly is contained in the largest contig, which is 417kbp._
_This explains also why the L50 is 1: you only need 1 contig to get at least 50% of the assembly.
</details>

<details>
<summary>In the "Genome statistics" section, there is a metric called NG50. How does this differ from N50? Which metric would you prefer to use to assess the quality of an assembly?</summary>
_NG50 is relative to the reference genome size, while N50 is relative to the assembly size. Using NG50 is preferred, because it avoids bias because of incomplete assemblies (see example below).

_Example: We have two assemblies of a genome:_

- _AssemblyA has 12 contigs of each 250kbp_
- _AssemblyB has 2 contig of each 500kbp, and 25 contigs of each 1kbp_

We know from previous studies that the genome of this species should be around 2.5 Mbp.
Using this information we can calculate the N50:

_The N50 of AssemblyA is 250kbp (half of the assembly is in contigs of size 250kbp or larger), the N50 of AssemblyB is 500kbp (half of the assembly is in contigs of size 500kbp or larger)_

_Based on N50 alone it looks like Assembly is better (higher value). However, let's calculate NG50 now. This means we have to see the size of the contigs that we need to have 50% of the reference genome (= 1.25 Mbp)_200x

_The NG50 of AssemblyA is still 250kbp (we can build an assembly of 1.25 Mbp using contigs of 250kbp or higher). However, the NG50 of AssemblyB is now 1kbp: We need all contigs to get to an assembly of 1.25 Mbp, thus the NG50 is equal to the smalles contig (1 kbp)_

_Thus, NG50 is a better metric to assess an assembly, but it requires you to have a reference genome, or to know how large the genome should be.
</details>


<details>
<summary>Which assembly covers the largest fraction of the reference genome?</summary>
_The Abyss_k75 assembly (99.289%). You can find this in the "Genome fraction" row under the "Genome statistics section"._
</details>

<details>
<summary>Which assembly shows the most mismatches and indels compared to the reference?</summary>
_The Abyss k75 assembly has the most mismatches and indels, the Abyss k31 the fewest._
</details>

<details>
<summary>Given that the read quality was very good, what could be the reason for observing mismatches and indels in our assemblies, compared to the reference?</summary>
_Because we are sequencing a different strain than the reference genome. Bacteria, and especially pathogenic bacteria, evolve very quickly. Thus, if we sequence a bacterium of a certain species, it is very unlikely that we find the exact same genome as the reference._
_The observed mismatches and indels are thus likely differences that evolved between the reference strain, and the strain that was sequenced in our dataset._
</details>

<details>
<summary>Which of the genome assemblies do you think is best, and why?</summary>
_Either one of the SPAdes assemblies, or the ABySS k75. The SPAdes assemblies have a higher NG50 and have a larger largest contig._ 
_However, ABySS k75 captures a higher fractino of the reference genome, an has fewer contigs in total._
_Personally, I would pick the SPAdes assembly, but filter the smaller contigs (< 1000bp)_
</details>

## BUSCO analysis

In this case we compared our assemblies to the reference genome to assess if we had a good assembly or not. 
Of course, when sequencing a new species, there will be no reference genome available to compare with, so it gets more difficult to properly assess assembly quality. 
In that case, tools like [BUSCO](https://busco.ezlab.org/) come in handy. BUSCO is a tool that will try to detect the presence of conserved genes in your assembly.
The BUSCO tool uses a database which contains sets of conserved genes (universal single-copy orthologs) which are present in the major of species of a certain lineage (e.g. plants).
By checking how much of the genes that should be in your assembly, you can have a rough idea of how complete your assembly is.

We will run BUSCO on the SPAdes-isolate assembly, the ABySS-k75 assembly, and the ABySS-k31 assembly. The most important step to running BUSCO is figuring out what lineage to use.
To get an overview of which lineages are available, you can run the following commands:

```
busco --list-datasets
```

This will give you a list of lineages that you can use. But how do we find the correct lineage? There are 3 options:

- You look up the taxonomy (e.g. on [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), and check which of the available groups is applicable
- Let BUSCO figure out the best lineage (using the `--auto-lineage` option. However, this will increase run time and doesn't always work out.
- Only use a general lineage (e.g. Bacteria or Eukaryota) 

Here, we will let BUSCO try to figure out what we are dealing with:

```
busco -m geno -i spades_isolate/scaffolds.fasta --auto-lineage-prok
```
> The `-m` option selects the running mode (here `geno` for genome, other options are transcriptome or protein).
> Since we know that we are dealing with a bacterium, we can constrain the automatic lineage selection to only search for prokaryotes (`--auto-lineage-prok`). 

Have a look at the results, which are handily printed to the screen. Have a look at the results. For this tutorial, we will only be looking at the proportion of complete BUSCOs (C).
If all goes well, it should have identified the genome as bacterial (the generic domain), and as belonging to the mycoplasmatales dataset.

<details>
<summary>Compare the results for the bacteria dataset with the results for the mycoplasmatales. Why could there be such a big difference?</summary>
_BUSCO only finds 55% of conserved bacterial genes, but 97% of mycoplasmatales genes. This is because Mycoplasmoides is very different from normal bacteria. It has a small genome (580 kbp) compared to most bacteria (who are mostl around 3-5 Mbp), and is pathogenic._
_This means that this bacterium is very specialized, and doesn't have many of the "general" genes that other bacteria have". However, if we only look at conserved genes within the Mycoplasmatales, it does have most genes that we expect._
_This is because most species in the Mycoplasmatales are specialized in the same way, and thus share a lot more genes._

_As a general rule, a BUSCO score >90% is considered good, >95% is considered very good._
</details>

Now try running BUSCO yourself on the two ABySS assemblies, but specifying that you only want to run it on the `mycoplasmatales_odb10` lineage. You can run `busco -h` to print a help message on how to run BUSCO and use the different options.

<details>
<summary>Are there large differences in BUSCO score between the assemblies?</summary>
_There are no large differences, but the BUSCO score is slightly lower in the ABySS-k31 assembly (96.6% vs 97.7%)._
</details>

## Long Read Assembly

Long reads sequencing is becoming more and more common, even for small genomes. The range of tools used for assembling long reads is different than the ones we have seen above.
Many of these tools take a bit longer to run, so some patience might be required. We'll also use a trick that will allow the process to run in the background so you can logout and pick up where you left off later.

We will assemble Nanopore and PacBio reads from _Mycoplasma ovipneumoniae_, using the 3 different tools: [Canu](https://github.com/marbl/canu), [Flye](https://github.com/mikolmogorov/Flye), and [HifiAsm](https://github.com/chhylp123/hifiasm).

First, we'll assemble the Nanopore reads using Canu. Canu needs an estimation of the genome size to be able to calculate expected coverage. Luckily, there are already [a lot of _M. ovipneumoniae_ genomes](https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=29562) sequenced.
Based on these genome assemblies, we can say that the expected genome size is roughly 1.1 Mbp.

Below you can find the commands for running the assemblies normally. If you are running out of time, or you want to just launch the assemblies and look at the results later, we can ask the server to run the command in the background.
To do this, we will be using a combination of acommand called `nohup` and the background command `&`. 
The `nohup` command will tell the server to not stop a command on "hangup" (i.e. when disconnecting from the server), and the `&` command will tell the server to run the process in the background so we can do other stuff in the meantime (like log out from the server).
To run a command using `nohup` and `&`, simply put `nohup` in front of the command, and `&` at the end of the command. The output of your command will be saved to a file called `nohup.out`.
> For example, the command `fastqc sample.fastq` would become `nohup fastqc sample.fastq &`

```
canu -p canu -d MycOvi_Canu genomeSize=1.1m maxThreads=4 -nanopore MycOvi_Nano.fastq
```
> `-p` determines the prefix (so in our case all outputfiles will start with "canu"); `-p` sets the name of the output directory

We can now find the assembly in `MycOvi_Canu/canu.contigs.fasta`.

Now let's assemble the same reads using Flye:

```
flye --nano-raw MycOvi_Nano.fastq --out-dir MycOvi_Flye --threads 2
``` 

The final assembly can be found in `MycOvi_Flye/assembly.fasta`

Lastly, we will assemble the PacBio HiFi reads using both Flye and HifiAsm:
> HifiAsm is actually designed for diploid organisms. It might not be the most suitable tool for assembling a bacterial sequence.

```
flye --pacbio-hifi MycOvi_HiFi.fastq --out-dir MycOvi_HiFi_Flye --threads 2
hifiasm -o hifiasm -t2 MycOvi_HiFi.fastq
awk '/^S/{print ">"$2;print $3}' hifiasm.bp.hap1.p_ctg.gfa > hifiasm.bp.hap1.p_ctg.fasta
mkdir -p MycOvi_HiFiAsm
mv hifiasm* MycOvi_HiFiAsm/
```
> Since hifiasm doesn't output `fasta` directly, only `.gfa`, we use the `awk` command to convert from one to the other format.

The assemblies can now be found in `MycOvi_HiFi_Flye/assembly.fasta` and `MycOvi_HiFiAsm/hifiasm.bp.hap1.p_ctg.fasta`

Now run Quast on the 4 assemblies and comapre their statistics (we will not compare to a reference this time).

<details>
<summary>Are these assemblies better or worse than the assemblies we have made using short reads? (the expected genomes size is +- 1 Gbp)</summary>
_A lot better. Most assemblers manage to form one big contig close to the expected genome size._
</details>

Run a BUSCO analysis on one assembly using the Nanopore reads, and one assembly using the PacBio HiFi reads

<details>
<summary>Do the BUSCO scores indicate a good assembly?</summary>
_For the HiFi reads: yes (99.4% completeness for mycoplasmatales)_

_For the Nanopore reads: no (23-26% completeness for mollicutes)_

_Thus, while both Nanopore and HiFi assemblies had very good statistics, the Nanopore ones have very low BUSCO completeness score, while the HiFi showed very good BUSCO scores._
_This is because the Nanopore reads used old Nanopore technology, that had error rates up to 15-20%. Assembling a genome using these reads alone leads to assemblies full of mitakes._
_Because of that, Nanopore was often combined with short reads to fix their mistakes. The most recent Nanopore technologies has a significant lower error rate (<5%), and can now create good assemblies on its own as well._
</details>
