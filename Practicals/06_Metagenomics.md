# BINF201 – Practical 6 – Metagenomics

In this practical we will be having a look at metagenomic data. We will perform metagenome assemblies and perform binning to try to recover full genomes.
In addition we will have a look at marker-gene based analysis (e.g. 16S) and how to analyse those.

## Software installation and data retrieval

In this tutorial, we will have a look at the following metagenomics and related software:

- [MegaHit](https://github.com/voutcn/megahit) - An ultra-fast and memory-efficient (meta-)genome assembler
- [MetaFlye](https://github.com/mikolmogorov/Flye) - An adaptation of the Flye long-read assembler for metagenomes
- [Concoct](https://github.com/BinPro/CONCOCT) - An unsupervised binner for metagenomic contigs
- [Kraken2](https://github.com/DerrickWood/kraken2) - A taxonomic k-mer based taxonomic classifier
- [Metaphlan](https://huttenhower.sph.harvard.edu/metaphlan/) - A metagenomic community profiling tool
- [CheckM](https://github.com/Ecogenomics/CheckM) - A quality assesment tool for mircrobial genomes
- [KrakenTools](https://github.com/jenniferlu717/KrakenTools) - Toolkit for anlysing Kraken/Bracken output
- [Krona](https://github.com/marbl/Krona/releases) - Taxonomy visualisation tool
- [Bracken](https://github.com/jenniferlu717/Bracken) - Statistical method for species abundance calculation

### For students using NREC

The data and software has been setup on the NREC server. Before starting the practical, make sure to activate the correct environment before each part of the tutorial 
(`QC` for the QC part, `Metagenome` for the metagenome part)

You can make a copy of the data you will be working on by running this command from your home directory:
> If your not sure if you are in the home folder, run `cd` or `cd ~` to go to your home directory.

```
mkdir -p Practical6
ln -s /storage/06_Transcriptome/* Practical6/
```
> `mkdir -p` creates a folder called Practical6. The "-p" options tells mkdir to create subdirectories if necessary, and to not give an error if the folder(s) already exist
> `ln -s` creates what we call a "symbolic link". This creates a small file that just says "Instead of this file, use the file that I'm liking to". 
> This allows you to "copy" files without actually having to make a physical copy.

Now, go to the newly created directory (by running `cd Practical6`), and you are ready to start!

### For students running on their own pc

You will first have to setup the correct environment with the necessary tools. 
See [the intro practical](00_IntroSetup.md) on how to install mamba, and how to create an enviroment for downloading the necessary data.

To create a new environment including the necessary tools for this practical, run the following command:

```
mamba create -n Metagenomics -c bioconda -c conda-forge spades kraken2 krona metaphlan checkm-genome krakentools hclust2 concoct smalt quast bracken flye minimap2 flash2 pylablib
```
> This will create a new environment called "Transcriptome" with the necessary tools.
> The ncbi-datasets-cli will allow us to download our reference genome.

This will print some information to the screen, telling us to do some additional setup:

- Get taxonomy information for Krona:

```
ktUpdateTaxonomy.sh
```

- Set CheckM root:

```
checkm data setRoot {path to CheckM data}
```
> Check the output on your screen for the exact comment

- Download necessary kraken2 databases:

```
kraken2-build --standard --threads 2 --db KrakenStandardDB
kraken2-build --SILVA --threads --db KrakenSILVADB
bracken-build -d KrakenSILVADB -t 2 -k 35 -l 300
```
>The first database can take up 100Gb of space. If you don't have that kind of space, only download the SILVA one

- Configure Autometa databases: follow the instructions [here](https://autometa.readthedocs.io/en/latest/databases.html).

Then, create a new folder (e.g. `Practical6`), go into it (`cd Practical6`), and then download the necessary data, by either:

- Download directly from the [Zenodo repository]():

```
wget link_to_zenodo_file(s)
gunzip *gz
```

- Download and filter the data manually using the commands below:

> Remember to activate the download environment (see [here](00_IntroSetup.md)).
> This will download the data, and create subsamples
> Note: Largest downloaded file: 7G, total size after subsampling: 7.5G

> Important: Part of the data-prep is mapping the reads to Chr15. It is not recommended to do this on a normal desktop computer.
> Running the commands as they are written (i.e. with only 2 cores) might take multiple hours to complete the mappings!
> If you don't have access to high-performance computing resources, we highly recommend downloading the prepared files from Zenodo (see above)

```
prefetch SRR27874410
fasterq-dump -p --outdir ./ --split-files SRR27874410/SRR27874410.sra
seqtk sample -s666 SRR27874410_1.fastq 2500000 > Plant_F.fastq
seqtk sample -s666 SRR27874410_2.fastq 2500000 > Plant_R.fastq
rm -r SRR27874410*

prefetch SRR23255609
fasterq-dump -p --outdir ./ --split-files SRR23255609/SRR23255609.sra
seqtk sample -s666 SRR23255609_1.fastq 2500000 > Stool_F.fastq
seqtk sample -s666 SRR23255609_2.fastq 2500000 > Stool_R.fastq
rm -r SRR23255609*

prefetch SRR25440994
fasterq-dump -p --outdir ./ --split-files SRR25440994/SRR25440994.sra
mv SRR25440994_1.fastq Soil_1F.fastq
mv SRR25440994_2.fastq Soil_1R.fastq
rm -r SRR25440994*

prefetch SRR25440995
fasterq-dump -p --outdir ./ --split-files SRR25440995/SRR25440995.sra
mv SRR25440995_1.fastq Soil_2F.fastq
mv SRR25440995_2.fastq Soil_2R.fastq
rm -r SRR25440995*

prefetch SRR25440996
fasterq-dump -p --outdir ./ --split-files SRR25440996/SRR25440996.sra
mv SRR25440996_1.fastq Soil_3F.fastq
mv SRR25440996_2.fastq Soil_3R.fastq
rm -r SRR25440996*

prefetch ERR12117715	
fasterq-dump -p --outdir ./ --split-files ERR12117715/ERR12117715.sra
seqtk sample -s666 ERR12117715.fastq 75000 > Sea_HiFi.fastq
rm -r ERR12117715*
```

## The Data

The data for this tutorial are different kinds of metagenomes:

- A plant metagenome (short paired-end reads)
- A stool (faeces) metagenome (short paired-end reads)
- An ocean metagenome (long PacBio HiFi reads)
- Three soil metagenomes (16S amplicon sequencing, paired Illumina reads)

Most samples have been downsampled to speed up their analysis.
For an overview of where the data comes from, see the [information sheet on Zenodo]().

## Shotgun metagenomics

In this first part, we will have a look at the shotgun metagenomics data. This data is generated by sequencing all DNA from a sample. We have two short-read samples, and one long-read one.

### Quality control

As is usual, we'll start by performing quality control of our sequences. Run FastQC on the Plant, Stool, and Sea samples, and have a look at the QC reports.
> Remember to activate the QC environment!

<details>
<summary>Do we need any trimming?</summary>

_Some trimming could be good. Adapter levels are generally low, but in some read sets the quality drops off near the end._
</details>

<details>
<summary>Look at the %GC content plots. Why do we observe some of these weird patterns?</summary>

_Because we are dealing with metagenome samples: the samples contain DNA from multiple organsims, which often have different %QC, leading to the erratic %GC plots._
</details>

We will do some basic trimming using cutadapt to get rid of low-quality bases and adapters in the short reads:
> We will override our read files with the trimmed versions. If for some reason you need the original reads again, rerun the link command at the beginning of the tutorial.

```
cutadapt -q 30,30 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tempF.fq -p tempR.fq Plant_F.fastq Plant_R.fastq
mv tempF.fq Plant_F.fastq
mv tempR.fq Plant_R.fastq
cutadapt -q 30,30 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tempF.fq -p tempR.fq Stool_F.fastq Stool_R.fastq
mv tempF.fq Stool_F.fastq
mv tempR.fq Stool_R.fastq
```

### Taxonomic profiling of reads

Now that we have QC’ed and filtered our reads, we could try a first classification on read level. 
We will first use [kraken2](https://ccb.jhu.edu/software/kraken2/) for this. 
Kraken extracts k-mers from our reads and compares them to a database to classify them. We will first concatenate our paired-end reads from one sample together:
> Here we will classify the trimmed reads. It is however debated wether it is better to classify trimmed or raw reads. Raw reads are more representative, but can lead to some errors due to low-quality bases.

```
cat Plant_F.fastq Plant_F.fastq > Plant_all.fastq
cat Stool_R.fastq Stool_R.fastq > Stool_all.fastq
```

Then we will Kraken, using a 16S database from [SILVA](https://www.arb-silva.de/):

```
kraken2 --db /storage/dbs/KrakenSILVADB/ --threads 2 --output Plant.kraken.out --report Plant.kraken.report Plant_all.fastq
kraken2 --db /storage/dbs/KrakenSILVADB/ --threads 2 --output Stool.kraken.out --report Stool.kraken.report Stool_all.fastq
kraken2 --db /storage/dbs/KrakenSILVADB/ --threads 2 --output Sea.kraken.out --report Sea.kraken.report Sea_HiFi.fastq
```

Take a look at the generated reports:

```
less *report
```
> This will open all files endin in "report" using the `less` command. To switch between files, you can press `:n` for next file, or `:p` for previous file.

In the first column, we can see what proportion is classified as a certain taxonomy (last column). However, this overview is quite long and complex.
We'll use krakentools and kronatools to make a visual overview.:

```
kreport2krona.py -r Plant.kraken.report -o Plant.krona
ktImportText Plant.krona -o Plant.krona.html
kreport2krona.py -r Stool.kraken.report -o Stool.krona
ktImportText Stool.krona -o Stool.krona.html
kreport2krona.py -r Sea.kraken.report -o Sea.krona
ktImportText Sea.krona -o Sea.krona.html
```
> `kreport2krona.py` will conver the kraken output to an output that Kronatools can read
> `ktImportText` will transform the modified output to a krona plot

Now download and look at the generated .html files. You will probably see that most of the plot is unclassified reads.
This is because we are only classifying using a 16S database. 
16S is only one of the many genes in bacteria so if we sequence all the DNA from a sample, only few reads will be part of the 16S.

Try double clicking on the "Bacteria" section, in the bottom right corner

<details>
<summary>What is the most abundant group in each sample?</summary>

_Plant: Actinobacteriota-Streptomyces_

_Stool: Actinobacteriota-Bifidobacterium_

_Sea: Proteobacteria-Klebsiella_
</details>

Look a bit at the different report. You can "zoom in" on certain slices by double clicking, and return to levels above using the circles in the top right.

In a real analysis, you would not perform this analysis on a database of 16S sequences, but on a database containing full genomes. 
However, as these databases can easily take up >100 Gb it is not possible to perform these as part of the practical. 
However, in your data folder we have provided you with the results of the Kraken output on a larger database containing full genome sequences. 
Take a look at the `Plant.full.krona.html`, `Stool.full.krona.html`, and `Sea.full.krona.html` files.

<details>
<summary>What is the most common species in each sample?</summary>

- _Plant: Streptomyces rutgersensis_
- _Stool: Bifidobacterium longum_
- _Sea: Difficult to say, but probably a Candidatus Pelagibacter species_ 
</details>

<details>
<summary>Are there a lot of eukaryotic reads in the sample? If so, from which species?</summary>

- _Plant: almost none, the few that are there are from a Poaceae species_
- _Stool: almost none, the few that are there are from different plant species mainly_
- _Sea: a significant part: 11%. Mostly different plant species as well_
</details>

<details>
<summary>How diverse are the samples? Which sample looks to be the most diverse?</summary>

_The plant species is not diverse at all, the stool sample has some more species, but the sea sample is very diverse._
</details>

Kraken can give a good overview of the diversity, and since it’s read based also of the abundance of each species in the sample.
> This however assumes that every copy of a genome in the metagenome sample has the same probability of being sequenced, which is not always the case!

Another neat tool for taxonomic profiling and quantification, is [Metaphlan](https://huttenhower.sph.harvard.edu/metaphlan).
Metaphlan will map the reads to a database of marker genes (genes typical for certain clades/species).
Since the Metaphlan database is quite large (>20Gb), and the software can be quite slow, we will have provided the output in the data for this tutorial.

> If you want to run these yourself, you can use following commands:
>```
>metaphlan Plant_all.fastq --nproc 2 --input_type fastq --unclassified_estimation --index /storage/dbs/Metaphlan -o Plant_metaphlan
>metaphlan Stool_all.fastq --nproc 2 --input_type fastq --unclassified_estimation --index /storage/dbs/Metaphlan -o Stool_metaphlan
>metaphlan Sea_HiFi.fastq --nproc 2 --input_type fastq --unclassified_estimation --index /storage/dbs/Metaphlan -o Sea_metaphlan
>```

Take a look at the provided files: `Plant.metaphlan.txt`, `Stool.metaphlan.txt`, and `Sea.metaphlan.txt` files (e.g. using `cat` or `less`).
As you can see, these files are not straightforward to interpret. We can get a rough estimation of the number of identified species using `grep`:

```
grep "s__" Plant.metaphlan.txt | grep -v "t__" | wc -l
grep "s__" Stool.metaphlan.txt | grep -v "t__" | wc -l
grep "s__" Sea.metaphlan.txt | grep -v "t__" | wc -l
```
> Here we select all lines wich have "s__" in them (= identified at species level), and then filter everything out that is classified at the strain level ("t__")

<details>
<summary>How many different species are identified in each sample?</summary>

_None in the sea sample, 1 in the plant sample, 36 in the stool sample._
> The 0 in the sea sample is likely because Metaphlan can't really handle long reads very well.
</details>

<details>
<summary>Are these results in line with the Kraken results?</summary>

_More or less. We found only one main species in the plant sample and a bit more specis in the stool sample._
_We can't comment on the sea sample, since it seems to have failed._
</details>

Metaphlan can also do some extra analysis on metagenomic samples, of which we will use some in the amplicon analysis.
See the [Metaphlan](https://huttenhower.sph.harvard.edu/metaphlan) page for a full overview of all functionalities.

###  Metagenome assembly

So we made a taxonomic profile of the reads, and have a good idea of what is in there. 
However, for a more functional analysis it would be good to also have the genomes of some of these species. 
To get to these, we will first do a metagenome assembly using metaSPAdes for the short reads, and MetaFlye for the long reads.
In metagenome assembly we're just creating contigs from the reads. This works similar to a normal genome assemblies, but takes into account that multiple species are present (and thus multilpe DBGgraphs will be formed).
Metagenome assembly is a bit more resource intensive on short reads. Using the quick and efficient [MegaHit](https://github.com/voutcn/megahit) assembler, the assemblies would still take 15m each on 2 cores, and 30m with spades (1h if error-correction is performed).
To save you some time, the metaSPAdes assemblies have been run beforehand, and the resulting assembly and assembly graph have been provided (`Plant_meta.fasta`, `Stool_meta.fasta`, `Plant_meta.gfa`, `Stool_meta.gfa`)
For your information, the commands used for running the assemblies can be found below:

>```
>metaspades.py -t 2 -1 Plant_F.fastq -2 Plant_R.fastq -o Plant_metaspades
>metaspades.py -t 2 -1 Stool_F.fastq -2 Stool_R.fastq -o Stool_metaspades
>```

The reason we provided the assembly graphs (`.gfa`) is because we want to take a look at it.
A good tool to visualise the structure of our assembly graphs is [Bandage](http://rrwick.github.io/Bandage/). 
Install the software, and then download the `assembly_graph.fastg` from both assemblies (you might have to rename to avoid overwriting them when downloading). 
Start the Bandage software, and open the `.gfa` files of both assemblies. Once opened, click on the “draw graph” button to drawh the bandage graph.
You can draw a frame over a collection of lines, and see how long they are on the right side of the screen.

<details>
<summary>How many large clusters do you see? What do they correspond to according to you?</summary>

_Stool sample: There is one big connected cluster (total size +- 14 Mbp). This is likely one bacterial genome (or two closely related). Otheriwse we have smaller components (100kb and lower). These are likely contigs from different other species which are present in lower amounts in the sample._
_Plant sample: There is also one big connected cluster (total size +- 9Mbp). Again, this is likely one bacterial genomes. The other components are similar, but there are less short pieces in this sample._
</details>

Now let's assemble the long reads, using the metagenomic mode of Flye:

```
flye -t 2 --meta --pacbio-hifi Sea_HiFi.fastq -o Sea_metaFlye
```

Download the assembly graph (`assembly_graph.gfa` in the `Sea_metaFlye` folder) and have a look at it.

<details>
<summary>How is the assembly graph similar/different from the short-read assemblies?</summary>

_We have a lot fewer components, and all components are one contig and non-connected. There is no large clusters, and the contigs are quite small._
</details>

### Metagenome binning

We now have our metagenome: a collection of contigs. However, we have no idea which contigs belongs to which species/genome.
So, we need to try to separate contigs from different genomes into so-called genome bins.
Here we will use Concoct for binning, but other good binning tools are [Autometa](https://autometa.readthedocs.io/en/latest/index.html) or [MetaBat](https://bitbucket.org/berkeleylab/metabat).
We will first create a new directory for binning and put our metagenome assemblies there.

```
mkdir binning
ln -s ${PWD}/Plant_meta.fasta binning/Plant_meta.fasta
ln -s ${PWD}/Stool_meta.fasta binning/Stool_meta.fasta
ln -s ${PWD}/Sea_metaFlye/assembly.fasta binning/Sea_meta.fasta 
cd binning
```

One type of information that can help bin genomes is coverage.
We expect different genomes to have different coverages, or better stated we expected contigs coming from the same genome to have simlar covarage.
To get coverage information, we’ll have to map the reads to the metagenome. 
We’ll be using Smalt for mapping the short reads, and Minimap2 for mapping the long reads:

```
smalt index -k 11 -s 4 PMindex Plant_meta.fasta
smalt index -k 11 -s 4 SMindex Stool_meta.fasta
mkdir mapping
smalt map -n 2 PMindex ../Plant_F.fastq ../Plant_R.fastq | samtools view -bS | samtools sort > mapping/Plant.bam
smalt map -n 2 SMindex ../Stool_F.fastq ../Stool_R.fastq | samtools view -bS | samtools sort > mapping/Stool.bam
minimap2 -t 2 -ax map-hifi Sea_meta.fasta ../Sea_HiFi.fastq | samtools view -bS | samtools sort > mapping/Sea.bam
samtools index mapping/Plant.bam
samtools index mapping/Stool.bam
samtools index mapping/Sea.bam
```

To increase binning performance, Concoct wants us to cut up our contigs in similarly-sized pieces. We’ll do this by running:

```
cut_up_fasta.py Plant_meta.fasta -c 10000 -o 0 --merge_last -b Plant_10k.bed > Plant_10k.fa
cut_up_fasta.py Stool_meta.fasta -c 10000 -o 0 --merge_last -b Stool_10k.bed > Stool_10k.fa
cut_up_fasta.py Sea_meta.fasta -c 10000 -o 0 --merge_last -b Sea_10k.bed > Sea_10k.fa
```
> We will tell Concoct to cut up the contigs in pieces of 10000 bp (`-c 10000`), without allowing overlap between pieces (`-o 0`), and we will concatenate the remaining part (< 10000bp) to the final contig (`--merge_last`)

Then Concoct can create a coverage table from the mapping and the created `.bed` files from the cutting process:

```
concoct_coverage_table.py Plant_10k.bed mapping/Plant.bam > Plant_cov.tsv
concoct_coverage_table.py Stool_10k.bed mapping/Stool.bam > Stool_cov.tsv
concoct_coverage_table.py Sea_10k.bed mapping/Sea.bam > Sea_cov.tsv
```

Then we can run concoct itself:

```
concoct --threads 2 --composition_file Plant_10k.fa --coverage_file Plant_cov.tsv -b Plant_binning
concoct --threads 2 --composition_file Stool_10k.fa --coverage_file Stool_cov.tsv -b Stool_binning
concoct --threads 2 --composition_file Sea_10k.fa --coverage_file Sea_cov.tsv -b Sea_binning
```

However, now we only clustered the cut-up pieces of our assembly. So, we need to get back to the original contigs:

```
merge_cutup_clustering.py Plant_binning_clustering_gt1000.csv > Plant_binning_clustering_merged.csv
merge_cutup_clustering.py Stool_binning_clustering_gt1000.csv > Stool_binning_clustering_merged.csv
merge_cutup_clustering.py Sea_binning_clustering_gt1000.csv > Sea_binning_clustering_merged.csv
```

Then finally, we can extract the clustered bins as fasta files:

```
mkdir Plant_bins Stool_bins Sea_bins
extract_fasta_bins.py Plant_meta.fasta Plant_binning_clustering_merged.csv --output_path Plant_bins
extract_fasta_bins.py Stool_meta.fasta Stool_binning_clustering_merged.csv --output_path Stool_bins
extract_fasta_bins.py Sea_meta.fasta Sea_binning_clustering_merged.csv --output_path Sea_bins
```

### Binning QC

We now have genomic bins for each sample. A genomic bin is a set of contigs that is believed to come from the same genome. 
We can thus treat it as a genome, and use genome QC tools to assess their quality.
We can quickly have an overview of the genome statistics using Quast:

```
cd Plant_bins && quast.py *fa
cd ../Stool_bins && quast.py *fa
cd ../Sea_bins && quast.py *fa
```

Find and download the quast reports, and have a look at them.
> You can use the `pwd` command to check in what directory you currently are

<details>
<summary>Given that most bacterial genomes are between 2 and 10 Mbp, how many full genomes do you think we recovered from every sample?</summary>

- _Plant: 1 large genome (+- 7Mbp)_
- _Stool: 1 large genome (+- 12.2 Mbp), possible some smaller or uncomplete genomes (3-5 genomes between 2Mbp and 4Mbp)
- _Sea: Only 1 genome >2Mbp, so possibly no complete genomes_
</details>

<details>
<summary>Is this more or less then expected based on the taxonomic analysis (kraken/metaphlan)?</summary>

_For plant and stool yes: Plant had 1 dominant species, Stool a couple._
_There was a lot of diversity in the sea sample, but it looks like we did not get any good genomes (based on just their size)._
</details>

<details>
<summary>Why do you think we could not recover all species' genomes?</summary>

_Mainly because we are lacking coverage (we only worked on a subset of reads), so we have not enough information to recover all genomes._
_Especially species which are present in low amounts are difficult to recover from metagenomic samples._
</details>

Lastly, we will be using [CheckM](https://github.com/Ecogenomics/CheckM/wiki) to check how good the bins are that we reconstructed. 
This software will try to identify marker genes in our bins, and use that to estimate the completeness, purity (~ 1-contamination), and tries to classify our bin.

```
checkm lineage_wf Plant_bins/ Plant_checkm -x .fa -t 2 --reduced_tree
checkm lineage_wf Stool_bins/ Stool_checkm -x .fa -t 2 --reduced_tree
checkm lineage_wf Sea_bins/ Sea_checkm -x .fa -t 2 --reduced_tree
```
> We tell CheckM to analyze files ending in `.fa` (`-x .fa`), using two threads (`-t 2`). 
> We also added the `--reduced_tree` option to lower the memory requirements (14G instead of 40G for the full tree). 
> This however slightly lowers the accuracy of the analysis.

Normally, the analyses will have printed the overview table to the screen. Have a look at them.

<details>
<summary>How many good quality genomes are recovered from the data?</summary>

- _Plant: 1 high-quality genome part of the Streptomycetaceae family. The other bins show no markers._
- _Stool: Many genomes are discovered. 1 high-quality (Bifidobacteriaceae), 1 near-complete but possibly contaminated (or 2 close strains; Bacteroidales), some half-complete genomes (50-80% completeness), and the rest is likely nothing._
- _Sea: There are some genomes detected, but none are fully complete. Most of them are also difficult to classify._
</details>

<details>
<summary>How is the quality of the bins? Do you think we could have recovered more genomes if we would have used the full dataset instead ofa subsample of reads?</summary>

_We could discover some high-quality bins in Plant and Stool._
_In the plant sample, it looks like there is only one species really present._
_In the stool sample, there are more genomes present, and having the full dataset could have helped us recover more complete genomes._
_In the sea sample, we didn't manage to recover a high-quality genome, but could detect some uncomplete genomes. Using all data could help complete some of them._
</details>

As a last note: binning analysis on long reads is not always recommended. 
Most binning tools use coverage as an important factor to cluster contigs together.
Since long reads tend to have fewer reads, and thus lower coverages, it is not always easy/possible to cluster contigs together using coverage.
If enough coverage is available, a meteagenome assembly alone could be sufficient to recover the full genomes of separate species.

##  Amplicon sequencing data

In the previous part we have used whole metagenome sequencing data to recover genomes of some of the members of the bacterial communities.
We have not looked at the exact community structure and diversity. 
In this part we will look at amplicion sequencing data (only the 16S gene is sequenced, instead of all the DNA in the sample), and how we can use it to compare microbial communities.
The analysis presented here for comparing the communities an also used on shotgun sequencing data, so have a look at their documentation if you want to do that.

Because of time restraints, we have only looked at shotgun metagenome analysis. 
If you are interested in analyzing some amplicon data, here are some small exercises here regarding taxonomic profiling and comparing samples using Kraken, Bracken, and Metaphlan. 
If you want to have an overview of metagenomic tools, here is a good review: [https://link.springer.com/article/10.1007/s13238-020-00724-8](https://link.springer.com/article/10.1007/s13238-020-00724-8)

We will analyze three sets of 16S amplicon sequencing data from soil metagenomes. 
You can find them in the same data folder as the other data we used, and are named Soil_1, Soil_2, and Soil_3.

First, run `fastqc` on the three sets of reads and have a look at the reports.

<details>
<summary>Why do we observe such a high level of duplication in the data?</summary>

_Because we are dealing with amplicon data: we sequenced an amplification of markers genes. We thus only have copies of the same gene in different bacteria._
_Since we amplified (= copied) these genes multiple times, they will have a lot of reads with exactly the same sequence._
</details>

There are a tiny bit of adapters present, and the quality drops quite low near the end of the reads.
We'll do a quick cutadapt filtering to remove most of these problems:

```
cutadapt -q 30,30 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tempF.fq -p tempR.fq Soil_1F.fastq Soil_1R.fastq
mv tempF.fq Soil_1F.fastq
mv tempR.fq Soil_1R.fastq

cutadapt -q 30,30 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tempF.fq -p tempR.fq Soil_2F.fastq Soil_2R.fastq
mv tempF.fq Soil_2F.fastq
mv tempR.fq Soil_2R.fastq

cutadapt -q 30,30 -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o tempF.fq -p tempR.fq Soil_3F.fastq Soil_3R.fastq
mv tempF.fq Soil_3F.fastq
mv tempR.fq Soil_3R.fastq
```

Since we sequenced PCR-fragments (amplicons), in many cases the insert size (size of the amplicon) will be shorter than 2x read length.
This means there will often be overlap between the read pairs. We'll be using [FLASH](https://github.com/dstreett/FLASH2) to merge the reads:
> It is debated wether it is best to merge first or do trimming first. Here we pick trimming first since the quality dropped quite low sometimes near the end of the reads.
> Doing merging first might be advisable if you have low amount of reads and want to keep as many as possible.

```
flash2 Soil_1F.fastq Soil_1R.fastq -M 200
mv out.extendedFrags.fastq Soil_1M.fastq
mv out.notCombined_1.fastq Soil_1F.fastq
mv out.notCombined_2.fastq Soil_1R.fastq

flash2 Soil_2F.fastq Soil_2R.fastq -M 200
mv out.extendedFrags.fastq Soil_2M.fastq
mv out.notCombined_1.fastq Soil_2F.fastq
mv out.notCombined_2.fastq Soil_2R.fastq

flash2 Soil_3F.fastq Soil_3R.fastq -M 200
mv out.extendedFrags.fastq Soil_3M.fastq
mv out.notCombined_1.fastq Soil_3F.fastq
mv out.notCombined_2.fastq Soil_3R.fastq
```
> The `-M` parameter increases the maximal allowed overlap, since we're dealing with relatively long reads (250bp)

Mergeing reads can also be performed to improve genome assembly, as longer reads will make assembly easier.
As an additional bonus, it will lower the number of total reads, leading to shorter running times.

Let's analyse our soil samples with Kraken first. 
We have to again put all reads in one file:

```
cat Soil_1*.fastq > Soil_1A.fastq
cat Soil_2*.fastq > Soil_2A.fastq
cat Soil_3*.fastq > Soil_3A.fastq
```

As we are dealing with 16S data, using the SILVA database should give better results then then we were working with shotgun data.

```
kraken2 --db /storage/dbs/KrakenSILVADB/ --threads 2 --output Soil_1.kraken.out --report Soil_1.kraken.report Soil_1A.fastq
kraken2 --db /storage/dbs/KrakenSILVADB/ --threads 2 --output Soil_2.kraken.out --report Soil_2.kraken.report Soil_2A.fastq
kraken2 --db /storage/dbs/KrakenSILVADB/ --threads 2 --output Soil_3.kraken.out --report Soil_3.kraken.report Soil_3A.fastq
```

We’ll again visualize the results in a krona plot:
> An alternative way to visualise the kraken reports, is on the webserver [Pavian](https://fbreitwieser.shinyapps.io/pavian/)

```
kreport2krona.py -r Soil_1.kraken.report -o Soil_1.krona && ktImportText Soil_1.krona -o Soil_1.krona.html
kreport2krona.py -r Soil_2.kraken.report -o Soil_2.krona && ktImportText Soil_2.krona -o Soil_2.krona.html
kreport2krona.py -r Soil_3.kraken.report -o Soil_3.krona && ktImportText Soil_3.krona -o Soil_3.krona.html
```

Now open the krona plots and have a look at the results.

<details>
<summary>Are the samples diverse, or are there only few dominating species?</summary>

_All 3 samples look very divers, with not really a single species dominating the community._
</details>

<details>
<summary>Do the 3 samples have similar species distribution on first glance?</summary>

_Soil1 and Soil2 have very similar distributions, Soil3 has the same large groups, but the proportions are a bit different._
</details>

Now we’ll use a tool called [Bracken](https://ccb.jhu.edu/software/bracken/) to calculate species abundance for all three samples based on the Kraken results.
Bracken re-estimates the read counts from kraken, to make the more robust and comparable between samples. 
Since the SILVA database only classifies on genus level, we will do the analysis on genes level. 
Furthermore, we will limit ourselves to species with at least 10 reads:

```
bracken -d /storage/dbs/KrakenSILVADB/ -i Soil_1.kraken.report -o Soil_1.bracken -r 300 -l "G" -t 10
bracken -d /storage/dbs/KrakenSILVADB/ -i Soil_2.kraken.report -o Soil_2.bracken -r 300 -l "G" -t 10
bracken -d /storage/dbs/KrakenSILVADB/ -i Soil_3.kraken.report -o Soil_3.bracken -r 300 -l "G" -t 10
```
> We specify the average read length (`-l 300`), and specify we want to classify on Genus level (`-l "G"`). `-t 10` filters out results with fewer than 10 mapped reads.
> We use a read length of 300 since most of our merged reads are that length (= size of the amplicon).

Have a look at the Bracken files (using `less`, `head`, or `cat`).
You will see for every discovered genus the assigned reads, added reads (from higher/lower-level taxonomies), and the new amounts of reads.

Let's try finding the most abundant genera:

```
sort -nrk6 Soil_1.bracken | head -n10
sort -nrk6 Soil_2.bracken | head -n10
sort -nrk6 Soil_3.bracken | head -n10
```
> We sort numerically (`-n`), in reverse order (`-r`), and use field 6 (corrected reads) (`k6`) to sort.
> Then, we take the first 10 lines using head.

<details>
<summary>What is the most abundant genus in each sample?</summary>

- _Soil1: Acidothermus_
- _Soil2: uncultured (followed by Acidothermus)_
- _Soil3: Gaiella_
</details>

Using these estimations, we can try to estimate the alpha diversity (using the [KrakenTools](https://github.com/jenniferlu717/KrakenTools) package). 
Many diversity measurements are implemented in the software, and here we will be using the Shannon alpha diversity:

```
alpha_diversity.py -f Soil_1.bracken -a Sh
alpha_diversity.py -f Soil_2.bracken -a Sh
alpha_diversity.py -f Soil_3.bracken -a Sh
```
> `-a Sh` specifies that we want to use the Shanon Index

<details>
<summary>Which of the samples has the highest diversity?</summary>

_Soil3 (4.769)_
</details>

We can also calculate the beta-diversity between the three samples. 
The `beta_diversity.py` script from KronaTools calculates the Bray-Curtis dissimilarity matrix, 
where a score of 0 means samples with exactly the same species distribution, and 1 are totally different samples.

```
beta_diversity.py -i *bracken --type bracken
```

<details>
<summary>Do the observed differences between samples are what you expected based on the Krona plots?</summary>

_Yes, Soil1 and Soil2 have a low dissimilarity (0.382), while Soil3 is more dissimilar to the other 2 (0.426-0.562)._
</details>

Now, we will also use Metaphlan to analyse the samples. 
Sadly, 16S sequences are not included in the Metaphlan marker database, so running metaphlan itself is no use here. 
However, using KronaTools, we can convert our Kraken2 output to a metaphlan output. 
As not all libraries contain the same number of reads, we will use percentages instead of read counts (which is also used by metaphlan itself):

```
kreport2mpa.py -r Soil_1.kraken.report -o Soil_1.metaphlan.txt --percentages
kreport2mpa.py -r Soil_2.kraken.report -o Soil_2.metaphlan.txt --percentages
kreport2mpa.py -r Soil_3.kraken.report -o Soil_3.metaphlan.txt --percentages
```

Our kraken analysis only goes to genus level, so we’ll try counting the number of genera in the samples:

```
grep "g__" Soil_1.metaphlan.txt | wc -l
grep "g__" Soil_2.metaphlan.txt | wc -l
grep "g__" Soil_3.metaphlan.txt | wc -l
```

It looks like there are many genera identified, many of the with only few reads. 
Let’s do some filtering. We will remove “uncultured” genera, and only keep genera with a higher proportion than 0.5%.
We'll do this using some build-in Linux functions:

```
grep "g__" Soil_1.metaphlan.txt | awk '{ if ($NF > 0.49) print}' | grep -v "uncultured" > Soil_1.metaphlan.reduced.txt
grep "g__" Soil_2.metaphlan.txt | awk '{ if ($NF > 0.49) print}' | grep -v "uncultured" > Soil_2.metaphlan.reduced.txt
grep "g__" Soil_3.metaphlan.txt | awk '{ if ($NF > 0.49) print}' | grep -v "uncultured" > Soil_3.metaphlan.reduced.txt
```
> In this command, we filter the lines that have identifications on genus level `grep "g__"`, 
> and then only retain lines where the last column has a value strictly higher than 0.49 (`awk '{ if ($NF > 0.49) print}'`. 
> Lastly, we use grep again to filter out lines which contain "uncultured" somewhere (`grep -v "uncultured"`)

Since we only filtered out the "genus" lines, we can simply count the number of discovered genera by counting the lines in the files:

```
wc -l Soil*reduced.txt
```

This reduces the number of genera significantly. We can merge the Metaphlan profiles, to create one big table:

```
combine_mpa.py -i *reduced.txt -o Soil_merged.tbl
```

<details>
<summary>The combine_mpa script tells us that it printed 26 classifications, but if we count the genera in all 3 samples, we have 47. Why is there a difference?</summary>

_Because there is overlap in genera between samples. There are thus 26 unique genera found, but some of them are present in multiple samples._
</details>

If you look at the merged table, the names of the genera are quite long, containing the whole taxonomy.
We’ll clean up the table a bit using `sed`, by removing everything before the last `|` character:

```
sed -i "s/^.*|//g" Soil_merged.tbl
```
> Here, we are using the substitute command of  `sed`: `sed s/pattern/replacement/g`
> The leading `s` tells `sed` to do a substitution, and the leading `g` tells `sed` to do it "globally" (so don't stop after the first replacement.
> We tell `sed` to substitute the pattern `^.*|` and replace it with nothing (there is nothing between the last forwards slashes: `//`) 
> The `^.*|` pattern roughly translates to "any number of any character (`.*`) after the start of the string (`^`) until you meet a `|`.  
> This will effectively remove everything before (and including) the last `|` character.


Lastly, we need to slightly change the header before we can visualize the result.
Open the merged table using nano (`nano Soil_merged.tbl`), and and change the header by removing the leading `#` in the first line, and change the sample headers (e.g. Sample #1 in Soil1).
Press `Ctrl-O` to save the file, and `Ctrl-X` to quit.
Now let's have a look at the table:

```
column -t Soil_merged.tbl
```
> The `column` command will nicely align your columns

<details>
<summary>Compare the most abundant genera from the heatmap to the those of the first Krona plot. Do they match?</summary>

_Soil1: Acidothermus; Soil2: Acidothermus; Soil3: Nocardioides_
_This matches the output from the Krona plots._
</details>

As you might notice, the results from both analyses differ from each other. 
Without going into too much detail, the main difference is that our first analysis we used the Bracken results, and in our Metaphlan-based analysis, we used the Kraken result. 
As seen in the lecture, Kraken classifies the reads using the LCA approach (Last Common Ancestor). 
As such, not all reads are classified will be classified to species/genus level. 
As such, our method of extracting just the rows that are classified to the genus level will only represent a subset of our total reads. 
That’s why it’s better to run Bracken after you run Kraken, as Bracken modifies the classification to make sure that reads are classified on the desired taxonomical level. 
If you want, you could try to run the Metaphlan based analysis on the bracken output rather than the kraken output and compare the results.
