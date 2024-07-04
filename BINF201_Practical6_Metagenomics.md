# BINF201 – Practical 6 – Metagenomics

For this practical, you will need to have installed and configured Conda. If this is not the case, see [here](https://docs.conda.io/en/latest/miniconda.html) for installation instructions for your system.

> For Windows users, it is preferable to use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install), and install conda on there using the Linux install instructions.

In this practical we will work with metagenomic data. We will analyze both shotgun sequencing data and amplicon sequencing data.

## Software setup & data retrieval

As usual, we will set up a Conda environment with the necessary tools. However, as we have many different tools to install, installation will run smoother using Mamba instead of Conda. If you have not already installed [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html), you can do this now by running:

	conda install mamba

Mamba is a Conda “replacement” that can resolve dependencies a lot faster than Conda. We will use it to install our metagenomic environment:

	mamba create -n metagenomics -c bioconda -c conda-forge spades kraken2 krona metaphlan fastqc fastp checkm-genome krakentools hclust2 concoct smalt quast bracken

To activate the environment, you still use Conda:

	conda activate metagenomics

You will probably get a message from Krona saying you need to manually update the taxonomy data. We’ll do this first by running:

	ktUpdateTaxonomy.sh

> This 	step might give an error. If this is the case, don't worry. The software should still rune fine. We only may not have access to the latest taxonomy information.

You will also get another message saying you set the CheckM root. Do this as well by running the advised command (`checkm data setRoot path/that/is/displayed`)

We will also have to download a database for classification using kraken. We will be using the [SILVA](https://www.arb-silva.de/) 16S bacterial and archaea database, as it is relatively small:

	kraken2-build -db SILVA --special silva

All the necessary data for this tutorial can be found on MittUiB: exercises -> data -> Practical6_Metagenomics.zip

Once you have downloaded and unzipped this archive, think about decompressing all `.gz` files (using gunzip, for example).

##  Quality control & trimming

In this part we will analyze two shotgun metagenomes. One metagenome (SRA-accession SRR23255609) is from a human stool sample from a patient who had a liver transplant. The other metagenome (ERR2969992) is shotgun sequencing from plant tissue containing a leaf symbiont. For convenience, these datasets have been reduced to 2x1M reads.

First, we will perform the QC and trimming with [fastp](https://github.com/OpenGene/fastp).

	fastp -i Stool_1.fastq -I Stool_2.fastq -o Stool_1.trimmed.fastq -O Stool_2.trimmed.fastq -h Stool.report.html
	fastp -i Plant_1.fastq -I Plant_2.fastq -o Plant_1.trimmed.fastq -O Plant_2.trimmed.fastq -h Plant.report.html

To have a better look at the data, we will also run FastQC on the trimmed data:

	fastqc *trimmed.fastq

Then look at the FastQC reports first, and look at the “Per sequence GC-content”

1.  Why does FastQC give an error?

2.  Is the observed distribution what we expect based on the type of data we are analyzing?

Now look at the fastp reports and look at the insert size distribution and the read length.

3.  What is the average insert size? (roughly)

4.  Given our read length, and the insert sizes we observe, do we expect to have overlap between the reads?

As there is potentially some overlap between reads, we will rerun fastp, but this time, we will add the `--merge` flag This option will try to merge paired reads that overlap (i.e. with small insert size) into one read, which we can use later to improve our assembly:

	fastp -i Stool_1.fastq -I Stool_2.fastq -o Stool_1.trimmed.fastq -O Stool_2.trimmed.fastq --merge --merged_out Stool_m.trimmed.fastq -h Stool.report.html
	fastp -i Plant_1.fastq -I Plant_2.fastq -o Plant_1.trimmed.fastq -O Plant_2.trimmed.fastq --merge --merged_out Plant_m.trimmed.fastq -h Plant.report.html

## Taxonomic profiling of reads

Now that we have QC’ed and filtered our reads, we could try a first classification on read level. We will first use [kraken2](https://ccb.jhu.edu/software/kraken2/) for this. Kraken extracts k-mers from our reads and compares them to a database to classify them. We will first concatenate all our reads from one sample together:

	cat Plant_1.trimmed.fastq Plant_2.trimmed.fastq Plant_m.trimmed.fastq > Plant_all.trimmed.fastq
	cat Stool_1.trimmed.fastq Stool_2.trimmed.fastq Stool_m.trimmed.fastq > Stool_all.trimmed.fastq

Then we will run kraken using the SILVA database downloaded before:

	kraken2 --db SILVA/ --threads 2 --output Plant.kraken.out --report Plant.kraken.report Plant_all.trimmed.fastq
	kraken2 --db SILVA/ --threads 2 --output Stool.kraken.out --report Stool.kraken.report Stool_all.trimmed.fasta

We will visualise the classification using kronatools:

	kreport2krona.py -r Plant.kraken.report -o Plant.krona
	ktImportText Plant.krona -o Plant.krona.html
	kreport2krona.py -r Stool.kraken.report -o Stool.krona
	ktImportText Stool.krona -o Stool.krona.html

> `kreport2krona.py` will conver the kraken output to an output that Kronatools can read
> `ktImportText` will transform the modified output to a krona plot

Now look at the generated .html files.

1.  Are there many reads classified? Why (not)?

Zoom in on the Bacteria clade (click on the Bacteria name on the bottom right)

2.  What is the most abundant order (`o__`), family (`f__`), and genus (`g__`) in the data? 
(hint: you can see the exact percentage and number of reads by clicking on/hovering over  a wedge)

3.  In the Plant sample, we have a large chunk of `o__Chloroplast` and also some `f__Mitochondria`. Why are these classified as Bacteria?

In a real analysis, you would not perform this analysis on a database of 16S sequences, but on a database containing full genomes. However, as these databases can easily take up >100 Gb it is not possible to perform these as part of the practical. However, in your data folder we have provided you with the results of the Kraken output on a larger database containing full genome sequences. Take a look at the `Plant_allReads.krona.html` and `Stool_allReads.krona.html`

4.  Can you identify which two species are the main species of the symbiosis in the Plant sample?

5.  Look up information about the bacterial species. Is the host plant correctly identified?

6.  Do you think it is necessary for further analysis to remove “contaminating” human DNA from the Stool sample?

7.  Are there big differences in bacterial diversity of the stool sample between this analysis and your analysis using a smaller database?

Kraken can give a good overview of the diversity, and since it’s read based also of the abundance of each species in the sample. However, for a proper analysis you need to use a large database, which is not always feasible. [Metaphlan](https://huttenhower.sph.harvard.edu/metaphlan) is a good alternative (or supplementary) tool. Metaphlan will map the reads to a database of marker genes (genes typical for certain clades/species). Sadly, the Metaphlan database is still quite large (>20 Gb), so for your convenience we have provided you with the output.

Open the files `Plant.metaphlan.txt` and `Stool.metaphlan.txt`

As you can see, these files are not straightforward to interpret. We’ll try get a rough estimation of the number of identified species using grep:

	grep "s__" Plant.metaphlan.txt | grep -v "t__" | wc -l
	grep "s__" Stool.metaphlan.txt | grep -v "t__" | wc -l

> Here we select all lines wich have "s__" in them (= identified at species level), and then filter everything out that is classified at the strain level ("t__")

8.  How many species are identified in each sample?

Metaphlan can do some extra analysis on metagenomic samples. If you’re interested, see the optional part of this practical, or see the [Metaphlan](https://huttenhower.sph.harvard.edu/metaphlan) page.

##  Metagenomic assembly & binning

So we made a taxonomic profile of the reads, and have a good idea of what is in there. However, for a more functional analysis it would be good to also have the genomes of some of these species. To get to these, we will first do a metagenome assembly using [SPAdes](https://github.com/ablab/spades), and then separate the contigs from different species (= binning). First, let’s do the metagenome assembly. Here we will also use our merged reads for increased assembly accuracy:

	metaspades.py -1 Plant_1.trimmed.fastq -2 Plant_2.trimmed.fastq --merged Plant_m.trimmed.fastq -o Plant_metaspades
	metaspades.py -1 Stool_1.trimmed.fastq -2 Stool_2.trimmed.fastq --merged Stool_m.trimmed.fastq -o Stool_metaspades

Let’s check our metagenome assembly using [Bandage](http://rrwick.github.io/Bandage/). Install the software if you haven’t already and launch it. Then, navigate to your working directory and look in the spades output folders. Start by opening the `assembly_graph.fastg` file of the Plant sample. Once opened, click on the “draw graph” button.

1.  How many contigs do we have in total?

2.  How many large clusters do you see? Is this what you expect?

Now open the one of the Stool assemblies

3.  Do you have more or fewer contigs?

4.  Do you have more or fewer clusters?

5.  Based on the graph, how many good genomes do you think we can recover?

A metagenome assembly is just an intermediate step in this analysis. With the binning, we hope to separate the contigs that belong to different species. We will bin the contigs using a software called [CONCOCT](https://github.com/BinPro/CONCOCT). We will first create a new directory for binning and put our metagenome assemblies there.

	mkdir binning
	cp Plant_metaspades/contigs.fasta binning/Plant_meta.fasta
	cp Stool_metaspades/contigs.fasta binning/Stool_meta.fasta
	cd binning

To get coverage information, we’ll have to map the reads to the metagenome. We’ll be using Smalt for this. We use the unfiltered reads, to make sure we have the best coverage information:

	smalt index -k 11 -s 4 PMindex Plant_meta.fasta
	smalt index -k 11 -s 4 SMindex Stool_meta.fasta
	mkdir mapping
	smalt map -n 2 PMindex ../Plant_1.fastq ../Plant_2.fastq | samtools view -bS | samtools sort > mapping/Plant.bam
	smalt map -n 2 SMindex ../Stool_1.fastq ../Stool_2.fastq | samtools view -bS | samtools sort > mapping/Stool.bam
	samtools index mapping/Plant.bam
	samtools index mapping/Stool.bam

To increase binning performance, Concoct wants us to cut up our contigs in similarly-sized pieces. We’ll do this by running:

	cut_up_fasta.py Plant_meta.fasta -c 10000 -o 0 --merge_last -b Plant_10k.bed > Plant_10k.fa
	cut_up_fasta.py Stool_meta.fasta -c 10000 -o 0 --merge_last -b Stool_10k.bed > Stool_10k.fa

> We will tell Concoct to cut up the contigs in pieces of 10000 bp (`-c 10000`), without allowing overlap between pieces (`-o 0`), and we will concatenate the remaining part (< 10000bp) to the final contig (`--merge_last`)

Then Concoct can create a coverage table from the mapping and the created `.bed` files:

	concoct_coverage_table.py Plant_10k.bed mapping/Plant.bam > Plant_cov.tsv
	concoct_coverage_table.py Stool_10k.bed mapping/Stool.bam > Stool_cov.tsv

Then we can run concoct itself:

	concoct --composition_file Plant_10k.fa --coverage_file Plant_cov.tsv -b Plant_binning
	concoct --composition_file Stool_10k.fa --coverage_file Stool_cov.tsv -b Stool_binning

> You might run into an error here. If not, good for you! If so, you’ll
> have to change something in a validation script of scikit-learn, a
> python package. Open the script in nano:
> 
> 		nano +1880 ~/miniconda3/envs/metagenomics/lib/python3.10/site-packages/sklearn/utils/validation.py
> 
> Then, search for the line `feature_names = np.asarray(X.columns,
> dtype=object)` Change this line to `feature_names =
> np.asarray(X.columns.astype(str), dtype=object)` Once, done save the
> file by tapping `Ctrl+O`, and quit by tapping `Ctrl+X.`
> This should fix the error.

However, now we only clustered the cut-up pieces of our assembly. So, we need to get back to the original contigs:

	merge_cutup_clustering.py Plant_binning_clustering_gt1000.csv > Plant_binning_clustering_merged.csv
	merge_cutup_clustering.py Stool_binning_clustering_gt1000.csv > Stool_binning_clustering_merged.csv

Then finally, we can extract the clustered bins as fasta files:

	mkdir Plant_bins Stool_bins
	extract_fasta_bins.py Plant_meta.fasta Plant_binning_clustering_merged.csv --output_path Plant_bins
	extract_fasta_bins.py Stool_meta.fasta Stool_binning_clustering_merged.csv --output_path Stool_bins

We can quickly have an overview of the created bins using quast:

	cd Plant_bins && quast.py *fa
	cd ../Stool_bins && quast.py *fa

Now have a look at the Quast reports for both samples.

6.  Given that most bacterial genomes are between 2 and 10 Mbp, how many full genomes do you think we recovered from every sample?

7.  Is this in agreement with what we found in the kraken analysis and the metagenome graph analysis?

Lastly, we will be using [CheckM](https://github.com/Ecogenomics/CheckM/wiki) to check how good the bins are that we reconstructed. This software will try to identify marker genes in our bins, and use that to estimate the completeness, purity (~1-contamination), and tries to classify our bin.

	checkm lineage_wf Plant_bins/ Plant_checkm -x .fa -t 2 --reduced-tree
	checkm lineage_wf Stool_bins/ Stool_checkm -x .fa -t 2 --reduced-tree

> We tell Checkm to analyze files ending in `.fa` (`-x .fa`), using two threads (`-t 2`). We also added the `--reduced-tree` option to lower the memory requirements (14G instead of 40G for the full tree). This however slightly lowers the accuracy of the analysis.

>This analysis requires 14Gb of memory, which you might not have on your laptop. If the analysis fails, you can use the pre-generated output found in the data-folder (`Plant.checkm.out` & `Stool.checkm.out`).

Normally, the analyses will have printed the overview table to the screen.

8.  How many good quality genomes are recovered from the data? Is this what you have expected?

> Here we worked with only a subset of reads (2x1M). If we would have worked with the full dataset, we could have maybe recovered more high-quality genomes.

## Optional:  Amplicon sequencing data

Because of time restraints, we have only looked at shotgun metagenome analysis. If you are interested in analyzing some amplicon data, here are some small exercises here regarding taxonomic profiling and comparing samples using Kraken, Bracken, and Metaphlan. If you want to have an overview of metagenomic tools, here is a good review: [https://link.springer.com/article/10.1007/s13238-020-00724-8](https://link.springer.com/article/10.1007/s13238-020-00724-8)

We will analyze three sets of 16S amplicon sequencing data from soil metagenomes (SRR25440994, SRR25441013, and SRR25440997). You can find them in the same data folder as the other data we used, and are named Soil1, Soil2, and Soil3.

First, we will start by running FastQC on the reads:

	fastqc Soil*fastq

Now look at the FastQC reports.

1.  Why do we observe very high duplication levels and overrepresented sequences?

As there seem to be some adapter sequences present, and the quality drops off a bit near the end of the reads, we’ll run fastp on the reads. Since we are dealing with 16S amplicons, which are not that long (+- 300 bp), we will again merge the reads to create full-length amplicons:

	fastp -i Soil1_1.fastq -I Soil1_2.fastq -o Soil1_1.trimmed.fastq -O Soil1_2.trimmed.fastq -h Soil1.fastp.html --merge --merged_out Soil1_m.trimmed.fastq
	fastp -i Soil2_1.fastq -I Soil2_2.fastq -o Soil2_1.trimmed.fastq -O Soil2_2.trimmed.fastq -h Soil2.fastp.html --merge --merged_out Soil2_m.trimmed.fastq
	fastp -i Soil3_1.fastq -I Soil3_2.fastq -o Soil3_1.trimmed.fastq -O Soil3_2.trimmed.fastq -h Soil3.fastp.html --merge --merged_out Soil3_m.trimmed.fastq

Then, we’ll again put all reads in one file for the kraken analysis:

	cat Soil1*trimmed.fastq > Soil1_all.trimmed.fastq
	cat Soil2*trimmed.fastq > Soil2_all.trimmed.fastq
	cat Soil3*trimmed.fastq > Soil3_all.trimmed.fastq

As we are dealing with 16S data, using the SILVA database should give better results then then we were working with shotgun data.

	kraken2 --db SILVA/ --threads 2 --output Soil1.kraken.out --report Soil1.kraken.report Soil1_all.trimmed.fastq
	kraken2 --db SILVA/ --threads 2 --output Soil2.kraken.out --report Soil2.kraken.report Soil2_all.trimmed.fastq
	kraken2 --db SILVA/ --threads 2 --output Soil3.kraken.out --report Soil3.kraken.report Soil3_all.trimmed.fastq

We’ll again visualize the results in a krona plot:

	kreport2krona.py -r Soil1.kraken.report -o Soil1.krona && ktImportText Soil1.krona -o Soil1.krona.html
	kreport2krona.py -r Soil2.kraken.report -o Soil2.krona && ktImportText Soil2.krona -o Soil2.krona.html
	kreport2krona.py -r Soil3.kraken.report -o Soil3.krona && ktImportText Soil3.krona -o Soil3.krona.html

Now open the krona plots and have a look at the results.

2.  Are the samples very diverse (i.e., are there many or only few different species?)

3. Would you say the three samples are similar or very different from each other (looking on class/phylum level)?

Now we’ll use a tool called [Bracken](https://ccb.jhu.edu/software/bracken/) to calculate species abundance for all three samples. To do this, we first need to prepare our database. We can do this by running the following command:

	bracken-build -d SILVA/ -t 2 -k 35 -l 250

> We do this using two threads (`-t 2`). We also have to specify the *k*-mer value used for searching (`-k 35`, as 35 is the standard *k*-mer length used in Kraken2), and the estimated read length (`-l 250`). We could also try using a read length of 300 here, as a lot of our reads are merged into ~300bp amplicons

Then we can run bracken. Since the SILVA database only classifies on genus level, we will do the analysis on genes level. Furthermore, we will limit ourselves to species with at least 10 reads:

	bracken -d SILVA/ -i Soil1.kraken.report -o Soil1.bracken -r 250 -l "G" -t 10
	bracken -d SILVA/ -i Soil2.kraken.report -o Soil2.bracken -r 250 -l "G" -t 10
	bracken -d SILVA/ -i Soil3.kraken.report -o Soil3.bracken -r 250 -l "G" -t 10

> Again, we specify the average read length (`-l 250`), and specify we want to classify on Genus level (`-l "G"`). `-t 10` filters out results with fewer than 10 mapped reads

Using these estimations, we can try to estimate the alpha diversity (using the [KrakenTools](https://github.com/jenniferlu717/KrakenTools) package). Many diversity measurements are implemented in the software, and here we will be using the Shannon alpha diversity:

	alpha_diversity.py -f Soil1.bracken -a Sh
	alpha_diversity.py -f Soil2.bracken -a Sh
	alpha_diversity.py -f Soil3.bracken -a Sh

> `-a Sh` specifies that we want to use the Shanon Index

4.  Which of the samples has the highest alpha diversity?

5.  Do the differences in alpha-diversity reflect what you saw the krona plots?

We can also calculate the beta-diversity between the three samples. The `beta_diversity.py` script from KronaTools calculates the Bray-Curtis dissimilarity matrix, where a score of 0 means samples with exactly the same species distribution, and 1 are totally different samples.

	beta_diversity.py -f *bracken --type bracken

Lastly, we can browse our bracken results online, using Pavian: [https://fbreitwieser.shinyapps.io/pavian/](https://fbreitwieser.shinyapps.io/pavian/)

Go to the website, and upload the three `.report` files generated by Kraken. Then, go to the “Sample” view and “Comparison” view on the left of the screen.

6.  Which genus is the most abundant in each sample (excluding "uncultured")?

Now, we will also use Metaphlan to analyse the samples. Sadly, 16S sequences are not included in the Metaphlan marker database, so running metaphlan itself is no use here. However, using KronaTools, we can convert our Kraken2 output to a metaphlan output. As not all libraries contain the same number of reads, we will use percentages instead of read counts (which is also used by metaphlan itself):

	kreport2mpa.py -r Soil1.kraken.report -o Soil1.metaphlan.txt --percentages
	kreport2mpa.py -r Soil2.kraken.report -o Soil2.metaphlan.txt --percentages
	kreport2mpa.py -r Soil3.kraken.report -o Soil3.metaphlan.txt --percentages

Again, our kraken analysis only goes to genus level, so we’ll try counting the number of genera in the samples:

	grep "g__" Soil1.metaphlan.txt | wc -l
	grep "g__" Soil2.metaphlan.txt | wc -l
	grep "g__" Soil3.metaphlan.txt | wc -l

7.  How many genera are identified in each sample?

It looks like there are many genera identified, many of the with only few reads. Let’s do some filtering. We will remove “uncultured” genera, and only keep genera with a higher proportion than 0.5%

	grep "g__" Soil1.metaphlan.txt | awk '{ if ($NF > 0.49) print}' | grep -v "uncultured" > Soil1.metaphlan.reduced.txt
	grep "g__" Soil2.metaphlan.txt | awk '{ if ($NF > 0.49) print}' | grep -v "uncultured" > Soil2.metaphlan.reduced.txt
	grep "g__" Soil3.metaphlan.txt | awk '{ if ($NF > 0.49) print}' | grep -v "uncultured" > Soil3.metaphlan.reduced.txt
 
> In this command, we filter the lines that have identifications on genus level `grep "g__"`, and then only retain lines where the last column has a value strictly higher than 0.49 (`awk '{ if ($NF > 0.49) print}'`. Lastly, we use grep again to filter out lines which contain "uncultured" somewhere (`grep -v "uncultured"`
Then you can count the number of identified genera using `wc -l FILE`

8.  How many genera are identified in each sample after filtering?

Now, we will merge the Metaphlan profiles, to create one big table:

	combine_mpa.py -i *reduced.txt -o Soil_merged.tbl

We’ll clean up the table a bit using `sed`, by removing everything before the last “|” character:

	sed -i "s/^.*|//g" Soil_merged.tbl

> Here, we are using the substitute command of  `sed`: `sed s/pattern/replacement/g`
> The leading `s` tells `sed` to do a substitution, and the leading `g` tells `sed` to do it "globally" (so don't stop after the first replacement.
> We tell `sed` to substitute the pattern `^.*|` and replace it with nothing (there is nothing between the last forwards slashes: `//`) 
> The `^.*|` pattern roughly translates to "any number of any character (`.*`) after the start of the string (`^`) until you meet a `|`.  This will effectively remove everything before (and including) the last `|` character.  
> 
Lastly, we need to slightly change the header before we can visualize the result. 
Go into a text editor and open `Soil_merged.tbl` and change the header by removing the leading “#”, and changing “Sample #X” by "SoilX".

Then we’ll create a heatmap of the results:

	hclust2.py -i Soil_merged.tbl -o Soil_metaphlan.png --f_dist_f correlation \
	--s_dist_f braycurtis --cell_aspect_ratio 0.5 --flabel_size 10 \
	--slabel_size 10 --max_flabel_len 100 --max_slabel_len 100 --dpi 300

> For more info on the hclust2 options, see https://github.com/SegataLab/hclust2

This will create a heatmap with all 3 samples and the identified genera (>0.5% proportion). Open the heatmap and have a look at it.

9.  What are the most abundant genera in the samples?

As you might notice, the results from both analyses differ from each other. Without going into too much detail, the main difference is that our first analysis we used the Bracken results, and in our Metaphlan-based analysis, we used the Kraken result. As seen in the lecture, Kraken classifies the reads using the LCA approach (Last Common Ancestor). As such, not all reads are classified will be classified to species/genus level. As such, our method of extracting just the rows that are classified to the genus level will only represent a subset of our total reads. That’s why it’s better to run Bracken after you run Kraken, as Bracken modifies the classification to make sure that reads are classified on the desired taxonomical level. If you want, you could try to run the Metaphlan based analysis on the bracken output rather than the kraken output and compare the results.
