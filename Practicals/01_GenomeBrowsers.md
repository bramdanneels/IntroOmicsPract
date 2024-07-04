#Practical 1 - Genome browsers
In this first practical exercise we will have a look at different genome browsers, and familiarize with them. The main aim is to make you comfortable using genome browsers, as we will use genome browsers throughout the practical exercises to visualize some of our results. 

## USCS Genome Browser

First, we'll have a look at one of the most detailed online genome browsers: the UCSC genome browser

Go to the [UCSC genome browser webpage](https://genome-euro.ucsc.edu/).

Go to the Genome Browser tool, and pick the hg38 version of the human genome (GRCh38/hg38).

Find the BRCA1 gene using the search bar.

<details>
<summary>The gene models (in blue) have 3 types of lines: thick lines, medium lines, and fine lines (sometimes with arrows). What do these 3 different lines represent?</summary>
-Thin lines: introns
-Thick lines: translated exons
-Medium lines: un-translated exons (e.g. UTRs)
</details>

Open the info panel for the BRCA1 gene by clicking on one of the transcripts. Answer the following questions regarding BRCA1:

<details>
<summary>What is the function of the BRCA1 gene in humans?</summary>
Maintaining genome stabolity and tumor suppression.
</details>

<details>
<summary>Where on the genome is BRCA1 located?</summary>
Chromosome 17, bp 43044295 - 43125402 (chr17:43,044,295-43,125,402)
</details>

<details>
<summary>Is BRCA1 located on the forward or reverse strand?</summary>
The reverse strand (Strand: -)
</details>

Go back to the genome browser view, and try finding the start codon (= translation start) of a transcript.
> The easiest way to zoom in on a specific region, is to drag the region of interest to the middle of the screen, and then click on one of the zoom buttons on the top of the search

<details>
<summary>In what colour do the start codons have on the amino acid sequence?</summary>
Green
</details>

Zoom in further, untill you can see the nucleotide sequence on the top of the screen.

<details>
<summary>What sequence does the start codon have? What amino acid does it code for?</summary>
ATG, coding for Methionine (M)
> If you guessed CAT or TAC, remember that the gene is on the reverse strand. This means that you have to read the sequence from right to left, and have to take the complement (= reverse complement) compared to the reference sequence to have the actual sequence
</details>

The genome browser has actually a lot more information than you can see at first sight. The data is added in so-called "tracks". Each tracks represents a type of information that is overlayed on the genomic sequence.
>You can collapse or expand tracks by clicking on the track title (the text in the middle of the screen above a certain track)

We will try looking at some extra information from other genome tracks. All tracks that can be added or removed to the display are located below the display. They are categorized by type, and can be selected from the drop-down menu or by clicking on the track type and configuring the track.

First, go to the ARTN gene.

We want to have a look at SNPs (Single Nucleotide Polymorphisms). SNPs are single nucleotide genomic variations that can be found in a population. Click on the "dbSNP" track from the "Variants" category to open it's configuration.

Select the "full" view for the "Common dbSNP(155)" track, and click submit on the top of the screen. You should now be taken back to the display, and have a new track called "Common (1000 Genomes ... dbSNP Release 155).

Different types of SNPs are visualised in different colours.

<details>
<summary>What are the black, red, and green SNPs?</summary>
- Black: Intronic SNPs (SNPs in introns)
- Red: Missense SNPs (SNPs in a coding region that alters the amino acid of the codon)
- Green: Synonymous SNPs (SNPs in a coding region that does not alter the amino acid of the codon)
</details>

<details>
<summary>Of the 3 type of SNPs discussed in the previous question, which one is most likely to cause disease?</summary>
The Missense SNP (red). These SNPs change the amino acid sequence of a protein, leading to a possible change in the resulting protein, causing disease.
</details>

Now, we will have a look at adding transcripts from another organism. Search for the "Other RefSeq" track in the "Genes & Gene predictions" category. Change the display mode to "full", and click submit.

You should now see an extra track with non-human RefSeq genes.
> You might have to zoom out to see the whole gene again

Click one the first non-human gene to open the detailed view.

<details>
<summary>From what organisms is the gene?</summary>
Mouse (*Mus musculus*)
</details>

<details>
<summary>How similar are the mouse sequence from the human sequence for this gene?</summary>
84.3%: see mRNA/Genomic Alignments section
</details>

Go back to the browser view and zoom out further to see the neighbouring genes of ARTN.

<details>
<summary>What are the upstream and downstream flanking genes</summary>
- Upstream: ST3GAL3
- Downstream: IPO13
> LINC02918 and ENSG00000288573 are here not counted as up- or downstream, since they are located on a different strand
</details>

<details>
<summary>Right next to ARTN lies the LINC02918 gene. Why is this gene coloured in green?</summary>
It is coloured in green because it is a non-coding gene (ncRNA or non-coding RNA), and does not produce a protein product.
</details>

"Synteny" is a concept in biology that means that the gene order is conserved between species. A "syntenous" region is a region which has the same layout in two different species.

<details>
<summary>Is the region containing ARTN and its two neighbours syntenous between human and mouse?</summary>
Yes. In mouse the ARTN gene is also flanked by ST3GAL3 upstream, and IPO13 downstream.
> Note: while you might come to this conclusion just by seeing mouse transcripts where the human IPO13 and ST3GAL3 are, this does not mean these transcripts are in the same order on the mouse genome. You have to go to the actual mouse genome browser and look for ARTN to confirm this.
</details>

Now, check out the KCNJ11 gene, and activate the OMIM Alleles track to "full" (can be found in the "Phenotype and Literature category"

<details>
<summary>What is the OMIM database? What kind of information can be found there?</summary>
OMIM ([Online Mendelian Inheritance in Man](https://www.omim.org/help/about)) is a databse of human genes and genetic phenotypes. It categorizes so-called Mendelian diseases, diseases caused by changes in a single gene. 
</details>

<details>
<summary>Which diseases are associated with the OMIM allelic variants in KCNJ11?</summary>
- Hyperinsulinemic hypoglycemia
- Maturity-onset diabetes of the young
- Diabetes mellitus
</details>

Genes are often regulated by regulatory elements close to the actual gene sequence. One way to find such regions, is by looking for "open chromatin". Normally, DNA is bound to histones and coiled up into nucleosomes.
However, for proteins to have access to the DNA (for example for replication, transcription, or regulation), the chromatin needs to be open, i.e. the DNA needs to be accessible.
Open regions in the chromatin can be detected by so-called DNAse essays. These assays use DNAse (a DNA-cutting enzyme) that can access this open chromatin, and cut the DNA there. 
In DNAse assays, we can detect which parts of the DNA is accessible by detecting where DNAse can bind to the DNA.

Go to the configuration of the "ENCODE regulation tracks" in the "Regulation category". Activate the 3 DNase tracks to full, and click submit. Zoom out a bit, so you seen the surroundings of the KCNJ11 gene.

<details>
<summary>How many "open" chromatin regions can be detected near or in KCNJ11?</summary>
3: One before, one after, and one in the middle of the gene
</details>

Regulatory elements often contain so-called CpG islands (regions with a lot of CG's). Activate the "CpG" track in the "Regulation" category to full as well

<details>
<summary>Which of the regions with high DNase activity also has CpG islands?</summary>
The region upstream (before) the start of the gene.
</details>

Have a look at the repeats in this region (if the repeatmasker track is not on by default, activate it in the "Repeats" category).

<details>
<summary>What kind of repeat do we find in the regulatory region with the CpG islands?</summary>
Repeat (CGCCCG)n
</details>

In the "Regulation" category, activate the ORegAnno track. [ORegAnno](https://www.bcgsc.ca/resources/software/oreganno) provides a curated resource of gene regulatory information. The track provides evidence of proteins binding to a certain region in the genome.

<details>
<summary>Are there any proteins that bind to our supsected regulatory region? If so, which ones?</summary>
Yes:
- EGR1 (transcription factor)
- CTCF (transcription factor)
> You can find this information by clicking on the ORegAnno boxes
</details>

## IGV 

UCSC is a usefull genome browser for browsing public reference genome and their associated data.
However, it does not allow to visualise your own data. To do this, we will be using the Integrated Genome Browser (IGV).

If you haven't alredy, [install the IGV genome browser](http://software.broadinstitute.org/software/igv/).

Open the genome browser, and select the human reference genome (GRCh38/hg38) from the top-left dropdown menu. This will download and load the genome for you. You should see the chromosomes (1->Y) on top, and the RefSeq genes on the bottom (blue bars).

Now we need some data to visualize. We will use public ChIP-Seq data from [ReMap](https://remap.univ-amu.fr/). ReMap is a database of ChIP-seq and ChIP-exo data for the most common reference species (human, mouse, fruitfly, arabidopsis).
ChIP-seq and ChIP-exo are methods use to find the binding of specific proteins to the genomes. Proteins are added to the genome, where they will bind to their specific binding sites. The proteins are then "glued" to the DNA, and the DNA that is stuck to the proteins is than sequenced.
This allows us to detect at which specific sequence certain proteins (e.g. transcription factors) bind.

On the ReMap website, go to "Download", and scroll down to the ReMap data by Target. In the search box, look for the GABPA transcription factor. Download the ***Merged peaks*** file for the human genome (hg38). This should download a `.bed.gz` file.
> `.bed` files are text files which contain genomic regions as coordinates. They are always specific to a certain reference sequence (in this case hg38). The `.gz` extension is just a compression extension (from [Gzip](https://www.gnu.org/software/gzip/)).

Now load the downloaded file into IGV. Go to file -> load from file, and select the `.bed.gz` file (no need to decompress it first). A new track should now be displayed (remap2022_GABPA_...).

Now in the searchbar, search for the BRCA1 gene. You should see the similar gene display as in the UCSC browser with genes in blue. 

<details>
<summary>Does GABPA bind to the BRCA1 gene sequence? If so, in what tissue?</summary>
Yes, in liver tissue (GABPA:liver block in the remap track)
</details>

Zoom out and scroll a bit left and right around the BRCA1 gene.

<details>
<summary>Does GABPA bind to many other genes? What does this mean?</summary>
Yes, GABPA seems to bind to many other genes. This could mean that GABPA is a very general transcription factor, or that there is was lot of off-target binding and/or background noise during the experiments.
</details>

We will use the IGV browser more in the following tutorials.

## Optional: Ensembl genome browser

If you’re interested in the Ensemble genome browser, you can find tutorials about it on [https://training.ensembl.org/exercises](https://training.ensembl.org/exercises)
