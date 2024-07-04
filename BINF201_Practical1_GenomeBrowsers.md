
# BINF201 - Practical 1 - Genome browsers
In this first practical exercise we will have a look at different genome browsers, and familiarize with them. The main aim is to make you comfortable using genome browsers, as we will use genome browsers throughout the practical exercises to visualize some of our results. 

## USCS Genome Browser

First, we'll have a look at one of the most detailed online genome browsers: the UCSC genome browser

Go to the [UCSC genome browser webpage](https://genome.ucsc.edu/)
Go to the latest version of the human genome (GRCh38/hg38)
Find the BRCA1 gene

1. The gene models (in dark blue) have 3 types of lines (fine with arrows, medium, large). What do these 3 different lines represent?

2.	Open the info panel for BRCA1 by clicking on one of the transcripts. What is the function of BRCA1?

3.	Where in the genome is BRCA1 located?

4.	BRCA1 located on the forward or reverse strand? (compared to the reference)

5. Find the start codon (= translation start) of the first transcript, and try zooming in on it. In what color are start codons depicted on the sequence?

6. Zoom in further until you can see the nucleotide sequence. What are the first three codons, and for which amino acids do they code?
(hint: remember that strand that the gene is on)

Find the ARTN gene
Add the “All SNPs” track from the newest dbSNP (v155), 

> Tracks can be added, removed, and customized below the browser window. You can find the dbSNP 155 options in the "Variants" section. If necessary, you can edit an opened track by clicking on the grey bar on the left. You can also expand or collapse tracks by clicking on the track name, which might be necessary to see more information. Once you made the changes, you need to click "Submit" to go back to the browser.

7.	What color do non-synonymous SNPs have?

8. What kind of SNP is rs2242637? 

> You can use the search tool if you don’t find it immediately in the list

9. What about rs12737332?

Add mouse transcripts (hint: check the Other Refseq track in "Genes & Gene predictions), and open a mouse transcript that aligns with the human ARTN gene

10. What is the similarity (in % identities) between this transcript and the human ARTN?
 
Go back to the browser and zoom out a bit.

11. What are the neighboring genes of ARTN? (disregarding genes displayed in green)

12. Look up what “Synteny” means. Is there synteny between the neighbors of ARTN in human and mouse?

Find the KCNJ11 gene, and expand the OMIM Alleles track (click on the track title)

13. What is OMIM?

14. With what diseases are the OMIM alleles associated with in KCNJ11? (Give 2 that you find)

DNAse is an enzyme that degrades DNA. It only degrades DNA that is not “protected”, i.e. DNA that is not protected by histones. DNAse assays can be used to survey the genome for “open” genome regions, which often correspond to promoter regions. 

Turn on the DNAse signal (full) in the “ENCODE Regulation” track.

15. Based on the DNAse signal, where is the promoter region for KCNJ11?

Another typical feature of promotors, are the presence of CpG islands. CpG islands are regions in the genomes with higher GC content, and a higher-than-expected presence of CG dinucleotides. These sites can be easily methylated, and as such can regulate gene expression. Enable the "CpG Islands" track under "Regulation"

16. Does the presumed promotor side of KCNJ11 contain CpG islands?

17. What type of repeat element(s) can be found in the presumed KCNJ11 promotor? (hint: expand the track if necessary)

Look at (or add if necessary) the ORegAnno track. ORegAnno provides a curated resource for gene regulatory information.

18. According to the ORegAnno information, what transcription factors bind to the presumed KCNJ11 promotor?

## IGV 

If not the case yet,  [install the IGV genome browser](http://software.broadinstitute.org/software/igv/).
Select the latest human genome reference (GRCh38/hg38)


Go to [https://remap.univ-amu.fr/](https://remap.univ-amu.fr/)

1.	What is the REMAP project, and what kind of information does it collect?

2.	What is Chip-Seq?

Go to the download page, and search for GABPA. Download the merged peaks file for hg38. This should download a .bed.gz file.

3.	What is a bed file, and what information does it contain?

Load the file into IGV
Search for the BRCA1 gene

4.	Does GABPA bind to the BRCA1 promotor?

5.	Can you find 2 other genes that are likely regulated by GABPA?

## Optional: Ensembl genome browser

If you’re interested in the Ensemble genome browser, you can find tutorials about it on [https://training.ensembl.org/exercises](https://training.ensembl.org/exercises)