# Long Read Assembly - Full analysis

In [Practical 3](../Practicals/03_GenomeAssembly.md) we have used both short and long reads for assembling a bacterial genome.
However, due to time and computational constraints, we could not properly analyse the long reads, or analyze only a subset of the data.
In this document you can find an example of analyses on the full read datasets for the HiFi reads, and using proper pre- and post-processing for the Nanopore reads.
We will only document the used code, and provide the most relevant results, so you can compare to the results from the practical yourself.

## Nanopore assembly

In the practical we used the whole set of reads, but did not do any long-read specific pre-processing.
Here we will do some pre- and post-processing of the Nanopore reads and the resulting assembly.
These steps are based on [this tutorial by Tim Kahlke](https://timkahlke.github.io/LongRead_tutorials/).

First, we will use a tool called [Porechop](https://github.com/rrwick/Porechop) to remove possible adapters from the reads.

```
porechop -i MycOvi_Nano.fastq -o MycOvi_Nano.chopped.fastq --threads 64
```

Then, we use the [NanoFilt tool](https://github.com/wdecoster/nanofilt) to filter out short reads and remove the first 10 bases (which are often very low quality).

```
NanoFilt -l 500 --headcrop 10 < MycOvi_Nano.chopped.fastq > MycOvi_Nano.filtered.fastq
```
> The `-l 500` option keeps only reads >500bp.
> The `--headcrop 10` option trims the first 10 bases of every read.

The FastQC report of the filtered reads can be found [here](MycOvi_Nano.filtered_fastqc.html).

Then, we assembled the reads using Flye (similar to the tutorials):

```
flye -t 64 --nano-raw MycOvi_Nano.filtered.fastq --genome-size 1m --out-dir Nano_Flye
```
> Similar to Canu, you can set the genome size as well, which can slightly improve assembly quality.

We will also do some post-assembly polishing. We will map the reads back to the assembly (using minimap2), and use [Racon](https://github.com/isovic/racon) to polish.

```
minimap2 -t64 -Nano_Flye/assembly.fasta MycOvi_Nano.filtered.fastq > Flye_Nano.paf
racon MycOvi_Nano.filtered.fastq Flye_Nano.paf Nano_Flye/assembly.fasta > MycOvi_Racon.fasta
```
> We use minimap2 in default mode, which maps long sequences to long sequences. This creates a special alignment format: `.paf`, which can be used in Racon.

Then we compare these assemblies to the same reference genome used in the practical.
We will compare both the Flye output and the Racon output to see if there is an observable difference.
The Quast report can be found [here](MycOvi_Nano_Quast.html).

Lastly, we will run BUSCO on both assemblies. The results are printed below the commands.

```
busco -c 64 -l mycoplasmatales -m geno -i MycOvi_Flye.fasta
    ---------------------------------------------------
    |Results from dataset mycoplasmatales_odb10        |
    ---------------------------------------------------
    |C:24.1%[S:24.1%,D:0.0%],F:37.4%,M:38.5%,n:174     |
    |42    Complete BUSCOs (C)                         |
    |42    Complete and single-copy BUSCOs (S)         |
    |0    Complete and duplicated BUSCOs (D)           |
    |65    Fragmented BUSCOs (F)                       |
    |67    Missing BUSCOs (M)                          |
    |174    Total BUSCO groups searched                |
    ---------------------------------------------------
busco -c 64 -l mycoplasmatales -m geno -i MycOvi_Racon.fasta
    ---------------------------------------------------
    |Results from dataset mycoplasmatales_odb10        |
    ---------------------------------------------------
    |C:39.7%[S:39.7%,D:0.0%],F:28.2%,M:32.1%,n:174     |
    |69    Complete BUSCOs (C)                         |
    |69    Complete and single-copy BUSCOs (S)         |
    |0    Complete and duplicated BUSCOs (D)           |
    |49    Fragmented BUSCOs (F)                       |
    |56    Missing BUSCOs (M)                          |
    |174    Total BUSCO groups searched                |
    ---------------------------------------------------
```

### Conclusion

The statistics are improved compared to the practical, but the BUSCO scores are still very low.
We do however see that both pre-processing and post-processing (polishing) has significantly improved the BUSCO scores.
Doing further rounds of polishing could further increase the assembly quality (to an extent).

These results show that the main problem with our Nanopore reads is the base quality.
If we had short reads from the same sample, we could use those to polish/error-correct the reads and improve our assembly.
Increasing sequencing coverage could also be an option, but we already have 100x coverage (116.3 Mbp of reads for 1Mbp genome).

Since Nanopore reads are so error-prone, it can actually be detrimental to have high coverage.
Higher coverage means more reads that map to a position, and more chance of a wrong base being mapped.
This makes it hard for the assembler to discover what the correct base is.
As an example, subsetting the filtered Nanopore reads to 15000 reads (instead of +- 45000) and assemblying using Flye (without polishing after),
yielded an assembly with already 30% BUSCO completeness.

## PacBio HiFi assembly

During the practical we only used a subset of the reads (50000 out of total reads)
Here we will check if using all the reads would increase the assembly quality.
In addition, we will use a HiFi-appropriate trimming tool for pre-processing.

First we can use [HiFiAdapterFilt](https://github.com/sheinasim/HiFiAdapterFilt) to pre-process the HiFi reads.
Since HiFi reads have very high base quality, there is no real need for quality trimming.


```
hifiadapterfilt.sh -p MycOvi_HiFi -t 64 -o HiFi_Filt
```

Then we will use hifiasm to assemble the reads, and convert the `.gfa` to `.fasta`:

```
hifiasm -t 64 -o MycOvi_HiFi HiFi_Filt/MycOvi_HiFi.filt.fastq.gz
awk '/^S/{print ">"$2;print $3}' MycOvi_HiFi/hifiasm.bp.hap1.p_ctg.gfa > MycOvi_HiFi.fasta
```

The Quast report of this assembly compared to the two Nanopore assemblies and the reference can be found [here](MycOvi_HiFi_Quast.html).

And lastly, we will run BUSCO on this assembly as well:

```
busco -c 64 -l mycoplasmatales -m geno -i MycOvi_HiFi.fasta
    ---------------------------------------------------
    |Results from dataset mycoplasmatales_odb10        |
    ---------------------------------------------------
    |C:99.4%[S:93.7%,D:5.7%],F:0.6%,M:0.0%,n:174       |
    |173    Complete BUSCOs (C)                        |
    |163    Complete and single-copy BUSCOs (S)        |
    |10    Complete and duplicated BUSCOs (D)          |
    |1    Fragmented BUSCOs (F)                        |
    |0    Missing BUSCOs (M)                           |
    |174    Total BUSCO groups searched                |
    ---------------------------------------------------
```

### Conclusion

The HiFi assembly using all data has again very nice assembly statistics and very good BUSCO score.
However, we see that only 58% of the reference genome could be recovered.
Given the high BUSCO score, this means likely that the sample from which the long reads are derived, is very divergent from the reference.
Analysis of the 16S sequence confirms the species as _Mesomycoplasma ovipneumoniae_, so at least it is not a wrong classification.

