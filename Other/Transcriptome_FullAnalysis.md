# Transcriptome - Full analysis

During the [transcriptome practical](../Practicals/05_Transcriptome.md) we have generated several transcriptomes from reads mapped to the human Chr15.
Here we want to show how the transcriptome analysis would have looked like if we had used all reads, mapping to the whole human genome.

## QC and trimming

The reads have been trimmed as follows:

- The Illumina RNA reads have been trimmed using cutadapt (`-q 30,30`; `-m 30`; `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`; `-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`; `--poly-a`)
- The Nanopore RNA reads have been adapter trimmed using Porechop, quality trimmed using NanoFilt (see [here](Assembly_FullAnalysis.md) for the commands and references)
- The HiFi RNA reads have been adapter trimmed using HiFiAdapterFilt (see [here](Assembly_FullAnalysis.md) for the commands and reference)
- Cutadapt was used to further remove poly-A tails after trimmming the long reads

## Reference assembly

The full set of short reads (2x6.1M reads) were mapped to the full human genome using HiSat2, and we used Stringtie to assemble the transcripts:

```
hisat2_extract_splice_sites.py human.gtf > human.ss
hisat2_extract_exons.py human.gtf > human.exons
hisat2-build -p 64 --exon human.exons --ss human.ss human.fasta human.index

hisat2 -p 64 --dta -x human.index --max-intronlen 5000 \
-1 RNA_PE_F.filt.fastq -2 RNA_PE_R.filt.fastq \
| samtools view -Sb \
| samtools sort -@32 > RNA_PE.bam

stringtie -p 64 -o PE_stringtie.gtf RNA_PE.bam

gffcompare -r human.gff3 -o PE_gffcompare.txt PE_stringtie.gtf && cat PE_gffcompare.txt
```

We also mapped the long Nanopore (3M reads) and HiFi (650k reads) reads using minimap2, and Stringtie to assembly in long-read mode:

```
paftools.js gff2bed human.gtf > human.bed
minimap2 -t 64 -ax splice -uf -k14 --junc-bed human.bed human.fasta RNA_Nano.filt.fastq \
| samtools view -b | samtools sort -@32 > RNA_Nano.bam
stringtie -p 64 -L -o Nano_stringtie.gtf RNA_Nano.bam

minimap2 -t 64 -ax splice:hq -uf --junc-bed human.bed human.fasta RNA_HiFi.filt.fastq \
| samtools view -b | samtools sort -@32 > RNA_HiFi.bam
stringtie -p 64 -L -o HiFi_stringtie.gtf RNA_HiFi.bam

gffcompare -r human.gff3 -o Long_Chr15_gffcompare.txt Nano_stringtie.gtf HiFi_stringtie.gtf
cat Long_Chr15_gffcompare.txt
```

