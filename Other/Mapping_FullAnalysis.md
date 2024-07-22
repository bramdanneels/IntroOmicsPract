# Mapping - Full analysis

During the [mapping practical](../Practicals/04_ReadMapping.md) we have mapped several reads to the _A. thaliana_ genome.
However, we have used only a subset of reads.
In this analysis we want to show you wha the analysis and mapping statistics would look like when the full datasets would have been used.
We will only document the used code, and provide the most relevant results, so you can compare to the results from the practical yourself.

## QC and trimming

The different read sets have been trimmed of adapters, poly-A tails, and low-quality bases:

- The Illumina DNA reads have been trimmed using cutadapt (`-q 30,30`; `-m 30`; `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`; `-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`)
- The Illumina RNA reads have been trimmed using cutadapt (`-q 30,30`; `-m 30`; `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`; `-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`; `--poly-a`)
- The Nanopore reads have been adapter trimmed using Porechop, and quality trimmed using NanoFilt (see [here](Assembly_FullAnalysis.md) for the commands and references)
- The HiFi reads have been adapter trimmed using HiFiAdapterFilt (see [here](Assembly_FullAnalysis.md) for the commands and reference)
- The long read RNA sequences have been additionaly poly-A trimmed using cutadapt.

## Short read mappers - timing

Here we timed the tree mappers on the full Illumina Paired-End dataset (2x8.7M instead of 2x1M reads).
We also timed the indexing part (separately) for your information. The same number of threads (64) were used for all 3 mappers.

```
time smalt index -k 15 -s 8 TAIR_smalt.index TAIR.fasta
real    0m11.376s
user    0m9.831s
sys     0m0.529s

time smalt map -n 64 TAIR_smalt.index Atha_DNA_PE_F.fastq.filt Atha_DNA_PE_R.fastq.filt > DNA_PE_smalt.sam
real    2m7.585s
user    103m58.518s
sys     1m36.876s

time bwa-mem2 index -p TAIR_bwa.index TAIR.fasta
real    0m47.933s
user    0m44.729s
sys     0m1.724s

time bwa-mem2 mem -t 64 TAIR_bwa.index Atha_DNA_PE_F.fastq.filt Atha_DNA_PE_R.fastq.filt > DNA_PE_bwa.sam
real    1m21.938s
user    57m29.970s
sys     0m58.755s

time hisat2-build TAIR.fasta TAIR_hisat.index
real    1m14.870s
user    1m12.231s
sys     0m0.772s

time hisat2 -p 64 --no-spliced-alignment -x TAIR_hisat.index -1 Atha_DNA_PE_F.fastq.filt -2 Atha_DNA_PE_R.fastq.filt -S DNA_PE_hisat.sam
real    3m11.345s
user    188m15.460s
sys     12m57.816s
```

The resulting `.sam` files were convered to `.bam`, and the mapping statistics were calculated (results printed below the commands):

```
samtools view -bS DNA_PE_smalt.sam | samtools sort -@32 > DNA_PE_smalt.bam && rm DNA_PE_smalt.sam
samtools view -bS DNA_PE_bwa.sam | samtools sort -@32 > DNA_PE_bwa.bam && rm DNA_PE_bwa.sam
samtools view -bS DNA_PE_hisat.sam | samtools sort -@32 > DNA_PE_hisat.bam && rm DNA_PE_hisat.sam

for map in *bam; do \
echo $map;
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}';
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}';
echo "% reads mapped:" && samtools flagstat $map | grep "mapped (";
echo ""
done

DNA_PE_bwa.bam
Average coverage:
11.369
% of genome covered:
93.4738
% reads mapped:
14520810 + 0 mapped (95.95% : N/A)
14492295 + 0 primary mapped (95.94% : N/A)

DNA_PE_hisat.bam
Average coverage:
10.2055
% of genome covered:
90.2865
% reads mapped:
15887580 + 0 mapped (88.10% : N/A)
12959696 + 0 primary mapped (85.79% : N/A)

DNA_PE_smalt.bam
Average coverage:
11.3978
% of genome covered:
93.7138
% reads mapped:
14561488 + 0 mapped (96.40% : N/A)
14561488 + 0 primary mapped (96.40% : N/A)
```
> The `-@32` parameter of samtools tells it to use 32 threads for sorting the bam files.

We notice that BWA-mem is the fastest mapper on larger datasets. The mapping statistics are better here: 10x coverage, and most of the genome covered.
Similar to the practical, BWA-mem and Smalt have similar statistics, while HiSat is a bit less good.

## ChIP Seq Mapping

Here we map the full ChIP Seq read data to the reference genome (2x7.3M reads instead of 2x1M), and calculate the mapping statistics:
>We use BWA-mem here because it is faster on large datasets.

```
bwa mem -t 64 TAIR_bwa.index Atha_ChIP_PE_F.fastq.filt Atha_ChIP_PE_R.fastq.filt \
| samtools view -bS | samtools sort -@32 > CHIP_smalt.bam

samtools depth -a CHIP_smalt.bam | awk '{c++;s+=$3}END{print s/c}' #Coverage
5.77621
samtools depth -a CHIP_smalt.bam | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}' #Proportion genome covered
92.1482
samtools flagstat CHIP_smalt.bam | grep "mapped (" #Mapped reads
13916178 + 0 mapped (97.63% : N/A)
13916152 + 0 primary mapped (97.63% : N/A)
```

We notice that the ChIP-seq mapping has lower coverage, but a similar genome coverage.
This is a bit unexpected given that we expected only a part of the genome to be covered.
However, whith ChIP-seq we tend to search for "overrepresented" regions in the genome: 
regions that have higher-than-average coverage compared to the rest of the genome.
In normal genome sequencing we expect that the coverage will be roughly the same over the whole genome.
During ChIP-seq we assume that regions where the protein of interest binds, will be enriched in the sample and we will higher coverage there.

## Long Read DNA Mapping

Here we map the full long read sets to the reference.
The Nanopore set has 80k reads instead of 10k reads,
and the HiFi set has 260k reads instead of 10k reads.

```
minimap2 -t 64 -ax map-hifi TAIR.fasta Atha_DNA_HiFi.fastq.filt | samtools view -bS | samtools sort -@32 > DNA_HiFi.bam
minimap2 -t 64 -ax map-ont TAIR.fasta Atha_DNA_Nano.fastq.filt | samtools view -bS | samtools sort -@32 > DNA_Nano.bam

for map in DNA_HiFi.bam DNA_Nano.bam; do \
echo $map;
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}';
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}';
echo "% reads mapped:" && samtools flagstat $map | grep "mapped (";
echo ""
done

DNA_HiFi.bam
Average coverage:
32.7811
% of genome covered:
93.3224
% reads mapped:
635706 + 0 mapped (99.96% : N/A)
255786 + 0 primary mapped (99.91% : N/A)

DNA_Nano.bam
Average coverage:
9.15059
% of genome covered:
99.5027
% reads mapped:
126457 + 0 mapped (93.90% : N/A)
69027 + 0 primary mapped (89.37% : N/A)
```

We can see high coverage in both long read sets. We see a lower mapping rate for Nanopore, but it covers more of the genome.
The lower mapping rate is because there are still a lot of errors in the reads.
The difference in % of the genoem covered could be a difference in the sequenced strain of _A. thaliana_.

## Short Read RNA Mapping

Here we map 2x14M short RNA reads to the reference instead of 2x1M in the practical:

- Using HiSat2

```
hisat2_extract_exons.py TAIR.gtf > TAIR.exons
hisat2_extract_splice_sites.py TAIR.gtf > TAIR.ss
hisat2-build -p 64 --exon TAIR.exons --ss TAIR.ss TAIR.fasta TAIR_hisat_splice.index

hisat2 -t 64 -x TAIR_hisat_splice.index -1 Atha_RNA_PE_F.fastq.filt -2 Atha_RNA_PE_R.fastq.filt -S RNA_PE_hisat_splice.sam
samtools view -bS RNA_PE_hisat_splice.sam | samtools sort -@32 > RNA_PE_hisat_splice.bam && rm RNA_PE_hisat_splice.sam

map=RNA_PE_hisat_splice.bam
echo $map
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}'
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
echo "% reads mapped:" && samtools flagstat $map | grep "mapped ("

RNA_PE_hisat_splice.bam
Average coverage:
21.1383
% of genome covered:
31.3839
% reads mapped:
26615503 + 0 mapped (96.77% : N/A)
25853388 + 0 primary mapped (96.68% : N/A)
```

- Using STAR

```
mkdir STAR_in && ln -s $(pwd)/TAIR.fasta STAR_in/
STAR --runThreadN 64 --runMode genomeGenerate --genomeDir ./STAR_in \
--genomeFastaFiles ./STAR_in/TAIR.fasta --sjdbGTFfile TAIR.gtf \
--sjdbOverhang 100 --genomeSAindexNbases 10

STAR --runThreadN 64 --genomeDir ./STAR_in --sjdbGTFfile  TAIR.gtf --sjdbOverhang 100 --readFilesIn Atha_RNA_PE_F.fastq.filt Atha_RNA_PE_R.fastq.filt
samtools view -bS Aligned.out.sam | samtools sort -@32 > RNA_PE_STAR.bam && rm Aligned.out.sam
rm -r Log.* _STARgenome SJ.out.tab

map=RNA_PE_STAR.bam
echo $map
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}'
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
echo "% reads mapped:" && samtools flagstat $map | grep "mapped ("

RNA_PE_STAR.bam
Average coverage:
21.3444
% of genome covered:
31.6077
% reads mapped:
26820348 + 0 mapped (100.00% : N/A)
26197937 + 0 primary mapped (100.00% : N/A)
```

We can see that the average coverage is quite high considering the low proportion of the genome that is covered.
As discussed, this is because we only have reads mapping to genes, which are a minority of the genome (compared to non-coding DNA).
Since we have an average coverage of 21x on only 31.5% of the genome, this means we have quite high coverage per gene.

## Long Read RNA Mapping

Here we map 1.2M Nanopore RNA reads instead of 100k reads.

```
paftools.js gff2bed TAIR.gtf > TAIR.bed
minimap2 -t 64 -ax splice -uf -k14 --junc-bed TAIR.bed TAIR.fasta Atha_RNA_Nano.fastq.filt | samtools view -bS | samtools sort -@32 > RNA_NanoPore.bam

map=RNA_NanoPore.bam
echo $map
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}'
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
echo "% reads mapped:" && samtools flagstat $map | grep "mapped ("

RNA_NanoPore.bam
Average coverage:
6.96772
% of genome covered:
27.9814
% reads mapped:
885333 + 0 mapped (99.87% : N/A)
860088 + 0 primary mapped (99.87% : N/A)
```

The coverage is of the long reads is lower (we had fewer total bp), and the proportion of the genome covered is also lower.
The lower proportion covered could be due to having not enough reads, but is more likely due to different experimental conditions, 
leading to different (and probably less) genes expressed.

## Wrong mappings

We will do some of the "wrong" mappings on the full dataset as well.

First, we'll map the short RNA-reads using a genomic mapper (bwa-mem2):

```
bwa-mem2 mem -t 64 TAIR_bwa.index Atha_RNA_PE_F.fastq.filt Atha_RNA_PE_R.fastq.filt \
| samtools view -bS | samtools sort -@32 > RNA_bwa.bam
map=RNA_bwa.bam
echo $map
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}'
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
echo "% reads mapped:" && samtools flagstat $map | grep "mapped ("

RNA_bwa.bam
Average coverage:
21.199
% of genome covered:
31.8838
% reads mapped:
30138125 + 0 mapped (99.38% : N/A)
26552549 + 0 primary mapped (99.30% : N/A)
```

We can see that the mapping statistics are very similar to using a splice-aware aligner.

And we will also have a look at mapping the wrong species as well (using bwa-mem2).
Here we use 2x3.3M instead of 2x1M reads.

```
bwa-mem2 mem -t 64 TAIR_bwa.index Bnap_DNA_PE_F.fastq.filt Bnap_DNA_PE_R.fastq.filt \
| samtools view -bS | samtools sort -@32 > Bnap_bwa.bam
map=Bnap_bwa.bam
echo $map
echo "Average coverage:" && samtools depth -a $map | awk '{c++;s+=$3}END{print s/c}'
echo "% of genome covered:" && samtools depth -a $map | awk '{c++; if ($3>0) total+=1}END{print (total/c)*100}'
echo "% reads mapped:" && samtools flagstat $map | grep "mapped ("

Bnap_bwa.bam
Average coverage:
2.51476
% of genome covered:
21.1592
% reads mapped:
2653374 + 0 mapped (39.93% : N/A)
2631666 + 0 primary mapped (39.74% : N/A)
```

As you can see the mapping is really bad, with only 40% of reads mapped, and only 21% of the genome covered.
Mapping with smalt gives a bit better results (3x Cov, 33% of genome, 50% reads mapped), but is still quite bad.

## Visualization

Here you can find some screenshot of the mappings in IGV for your reference.
All mappings are focussed on the same region on Chr2 of _A. thaliana_.

DNA short read mapping using 3 differen mappers:

![DNA_PE_Map.PNG]

DNA long reads mapping using minimap2

