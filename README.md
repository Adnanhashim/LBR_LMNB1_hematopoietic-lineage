# LBR_LMNB1_hematopoietic-lineage

## Background: 
The nuclear envelope and lamina define the nuclear periphery and are crucial for the spatial organization of chromatin and genome function. Part of the nuclear lamina is the lamin B receptor (LBR), which binds lamin B1 (LMNB1) and tethers heterochromatin to the inner nuclear membrane (INM).

## Results:
Using primary human CD34+ progenitors, as well as mature monocytes and granulocytes, we generated the first global maps of LBR-chromatin interaction. We compared those with LMNB1 chromatin domains and show that 98% of LMNB1 LADs overlapped with LBR-associated chromatin domains (LCDs), while a significant number of LCDs (50 %) exhibit unique, LMNB1-independent chromatin interactions. In all three cell types analyzed, we found that LBR LCDs are more gene rich than common LCD-LAD domains, and that they exhibit a prominent enrichment in 5-hydroxy methylation, while overall representing a mostly transcriptionally chromatin environment.   

## Conclusions:
Taken together, our study provides an important starting point for the systematic characterization of LBR-chromatin interactions.

----------------------------------

## Data Analysis:

### ChIP seq Analysis:

#### Quality control: FASTQC

FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

```bash
fastqc *.fastq -o 1_fastQC
```

#### Alignment to the Reference Genome
Alignment of raw reads to a reference genome is one of the key steps in ChIP seq data analysis. Here I am using Bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment and samtools (http://www.htslib.org/) for reading, editing, sorting and indexing SAM/BAM formats. 

```bash
#!/bin/bash
nohup sh -c 'for fq1 in *R1_001.fastq.gz; do
fq2=${fq1/_R1_/_R2_}
out=`basename $fq1 _R1_001.fastq.gz`
echo $fq1, $fq2, $out

if [[ ! -f "$fq1" ]]; then echo "WRONG"; break; fi
if [[ ! -f "$fq2" ]]; then echo "WRONG"; break; fi

Bowtie2_INDEX=/lsc/common/Adnantools/bowtie_index/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

##alignment
bowtie2 -p8 -x $Bowtie2Index -1 $fq1 -2 $fq2 > ${out}.sam &&

##sam to sortedbam
samtools view -bS ${out}.sam | samtools sort -@ 16 - -o ${out}_sorted.bam 
rm ${out}.sam ; done' > nohup_bowtie2_alignment_sorted_bam.txt
```

#### Removal of duplicate reads

Duplicate reads were removed using samtools v0.1.19 rmdup tool (keeping duplicate reads doesn’t significantly affect LAD detection because of the large size of LADs). MarkDuplicates (Picard) (https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-) can also be used for similar process.
```bash
nohup sh -c 'for i in *sorted.bam; do
echo $i
/lsc/oldsamtools/samtools rmdup -s $i ${i%%sorted.bam}_rmdup.bam; done' > nohup_BAM_to_rmdupBAM.txt 
```

#### ChIP enrichment

plotFingerprint https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html is used for the cumulative enrichment plots that showed good sample quality, enrichment in ChIP samples and no indication of enrichment in input samples.

```bash
#example of Monocytes replicate2
plotFingerprint -b MonocytesLBR_Rep2_chip_rmdup.bam MonocytesLBR_Rep2_input_rmdup.bam-plot MonocytesLBR_rep2_ChIP_enrichment_fingerprint.png
```

![ChIP enrichment Fingerprint](https://github.com/Adnanhashim/LBR_LMNB1_hematopoietic-lineages/blob/master/Screenshot%202020-07-31%2011.36.22.png)

#### Peak Calling

Enriched Domain Detector (EDD) (https://github.com/CollasLab/edd) is used as ChIP-seq peak caller for detection of LBR and LMNB1 megabase domains of enrichment. EDD has a algorithm for analysis of broad (megabases) enrichment domains from ChIP-seq data. EDD is able to identify genomic domains interacting with broadly distributed proteins i.e. Lamin B, Lamin A and Lamin B Receptor (LBR).

```bash
#example of CD34+ replicate1
chr_size=/lsc/common/Adnantools/hg38_sizes.txt
unaligned_regions=/lsc/common/Adnantools/unalignable_regions.bed 
chip= ./CD34_Rep1_chip_rmdup.bam
input=./CD34_Rep1_input_rmdup.bam

edd --write-log-ratios --write-bin-scores $chr_size $unaligned_regions $chip $input outputDir
```

#### Jaccard Statistics

The Jaccard statistic (https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html#:~:text=The%20Jaccard%20statistic%20is%20used,union%20of%20the%20two%20sets.&text=The%20bedtools%20jaccard%20tool%20implements,the%20length%20of%20the%20union.) is used to find the similarity of the two sets based on the intersections between them.

```bash
#example of CD34+ 

bedtools jaccard -a CD34_Rep1_peaks.bed -b CD34_Rep2_peaks.bed 
```

#### Consensus Peaks

Duplicates of each lineage were converted to consensus peaks/LCDs with bedtools intersect 

```bash
#example of CD34+ 
intersectBed -a CD34_Rep1_peaks.bed -b CD34_Rep2_peaks.bed > CD34_consensus_peaks.bed
```

Consensus peaks (cLADs/cLCDs) within hematopiotic lineages.

```bash
multiIntersectBed -header -i Granulocytes_consensus_peaks.bed Monocytes_consensus_peaks.bed CD34_consensus_peaks.bed -names Granulocytes Monocytes CD34 -empty -g /lsc/common/adnantools/hg19_sizes.txt  > multiintersect_consenus_all_with_empty.bed
```

#### Clustering heatmap of jaccard statistic to measure the similarities of LBR and LMNB1 ChIP sequencing data sets

we need to install a tiny script from bedtools for this analysis.
```
curl -O https://s3.amazonaws.com/bedtools-tutorials/web/make-matrix.py
```
To obtain pairwise Jaccard measurements of all bed files (CD34, Granulocytes and Monocytes)

```
parallel "bedtools jaccard -a {1} -b {2} \
         | awk 'NR>1' \
          | cut -f 3 \
          > {1}.{2}.jaccard" \
          ::: `ls *.bed` ::: `ls *.bed`
          
find . | grep jaccard | xargs grep "" | sed -e s"/\.\///" | perl -pi -e "s/.bed./.bed\t/" | perl -pi -e "s/.jaccard:/\t/" > pairwise.jaccardstats.txt
awk 'NF==3' pairwise.jaccardstats.txt | python make-matrix.py > hematopiotic.distance.matrix
cut -f 1 hematopiotic.distance.matrix | cut -f 1 -d "-" | cut -f 1 -d "_" > labels.txt

```

Create heatmap in R
```
library(ggplot2)
library(RColorBrewer)
library(gplots)

distance_matrix <- read.table('hematopiotic.distance.matrix')
labels <- read.table('labels.txt')
ngroups <- length(unique(labels))
jaccard_table <- distance_matrix
jaccard_matrix <- as.matrix(jaccard_table)

pdf('Biclustering_heatmap_LBR.pdf')
heatmap.2(jaccard_matrix, col = brewer.pal(9,"Blues"), margins = c(14, 14), density.info = "none", lhei=c(2, 8), trace= "none", sepcolor="black",colsep=1:ncol(jaccard_matrix),rowsep=1:nrow(jaccard_matrix), sepwidth=c(0.003,0.003))
dev.off()
```
reference: http://quinlanlab.org/tutorials/bedtools/bedtools.html

![Enrichmentoverchromosomes](https://github.com/Adnanhashim/LBR_LMNB1_hematopoietic-lineages/blob/master/Biclustering_heatmap_LBR.png)


#### Overview of cLCDs and cLADs at the chromosomal level

The location of common cLCDs(cLCDs overlapping with cLADs) and unique cLCDs (cLCDs not overlapping to cLADs and unique) binding regions was mapped over the whole genome and illustrated for each chromosome by using chipseeqer (https://www.bioconductor.org/packages/release/bioc/html/ChIPseeker.html) in R . Both, common LCDs and unique cLCDs could be identified for all chromosomes. 

```R
library("ChIPseeker")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
mylist <- list("uniq"= "unique_cLCDs_306_length_MB.bed", "common"="common_cLCDs_cLADs_311_length_MB.bed")
tiff("chromosomal_coverage_common_unique_cLCDs.tiff")
covplot(mylist, weightCol="V6")
dev.off()
```

![Enrichmentoverchromosomes](https://github.com/Adnanhashim/LBR_LMNB1_hematopoietic-lineages/blob/master/LMNB1_LBR_Figure_3_temp.png)


