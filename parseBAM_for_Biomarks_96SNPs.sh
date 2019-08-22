#!/usr/bin/sh

### $1: NA12878_S10_R1_001.bam
bam = $1
samtools view -bS -L Biomarks_96SNPs_hg38.bed $bam > ${bam/.bam/_Biomarks_96SNPs.bam}
