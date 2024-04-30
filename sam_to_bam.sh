#!/bin/bash

# Define the path to samtools
samtools_path="/usr/bin/samtools"

# Define the path to the directory containing SAM files
sam_dir="/student/nchenchu/courses/bcb5250/RNAseq_protocol/chrX_data/hisat_out"

# Change directory to the directory containing SAM files
cd "$sam_dir" || exit

# Sort and convert each SAM file to BAM
$samtools_path sort -@ 8 -o ERR188044_chrX.bam ERR188044_chrX.sam
$samtools_path sort -@ 8 -o ERR188104_chrX.bam ERR188104_chrX.sam
$samtools_path sort -@ 8 -o ERR188234_chrX.bam ERR188234_chrX.sam
$samtools_path sort -@ 8 -o ERR188245_chrX.bam ERR188245_chrX.sam
$samtools_path sort -@ 8 -o ERR188257_chrX.bam ERR188257_chrX.sam
$samtools_path sort -@ 8 -o ERR188273_chrX.bam ERR188273_chrX.sam
$samtools_path sort -@ 8 -o ERR188337_chrX.bam ERR188337_chrX.sam
$samtools_path sort -@ 8 -o ERR188383_chrX.bam ERR188383_chrX.sam
$samtools_path sort -@ 8 -o ERR188401_chrX.bam ERR188401_chrX.sam
$samtools_path sort -@ 8 -o ERR188428_chrX.bam ERR188428_chrX.sam
$samtools_path sort -@ 8 -o ERR188454_chrX.bam ERR188454_chrX.sam
$samtools_path sort -@ 8 -o ERR204916_chrX.bam ERR204916_chrX.sam


