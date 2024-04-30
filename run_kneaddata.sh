#!/bin/bash
dir="/student/nchenchu/courses/bcb5250/metagenomics_wgs/raw_data"
files="$dir"/*_R1.fastq

echo "Run script ========================="
for FILE in $files; do
echo "$FILE"
# Extract the base name of the file
basename "$FILE"
FILENAME="$(basename -- $FILE)"
echo $FILENAME
# Extract the sample ID
ID=$(echo "${FILENAME}" | sed 's/_R1.fastq//')
echo $ID
# Run kneaddata for the current sample
kneaddata -i1 "$dir"/"$ID"_R1.fastq -i2 "$dir"/"$ID"_R2.fastq -o kneaddata_out -db /public/ahnt/courses/bcb5250/metagenomics_wgs/bowtie2db/hg37dec_v0.1 --trimmomatic /opt/Trimmomatic-0.39 -t 8 --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" --bowtie2 /opt/bowtie2-2.4.2-linux-x86_64 --bowtie2-options "--very-sensitive --dovetail" --remove-intermediate-output --trf /public/ahnt/software/
  done
