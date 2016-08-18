#!/bin/bash

#PBS -l nodes=1:ppn=1,walltime=10:00:00

#PBS -N reduced_genome

#PBS -l mem=20GB

module load bedtools
module load kent
cd folder_for_analysis
#fragment length
fl=30
#primary enzyme for file name
enzyme=hindiii
#genome for file name
genome=mm10

#get coordinates for all RE sites in the genome
oligoMatch ${enzyme}.fa ${genome}.fa ${genome}_${enzyme}_restriction_sites_oligomatch.bed
#get coordinates of upstream fragments
awk -v fl=$fl '{print $1"\t"$2-fl"\t"$2}' ${genome}_${enzyme}_restriction_sites_oligomatch.bed >up.txt
#get coordinates of downstream fragments
awk -v fl=$fl '{print $1"\t"$3"\t"$3+fl}' ${genome}_${enzyme}_restriction_sites_oligomatch.bed > down.txt
#combine up and downstream fragments
cat up.txt down.txt > ${genome}_${enzyme}_flanking_sites_${fl}_2.bed
#remove any fragments with negative coordinates
awk '{if($2 >= 0 && $3 >=0) print $0}' ${genome}_${enzyme}_flanking_sites_${fl}_2.bed | grep -v -E 'random|JH|GL' - | sort -k1,1 -k2,2n | uniq  > ${genome}_${enzyme}_flanking_sites_${fl}_unique_2.bed
#get the sequence of unique flanking coordinates
fastaFromBed -fi ${genome}.fa -bed ${genome}_${enzyme}_flanking_sites_${fl}_unique_2.bed -fo ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa
#get only sequences from FASTA file
grep -v '^>' ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa | sort | uniq -i -u | grep -xF -f - -B 1 ${genome}_${enzyme}_flanking_sequences_${fl}_unique_2.fa | grep -v '^--' > ${genome}_${enzyme}_flanking_sequences_${fl}_unique.fa

#remove unwanted intermediate files
rm up.txt
rm down.txt

#make a BED file of unqiue sequences
grep '^>' ${genome}_${enzyme}_flanking_sequences_${fl}_unique.fa > ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed
sed -i 's/>//g' ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed
sed -i 's/:\|-/\t/g' ${genome}_${enzyme}_flanking_sites_${fl}_unique.bed

exit 0;
