#!/usr/bin/bash

mkdir results
mkdir results/fastqc_reports
fastqc SRR_1.fastq SRR_2.fastq --outdir=results/fastqc_reports/

mkdir temp
echo "BWA INDEX"
bwa index reference.fna
echo "FINISHED BWA INDEX"
echo "BWA MEM"
bwa mem reference.fna SRR_1.fastq SRR_2.fastq -o sample.sam
echo "FINISHED BWA MEM"
echo "SAMTOOLS VIEW"
samtools view -b -o sample.bam sample.sam
echo "FINISHED SAMTOOLS VIEW"
echo "SAMTOOLS FLAGSTAT"
samtools flagstat sample.bam > results/flagstat_result.txt

n=1
until [ "$n" -eq 8 ]
do
read line 
n=$((n+1))
done < 'results/flagstat_result.txt'

symbol=d
p1=$(expr index "$line" $symbol)
p2=$[p1+2]
numb=${line:p2:1}

if [ "$numb" -eq 7 ]
then 
echo "OK"
echo "SAMTOOLS SORT AND INDEX"
samtools sort -o sorted.bam sample.bam
samtools index sorted.bam
echo "FINISHED SAMTOOLS SORT AND INDEX"
echo "FREEBAYES"
freebayes -f reference.fna -b sorted.bam --vcf results/sample.vcf
else echo "NOT OK"
fi

