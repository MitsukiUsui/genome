#!/bin/bash
#$ -S /bin/bash
#$ -N shuffle
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=1G

taxid=${1}
seqFilepath=/data/mitsuki/out/altorf/genome/fasta/${taxid}_chr.fna
gffFilepath=/data/mitsuki/data/refseq/genomic_gff/${taxid}_genomic.gff
outFilepath=/data/mitsuki/out/altorf/genome/fasta/${taxid}_sim.fna
shuffleMode1=0
shuffleMode2=3

time ../../build/shuffle_genome \
     ${seqFilepath} ${gffFilepath} ${outFilepath} ${shuffleMode1} ${shuffleMode2}

