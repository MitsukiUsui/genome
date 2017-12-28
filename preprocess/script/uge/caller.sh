#!/bin/bash
#$ -S /bin/bash
#$ -N split
#$ -q standard.q
#$ -cwd
#$ -v PATH
#$ -o ./log/$JOB_ID.out
#$ -e ./log/$JOB_ID.err
#$ -l mem_free=1G

taxid=${1}
seqFilepath=/data/mitsuki/data/refseq/genomic_fna/${taxid}_genomic.fna
chrFilepath=/data/mitsuki/out/altorf/genome/fasta/${taxid}_chr.fna
plsFilepath=/data/mitsuki/out/altorf/genome/fasta/${taxid}_pls.fna
logFilepath=/data/mitsuki/out/altorf/genome/preprocess/split/${taxid}.log

time ../../build/fasta_split -c 2 ${seqFilepath} ${chrFilepath} ${plsFilepath} ${logFilepath}

