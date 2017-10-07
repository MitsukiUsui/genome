taxid=${1}
seqFilepath=/data/mitsuki/data/refseq/genomic_fna/${taxid}_genomic.fna
chrFilepath=/data/mitsuki/out/altorf/genome/fasta/${taxid}_chr.fna
plsFilepath=/data/mitsuki/out/altorf/genome/fasta/${taxid}_pls.fna
logFilepath=/data/mitsuki/out/altorf/genome/preprocess/split/${taxid}.log

time /home/mitsuki/altorf/genome/preprocess/build/fasta_split ${seqFilepath} ${chrFilepath} ${plsFilepath} ${logFilepath}

