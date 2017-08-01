seqDir="/data/mitsuki/out/altorf/genome/fasta"
outDir="/data/mitsuki/out/altorf/genome/property/overlap"

basename=${1}
seqFilepath=${seqDir}/${basename}_chromosome.fna
outFilepath=${outDir}/${basename}_chromosome_overlap.csv

~/altorf/genome/property/orf_overlap.py ${seqFilepath} ${outFilepath}

seqFilepath=${seqDir}/${basename}_chromosome_simu03_0.fna
outFilepath=${outDir}/${basename}_chromosome_sim_overlap.csv

~/altorf/genome/property/orf_overlap.py ${seqFilepath} ${outFilepath}

