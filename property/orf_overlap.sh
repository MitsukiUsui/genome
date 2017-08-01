seqDir="/data/mitsuki/out/altorf/genome/fasta"
outDir="/data/mitsuki/out/altorf/genome/patternanalyze/overlap"

basename=${1}
seqFilepath=${seqDir}/${basename}_chromosome.fna
outFilepath=${outDir}/${basename}_chromosome_overlap.csv

~/altorf/genome/patternanalyze/orf_overlap.py ${seqFilepath} ${outFilepath}

seqFilepath=${seqDir}/${basename}_chromosome_simu03_0.fna
outFilepath=${outDir}/${basename}_chromosome_sim_overlap.csv

~/altorf/genome/patternanalyze/orf_overlap.py ${seqFilepath} ${outFilepath}

