seqDir="/data/mitsuki/out/altorf/genome/fasta"
gffDir="/data/mitsuki/data/refseq/genomic_gff"
outDir="/data/mitsuki/out/altorf/genome/property/relative"

basename=${1}
seqFilepath=${seqDir}/${basename}_chromosome.fna
gffFilepath=${gffDir}/${basename}_genomic.gff
outFilepath=${outDir}/${basename}_chromosome_relative.csv

~/altorf/genome/property/relative_overlap.py ${seqFilepath} ${gffFilepath} ${outFilepath}

seqFilepath=${seqDir}/${basename}_chromosome_simu03_0.fna
outFilepath=${outDir}/${basename}_chromosome_sim_relative.csv

~/altorf/genome/property/relative_overlap.py ${seqFilepath} ${gffFilepath} ${outFilepath}

