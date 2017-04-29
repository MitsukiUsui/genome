#input="../speciespick/head_picked.txt"
input="../speciespick/picked_assembly_summary.csv"

array=(`awk -F, 'NR>1 {print $4}' ${input}`) 

indir="/home/mitsuki/data/refseq/genomic_fna/"
outdir="/home/mitsuki/out/altorf/genome/"
for basename in ${array[@]}
do
    seqFilepath=${indir}${basename}_genomic.fna
    chrFilepath=${outdir}fasta/${basename}_chromosome.fna
    plsFilepath=${outdir}fasta/${basename}_plasmid.fna
    logFilepath=${outdir}preprocess/split/${basename}_split.log
    echo ${basename}
    ./fasta_split ${basename} ${seqFilepath} ${chrFilepath} ${plsFilepath} ${logFilepath} 1
#	echo "./fasta_split ${basename} ${seqFilepath} ${chrFilepath} ${plsFilepath} ${logFilepath}"
done
