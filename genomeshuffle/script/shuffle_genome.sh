input="../../speciespick/picked_assembly_summary.csv"
#input="../../speciespick/head_picked.txt"

array=(`awk -F, 'NR>1 {print $4}' ${input}`) 

seqDir="/home/mitsuki/out/altorf/genome/fasta/"
gffDir="/home/mitsuki/data/refseq/genomic_gff/"
shuffleMode1=0
shuffleMode2=3

for basename in ${array[@]}
do
    seqFilepath=${seqDir}${basename}_chromosome.fna
	gffFilepath=${gffDir}${basename}_genomic.gff
	outFilepath=${seqDir}${basename}_chromosome_simu${shuffleMode1}${shuffleMode2}_0.fna
	echo ${basename}
	../build/shuffle_genome ${seqFilepath} ${gffFilepath} ${outFilepath} ${shuffleMode1} ${shuffleMode2}
	echo
done
