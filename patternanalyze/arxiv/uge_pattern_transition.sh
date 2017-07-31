input="../speciespick/picked_assembly_summary.csv"
#input="../speciespick/head_picked.txt"

array=(`awk -F, 'NR>1 {print $4}' ${input}`) 

logDir="/home/mitsuki/out/altorf/genome/patternanalyze/log"

for basename in ${array[@]}
do
	outlog=${logDir}/sge_${basename}.out
	errlog=${logDir}/sge_${basename}.err
	qsub -S /bin/bash -q all.q -o ${outlog} -e ${errlog} ./pattern_transition.sh ${basename}
done
