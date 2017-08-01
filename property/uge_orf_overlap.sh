input="../speciespick/picked_assembly_summary.csv"
#input="../speciespick/head_picked.txt"

array=(`awk -F, 'NR>1 {print $4}' ${input}`) 

logDir="/data/mitsuki/out/altorf/genome/property/log"

for basename in ${array[@]}
do
	outlog=${logDir}/${basename}.sgeout
	errlog=${logDir}/${basename}.sgeerr
    cmd=./relative_overlap.sh
	qsub -S /bin/bash -q standard.q -o ${outlog} -e ${errlog} ${cmd} ${basename}
done
