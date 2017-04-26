input="summary_split.csv"
array=(`awk -F, 'NR>1 {print $2}' ${input}`) 

outdir="/home/mitsuki/out/altorf/genome/preprocess/comp/"
for filepath in ${array[@]}
do
    basename=`basename $filepath .fna`
    outFilepath=${outdir}${basename}.comp
    echo ""${basename}
    fatt composition ${filepath} > ${outFilepath}
done
