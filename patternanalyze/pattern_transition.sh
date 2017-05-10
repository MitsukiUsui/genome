seqDir="/home/mitsuki/out/altorf/genome/fasta/"
outDir="/home/mitsuki/out/altorf/genome/patternanalyze/"
shuffleMode1=0
shuffleMode2=3

basename=${1}

realFilepath=${seqDir}${basename}_chromosome.fna
simuFilepath=${seqDir}${basename}_chromosome_simu${shuffleMode1}${shuffleMode2}_0.fna
windowSize_bp=500
countFilepath=${outDir}count/${basename}_chromosome_simu${shuffleMode}${shuffleMode2}_count.csv
transFilepath=${outDir}trans/${basename}_chromosome_simu${shuffleMode}${shuffleMode2}_trans.csv
transdegFilepath=${outDir}transdeg/${basename}_chromosome_simu${shuffleMode}${shuffleMode2}_transdeg.csv

echo ${basename}
~/altorf/genome/patternanalyze/pattern_transition.py ${realFilepath} ${simuFilepath} ${windowSize_bp} ${countFilepath} ${transFilepath} ${transdegFilepath}
