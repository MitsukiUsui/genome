import numpy as np
import pandas as pd
import sys
from Bio import SeqIO

def codon_dir(codon):
	codon=codon.upper()
	if codon in {"TAA","TAG","TGA"}:#stop codon for foward lane
		return 1
	elif codon in {"TTA","CTA","TCA"}:#stop codon for backward lane
		return -1
	else:
		return 0
	
def main(inFilepath,outFilepath):
	lane_lst=["+1","+2","+3","-1","-2","-3"]
	dct_lst=[]
	for seqRecord in SeqIO.parse(inFilepath, "fasta"):
		for i in range(2,len(seqRecord)):
			codon=seqRecord[i-2:i+1].seq
			codonDir=codon_dir(str(codon))

			if codonDir!=0:
				dct={}
				dct["Chromosome"]=seqRecord.name
				dct["Start"]=i-2
				dct["End"]=i+1
				for lane in lane_lst:
					dct[lane]=0

				if codonDir==1:
					dct["Feature"]=str(codon)
					dct['+'+str((i-2)%3+1)]=1
				else:
					dct["Feature"]=str(codon.reverse_complement())
					dct['-'+str((i-2)%3+1)]=1
				dct_lst.append(dct)

	out_df=pd.DataFrame(dct_lst)
	out_df=out_df[["Chromosome","Start","End","Feature"]+lane_lst]
	out_df.to_csv(outFilepath, sep='\t', index=False)

if __name__=="__main__":
	#filepath="/home/mitsuki/out/altorf/genome/fasta/GCF_000265365.1_ASM26536v1_chromosome.fna"
	main(sys.argv[1],sys.argv[2])
