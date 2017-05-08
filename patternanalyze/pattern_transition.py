#!/home/mitsuki/.pyenv/versions/anaconda3-2.5.0/bin/python3.5

import numpy as np
import pandas as pd
import sys
from Bio import SeqIO

##============================================================
##update orfCounter given codon on position i
##============================================================
def update_orfCounter(orfCounter,codon,i):	#{{{
	def get_lane_to_be_updated(i):
		if i%3==0:
			return (1,4)
		elif i%3==1:
			return (2,5)
		elif i%3==2:
			return (0,3)
	
	def codon_dir(codon):
		if codon in {"TAA","TAG","TGA"}:#stop codon for foward lane
			return 1
		elif codon in {"TTA","CTA","TCA"}:#stop codon for backward lane
			return -1
		else:
			return 0
	
	fl,bl=get_lane_to_be_updated(i)#forward lane, backward lane
	codonDir=codon_dir(codon)
	
	if(codonDir==0):
		orfCounter[fl]+=1
		orfCounter[bl]+=1
	elif(codonDir==1):
		orfCounter[fl]=0
		orfCounter[bl]+=1
	elif(codonDir==-1):
		orfCounter[fl]+=1
		orfCounter[bl]=0

	return orfCounter#}}}


##============================================================
##Given orfCounter, return  Naive Pattern of NSRF
##============================================================
def get_pattern(orfCounter,windowSize_aa):#{{{
	pattern=""
	for i in range(6):
		if orfCounter[i]>=windowSize_aa:
			pattern+='1'
		else:
			pattern+='0'
	return pattern#}}}


##============================================================
##Given Naive Pattern of NSRF, return Relative Pattern of NSRF
##============================================================
def degenerate_pattern(pattern):#{{{
	def rotate(pattern):
		newPattern=""
		for i in [2,0,1,5,3,4]:
			newPattern+=pattern[i]
		return newPattern
	
	def reverse(pattern):
		newPattern=""
		for i in reversed(range(6)):
			newPattern+=pattern[i]
		return newPattern
				
	def score(pattern):
		return int(pattern,2)
			
	maxVal=-1
	for i in range(6):
		if score(pattern)>maxVal:
			maxVal=score(pattern)
		pattern=rotate(pattern)
		if i==2:
			pattern=reverse(pattern)
		
	return format(maxVal,'06b')#}}}


def get_score(pattern):
	assert len(pattern)==6
	score=0
	for i in range(6):
		if pattern[i]=='1':
			score+=1
	return score


patterns_lst=['000000',
			 '100000','010000','001000','000100','000010','000001',
			 '110000','011000','101000','000110','000011','000101',
			 '100100','010010','001001','100010','010001','001100','100001','010100','001010',
			 '111000','000111',
			 '110100','011010','101001','110010','011001','101100','110001','011100','101010',
			 '100110','010011','001101','010110','001011','100101','001110','100011','010101',
			 '111100','111010','111001','100111','010111','001111',
			 '110110','011011','101101','110011','011101','101110','110101','011110','101011',
			 '111110','111011','111101','110111','011111','101111',
			 '111111']
patternEncode={patterns_lst[i]:i for i in range(64)}


patternDegs_lst=['000000',
				 '100000',
				 '110000',
				 '100100','100010','100001',
				 '111000',
				 '110100','110010','110001',
				 '111100',
				 '110110','110011','110101',
				 '111110',
				 '111111']

patternDegEncode_basic={patternDegs_lst[i]:i for i in range(16)}
patternDegEncode={pattern:patternDegEncode_basic[degenerate_pattern(pattern)] for pattern in patterns_lst}
patternDegScore=[get_score(pattern) for pattern in patternDegs_lst]


##============================================================
##Compare Naive/relative Pattern of NSRF between real genome and simulated genome
##	
##	prefix _p stands for pair, which means _p[real] belongs to real genome and vice versa
##		filepath_p, seqRecords_p, seqRecord_p, orfCounter_p
##============================================================
def main(realFilepath, simuFilepath, windowSize_bp, countFilepath, transitionFilepath, transitiondegenerateFilepath):

	#Configuration
	assert windowSize_bp % 3 == 2
	windowSize_aa=int((windowSize_bp-2)/3)
	print("BEGIN processing...")
	print("\tREAL        : "+realFilepath)
	print("\tSIMU        : "+simuFilepath)
	print("\tWindow size : "+str(windowSize_bp)+" ("+str(windowSize_aa)+")")
	print("\tOUT COUNT   : "+countFilepath)
	print("\tOUT TRANS   : "+transitionFilepath)
	print("\tOUT TRANSDEG: "+transitiondegenerateFilepath)
	
	#define index for real genome and simulated genome
	real = 0
	simu = 1

	#define filepath_p
	filepath_p=[realFilepath, simuFilepath]
	
	#define seqRecords_p, numOfSeqs
	seqRecords_p=[]
	for idx in range(2):
		seqRecords=[]
		for seqRecord in SeqIO.parse(filepath_p[idx], "fasta"):
			seqRecords.append(seqRecord)
		seqRecords_p.append(seqRecords)
	assert len(seqRecords_p[real])==len(seqRecords_p[simu])
	numOfSeqs=len(seqRecords_p[real])


	#define 
	orfCounter_p=[np.zeros(6),np.zeros(6)]
	transitionMat=np.zeros((64,64)).astype(int)
	transitionDegMat=np.zeros((16,16)).astype(int)

	for seqIdx in range(numOfSeqs):
		seqRecord_p=[]
		seqRecord_p.append(seqRecords_p[real][seqIdx])
		seqRecord_p.append(seqRecords_p[simu][seqIdx])
		assert len(seqRecord_p[real])==len(seqRecord_p[simu]), "the lengths of the sequences is not the same ("+str(len(seqRecord_p[real]))+", "+str(len(seqRecord_p[simu]))+")"
		
		#begin processing by 1 bp	
		for i in range(2, len(seqRecord_p[real])):
			pid_p=[]#pid for pattern id
			pdid_p=[]#pdid for pattern degenerate id

			for idx in range(2):
				codon=str(seqRecord_p[idx].seq[i-2:i+1])
				orfCounter_p[idx]=update_orfCounter(orfCounter_p[idx], codon, i)
				pattern=get_pattern(orfCounter_p[idx], windowSize_aa)

				pid_p.append(patternEncode[pattern])
				pdid_p.append(patternDegEncode[pattern])

			if(i>=windowSize_bp-1):#when window can be defined
				transitionMat[pid_p[real]][pid_p[simu]]+=1
				transitionDegMat[pdid_p[real]][pdid_p[simu]]+=1
		

	distDegArr_p=[np.sum(transitionDegMat, axis=simu), np.sum(transitionDegMat, axis=real)]
	countArr_p=[np.zeros(7).astype(int), np.zeros(7).astype(int)]
	for idx in range(2):
		for i, count in enumerate(distDegArr_p[idx]):
			countArr_p[idx][patternDegScore[i]]+=count
	
	count_df=pd.DataFrame(index=range(7))
	count_df["real"]=countArr_p[real]
	count_df["simu"]=countArr_p[simu]
	count_df.to_csv(countFilepath)
	df=pd.DataFrame(transitionMat,columns=patterns_lst,index=patterns_lst)
	df.index.name = 'idx'
	df.to_csv(transitionFilepath)
	deg_df=pd.DataFrame(transitionDegMat,columns=patternDegs_lst,index=patternDegs_lst)
	deg_df.index.name = 'idx'
	deg_df.to_csv(transitiondegenerateFilepath)
	
if __name__=="__main__":
	main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6])
