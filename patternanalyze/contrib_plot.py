import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def get_num_of_nsrf(pattern):
	assert len(pattern)==6	  
	retVal=0
	for i in range(6):
		if pattern[i]=='1':
			retVal+=1
	return retVal


lookupFilepath="../speciespick/picked_assembly_summary_code.csv"
lookup_df=pd.read_csv(lookupFilepath)
print(lookup_df.shape)

compFilepath="../preprocess/out/summary_composition.csv"
comp_df=pd.read_csv(compFilepath)
comp_df=comp_df[comp_df["dna_type"]=="chromosome"]#exclude prasmid
print(comp_df.shape)

tmp_df=pd.merge(lookup_df, comp_df, on="ftp_basename",how="left")
tmp_df=tmp_df[tmp_df["genetic_code"]==11]#!!!TBI!!! failed to handle genetic code=4 for now
print(tmp_df.shape)

dct_lst=[]

for key, row in tmp_df.iterrows():
	filedir="/data/mitsuki/out/altorf/genome/patternanalyze/transdeg"
	filepath="{0}/{1}_chromosome_simu3_transdeg.csv".format(filedir, row["ftp_basename"])
	print("Processing {0}/{1}".format(key+1, tmp_df.shape[0]))
	
	
	df=pd.read_csv(filepath,dtype={'idx':str})
	df=df.set_index('idx') 
	
	total=np.sum(np.sum(df))
	for idx in df.index:
		for clm in df.columns:
			dif=get_num_of_nsrf(idx)-get_num_of_nsrf(clm)
			df.loc[idx,clm]=df.loc[idx,clm]*dif/total
	countDiff=np.sum(np.sum(df))

	sym_df=pd.DataFrame(0,index=df.index, columns=df.columns)
	for i in range(df.shape[0]):
		for j in range(i):
			sym_df.iloc[i,j]=df.iloc[i,j]+df.iloc[j, i]
	
	dct={}
	dct["ftp_basename"]=row["ftp_basename"]
	dct["G+C"]=row["G+C"]
	dct["taxid"]=row["taxid"]
	dct["count_diff"]=countDiff
	for i in range(df.shape[0]):
		for j in range(df.shape[1]):
			name="{0}->{1}".format(i,j)
			dct[name]=sym_df.iloc[i][j]
	dct_lst.append(dct)
	

#for _,basename in enumerate(tmp_df["ftp_basename"]):
#	filepath="{0}/{1}_chromosome_simu3_transdeg.csv".format(filedir, basename)
#	print("Processing {0}/{1}".format(_+1, tmp_df.shape[0]))
#	df=pd.read_csv(filepath,dtype={'idx':str})
#	df=df.set_index('idx') 
#	
#	total=np.sum(np.sum(df))
#	for idx in df.index:
#		for clm in df.columns:
#			dif=get_num_of_nsrf(idx)-get_num_of_nsrf(clm)
#			df.loc[idx,clm]=df.loc[idx,clm]*dif/total
#
#	sym_df=pd.DataFrame(0,index=df.index, columns=df.columns)
#	for i in range(df.shape[0]):
#		for j in range(i):
#			sym_df.iloc[i,j]=df.iloc[i,j]+df.iloc[j, i]
#	
#	dct={}
#
#	for i in range(df.shape[0]):
#		for j in range(df.shape[1]):
#			name="{0}->{1}".format(i,j)
#			dct[name]=sym_df.iloc[i][j]
#	dct_lst.append(dct)

columns_lst=["ftp_basename","G+C","taxid","count_diff"]
for i in range(df.shape[0]):
	for j in range(df.shape[1]):
		columns_lst.append("{}->{}".format(i,j))
out_df=pd.DataFrame(dct_lst)
out_df=out_df[columns_lst]
out_df.to_csv("out.csv", index=False)
