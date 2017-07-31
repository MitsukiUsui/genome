import pandas as pd
import numpy as np

def main(lookupFilepath, outFilepath):
	lookup_df=pd.read_csv(lookupFilepath)
	
	dct_lst=[]
	for _, row in lookup_df.iterrows():
		print(row["ftp_basename"])

		direc="/data/mitsuki/out/altorf/genome/patternanalyze/overlap"
		filepath_lst=[]
		filepath_lst.append(direc+"/{}_chromosome_overlap.csv".format(row["ftp_basename"]))
		filepath_lst.append(direc+"/{}_chromosome_sim_overlap.csv".format(row["ftp_basename"]))
		ident_lst=['r','s']

		dct={}
		dct["taxid"]=row["taxid"]

		for i in range(2):
			filepath=filepath_lst[i]
			ident=ident_lst[i]

			df=pd.read_csv(filepath, index_col=0)
			for thres in range(50, 1000+1, 50):
				aveOverlap=np.dot(df[str(thres)], np.arange(7))/np.sum(df[str(thres)])
				dct[ident+str(thres)]=aveOverlap
		dct_lst.append(dct)
	out_df=pd.DataFrame(dct_lst)
	column_lst=["taxid"]
	for ident in ['r','s']:
		for thres in range(50,1000+1,50):
			column_lst.append(ident+str(thres))
	out_df=out_df[column_lst]
	out_df.to_csv(outFilepath, index=False)

if __name__=="__main__":
	lookupFilepath="../speciespick/picked_assembly_summary_code.csv"
	outFilepath="summary_overlap.csv"
	main(lookupFilepath, outFilepath)
