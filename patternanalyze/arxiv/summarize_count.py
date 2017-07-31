import pandas as pd
import numpy as np
import sys


def main(outFilepath):
    lookupFilepath="../speciespick/picked_assembly_summary_code.csv"
    lookup_df=pd.read_csv(lookupFilepath)

    real_lst=[]
    sim_lst=[]
    for key,row in lookup_df.iterrows():
        countDir="/home/mitsuki/out/altorf/genome/patternanalyze/count"
        countFilepath="{0}/{1}_{2}_{3}_count.csv".format(countDir,row["ftp_basename"],"chromosome","simu3")
        count_df=pd.read_csv(countFilepath,index_col=0)

        assert count_df["real"].sum()==count_df["simu"].sum()
        totalWindow=count_df["real"].sum()
        real_lst.append(np.dot(np.array(count_df["real"]),np.arange(7))/totalWindow)
        sim_lst.append(np.dot(np.array(count_df["simu"]),np.arange(7))/totalWindow)


    out_df=pd.DataFrame()
    out_df["taxid"]=lookup_df["taxid"]
    out_df["count_real"]=real_lst
    out_df["count_sim"]=sim_lst
    out_df.to_csv(outFilepath,index=False)
    
if __name__=="__main__":
    main(sys.argv[1])
