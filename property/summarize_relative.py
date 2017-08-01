import pandas as pd

def main(lookupFilepath, outFilpeath):
    lookup_df=pd.read_csv(lookupFilepath)
    out_df=pd.DataFrame(columns=["ftp_basename","seq_type","thres","+1","+2","+3","-1","-2","-3"])
    
    for basename in lookup_df["ftp_basename"]:
        print(basename)
        direc="/data/mitsuki/out/altorf/genome/property/relative"
        realFilepath=direc+"/{}_chromosome_relative.csv".format(basename)
        simFilepath=direc+"/{}_chromosome_sim_relative.csv".format(basename)
        
        real_df=pd.read_csv(realFilepath)
        real_df["ftp_basename"]=basename
        real_df["seq_type"]="real"
        out_df=pd.concat([out_df, real_df])
        
        sim_df=pd.read_csv(simFilepath)
        sim_df["ftp_basename"]=basename
        sim_df["seq_type"]="sim"
        out_df=pd.concat([out_df, sim_df])
   
    int_lst=["thres","+1","+2","+3","-1","-2","-3"]
    out_df[int_lst]=out_df[int_lst].astype(int)
    out_df=out_df[["ftp_basename", "seq_type"]+int_lst]
    out_df.to_csv(outFilepath, index=False)
    print("OUTPUT to {}".format(outFilepath))
        

if __name__=="__main__":
    lookupFilepath="../speciespick/picked_assembly_summary_code.csv"
    outFilepath="summary_relative.csv"
    main(lookupFilepath, outFilepath)
