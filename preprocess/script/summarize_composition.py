import pandas as pd

catalogDf=pd.read_csv("summary_split.csv")
catalogDf.head()

#read each output of fatt composition and summarize into out_df
basedir="/data/mitsuki/out/altorf/genome/preprocess/comp/"
dct_lst=[]
for idx,record in catalogDf[catalogDf["num_of_seqs"]>0].iterrows():#exclude empty files
    dct={}
    dct["ftp_basename"]=record["ftp_basename"]
    dct["dna_type"]=record["dna_type"]
    
    filepath=basedir+record["ftp_basename"]+"_"+record["dna_type"]+".comp"
    with open(filepath,'r') as f:
        f.readline()#skip Total # n-mers header
        f.readline()#skip length

        for line in f:
            line=line.strip()
            if line in ["1-mer stats","2-mer stats","3-mer stats"]:
                continue
            else:
                line_splitted=line.split('\t')
                dct[line_splitted[0]]=line_splitted[1]
    dct_lst.append(dct)
out_df=pd.DataFrame(dct_lst)



# define columns to output(columns_sorted)
basicBase_tpl=('A','C','G','T')
columns_sorted=list(basicBase_tpl)
#add other single base
ambiguousBase_lst=[]
for column in out_df.columns:
    if(len(column)==1) and not(column in columns_sorted):
        ambiguousBase_lst.append(column)
ambiguousBase_lst.sort()
columns_sorted+=ambiguousBase_lst
#add 2-mer
for i in basicBase_tpl:
    for j in basicBase_tpl:
        columns_sorted.append(i+j)
#add 3-mer
for i in basicBase_tpl:
    for j in basicBase_tpl:
        for k in basicBase_tpl:
            columns_sorted.append(i+j+k)
            

#calculate important measures according to the output of fatt composition
out_df["num_basic"]=out_df[list(basicBase_tpl)].astype(float).sum(axis=1)
out_df["num_ambiguous"]=(out_df[ambiguousBase_lst].astype(float).sum(axis=1)).astype(int)
out_df["per_ambiguous"]=out_df["num_ambiguous"]*100/(out_df["num_basic"]+out_df["num_ambiguous"])
out_df["G+C"]=(out_df[['C','G']].astype(float).sum(axis=1))/out_df["num_basic"]

#output designated measures into summary_composition.csv
out_df=out_df[["ftp_basename","dna_type","G+C","num_ambiguous","per_ambiguous"]+columns_sorted]
out_df.to_csv("summary_composition.csv",index=False)
