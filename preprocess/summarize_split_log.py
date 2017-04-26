import pandas as pd

catalogDf=pd.read_csv("/home/mitsuki/altorf/genome/speciespick/picked_assembly_summary.csv")

dct_lst=[]
basedir="/home/mitsuki/out/altorf/genome/preprocess/split/"
for basename in catalogDf["ftp_basename"]:
    filepath=basedir+basename+"_split.log"
    columns=["ftp_basename","dna_type","filepath","num_of_seqs","length_of_seqs"]
    with open(filepath,'r') as f:
        for line in f:
            dct={}
            line_splitted=line.strip().split(',')
            assert len(line_splitted)==len(columns)
            for i,column in enumerate(columns):
                dct[column]=line_splitted[i]
            dct_lst.append(dct)
    print("\tDONE processing "+basename)


summaryDf=pd.DataFrame(dct_lst)
summaryDf.to_csv("summary_split.csv",index=False)
print("DONE summarizing to summary_split.csv")