#!/home/mitsuki/.pyenv/versions/anaconda3-4.2.0/bin/python

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO

def get_orf_df(seqRec):

    dct_lst=[]
    for lane in range(6):
        gap=lane%3
        start=gap
        pos=start+3
        print("\t\tlane: {}, gap: {}, start: {}, pos: {}".format(lane, gap, start, pos))

        while pos<=len(seqRec):
            codon=seqRec.seq[pos-3:pos] 
            if lane<3:
                if str(codon).upper() in ["TAA","TAG", "TGA"]:
                    dct={"lane":lane, "start":start, "end":pos}
                    dct_lst.append(dct)
                    start=pos
            else:
                if str(codon).upper() in ["TTA","CTA", "TCA"]:
                    dct={"lane":lane, "start":start, "end":pos-3}
                    dct_lst.append(dct)
                    start=pos-3
            pos+=3

        dct={"lane":lane, "start":start, "end":pos-3}  #最後のorfの処理
        dct_lst.append(dct)
        
    df=pd.DataFrame(dct_lst)
    df["length"]=df["end"]-df["start"]
    return df

def count_overlap(df, length):
    pos_lst=[]
    for _,row in df.iterrows():
        pos_lst.append((row["start"], row["lane"]))
        pos_lst.append((row["end"], row["lane"]))
            
    pos_lst=sorted(pos_lst, key=lambda x: x[0])
    
    overlapCounter=np.zeros(7, dtype=int)
    isOpen=np.zeros(6, dtype=bool)
    lastPos=0
    for pos, lane in pos_lst:
        overlapCounter[np.sum(isOpen)]+=pos-lastPos
        isOpen[lane]=~(isOpen[lane])
        lastPos=pos
    overlapCounter[np.sum(isOpen)]+=length-lastPos
    return overlapCounter
    
def main(seqFilepath, outFilepath):
    thres_lst=list(range(50, 1000+1, 50))
    out_mat=np.zeros((len(thres_lst), 7)).astype(int)
    
    seqRec_lst=[]
    for seqRec in SeqIO.parse(seqFilepath, "fasta"):
        seqRec_lst.append(seqRec)
    for seqRec in seqRec_lst:
        print("PROCESSING {}".format(seqRec.id))

        orf_df=get_orf_df(seqRec)
        print("\tfound {} orfs".format(orf_df.shape[0]))

        for i, thres in enumerate(thres_lst):
            filtered_df=orf_df[orf_df["length"]>=thres]
            overlapCounter=count_overlap(filtered_df, len(seqRec))
            out_mat[i]+=overlapCounter
            aveOverlap=np.dot(overlapCounter, np.arange(7))/len(seqRec)
            print("\tthres = {}, average overlap = {}, {} orfs".format(thres, aveOverlap, filtered_df.shape[0]))

    out_df=pd.DataFrame(out_mat)
    out_df["thres"]=thres_lst
    out_df=out_df[["thres", 0, 1, 2, 3, 4, 5, 6]]
    out_df.to_csv(outFilepath, index=False)

if __name__=="__main__":
    seqFilepath=sys.argv[1]
    outFilepath=sys.argv[2]
    main(seqFilepath, outFilepath)
