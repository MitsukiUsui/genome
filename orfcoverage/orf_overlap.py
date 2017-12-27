#!/home/mitsuki/.pyenv/versions/anaconda3-4.2.0/bin/python

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO

def create_orf_df(seq):
    """
    :param seq:
    :return: dataframe of all stop to stop orfs
    """

    dct_lst=[]
    for lane in range(6):
        start = lane % 3
        pos = start+3

        while pos <= len(seq):
            codon = seq[pos-3 : pos]
            if lane in (0, 1, 2): #forward lane
                if str(codon).upper() in ["TAA", "TAG", "TGA"]:
                    dct={"lane":lane, "start":start, "end":pos}
                    dct_lst.append(dct)
                    start=pos
            elif lane in (3, 4, 5): #backward lane
                if str(codon).upper() in ["TTA", "CTA", "TCA"]:
                    dct={"lane":lane, "start":start, "end":pos-3}
                    dct_lst.append(dct)
                    start=pos-3
            pos+=3

        dct={"lane":lane, "start":start, "end":pos-3} #process last orf
        dct_lst.append(dct)
        
    df=pd.DataFrame(dct_lst)
    df["length"]=df["end"]-df["start"]
    df = df[["lane", "start", "end", "length"]]
    return df

def count_overlap(orf_df, seqLength):
    """
    :param orf_df: should contains lane, start, end
    :param seqLength: length of seq. needed for calc last overlap length
    :return: overlapCount_arr (bp)
    """

    stateSwitch_lst=[] #switch for state (open/close)
    for _,row in orf_df.iterrows():
        stateSwitch_lst.append((row["start"], row["lane"]))
        stateSwitch_lst.append((row["end"], row["lane"]))
    stateSwitch_lst=sorted(stateSwitch_lst, key=lambda x: x[0])
    
    overlapCount_arr=np.zeros(7, dtype=int)
    isOpenLane=np.zeros(6, dtype=bool)
    lastPos=0
    for pos, lane in stateSwitch_lst:
        overlapCount_arr[np.sum(isOpenLane)] += pos - lastPos
        isOpenLane[lane] = ~(isOpenLane[lane])
        lastPos = pos
    overlapCount_arr[np.sum(isOpenLane)] += seqLength - lastPos
    return overlapCount_arr

def main(seqFilepath, outFilepath):
    seqRec_lst=[]
    for seqRec in SeqIO.parse(seqFilepath, "fasta"):
        seqRec_lst.append(seqRec)

    dct_lst = []
    thres_lst=list(range(50, 1000+1, 50)) #threshold for orf length

    for seqRec in seqRec_lst:
        print("PROCESSING {}".format(seqRec.id))

        orf_df=create_orf_df(seqRec.seq)
        print("\tfound {} orfs".format(orf_df.shape[0]))

        for thres in thres_lst:
            filtered_df = orf_df[orf_df["length"]>=thres]
            overlapCount_arr = count_overlap(filtered_df, len(seqRec))

            dct = {}
            dct["thres"] = thres
            dct["id"] = seqRec.id
            for i in range(7):
                dct[i]=overlapCount_arr[i]
            dct_lst.append(dct)

            orfCoverage=np.dot(overlapCount_arr, np.arange(7))/len(seqRec)
            print("\tthres = {}, orfCoverage = {}, {} orfs".format(thres, orfCoverage, filtered_df.shape[0]))

    out_df = pd.DataFrame(dct_lst)
    out_df.to_csv(outFilepath, index=False)

if __name__=="__main__":
    seqFilepath=sys.argv[1]
    outFilepath=sys.argv[2]
    main(seqFilepath, outFilepath)
