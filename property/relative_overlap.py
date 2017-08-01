#!/home/mitsuki/.pyenv/versions/anaconda3-4.2.0/bin/python

import sys
from BCBio import GFF
from Bio import SeqIO
import Bio
import pandas as pd
import numpy as np

from orf_overlap import get_orf_df

def get_pos_lst(cds_lst, orf_df):
    """
    given cds_lst and orf_df, return sorted pos_lst
    format: (pos, isStart, type, name, lane)
        type:0 for cds, 1 for orf
    """
    pos_lst=[] 

    #update pos_lst with cds_lst
    for cdsId,cds in enumerate(cds_lst):
        assert len(cds.features)==1
        start=int(cds.features[0].location.start)
        end=int(cds.features[0].location.end)
        lane=start%3
        if cds.features[0].location.strand==-1:
            lane+=3
                    
        length=end-start
        if length%3==0 and length>=6:  #if "typical cds"
            pos_lst.append((start, True, 0, cdsId, lane))
            pos_lst.append((end, False, 0, cdsId, lane))

    #update pos_lst with df
    for orfId,row in orf_df.iterrows():
        pos_lst.append((row["start"], True, 1, orfId, row["lane"]))
        pos_lst.append((row["end"], False, 1, orfId, row["lane"]))
                
    pos_lst=sorted(pos_lst, key=lambda x: x[0])
    return pos_lst
    

def get_overlap_dctdct(pos_lst):
    """
    given sorted pos_lst, return overlap_dctdct
    format: key = (orfId, cdsId), value = {ostart, oend, relLane}
    
    """

    def get_relative_lane(cdsLane, orfLane):
        dct={}
        dct[0]=["+1", "+2", "+3", "-1", "-2", "-3"]
        dct[1]=["+3", "+1", "+2", "-3", "-1", "-2"]
        dct[2]=["+2", "+3", "+1", "-2", "-3", "-1"]
        dct[3]=["-1", "-3", "-2", "+1", "+3", "+2"]
        dct[4]=["-2", "-1", "-3", "+2", "+1", "+3"]
        dct[5]=["-3", "-2", "-1", "+3", "+2", "+1"]
        return dct[cdsLane][orfLane]
    
    overlap_dctdct={}
    proCds_dct= {} #dictionary of cdss in process (key: cdsId, value: naive lane)
    proOrf_dct={}  #dictionary of orfs in process (key: orfId, value: naive lane)
    for pos in pos_lst:
        if pos[2]==0:#if CDS
            cdsId=pos[3]
            cdsLane=pos[4]
            if pos[1]:#if start
                proCds_dct[cdsId]=cdsLane  #update processing CDS dictionary
                for orfId, orfLane in proOrf_dct.items():  #add new overlap object with orfs in proOrf_dct
                    key=(orfId, cdsId)
                    overlap_dctdct[key]={}
                    overlap_dctdct[key]["ostart"]=pos[0]
                    overlap_dctdct[key]["relLane"]=get_relative_lane(cdsLane, orfLane)                   
            else:#if last
                del proCds_dct[cdsId]  #delete from processing CDS dicrionary
                for orfId, _ in proOrf_dct.items():
                    key=(orfId, cdsId)
                    overlap_dctdct[key]["oend"]=pos[0]

        elif pos[2]==1:  #if orf
            orfId=pos[3]
            orfLane=pos[4]
            if pos[1]:  #if start
                proOrf_dct[orfId]=orfLane
                for cdsId, cdsLane in proCds_dct.items():
                    key=(orfId, cdsId)
                    overlap_dctdct[key]={}
                    overlap_dctdct[key]["ostart"]=pos[0]
                    overlap_dctdct[key]["relLane"]=get_relative_lane(cdsLane, orfLane)
            else:  #if last
                del proOrf_dct[orfId]  #delete from processing orf dicrionary
                for cdsId, _ in proCds_dct.items():
                    key=(orfId, cdsId)
                    overlap_dctdct[key]["oend"]=pos[0]
    return overlap_dctdct
    
    
    
def main(seqFilepath, gffFilepath, outFilepath):
    # load fasta
    seqRec_lst=[]
    seqName_lst=[]
    for seqRec in SeqIO.parse(seqFilepath, "fasta"):
        seqRec_lst.append(seqRec)
        seqName_lst.append(seqRec.id)
    print("LOADED {} seqs from {}".format(len(seqRec_lst), seqFilepath))

    # load gff and distribute CDS
    cds_lstlst=[[] for _ in range(len(seqName_lst))]
    with open(gffFilepath) as f:
        for rec in GFF.parse(f,target_lines=1):
            assert len(rec.features)==1
            if rec.features[0].type=="CDS":
                try:
                    idx=seqName_lst.index(rec.id)
                    cds_lstlst[idx].append(rec)
                except ValueError:
                    pass # corresponding sequence does not exists in seqRec_lst

    
    for idx, seqName in enumerate(seqName_lst):
        print("\tLOADED {0} CDSs in {1}".format(len(cds_lstlst[idx]), seqName))

    
    thres_lst=list(range(50, 1000+1, 50))
    columns=["+1", "+2", "+3", "-1", "-2", "-3"]
    out_mat=np.zeros((len(thres_lst), len(columns))).astype(int)

    for seqRec, cds_lst, seqName in zip(seqRec_lst, cds_lstlst, seqName_lst):
        orf_df=get_orf_df(seqRec)
        
        for i, thres in enumerate(thres_lst):
            filtered_df=orf_df[orf_df["length"]>=thres]
            pos_lst=get_pos_lst(cds_lst, filtered_df)
            overlap_dctdct=get_overlap_dctdct(pos_lst)
            
            for _, dct in overlap_dctdct.items():
                out_mat[i, columns.index(dct["relLane"])] += (dct["oend"]-dct["ostart"])
        print("\tDONE with {}".format(seqName))
            

    out_df=pd.DataFrame(out_mat, columns=columns)
    out_df["thres"]=thres_lst
    out_df=out_df[["thres"]+columns]
    out_df.to_csv(outFilepath, index=False)
    print("OUTPUT to {}".format(outFilepath))


if __name__=="__main__":
    seqFilepath=sys.argv[1]
    gffFilepath=sys.argv[2]
    outFilepath=sys.argv[3] 
    main(seqFilepath, gffFilepath, outFilepath)
