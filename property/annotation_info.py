from BCBio import GFF
from Bio import SeqIO
import pandas as pd
import numpy as np

def main(lookupFilepath):
    lookup_df=pd.read_csv(lookupFilepath)
    
    dct_lst=[]
    for basename in lookup_df["ftp_basename"]:
        print("PROCESSING {}".format(basename))
        
        seqFilepath="/data/mitsuki/out/altorf/genome/fasta/{}_chromosome.fna".format(basename)
        gffFilepath="/data/mitsuki/data/refseq/genomic_gff/{}_genomic.gff".format(basename)

        # load fasta
        seqRec_lst=[]
        seqName_lst=[]
        for seqRec in SeqIO.parse(seqFilepath, "fasta"):
            seqRec_lst.append(seqRec)
            seqName_lst.append(seqRec.id)

        # load gff and distribute CDS
        cds_lstlst=[[] for _ in range(len(seqName_lst))]
        with open(gffFilepath) as f:
            for rec in GFF.parse(f,target_lines=1):
                idx=-1
                for i, seqName in enumerate(seqName_lst):
                    if rec.id==seqName:
                        idx=i
                        break
                if idx>=0:
                    assert len(rec.features)==1
                    if rec.features[0].type=="CDS":
                        cds_lstlst[idx].append(rec)
        
        # colect information
        for seqRec, cds_lst, seqName in zip(seqRec_lst, cds_lstlst, seqName_lst):
            dct={}
            dct["ftp_basename"]=basename
            dct["seq_name"]=seqName
            dct["seq_length"]=len(seqRec)
            dct["num_cds"]=len(cds_lst)

            countTypical=0
            cdsLength=0
            typicalLength=0
            for cds in cds_lst:
                start=int(cds.features[0].location.start)
                end=int(cds.features[0].location.end)
                length=end-start
                cdsLength+=length
                if length%3==0 and length>=6:
                    countTypical+=1
                    typicalLength+=length
            dct["num_typical"]=countTypical
            dct["cds_length"]=cdsLength
            dct["typical_length"]=typicalLength
            dct_lst.append(dct)


if __name__=="__main__":
    lookupFilepath="../speciespick/picked_assembly_summary_code.csv"
    main(lookupFilepath)