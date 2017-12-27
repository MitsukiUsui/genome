#!/home/mitsuki/.pyenv/versions/anaconda3-4.2.0/bin/python

import sys
from BCBio import GFF
from Bio import SeqIO
import pandas as pd

from orfcoverage.orf_overlap import create_orf_df
from orfcoverage.phase import PhaseController

class Interval:
    def __init__(self, intervalId, start, end, phase):
        assert phase in PhaseController.phase_lst
        assert (end - start) % 3 == 0

        self.intervalId = intervalId
        self.start = start
        self.end = end
        self.phase = phase

def create_orf_interval_lst(orf_df, thres):
    filtered_df = orf_df[orf_df["length"] >= thres]
    orfInterval_lst = []
    for orfId, row in filtered_df.iterrows():
        phase = PhaseController.phase_lst[row["lane"]]
        interval = Interval(orfId, row["start"], row["end"], phase)
        orfInterval_lst.append(interval)
        
def create_cds_interval_lst(gffFilepath, seqId):
    cdsId = 0
    cdsInterval_lst = []
    with open(gffFilepath) as f:
        for rec in GFF.parse(f, target_lines=1):
            assert len(rec.features)==1
            feature = rec.fetures[0]
            if (rec.id == seqId) and (feature.type == "CDS"): #if this is targeted cds
                start = feature.location.start
                end = feature.location.end
                length = end - start
                if (length >= 6) and (length % 3 == 0): # if this is targeted typical cds
                    lane = (start % 3) if (feature.location.strand == 1) else (3 + (start % 3))
                    phase = PhaseController.phase_lst[lane]
                    interval = Interval(cdsId, start, end, phase)
                    cdsInterval_lst.append(interval)
                    cdsId += 1

def count_relative_overlap(cdsInterval_lst, orfInterval_lst):
    """
    :param orfInterval_lst:
    :param cdsInterval_lst:
    :return: relativeOverlap_dct: keys(phase_lst+"intergenic"), val(count in bp)
    """

    class Switch:
        def __init__(self, interval, isStart, type):
            self.intervalId = interval.intervalId
            self.isStart = isStart
            self.position = (interval.start) if (isStart) else (interval.end)
            self.phase = interval.phase
            self.type = type

    class Overlap:
        def __init__(self, start, end, relativePhase):
            assert relativePhase in (PhaseController.phase_lst + ["intergenic"])

            self.start = start
            self.end = end
            self.relativePhase = relativePhase

    #initialize switch
    switch_lst=[]
    for interval_lst, type in zip((cdsInterval_lst, orfInterval_lst), ("cds", "orf")):
        for interval in interval_lst:
            startSwitch = Switch(interval, isStart = True, type = type)
            endSwitch = Switch(interval, isStart = False, type = type)
            switch_lst.append(startSwitch)
            switch_lst.append(endSwitch)
    switch_lst = sorted(switch_lst, key = lambda x: x.position)

    inProcessCds_dct = {} #(key: id, value: phase)
    inProcessOrf_dct = {} #(key: id, value: phase)
    overlap_dct = {} #(key: (cdsId, orfId), value: overlap)
    phaseController = PhaseController()

    for switch in switch_lst:
        if switch.type == "cds":
            cdsId = switch.intervalId
            cdsPhase = switch.phase
            if switch.isStart:
                inProcessCds_dct[cdsId] = cdsPhase
                for orfId, orfPhase in inProcessOrf_dct.items():
                    key = (cdsId, orfId)
                    relativePhase = phaseController.relative(cdsPhase, orfPhase)
                    overlap = Overlap(start = switch.position, end = switch.position, relativePhase = relativePhase)
                    overlap_dct[key] = overlap
            else:
                del inProcessCds_dct[cdsId]
                for orfId, _ in inProcessOrf_dct.items():
                    key = (cdsId, orfId)
                    overlap_dct[key].end = switch.position
        elif switch.type == "orf":
            orfId = switch.intervalId
            orfPhase = switch.phase
            if switch.isStart:
                inProcessOrf_dct[orfId] = orfPhase
                for cdsId, cdsPhase in inProcessCds_dct.items():
                    key = (cdsId, orfId)
                    relativePhase = phaseController.relative(cdsPhase, orfPhase)
                    overlap = Overlap(start = switch.position, end = switch.position, relativePhase = relativePhase)
                    overlap_dct[key] = overlap
            else:
                del inProcessOrf_dct[orfId]
                for cdsId, _ in inProcessCds_dct.items():
                    key = (cdsId, orfId)
                    overlap_dct[key].end = switch.position

    #return overlapCounter_dct
    overlapCount_dct = {"intergenic" : 0}
    for phase in PhaseController.phase_lst:
        overlapCount_dct[phase] = 0
    for _, overlap in overlap_dct.items():
        overlapCount_dct[overlap.relativePhase] = overlap.end - overlap.start
    return overlapCount_dct


def main(seqFilepath, gffFilepath, outFilepath):
    seqRec_lst=[]
    for seqRec in SeqIO.parse(seqFilepath, "fasta"):
        seqRec_lst.append(seqRec)
    print("DONE: load {} seqs from {}".format(len(seqRec_lst), seqFilepath))

    dct_lst = []
    thres_lst=list(range(50, 1000+1, 50))
    for seqRec in seqRec_lst:
        print("START: process {}".format(seqRec.id))
        orf_df = create_orf_df(seqRec.seq)
        cdsInterval_lst = create_cds_interval_lst(gffFilepath, seqRec.id)

        for thres in thres_lst:
            orfInterval_lst = create_orf_interval_lst(orf_df, thres)
            overlapCount_dct = count_relative_overlap(cdsInterval_lst, orfInterval_lst)
            overlapCount_dct["thres"] = thres
            overlapCount_dct["id"] = seqRec.id
            dct_lst.append(overlapCount_dct)
            print("\nDONE: process with thres = {}".format(thres))

    out_df = pd.DataFrame(dct_lst)
    out_df.to_csv(outFilepath, index=False)

if __name__=="__main__":
    seqFilepath=sys.argv[1]
    gffFilepath=sys.argv[2]
    outFilepath=sys.argv[3] 
    main(seqFilepath, gffFilepath, outFilepath)
