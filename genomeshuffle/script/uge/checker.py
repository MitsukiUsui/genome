#!/usr/bin/env python3

import sys
import os

sys.path.append("/home/mitsuki/altorf/genome/helper")
from dbcontroller import DbController

def main(dbFilepath):
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("shuffle")
    for taxid in taxid_lst:
        simFilepath="/data/mitsuki/out/altorf/genome/fasta/{}_sim.fna".format(taxid)
        bedFilepath="/data/mitsuki/out/altorf/genome/genomeshuffle/bed/{}.bed".format(taxid)
        if os.path.isfile(simFilepath) and os.path.isfile(bedFilepath):
            pass
        else:
#            print("ERROR: {} not finished".format(taxid))
            print(taxid)
    #dc.mark_as_done(taxid_lst, "shuffle")
    #print("DONE: mark {} taxids as done".format(len(taxid_lst)))


if __name__ == "__main__":
    dbFilepath = "/home/mitsuki/altorf/genome/db/altorf.db"
    main(dbFilepath)
