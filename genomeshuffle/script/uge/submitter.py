#!/usr/bin/env python3

import sys

sys.path.append("/home/mitsuki/altorf/genome/helper")
from dbcontroller import DbController
from myutil import myrun

def main(dbFilepath):
    dc = DbController(dbFilepath)
#    taxid_lst = dc.get_target("shuffle")
    with open("rerun.list", "r") as f:
        taxid_lst=[line.strip() for line in f]

    print("START: submit {} jobs".format(len(taxid_lst)))
    for taxid in taxid_lst:
        cmd = "qsub caller.sh {}".format(taxid)
#        print(cmd)
        myrun(cmd)

if __name__ == "__main__":
    dbFilepath = "/home/mitsuki/altorf/genome/db/altorf.db"
    main(dbFilepath)
