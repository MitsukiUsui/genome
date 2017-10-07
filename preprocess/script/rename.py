import sys
import subprocess

sys.path.append("/home/mitsuki/altorf/genome/helper")
from dbcontroller import DbController
from myutil import myrun

"""
!!! UNTESTED with myrun !!!
"""

def main():
    directory_lst = ["/data/mitsuki/data/refseq/genomic_fna",
                     "/data/mitsuki/data/refseq/genomic_gff",
                     "/data/mitsuki/data/refseq/cds_from_genomic"]
    suffix_lst = ["_genomic.fna",
                  "_genomic.gff",
                  "_cds_from_genomic.fna"]

    dbFilepath = "../db/altorf.db"
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("preprocess")

    for _, taxid in enumerate(taxid_lst):
        row = dc.get_row(taxid)

        for directory, suffix in zip(directory_lst, suffix_lst):
            realFilepath = directory + "/" + row["ftp_basename"] + suffix
            gzFilepath = realFilepath + ".gz"
            linkFilepath = directory + "/" + str(taxid) + suffix

            # gunzip
            cmd = "gunzip {}".format(gzFilepath)
            success = myrun(cmd)
            if not(success):
                print("ERROR: {}".format(cmd))

            # create symbolic link
            cmd = "ln -s {} {}".format(realFilepath, linkFilepath)
            success = myrun(cmd)
            if not(success):
                print("ERROR: {}".format(cmd))


if __name__ == "__main__":
    main()
