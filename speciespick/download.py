import sys

sys.path.append("/home/mitsuki/altorf/genome/helper/")
from dbcontroller import DbController
from myutil import myrun

"""
!!!UNTESTED with myrun function!!!
"""


def main(dbFilepath):
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("download")
    print("START: download {} species from RefSeq".format(len(taxid_lst)))

    for _, taxid in enumerate(taxid_lst):
        print("START: taxid={} ({}/{})".format(taxid, _ + 1, len(taxid_lst)))
        row = dc.get_row(taxid)

        directory_lst = ["/data/mitsuki/data/refseq/genomic_fna",
                         "/data/mitsuki/data/refseq/genomic_gff",
                         "/data/mitsuki/data/refseq/cds_from_genomic"]
        suffix_lst = ["_genomic.fna.gz",
                      "_genomic.gff.gz",
                      "_cds_from_genomic.fna.gz"]

        for directory, suffix in zip(directory_lst, suffix_lst):
            filename = row["ftp_basename"] + suffix
            outFilepath = directory + "/" + filename
            ftpFilepath = row["ftp_path"] + "/" + filename
            cmd = "wget -q -O {} {}".format(outFilepath, ftpFilepath)
            success = myrun(cmd)
            if (not(success)):
                break

        if success:
            dc.mark_as_done([taxid], "download")
        else:
            print("ERROR: {}".format(cmd))


if __name__ == "__main__":
    dbFilepath = "../db/altorf.db"
    main(dbFilepath)
