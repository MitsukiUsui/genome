import subprocess
import sys

sys.path.append("/home/mitsuki/altorf/genome/db")
from dbcontroller import DbController


def wget_helper(outFilepath, ftpFilepath):
    cmd = "wget -q -O {} {}".format(outFilepath, ftpFilepath)
    try:
        subprocess.run(cmd.split(), check=True)
    except subprocess.CalledProcessError as e:
        print(e)
        return False
    return True


def main(dbFilepath):
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("download")
    print("START: download {} species from RefSeq".format(len(taxid_lst)))

    prefix_lst = ["_genomic.fna.gz", "_genomic.gff.gz", "_cds_from_genomic.fna.gz"]
    for _, taxid in enumerate(taxid_lst):
        print("START: taxid={} ({}/{})".format(taxid, _ + 1, len(taxid_lst)))

        successCount = 0
        row = dc.get_row(taxid)

        # genomic.fna
        filename = row["ftp_basename"] + "_genomic.fna.gz"
        outFilepath = "/data/mitsuki/data/refseq/genomic_fna/" + filename
        ftpFilepath = row["ftp_path"] + "/" + filename
        successCount += wget_helper(outFilepath, ftpFilepath)

        # genomic.gff
        filename = row["ftp_basename"] + "_genomic.gff.gz"
        outFilepath = "/data/mitsuki/data/refseq/genomic_gff/" + filename
        ftpFilepath = row["ftp_path"] + "/" + filename
        successCount += wget_helper(outFilepath, ftpFilepath)

        # genomic.gff
        filename = row["ftp_basename"] + "_cds_from_genomic.fna.gz"
        outFilepath = "/data/mitsuki/data/refseq/cds_from_genomic/" + filename
        ftpFilepath = row["ftp_path"] + "/" + filename
        successCount += wget_helper(outFilepath, ftpFilepath)

        if successCount == 3:
            dc.mark_as_done(taxid, "download")
        else:
            print("ERROR: download only {} files".format(successCount))


if __name__ == "__main__":
    dbFilepath = "../db/altorf.db"
    main(dbFilepath)
