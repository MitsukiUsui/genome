import sys

sys.path.append("/home/mitsuki/altorf/genome/helper")
from dbcontroller import DbController
from myutil import myrun


def main(dbFilepath):
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("preprocess")
    print("START: submit {} jobs".format(len(taxid_lst)))
    for taxid in taxid_lst:
        que = "s8.q"
        out = "/home/mitsuki/altorf/genome/preprocess/script/uge/log/{}.out".format(taxid)
        err = "/home/mitsuki/altorf/genome/preprocess/script/uge/log/{}.err".format(taxid)
        cmd = "/home/mitsuki/altorf/genome/preprocess/script/uge/caller.sh"

        ugecmd = "qsub -S /bin/bash -q {0} -o {1} -e {2} {3} {4}".format(que, out, err, cmd, taxid)
        myrun(ugecmd)


if __name__ == "__main__":
    dbFilepath = "/home/mitsuki/altorf/genome/db/altorf.db"
    main(dbFilepath)
