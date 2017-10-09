import sys

sys.path.append("/home/mitsuki/altorf/genome/helper")
from dbcontroller import DbController
from myutil import myrun


def main(dbFilepath):
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("shuffle")
    print("START: submit {} jobs".format(len(taxid_lst)))
    for taxid in taxid_lst:
        directory = "/home/mitsuki/altorf/genome/genomeshuffle/script/uge"

        que = "standard.q"
        out = "{}/log/{}.out".format(directory, taxid)
        err = "{}/log/{}.err".format(directory, taxid)
        cmd = "{}/caller.sh".format(directory)

        ugecmd = "qsub -S /bin/bash -q {0} -o {1} -e {2} {3} {4}".format(que, out, err, cmd, taxid)
        myrun(ugecmd)


if __name__ == "__main__":
    dbFilepath = "/home/mitsuki/altorf/genome/db/altorf.db"
    main(dbFilepath)
