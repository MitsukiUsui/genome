import sys

sys.path.append("/home/mitsuki/altorf/genome/helper")
from dbcontroller import DbController

def main(dbFilepath):
    dc = DbController(dbFilepath)
    taxid_lst = dc.get_target("preprocess")
    dc.mark_as_done(taxid_lst, "preprocess")
    print("DONE: mark {} taxids as done".format(len(taxid_lst)))


if __name__ == "__main__":
    dbFilepath = "/home/mitsuki/altorf/genome/db/altorf.db"
    main(dbFilepath)
