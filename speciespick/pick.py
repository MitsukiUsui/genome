import pandas as pd
import numpy as np

import sys
sys.path.append("/home/mitsuki/altorf/genome/db")
from dbcontroller import DbController


def get_df(filepath_lst, domain_lst):
    assert len(filepath_lst) == len(domain_lst)

    dct_lst = []
    for filepath, domain in zip(filepath_lst, domain_lst):
        df = pd.read_csv(filepath, delimiter='\t', skiprows=1, dtype="object")
        print("DONE : load {} species from {}".format(df.shape[0], filepath))

        filtered_df = df[(df["refseq_category"] == "representative genome") & (df["assembly_level"] == "Complete Genome")]
        for _, row in filtered_df.iterrows():
            dct = {}
            column_lst = ["taxid", "organism_name", "ftp_path"]
            for column in column_lst:
                dct[column] = row[column]
            dct["ftp_basename"] = row["ftp_path"].split('/')[-1]
            dct["domain"] = domain
            dct_lst.append(dct)
        print("DONE : pick {} species".format(filtered_df.shape[0]))

    return pd.DataFrame(dct_lst)


def add_genetic_code(df, taxFilepath):
    names_lst = [str(i) for i in range(14)]
    names_lst[0] = "taxid"
    names_lst[6] = "genetic_code"
    tax_df = pd.read_csv(taxFilepath, names=names_lst)

    tax_df["taxid"] = tax_df["taxid"].astype(int)
    df["taxid"] = df["taxid"].astype(int)
    ret_df = pd.merge(df, tax_df[["taxid", "genetic_code"]], on="taxid", how="left")

    assert ret_df.shape[0] == df.shape[0]  # no duplicate row in tax_df
    assert ret_df["genetic_code"].isnull().sum() == 0, "genetic code undifined in {} rows".format(ret_df["genetic_code"].isnull().sum())

    print("DONE : add genetic code")
    unique, counts = np.unique(ret_df["genetic_code"], return_counts=True)
    print("\t", dict(zip(unique, counts)))

    return ret_df


def update_db(df, dbFilepath):
    print("START: update {}".format(dbFilepath))
    dc = DbController(dbFilepath)
    insertCount = 0
    for _, row in df.iterrows():
        success = dc.insert_species(row["taxid"], row["domain"], row["organism_name"], row["ftp_path"], row["ftp_basename"], row["genetic_code"])
        if success:
            insertCount += 1
    print("DONE:  insert {}/{} rows".format(insertCount, df.shape[0]))


def main(filepath_lst, domain_lst, taxFilepath, dbFilepath):
    df = get_df(filepath_lst, domain_lst)
    df = add_genetic_code(df, taxFilepath)
    update_db(df, dbFilepath)


if __name__ == "__main__":
    filepath_lst = ["./data/archaea_assembly_summary.txt", "./data/bacteria_assembly_summary.txt"]
    domain_lst = ["archaea", "bacteria"]
    taxFilepath = "/data/mitsuki/data/taxonomy/nodes.processed"
    dbFilepath = "../db/altorf.db"
    main(filepath_lst, domain_lst, taxFilepath, dbFilepath)
