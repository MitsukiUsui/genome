import pandas as pd
import numpy as np


def get_df(filepath_lst, domain_lst):
    assert len(filepath_lst) == len(domain_lst)

    dct_lst = []
    for filepath, domain in zip(filepath_lst, domain_lst):
        df = pd.read_csv(filepath, delimiter='\t', skiprows=1, dtype="object")
        print("DONE : load {} species from {}".format(df.shape[0], filepath))

        # filtering
        filtered_df = df[(df["refseq_category"] == "representative genome") & (df["assembly_level"] == "Complete Genome")]
        taxid_set = set(filtered_df["taxid"]) # choose only one from one taxid

        for taxid in taxid_set:
            dct = {}
            record = (filtered_df[filtered_df["taxid"] == taxid]).iloc[0, :]  # 先頭のレコードのみを取得

            column_lst = ["# assembly_accession", "taxid", "organism_name", "asm_name", "ftp_path"]  # 残すカラムを選択
            for column in column_lst:
                dct[column] = record[column]
            dct["ftp_basename"] = record["ftp_path"].split('/')[-1]
            dct["domain"] = domain
            dct_lst.append(dct)
        print("DONE : pick {} species".format(len(taxid_set)))

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


def main(filepath_lst, domain_lst, taxFilepath, outFilepath):
    df = get_df(filepath_lst, domain_lst)
    df = add_genetic_code(df, taxFilepath)
    df.to_csv(outFilepath, index=False)
    print("DONE : output csv to {}".format(outFilepath))


if __name__ == "__main__":
    filepath_lst = ["./data/archaea_assembly_summary.txt", "./data/bacteria_assembly_summary.txt"]
    domain_lst = ["archaea", "bacteria"]
    taxFilepath = "/data/mitsuki/data/taxonomy/nodes.processed"
    outFilepath = "../target.csv"
    main(filepath_lst, domain_lst, taxFilepath, outFilepath)
