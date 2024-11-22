import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.MolStandardize import rdMolStandardize


def standardize_mol(mol: Chem.Mol) -> Chem.Mol:
    try:
        _m = rdMolStandardize.ChargeParent(mol)
        return _m
    except ValueError as e:
        print(e)
        return None


def main() -> None:
    input_sdf = "../../data/Enamine_FDA_approved_Drugs_1123cmpds_20231109.sdf"
    df = PandasTools.LoadSDF(input_sdf)
    df = df.rename(columns={"Catalog ID": "Znumber"})
    print(df.shape)

    # 金属が入っている分子のZnumberを読み込む
    removal_znumbers = []
    with open("../../data/enamine_fda_removal_znumbers.txt", "r") as f:
        for r in f.readlines():
            r = r.rstrip().split(',')
            removal_znumbers.append(r[0])
    # 金属が入っている分子を除去
    df = df[~df.Znumber.isin(removal_znumbers)]

    df["smiles"] = df["ROMol"].apply(lambda x: Chem.MolToSmiles(x))
    print(df.shape)
    df["StMol"] = df["ROMol"].apply(lambda x: standardize_mol(x))
    print(df.shape)
    df["st_smiles"] = df["StMol"].apply(lambda x: Chem.MolToSmiles(x))
    # 前処理が失敗した行を削除
    df = df.dropna(subset="st_smiles").reset_index()
    print(df.shape)
    df.loc[:, ["Znumber", "st_smiles"]].rename(columns={"st_smiles": "smiles"}).to_csv(
        "../../data/enamine_fda_processed.csv", index=False
    )

if __name__ == "__main__":
    main()
