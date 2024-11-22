import sys

from pathlib import Path

import pandas as pd


def main() -> None:
    raw_csv_file = sys.argv[1]
    p_path = sys.argv[2]
    input_file_path = sys.argv[3]
    # 絶対パスに変換する
    p_path = str(Path(p_path).resolve())

    mol_df = pd.read_csv(raw_csv_file)
    lig_info = mol_df.smiles.tolist()
    c_name = mol_df.Znumber.values
    df = pd.DataFrame(
        {
            "complex_name": c_name,
            "protein_path": [p_path] * len(lig_info) ,
            "ligand_description": lig_info,
            "protein_sequence": [""] * len(lig_info),
        }
    )
    df.to_csv(input_file_path, index=False)


if __name__ == "__main__":
    main()

