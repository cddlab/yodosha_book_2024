import argparse
import os

from pathlib import Path
from typing import List

import pandas as pd

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="smi2pdbqt.py",
        usage="preparing pdbqt files from csv containing smiles column",
    )
    parser.add_argument("--input", "-i", help="input csv path")
    parser.add_argument("--mol_id", "-m", help="compound id in input csv")
    parser.add_argument("--out_dir", "-d", help="output mol directory")
    args = parser.parse_args()
    if not Path(args.out_dir).exists():
        Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    smiles = df.smiles.values
    c_ids = df[args.mol_id].values
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xF00D
    for c, s in tqdm(zip(c_ids, smiles), total=len(smiles)):
        try:
            m0 = Chem.MolFromSmiles(s)
            if m0 is not None:
                m_h = Chem.AddHs(m0)
                AllChem.EmbedMolecule(m_h, ps)
                preparator = MoleculePreparation(flexible_amides=True)
                mol_setups: List[meeko.molsetup.RDKitMoleculeSetup] = preparator.prepare(m_h)
                ligand_pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
                if is_ok:
                    ligand_pdbqt_string = preparator.write_pdbqt_string()
                else:
                    print(error_msg)
                lig_input = os.path.join(args.out_dir, f"{c}.pdbqt")

                with open(lig_input, "w") as f:
                    f.write(ligand_pdbqt_string)
        except ValueError as e:
            print(f"{s}: {e}")


if __name__ == "__main__":
    main()
