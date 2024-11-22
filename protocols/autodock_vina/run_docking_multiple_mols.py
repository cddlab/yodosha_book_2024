import argparse
import os
import subprocess

from pathlib import Path
from typing import List

from tempfile import TemporaryDirectory

import pandas as pd

from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm



def main() -> None:
    parser = argparse.ArgumentParser(
        usage="vina docking with executable binary",
    )
    parser.add_argument("--config", "-c", help="vina config")
    parser.add_argument("--input", "-i", help="input csv path")
    parser.add_argument("--mol_id", "-m", help="compound id in input csv")
    parser.add_argument("--out_dir", "-d", help="output mol directory")
    parser.add_argument("--output", "-o", help="output csv path")
    parser.add_argument("--pdbqt", "-p", help="protein pdbqt")
    parser.add_argument("--num_confs", "-n", default=-1, type=int, help="number of conformers")
    args = parser.parse_args()
    if not Path(args.out_dir).exists():
        Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    smiles = df.smiles.values
    # 化合物IDを取得し，ドッキング結果のPDBQTファイル名に用いる
    c_ids = df[args.mol_id].values
    vina_scores = []
    ps = AllChem.ETKDGv3()
    # 再現性のためにシード固定
    ps.randomSeed = 0xF00D
    num_conf = args.num_confs
    for c, s in tqdm(zip(c_ids, smiles), total=len(smiles)):
        try:
            score: float = 10.0
            m0 = Chem.MolFromSmiles(s)
            if m0 is not None:
                m_h = Chem.AddHs(m0)
                if num_conf > 0:
                    AllChem.EmbedMultipleConfs(m_h, num_conf, ps)
                else:
                    AllChem.EmbedMolecule(m_h, ps)
                preparator = MoleculePreparation(flexible_amides=True)
                mol_setups: List[meeko.molsetup.RDKitMoleculeSetup] = preparator.prepare(m_h)
                ligand_pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setups[0])
                if is_ok:
                    ...
                else:
                    print(error_msg)

                with TemporaryDirectory() as tmpdir:
                    lig_input = os.path.join(tmpdir, "tmp.pdbqt")
                    lig_output = os.path.join(args.out_dir, f"{c}.pdbqt")
                    with open(lig_input, "w") as f:
                        f.write(ligand_pdbqt_string)
                    cmd = [
                        "vina",
                        f"--config={args.config}",
                        f"--ligand={lig_input}",
                        f"--receptor={args.pdbqt}",
                        f"--out={lig_output}",
                        f"--verbosity=1",
                    ]
                    p = subprocess.run(cmd, capture_output=True, timeout=60)
                    if p.returncode > 0:
                        raise RuntimeError(p.returncode, p.args, p.stderr.decode())
                    pdbqt_mol = PDBQTMolecule.from_file(lig_output)
                    score = [m.score for m in pdbqt_mol][0]

        except:
            # ドッキングに失敗した分子にスコアをいれる
            score = 10
        vina_scores.append(score)

    df["score"] = vina_scores
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
