import argparse
import warnings

from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from posebusters import PoseBusters
from rdkit import Chem
from tqdm import tqdm

import sys
sys.path.append("../autodock_vina/")
from get_ligand_coordinate import calc_grid_box_info

warnings.filterwarnings('ignore')



def calc_center_dist(mol: Chem.Mol) -> float:
    # AutoDock Vinaのconfig.txtで指定したcenter座標
    CENTER_COORDINATES: np.ndarray = np.array([2.94935341, -10.73860162, -7.62855921])
    mol_noH = Chem.RemoveAllHs(mol)
    cen, _ = calc_grid_box_info(mol_noH)
    dist = np.linalg.norm(cen - CENTER_COORDINATES)
    return dist


def run_posebusters(
        sdf_file_list: List[str | Path], receptor_file: str | Path | None=None
    ) -> pd.DataFrame:
    bust_df = pd.DataFrame()
    if receptor_file is not None:
        config = "dock"
    else:
        config = "mol"

    buster = PoseBusters(config=config)
    for f in tqdm(sdf_file_list):
        try:
            _df = buster.bust([f], None, receptor_file, full_report=False)
            bust_df = pd.concat([bust_df, _df], axis=0)
        except ValueError as e:
            print(f"PoseBusters failed with {f}: {e}")
    bust_df['all_pass'] = bust_df.all(axis='columns')
    bust_df = bust_df.reset_index()
    int_cols = [
        'mol_pred_loaded',
        'sanitization',
        'inchi_convertible',
        'all_atoms_connected',
        'bond_lengths',
        'bond_angles',
        'internal_steric_clash',
        'aromatic_ring_flatness',
        'double_bond_flatness',
        'internal_energy',
        'all_pass'
    ]
    int_col_dict = {c: 'i4' for c in int_cols}
    bust_df = bust_df.astype(int_col_dict)
    del bust_df["molecule"]
    return bust_df


def get_pose_properties(sdf_file_path: str | Path) -> Dict[str, int | float | None]:
    f_stem = sdf_file_path.stem
    if len(f_stem.split("_")) > 1:
        rank = int(f_stem.split("_")[0].replace("rank", ""))
        confidence = 100
        if "confidence-" in f_stem.split("_")[1]:
            confidence = float(f_stem.split("_")[1].replace("confidence-", ""))
        else:
            confidence = float(f_stem.split("_")[1].replace("confidence", ""))

        m = Chem.SDMolSupplier(str(sdf_file_path))[0]
        if m is not None:
            m_noh = Chem.RemoveAllHs(m)
            dist = calc_center_dist(m_noh)
        else:
            dist = None
        return {
            "Znumber": sdf_file_path.parent.stem.split('_')[0],
            "rank_id": rank,
            "confidence": confidence,
            "center_distance": dist,
        }


def main() -> None:
    parser = argparse.ArgumentParser(
        usage="get and calculate diffdock results",
    )
    parser.add_argument("--input_dir", "-i", help="input sdf directory path")
    parser.add_argument("--output", "-o", help="output csv path")
    parser.add_argument("--pdb", "-p", default=None, help="protein pdb")
    args = parser.parse_args()

    sdf_dir = Path(args.input_dir)

    sdf_files = sorted(list(sdf_dir.glob("./*/*.sdf")))
    sdf_files = [f for f in sdf_files if f.stem != "rank1"]
    properties = [get_pose_properties(f) for f in sdf_files]
    df = pd.json_normalize(properties).dropna()
    bust_df = run_posebusters(sdf_files)
    bust_df = pd.concat(
        [
            bust_df,
            bust_df.file.str.split('/', expand=True).rename(
                columns={
                    0: 'directory',
                    1: 'Znumber',
                    2: 'rank_id',
                }
            ).drop(["directory"], axis=1)
        ], axis=1,
    )
    bust_df = bust_df[bust_df['rank_id'].str.contains("_")]
    bust_df['rank_id'] = bust_df['rank_id'].str.split('_', expand=True).rename(
        columns={
            0: 'rank_id',
            1: 'conf',
        }
    ).drop('conf', axis=1)
    bust_df['rank_id'] = bust_df['rank_id'].str.replace("rank", "").astype('i4')
    df = df.merge(bust_df, on=['Znumber', 'rank_id'])
    df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
