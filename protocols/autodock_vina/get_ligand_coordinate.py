import sys

import numpy as np

from rdkit import Chem


def get_coords(mol: Chem.Mol) -> np.ndarray:
    """重原子の座標を取得する関数"""
    mol = Chem.RemoveHs(mol)
    conf = mol.GetConformer(0)
    if conf is None:
        raise RuntimeError("cannot get conformer from mol")
    coords = np.array(conf.GetPositions())
    return coords


def calc_grid_box_info(mol) -> tuple[np.ndarray, np.ndarray]:
    """分子の重心・分子を包含するgrid boxサイズを求める関数"""
    coords = get_coords(mol)
    center = np.average(coords, axis=0)
    coord_min = np.min(coords, axis=0)
    coord_max = np.max(coords, axis=0)
    ligand_size = coord_max - coord_min
    return center, ligand_size


if __name__ == "__main__":

    mol_file_path = list(sys.argv[1:])
    mols = [Chem.SDMolSupplier(f)[0] for f in mol_file_path]
    centers = []
    grid_sizes = []
    for mol in mols:
        c, g = calc_grid_box_info(mol)
        centers.append(c)
        grid_sizes.append(g)

    print("center coordinates")
    print(np.average(centers, axis=0))
    print("ligand size")
    print(np.max(grid_sizes, axis=0))
