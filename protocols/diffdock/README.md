# [DiffDock](https://github.com/gcorso/DiffDock)

## 実行環境
- NVIDIA RTX A4000（VRAM: 16GB）
- CUDA 12.0
- NVIDIA Driver 525.116.04

## Installation
公式リポジトリに従い，conda（またはmamba）でインストール:
```bash
git clone https://github.com/gcorso/DiffDock.git
cd DiffDock
# 再現性のため特定のコミットを使用可能
git checkout 561f70a8976e196a5caa1907d5437dd0de7cdba3
mamba env create --file environment.yml
mamba activate diffdock
```

openfoldのインストールに失敗する場合：
1. environment.ymlの以下をコメントアウト
```txt
# - pip:
# - openfold @ git+https://github.com/aqlaboratory/openfold.
git@4b41059694619831a7db195b7e0988fc4ff3a307
```
2. 別途インストール
```bash
mamba env create --file environment.yml
mamba activate diffdock
pip install openfold@git+https://github.com/aqlaboratory/openfold.git@4b41059694619831a7db195b7e0988fc4ff3a307
```


## Protocol
### 1. データ準備
DiffDockでは以下が必要：
- タンパク質: PDBファイルまたはアミノ酸配列
- リガンド: SDFファイルまたはSMILES
- 入力情報を記載したCSV
注：リガンドの3D構造は提供しても無視され、RDKitで再計算されます.

CSVファイル作成例:
```python
mol_df = pd.read_csv("/path/to/csv_file")
lig_info = mol_df.smiles.tolist()
complex_name = mol_df.Znumber.tolist()
protein_path = /path/to/receptor.pdb
df = pd.DataFrame(
    {
        "complex_name": complex_name,
        "protein_path": [protein_path] * len(lig_info),
        "ligand_description": lig_info,
        "protein_sequence": [""] * len(lig_info),
    }
)
df.to_csv("processed_diffdock_data.csv", index=False)
```

前処理済みリガンドデータからのCSV作成：
```bash
DATA_DIR="../../data"
# IPAB化合物
python make_input_csv.py \
    $DATA_DIR/ipab_processed.csv \
    $DATA_DIR/receptor.pdb \
    ipab_input.csv

# FDA化合物
python make_input_csv.py \
    $DATA_DIR/enamine_fda_processed.csv \
    $DATA_DIR/receptor.pdb \
    enamine_fda_input.csv
```

### 2. ドッキング実行
DiffDockリポジトリで以下を実行:
```bash
WORK_DIR="/path/to/yodosha_book_2024/protocols/diffdock"
cd /path/to/DiffDock_repository
mamba activate diffdock
python -m inference \
    --protein_ligand_csv $WORK_DIR/ipab_input.csv \
    --out_dir $WORK_DIR/results_ipab \
    --inference_steps 20 \
    --samples_per_complex 40 \
    --batch_size 10 \
    --actual_steps 18 \
    --no_final_step_noise
```
FDA化合物の場合は`ipab_input.csv`→`enamine_fda_input.csv`、`results_ipab`→`results_fda`に変更

結果は各リガンドのディレクトリに保存:
```bash
tree results_fda/ | head -n 13
results_diffdock/
├── Z1020762628
│   ├── rank1.sdf
│   ├── rank10_confidence-2.29.sdf
│   ├── rank1_confidence-0.63.sdf
│   ├── rank2_confidence-1.28.sdf
│   ├── rank3_confidence-1.33.sdf
│   ├── rank4_confidence-1.37.sdf
│   ├── rank5_confidence-1.59.sdf
│   ├── rank6_confidence-1.74.sdf
│   ├── rank7_confidence-1.78.sdf
│   ├── rank8_confidence-1.85.sdf
│   └── rank9_confidence-1.92.sdf
```

### 解析例
PoseBustersによる評価を実行する場合:
```
mamba activate diffdock
pip install posebusters
python calc_results.py -i results_ipab -o docked_ipab.csv
python calc_results.py -i results_fda -o docked_fda.csv
```

[analyze_results.ipynb](analyze_results.ipynb)に、AutoDock Vinaとの比較を含む解析例を示す．
