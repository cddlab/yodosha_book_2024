# [KarmaDock](https://github.com/schrojunzhang/KarmaDock)

## 実行環境
- NVIDIA RTX A4000
- CUDA 12.0
- NVIDIA Driver 525.116.04

## Installation
```bash
git clone https://github.com/schrojunzhang/KarmaDock.git
cd KarmaDock
mamba env create --file karmadock_env.yml
mamba activate karmadock
```

## Protocol

### 1. データ準備
KarmaDockでは分子IDとSMILESを記載したTSVファイルをインプットとして使用する

CSVからTSVフォーマットへの変換:
```bash
DATA_DIR="../../data"
# IPAB化合物
python csv2tsv.py \
    $DATA_DIR/ipab_processed.csv \
    ipab_input.csv

# FDA化合物
python csv2tsv.py \
    $DATA_DIR/enamine_fda_processed.csv \
    enamine_fda_input.csv
```

結合部位指定用のMOL2ファイル作成:
```bash
cd data/
mamba install -c openbabel openbabel
obabel -ipdbqt lig_1.pdbqt -omol2 -Olig_1.mol2

```

### 2. ドッキング実行
#### グラフデータの生成
```bash
REPO_DIR="/path/to/KarmaDock_repository"
cd /path/to/KarmaDock_repository

mamba activate karmadock
python -u $REPO_DIR/utils/virtual_screening_pipeline.py \
    --mode generate_graph \
    --ligand_smi ipab_input.smi \
    --protein_file ../../data/receptor.pdb \
    --crystal_ligand_file ../../data/lig_1.mol2 \
    --graph_dir graph_dir_ipab \
    --random_seed 2023
```
#### ドッキング計算の実行
```bash
python -u $REPO_DIR/utils/virtual_screening_pipeline.py \
    --mode vs \
    --protein_file ../../data/receptor.pdb \
    --crystal_ligand_file ../../data/lig_1.mol2 \
    --graph_dir graph_dir_ipab \
    --out_dir results_ipab \
    --score_threshold 1 \
    --batch_size 64 \
    --random_seed 2023 \
    --out_uncoorected \
    --out_corrected
```
注：`score_threshold`は全化合物の結果を得るため、デフォルトより低く設定している．

FDA化合物の場合は以下を変更：
- `ipab_input.smi` → `enamine_fda_input.smi`
- `graph_dir_ipab` → `graph_dir_fda`
- `results_ipab` → `results_fda`

#### 結果ファイル
各リガンドについて以下が生成される：
- `*_uncorrected.sdf`: 初期推論ポーズ
- `*_ff_corrected.sdf`: 分子力場による補正後
- `*_align_corrected.sdf`: RDKitアラインメントによる補正後
- `score.csv`：全リガンドのKarmaDockスコア

```bash
tree results_fda/ | head -n 4
results_fda/
├── Z1020762628_uncorrected.sdf
├── Z1020762628_ff_corrected.sdf
├── Z1020762628_align_corrected.sdf
```

### 解析例
[analyze_results.ipynb](analyze_results.ipynb)に、AutoDock Vinaとの比較を含む解析例を示す．
