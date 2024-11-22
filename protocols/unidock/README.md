# [Uni-Dock](https://github.com/dptech-corp/Uni-Dock)

## 実行環境
- NVIDIA RTX A4000
- CUDA 12.0
- NVIDIA Driver 525.116.04

## Installation
```bash
mamba create -n unidock_env unidock -c conda-forge
mamba activate unidock_env
```
注：本プロトコールではv1.1.2を使用


## Protocol

### 1. データ準備
必要なファイル：
- リガンド：PDBQTファイル（Meekoで作成）とそのパスリスト
- タンパク質：PDBQTファイル（AutoDock Vinaと同様）

```bash
# Meeko, RDKitがインストールされていなければインストールして以下を実行
mamba activate unidock_env
mamba install -c conda-forge meeko rdkit

DATA_DIR="../../data"
# IPAB化合物
python smi2pdbqt.py \
    -i $DATA_DIR/ipab_processed.csv \
    -m Znumber \
    -d ipab_input
find ./ipab_input/ -type f -name "*.pdbqt" > ipab_indata.txt

# FDA化合物
python smi2pdbqt.py \
    $DATA_DIR/enamine_fda_processed.csv \
    -m Znumber \
    -d enamine_fda_input
find ./enamine_fda_input/ -type f -name "*.pdbqt" > enamine_fda_indata.txt
```


### 2. ドッキング実行
Uni-Dockがインストールされた環境をアクティベートして，以下のコマンドを実行する.

```bash
mamba activate unidock_env

# IPAB化合物
mkdir results_ipab
DATA_DIR="../../data"
unidock \
    --receptor $DATA_DIR/receptor.pdbqt \
    --ligand_index indata_ipab.txt \
    --dir results_ipab \
    --config config.txt

# IPAB化合物
mkdir results_fda
DATA_DIR="../../data"
unidock \
    --receptor $DATA_DIR/receptor.pdbqt \
    --ligand_index indata_enamine_fda.txt \
    --dir results_fda \
    --config config.txt
```

結果は各ディレクトリにPDBQTファイルとして保存される．


### 解析例
[analyze_results.ipynb](analyze_results.ipynb)に、AutoDock Vinaとの比較を含む解析例を示す．
