# yodosha_book_2024

このリポジトリは [AlphaFold 時代の構造インフォマティクス実践ガイド (羊土社)](https://www.yodosha.co.jp/yodobook/book/9784758122764/)の「第3章 立体構造によるタンパク質の機能推定 2.分子ドッキング法」で紹介したプロトコールのコードを掲載している．なお，書籍本文と一部変更している．

## Table of Contents
1. [File Structure](#file-structure)
2. [Datasets and Experimental Setting](#datasets-and-experimental-setting)
3. [Protocols](#protocols)

## File Structure
```
.
├── LICENSE
├── README.md
├── data
└── protocols
    ├── autodock_vina
    ├── diffdock
    ├── karmadock
    └── unidock
```

## Datasets and Experimental Setting

本プロトコールでは，[第2回IPABコンテスト](http://www.ipab.org/eventschedule/contest/contest2) ([論文](https://www.nature.com/articles/s41598-017-10275-4))で題材となったcYesに対する化合物を例に各種ドッキングの実行方法を紹介する．

使用したデータセット詳細 ([data](data)):
- タンパク質: 参加チームGroup10がホモロジーモデリングで作成したモデル
- リガンド(lig_1, lig_2, lig_3): 同チームが選出した代表化合物
- 複数化合物ドッキング例: IPAB参加者が提出した化合物から381化合物を選択
- デコイ: Enamine社 FDA approved Drugs Collections（[version 9，November 2023](https://enamine.net/compound-libraries/bioactive-libraries/fda-approved-drugs-collection)，以下FDA 化合物）の1123 化合物を選択．金属を含有する化合物を除いた1118化合物を使用．

通常のバーチャルスクリーニングは，ある化合物群に対し一定の活性閾値（例：IC50 <1 μM）を設定して，ドッキングスコアなどを指標に活性の有無を判別する．本プロトコールでは化合物数/計算時間を考慮して，以下の擬似的な正例・負例を判別するタスクを設定した:
- 正例: IPAB提出化合物 (コンペティション参加者が活性を示すことを期待して選抜した化合物)
- 負例: FDA化合物 (多様な構造を有する化合物群であり，多くが活性を示す可能性は低いと仮定)

解析例では，各ドッキング手法の正例 (IPAB化合物)および負例 (FDA化合物)の判別性能評価の一例を紹介している．

## Protocols
- [AutoDock Vina](protocols/autodock_vina/README.md)
- [DiffDock](protocols/diffdock/README.md)
- [KarmaDock](protocols/karmadock/README.md)
- [Uni-Dock](protocols/unidock/README.md)
