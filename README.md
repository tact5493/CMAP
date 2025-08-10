# CMAP

## Instruction
C-MAP (CHI Magnetic Analysis Program) はCHIプラズマの磁気面再構成に使用される．
Taylor 緩和後($\lambda=const.$) を仮定している．
仮定条件については T.Motoki, Master thesis (2024). を参照．

### Recommended environment
OS : mac os 
Version : Python 3.9.6
Editor : VS code
> CHI 実験結果を利用するため，QUESTサーバーにマウントを推奨

### File configuration
```
CMAP/
├── README.md
├── README.md
├── src/
│   ├── scripts/
│   │   ├── get_data.py
│   │   ├── main.py
│   │   ├── plot.py
│   │   └── ...
│   └── utils/
│       ├── get_coordinate_params.py
│       ├── get_physical_constants.py
│       ├── load_experiment_data.py
│       ├── mfield_sub.py
│       ├── Parameter.py
│       └── ...
├── data/
│   ├── .gitkeep
│   ├── ele_posi.csv
│   ├── electrode.csv
│   ├── PFdata.csv
│   ├── quest_lines.txt
│   └── Qvessel.png
├── Image/
│   ├── .gitkeep
│   ├── psi_cont_19.450_small.png
│   └── psi_cont_19.450.png
└── ...
```

### Input
- 磁気コイル
- 入射電流
- トロイダル電流
- TF, PF コイル電流

### Output
- 2次元ポロイダル磁束分布（磁気面）
- 2次元ポロイダル電流分布


## Mechanism


## How to use


## Example output

![Output_example](/Image/psi_cont_19.450_small.png) 


## Contact
Takuto MOTOKI: takuto.505521@gmail.com