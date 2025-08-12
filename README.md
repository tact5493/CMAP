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

## Branch Strategy
### main:
稼働を確認．（2025.03.--）
#### File configuration
```
CMAP/
├── README.md
├── requirement.txt
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
└── src/
```

### experiment:
リファクタリング中
コードの細分化

#### File configuration
```
CMAP/
├── README.md
├── requirement.txt
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
└── src/
    ├── scripts/
    │   ├── get_data.py
    │   ├── main.py
    │   └── plot.py
    └── utils/
        ├── cal_vacuum_field.py
        ├── decide_filament_position.py
        ├── get_coordinate_params.py
        ├── get_physical_constants.py
        ├── load_experiment_data.py
        ├── log_analysis_result.py
        ├── mfield_sub.py
        ├── Parameter.py
        ├── prepare_time_arrays.py
        ├── process_pickup_coil_data.py
        ├── read_or_generate_A0.py
        ├── save_and_plot_results.py
        ├── save_best_iteration_files.py
        ├── set_elf0.py
        ├── setup_coil_and_injection.py
        ├── setup_directories.py
        └── setup_pf_coil.py
```

### feauture:
新機能の開発

## Contact
Takuto MOTOKI: takuto.505521@gmail.com