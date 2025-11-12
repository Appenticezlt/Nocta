# Nocta
一个面向蛋白质序列与结构的分析管线，通过整合各种工具提供从结构比对、层级聚类、能量分析到序列位点香农熵计算的完整流程。该项目支持大规模 PDB 文件的批量处理，并可以选择生成可视化结果，包括热图、聚类树、能量分布和熵柱状图等。
## 1.目录结构

```text
Nocta/
├── data/                  # 待分析 PDB 文件目录
├── results/               # 默认输出目录
├── main.py                # 主入口脚本
├── TMalign.py             # TMalign 批量比对模块
├── energy.py              # Rosetta 能量计算模块
├── shannon.py             # Shannon entropy 计算模块
├── cluster.py             # 层级聚类模块
├── visualize.py           # 可视化模块
├── README.md              # 使用说明
└── requirements.txt       # 需要安装的依赖包
└── setup.sh               #自动化部署环境脚本 
```

## 2.部署说明

### 2.1 脚本安装

使用 Conda 管理 Python 环境，Nocta 提供了自动化安装脚本`setup_env.sh`：
P.S 如果遇到syntax error等错误，可能是因为文件在 Windows 下编辑过 ，可能含有 **CRLF 换行符**，在 Linux 下会报错，运行以下命令可以解决：

```
dos2unix setup.sh
```
之后运行：

```bash
bash setup_env.sh
```

### 2.2手动方式（可选）

如果不使用脚本，也可手动安装：

1.创建 Conda 环境：

```
conda create -n analysis python=3.10 -y
conda activate analysis
```

2.安装 Python 依赖：

```
pip install -r requirements.txt
```
3.安装USalign和TMalign
```
conda install bioconda::tmalign
conda install bioconda::usalign
```
4.安装 PyRosetta：

```
# 在 Python 交互式环境中运行，其中pyrosetta在requirements.txt中已经安装
import pyrosetta_installer
pyrosetta_installer.install_pyrosetta()
```

5.验证 PyRosetta 是否正常：

```
import pyrosetta
pyrosetta.init("-mute all")
print("PyRosetta successfully initialized!")
```

安装完成后，你就可以运行 `main.py` 进行结构分析、能量计算和 Shannon 熵计算。

## 各模块文件介绍

- **TMalign.py**：利用TMalign和USalign进行批量 PDB 比对，输出 RMSD 和 TM-score 矩阵。
- **energy.py**：使用 PyRosetta 中的计算单个 PDB 或批量 PDB 的总能量及结合能。
- **shannon.py**：计算每个位点的 Shannon entropy，可自动提取 PDB 序列。
- **cluster.py**：基于 RMSD 或 TM-score 的层级聚类，生成簇标签。
- **visualize.py**：提供可视化函数，包括热图、聚类树、能量分布和 entropy 图表。


### 参数说明

| 参数                    | 说明                                     |
| ----------------------- | ---------------------------------------- |
| `--input`               | 输入 PDB 文件目录                        |
| `--output`              | 输出目录                                 |
| `--useTMalign`          | 是否运行结构比对                         |
| `--alignnum_workers`    | 结构比对并行线程数                       |
| `--savetmalign_results` | 是否保存 TMalign 输出结果                |
| `--doClusterRMSD`       | 是否进行基于 RMSD 的层级聚类             |
| `--doClusterTMscore`    | 是否进行基于 TM-score 的层级聚类         |
| `--cluster_threshold`   | 层级聚类阈值（距离或 1-TM-score）        |
| `--cluster_method`      | 聚类方法（single/complete/average/ward） |
| `--useEnergy`           | 是否运行 Rosetta 能量计算                |
| `--saveenergy_results`  | 是否保存能量计算结果                     |
| `--doentropy`           | 是否计算 Shannon 熵                      |
| `--visualize_*`         | 可视化开关，根据模块生成对应图表         |
| `--visualize_combined`  | 汇总各模块可视化图表，生成 combined_plot |

### ￥_WORKERS 选择建议

| num_workers | CPU 使用情况        | 并行度   | 建议             |
| ----------- | ------------------- | -------- | ---------------- |
| 1           | 顺序执行，最慢      | 无并行   | 调试用           |
| 4           | 同时跑 4 个 USalign | 4× 加速  | 轻量任务         |
| 8           | 同时跑 8 个         | 高并行度 | 16核机器推荐     |
| 16+         | CPU 饱和            | 最高并发 | 超算或多核服务器 |

## 5.使用指南

### 5.1 仅运行TMalign 比对（不输出结果）

```bash
python main.py --input data/ --output results/ --useTMalign
```

### 5.2.仅运行 Rosetta 能量计算

```
python main.py --input data/ --output results/ --useEnergy
```

可选：

- `--energynum_workers 8`：多线程加速能量计算
- `--saveenergy_results`：保存能量计算结果

------

### 3. 仅计算 Shannon 熵

```
python main.py --input data/ --output results/ --doentropy
```

- 自动读取 PDB 文件提取序列，如果已存在 `sequences.fasta`，将直接使用
- 结果保存为 Excel 文件 `Shannon_entropy.xlsx`

------

### 4.层级聚类（可选）

- 在运行过 TMalign 后，可以进行聚类分析：

```
python main.py --input data/ --output results/ --useTMalign \
               --doClusterRMSD --doClusterTMscore \
               --cluster_method average --cluster_threshold 2.0
```

- `--doClusterRMSD`：基于 RMSD 的层级聚类
- `--doClusterTMscore`：基于 TM-score 的层级聚类
- `--cluster_method`：聚类方法（single / complete / average / ward）
- `--cluster_threshold`：距离阈值（RMSD 或 1-TM-score）

------

### 5.可视化模块（需先生成数据）

- TMalign 热图：

```
python main.py --visualize_tmalign
```

- 聚类树状图：

```
python main.py --visualize_cluster
```

- 能量分布图：

```
python main.py --visualize_energy
```

- Shannon 熵柱状图：

```
python main.py --visualize_entropy
```

- 汇总可视化：

```
python main.py --visualize_combined
```

> 注意：可视化模块需要对应计算结果已经存在，否则会报错。

------

### 6.一条命令运行全部流程

```
python main.py --input data/ --output results/ \
               --useTMalign --savetmalign_results \
               --doClusterRMSD --doClusterTMscore \
               --useEnergy --saveenergy_results \
               --doentropy \
               --visualize_tmalign --visualize_cluster \
               --visualize_energy --visualize_entropy \
               --visualize_combined
```

- 这条命令会执行：
  1. 结构比对 (TMalign / USalign)
  2. 层级聚类
  3. 能量计算
  4. Shannon 熵计算
  5. 可视化各模块并生成汇总图

## 🔹 作者

未烁寿 (博士生研究方向：我也不知道最后会是什么方向)

GitHub: https://github.com/Appenticezlt/Nocta
