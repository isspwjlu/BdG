# BdG DOS — 基于全 Wannier 轨道的 Bogoliubov-de Gennes 态密度计算

## 目录

1. [概述](#1-概述)
2. [理论原理](#2-理论原理)
3. [数值方法](#3-数值方法)
4. [程序架构与实现](#4-程序架构与实现)
5. [使用指南](#5-使用指南)
6. [技术细节](#6-技术细节)
7. [参考文献](#7-参考文献)

---

## 1. 概述

本程序通过构建并对角化全 Wannier 轨道基的 Bogoliubov-de Gennes (BdG) 哈密顿量，计算晶体材料在超导态下的态密度 (Density of States, DOS)。程序基于 Wannier90 输出的紧束缚参数 (`hr.dat`)，自动探测 Wannier 函数数目，可应用于任意晶体材料。

主要功能：

- 在全 Wannier 轨道-自旋空间中构造 BdG 哈密顿量（无轨道对称性约化）
- 支持面外 (z) 和面内 (x) 两个方向的 Zeeman 磁场
- 面内 Zeeman 效应可配置增强因子 `G_EFFECT`，适用于具有强自旋-轨道耦合的材料
- 基于 Python 多进程的 k 点并行计算
- 结果输出为数据文件 (`.dat`) 和图像文件 (`.png`)

## 2. 理论原理

### 2.1 BdG 哈密顿量

Bogoliubov-de Gennes 形式是处理超导问题的平均场框架，它将超导体的准粒子激发描述为电子和空穴的混合。在 Nambu 表象中，BdG 哈密顿量为：

```
H_BdG(k) = ⎡  A(k)         D    ⎤
           ⎣  -D      -conj(A(-k)) ⎦
```

其中电子块 `A(k) = H(k) + H_Z − μ·I` 包含正常态哈密顿量、Zeeman 项和化学势；配对块 `D` 描述超导关联。

### 2.2 Nambu 表象

Nambu 基矢将产生和湮灭算符组合为一个扩展的旋量：

```
Ψ(k) = (c_{k,1}, c_{k,2}, ..., c_{k,N}, c^{†}_{-k,1}, ..., c^{†}_{-k,N})^T
```

其中 N = nwannier 是 Wannier 函数数目。在该基下，BdG 哈密顿量是 2N × 2N 的厄密矩阵。前 N 维对应电子（粒子）空间，后 N 维对应空穴空间。这种构造自动保证了谱的粒子-空穴对称性：若 E 是本征值，则 −E 也必是本征值。

### 2.3 正常态：Wannier 插值

正常态哈密顿量 H(k) 通过 Wannier90 输出的紧束缚参数进行傅里叶插值：

```
H(k) = Σ_R H(R) · exp(2πi · k · R)
```

其中 R 是实空间格矢，H(R) 是实空间跳跃积分矩阵。该插值在 DFT 计算使用的网格上精确重现第一性原理能带，并通过平滑傅里叶变换获得任意 k 点的哈密顿量。

### 2.4 超导配对项

程序采用 s 波单态配对各向同性能隙：

```
D = Δ · (iσ_y) ⊗ I_orb
```

其中 Δ 是超导能隙大小，iσ_y 确保自旋单态配对对称性，I_orb 是轨道空间单位矩阵。在 BdG 矩阵中配对项出现在非对角块，空穴部分为 −D 以保证粒子-空穴对称性。

### 2.5 Zeeman 耦合

外加磁场的 Zeeman 项为：

```
H_Z = (g·μ_B/2) · [G_x·σ_x·B_x + G_y·σ_y·B_y + σ_z·B_z] ⊗ I_orb
```

其中 g = 2.0 是电子自旋 g 因子，μ_B = 5.788×10⁻⁵ eV/T 是玻尔磁子。面内分量 (Bx, By) 乘以增强因子 G_x = G_y = G_EFFECT，面外分量 (Bz) 使用裸 g 因子。该构造通过 Kronecker 积 (Pauli) ⊗ I_orb 将自旋 Zeeman 分裂投影到所有轨道上。

### 2.6 自旋-轨道耦合的有效增强

对于具有强自旋-轨道耦合的体系，自旋被锁定在特定方向（如 Ising SOC 中锁定在面外），这会修改有效 g 因子。本程序通过 G_EFFECT 参数建模：

- **G_EFFECT = 1**：标准 Zeeman 效应（裸电子 g = 2），适用于常规体系
- **G_EFFECT > 1**：有效增强的面内 Zeeman，适用于强 SOC 体系

该参数仅作用于面内方向，面外方向不受影响。

### 2.7 面外与面内 Zeeman 的物理差异

#### 2.7.1 DOS 相干峰劈裂行为：Bz vs Bx

面外 (Bz) 和面内 (Bx) 磁场对 BdG 能谱的影响有本质区别，这反映在 DOS 的相干峰行为上。

**Bz（面外）：平移 + 劈裂**

σ_z 在自旋基下是对角的。`H_Z ∝ σ_z·Bz` 将自旋↑和自旋↓通道的化学势反向移动，准粒子激发谱变为：

```
E ≈ ±√(ξ² + Δ²) ± h_z,   h_z = (g·μ_B/2)·Bz
```

两支相干峰对称劈裂为四支，只要 h_z < Δ 就能观察到清晰的劈裂。

**Bx（面内）：退配对 + 压缩**

σ_x 是**非对角**的，且与配对项 iσ_y 满足**反对易关系** `{σ_x, iσ_y} = 0`。这意味着面内磁场产生**Pauli 顺磁退配对效应**：

- 不是能级平移，而是超导能隙本身被压缩：`Δ_eff(B) < Δ₀`
- 当面内 Zeeman 能超过 Pauli 极限 `h_x > Δ₀/√2` 时，能隙完全关闭
- 表现为相干峰的**坍缩**而非劈裂

**能标对比（G_EFFECT = 100）：**

| 场 | 有效 Zeeman 能 | h_eff/Δ |
|----|--------------|---------|
| Bz = 5 T | h_z = 0.29 meV | 1.45× |
| Bx = 10 T | h_x = 57.9 meV | **289×** |
| Bx = 50 T | h_x = 289.4 meV | **1447×** |

Bz=5T 时 h_z ~ Δ，体系仍处于超导态，可观察到劈裂。Bx 即使取最小值 10T，因 G_EFFECT=100 放大后 h_x >> Δ，超导配对已被完全压制，体系进入正常态，相干峰已不存在。

若设 G_EFFECT=1，Bx 和 Bz 在相同场强下的 DOS 行为完全对称。

#### 2.7.2 面内磁场方向：Bx vs By vs 对角线

**结论：计算结果完全相同。**

原因在于 Zeeman 项的构造：

```python
sigma_B = SIGMAX * (g_effect * Bx) + SIGMAY * (g_effect * By) + SIGMAZ * Bz
```

- σ_x 和 σ_y 具有相同的本征值（±1）
- 面内分量 By 和 Bx 使用相同的增强因子 G_EFFECT
- 对于相同大小的面内磁场，无论沿 x、y 还是 xy 对角线方向，有效 Zeeman 分裂大小相同：
  ```
  |h_eff| = G_EFFECT · (g·μ_B/2) · √(Bx² + By²)
  ```
- σ_x 和 σ_y 都与 iσ_y 反对易，物理机制相同

区别仅在于自旋量子化轴的方向，但 BdG 谱具有旋转不变性，不影响 DOS 计算结果。

## 3. 数值方法

### 3.1 Monkhorst-Pack k 点网格

布里渊区积分采用 Monkhorst-Pack 方法：

```
k_i = (n_i + 0.5) / N_i，n_i = 0, 1, ..., N_i−1
```

网格点偏移半格点位置，避免落在高对称点上，从而提高积分收敛效率。网格密度由 `K_MESH = [N_kx, N_ky, N_kz]` 控制，应根据材料的倒格矢长度确定。

### 3.2 H(k) 的傅里叶插值

预计算所有 k 点与 R 格矢之间的相位因子：

```
φ_{k,R} = exp(2πi · k · R)
```

然后对每个 k 点求和得到 H(k)：

```
H(k) = Σ_R H(R) · φ_{k,R}
```

H(R) 以稀疏矩阵 (CSR) 格式存储，仅在求和时转为稠密矩阵，以平衡存储效率和计算速度。

### 3.3 BdG 矩阵构建

构建完整的 2N × 2N BdG 矩阵的步骤：

1. **厄密对称化**：H_sym(k) = (H(k) + H†(k)) / 2，消除数值误差
2. **Zeeman 项**：按磁场分量分别计算面外（裸 g）和面内（增强因子）贡献
3. **电子块**：A = H_sym + H_Z − μ·I
4. **空穴块**：−conj(H(−k) + H_Z) + μ·I
5. **配对项**：布置在非对角块

### 3.4 矩阵对角化

使用 SciPy 的 `linalg.eigvalsh` 求解特征值：

- 专为厄密矩阵优化的 LAPACK 求解器 (ZHEEVR)
- 仅计算特征值不计算特征向量，节省内存和计算时间
- 对每个 k 点返回 2N 个本征值（准粒子激发谱）

### 3.5 态密度的 Lorentzian 展宽

离散的能谱通过 Lorentzian 函数展宽为连续 DOS：

```
N(E) = (1/N_k) · Σ_{k,n} η / [π · ((E − E_{k,n})² + η²)]
```

利用 NumPy 广播实现向量化：

```
diff = eigenvalues[:, :, None] − energy_grid[None, None, :]
lorentz = η / (diff² + η²)
dos = sum(lorentz, axis=(0,1)) / (π · N_k)
```

展宽参数 η 控制分辨率：η 越大谱线越平滑但细节损失越多，η 越小特征越尖锐但需要更密的 k 网格以消除涨落噪声。典型值取 Δ 的百分之几。

### 3.6 DOS 高能区归一化

由于 Lorentzian 展宽在有限能量范围（E_RANGE）处被截断，不同磁场下 Zeeman 分裂的大小不同，导致截断处 Lorentzian 尾巴的贡献不一致。这反映为不同 B 场的 DOS 谱在高能区存在垂直偏移。

为解决此问题，程序实现了高能尾均值归一化：

1. 取完整对称 DOS 数组正能侧末尾一定比例（默认 10%）的能量区间
2. 计算该区间内的平均 DOS 值，作为该 B 场下的归一化基准
3. 将整条 DOS 谱除以该基准值

归一化后 DOS 为无量纲量，高能区趋近于 1，各 B 场曲线在高能区对齐。

## 4. 程序架构与实现

### 4.1 模块结构

```
bdg_dos/
├── __init__.py          # 包标记
├── config.py            # 物理与数值参数配置文件
├── hr_io.py             # Wannier90 hr.dat 文件读取
├── kpoints.py           # Monkhorst-Pack k 点网格生成
├── bdg.py               # BdG 哈密顿量构建与对角化
├── dos.py               # 态密度计算、绘图、数据导出
├── main.py              # 主程序入口（并行调度）
├── README.md            # 英文说明文档
└── README_CN.md         # 中文技术文档（本文档）
```

### 4.2 模块详解

**config.py**：集中管理所有可调参数

```python
# 物理参数
DELTA = 0.2           # meV，超导能隙
FERMI_LEVEL = -2.0697 # eV，化学势
G_FACTOR = 2.0        # 电子自旋 g 因子
G_EFFECT = 100        # 面内 Zeeman 增强因子

# 数值参数
ETA = 0.05 * DELTA    # meV，Lorentzian 展宽
E_RANGE = 3 * DELTA   # meV，能量扫描范围
N_ENERGY = 2000       # 能量网格点数

# k 网格
K_MESH = [1000, 1000, 1]

# 磁场参数
B_FIELDS_Z = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
B_FIELDS_X = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]

# DOS 归一化
NORMALIZE_BACKGROUND = True    # 高能尾均值归一化
NORMALIZE_FRACTION = 0.1       # 参考区间比例
```

**hr_io.py**：解析 Wannier90 hr.dat 文件，返回稀疏 H(R) 矩阵和格矢列表。

**kpoints.py**：生成 Monkhorst-Pack 网格并提供分批迭代器，用于大网格的内存管理。

**bdg.py**：核心计算模块

| 函数 | 功能 |
|------|------|
| `precompute_phases()` | 预计算 exp(2πi·k·R) 相位因子 |
| `compute_Hk_from_phases()` | 单 k 点 H(k) 傅里叶插值 |
| `build_BdG_matrix_full()` | 构建完整的 2N×2N BdG 矩阵 |
| `diagonalize_BdG()` | 厄密矩阵特征值求解 |

**dos.py**：DOS 处理与输出模块

| 函数 | 功能 |
|------|------|
| `compute_dos_vectorized()` | 向量化 DOS 累计 |
| `normalize_dos_high_energy()` | 高能尾均值归一化 |
| `plot_dos()` | 多磁场 DOS 对比图 |
| `plot_dos_single()` | 单磁场 DOS 图 |
| `save_dos_dat()` | 数据导出为文本文件 |

**main.py**：工作流调度器

- 读取配置 → 加载 hr.dat → 生成 k 网格 → 并行分派 → 合并 → 输出

### 4.3 数据流

```
hr.dat → hr_io.read_hr_dat_full()
           ├── nwannier (自动探测)
           ├── r_vectors (格矢列表)
           └── hr_data (稀疏 H(R) 矩阵)

config.py → 参数（Δ, μ, G_EFFECT, 网格, 磁场...）

kpoints.generate_mp_grid_weighted() → kpoints 数组

                        ▼
      k 点分片 → ProcessPoolExecutor
                        │
            ┌───────────┼───────────┐
            ▼           ▼           ▼
        Worker 0    Worker 1    Worker 2 ...
            │           │           │
        相位预计算    相位预计算    相位预计算
        H(k) 插值     H(k) 插值     H(k) 插值
        BdG 构建      BdG 构建      BdG 构建
        对角化        对角化        对角化
        DOS 累加      DOS 累加      DOS 累加
            │           │           │
            └───────────┼───────────┘
                        ▼
                 DOS 合并 (主进程)
                        ▼
             高能尾归一化（可选）
                        ▼
             dos.py → .dat + .png
```

### 4.4 并行策略

**第一层：磁场顺序处理**。各磁场值按列表顺序依次计算，结果显示进度条便于跟踪。

**第二层：k 点多进程并行**。对每个磁场，将 k 点均匀切分为 N_workers 个切片，通过 `ProcessPoolExecutor` 分配到各工作进程。每个进程将切片再分为 `BATCH_SIZE` 大小的子批次以控制内存，逐 k 点构建 BdG 矩阵并对角化，累加 DOS 后返回结果。

k 点之间计算天然独立，可达到近乎线性的加速比。子批次机制避免了大规模 k 网格下的内存溢出。

## 5. 使用指南

### 5.1 环境配置

**要求**：Python ≥ 3.9

**安装依赖**：

```bash
pip install numpy scipy matplotlib tqdm
```

### 5.2 输入文件准备

需要 Wannier90 输出的 `hr.dat` 文件。默认路径为 `../NbSe2_hr.dat`，可通过 `config.py` 中的 `HR_PATH` 参数修改。

hr.dat 格式说明：

- 第 1 行：注释
- 第 2 行：nwannier（Wannier 函数数，程序自动读取，无需用户指定）
- 第 3 行：nrpts（R 矢量数）
- 退简并度行
- 数据行：每行包含 `rx ry rz i j Re(H) Im(H)`

### 5.3 参数配置

编辑 `config.py` 设置计算参数。必须根据目标材料调整的参数：

| 参数 | 说明 | 调整依据 |
|------|------|----------|
| `DELTA` | 超导能隙 (meV) | 实验值或理论估算 |
| `FERMI_LEVEL` | 化学势 (eV) | DFT 计算结果 |
| `G_EFFECT` | 面内 Zeeman 增强因子 | 常规体系 = 1，Ising SOC 体系 > 1 |
| `K_MESH` | k 网格密度 | 根据倒格矢长度和收敛测试确定 |
| `B_FIELDS_Z/X` | 磁场列表 (T) | 根据研究需要设定 |
| `NORMALIZE_BACKGROUND` | 是否归一化（除以高能区均值） | `True` 开启，`False` 关闭 |
| `NORMALIZE_FRACTION` | 参考区间占正能侧比例 | 默认 0.1（尾部 10%） |

### 5.4 运行程序

```bash
cd bdg_dos/
python main.py
```

指定进程数：

```bash
python main.py --workers 20
```

不指定时默认使用 CPU 核心数减 1。

运行将依次：
1. 输出体系信息（nwannier, nrpts, BdG 矩阵维度）
2. 报告 k 网格信息
3. 依次处理面外 (z) 和面内 (x) 方向各磁场值
4. 每个磁场显示 tqdm 进度条

### 5.5 输出文件

```
output/
├── z/                        # 面外方向
│   ├── DOS_z_B0.0T.dat
│   ├── DOS_z_B0.0T.png
│   ├── DOS_z_B1.0T.dat
│   └── ...
├── x/                        # 面内方向
│   ├── DOS_x_B0.0T.dat
│   ├── DOS_x_B0.0T.png
│   └── ...
├── DOS_z_direction.png       # 所有 z 场叠加图
└── DOS_x_direction.png       # 所有 x 场叠加图
```

**数据文件格式**：

```text
# DOS calculation
# B_direction = z
# Delta = 0.2 meV, Fermi level = -2.0697 eV
# Energy(meV)    DOS_B=0T    DOS_B=1T    ...
  -0.600000      0.000123    0.000145    ...
  -0.599400      0.000130    0.000152    ...
```

各列：能量 (meV) | 各磁场下的 DOS (states/eV)。

## 6. 技术细节

### 6.1 Wannier90 hr.dat 格式解析

解析流程：

1. 读取文件第 2 行获取 nwannier
2. 读取第 3 行获取 nrpts
3. 跳过退简并度行（共 (nrpts + 4) // 5 行，每行最多 5 个整数）
4. 逐行解析数据，按 (rx, ry, rz) 键值分组
5. 对每个 R 构建 CSR 稀疏矩阵

```python
# 核心数据结构
hr_dict = {(rx, ry, rz): (rows_list, cols_list, data_list)}
hr_data = [csr_matrix((data, (rows, cols)), shape=(N, N)) for each R]
```

### 6.2 k 点分片与并行调度

分片保证负载均衡：

```python
chunk_size = (N_k + N_w - 1) // N_w
slices = [(i * chunk_size, min((i+1) * chunk_size, N_k))
          for i in range(N_w) if i * chunk_size < N_k]
```

每个工作进程通过 `_process_k_slice()` 函数独立处理分配的 k 点，内部进一步按 `BATCH_SIZE` 分批次。由于数据只读共享，无需进程间通信或加锁。

### 6.3 内存优化

- **子批次处理**：防止大规模相位因子矩阵一次性占用过多内存
- **及时释放**：每批次后用 `del` 释放中间变量
- **稀疏存储**：H(R) 以 CSR 格式保存，仅在求和时转为稠密
- **特征值模式**：使用 `eigvalsh` 避免计算特征向量

## 7. 常见问题与修复

### 7.1 多能带体系 V 型 DOS 问题

#### 问题描述

对于多能带材料（如 InNbSe₂），即使设置 s 波配对（各向同性超导能隙），DOS 在零能附近呈现 **V 形**（非零 DOS ~0.5 states/eV），而非理论预期的 **U 形**（能隙完全打开，零能 DOS = 0）。

#### 根本原因：时间反演对称性 (TRS) 数值破坏

BdG 哈密顿量的空穴块为：

```
hole_block = -conj[H(-k) + H_Z] + μ·I
```

传统实现中，`H(-k)` 通过对 Wannier90 相位因子取复共轭后重新进行傅里叶求和得到：

```python
H_minusk = compute_Hk_from_phases(hr_data, np.conj(phase_k), nwannier)
```

这种做法在**多能带体系**中会严重破坏时间反演对称性：

```
T H(k) T^{-1} = H(-k),  其中 T = -iσ_y·K  (自旋 1/2 时间反演算符)
```

数值验证显示 `H(-k)` 与 `iσ_y·H(k)^*·iσ_y` 的最大偏差可达 **1.8 eV**，远超超导能隙的 0.2 meV。这直接导致粒子块与空穴块的能带结构不匹配，超导能隙无法正确打开。

#### 单能带 vs 多能带差异

| | 单能带（如 NbSe₂） | 多能带（如 InNbSe₂） |
|--|--|--|
| Wannier 函数数 | ~22 (11 轨道 × 2 自旋) | ~28 (14 轨道 × 2 自旋) |
| 能带结构 | 仅 1 条带在 μ 附近 | 多条带在 μ 附近 |
| TRS 破坏程度 | 可忽略（ED 与 TRS 路径等价） | 严重（差值达 1.8 eV） |
| DOS 结果 | U 形（正确） | V 形（错误） |

单能带体系 `H(-k) ≈ H(k)` 近似成立，因此传统方法碰巧正确。多能带体系中不同轨道间存在复杂杂化，导致 `H(-k)` 的数值 TRS 破坏被急剧放大。

#### 修复方案：基于 TRS 的空穴块构造

TRS 确保：

```
H(-k) = -iσ_y · H(k)^* · iσ_y
conj[H(-k)] = -iσ_y · H(k) · iσ_y
```

因此空穴块可以绕过 `H(-k)` 的计算，直接由 `H(k)` 通过 TRS 变换得到：

```
hole_block = iσ_y · H(k) · iσ_y - conj(H_Z) + μ·I
```

#### 实现变更

在 `bdg.py` 的 `build_BdG_matrix_full()` 函数中：

```python
# 旧：通过独立傅里叶求和计算 H(-k)（数值 TRS 破坏）
Hmk_sym = (H_minusk + H_minusk.conj().T) / 2.0
hole_block = -np.conj(Hmk_sym + H_Z) + mu * np.eye(nwannier, dtype=complex)

# 新：利用 TRS 从 H(k) 直接构造空穴块（始终精确）
i_sigma_y_full = np.kron(i_sigma_y, I_norb)
hole_Hk = i_sigma_y_full @ Hk_sym @ i_sigma_y_full
hole_block = hole_Hk - np.conj(H_Z) + mu * np.eye(nwannier, dtype=complex)
```

注意此处 `i_sigma_y = [[0, 1], [-1, 0]]`，即 `i·σ_y` 矩阵。

#### 验证

- 在 InNbSe₂ 的 FS 最近 k 点（|E−μ| = 0.0006 meV）上，修复后 `min|E| = 0.200001 meV = Δ`，与理论一致
- Kramer 配对量 ⟨v|D|Tv⟩/Δ = −1.000000 对所有能带精确成立
- 对角配对的能带基 BdG 也给出正确能隙，验证了配对矩阵本身没有问题

#### 关键结论

**配对矩阵本身是正确的**：D_pair = Δ·(iσ_y)⊗I_orb，对所有能带的 Kramer 配对量都是精确的 −Δ。问题的根源在于**空穴块的构造方式破坏了 TRS**，而非配对矩阵的形式。修复后多能带体系也能正确呈现 U 形 DOS。

---

## 8. 参考文献

1. P. G. de Gennes, *Superconductivity of Metals and Alloys* (Benjamin, 1966)
2. A. A. Mostofi et al., *Comput. Phys. Commun.* **185**, 2309 (2014) — Wannier90
3. H. J. Monkhorst and J. D. Pack, *Phys. Rev. B* **13**, 5188 (1976)
4. N. Marzari et al., *Rev. Mod. Phys.* **84**, 1419 (2012)
