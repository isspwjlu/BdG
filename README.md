# Superconducting density of states (DOS) calculation by solving BdG equation under magnetic field

Density of states from a full Wannier-basis Bogoliubov–de Gennes calculation, with out-of-plane and in-plane Zeeman fields.

[![Python](https://img.shields.io/badge/python-%3E%3D3.9-blue)](https://python.org)
[![NumPy](https://img.shields.io/badge/numpy-%3E%3D1.21-blueviolet)](https://numpy.org)
[![SciPy](https://img.shields.io/badge/scipy-%3E%3D1.7-green)](https://scipy.org)
[![License](https://img.shields.io/badge/license-MIT-lightgrey)](LICENSE)

## Overview

This code computes the superconducting density of states of a crystalline material by diagonalising a Bogoliubov–de Gennes (BdG) Hamiltonian built from a Wannier90 tight-binding model. The Zeeman coupling includes a configurable enhancement factor for in-plane fields, suitable for modelling systems with spin–orbit coupling.

The BdG matrix is constructed in the full Wannier orbital–spin basis — **2N × 2N** in Nambu space, where N is the number of Wannier functions — and diagonalised on a dense Monkhorst–Pack k-grid. The density of states is accumulated via Lorentzian broadening and exported as data files and plots.

## Features

- Full Wannier-orbital BdG formalism (no orbital-symmetry reduction)
- Out-of-plane Zeeman (standard g = 2) and in-plane Zeeman with tunable enhancement factor
- Parallel execution over k-points via Python multiprocessing
- Monkhorst–Pack grid with configurable mesh density
- Multi-field DOS curves overlaid on a single plot, plus per-field output files

## Theoretical background

The BdG Hamiltonian in Nambu space takes the form

```
H_BdG(k) = ⎡ A(k)       D   ⎤
           ⎣ -D    -conj(A(-k)) ⎦
```

where:

- **Normal-state block:** `A(k) = H(k) + H_Z − μ·I`.`H(k) = Σ_R H(R) exp(2πi·k·R)` is the Wannier-interpolated Hamiltonian from Wannier90.
- **Pairing block:** `D = Δ · (iσ_y) ⊗ I_orb` describes s-wave singlet pairing.
- **Zeeman coupling:**

  ```
  H_Z = (g·μ_B/2) · [σ_x·(g_eff·B_x) + σ_y·(g_eff·B_y) + σ_z·B_z] ⊗ I_orb
  ```

  Out-of-plane (Bz) uses the bare electron g-factor (g = 2). In-plane (Bx, By) is scaled by the configurable parameter `G_EFFECT` to model enhancement from spin–orbit coupling. `G_EFFECT = 1` corresponds to the bare Zeeman term; values > 1 represent effective enhancement (e.g. from Ising SOC in monolayer superconductors).

## Project structure

```
bdg_dos/
├── __init__.py       # Package marker
├── config.py         # Physical and numerical parameters
├── hr_io.py          # Wannier90 hr.dat reader
├── kpoints.py        # Monkhorst–Pack k-point grid generator
├── bdg.py            # BdG Hamiltonian: construction, Zeeman, diagonalisation
├── dos.py            # DOS accumulation, plotting, data export
├── main.py           # Entry point — parallel execution driver
└── README.md         # This file
```

| Module         | Purpose                                                                                                                                                                                                                                                                                 |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `config.py`  | All tunable parameters: superconducting gap Δ, chemical potential μ, G_EFFECT, k-mesh density, field lists, Lorentzian broadening η, energy range, and batch size for parallelisation.                                                                                               |
| `hr_io.py`   | Parses the Wannier90 `hr.dat` file. Auto-detects the number of Wannier functions N from the file header. Returns sparse H(R) matrices and the list of R vectors.                                                                                                                      |
| `kpoints.py` | Generates Monkhorst–Pack k-point grids with optional shifts and uniform weights.                                                                                                                                                                                                       |
| `bdg.py`     | Core engine:`precompute_phases()` computes exp(2πi·k·R); `compute_Hk_from_phases()` performs the Fourier interpolation to obtain H(k); `build_BdG_matrix_full()` assembles the 2N×2N BdG matrix; `diagonalize_BdG()` solves the eigenproblem via `scipy.linalg.eigvalsh`. |
| `dos.py`     | Accumulates the density of states with Lorentzian broadening, produces plots (single-field and multi-field) via matplotlib, and exports data in plain-text `.dat` format.                                                                                                             |
| `main.py`    | Orchestrates the full workflow: reads `hr.dat`, builds the k-mesh, dispatches slices of k-points across a `ProcessPoolExecutor`, and drives the calculation for both out-of-plane (z) and in-plane (x) field directions.                                                            |

## Getting started

### Prerequisites

- Python ≥ 3.9
- [NumPy](https://numpy.org), [SciPy](https://scipy.org), [matplotlib](https://matplotlib.org), [tqdm](https://tqdm.github.io)

Install with:

```bash
pip install numpy scipy matplotlib tqdm
```

### Running the code

Place your Wannier90 `hr.dat` file one directory above the `bdg_dos/` folder (or adjust the path in `main.py`), then run:

```bash
cd bdg_dos/
python main.py
```

Specify the number of worker processes (defaults to CPU count − 1):

```bash
python main.py --workers 20
```

### What happens

1. N (number of Wannier functions) is auto-detected from `hr.dat`
2. A Monkhorst–Pack k-grid is generated (configurable via `K_MESH` in `config.py`)
3. For each magnetic field in `B_FIELDS_Z` and `B_FIELDS_X`, the BdG Hamiltonian is diagonalised across all k-points in parallel
4. DOS curves are accumulated, plotted, and saved to the `output/` directory

## Configuration

Edit `config.py` to adjust calculation parameters:

| Parameter       | Default                   | Description                                  |
| --------------- | ------------------------- | -------------------------------------------- |
| `DELTA`       | 0.2 meV                   | Superconducting gap                          |
| `FERMI_LEVEL` | −2.0697 eV               | Chemical potential μ (from DFT)             |
| `G_FACTOR`    | 2.0                       | Electron spin g-factor                       |
| `G_EFFECT`    | 1                       | In-plane Zeeman enhancement factor           |
| `ETA`         | 0.05·Δ                  | Lorentzian broadening width                  |
| `E_RANGE`     | 3·Δ                     | Energy half-range for DOS                    |
| `N_ENERGY`    | 2000                      | Number of energy grid points (negative half) |
| `K_MESH`      | [1000, 1000, 1]           | k-point mesh along [kx, ky, kz]              |
| `B_FIELDS_Z`  | [0, 1, 2, 3, 4, 5] T      | Magnetic field values, out-of-plane          |
| `B_FIELDS_X`  | [0, 10, 20, 30, 40, 50] T | Magnetic field values, in-plane              |
| `BATCH_SIZE`  | 100000                    | K-points per sub-batch per worker process    |

## Output

```
output/
├── z/                              # Per-field results — out-of-plane
│   ├── DOS_z_B0.0T.dat
│   ├── DOS_z_B0.0T.png
│   ├── DOS_z_B1.0T.dat
│   └── ...
├── x/                              # Per-field results — in-plane
│   ├── DOS_x_B0.0T.dat
│   ├── DOS_x_B0.0T.png
│   └── ...
├── DOS_z_direction.png            # Combined plot, all z-fields
└── DOS_x_direction.png            # Combined plot, all x-fields
```

### Data format

`.dat` files are tab-separated plain text:

```text
# DOS calculation
# B_direction = z
# Delta = 0.2 meV, Fermi level = -2.0697 eV
# Energy(meV)    DOS_B=0T    DOS_B=1T    DOS_B=2T
  -0.600000      0.000123    0.000145    0.000167
  -0.599400      0.000130    0.000152    0.000175
  ...
```

- Column 1: energy (meV) relative to the Fermi level
- Columns 2+: density of states (states/eV) at each magnetic field

## License

This project is distributed under the MIT License. See `LICENSE` for details.
