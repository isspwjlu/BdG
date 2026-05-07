"""Configuration for NbSe2 BdG DOS — Full Wannier-basis BdG (44×44 matrix).

Units: eV (energy), Tesla (magnetic field), unless noted.
"""
import numpy as np


# Wannier90 hr.dat file path
HR_PATH = "../NbSe2_hr.dat"

# Output
OUTPUT_DIR = "output"

# Physical constants
MU_B = 5.788e-5       # eV/T, Bohr magneton
K_B = 8.617e-5        # eV/K, Boltzmann constant

# Material parameters (NbSe2 monolayer)
# NWANNIER is auto-detected from hr.dat header — not set here.
DELTA = 0.2           # meV, superconducting gap
FERMI_LEVEL = -2.0697 # eV, chemical potential μ from DFT
G_FACTOR = 2.0        # electron spin g-factor
TEMP = 20.0           # K (for plot labels only)

# Effective in-plane Zeeman enhancement factor (Ising SOC regime).
# Ising SOC (~100 meV) in NbSe₂ locks spins out-of-plane, which modifies
# the effective in-plane g-factor. G_EFFECT captures this enhancement:
#   H_Z_x = G_EFFECT * (g·μ_B/2) · σ_x · Bx
#   = 1   → standard Zeeman (bare electron g=2)
#   > 1  → effective enhancement due to Ising SOC mixing
# Default 100 corresponds to a 100× effective enhancement of in-plane
# Zeeman, consistent with the large in-plane critical field observed
# in monolayer NbSe₂ (B_c∥ ~ 30 T).
G_EFFECT = 1

# DOS calculation
ETA = 0.05 * DELTA    # meV, Lorentzian broadening
E_RANGE = 3 * DELTA   # meV, energy half-range (3×gap = 0.6 meV)
N_ENERGY = 2000       # energy grid points (negative half)

# k-point mesh  [nkx, nky, nkz]
# Mesh density should be determined from the reciprocal lattice vectors.
K_MESH = [1000, 1000, 1]

# Magnetic fields (Tesla)
B_FIELDS_Z = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
B_FIELDS_X = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]

# DOS normalization by high-energy tail average.
# When enabled, the entire DOS is divided by the average value in the
# high-energy tail (top NORMALIZE_FRACTION of the energy range), so that
# the high-energy DOS is normalized to ~1 for all B-fields. This corrects
# the vertical baseline shift caused by truncation of Lorentzian tails at
# finite energy range.
NORMALIZE_BACKGROUND = True
NORMALIZE_FRACTION = 0.1   # fraction of positive-energy side used as reference

# Performance
BATCH_SIZE = 100000   # k-points per sub-batch in worker


