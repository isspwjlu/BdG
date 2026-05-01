"""
Density of States (DOS) calculation and plotting for NbSe2 BdG Hamiltonian.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt


def compute_dos_grid(eigenvalues, energy_grid, eta):
    """
    Compute DOS on an energy grid from eigenvalues using Lorentzian broadening.

    Args:
        eigenvalues: (Nk, 2*nwannier) array of all eigenvalues in eV
        energy_grid: (Ne,) array of energy values in eV
        eta: broadening parameter in eV

    Returns:
        dos: (Ne,) array of DOS values

    Formula:
        N(E) = (1/Nk) * sum_{k,n} eta / ((E - E_{k,n})^2 + eta^2)
    """
    Nk = eigenvalues.shape[0]
    dos = np.zeros(len(energy_grid))

    for k in range(Nk):
        for band in range(eigenvalues.shape[1]):
            e = eigenvalues[k, band]
            diff = energy_grid - e
            lorentz = eta / (diff**2 + eta**2)
            dos += lorentz

    dos /= (np.pi * Nk)
    return dos


def compute_dos_vectorized(eigenvalues, energy_grid, eta):
    """
    Vectorized DOS calculation using broadcasting.
    eigenvalues: (Nk, 2*nwannier) in eV
    energy_grid: (Ne,) in eV
    eta: in eV
    Returns: (Ne,) DOS
    """
    # Broadcasting: eigenvalues[:, :, None] - energy_grid[None, None, :]
    diff = eigenvalues[:, :, None] - energy_grid[None, None, :]
    lorentz = eta / (diff**2 + eta**2)
    dos = np.sum(lorentz, axis=(0, 1)) / (np.pi * eigenvalues.shape[0])
    return dos


def normalize_dos_high_energy(dos_full, fraction=0.1):
    """
    Normalize DOS by the average value in the high-energy tail.

    The high-energy tail of the DOS should approach a constant (normal-state)
    value. Due to truncation of Lorentzian tails at finite energy range,
    different magnetic fields can exhibit different baselines. This function
    divides the entire DOS by the average in the high-energy tail, so that
    the high-energy region is normalized to ~1 for all B-fields.

    Args:
        dos_full: (Ne,) full symmetric DOS array (negative then positive E).
        fraction: fraction of the positive-energy side used as background.

    Returns:
        dos_normalized: (Ne,) normalized DOS (dimensionless).
    """
    n_pos = len(dos_full) // 2
    n_bg = max(1, int(n_pos * fraction))
    bg = np.mean(dos_full[-n_bg:])
    if bg <= 0:
        return dos_full  # fallback: skip if background is non-positive
    return dos_full / bg


def plot_dos(energy_grid, dos_data, B_fields, direction, output_path, normalized=False):
    """
    Plot DOS vs energy for different magnetic fields.

    Args:
        energy_grid: (Ne,) energy values in meV
        dos_data: dict mapping B_value -> dos array
        B_fields: list of B values
        direction: 'z' or 'x'
        output_path: path to save figure
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Color map for different B fields
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(B_fields)))

    for i, B in enumerate(B_fields):
        dos = dos_data[B]
        ax.plot(energy_grid, dos, color=colors[i], label=f'B = {B} T',
                linewidth=1.5)

    ax.set_xlabel('Energy (meV)', fontsize=12)
    ax.set_ylabel('DOS (states/eV)', fontsize=12)
    label = 'NbSe$_2$ BdG DOS'
    if normalized:
        label += ' (normalized)'
    ax.set_title(f'{label} - {direction}-direction', fontsize=14)
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([energy_grid.min(), energy_grid.max()])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_dos_single(energy_grid, dos, B, direction, output_path, normalized=False):
    """
    Plot DOS vs energy for a single magnetic field.

    Args:
        energy_grid: (Ne,) energy values in meV
        dos: (Ne,) DOS array
        B: magnetic field value (float)
        direction: 'z' or 'x'
        output_path: path to save figure
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(energy_grid, dos, color='navy', linewidth=1.5)

    ax.set_xlabel('Energy (meV)', fontsize=12)
    ax.set_ylabel('DOS (states/eV)', fontsize=12)
    label = 'NbSe$_2$ BdG DOS'
    if normalized:
        label += ' (normalized)'
    ax.set_title(f'{label} - B = {B} T ({direction}-direction)', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([energy_grid.min(), energy_grid.max()])

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def save_dos_dat(energy_grid, dos_data, B_fields, direction, output_path, normalized=False):
    """
    Save DOS data to .dat file.

    Format:
        # NbSe2 DOS calculation
        # B_direction = z
        # Delta = 0.3 meV, Fermi level = -2.057 eV
        # Energy(meV)    DOS_B=0T    DOS_B=2T    ...
    """
    from config import DELTA, FERMI_LEVEL

    with open(output_path, 'w') as f:
        f.write('# NbSe2 DOS calculation\n')
        f.write(f'# B_direction = {direction}\n')
        f.write(f'# Delta = {DELTA} meV, Fermi level = {FERMI_LEVEL} eV\n')
        if normalized:
            f.write('# Normalized by high-energy tail average\n')
        f.write('# Energy(meV)\t' + '\t'.join([f'DOS_B={B}T' for B in B_fields]) + '\n')

        for i in range(len(energy_grid)):
            line = f'{energy_grid[i]:12.6f}'
            for B in B_fields:
                line += f'\t{dos_data[B][i]:12.6f}'
            f.write(line + '\n')
