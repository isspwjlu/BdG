"""Main program for NbSe2 BdG DOS — Full Wannier-basis BdG (44×44).

The BdG Hamiltonian is constructed in the full Wannier orbital-spin basis
(2·nwannier = 44-dimensional Nambu space). Zeeman term in the Wannier basis:

    H_Z = (g·μ_B/2)·σ·B ⊗ I_orb

In-plane Zeeman (Bx, By) is scaled by G_EFFECT (from config): G_EFFECT=1
is standard Zeeman, G_EFFECT>1 is Ising SOC enhancement.
"""

import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from tqdm import tqdm

from config import (
    DELTA, FERMI_LEVEL, G_EFFECT, ETA,
    E_RANGE, N_ENERGY, K_MESH, B_FIELDS_Z, B_FIELDS_X, BATCH_SIZE, OUTPUT_DIR,
    NORMALIZE_BACKGROUND, NORMALIZE_FRACTION,
)
from hr_io import read_hr_dat_full
from kpoints import generate_mp_grid_weighted
from bdg import (
    compute_Hk_from_phases, precompute_phases,
    build_BdG_matrix_full, diagonalize_BdG,
)
from dos import plot_dos, plot_dos_single, save_dos_dat, normalize_dos_high_energy


def _process_k_slice(args):
    """Process a k-point slice: build BdG, diagonalise, accumulate DOS."""
    (k_start, k_end, kpoints, r_vectors, hr_data, nwannier,
     delta_eV, B_vec, mu, energy_grid, eta_eV, nk) = args

    n_energy = len(energy_grid)
    dos_accum = np.zeros(n_energy, dtype=np.float64)
    k_slice = kpoints[k_start:k_end]

    sub_batch_size = BATCH_SIZE
    n_sub = (len(k_slice) + sub_batch_size - 1) // sub_batch_size

    for sub_idx in range(n_sub):
        sb_start = sub_idx * sub_batch_size
        sb_end = min((sub_idx + 1) * sub_batch_size, len(k_slice))
        k_sub = k_slice[sb_start:sb_end]

        phases_sub = precompute_phases(k_sub, r_vectors)

        for i in range(sb_end - sb_start):
            Hk = compute_Hk_from_phases(hr_data, phases_sub[i], nwannier)
            H_minusk = compute_Hk_from_phases(
                hr_data, np.conj(phases_sub[i]), nwannier)

            H_bdg = build_BdG_matrix_full(
                Hk, delta_eV, B_vec, nwannier, mu,
                H_minusk=H_minusk, g_effect=G_EFFECT)

            eigs = diagonalize_BdG(H_bdg)

            diff = eigs[:, None] - energy_grid[None, :]
            lorentz = eta_eV / (diff**2 + eta_eV**2)
            dos_accum += np.sum(lorentz, axis=0)

        del phases_sub

    return dos_accum


def _run_bdg_dos_direction(kpoints, k_slices, n_slices, hr_data, r_vectors,
                            nwannier, delta_eV, energy_grid, eta_eV, nk,
                            B_list, direction, output_dir, executor):
    """Process all B-fields for one direction."""
    results = {}

    for B in B_list:
        print(f"\n=== Processing B = {B} T ({direction}-direction) ===")

        B_vec = np.array([B, 0.0, 0.0] if direction == 'x'
                         else [0.0, 0.0, B])

        slice_args = [
            (k_start, k_end, kpoints, r_vectors, hr_data, nwannier,
             delta_eV, B_vec, FERMI_LEVEL, energy_grid, eta_eV, nk)
            for k_start, k_end in k_slices
        ]

        dos_accum = np.zeros(len(energy_grid))

        for ds in tqdm(executor.map(_process_k_slice, slice_args),
                       total=n_slices, desc=f"B={B}T"):
            dos_accum += ds

        dos_final = dos_accum / (np.pi * nk)
        energy_full = np.concatenate([energy_grid, -energy_grid[::-1]])
        dos_full = np.concatenate([dos_final, dos_final[::-1]])

        if NORMALIZE_BACKGROUND:
            dos_full = normalize_dos_high_energy(dos_full, fraction=NORMALIZE_FRACTION)

        results[B] = dos_full

        dos_file = os.path.join(output_dir, f'NbSe2_DOS_{direction}_B{B:.1f}T.dat')
        save_dos_dat(energy_full * 1000, {B: dos_full}, [B], direction,
                     dos_file, normalized=NORMALIZE_BACKGROUND)
        plot_file = os.path.join(output_dir, f'NbSe2_DOS_{direction}_B{B:.1f}T.png')
        plot_dos_single(energy_full * 1000, dos_full, B, direction, plot_file,
                        normalized=NORMALIZE_BACKGROUND)

    return energy_full * 1000, results


def run_bdg_dos_parallel(hr_path, k_mesh, B_list, direction, output_dir,
                         n_workers=None, _executor=None):
    """Parallel full BdG DOS calculation.

    Args:
        hr_path: path to Wannier90 hr.dat
        k_mesh: [nkx, nky, nkz] mesh
        B_list: field values (Tesla)
        direction: 'z' or 'x'
        output_dir: output directory
        n_workers: number of processes (auto if None)
        _executor: shared ProcessPoolExecutor
    """
    if n_workers is None:
        n_workers = max(1, multiprocessing.cpu_count() - 1)

    os.makedirs(output_dir, exist_ok=True)

    delta_eV = DELTA / 1000.0
    eta_eV = ETA / 1000.0
    e_range_eV = E_RANGE / 1000.0
    energy_grid = np.linspace(-e_range_eV, 0, N_ENERGY)

    print(f"Loading {hr_path}...")
    nwannier, r_vectors, hr_data = read_hr_dat_full(hr_path)
    print(f"nwannier={nwannier}, nrpts={len(r_vectors)}")
    print(f"BdG matrix: {2 * nwannier}×{2 * nwannier}")
    print(f"G_EFFECT = {G_EFFECT} (in-plane Zeeman scaling)")

    print(f"Generating {k_mesh[0]}x{k_mesh[1]}x{k_mesh[2]} k-mesh...")
    kpoints, _ = generate_mp_grid_weighted(k_mesh)
    nk = len(kpoints)
    print(f"Total k-points: {nk}")

    chunk_size = (nk + n_workers - 1) // n_workers
    k_slices = []
    for w in range(n_workers):
        ks = w * chunk_size
        ke = min((w + 1) * chunk_size, nk)
        if ks < nk:
            k_slices.append((ks, ke))
    print(f"{len(k_slices)} slices of ~{chunk_size} k-points, {n_workers} workers")

    own_pool = False
    if _executor is None:
        _executor = ProcessPoolExecutor(max_workers=n_workers)
        own_pool = True

    try:
        energy_full, results = _run_bdg_dos_direction(
            kpoints, k_slices, len(k_slices), hr_data, r_vectors, nwannier,
            delta_eV, energy_grid, eta_eV, nk,
            B_list, direction, output_dir, _executor)
    finally:
        if own_pool:
            _executor.shutdown()

    return energy_full, results


def main(use_parallel=True, n_workers=None):
    """Main entry point."""
    hr_path = '../NbSe2_hr.dat'
    output_dir = OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("NbSe2 BdG DOS — Full Wannier-basis BdG")
    nw = n_workers or (multiprocessing.cpu_count() - 1)
    if use_parallel:
        print(f"Mode: parallel, {nw} workers")
    print(f"G_EFFECT = {G_EFFECT} (in-plane Zeeman enhancement factor)")
    print(f"K_MESH = {K_MESH}")
    print("=" * 60)

    if use_parallel:
        with ProcessPoolExecutor(max_workers=nw) as executor:
            print("\n>>> z-direction magnetic fields")
            energy_z, results_z = run_bdg_dos_parallel(
                hr_path, K_MESH, B_FIELDS_Z, 'z',
                os.path.join(output_dir, 'z'),
                n_workers=nw, _executor=executor)

            print("\n>>> x-direction magnetic fields")
            energy_x, results_x = run_bdg_dos_parallel(
                hr_path, K_MESH, B_FIELDS_X, 'x',
                os.path.join(output_dir, 'x'),
                n_workers=nw, _executor=executor)
    else:
        energy_z, results_z = run_bdg_dos_parallel(
            hr_path, K_MESH, B_FIELDS_Z, 'z',
            os.path.join(output_dir, 'z'), n_workers=1)
        energy_x, results_x = run_bdg_dos_parallel(
            hr_path, K_MESH, B_FIELDS_X, 'x',
            os.path.join(output_dir, 'x'), n_workers=1)

    plot_dos(energy_z, results_z, B_FIELDS_Z, 'z',
             os.path.join(output_dir, 'NbSe2_DOS_z_direction.png'),
             normalized=NORMALIZE_BACKGROUND)
    plot_dos(energy_x, results_x, B_FIELDS_X, 'x',
             os.path.join(output_dir, 'NbSe2_DOS_x_direction.png'),
             normalized=NORMALIZE_BACKGROUND)

    print("\nDone!")


if __name__ == '__main__':
    import sys
    nw = None
    if len(sys.argv) > 1:
        for i, a in enumerate(sys.argv[1:]):
            if a in ('-w', '--workers') and i + 2 < len(sys.argv):
                nw = int(sys.argv[i + 2])
    main(n_workers=nw)
