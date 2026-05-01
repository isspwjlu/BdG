"""Monkhorst-Pack k-point grid generation for NbSe2 BdG DOS calculations."""

import numpy as np


def generate_mp_grid(mesh, shift=None):
    """
    Generate Monkhorst-Pack k-point grid.

    Args:
        mesh: [nkx, nky, nkz] - number of k-points along each direction
        shift: optional [sx, sy, sz] - shift from origin (usually 0 or 0.5)

    Returns:
        kpoints: (N, 3) numpy array of k-point coordinates in reciprocal lattice units
    """
    if shift is None:
        shift = [0.0, 0.0, 0.0]

    kpoints = []
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            for k in range(mesh[2]):
                kx = (i + 0.5) / mesh[0] - 0.5 * shift[0]
                ky = (j + 0.5) / mesh[1] - 0.5 * shift[1]
                kz = (k + 0.5) / mesh[2] - 0.5 * shift[2]
                kpoints.append([kx, ky, kz])

    return np.array(kpoints)


def generate_mp_grid_weighted(mesh, shift=None):
    """
    Generate Monkhorst-Pack k-point grid with weights.

    Args:
        mesh: [nkx, nky, nkz] - number of k-points along each direction
        shift: optional [sx, sy, sz] - shift from origin (usually 0 or 0.5)

    Returns:
        kpoints: (N, 3) array of k-point coordinates
        weights: (N,) array of weights (all 1/Nk)
    """
    kpoints = generate_mp_grid(mesh, shift)
    nk = len(kpoints)
    weights = np.ones(nk) / nk
    return kpoints, weights


def kpoint_batches(kpoints, batch_size):
    """
    Split kpoints into batches for memory-efficient processing.

    Args:
        kpoints: (N, 3) array of k-point coordinates
        batch_size: number of k-points per batch

    Yields:
        batches of k-point indices (numpy arrays)
    """
    nk = len(kpoints)
    indices = np.arange(nk)

    for start in range(0, nk, batch_size):
        end = min(start + batch_size, nk)
        yield indices[start:end]
