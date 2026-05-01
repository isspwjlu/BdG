"""BdG Hamiltonian construction and diagonalization for NbSe2.

The Wannier Hamiltonian H(k) from Wannier90 already includes orbital and spin
degrees of freedom. The BdG matrix is built in the FULL Wannier orbital-spin
basis (nwannier × nwannier) and then doubled to (2·nwannier × 2·nwannier)
in Nambu space.

Only the full BdG method (44×44) is used. Fermi-surface projection and
orbital Peierls methods have been removed.

    H_BdG(k) = [ A(k)       D      ]
               [ -D     -conj(A(-k)) ]

    where:
    - A(k) = H(k) + H_Z - μ·I       (electron sector)
    - D = Δ·(iσ_y)⊗I_orb            (s-wave singlet pairing)
    - H_Z = (g·μ_B/2)·σ·B ⊗ I_orb   (Zeeman, with in-plane g_effect factor)
"""

import numpy as np
from scipy import linalg

# Pauli matrices
SIGMA0 = np.eye(2, dtype=complex)
SIGMAX = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMAY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMAZ = np.array([[1, 0], [0, -1]], dtype=complex)

# Physical constants
MU_B = 5.788e-5  # eV/T, Bohr magneton
G_FACTOR = 2.0   # electron g-factor


def precompute_phases(kpoints, r_vectors):
    """
    Precompute exp(2πi·k·R) for all k-points and R vectors.

    Args:
        kpoints: (nk, 3) k-points in reciprocal lattice units
        r_vectors: list of R-vector integer triples

    Returns:
        phases: (nk, nr) complex array exp(2πi·k_i·R_j)
    """
    nk = len(kpoints)
    nr = len(r_vectors)
    phases = np.zeros((nk, nr), dtype=complex)
    k_arr = np.array(kpoints)
    for i, R in enumerate(r_vectors):
        phases[:, i] = np.exp(2j * np.pi * (
            k_arr[:, 0] * R[0] +
            k_arr[:, 1] * R[1] +
            k_arr[:, 2] * R[2]))
    return phases


def compute_Hk_from_phases(hr_data, phases_k, nwannier):
    """
    Compute H(k) = Σ_R H(R) exp(2πi·k·R).

    Args:
        hr_data: list of sparse H(R) matrices
        phases_k: (nr,) complex phase factors for one k-point
        nwannier: number of Wannier functions

    Returns:
        Hk: (nwannier, nwannier) dense complex H(k)
    """
    Hk = np.zeros((nwannier, nwannier), dtype=complex)
    for H_R, phase in zip(hr_data, phases_k):
        Hk += H_R.toarray() * phase
    return Hk


def build_BdG_matrix_full(Hk, delta_eV, B_vec, nwannier, mu,
                          H_minusk=None, g_effect=1.0):
    """
    Build the full (2·nwannier)-dimensional BdG matrix.

    Args:
        Hk: (nwannier, nwannier) H(k)
        delta_eV: superconducting gap (eV)
        B_vec: (3,) field (Bx, By, Bz) in Tesla
        nwannier: number of Wannier functions
        mu: chemical potential (eV)
        H_minusk: H(-k); if None uses conj(H(k))
        g_effect: scaling for in-plane (Bx, By) Zeeman

    Returns:
        H_bdg: (2·nwannier, 2·nwannier) Hermitian BdG matrix
    """
    norb = nwannier // 2
    I_norb = np.eye(norb, dtype=complex)

    # Zeeman: in-plane scaled by g_effect (Ising SOC suppression)
    Bx, By, Bz = B_vec
    sigma_B = (SIGMAX * (g_effect * Bx) +
               SIGMAY * (g_effect * By) +
               SIGMAZ * Bz)
    H_Z = (0.5 * G_FACTOR * MU_B) * np.kron(sigma_B, I_norb)

    # Pairing: Δ·(iσ_y)⊗I_orb  (s-wave singlet)
    i_sigma_y = np.array([[0, 1], [-1, 0]], dtype=complex)
    D_pair = delta_eV * np.kron(i_sigma_y, I_norb)

    # Hermiticity symmetrisation
    Hk_sym = (Hk + Hk.conj().T) / 2.0

    # Electron block
    A = Hk_sym + H_Z - mu * np.eye(nwannier, dtype=complex)

    # Hole block
    if H_minusk is not None:
        Hmk_sym = (H_minusk + H_minusk.conj().T) / 2.0
        hole_block = -np.conj(Hmk_sym + H_Z) + mu * np.eye(nwannier, dtype=complex)
    else:
        hole_block = -np.conj(A)

    dim = 2 * nwannier
    H_bdg = np.zeros((dim, dim), dtype=complex)
    H_bdg[:nwannier, :nwannier] = A
    H_bdg[:nwannier, nwannier:] = D_pair
    H_bdg[nwannier:, :nwannier] = -D_pair
    H_bdg[nwannier:, nwannier:] = hole_block
    return H_bdg


def diagonalize_BdG(H_bdg):
    """
    Diagonalise BdG matrix and return eigenvalues.

    Args:
        H_bdg: (2·nwannier, 2·nwannier) Hermitian matrix, or None

    Returns:
        eigenvalues: (2·nwannier,) array in eV, or empty array
    """
    if H_bdg is None:
        return np.array([])
    return linalg.eigvalsh(H_bdg)
