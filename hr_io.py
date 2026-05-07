"""Read Wannier90 hr.dat file and construct H(R) sparse matrices.

Wannier90 v1.x and v3.x both use block spin-orbital ordering:
    [up1, up2, ..., upN, down1, down2, ..., downN]
i.e., all spin-up Wannier functions first, then all spin-down.

The only difference is that v1.x (without SOC) may output separate
_up_hr.dat / _dn_hr.dat files, while v3.x (with SOC, spinors=True)
outputs a single combined _hr.dat file. The index ordering within
the combined matrix is the same for both versions.
"""

import numpy as np
from scipy import sparse


def read_hr_dat_full(filepath, version="v3"):
    """
    Read Wannier90 hr.dat file.

    Both v1 and v3 versions use the same block spin ordering
    (all up first, then all down), so the version parameter is
    for documentation / future extension only.

    Args:
        filepath: path to hr.dat file
        version: 'v1' or 'v3' (default). Currently treated identically.

    Returns:
        nwannier: number of Wannier functions
        r_vectors: list of R vector tuples (rx, ry, rz)
        hr_data: list of sparse matrices H(R), one per R vector
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Parse header
    nwannier = int(lines[1].strip())
    nrpts = int(lines[2].strip())

    # Find where hopping data starts
    degen_lines = (nrpts + 4) // 5
    data_start = 3 + degen_lines

    # Group hopping data by R vector
    hr_dict = {}
    for line in lines[data_start:]:
        parts = line.split()
        if len(parts) < 7:
            continue
        rx, ry, rz = int(parts[0]), int(parts[1]), int(parts[2])
        i = int(parts[3]) - 1  # 0-indexed
        j = int(parts[4]) - 1  # 0-indexed
        re_h = float(parts[5])
        im_h = float(parts[6])
        h_val = complex(re_h, im_h)

        key = (rx, ry, rz)
        if key not in hr_dict:
            hr_dict[key] = ([], [], [])
        hr_dict[key][0].append(i)
        hr_dict[key][1].append(j)
        hr_dict[key][2].append(h_val)

    # Build list of R vectors and sparse matrices
    r_vectors = list(hr_dict.keys())
    hr_data = []
    for R in r_vectors:
        rows, cols, data = hr_dict[R]
        H_R = sparse.csr_matrix(
            (data, (rows, cols)),
            shape=(nwannier, nwannier),
            dtype=complex
        )
        hr_data.append(H_R)

    return nwannier, r_vectors, hr_data
