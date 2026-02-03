"""Visualize a triangular lattice with bond coloring and failed nodes.

This script loads lattice positions and bond data and produces a publication-ready
figure with a consistent font and clean axes.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc

from read_output import *


# ----- Publication style defaults -----
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Times New Roman"
plt.rcParams["mathtext.it"] = "Times New Roman:italic"
plt.rcParams["mathtext.bf"] = "Times New Roman:bold"


# ----- Configuration -----
UC_SIZE = 4.0 # Unit cell size in plotting units.
DATA_DIR = Path("../src/IO/triang/")
RUN_ID = None  # Set to a specific run name if needed.
FIGSIZE = (6, 3)
NX, NY = 60, 30
BOND_LINEWIDTH = 1.0
CRACK_NODE_SIZE = 20


def _normalize(arr: np.ndarray) -> np.ndarray:
    """Normalize an array safely to [0, 1]."""
    max_val = np.max(arr)
    if max_val == 0:
        return np.zeros_like(arr)
    return arr / max_val


def plot_frame(
    positions: np.ndarray,
    bonds: np.ndarray,
    n_init: int,
    *,
    nx: int = NX,
    ny: int = NY,
    uc_size: float = UC_SIZE,
    figsize: tuple[float, float] = FIGSIZE,
) -> None:
    """Plot a single frame of the lattice."""
    bond_ends = [
        [
            (positions[int(bond[0])][0], positions[int(bond[0])][1]),
            (positions[int(bond[1])][0], positions[int(bond[1])][1]),
        ]
        for bond in bonds
    ]

    crack_node_pos = positions[n_init:]

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-uc_size, nx * uc_size + uc_size)
    ax.set_ylim(-uc_size, ny * uc_size + uc_size)
    ax.set_aspect("equal")
    ax.tick_params(axis="both", which="both", length=0, labelbottom=False, labelleft=False)

    colorray = _normalize(np.abs(bonds[:, 2]) + np.abs(bonds[:, 4]))
    line_collection = mc.LineCollection(bond_ends, array=colorray, cmap=plt.cm.Reds, linewidths=2)
    ax.add_collection(line_collection)
    line_collection.set_linewidth(BOND_LINEWIDTH)

    # Put failed nodes on top (make them visible if inside the lattice)
    ax.scatter(crack_node_pos[:, 0], crack_node_pos[:, 1], c="r", zorder=3, s=CRACK_NODE_SIZE)

    plt.tight_layout()
    plt.show()


def main() -> None:
    """Entry point for lattice visualization."""
    P_array, n_init = read_pos(str(DATA_DIR), RUN_ID)
    bondframes = read_bonds(str(DATA_DIR), RUN_ID)

    # Boundary indices are loaded for completeness; unused in this visualization.
    _boundary_inds_bot, _boundary_inds_top = read_boundary_indeces(str(DATA_DIR))

    frame_no = np.arange(1, len(bondframes), 2)
    frame_no = np.append(frame_no, len(bondframes) - 1)

    for frame in frame_no[-1:]:
        posxy = P_array[frame]
        ind_bonds = bondframes[frame]
        plot_frame(posxy, ind_bonds, n_init)


if __name__ == "__main__":
    main()