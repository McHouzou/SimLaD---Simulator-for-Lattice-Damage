"""Create input files for a triangular lattice simulation.

This script builds a triangular lattice, applies optional disorder/cracks, and
writes node, bond, bend, and simulation input files in a reproducible format.
"""

from __future__ import annotations

import argparse
from typing import Tuple

import numpy as np

# TODO: Add new boundary-condition/flag functionality in `make_bnfiles` (see kagome input creator).

# ---- Output formatting ----
DECIMAL_TOL = 5  # number of decimals in output files

def dist_nodes(ra: np.ndarray, rb: np.ndarray) -> float:
    """
    Return the Euclidean distance between two nodes.

    Parameters
    ----------
    ra, rb
        Position vectors of the two nodes.
    """
    rab = rb - ra
    return np.sqrt(np.dot(rab, rab))


def bond_CoM(ra: np.ndarray, rb: np.ndarray) -> np.ndarray:
    """
    Return the midpoint (center of mass) between two nodes.
    """
    return 0.5 * (ra + rb)



class Triangular_Lattice:

    def __init__(self, dir: str, Nx: int, Ny: int) -> None:
        """
        Parameters
        ----------
        dir
            Output directory for input files.
        Nx, Ny
            Number of unit cells in the x and y direction.
        """
        self.dir = dir
        self.Nx = Nx
        self.Ny = Ny
        self.Nnodes = self.Nx * self.Ny  # number of unit cells


    def makeNodes_AddBondInd(self) -> None:
        """
        Construct node coordinates and bond connectivity arrays.
        """
        self.rnodes = []
        self.ind_bonds = []

        # Basis vectors for the triangular lattice (two lattice directions)
        self.a = UC_SIZE * np.array([[1, 0], [1 / 2, np.sqrt(3) / 2]])
        #basis displacements for nodes
        #self.rm = np.array([r0[i] - np.sqrt(3)*self.X[i]*self.p[i] + self.X[(i-1)%3]*self.a[(i+1)%3] for i in range(3)])

        # print(self.rm - r0)
        for j in range(self.Ny):
            for i in range(self.Nx):
                # Node position for unit cell (i, j) with staggered rows
                self.rnodes.extend([(i - (j + 1) // 2) * self.a[0] + j * self.a[1]])

                # Index of node in this unit cell
                uc_ind_0 = i + j * self.Nx
                # External bonds to neighbors (row parity matters)
                if j%2 == 0:
                    uc_ind_0_up_left = uc_ind_0 + self.Nx
                    uc_ind_0_up_right = uc_ind_0 + self.Nx + 1
                else:
                    uc_ind_0_up_left = uc_ind_0 + self.Nx - 1
                    uc_ind_0_up_right = uc_ind_0 + self.Nx
                self.ind_bonds.extend(
                    [
                        # Horizontal bond
                        (uc_ind_0, (uc_ind_0 + 1) % self.Nnodes),
                        # Two upward bonds
                        (uc_ind_0, uc_ind_0_up_left % self.Nnodes),
                        (uc_ind_0, uc_ind_0_up_right % self.Nnodes),
                    ]
                )


        self.rnodes = np.array(self.rnodes)
        self.ind_bonds = np.array(self.ind_bonds)

        # Delete bonds that wrap around or exceed a unit cell length
        bonds_to_del = []
        for i in range(self.ind_bonds.shape[0]):
            ra = self.rnodes[self.ind_bonds[i][0]]
            rb = self.rnodes[self.ind_bonds[i][1]]
            if dist_nodes(ra, rb) > UC_SIZE * 1.01:  # allow for small numerical errors
                bonds_to_del.append(i)
        self.ind_bonds = np.delete(self.ind_bonds, bonds_to_del, 0)


    def add_randomness_nodes(self, s: float = 1 / 24, m: float = 0) -> None:
        """
        Add Gaussian noise to interior node positions (boundary nodes fixed).

        Parameters
        ----------
        s
            Noise amplitude (non-dimensional, scaled by UC_SIZE).
        m
            Reserved for mean shift (currently unused).
        """
        self.noise_amp = s  # non-dimensional
        mean_noise = 0
        stdev = s * UC_SIZE  # mm
        # Identify top/bottom boundary nodes by y-coordinate
        ind_max = np.where(self.rnodes == np.max(self.rnodes[:, 1]))[0]
        ind_min = np.where(self.rnodes == np.min(self.rnodes[:, 1]))[0]

        rand_gitter = np.random.normal(mean_noise, stdev, (self.rnodes.shape[0], 2))
        rand_gitter[ind_max] = 0
        rand_gitter[ind_min] = 0

        self.rnodes += rand_gitter


    def make_Bond_End_Array(self) -> None:
        """
        Create bond endpoint pairs for plotting.
        """
        self.bond_ends = []
        for i in range(self.ind_bonds.shape[0]):
            self.bond_ends.append([(self.rnodes[self.ind_bonds[i][0]][0], self.rnodes[self.ind_bonds[i][0]][1]), (self.rnodes[self.ind_bonds[i][1]][0], self.rnodes[self.ind_bonds[i][1]][1])])


    def strengthen_corners(self, factor: float = 10, box_size: int = 4) -> None:
        """
        Strengthen the corner bonds by a given factor. This is to avoid boundary failures/delamination due to size effects.
        It is important to note that this is done for a region and not just the single corner bonds.
        
        factor: by how much to strengthen the bonds.
        self.toughened_bonds: array that stores the indeces of the bonds that should be toughened.
        box_size: size of the square region (in unit cells) at each corner to be toughened.
        """

        self.toughened_bonds = []
        # Find lattice bounds (lattice is not centered at origin)
        min_x = np.min(self.rnodes[:,0])
        min_y = np.min(self.rnodes[:,1])
        max_x = np.max(self.rnodes[:,0])
        max_y = np.max(self.rnodes[:,1])

        for i in range(self.ind_bonds.shape[0]):
            r_com = bond_CoM(self.rnodes[self.ind_bonds[i][0]], self.rnodes[self.ind_bonds[i][1]])

            """
            # Exclude corners by checking if CoM is in corner boxes.
            if (r_com[0] < min_x + box_size*uc_size and r_com[1] < min_y + box_size*uc_size) or \
               (r_com[0] > max_x - box_size*uc_size and r_com[1] < min_y + box_size*uc_size) or \
               (r_com[0] < min_x + box_size*uc_size and r_com[1] > max_y - box_size*uc_size) or \
               (r_com[0] > max_x - box_size*uc_size and r_com[1] > max_y - box_size*uc_size):
                self.toughened_bonds.append((i, factor))
            """

                # Exclude edges (x > max_x - box_size*uc_size etc.) if wanted.
            if r_com[0] > max_x - box_size * UC_SIZE:
                self.toughened_bonds.append((i, factor))

        self.toughened_bonds = np.array(self.toughened_bonds)


    def make_bnfiles(
        self,
        xpins_tb: int = 1,
        xroller_lr: int = 0,
        n: float | str = 10.0,
        toughen_corners: bool = False,
    ) -> None:
        """
        Write nodes.inp and bonds.inp for the simulator.

        Include random noise in failure conditions if wanted. This is controlled through
        Weibull sampling for famax/fbmax.

        pins: boolean for applied displacement on top boundary nodes. 1 - applied displacement, 0 - fixed.
        n: weibull modulus for failure stresses. If n = "inf", then all bonds have the same failure stress.
        xpins_tb: boolean for top and bottom boundary nodes in x direction. 0 - roller, 1 - fixed. (as y always fixed or applied displacement)
        xpins_lr: boolean for left and right boundary nodes in x direction. 0 - free, 1 - rollers. (as y always free)
        """
        node_str_out = []
        bonds_str_out = []

        # The two integers at the end signify the type of constraint for each DOF:
        # 0 - roller, 1 - fixed, 2 - applied displacement.
        for i in range(self.rnodes.shape[0]):
            node_str_out.append([i, self.rnodes[i][0], self.rnodes[i][1], 0, 0])
        node_str_out = np.array(node_str_out)

        # Make sure that the nodes are not at a zero coordinate.
        #node_str_out[:,1] += 0.1
        node_str_out[:,2] += 0.1

        # Add top/bottom edge nodes to boundary nodes. Applied displacement to top boundary.
        rows, cols = np.where(node_str_out == np.amin(node_str_out[:, 2]))
        cols += 1
        node_str_out[rows, cols] = xpins_tb  # x_dof boolean
        cols += 1
        node_str_out[rows, cols] = 1  # y_dof boolean
        rows, cols = np.where(node_str_out == np.amax(node_str_out[:, 2]))
        cols += 1
        node_str_out[rows, cols] = xpins_tb  # x_dof boolean
        cols += 1
        node_str_out[rows, cols] = 2  # y_dof boolean - applied displacement

        # Apply roller boundary conditions to left and right edges.
        rows, cols = np.where(node_str_out == np.amin(node_str_out[:, 1]))
        cols += 2
        node_str_out[rows, cols] = xroller_lr  # x_dof boolean
        rows, cols = np.where(node_str_out == np.amax(node_str_out[:, 1]))
        cols += 2
        node_str_out[rows, cols] = xroller_lr  # x_dof boolean


        # node_str_out[:,1] -= np.amin(node_str_out[:,1])
        node_str_out[:,2] -= np.amin(node_str_out[:,2])

        node_file = open(self.dir + "/nodes.inp", "w")
        bond_file = open(self.dir + "/bonds.inp", "w")
        node_file.write(f"{node_str_out.shape[0]} NumNodes"+"\n")
        bond_file.write(f"{self.ind_bonds.shape[0]} NumBonds"+"\n")

        for i in range(node_str_out.shape[0]):
            node_file.write(
                f"{round(node_str_out[i][1], DECIMAL_TOL)} {round(node_str_out[i][2], DECIMAL_TOL)} "
                f"{int(node_str_out[i][3])} {int(node_str_out[i][4])} Positionxy_boundarybool\n"
            )

        if n != "inf":
            # Divide by 3 to get same noise for the three discretized elements of a bond.
            fa_weibull = np.random.weibull(n, size=self.ind_bonds.shape[0] // 3)
        else:
            fa_weibull = np.ones(self.ind_bonds.shape[0] // 3)
        
        #fb_noises = np.random.normal(loc=0.0, scale=sigma_b, size=self.ind_bonds.shape[0]//3)
        fa_noises = np.hstack((fa_weibull, np.repeat(fa_weibull, 2)))
        #fb_noises = np.hstack((fb_noises, np.repeat(fb_noises, 2)))

        for i in range(self.ind_bonds.shape[0]):
            dist = dist_nodes(self.rnodes[self.ind_bonds[i][0]], self.rnodes[self.ind_bonds[i][1]])
            # Make toughened corners if wanted, by setting a factor to 1 or a given value.
            factor = 1
            if toughen_corners:
                for j in range(self.toughened_bonds.shape[0]):
                    if i == self.toughened_bonds[j][0]:
                        factor = self.toughened_bonds[j][1]
                        fa_noises[i] = 1.0  # reset noise to 1 for toughened bonds
                        break

            # Write bond information to file.
            bond_file.write(
                f"{self.ind_bonds[i][0]} {self.ind_bonds[i][1]} {round(k / dist, DECIMAL_TOL)} "
                f"{round(dist, DECIMAL_TOL)} {famax * fa_noises[i] * factor} {fbmax * fa_noises[i] * factor} "
                "Node1_Node2_K_R0_FMAX_FBMAX\n"
            )
        
        node_file.close()
        bond_file.close()
        
        """
        #Plot histogram of failure stresses
        plt.hist(fa_noises, bins = 30)
        plt.xlabel("Failure Stress")
        plt.ylabel("Frequency")
        plt.title("Histogram of Failure Stresses")
        plt.show()
        """


    def cleanup_arrays(self) -> None:
        """
        After "deleting" nodes and bonds, you have to actually delete them from the arrays
        and re-index the nodes and bonds.
        """
        cnt = 0
        tot_cnt = 0
        for i in range(self.rnodes.shape[0]):
            if np.all(self.rnodes[i]) == 0:
                cnt += 1
                tot_cnt += 1
            else:
                self.ind_bonds[self.ind_bonds >= i - tot_cnt] -= cnt
                cnt = 0

        self.rnodes = np.delete(self.rnodes, np.any(self.rnodes == 0.0, axis=1), axis=0)[self.Nx :]
        self.ind_bonds = np.delete(self.ind_bonds, np.any(self.ind_bonds == -2, axis=1), axis=0) - self.Nx


    def add_boundary_bonds(self) -> None:
        """
        Add boundary bonds to be able to simulate de-lamination if need be.
        """
        topnodes_ind = np.where(self.rnodes[:,1] == np.max(self.rnodes[:,1]))[0]
        botnodes_ind = np.where(self.rnodes[:,1] == np.min(self.rnodes[:,1]))[0]

        for i in range(topnodes_ind.shape[0] - 1):
            self.ind_bonds = np.vstack((self.ind_bonds, np.array([topnodes_ind[i], topnodes_ind[i+1]])))
        for i in range(botnodes_ind.shape[0] - 1):
            self.ind_bonds = np.vstack((self.ind_bonds, np.array([botnodes_ind[i], botnodes_ind[i+1]])))


    def discretize_lattice(self) -> None:
        """
        Break all bonds in three parts and replace them with extra bonds. Add nodes in between appropriately.
        (Starting with *---* bonds and end up with *-*-*-*)
        """
        bond_array_size = self.ind_bonds.shape[0]
        discr_el = 3

        # Go over all bonds and add intermediate nodes with appropriate bond indexing.
        for i in range(bond_array_size):
            rab = self.rnodes[self.ind_bonds[i][1]] - self.rnodes[self.ind_bonds[i][0]]
            rk = []
            rk.append(self.rnodes[self.ind_bonds[i][0]] + 1 / discr_el * rab)
            for j in range(1, discr_el - 1):
                rk.append(rk[-1] + 1 / discr_el * rab)

            # Add the new nodes to array
            self.rnodes = np.append(self.rnodes, np.array(rk), axis = 0)

            # Add the new bonds to array
            Nnodes = self.rnodes.shape[0]
            self.ind_bonds = np.append(
                self.ind_bonds,
                np.array([[self.ind_bonds[i][0], Nnodes - discr_el + 1]]),
                axis=0,
            )
            for j in range(1, discr_el - 1):
                self.ind_bonds = np.append(
                    self.ind_bonds,
                    np.array([[Nnodes - discr_el + j, Nnodes - discr_el + 1 + j]]),
                    axis=0,
                )

            # Change existing bond to point to closest node in line.
            self.ind_bonds[i][0] = Nnodes - 1


    def make_bends(self) -> None:
        """
        Go over arrays ind_bonds and rnodes, and find all adjacent triplets, and register them as bends.
        """

        self.bends = []  # store indices for bend + equilibrium angle
        for i in range(self.rnodes.shape[0]):
            rows, cols = np.where(self.ind_bonds == i)

            newbends = []
            for j in range(rows.shape[0]):
                ind1 = self.ind_bonds[rows[j]][(cols[j] + 1)%2]
                newbendstemp = []
                for k in range(cols.shape[0]):
                    if k != j:
                        ind3 = self.ind_bonds[rows[k]][(cols[k] + 1)%2]

                        r12 = self.rnodes[i] - self.rnodes[ind1]
                        r23 = self.rnodes[ind3] - self.rnodes[i]
                        r12n = r12 / np.linalg.norm(r12)
                        r23n = r23 / np.linalg.norm(r23)
                        r23northo = np.array([r23n[1], -1 * r23n[0]])

                        c = np.dot(r12n, r23n)
                        s = np.dot(r12n, r23northo)
                        angle = np.arctan2(s, c)

                        # Last argument tells if bend is double or single
                        newbendstemp.append([ind1, i, ind3, angle * 180 / np.pi, 1])

                # Keep only largest and smallest angle for each central node
                newbendstemp = np.array(newbendstemp)
                if newbendstemp.shape[0] == 1:
                    newbends.extend(newbendstemp)
                else:
                    indeces = [np.argmin(newbendstemp[:, 3]), np.argmax(newbendstemp[:, 3])]
                    newbends.extend(newbendstemp[indeces])

            # Remove duplicates unless only two neighbors exist
            newbends = np.array(newbends)
            if newbends.shape[0] == 2:
                newbends[0][4] = 2
                self.bends.append(newbends[0])
            else:
                # Ensure bends are written with edge nodes in ascending order
                for j in range(newbends.shape[0]):
                    if newbends[j][0] > newbends[j][2]:
                        temp = newbends[j][0]
                        newbends[j][0] = newbends[j][2]
                        newbends[j][2] = temp
                        newbends[j][3] *= -1

                newbends = np.unique(newbends, axis = 0)
                self.bends.extend(newbends)


        self.bends = np.array(self.bends)

        # Write bends file
        bend_file = open(self.dir + "/bends.inp", "w")
        bend_file.write(f"{self.bends.shape[0]} NumBends"+"\n")
        for i in range(self.bends.shape[0]):
            dist1 = dist_nodes(self.rnodes[int(self.bends[i][0])], self.rnodes[int(self.bends[i][1])])
            dist2 = dist_nodes(self.rnodes[int(self.bends[i][1])], self.rnodes[int(self.bends[i][2])])
            bend_file.write(
                f"{int(self.bends[i][0])} {int(self.bends[i][1])} {int(self.bends[i][2])} "
                f"{round(kbe / min(dist1, dist2), DECIMAL_TOL)} {round(self.bends[i][3], DECIMAL_TOL)} "
                "Node1_Node2_Node3_K_A0\n"
            )
        bend_file.close()


    def make_simfile(self, dx: float = 0.1) -> None:
        """
        Store simulation parameters in a sim.inp file for the simulator.
        """
        sim_file = open(self.dir + "/sim.inp", "w")
        sim_file.write(f"{dx} strainStep" + "\n")
        sim_file.write(f"{kbe:.5} bendStiff" + "\n")
        sim_file.write(f"{yy:.5} beamThickness" + "\n")
        sim_file.write(f"{ii:.5} secondMomentOfArea" + "\n")
        if self.noise_amp:
            sim_file.write(f"{round(1 / self.noise_amp, 1)} inverseNoise" + "\n")
        else:
            sim_file.write(f"0.0 inverseNoise" + "\n")
        sim_file.close()


    def make_crack(self, length: int, side: str = "r") -> None:
        """
        Introduce a crack in the lattice. This should be done before discretizing and adding bends.

        length: in unit cells.
        side: "l" or "r" or "m" - left, right or middle (notch) of the lattice.

        Crack is at the midpoint of y by default.
        """
        # Identify bonds that cross the mid-height and remove those within the crack region.
        mid_height = 0.5 * (np.min(self.rnodes[:, 1]) + np.max(self.rnodes[:, 1]))

        bonds_to_del = []
        for i in range(self.ind_bonds.shape[0]):
            if (self.rnodes[self.ind_bonds[i][0]][1] - mid_height) * (
                self.rnodes[self.ind_bonds[i][1]][1] - mid_height
            ) < 0:  # checks whether dy1*dy2<0, where dyi is taken wrt. mid_height
                if side == "r":
                    if self.rnodes[self.ind_bonds[i][0]][0] > UC_SIZE * (self.Nx - length + 1):
                        bonds_to_del.append(i)
                elif side == "l":
                    if self.rnodes[self.ind_bonds[i][0]][0] < UC_SIZE * (length - 1):
                        bonds_to_del.append(i)
                elif side == "m":
                    if (
                        (self.rnodes[self.ind_bonds[i][0]][0] + self.rnodes[self.ind_bonds[i][1]][0]) / 2
                        > UC_SIZE * (self.Nx / 2 - length / 2)
                        and (self.rnodes[self.ind_bonds[i][0]][0] + self.rnodes[self.ind_bonds[i][1]][0]) / 2
                        < UC_SIZE * (self.Nx / 2 + length / 2)
                    ):
                        bonds_to_del.append(i)
                else:
                    raise ValueError("Incorrect argument for variable side. Check documentation.")

        #Delete the bonds that are to be deleted.
        self.ind_bonds = np.delete(self.ind_bonds, bonds_to_del, 0)




def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate input files for a triangular lattice.")
    parser.add_argument("SR", type=float, help="Slenderness ratio (unitless).")
    parser.add_argument("nmod", type=float, help="Weibull modulus (use a large value for ~deterministic).")
    return parser.parse_args()


def _compute_material_params(sr: float) -> Tuple[float, float, float, float, float, float]:
    """Compute material parameters from the slenderness ratio."""
    yy_local = UC_SIZE / sr  # Thickness of beams (mm)
    depth_local = AA / yy_local  # Depth of cross section (mm)
    ii_local = (depth_local * yy_local**3) / 12  # Second moment of area (mm^4)
    e_local = 3120  # Young's modulus (MPa)
    k_local = e_local * AA  # Bond stiffness
    kbe_local = e_local * ii_local  # Bend stiffness
    return yy_local, ii_local, k_local, kbe_local, FAMAX, FBMAX


def main() -> None:
    args = _parse_args()

    # CALCULATED/DEFINED MATERIAL QUANTITIES
    global yy, ii, k, kbe, famax, fbmax
    yy, ii, k, kbe, famax, fbmax = _compute_material_params(args.SR)

    # Make lattice and files.
    lattice = Triangular_Lattice(OUTPUT_DIR, NX, NY)

    print("Preparing geometry...")
    lattice.makeNodes_AddBondInd()

    lattice.make_crack(NX // 2, side="l")  # make a notch in the middle of the lattice
    lattice.add_randomness_nodes(NOISE / UC_SIZE)
    lattice.discretize_lattice()
    # lattice.strengthen_corners(factor=100, box_size=4)

    print("Making files...")
    lattice.make_bnfiles(xpins_tb=1, xroller_lr=1, n=args.nmod, toughen_corners=False)
    lattice.make_bends()
    lattice.make_simfile(0.125)


# ---- Configuration ----
NX = 60
NY = 30
UC_SIZE = 4.0  # unit cell size (mm)
AA = 3.8  # area of cross section (mm^2)
OUTPUT_DIR = "../src/IO/triang"
FAMAX = 45  # maximum axial failure stress (MPa)
FBMAX = 85  # maximum bending failure stress (MPa)
NOISE = 0.0  # in mm


if __name__ == "__main__":
    main()


"""
#Plotting
print("Plotting...")
T.make_Bond_End_Array()

fig, ax = plt.subplots(figsize=(5, 5))
# set axes limits manually because Collections do not take part in autoscaling
ax.set_xlim(-2*uc_size, (Nx+1)*uc_size)
ax.set_ylim(-2*uc_size, (Ny+1)*uc_size)
ax.set_aspect("equal")
# create a LineCollection with the half-circles
# its properties can be set per line by passing a sequence (here used for *colors*)
# or they can be set for all lines by passing a scalar (here used for *linewidths*)
line_collection = mc.LineCollection(T.bond_ends, linewidths=2)
ax.add_collection(line_collection)

#Change colour of lines to grey and width of lines to 0.5
line_collection.set_linewidth(0.5)
line_collection.set_color('grey')


#ax.scatter(T.rnodes[:,0], T.rnodes[:,1], c = "r")
#for i in range(T.rnodes.shape[0]//3):
#    ax.annotate(f"{i}", (T.rnodes[:,0][i], T.rnodes[:,1][i]))


plt.show()
"""