"""Generate kagome lattice input files.

This module builds a kagome lattice, applies optional clean-up/crack/disorder
operations, and writes node, bond, bend, and simulation input files.

It can also be executed directly to generate inputs (see `main()`).
"""

from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections as mc

# ---- Formatting / numeric precision ----
decimal_tol = 5

# ---- Geometric constants (mm) ----
uc_size = 12.0  # unit cell size
t = 0.8  # beam thickness
b = 6.0  # in-plane depth of beams

# ---- Material constants ----
E = 2950
k = 14160
kbe = 755.2
famax = 35
fbmax = 108


def dist_nodes(ra, rb):
    """
    Return the Euclidean distance between two nodes,
    given their position vectors ra and rb.
    """
    rab = rb - ra
    return np.sqrt(np.dot(rab, rab))


class Kagome_lattice(object):


    def __init__(self, dir, Nx, Ny, X):
        """
        - dir: Output directory for input files.
        - Nx, Ny: Number of unit cells in the x and y direction.
        - X: Distortion parameters {x1, x2, x3, z} (kagome chirality).
        """
        self.dir = dir
        self.Nx = Nx + 1
        self.Ny = Ny + 2
        self.Nnodes = self.Nx*self.Ny*3
        self.X = X


    def makeNodes_AddBondInd(self):
        """
        Construct two arrays: node positions and bond indices.
        """
        self.rnodes = []
        self.ind_bonds = []

        # Basis vectors for the kagome lattice
        self.a = uc_size * np.array([[np.cos(2 * i * np.pi / 3), np.sin(2 * i * np.pi / 3)] for i in range(3)])
        self.p = uc_size * np.array([[-1 * np.sin(2 * i * np.pi / 3), np.cos(2 * i * np.pi / 3)] for i in range(3)])
        r0 = np.array([self.a[0]/2, self.a[0]/2 + self.a[1]/2, [0,0]])

        self.rm = np.array([r0[i] - np.sqrt(3)*self.X[i]*self.p[i] + self.X[(i-1)%3]*self.a[(i+1)%3] for i in range(3)])

        # print(self.rm - r0)
        for j in range(self.Ny):
            for i in range(self.Nx):
                self.rnodes.extend(self.rm + (i - j//2) * self.a[0] - j * self.a[2])

                # Internal bonds in the unit cell
                uc_ind_0 = 3 * (i + self.Nx * j)  # index of 0 node
                self.ind_bonds.extend([(uc_ind_0, 2 + uc_ind_0), (uc_ind_0, 1 + uc_ind_0), (1 + uc_ind_0, 2 + uc_ind_0)])

                # External bonds to neighboring unit cells
                if j%2 != 0:
                    uc_ind_0_up = uc_ind_0 + self.Nx*3
                    uc_ind_0_down = uc_ind_0 - self.Nx*3
                else:
                    uc_ind_0_up = uc_ind_0 + self.Nx*3 - 3
                    uc_ind_0_down = uc_ind_0 - self.Nx*3 - 3
                self.ind_bonds.extend([(uc_ind_0, 2 + uc_ind_0 + 3), (uc_ind_0 + 1, uc_ind_0_up), (2 + uc_ind_0, 1 + uc_ind_0_down)])


        self.rnodes = np.array(self.rnodes)
        self.ind_bonds = np.array(self.ind_bonds)

        # Delete bonds that wrap around or are otherwise unphysical
        bonds_to_del = []
        for i in range(self.ind_bonds.shape[0]):
            ra = self.rnodes[self.ind_bonds[i][0] % self.Nnodes]
            rb = self.rnodes[self.ind_bonds[i][1] % self.Nnodes]
            if dist_nodes(ra, rb) > uc_size:
                bonds_to_del.append(i)
        self.ind_bonds = np.delete(self.ind_bonds, bonds_to_del, 0)


    def make_Bond_End_Array(self):
        """
        Create an array with bond ends in x-y space for visualization.
        """
        self.bond_ends = []
        for i in range(self.ind_bonds.shape[0]):
            self.bond_ends.append([(self.rnodes[self.ind_bonds[i][0]][0], self.rnodes[self.ind_bonds[i][0]][1]), (self.rnodes[self.ind_bonds[i][1]][0], self.rnodes[self.ind_bonds[i][1]][1])])


    def cleanup_lattice(self):
        """
        Delete bonds/nodes outside a rectangular domain.
        # domain = [(x0,y0), (x1,y1)] - lower left and upper right corner of rectangle
        """
        self.rnodes[:,1] += 1.1
        # Cleanup lower boundary
        # ymin = (self.a[0] + self.a[1] + self.rm[1])[1] - 0.1
        tol = 0.01
        ymin = self.rnodes[self.Nx*3+1][1]
        xmin = min(self.rnodes[self.Nx*3+2][0], self.rnodes[self.Nx*3+1][0], self.rnodes[2*self.Nx*3][0])
        xmax = max(self.rnodes[self.Nx*3-3][0], self.rnodes[self.Nx*3-3][0], self.rnodes[2*self.Nx*3-1][0])
        #xmin = self.rnodes[self.Nx*3+2][0]
        ind_node_to_del = []
        ind_bond_to_del = []

        for i in range(self.ind_bonds.shape[0]):
            b = self.ind_bonds[i]
            if self.rnodes[b[0]][1] < ymin-tol:
                ind_node_to_del.append(b[0])
                ind_bond_to_del.append(i)
            elif self.rnodes[b[1]][1] < ymin-tol:
                ind_node_to_del.append(b[1])
                ind_bond_to_del.append(i)

        ind_node_to_del = np.unique(np.array(ind_node_to_del))
        ind_bond_to_del = np.array(ind_bond_to_del)
        self.ind_bonds[ind_bond_to_del] = 0
        self.rnodes[ind_node_to_del] = 0


        # Cleanup leftmost boundary
        ind_node_to_del = []
        ind_bond_to_del = []

        for i in range(self.ind_bonds.shape[0]):
            b = self.ind_bonds[i]
            if (self.rnodes[b[0]][0] < xmin):
                ind_node_to_del.append(b[0])
                ind_bond_to_del.append(i)
            elif self.rnodes[b[1]][0] < xmin:
                ind_node_to_del.append(b[1])
                ind_bond_to_del.append(i)

        ind_node_to_del = np.unique(np.array(ind_node_to_del))
        ind_bond_to_del = np.array(ind_bond_to_del)
        self.ind_bonds[ind_bond_to_del] = 0
        self.rnodes[ind_node_to_del] = 0


        # Cleanup rightmost boundary
        # xmax = (self.a[0]*(self.Nx-1))[0]
        ind_node_to_del = []
        ind_bond_to_del = []

        for i in range(self.ind_bonds.shape[0]):
            b = self.ind_bonds[i]
            if self.rnodes[b[0]][0] > xmax+tol:
                ind_node_to_del.append(b[0])
                ind_bond_to_del.append(i)
            elif self.rnodes[b[1]][0] > xmax+tol:
                ind_node_to_del.append(b[1])
                ind_bond_to_del.append(i)

        ind_node_to_del = np.unique(np.array(ind_node_to_del))
        ind_bond_to_del = np.array(ind_bond_to_del)
        self.ind_bonds[ind_bond_to_del] = 0
        self.rnodes[ind_node_to_del] = 0


        # Cleanup some default bonds at the ends of the trimmed region
        ind_node_to_del = []
        ind_bond_to_del = []
        for j in range(self.Ny//2):
            indbl = 3*self.Nx*j*2 + 1
            indbr = 3*(2*self.Nx*(j+1)-1) + 1
            row,col = np.where(self.ind_bonds == indbr)
            ind_bond_to_del.extend(row)
            ind_node_to_del.append(indbr)

            row,col = np.where(self.ind_bonds == indbl)
            ind_bond_to_del.extend(row)
            ind_node_to_del.append(indbl)

        ind_node_to_del = np.unique(np.array(ind_node_to_del))
        ind_bond_to_del = np.array(ind_bond_to_del)
        if ind_bond_to_del.any():
            self.ind_bonds[ind_bond_to_del] = 0
        if ind_node_to_del.any():
            self.rnodes[ind_node_to_del] = 0


    def make_bnfiles(self, xpins = 1):
        """
        Write input files for bonds and nodes.
        xpins: x-constraint for top/bottom nodes (0 free/roller, 1 fixed).
        """
        node_str_out = []
        bonds_str_out = []

        # The two integers at the end signify the constraint type for each DOF:
        # 0 - free, 1 - fixed, 2 - applied displacement
        for i in range(self.rnodes.shape[0]):
            node_str_out.append([i, self.rnodes[i][0], self.rnodes[i][1], 0, 0])
        node_str_out = np.array(node_str_out)

        # Add edge nodes to boundary nodes. Applied displacement on top boundary.
        rows, cols = np.where(node_str_out == np.amin(node_str_out[:,2])); 
        cols += 1; node_str_out[rows, cols] = xpins; #x_dof boolean
        cols += 1; node_str_out[rows, cols] = 1; #y_dof boolean
        rows, cols = np.where(node_str_out == np.amax(node_str_out[:,2])); 
        cols += 1; node_str_out[rows, cols] = xpins; #x_dof boolean
        cols += 1; node_str_out[rows, cols] = 2; #y_dof boolean - applied displacement


        # node_str_out[:,1] -= np.amin(node_str_out[:,1])
        # node_str_out[:,2] -= np.amin(node_str_out[:,2])

        node_file = open(self.dir + "/nodes.inp", "w")
        bond_file = open(self.dir + "/bonds.inp", "w")
        node_file.write(f"{node_str_out.shape[0]} NumNodes"+"\n")
        bond_file.write(f"{self.ind_bonds.shape[0]} NumBonds"+"\n")

        for i in range(node_str_out.shape[0]):
            node_file.write(f"{round(node_str_out[i][1], decimal_tol)} {round(node_str_out[i][2], decimal_tol)} {int(node_str_out[i][3])} {int(node_str_out[i][4])} Positionxy_boundaryboolx_y"+'\n')

        for i in range(self.ind_bonds.shape[0]):
            dist = dist_nodes(self.rnodes[self.ind_bonds[i][0]], self.rnodes[self.ind_bonds[i][1]])
            bond_file.write(f"{self.ind_bonds[i][0]} {self.ind_bonds[i][1]} {round(k/dist, decimal_tol)} {round(dist, decimal_tol)} {famax} {fbmax} Node1_Node2_K_R0_FMAX_FBMAX"+'\n')
        node_file.close()
        bond_file.close()


    def cleanup_arrays(self):
        """
        After "deleting" nodes and bonds, remove them and re-index arrays.
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

        self.rnodes = np.delete(self.rnodes, np.any(self.rnodes == 0.0, axis = 1), axis = 0)[self.Nx:]
        self.ind_bonds = np.delete(self.ind_bonds, np.any(self.ind_bonds == -2, axis = 1), axis = 0) - self.Nx


    def add_boundary_bonds(self):
        """
        Add boundary bonds to enable delamination if needed.
        """
        topnodes_ind = np.where(self.rnodes[:,1] == np.max(self.rnodes[:,1]))[0]
        botnodes_ind = np.where(self.rnodes[:,1] == np.min(self.rnodes[:,1]))[0]

        for i in range(topnodes_ind.shape[0] - 1):
            self.ind_bonds = np.vstack((self.ind_bonds, np.array([topnodes_ind[i], topnodes_ind[i+1]])))
        for i in range(botnodes_ind.shape[0] - 1):
            self.ind_bonds = np.vstack((self.ind_bonds, np.array([botnodes_ind[i], botnodes_ind[i+1]])))


    def discretize_lattice(self):
        """
        Split each bond into three segments by adding intermediate nodes.
        (Starting with *---* bonds and ending with *-*-*-*)
        """
        bond_array_size = self.ind_bonds.shape[0]
        discr_el = 3

        # Go over all bonds and add intermediate nodes with proper indexing
        for i in range(bond_array_size):
            rab = self.rnodes[self.ind_bonds[i][1]] - self.rnodes[self.ind_bonds[i][0]]
            rk = []
            rk.append(self.rnodes[self.ind_bonds[i][0]] + 1/discr_el * rab)
            for j in range(1, discr_el-1):
                rk.append(rk[-1] + 1/discr_el * rab)

            # Add the new nodes to array
            self.rnodes = np.append(self.rnodes, np.array(rk), axis = 0)

            # Add the new bonds to array
            Nnodes = self.rnodes.shape[0]
            self.ind_bonds = np.append(self.ind_bonds, np.array([[self.ind_bonds[i][0], Nnodes - discr_el + 1]]), axis = 0)
            for j in range (1, discr_el - 1):
                self.ind_bonds = np.append(self.ind_bonds, np.array([[Nnodes - discr_el + j, Nnodes - discr_el + 1 + j]]), axis = 0)

            # Update existing bond to point to closest node in line
            self.ind_bonds[i][0] = Nnodes - 1


    def make_bends(self):
        """
        Find adjacent triplets (i-j-k) and register them as bending elements.
        """

        self.bends = []  # store indices for bend + equilibrium angle
        for i in range(self.rnodes.shape[0]):
            rows,cols = np.where(self.ind_bonds == i)

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
                        r23northo = np.array([r23n[1], -1*r23n[0]])

                        c = np.dot(r12n, r23n)
                        s = np.dot(r12n, r23northo)
                        angle = np.arctan2(s,c)

                        # Last argument indicates single vs. double bend
                        newbendstemp.append([ind1, i, ind3, angle*180/np.pi, 1])

                #Make sure there are no bends that shouldn't be there (keep only largest and smallest angle)
                newbendstemp = np.array(newbendstemp)
                if newbendstemp.shape[0] == 1:
                    newbends.extend(newbendstemp)
                else:
                    indeces = [np.argmin(newbendstemp[:,3]), np.argmax(newbendstemp[:,3])]
                    newbends.extend(newbendstemp[indeces])

            # Remove duplicates (unless only two neighbors exist)
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
            bend_file.write(f"{int(self.bends[i][0])} {int(self.bends[i][1])} {int(self.bends[i][2])} {round(kbe/min(dist1, dist2), decimal_tol)} {round(self.bends[i][3], decimal_tol)} Node1_Node2_Node3_K_A0"+'\n')
        bend_file.close()


    def make_crack(self, length, side = "r"):
        """
        Introduce a crack in the lattice (before discretization and bends).

        length: crack length in unit cells.
        side: "l" or "r" - left or right side of lattice.

        Crack is at the midpoint by default.
        """
        # Find bonds crossing the mid-height and split them by duplicating nodes.
        mid_height = 1/2 * (np.min(self.rnodes[:,1]) + np.max(self.rnodes[:,1]))
        #crack_length = 1/2
        tol = 0.01

        if side == "r":
            crack_ind = np.where((self.rnodes[:,1] >= mid_height-tol) & (self.rnodes[:,1] <= mid_height+tol) & (self.rnodes[:,0] > uc_size*(self.Nx - length+1)))[0]
        elif side == "l":
            crack_ind = np.where((self.rnodes[:,1] >= mid_height-tol) & (self.rnodes[:,1] <= mid_height+tol) & (self.rnodes[:,0] < uc_size*(length-1)))[0]
        else:
            raise ValueError("Incorrect argument for variable side. Check documentation.")

        for i in crack_ind:
            rows, cols = np.where(self.ind_bonds == i)
            cols = (cols + 1)%2  # take the other node - not the one on the midplane
            abovecrack = self.rnodes[self.ind_bonds[rows, cols]][:,1] > mid_height
            self.rnodes = np.vstack([self.rnodes, self.rnodes[i]])  # add another node at same location
            for j in range(len(abovecrack)):
                if abovecrack[j] == False:
                    self.ind_bonds[rows[j], (cols[j]+1)%2] = self.rnodes.shape[0] - 1


    def add_randomness(self, s = 1/24, m = 0):
        """
        Add randomness to node positions (excluding boundary nodes).
        """
        self.noise_amp = s #non-dimensional
        mean_noise = 0
        stdev = s * uc_size # mm
        ind_max = np.where(self.rnodes == np.max(self.rnodes[:,1]))[0]
        ind_min = np.where(self.rnodes == np.min(self.rnodes[:,1]))[0]

        rand_gitter = np.random.normal(mean_noise, stdev, (self.rnodes.shape[0], 2))
        rand_gitter[ind_max] = 0
        rand_gitter[ind_min] = 0

        self.rnodes += rand_gitter


    def rotate180(self):
        """
        Mirror lattice around x and y (180° rotation).
        Useful after cleanup_lattice for alternative unit-cell cuts.
        """
        self.rnodes *= -1
        self.rnodes[:,0] -= np.min(self.rnodes[:,0]) - 2
        self.rnodes[:,1] -= np.min(self.rnodes[:,1]) - 2


    def make_simfile(self, dx = 0.1):
        """
        Write simulation parameters to a file for the simulator.
        """
        sim_file = open(self.dir + "/sim.inp", "w")
        sim_file.write(f"{dx} strainStep"+"\n")
        sim_file.write(f"{kbe} bendStiff" + "\n")
        #sim_file.write(f"{E} Youngs_Modulus")
        sim_file.write(f"{t} Beam_Thickness\n")
        sim_file.write(f"{b} Beam_Depth\n")

        if self.noise_amp:
            sim_file.write(f"{round(1/self.noise_amp, 1)} inverseNoise" + "\n")
        else:
            sim_file.write(f"0.0 inverseNoise" + "\n")
        sim_file.close()


def _parse_args() -> argparse.Namespace:
    """Parse command-line arguments for standalone runs."""
    parser = argparse.ArgumentParser(description="Generate kagome lattice input files.")
    parser.add_argument("x_amp", type=float, help="Right-polarized distortion amplitude.")
    parser.add_argument("noise", type=float, help="Noise amplitude in mm.")
    parser.add_argument("--nx", type=int, default=10, help="Number of unit cells in x direction.")
    parser.add_argument("--ny", type=int, default=10, help="Number of unit cells in y direction.")
    parser.add_argument(
        "--out-dir",
        type=str,
        default="../src/IO/kagome",
        help="Directory where input files will be saved.",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot the generated lattice after discretization.",
    )
    return parser.parse_args()


def main() -> None:
    """Run kagome input generation from the command line."""
    args = _parse_args()

    # Distortion parameters (Right Polarized by default; Kane-Lubensky 2014 Nat. Phys)
    X = [args.x_amp, -args.x_amp, -args.x_amp]

    # Create lattice object
    K = Kagome_lattice(args.out_dir, args.nx, args.ny, X)

    print("Preparing geometry...")
    K.makeNodes_AddBondInd()
    K.cleanup_lattice()
    K.cleanup_arrays()
    # K.rotate180()

    print("Discretizing and registering bends...")
    # K.make_crack(20)
    K.add_randomness(args.noise / uc_size)
    # K.add_boundary_bonds()
    K.discretize_lattice()
    K.make_Bond_End_Array()

    print("Making files...")
    K.make_bends()
    K.make_bnfiles(xpins=1)  # xpins: x-constraints for boundary nodes
    K.make_simfile(0.2)

    if args.plot:
        # Plotting for quick visual verification
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_xlim(0, (args.nx + 1) * uc_size)
        ax.set_ylim(0, (args.ny + 1) * uc_size)
        ax.set_aspect("equal")
        line_collection = mc.LineCollection(K.bond_ends, linewidths=2)
        ax.add_collection(line_collection)
        plt.show()


if __name__ == "__main__":
    main()
