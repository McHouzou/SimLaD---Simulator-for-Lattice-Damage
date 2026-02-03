# SimLaD — Simulator for Lattice Damage

SimLaD is a C++ simulator for fracture and damage in 2D beam‑lattices, accompanied by Python scripts for input generation and post‑processing. It supports both underconstrained, critical (Maxwell) and overconstrained lattices, includes bending and axial contributions, and uses a FIRE‑based relaxation scheme to resolve quasi‑static fracture events within a user defined tolerance.

## Highlights
- Beam‑lattice mechanics 
with axial and bending contributions.
- Quasi‑static fracture via strain stepping and event‑driven topology updates.
- Input generators for triangular and kagome lattices.
- Post‑processing utilities for visualization and analysis.

## Repository structure
- src/lattice_fire.cpp — Main C++ simulator (FIRE relaxation + fracture/topology updates).
- analysis/ — Python tools for input generation and analysis.
	- InputCreator_triang.py — Triangular lattice input generator.
	- InputCreator_kagome.py — Kagome lattice input generator (CLI included).
	- read_output.py — Output readers for positions, bonds, and boundaries.
	- Vis_lattice.py — Visualization script for lattice state.

## Build
Build the simulator with a C++ compiler (e.g. clang++, g++), using optimization flags (e.g. -O3 or -O2). Example build command:
clang++ -O3 -march=native -o lf_exec lattice_fire.cpp

## Run
The simulator expects input files in a directory (nodes.inp, bonds.inp, bends.inp, sim.inp). Provide input and output directories using flags:

- -i path/to/input/
- -o path/to/output/

Example:
./lf_exec -i IO/triang/ -o IO/triang/

## Input generators
Generate triangular lattice inputs:

python3 analysis/InputCreator_triang.py <SR> <nmod>

Generate kagome lattice inputs:

python3 analysis/InputCreator_kagome.py <x_amp> <noise> --nx 10 --ny 10 --out-dir ../src/IO/kagome --plot

## Output files
The simulator writes:
- pos.txt — node positions per frame.
- bondnew.txt — bond stresses per frame.
- energy.txt — axial and bending energy per frame.
- bforce.txt — boundary forces and current lattice height.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Citation
If you use this code or the datasets in your work, please cite the associated publication.
