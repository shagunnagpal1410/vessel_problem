# 🚀 Blood Vessel Simulation

This C++ project simulates time-dependent blood velocity and pressure fields inside a 3D vessel and its branches using meshfree numerical methods (LSMPS) and Chorin's projection method to solve incompressible Navier-Stokes equations. It reads `vessel_points.csv` and `branch_points.csv`, computes velocity fields, and updates pressure and velocity over time. Sparse linear systems for pressure correction are solved using Eigen's SparseLU.

## 📂 Project Structure

.
├── 7july.cpp              # Main simulation code (velocity & pressure solver)
├── eigen/                 # Eigen library folder (already included)
├── vessel_points.csv      # CSV file for vessel geometry & point types
├── branch_points.csv      # CSV file for branch geometry & point types
├── .vscode/
│   └── tasks.json         # VS Code build task with Eigen included
├── Velocity1.csv          # Example output file for velocity at timestep 1
├── non_zero_entries.csv   # CSV reporting non-zero matrix entries

## ⚙️ How to Build & Run

Open this folder in Visual Studio Code. Press `Ctrl + Shift + B` to compile (runs the build task configured in `tasks.json` which already includes Eigen). Run the executable. Enter:
- Maximum simulation time (e.g. `0.001`)
- Time step (e.g. `0.001`)
The simulation outputs velocity CSV files for each time step.

## 📥 Input Files

### vessel_points.csv and branch_points.csv

These CSV files store 3D point coordinates and point types:

x, y, z, type

- `x, y, z` → 3D coordinates (in meters or desired units)
- `type` → point classification:
    - `0` → Fluid point
    - `1` → Wall (stationary)
    - `2` → Inlet (positive z-direction velocity)
    - `3` → Inlet (positive x-direction velocity)
    - `4` → Outlet (non-zero z-direction velocity)

## 🔬 What Does This Code Do?

✅ Reads vessel and branch geometry from CSV files.  
✅ Partitions points into voxels for fast neighbor searching.  
✅ Calculates spatial derivatives using Least-Squares Moving Particle Semi-implicit (LSMPS).  
✅ Applies Chorin’s projection method to:
- Solve for intermediate velocity field v*
- Compute pressure correction
- Project velocity to satisfy incompressibility
✅ Builds and solves a sparse linear system for pressure using Eigen::SparseLU.  
✅ Outputs:
- Velocity CSVs at each time step
- A diagnostic CSV showing non-zero entries in the sparse pressure matrix

## 📝 Output Files

- `Velocity1.csv` → Velocities at time step 1
- `Velocity2.csv` → Velocities at time step 2
- `non_zero_entries.csv` → Number of non-zero entries per row in the sparse matrix for debugging

Each velocity CSV contains:
x, y, z, vx, vy, vz

## 🛠 Dependencies

- Eigen (already included in the eigen/ folder)
- A C++ Compiler (e.g. g++, MSVC, clang)


