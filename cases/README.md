# Examples

Build OpenGLSA from the project root as described in the [README](../README.md#build), then follow the individual case instructions. Run each calculation from the directory containing its `control.nml`. The commands use four MPI ranks; adjust this for the machine and mesh partition.

| Case | Configuration | Included workflows |
| --- | --- | --- |
| [CYL](CYL/README.md) | Cylinder, Ma=0.1, Re=100, 315 x 301 | NS with SFD; global stability at zero out-of-plane wavenumber |
| [RAMP](RAMP/README.md) | Compression ramp, Ma=7.7, Re=4200, 1601 x 321 | NS with a positivity-preserving limiter; global stability at wavenumber 0.95 |
| [SPH](SPH/README.md) | Sphere, Ma=0.1, Re=213, 241 x 361 | Cylindrical-coordinate NS; global stability at azimuthal order 1 |
| [LDC](LDC/README.md) | Lid-driven cavity: 2D, Ma=0.1, Re=900, 201 x 201; 3D, Ma=0.2, Re=200, 101 x 101 x 101 | NS from a quiescent state; global stability with supplied base flows (2D wavenumber 7.35); optional 3D Arnoldi starting-vector generation |
| [Isentropic Vortex](Isentropic%20Vortex/README.md) | Periodic inviscid vortex, 128 x 128 | Nonlinear transport and optional initial-field generation |

The CYL, RAMP and SPH cases contain:

- `mesh/`: ASCII mesh and a binary `grid.dat`.
- `baseflow/`: nonlinear controls, mesh and a supplied flow restart.
- `GSA/`: global-stability controls, mesh and a supplied frozen base flow.

LDC has separate `2D/` and `3D/` directories, each containing `baseflow/` and `GSA/`. Its `baseflow/` inputs initialize a quiescent flow without a restart file. The 3D inputs use `grid3d.dat` and `flow3d.dat` and run with `openglsa3d`; `3D/init_vector.f90` optionally generates a prescribed Arnoldi starting vector.

The supplied `GSA/flow.dat` files (or `GSA/flow3d.dat` for 3D LDC) allow the stability runs to start without first rerunning `baseflow/`. Treat NS and GSA as separate workflows. Do not replace a supplied GSA base flow with an arbitrary NS snapshot; first establish that the replacement is the intended converged base flow.

Solver snapshots are Fortran binary files. Convert them to Tecplot ASCII with `op2tec` (2D) or `op2tec3d` (3D); consult the [User Manual](../docs/USER_MANUAL.md) for record layouts, complex-mode pairing and eigenvalue interpretation. Histories are appended and matching snapshot names may be overwritten, so use a working copy for reruns. Each `GSA/` directory contains `Eigenvalues.ref.dat`, the reference eigenvalues for the supplied inputs. Runs write `Eigenvalues.dat`, so the reference file is not overwritten.
