# Lid-driven cavity

Lid-driven cavity (LDC) examples with separate 2D and fully 3D workflows. Both use LINEAR5/CD6 and `BC_Lid_Driven_Cavity`. Build the project as described in the [project README](../../README.md#build).

## Files and settings

| Directory | Inputs and active controls |
| --- | --- |
| `2D/baseflow/` | `control.nml`, `grid.dat`; 201 x 201, Ma=0.1, Re=900, `RunMode=0`, `IF_init=2`, `dt=0.00015`, `End_time=400` |
| `2D/GSA/` | `control.nml`, `grid.dat`, `flow.dat`, `Eigenvalues.ref.dat`; same mesh and flow parameters, `RunMode=1`, `IF_init=0`, `dt=0.0001`, propagator horizon 0.5, `wavenumber=7.35`, `nev=20`, `tol=1e-8`, `max_iter=100` |
| `3D/baseflow/` | `control.nml`, `grid3d.dat`; 101 x 101 x 101, Ma=0.2, Re=200, `RunMode=0`, `IF_init=2`, `dt=0.0005`, `End_time=100` |
| `3D/GSA/` | `control.nml`, `grid3d.dat`, `flow3d.dat`, `Eigenvalues.ref.dat`; same mesh and flow parameters, `RunMode=1`, `IF_init=0`, `dt=0.0005`, propagator horizon 0.5, `nev=10`, `tol=1e-8`, `max_iter=100` |

The nonlinear runs start from a quiescent state (`IF_init=2`), so no flow restart is required. GSA runs use the supplied frozen base flows and a random Arnoldi starting vector by default. The 2D GSA calculation includes out-of-plane perturbations at the specified wavenumber; the 3D calculation resolves all three spatial directions. Controls in the inactive namelist group are not used.

## Run

Run each command from the indicated directory relative to `cases/LDC/`. The four MPI ranks are an example; adjust the count for the machine and mesh partition.

| Working directory | Command |
| --- | --- |
| `2D/baseflow/` | `mpiexec -n 4 ../../../../build/openglsa > run.log 2>&1` |
| `2D/GSA/` | `mpiexec -n 4 ../../../../build/openglsa > run.log 2>&1` |
| `3D/baseflow/` | `mpiexec -n 4 ../../../../build/openglsa3d > run.log 2>&1` |
| `3D/GSA/` | `mpiexec -n 4 ../../../../build/openglsa3d > run.log 2>&1` |

NS and GSA are independent workflows: the supplied GSA base flows allow stability analysis without rerunning NS. Check `resid.history` before using a newly computed NS field as a replacement base flow; reaching `End_time` alone does not establish convergence.

For GSA, inspect the convergence count and PARPACK status in `run.log`. Successful extraction writes `Eigenvalues.dat` and `Mode-*.dat`; compare eigenvalues with `Eigenvalues.ref.dat`. The requested eigenvalue count is not a guarantee of convergence. Use a working copy for reruns because histories are appended and matching output files may be overwritten.

## Optional 3D Arnoldi starting vector

From `3D/GSA/`, compile and run the supplied generator:

```bash
gfortran -O2 -o init_vector ../init_vector.f90
./init_vector 101 101 101 1.0
```

It reads `grid3d.dat` and writes (or replaces) `fluc3d.dat`, with deterministic velocity perturbations that vanish at the boundary and zero density and temperature perturbations. The dimensions must match `control.nml`; use the same Fortran unformatted-record conventions as the solver. Set `IF_init=1` in `control_gsa` to use this file as the Arnoldi starting vector. The supplied setting `IF_init=0` does not use it.

## Post-processing

From `2D/GSA/`, convert an exported mode column to Tecplot ASCII:

```bash
../../../../build/op2tec Mode-0001.dat 201 201
```

From `3D/GSA/`, use:

```bash
../../../../build/op2tec3d Mode-0001.dat 101 101 101
```

The same converters work on NS snapshots from the corresponding `baseflow/` directory (`flow-*.dat` in 2D and `flow3d-*.dat` in 3D). Complex modes use paired real/imaginary column files; see the [User Manual](../../docs/USER_MANUAL.md#63-complex-mode-files) for interpretation.
