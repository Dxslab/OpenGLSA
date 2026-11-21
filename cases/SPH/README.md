# Sphere

Sphere flow at Ma=0.1 and Re=213 on a 241 x 361 mesh, using the cylindrical-coordinate formulation and `BC_Sphere`. The two i sides have symmetry flags enabled. The GSA calculation uses azimuthal mode order `wavenumber=1`.

## Files and settings

| Directory | Inputs and active controls |
| --- | --- |
| `mesh/` | `sphere.x` and `grid.dat`, both 241 x 361 |
| `basicflow/` | `control.nml`, `grid.dat`, `flow.dat`; `RunMode=0`, `IF_init=1`, `dt=0.0001`, `End_time=100` |
| `GSA/` | `control.nml`, `grid.dat`, `Eigenvalues.ref.dat`, `flow.dat`; `RunMode=1`, `dt=0.0001`, propagator horizon 1, `nev=1`, `ncv=40`, `tol=1e-8`, `max_iter=20` |

Both workflows use LINEAR5/CD6 and temperature-dependent viscosity with reference temperature 288.15 K. SFD is off by default. The supplied flow headers start at step/time zero. The `control_gsa2d` settings in the basic-flow file are inactive while `RunMode=0`; use `GSA/control.nml` for the stated stability configuration.

## Run

After building the project, run from this directory:

```bash
cd basicflow
mpiexec -n 4 ../../../build/openglsa > run.log 2>&1
```

Outputs are `flow-*.dat`, appended `resid.history` and overwritten `resid.dat`. Inspect the flow and residuals before treating an evolved snapshot as a new base flow; termination at `End_time` is not a convergence criterion.

For stability analysis with the supplied `GSA/flow.dat`:

```bash
cd ../GSA
mpiexec -n 4 ../../../build/openglsa > run.log 2>&1
```

Starting from the case root instead, use `cd GSA`. Check PARPACK status and the actual number of converged eigenvalues in `run.log`. Successful extraction writes `Eigenvalues.dat` and `Mode-*.dat`.

## Post-processing

From `basicflow/`:

```bash
python3 ../../../utils/historyPlot.py resid.history
../../../build/op2tec flow-00000000.dat 241 361
```

Replace the initial snapshot name with a later saved snapshot as needed. From `GSA/`:

```bash
../../../build/op2tec Mode-0001.dat 241 361
```

The converter writes Tecplot ASCII `.tec` files. Complex modes require paired real/imaginary column interpretation; see the [manual](../../docs/USER_MANUAL.md#63-complex-mode-files).

Mesh regeneration is optional: from `mesh/`, run `../../../build/x2op sphere.x`. The existing run directories already contain binary meshes; preserve their matching dimensions and coordinate conventions.
