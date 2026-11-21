# Cylinder

Cylinder flow at Ma=0.1 and Re=100 on a 315 x 301 mesh. The computational i direction is periodic. `basicflow/` uses NS with SFD; `GSA/` uses IRAM about the supplied frozen base flow with `wavenumber=0`.

## Files and settings

| Directory | Inputs and active controls |
| --- | --- |
| `mesh/` | `cylinder.x` (316 x 301 ASCII mesh) and `grid.dat` (315 x 301) |
| `basicflow/` | `control.nml`, `grid.dat`, `flow.dat`; `RunMode=0`, `IF_init=1`, `dt=0.00025`, `End_time=100` |
| `GSA/` | `control.nml`, `grid.dat`, `Eigenvalues.ref.dat`, `flow.dat`; `RunMode=1`, `dt=0.00025`, propagator horizon 1, `nev=2`, `ncv=20`, `tol=1e-8`, `max_iter=50` |

Both workflows use LINEAR5/CD6 and `BC_Cylinder`. SFD is active only in the NS workflow, with `Chi=0.4273` and `Delta=3.2650`. The supplied flow headers start at step/time zero. The controls in the inactive namelist group are not used.

## Run

After building the project, run from this directory:

```bash
cd basicflow
mpiexec -n 4 ../../../build/openglsa > run.log 2>&1
```

This writes `flow-*.dat`, `resid.history`, `resid-sfd.history` and the overwritten residual field `resid.dat`. Inspect both histories when assessing base-flow convergence; reaching `End_time` alone does not establish convergence.

To run stability analysis with the supplied `GSA/flow.dat`, independently of the preceding NS calculation:

```bash
cd ../GSA
mpiexec -n 4 ../../../build/openglsa > run.log 2>&1
```

Starting from the case root instead, use `cd GSA`. Check the convergence count and PARPACK status in `run.log`. Successful extraction writes `Eigenvalues.dat` and `Mode-*.dat`; the requested eigenvalue count is not a guarantee of convergence.

## Post-processing

From `basicflow/`, for example:

```bash
python3 ../../../utils/historyPlot.py resid.history
python3 ../../../utils/historyPlot.py resid-sfd.history
../../../build/op2tec flow-00000000.dat 315 301 omesh
```

Replace the initial snapshot name with a saved snapshot to view the evolved flow. From `GSA/`, convert an exported column with:

```bash
../../../build/op2tec Mode-0001.dat 315 301 omesh
```

`omesh` closes the visualization mesh by repeating its first i-column. Output is Tecplot ASCII (`.tec`). Complex modes use paired real/imaginary column files; see the [manual](../../docs/USER_MANUAL.md#63-complex-mode-files).

The supplied binary meshes can be used directly. To regenerate the mesh from `mesh/`, run `../../../build/x2op cylinder.x omesh`. This drops the duplicate final i-column, giving 315 points. Copy the resulting `grid.dat` into a run directory only when deliberately replacing its mesh; keep the matching flow dimensions.
