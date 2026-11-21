# Compression Ramp

Compressible ramp flow at Ma=7.7, Re=4200 and reference temperature 125 K on a 1601 x 321 mesh. The wall follows the supplied geometry and uses `BC_Flat_Plate`, with `BC_Para(1)=0` and nondimensional wall temperature `BC_Para(2)=2.344`.

## Files and settings

| Directory | Inputs and active controls |
| --- | --- |
| `mesh/` | `RAMP.x` and `grid.dat`, both 1601 x 321 |
| `basicflow/` | `control.nml`, `grid.dat`, `flow.dat`; `RunMode=0`, `IF_init=1`, `dt=0.0005`, `End_time=10` |
| `GSA/` | `control.nml`, `grid.dat`, `Eigenvalues.ref.dat`, an independent supplied `flow.dat`; `RunMode=1`, `dt=0.0005`, propagator horizon 2, `wavenumber=0.95`, `nev=41`, `ncv=101`, `tol=1e-8`, `max_iter=100` |

The nonlinear workflow uses WENOMR5, characteristic LLF reconstruction (`IF_Character=2`), `pp_limiter=1` and CD6. SFD is off by default. In GSA, the input WENOMR5 setting is overridden by LINEAR5, the linearized LLF branch is used, and the nonlinear positivity limiter is not applied.

The supplied `basicflow/flow.dat` has its header reset to step/time zero. `GSA/flow.dat` is a different supplied flow with stored time approximately 10000; GSA resets base-flow time internally, so no header reset is needed. Do not overwrite it with the output of the short `basicflow/` run unless that replacement has been established as the intended converged base flow.

## Run

After building the project, run from this directory:

```bash
cd basicflow
mpiexec -n 4 ../../../build/openglsa > run.log 2>&1
```

Outputs are `flow-*.dat`, appended `resid.history`, overwritten `resid.dat` and `pplim.txt`. The latter contains the latest saved limiter activation map, not a cumulative map or one file per snapshot.

For stability analysis with the supplied base flow:

```bash
cd ../GSA
mpiexec -n 4 ../../../build/openglsa > run.log 2>&1
```

Starting from the case root instead, use `cd GSA`. Check the PARPACK status and actual convergence count in `run.log`. Successful extraction writes `Eigenvalues.dat` and `Mode-*.dat`. The 1601 x 321 mesh and Arnoldi basis require more resources than the smaller examples; select MPI ranks to suit the machine.

## Post-processing

From `basicflow/`:

```bash
python3 ../../../utils/historyPlot.py resid.history
../../../build/op2tec flow-00020000.dat 1601 321
```

`flow-00020000.dat` is the expected final snapshot name for the supplied zero-time restart, `dt=0.0005` and `End_time=10`. From `GSA/`:

```bash
../../../build/op2tec Mode-0001.dat 1601 321
```

The converter writes Tecplot ASCII `.tec` files. Interpret complex modes using the paired real/imaginary columns and frequencies using the actual propagator interval, as described in the [manual](../../docs/USER_MANUAL.md#6-binary-layouts-and-eigenvalue-interpretation).

Mesh regeneration is optional: from `mesh/`, run `../../../build/x2op RAMP.x`. This creates `grid.dat` and `grid2d.tec` but no flow restart. Existing run directories already contain their matching binary meshes.
