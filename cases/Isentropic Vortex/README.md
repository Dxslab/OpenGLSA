# Isentropic Vortex

Periodic inviscid vortex transport on a deformed 128 x 128 mesh spanning 10 x 10. The supplied controls use `RunMode=0`, characteristic WENOMR5, `Ma=0.84515425472`, `dt=0.005` and `End_time=100`.

## Run

After building OpenGLSA, run from this directory using the supplied `grid.dat` and `flow.dat`:

```bash
mpiexec -n 4 ../../build/openglsa > run.log 2>&1
```

To regenerate the initial data first, compile with a compatible Fortran compiler:

```bash
gfortran -ffree-line-length-none initVortex.f90 -o initVortex
./initVortex
```

This overwrites `grid.dat` and `flow.dat` and creates the ASCII preview `vortex.tec`. Keep the generator's `nx/ny` and `gamma/Minf` consistent with `control.nml`; the mesh is periodic in both directions.

## Output and visualization

The solver writes `flow-*.dat`, `resid.history` and overwritten `resid.dat`. With the supplied zero-time restart, the expected final snapshot is `flow-00020000.dat`. Convert it and plot the history with:

```bash
../../build/op2tec flow-00020000.dat 128 128
python3 ../../utils/historyPlot.py resid.history
```

Load `flow-00020000.tec` in Tecplot. It contains primitive variables; the generator's `vortex.tec` instead contains density, momenta and total energy. Use corresponding quantities when comparing the initial and transported fields.
