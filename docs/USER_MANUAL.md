# OpenGLSA 1.0 user manual

The controls, case preparation, outputs and post-processing here cover the two-dimensional solver. The three-dimensional solver follows the same workflow; [section 8](#8-three-dimensional-solver) describes the differences. See the [README](../README.md#build) for dependencies and compilation.

Use [control.nml.example](../control.nml.example) as a configuration reference. Each case also needs a compatible mesh and, when required, base-flow or restart data. The [included examples](../cases/README.md) contain configurations, data and commands for specific cases.

## 1. Execution model and working directory

The program reads files from its working directory. Create a directory for each calculation, put `control.nml` and the required data there, and run the executable by its path:

```bash
cd /absolute/path/to/my-case
mpiexec -n 4 /absolute/path/to/OpenGLSA/build/openglsa > run.log 2>&1
```

The solver partitions the mesh in two dimensions based on the MPI rank count. Each local block must have at least `max(LAP,5)` points in each computational direction. Using too many ranks on a small mesh can produce blocks below this minimum.

| `RunMode` | Operation | Required input data beyond `control.nml` |
| --- | --- | --- |
| `0` | Nonlinear NS, optionally with SFD | `grid.dat`; also `flow.dat` when `IF_init=1` |
| `1` | IRAM global stability about a frozen base flow | `grid.dat`, `flow.dat`; also `fluc.dat` when `IF_init=1` |
| `2` | Direct LNS about a frozen base flow | `grid.dat`, `flow.dat`; also `fluc.dat` when `IF_init=1` |

LNS/IRAM retain five perturbation components `(d,u,v,w,T)`, including the out-of-plane velocity, and use a prescribed wavenumber. Do not mix flow files, coordinates or nondimensional parameters from different cases. Reruns append to flow history files and can overwrite snapshots.

## 2. Reading the controls

Copy the template to `control.nml`. Fortran namelists start with `&group_name` and end with `/`; `!` introduces a comment. Quote strings, as in this example:

```fortran
&control_common
   RunMode = 0
   IF_Character = 0
   Scheme_Invis = 'LINEAR5'
   Scheme_Vis = 'CD6'
   ! Remaining case controls go here.
/
&control_cfd2d
   IF_init = 0
   dt = 0.001D0
   End_time = 100.0D0
/
```

Use the complete template to set mesh dimensions, boundaries and physical parameters. The solver reads `control_common`, then `control_cfd2d` for mode 0 or `control_gsa2d` for modes 1/2. The template contains all three groups; the solver ignores the inactive mode group. Keep `IF_Character` in `control_common`.

For omitted parameters, `src/Read_Para.f90` supplies the defaults listed below. Some values in the template differ from these defaults. Check that the parameter combinations suit your case even if the parser accepts them.

### 2.1 `control_common`: modes, transport and numerical schemes

| Parameter | Parser default | Meaning and choices |
| --- | --- | --- |
| `RunMode` | `0` | `0` NS; `1` IRAM; `2` LNS. |
| `IF_Viscous` | `1` | `0` inviscid; `1` include viscous terms. Check boundary implementations before choosing inviscid operation; some are written for viscous problems. |
| `IF_Cylindrical_Coordinate` | `0` | `0` Cartesian; `1` cylindrical-coordinate branch. Mesh geometry, boundaries and wavenumber must agree with the selected branch. |
| `IF_Character` | `0` | NS: `0` component-wise flux reconstruction, `1` characteristic reconstruction, `2` characteristic LLF path. LNS/IRAM: `2` selects linearized LLF splitting; other values select the alternative linearized splitting. |
| `pp_limiter` | `0` | `0` off; use `1` to enable the nonlinear positivity-preserving limiter in NS with `IF_Character=2`. It is not applied to LNS/IRAM perturbations. |
| `Scheme_Invis` | `'LINEAR5'` | NS accepts `'LINEAR5'`, `'WENOZ5'`, `'WENOZQ5'`, `'WENOMR5'`. LNS/IRAM force `LINEAR5` after input parsing. |
| `Scheme_Vis` | `'CD6'` | `'CD4'` or `'CD6'` central viscous differences. Boundary closures can be lower order. |

Choose the NS scheme for base-flow preparation; LNS/IRAM always use LINEAR5.

### 2.2 `control_common`: mesh and physical parameters

| Parameter | Parser default | Meaning |
| --- | --- | --- |
| `nx_global` | `256` | Global point count in computational i; must match `grid.dat`. |
| `ny_global` | `256` | Global point count in computational j; must match `grid.dat`. |
| `LAP` | `4` | Halo width. Retain 4 for the supplied stencil configuration; reducing it requires checking all stencil accesses. |
| `Ma` | `0.1D0` | Reference Mach number; positive. |
| `Re` | `100.0D0` | Reference Reynolds number; positive when viscosity is used. |
| `gamma` | `1.4D0` | Specific heat ratio; greater than 1. |
| `Pr` | `0.72D0` | Prandtl number; positive. |
| `T_inf` | `288.0D0` | Reference temperature in kelvin for the temperature-dependent viscosity law; positive. |

The solver uses nondimensional primitive variables and time. Pressure is `p=d*T/(gamma*Ma**2)`, with `Cv=1/(gamma*(gamma-1)*Ma**2)`. For variable viscosity, it uses `S=110.4/T_inf` and `mu=(1/Re)*(1+S)*T**1.5/(S+T)`.

Record the reference length, velocity, density and temperature for each case. Snapshots do not store these dimensional scales. To convert time or frequency to dimensional units, use the case's reference time.

### 2.3 `control_common`: output and domain boundaries

| Parameter | Parser default | Meaning |
| --- | --- | --- |
| `Istep_Show` | `100` | Positive interval in time steps for residual reporting. IRAM overrides it as described below. |
| `Istep_Save` | `100` | Positive interval in time steps for NS/LNS snapshots. IRAM saves modes after eigenvector extraction, rather than time snapshots. |
| `Isymmetric_X` | `0,0` | Symmetry flags on lower/upper computational i sides: `0` off, `1` on. |
| `Isymmetric_Y` | `0,0` | Symmetry flags on lower/upper computational j sides: `0` off, `1` on. |
| `Iperiodic_X` | `0` | Periodicity in computational i: `0` off, `1` on. |
| `Iperiodic_Y` | `0` | Periodicity in computational j: `0` off, `1` on. |
| `Periodic_ISpan` | `0.0D0,1.0D0` | Physical coordinate translation `(Delta x,Delta y)` for i-periodic mesh halos. Set for the actual mesh. |
| `Periodic_JSpan` | `0.0D0,1.0D0` | Physical coordinate translation `(Delta x,Delta y)` for j-periodic mesh halos. Set for the actual mesh. |
| `BoundaryCondition` | `'BC_None'` | Problem-specific boundary keyword; see section 3. |
| `BC_Para` | All 100 entries zero | Boundary-specific real array. Current flat-plate implementation reads entries 1 and 2. |

The X/Y labels on symmetry and periodicity flags refer to the i/j partition directions. On a curved mesh, these may differ from the Cartesian wall normals. Set the periodic coordinate translations to match the mesh and the physical flow boundary conditions; the default vectors may not suit your mesh. Use indexed namelist assignments such as `BC_Para(1)=0.0D0` and `BC_Para(2)=1.0D0` to set individual entries.

### 2.4 `control_cfd2d`: nonlinear NS and SFD

| Parameter | Parser default | Meaning |
| --- | --- | --- |
| `IF_init` | `0` | `0` uniform freestream initialization (u = 1); `1` read `flow.dat`, retaining its step and time; `2` quiescent initialization (u = 0), intended for closed flows such as the lid-driven cavity. |
| `dt` | `0.001D0` | Positive fixed nondimensional time step. Reported CFL values do not automatically adjust it. |
| `End_time` | `10.0D0` | Final target simulation time. On restart this is an absolute target, not an additional duration. Whole steps can pass the target slightly. |
| `IF_SFD` | `0` | `0` off; `1` apply Selective Frequency Damping during NS integration. |
| `Chi` | `0.0D0` | SFD feedback coefficient; choose a physically appropriate nonnegative value. A zero value supplies no feedback. |
| `Delta` | `1.0D0` | SFD filter time scale; strictly positive when SFD is active. |

SFD can prepare a base flow during NS integration. The stopping target is `End_time`; integration does not stop automatically at a residual threshold. On restart, SFD initializes its filter field from the restarted flow; the four-field restart file does not store the filter field.

### 2.5 `control_gsa2d`: direct LNS and IRAM

| Parameter | Parser default | Meaning and mode applicability |
| --- | --- | --- |
| `IF_init` | `0` | LNS: `0` random perturbations, `1` read `fluc.dat`. IRAM: `1` reads `fluc.dat` as the Arnoldi starting vector and stops with an error if the file is missing; other values let PARPACK generate a random starting vector. Both modes always read the base flow from `flow.dat`. |
| `dt` | `0.001D0` | Positive fixed LNS time step, also used for each IRAM propagator application. |
| `End_time` | `10.0D0` | IRAM: requested horizon for each propagator application, starting from zero. LNS: final target time; a perturbation restart retains its stored time. |
| `tol` | `1.0D-6` | IRAM eigenvalue convergence tolerance, not a time-integration error tolerance. Not used to stop direct LNS. |
| `nev` | `10` | Requested IRAM eigenvalue count; must be positive and compatible with the problem dimension. |
| `ncv` | `30` | Requested Arnoldi subspace dimension. The code uses `max(2*nev+1,ncv)`; also check PARPACK dimension constraints. |
| `max_iter` | `10` | Maximum Arnoldi update cycles; not a count of NS/LNS time steps. |
| `wavenumber` | `0.0D0` | Used in both modes. Cartesian: out-of-plane wavenumber in the selected nondimensional coordinates. Cylindrical: azimuthal mode order. |

For cylindrical calculations, set `wavenumber` to a nonnegative integer azimuthal mode order and use the matching axis/symmetry settings.

For IRAM, use a propagator horizon of at least ten time steps. Set `dt`, `Istep_Show` and `Istep_Save` to positive values.

Use `IF_init=0` to generate normalized random LNS perturbations, or supply `fluc.dat` with `IF_init=1` for a prescribed initial perturbation.

In IRAM, `IF_init=1` uses `fluc.dat` as the Arnoldi starting vector. IRAM ignores its header step/time, and the field need not be normalized. You can copy a previously computed `Mode-*.dat` or an LNS `fluc-*.dat` snapshot to `fluc.dat` for this purpose. The starting vector must not be identically zero.

## 3. Problem-specific boundaries

| Keyword | Current behavior and assumptions |
| --- | --- |
| `BC_None` | No problem-specific boundary treatment. Configure other settings, such as periodicity, to provide the required domain treatment. |
| `BC_Flat_Plate` | Flat-plate NS and LNS routines. `BC_Para(1)=x0` marks the start of the no-slip wall on the lower j boundary; `BC_Para(2)=Tw` chooses wall temperature. `Tw<=0` selects the implemented adiabatic extrapolation; `Tw>0` selects an isothermal wall at nondimensional `Tw`. |
| `BC_Lid_Driven_Cavity` | Nonlinear and linearized routines. Upper-j lid velocity is `(1-(2*x-1)**18)**2`, other wall velocities are zero and wall temperature is 1. The lid profile assumes the intended nondimensional x interval. These values are fixed in the code and cannot be set through `BC_Para`. |
| `BC_Sphere` | Implemented geometry-specific wall, outer-field and symmetry/axis treatment. Requires a matching mesh and coordinate convention. |
| `BC_Cylinder` | Implemented geometry-specific treatment. Requires a matching mesh and boundary placement. |
| `BC_User_Def` | Empty NS/LNS extension hooks; implement the needed treatment before use. |

These routines assume specific geometries and boundary locations, so they cannot be assigned to arbitrary mesh faces. For a new geometry, review `src/BoundaryConditions.f90`, `src/Apply_BC.f90` and the associated coordinate and halo handling. The sphere, cylinder and cavity routines do not read `BC_Para`. Check that an existing boundary routine fits your setup before using it.

## 4. Workflows

### 4.1 Nonlinear flow and base-flow preparation (`RunMode=0`)

1. Supply `grid.dat` and matching dimensions, physical parameters and boundaries.
2. Select `IF_init=0` for a uniform start, or put a compatible snapshot in `flow.dat` and select `IF_init=1`.
3. Set `dt`, final `End_time`, output intervals and the NS reconstruction.
4. Enable and configure SFD if needed. Inspect both NS and SFD histories; choose case-appropriate convergence criteria.
5. After checking base-flow convergence, copy the selected binary snapshot to `flow.dat` in a separate LNS/IRAM case directory, together with the same mesh.

To prepare an LNS/IRAM case from a snapshot:

```bash
mkdir -p ../gsa-case
cp grid.dat ../gsa-case/grid.dat
cp flow-00010000.dat ../gsa-case/flow.dat
cp control.nml ../gsa-case/control.nml
# Edit the copied control.nml for RunMode=1 or 2 and the matching controls.
```

This example prepares input files; it is not a supplied benchmark. LNS/IRAM reset the base-flow header time and step internally.

### 4.2 Direct perturbation evolution (`RunMode=2`)

1. Supply `grid.dat` and the frozen base flow `flow.dat`.
2. Set `RunMode=2` and configure `control_gsa2d`.
3. Select random initialization with `IF_init=0`, or supply a five-field `fluc.dat` with `IF_init=1`.
4. Set the wavenumber, time step, final time and output intervals, and run.
5. Inspect `resid.history` and convert `fluc-*.dat` for visualization.

The base flow stays frozen as the perturbations advance. Snapshots contain signed primitive perturbations rather than total nonlinear fields. A perturbation restart retains its step and time. To continue integration, set `End_time` above the restart time.

### 4.3 Global stability (`RunMode=1`)

1. Supply the same mesh and a suitable base flow in `flow.dat`.
2. Set `RunMode=1`, `dt`, propagator horizon `End_time`, wavenumber and IRAM controls.
3. Keep `IF_init=0` for a random Arnoldi starting vector, or supply a five-field `fluc.dat` with `IF_init=1` to start from a prescribed perturbation.
4. Run and inspect `run.log` for the requested/actual Arnoldi dimension, `PDNAUPD`/`PDNEUPD` status and the number of converged Ritz values.
5. Read `Eigenvalues.dat` and convert the required `Mode-*.dat` files.

IRAM repeatedly advances LNS over the same interval and resets perturbation time and step before each application. It solves a real nonsymmetric PARPACK problem with `BMAT='I'`, `WHICH='LM'` and standard operator mode. It does not assemble a global matrix. Check `run.log` for progress and convergence status, then use the extracted eigenvalues and mode files for post-processing.

## 5. Files produced by each mode

### 5.1 NS (`RunMode=0`)

| File | Format | When and meaning |
| --- | --- | --- |
| `flow-NNNNNNNN.dat` | Binary | Written during initialization, every `Istep_Save` steps and at exit; density, two velocity components and temperature. |
| `resid.history` | ASCII, 6 columns | Appended every `Istep_Show` steps: time, four conservative-variable residuals, their mean. |
| `resid-sfd.history` | ASCII, 6 columns | With SFD: time, four flow/filter differences, their mean. |
| `resid.dat` | Tecplot ASCII | Residual field written at snapshot intervals and exit; overwritten, not a numbered series. |
| `pplim.txt` | ASCII | With `IF_Character=2` and `pp_limiter/=0`, accompanies flow writes; overwritten with the current limiter map. |

For NS history residuals, the solver sums absolute conservative-variable changes divided by `dt`, excluding certain physical boundary points. It then divides by the full global point count. These are not RMS residuals, and the solver evaluates them before the SFD filter update. `resid.dat` contains `X Y R_mean R_rho R_rhou R_rhov R_rhoE` and uses the full field of absolute changes. Its averaging and boundary treatment may differ from the history statistic.

SFD differences measure the separation between the flow and filter fields without division by `dt`. They therefore differ from the NS time-change residual. `pplim.txt` records global indices, coordinates and direction-specific density and pressure limiter hit counts from the last RK sub-stage. The counts are not cumulative. Each write replaces the previous snapshot's map.

### 5.2 Direct LNS (`RunMode=2`)

| File | Format | When and meaning |
| --- | --- | --- |
| `fluc-NNNNNNNN.dat` | Binary | Initialization, every `Istep_Save` steps and exit; five primitive perturbation fields. |
| `resid.history` | ASCII, 7 columns | Time, five conservative-perturbation RMS change rates, combined RMS. |

For LNS residuals, the solver sums squared conservative-variable changes divided by `dt**2`, excluding certain boundary points. It divides by the full global point count before taking square roots. These values are not directly comparable to NS absolute-change residuals. Perturbation growth can indicate physical instability, so a growing LNS residual alone cannot serve as a nonlinear steady-state convergence test.

### 5.3 IRAM (`RunMode=1`)

| File | Format | Meaning |
| --- | --- | --- |
| `Eigenvalues.dat` | ASCII, header plus 6 columns | Mode number, propagator magnitude, phase angle, growth, angular frequency, frequency divided by `2*pi`. |
| `Mode-NNNN.dat` | Binary | Exported real PARPACK vector columns, converted to five primitive perturbation fields. Complex modes require pairing as explained below. |

After successful eigenvector extraction, the solver writes these files for the reported converged values. Run in a separate clean directory to keep earlier results out of the current output.

## 6. Binary layouts and eigenvalue interpretation

### 6.1 Binary records

Mesh and snapshot files use Fortran sequential unformatted records. Each global field occupies one record of `real64` values with shape `(nx_global,ny_global)`, i fastest, excluding halos. The files do not store array dimensions. Record markers, endianness and default integer size must match the toolchain that wrote the file. Renaming a binary file to `.tec` does not make it a text file for visualization.

| File type | Records in order |
| --- | --- |
| `grid.dat` | X-coordinate array; Y-coordinate array. |
| `flow.dat`, `flow-*.dat` | Header `(Istep,tt)`; `d`; `u`; `v`; `T`. |
| `fluc.dat`, `fluc-*.dat` | Header `(Istep,tt)`; `d_fluc`; `u_fluc`; `v_fluc`; `w_fluc`; `T_fluc`. |
| `Mode-*.dat` | Header `(mode_index,0.0_real64)`; five primitive vector-column fields in the same order as perturbations. |

Each header is a single mixed-type record: a default integer followed by `real64`. The supplied reset and interpolation tools assume a 4-byte integer and compatible record markers. MPI helpers assemble and distribute mesh and flow records through one global file shared by all ranks.

### 6.2 Growth and frequency

For a propagator eigenvalue `mu` and actual interval `tau=Istep*dt`, the code writes:

```text
Magnitude       = abs(mu)
Angle           = atan2(imag(mu),real(mu))
Growth          = log(abs(mu))/tau
Frequency       = Angle/tau
Strouhal Number = Angle/(2*pi*tau)
```

`Frequency` is angular frequency. The last column gives cyclic frequency in the chosen nondimensional time units. Interpreting it as a physical Strouhal number requires compatible reference scaling. Phase is the principal angle, so the reported frequency depends on this logarithm branch.

The solver prints the actual `tau`, which can differ from the requested horizon because the time step is fixed. The stored magnitudes and phases describe the propagator; they are not directly the continuous-time growth-rate eigenvalues.

### 6.3 Complex mode files

PARPACK's real nonsymmetric extraction stores a complex eigenvector in two consecutive real columns: the real and imaginary parts for the positive-imaginary eigenvalue. Its conjugate uses the same real part and the negative imaginary part. The writer exports each real column to a separate `Mode-NNNN.dat` file. Each file therefore holds one real column, not a full complex vector.

Pair the columns using the conjugate eigenvalues and PARPACK ordering before computing complex amplitude or phase, or reconstructing an oscillatory mode. `op2tec` converts one real column at a time and does not pair them.

Vector normalization and phase are arbitrary, so an exported mode does not prescribe the amplitude of a nonlinear disturbance. The storage convention follows the [arpack-ng 3.9.1 PDNEUPD interface](https://github.com/opencollab/arpack-ng/blob/3.9.1/PARPACK/SRC/MPI/pdneupd.f).

## 7. Post-processing

### 7.1 Mesh conversion with `x2op`

From a directory containing an ASCII structured mesh:

```bash
/absolute/path/to/build/x2op mesh.xyz
```

The ASCII layout is `nx ny`, followed by all X coordinates and then all Y coordinates, with i varying fastest. The utility writes binary `grid.dat` and ASCII `grid2d.tec`. The ASCII output includes uniform placeholder density, velocity and temperature for mesh preview only. The utility does not produce `flow.dat` and can overwrite existing output files.

For an O-mesh whose last i-column duplicates the first:

```bash
/absolute/path/to/build/x2op mesh.xyz omesh
```

The utility drops the duplicate final column. Set `nx_global` to the resulting point count, excluding that column. Configure the solver's periodicity and translations separately.

### 7.2 Tecplot ASCII conversion with `op2tec`

Run from the directory containing `grid.dat`:

```bash
/absolute/path/to/build/op2tec flow-00010000.dat 361 641
/absolute/path/to/build/op2tec fluc-00010000.dat 361 641
/absolute/path/to/build/op2tec Mode-0001.dat 361 641
```

Replace the example dimensions with those of the binary mesh. Each command writes an ASCII `.tec` file beside its input snapshot, for example `flow-00010000.tec`. The file contains a structured `I/J` zone with `DATAPACKING=POINT`, `SOLUTIONTIME` and i-fastest rows. Conversion does not require TecIO. Load the `.tec` file as a Tecplot ASCII dataset. The variables are:

- Flow: `X Y Density U V Temperature`.
- Perturbations/mode columns: `X Y Density U V W Temperature`.

For perturbations, these labels refer to signed perturbation quantities. The converter does not export pressure. For a total flow, reconstruct it as `p=d*T/(gamma*Ma**2)`. For perturbation pressure, use the linearized relation and the base flow; applying the nonlinear formula to perturbations alone gives the wrong quantity.

For an O-mesh, append `omesh` to duplicate the first i-column in the visualization:

```bash
/absolute/path/to/build/op2tec flow-00010000.dat 361 641 omesh
```

The converter identifies five-field input by the basename prefix `fluc`/`FLUC` or `Mode`/`MODE`. If you rename a five-field file with an unrelated prefix, it will be read as four-field flow data. This rule also applies to snapshots supplied by path; the converter always reads `grid.dat` from the working directory.

Mode files have zero time in their header and do not form a temporal animation series. Conversion leaves the binary input unchanged but can overwrite an existing `.tec`.

To convert a batch of snapshots:

```bash
for snapshot in flow-*.dat; do
    [ -f "$snapshot" ] || continue
    /absolute/path/to/build/op2tec "$snapshot" 361 641 || break
done
```

For 3D snapshots, run `op2tec3d` from the directory containing `grid3d.dat`:

```bash
/absolute/path/to/build/op2tec3d flow3d-00010000.dat 101 81 41
/absolute/path/to/build/op2tec3d Mode-0001.dat 101 81 41 omesh
```

Use the actual mesh dimensions. All 3D snapshots contain five fields (`d,u,v,w,T`), regardless of filename. The converter writes Tecplot ASCII `.tec` files with variables `X Y Z Density U V W Temperature`, a structured `I/J/K` zone, `DATAPACKING=POINT`, `STRANDID=1` and the snapshot's `SOLUTIONTIME`. Rows vary fastest in i, then j, then k.

The optional `omesh` argument appends a copy of the first i-column, including coordinates and all fields. Conversion does not require TecIO. It leaves the binary input unchanged and can overwrite an existing output `.tec`.

### 7.3 Residual-history plotting

Install NumPy and Matplotlib as described in the README, then run:

```bash
python3 /absolute/path/to/OpenGLSA/utils/historyPlot.py resid.history --refresh 10
python3 /absolute/path/to/OpenGLSA/utils/historyPlot.py resid-sfd.history --refresh 10
```

The optional first argument defaults to `resid.history`. `--refresh` sets the refresh interval in seconds; `--linewidth` sets the line width. Both values must be positive. The utility plots six- and seven-column histories against simulation time on a logarithmic vertical scale. It omits nonpositive and nonfinite values.

The interactive window requires a graphical plotting backend. Close the window to stop. IRAM does not write a direct-LNS history for this utility.

### 7.4 Resetting restart headers

```bash
cp flow.dat flow.dat.backup
/absolute/path/to/build/reset_step flow.dat
```

`reset_step FILE [MORE_FILES ...]` sets the step and time to zero in place and leaves field records unchanged. It expects 4-byte record markers and a 12-byte header payload (4-byte integer plus 8-byte real). Check its diagnostics for files skipped because of an unexpected layout. Do not apply it to `grid.dat`.

NS/LNS continuation normally retains the original header. Reset a working copy only to start from a new time origin. The utility does not rename files. LNS/IRAM reset frozen base-flow headers internally.

### 7.5 Interpolating a base flow onto another mesh

`utils/interp_flow.py` requires NumPy, Matplotlib and SciPy. Its defaults read `grid.dat.old`, `flow.dat.old` and the new `grid.dat`, then write `flow.dat`. For explicit paths and sizes:

```bash
python3 /absolute/path/to/OpenGLSA/utils/interp_flow.py \
  --old-grid grid.dat.old --old-flow flow.dat.old \
  --new-grid grid.dat --new-flow flow.dat \
  --old-shape 361x641 --new-shape 721x641 --method linear
```

Other options include `--control` (default `control.nml`), `--method cubic`, `--istep` (default 0), `--time` (default 0) and `--force` to overwrite an existing output. The script can read new dimensions from `control.nml`; otherwise it may guess them. Specify shapes explicitly for reproducible results.

The utility handles four-field base flows. It does not support five-field perturbations or complex eigenmodes. It interpolates on triangles formed from the original structured mesh cells and reports projections for points outside the old mesh. Inspect the projection distances, density, temperature and wall values. Before running a stability analysis, check the base flow's consistency and converge it again on the new mesh.

## 8. Three-dimensional solver

Run the 3D solver `openglsa3d` in the same way as the 2D solver:

```bash
mpiexec -n 8 /absolute/path/to/OpenGLSA/build/openglsa3d > run.log 2>&1
```

Start from [control.nml.example.3d](../control.nml.example.3d). Mode groups are named `control_cfd` and `control_gsa` instead of `control_cfd2d` and `control_gsa2d`. `control_common` adds `nz_global`, `Isymmetric_Z`, `Iperiodic_Z` and `Periodic_KSpan`; all span vectors have three components `(Delta x,Delta y,Delta z)`.

There is no `wavenumber` or `IF_Cylindrical_Coordinate`. Supported boundaries are `BC_None`, `BC_Flat_Plate`, `BC_Lid_Driven_Cavity` and `BC_User_Def`. MPI ranks are split in all three directions.

Input and snapshot filenames include `3d`; their contents have the same meaning as in 2D.

| 2D | 3D |
| --- | --- |
| `grid.dat` | `grid3d.dat` (X, Y, Z) |
| `flow.dat`, `flow-*.dat` | `flow3d.dat`, `flow3d-*.dat` |
| `fluc.dat`, `fluc-*.dat` | `fluc3d.dat`, `fluc3d-*.dat` |
| `Mode-*.dat`, `Eigenvalues.dat` | unchanged |

`IF_init=1` reads `flow3d.dat` in NS, and `fluc3d.dat` in LNS and IRAM. All 3D flow files hold five fields `d,u,v,w,T`, so `resid.history` and `resid-sfd.history` have seven columns. In binary files each field is written as `nz_global` records, one `(nx_global,ny_global)` plane per k.

`x2op3d mesh.xyz [omesh]` reads ASCII `nx ny nz` followed by all X, Y, Z coordinates, and writes `grid3d.dat`, a uniform `flow3d.dat` and `grid3d.tec`. Use `op2tec3d` (section 7.2) to convert snapshots and modes. `interp_flow.py` handles 2D flows only.
