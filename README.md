# OpenGLSA 1.0

OpenGLSA is an MPI-parallel Fortran solver for compressible Navier-Stokes flow on structured meshes. It supports nonlinear time integration with optional Selective Frequency Damping (SFD) and direct linearized Navier-Stokes (LNS) integration. For global stability analysis, it uses matrix-free IRAM with PARPACK.

See the [User Manual](docs/USER_MANUAL.md) for controls, workflows, input/output formats, boundary conditions, eigenvalue interpretation and post-processing. The [examples](cases/README.md) include inputs and run instructions for each case.

## Build

OpenGLSA requires CMake 3.18+, a Fortran compiler, MPI with Fortran support, and arpack-ng built with PARPACK. Follow the arpack-ng [installation instructions](https://github.com/opencollab/arpack-ng#readme) to prepare its build environment.

Load the compiler, MPI and library environment, and set `CC`, `CXX` and `FC`
for that toolchain. Use the same environment for dependencies and OpenGLSA.
The dependency script accepts `gnu` or `intel`:

```bash
bash install_deps.sh gnu
cmake -S . -B build
cmake --build build --parallel 4
```

The script installs dependencies in `.deps/`. The build produces `openglsa`,
`openglsa3d`, `x2op`, `x2op3d`, `reset_step`, `op2tec`, `op2tec3d`, `op2plt` and
`op2plt3d`.

The script builds the serial TecIO static library from
[`Izumiko/tecio`](https://github.com/Izumiko/tecio) at commit
`732bcc7f298a2c38b70b65865a2b2db488528649`, using `teciosrc/`.
TecIO requires a C++14 compiler and Boost headers >= 1.69.0. It does not need
compiled Boost libraries. The script downloads Boost 1.69.0 and installs the
headers in `.deps/include`. Downloads and builds stay in `.deps-build/`, and
repeat runs reuse the archives. Set `JOBS` to control TecIO build parallelism
(default: 4).

The tools `op2plt` and `op2plt3d` write binary Tecplot `.plt` files:

```bash
# Run from the case directory containing grid.dat or grid3d.dat.
/path/to/build/op2plt flow.dat  nx ny       # 2D
/path/to/build/op2plt3d flow3d.dat nx ny nz  # 3D
```

Replace the dimensions with the grid sizes. Add `omesh` as the final argument
to duplicate the first i-plane and close an O-mesh. For filenames beginning
with `fluc` or `Mode` (also `FLUC`/`MODE`), the 2D tool reads five fields.
For other filenames, it reads four. The 3D tool always reads five fields,
slab by slab.

Both tools write the output beside the input, replacing its extension with `.plt`
and overwriting any existing output. Only these two tools link to TecIO.
`op2tec` and `op2tec3d` write ASCII `.tec` files.

To build without TecIO or a C++ compiler, configure with
`cmake -S . -B build -DBUILD_TECIO_UTILS=OFF`. This skips both binary converters.
The solver still requires ARPACK/PARPACK. Dependency licenses are in
`.deps/share/licenses/`.

## Package contents

| Path | Contents |
| --- | --- |
| `src/` | Source files |
| `utils/` | Source files for pre- and post-processing tools |
| `cases/` | Examples |
| `control.nml.example` | Runtime configuration template |
| `install_deps.sh` | Dependency installation script |
| `docs/USER_MANUAL.md` | Control options, workflows, outputs and post-processing |
| `LICENSE`, `NOTICE` | Project license and upstream notices |

## License

Copyright (c) 2026 Ruiqing Yang and Baofeng Ma.

OpenGLSA uses [Apache-2.0](LICENSE) and includes portions adapted from [OpenCFD-SC](https://github.com/OpenCFD-IMECH/OpenCFD-SC), Copyright 2021 Xinliang Li. Preserve [LICENSE](LICENSE), [NOTICE](NOTICE) and applicable source notices when redistributing. External dependencies retain their own licenses.
