# Install RASPA3–MBX on Linux x86-64

Use the `mbx` or `main` branch for rigid MBX and two-component TMMC. For the
standalone contribution based on RASPA3 v3.0.30, use the
[2d-tmmc installation guide](https://github.com/Chung-Research-Group/RASPA3-MBX/blob/2d-tmmc/INSTALLATION.md).
Official iRASPA packages follow the [upstream installation instructions](https://github.com/iRASPA/RASPA3/blob/main/INSTALLATION.md).

## Prepare the environment

```bash
git clone --branch mbx https://github.com/Chung-Research-Group/RASPA3-MBX.git
cd RASPA3-MBX
conda env create -f env.yml
conda activate raspa3-mbx
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cmake --version
```

The environment pins CMake 4.0.3, Clang 20.1.6, GCC 15 libstdc++, Ninja,
LLVM OpenMP and LP64 OpenBLAS. The Linux sources use C++26 standard-library
modules. The preset uses the environment’s HDF5, zlib and OpenCL loader.
OpenCL calculations additionally require a device implementation on the host.

## Build public MBX for the rigid-guest model

For a force-field-only installation, go directly to **Build and install RASPA3**.
For rigid MBX, build the pinned public MBX 1.4.0 library in the active
environment. If it is already installed, set `MBX_ROOT` to its installation
prefix. Use matching headers, library and LLVM OpenMP runtime.

```bash
export MBX_ROOT="$HOME/.local/mbx-1.4.0"
git clone https://github.com/paesanilab/MBX.git ../MBX-1.4.0
(
  cd ../MBX-1.4.0
  git checkout --detach 0e01b75b47611d7d51f27a34b112bdc5e2090a50
  export FFTW_HOME="$CONDA_PREFIX"
  autoreconf -fi &&
    CC="$CONDA_PREFIX/bin/clang" CXX="$CONDA_PREFIX/bin/clang++" \
    LDFLAGS="-L$CONDA_PREFIX/lib" ./configure --prefix="$MBX_ROOT" &&
    make -j4 && make check && make install
)
ldd "$MBX_ROOT/bin/single_point" | grep -E 'fftw|[ig]?omp'
```

The dependency suite contains 34 tests. The dependency and application should
resolve the same LLVM OpenMP runtime. The [public MBX audit](docs/mbx-public-migration.md)
records the version, model and licensing details.

## Build and install RASPA3

From the RASPA3-MBX checkout, configure the application you want.

**Rigid MBX and force-field calculations**

```bash
cmake --preset raspa-mbx
build=build-mbx-app
```

**Force-field calculations only** (no MBX dependency)

```bash
cmake --preset raspa-classical
build=build-classical
```

**Install the selected build**

```bash
export RASPA_DIR="$HOME/.local/$build"
cmake --build "$build" --target raspa3 --parallel 4
cmake --install "$build" --prefix "$RASPA_DIR" --component app
cmake --install "$build" --prefix "$RASPA_DIR" --component examples
```

## Verify the installed executable

In each new shell, reactivate the Conda environment and set `RASPA_DIR` to the application prefix used above. The executable
reads `simulation.json` from the working directory. Copy the installed example
into a writable directory so different users can run it independently.

```bash
export PATH="$RASPA_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
work=$(mktemp -d)
cp -R "$RASPA_DIR/share/raspa3/examples/advanced/2d_nvt_w_calf20" "$work/calf20"
(cd "$work/calf20" && OMP_NUM_THREADS=1 ./run)
cat "$work/calf20/tmmc/tmmc_statistics.txt"
```

The force-field run fixes four waters and four CO₂ molecules in CALF-20. Its
single production row should have populations `4 4`, four positive outward
exchange weights, and five collection weights that sum to the reported trial
count within floating-point precision.

For an MBX-enabled install, also run the supplied rigid CO₂ snapshot calculation.
This exercises the MBX backend directly.

```bash
cp -R "$RASPA_DIR/share/raspa3/examples/raspa3_mbx/3_snapshot_energy_evaluation" "$work/snapshot"
(cd "$work/snapshot" && ./run.sh)
```

## Additional build profiles

The tracked application presets work without a user preset file. For a focused
MBX test build, start with the supplied template in a fresh checkout. If you
already have `CMakeUserPresets.json`, keep your local settings and add the
profiles you need from the template.

```bash
cp -n CMakeUserPresets.json.example CMakeUserPresets.json
cmake --preset mbx-local
cmake --build --preset mbx-local --target raspa3 unit_tests_mbx --parallel 4
ctest --preset mbx-local --output-on-failure
```

The template also provides `mbx-server` (portable x86-64) and
`mbx-release-avx2` (AVX2 CPUs). For multi-user installations, point the
[modulefile template](packaging/modulefiles/raspa3-mbx.lua.example) at the
application prefix and the dependency environment. The [runtime instructions](README.md#install-and-expose-raspa3-mbx-as-a-server-module)
include a Slurm example. Changing compiler, CMake module version or MBX headers
requires a fresh build directory.

The combined CMake configuration recognizes the experimental `import std` gates
from 4.0.3 through 4.4.x. If the compiler reports missing standard modules,
check `command -v cmake`, `cmake --version` and the active environment. If MBX
lives outside the usual include/lib layout, pass `-DMBX_INCLUDE_DIR=...`,
`-DMBX_CONFIG_INCLUDE_DIR=...` and `-DMBX_LIBRARY=...` on the configure command.
Use headers and the library from the same installation. MBX licensing terms
are supplied with the public MBX source.
