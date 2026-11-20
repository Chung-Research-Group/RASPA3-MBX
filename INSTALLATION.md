# Install standalone 2D-TMMC on Linux x86-64

This branch adds component-resolved particle-number windows to the conventional
force-field RASPA3 v3.0.30 executable. For the combined rigid-MBX implementation,
use the [mbx branch guide](https://github.com/Chung-Research-Group/RASPA3-MBX/blob/mbx/INSTALLATION.md).
The published upstream packages contain the upstream release.

## Build and install

```bash
git clone --branch 2d-tmmc https://github.com/Chung-Research-Group/RASPA3-MBX.git
cd RASPA3-MBX
conda env create -f env-linux-2d-tmmc.yml
conda activate raspa-2d-tmmc
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cmake --version
cmake --preset raspa-2d-tmmc
cmake --build build-2d-tmmc --target raspa3 --parallel 4
export RASPA_DIR="$HOME/.local/raspa-2d-tmmc"
cmake --install build-2d-tmmc --prefix "$RASPA_DIR" --component app
cmake --install build-2d-tmmc --prefix "$RASPA_DIR" --component examples
```

The Linux preset uses CMake 4.3.x, Clang 20.1.6, GCC 15 libstdc++, Ninja,
LLVM OpenMP, LP64 OpenBLAS, zlib and an OpenCL loader. It builds the native
C++26 executable. Initial configuration fetches HDF5 2.1.1 and Voro++ from
GitHub. For an offline build, supply local source checkouts to CMake with
`-DFETCHCONTENT_SOURCE_DIR_HDF5=/path/to/hdf5` and
`-DFETCHCONTENT_SOURCE_DIR_VORO++=/path/to/voro`. Use the revisions in
the top-level CMakeLists.txt. The application and this preset have no MBX
dependency. The build-tree executable is `build-2d-tmmc/app/raspa3`.

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

For further checks, run the executable regression cases from the checkout with
`python tests/native_tmmc.py --binary "$RASPA_DIR/bin/raspa3" --source . --output "$work/native"`.
The [input reference](docs/manual/commands.md#joint-particle-number-windows-for-transition-matrix-monte-carlo)
explains how system and component macrostate bounds interact.
