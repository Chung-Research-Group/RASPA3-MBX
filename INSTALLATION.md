# RASPA 3 – Installation Guide

RASPA 3 is a modern Monte‑Carlo simulation package that uses C++ modules. The source build selects C++23 on Windows and C++26 on Linux and macOS.
Pre‑built binaries are **statically linked** and therefore do not require any external runtime libraries.
If you only want to *use* RASPA 3, download a pre‑built package.
If you want to *develop* or *debug* RASPA 3, follow the build‑from‑source instructions.

Use a compiler and standard-library combination that provides the C++ `std` module.

---

## 1. Quick Start (Pre‑built Packages)

All official releases are published on the **GitHub “Releases” tab**.
Choose the file that matches your operating system, CPU architecture and AVX capability.

| Platform        | File suffix | Example                                           |
| --------------- | ----------- | ------------------------------------------------- |
| Conda           | –           | install with `conda install raspa3 raspalib`      |
| Debian / Ubuntu | `.deb`      | `raspa‑{VERSION}‑{OS}‑{ARCH}‑{AVX}.deb`      |
| RPM‑based Linux | `.rpm`      | `raspa‑{VERSION}‑{OS}‑{ARCH}‑{AVX}.rpm` |
| macOS           | `.pkg`      | `raspa‑{VERSION}‑mac‑{ARCH}.pkg`                      |
| Windows         | `.exe`      | `raspa‑{VERSION}‑windows‑{ARCH}‑{AVX}.exe`       |

### 1.1 Conda (all platforms)

```bash
conda install -c conda-forge raspa3 raspalib
```

### 1.2 Debian / Ubuntu (.deb)

```bash
apt install ./raspa3-<VERSION>-ubuntu<20|22|24>-<ARCH>-<AVX>.deb
```

### 1.3 RPM‑based Linux (.rpm)

```bash
rpm -i raspa3-<VERSION>-<OS>-<ARCH>-<AVX>.rpm
```

*To unpack without installing:* `rpm2cpio raspa3-*.rpm | cpio -idmv`

### 1.4 macOS (.pkg)

Pick either the arm64 version (Apple silicon) or the x86_64 version (intel chips). Double‑click the `.pkg` and follow the wizard.

### 1.5 Windows (.exe)

Run the `.exe`; the wizard adds `%RASPA3_HOME%\bin` to your `%PATH%`.

### 1.6 Choosing AVX2 vs AVX‑512

| Binary         | Use when                                             |
| -------------- | ---------------------------------------------------- |
| core‑avx2      | Intel ≥ Haswell (2013) or any AMD ≥ Excavator (2015) |
| skylake‑avx512 | Intel CPUs with AVX‑512 (≈ 25 % faster)              |

---

## 2. Building from Source

### 2.1 Prerequisites

| Package                | Min. version |
| ---------------------- | ------------ |
| clang/clang++          | 18           |
| cmake                  | 4.3.0        |
| ninja                  | 1.11         |
| python                 | 3.11         |
| pybind11               | 2.12         |
| BLAS + LAPACK          | 64‑bit ints  |
| HDF5                   | 1.12         |
| OpenMP · OpenCL        | latest       |

### 2.2 Build with Conda (recommended)

```bash
git clone https://github.com/iRASPA/RASPA3.git
cd RASPA3
conda env create -f env.yml
conda activate raspa3
cmake --preset=linux_conda  # or mac_conda / windows_conda_raspa3
ninja -C build
ninja -C build install      # optional
```

**Note**: remove the package ocl-icd from env.yml for Mac and Windows installations, as OpenCL is only available for Linux. For MacOS using the system OpenCL is recommended.

**Note**: it is often necessary to add the conda installed libraries raspa3 links against to your dynamic library path using `export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH`.

### Optional MBX support

`BUILD_MBX` defaults to `OFF`. A standard build needs no MBX headers or library. For many-body Monte Carlo, install [MBX](https://github.com/paesanilab/MBX) with the same C++ standard library and ABI used by RASPA. The integration uses the public MBX 1.4.0 API, double-precision FFTW and OpenMP. For the Linux Conda toolchain, configure and build in a separate directory.

```sh
cmake --preset linux_conda -B build-mbx -DBUILD_APP=ON -DBUILD_MBX=ON -DMBX_ROOT=/path/to/installed/mbx
cmake --build build-mbx --target raspa3
```

With another RASPA toolchain, pass `-DBUILD_MBX=ON` and `-DMBX_ROOT=/path/to/installed/mbx` to its usual CMake configuration. The MBX prefix should contain `include/bblock/system.h`, `include/mbx_version.h` and `lib/libmbx` with the platform's library suffix.

In `simulation.json`, set `"UseMBX": true` and `"MBXSettingsFile": "mbx.json"` on the chosen system. Other systems use the standard RASPA force field by default. The [CO₂ adsorption example](examples/advanced/8_mbx_co2_in_mfi) includes both inputs.

The Monte Carlo target includes MBX two-, three- and four-body interactions, dispersion, permanent electrostatics relative to the bare host, and induction. RASPA supplies the short-range framework–guest potential, its tail correction and the external field. Rigid CO₂ uses C,O,O site order and water uses O,H,H. Supported moves include translation, rotation, conventional and CBMC exchange, CBMC reinsertion, Widom insertion, and isotropic volume moves in a host-free box without an external field. The upstream single-box transition-matrix driver can use the same fixed-volume exchange moves. The final Monte Carlo JSON records the seven MBX terms in kcal/mol, including the one-body energy for reference. Reported pressure tensors use the classical virial. Binary restarts preserve the MBX settings from the checkpoint and require the same `UseMBX` choice.

### 2.3 Linux examples

#### Ubuntu 24.04

```bash
sudo apt install -y git ca-certificates cmake ninja-build \
  llvm lld clang clang-tidy libc++-dev libc++abi-dev \
  libomp-dev libclang-rt-dev python3 python3-dev pybind11-dev \
  liblapack64-dev libblas64-dev

git clone https://github.com/iRASPA/RASPA3.git
cd RASPA3
cmake -B build --preset=linux-x86_64-core-avx2-ubuntu-24
ninja -C build
ninja -C build install
```

#### Fedora 40

```bash
sudo dnf install -y git llvm lld cmake clang clang-tools-extra ninja-build \
  libomp-devel libcxx libcxx-devel libcxxabi libcxxabi-devel \
  lapack64 blas64 python3 python3-devel pybind11-devel

git clone https://github.com/iRASPA/RASPA3.git
cd RASPA3
cmake -B build --preset=linux-x86_64-core-avx2-fedora-40
ninja -C build
ninja -C build install
```

### 2.4 macOS (Homebrew)

```zsh
brew install llvm lld libomp hdf5 libaec ninja cmake doxygen graphviz lapack pybind11
export PATH="$(brew --prefix llvm)/bin:$PATH"  # clang‑18
export HDF5_ROOT="$(brew --prefix hdf5)"
git clone https://github.com/iRASPA/RASPA3.git
cd RASPA3
cmake -B build --preset=macos-apple-silicon  # or macos-x64-core-avx2
ninja -C build
ninja -C build install
```

---

## 3. CMake Presets

Run `cmake --list-presets` to see every preset.
Common shortcuts:

| Preset                                 | Target                |
| -------------------------------------- | --------------------- |
| linux\_conda                           | Generic Linux + Conda |
| linux-x86\_64-core-avx2-ubuntu-24      | Ubuntu 24, AVX2       |
| linux-x86\_64-skylake-avx512-fedora-40 | Fedora 40, AVX‑512    |
| mac\_conda                             | macOS + Conda         |
| macos-apple-silicon                    | macOS Apple Silicon   |
| windows\_conda\_raspa3                 | Windows + Conda       |

---

Need help?  Open an issue on [GitHub](https://www.github.com/iRASPA/RASPA3/issues) or ask on the discussion forum.
