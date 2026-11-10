# mmkit — Formula, Compound, and Adduct

A C++17 molecular formula library supporting element counts, charge, parsing, Hill notation, exact mass, arithmetic, integer scaling, and containment checks. Negative element counts are supported. Masses exclude electron mass corrections.

## Install a specific release

Versions use `MAJOR.MINOR.PATCH` (for example `1.0.0`), with Git tags such as `v1.0.0`. Install from the assets of the selected [GitHub Release](https://github.com/MotohiroOGAWA/mmkit/releases), using a fixed version rather than a latest-release URL. The commands below are examples for version `1.0.0` on Linux x86_64 and require that release's assets to have been published.

| Language | Install from release `v1.0.0` | Guide |
| --- | --- | --- |
| Python | `python -m pip install "https://github.com/MotohiroOGAWA/mmkit/releases/download/v1.0.0/mmkit-1.0.0.tar.gz"` | [Python](docs/python.md) |
| Java | `curl -fLO "https://github.com/MotohiroOGAWA/mmkit/releases/download/v1.0.0/mmkit-1.0.0-java-linux-x86_64.tar.gz"` then extract and add the JAR and JNI library to the application | [Java](docs/java.md) |
| C# | `curl -fL "https://github.com/MotohiroOGAWA/mmkit/releases/download/v1.0.0/Mmkit.1.0.0.linux-x86_64.nupkg" -o Mmkit.1.0.0.nupkg` then `dotnet add package Mmkit --version 1.0.0 --source "$PWD"` in the application directory | [C#](docs/csharp.md) |
| C++ | `curl -fLO "https://github.com/MotohiroOGAWA/mmkit/releases/download/v1.0.0/mmkit-1.0.0-cpp-linux-x86_64.tar.gz"` then extract and use its CMake package | [C++ installation](#c-installation) |

Python's source archive builds locally and needs the C++ build dependencies. A matching release wheel can instead be installed directly by its asset URL. All native packages require compatible RDKit shared libraries. Choose assets for your OS, architecture, and Python ABI; these commands do not install from PyPI, NuGet.org, or Maven Central.

To prepare a new release, run `./release.sh` and select `major`, `minor`, or `patch`. Use `./release.sh patch --dry-run` to preview the next version. With authenticated GitHub CLI and a clean committed working tree, the script updates versions, tests all bindings, builds all language assets, commits the version changes, pushes the commit and tag, and publishes the GitHub Release. Use `--local` to build without publishing or `--publish-only` to publish previously generated assets.

For another release, change both the tag and every filename/package version. See [RELEASE.md](RELEASE.md) for preparing and publishing all language artifacts.

## C++ installation

Requires the RDKit C++ headers and libraries and CMake 3.20 or later.

```sh
./run_tests.sh
```

This command configures, builds, and runs Formula, Compound, and Adduct tests together. It works from any working directory and exits with a nonzero status if configuration, compilation, or a test fails. CTest prints a summary and failure details.

To also run Python, Java, and C# integration tests:

```sh
./run_tests.sh --bindings
```

Set `PYTHON` to the Python executable with pybind11, RDKit, and bidict installed. Optional environment variables are `BUILD_TYPE` (default `Debug`), `JOBS` (default `2`), and `BUILD_DIR`. Core and binding tests use separate build directories under `build/`.

See [RELEASE.md](RELEASE.md) for release builds, package verification, and distribution steps.

```cpp
#include <mmkit/Formula.hpp>

auto neutral = mmkit::Formula::parse("C6H12O6");
auto ion = neutral.add(mmkit::Formula::parse("H")).with_charge(1);
auto mass = ion.exact_mass();
```

Other C++ classes can use `mmkit::Formula` directly as a value type. `mmkit::Compound` provides canonical SMILES, molecular formulas, charge, exact mass, atom maps, and RDKit binary conversion. `mmkit::Adduct` parses ion notation, calculates formula and mass shifts, applies adducts to formulas or masses, and splits terms or matches references.

```cpp
#include <mmkit/Compound.hpp>

auto compound = mmkit::Compound::from_smiles("CCO");
auto formula = compound.formula(); // C2H6O
```

Compound preserves existing atom maps and fills missing maps by default. Pass `true` as the second constructor or `from_smiles` argument to regenerate maps. `copy()` also regenerates maps. Binary serialization is an RDKit interoperability format, not a portable storage format.

```sh
cmake --install build/tests-off --prefix /your/install/path
```

After installation, use `find_package(mmkit CONFIG REQUIRED)` and `target_link_libraries(your_target PRIVATE mmkit::mmkit)` in consuming projects.

The default build and core tests use only C++. Language bindings are optional.

## Adduct

```cpp
#include <mmkit/mmkit.hpp>

auto compound = mmkit::Compound::from_smiles("CCO");
auto adduct = mmkit::Adduct::parse("[M+H]+");
auto ion = adduct.apply_to_formula(compound.formula()); // C2H7O+
auto mz = adduct.apply_to_mz(compound.exact_mass());
```

Both `[M+2H]+2` and `[M+2H]2+` charge suffixes are accepted. Mass shifts use neutral atomic masses without electron corrections. Reference matching requires equal formula counts, matching ion type and molecule count, and the same charge polarity. Matched terms are consumed in reference order; residuals have one molecule and zero charge.

## Language bindings

The main implementation lives in `mmkit/*.cpp`. Python, Java, and C# bindings live under `bindings/` and use the same C++ classes.

| Language | Usage guide | Integration tests |
| --- | --- | --- |
| Python | [Python usage](docs/python.md) | [tests/python/](tests/python/) |
| Java | [Java usage](docs/java.md) | [tests/java/](tests/java/) |
| C# | [C# usage](docs/csharp.md) | [tests/csharp/](tests/csharp/) |

Each guide covers prerequisites, building, using Formula, Compound, and Adduct, and running integration tests.
