# Python usage

Source-build commands run from the repository root. Release installation does not require a repository checkout. The binding calls the shared C++ implementation of Formula, Compound, and Adduct.

## Install a specific release

Install the versioned source distribution from GitHub Releases:

```sh
python -m pip install "https://github.com/MotohiroOGAWA/mmkit/releases/download/v0.2.0/mmkit-0.2.0.tar.gz"
```

This builds the extension locally, so the build requirements below apply. To avoid compiling, select a compatible wheel on the release page and use `python -m pip install "<exact-wheel-asset-URL>"`. Wheels must match your Python ABI and platform and still require compatible RDKit runtime libraries. The `rdkit` and `bidict` dependencies are resolved by pip from your configured package index.

Change both `v0.2.0` and the filename version to select another release. This command works after the corresponding assets have been published.

## Source build requirements

Python 3.10+, a C++17 compiler, CMake 3.20+, and RDKit C++ development headers and libraries are required. Package installation builds the extension using pybind11 and scikit-build-core. RDKit and bidict are Python runtime dependencies.

## Build and install

```sh
python -m pip install .
```

Use compatible Python and C++ RDKit versions when passing RDKit Mol objects. Reinstall an existing editable installation after moving the package sources.

## Use the classes

```python
from mmkit import Formula, Compound, Adduct

water = Formula.parse("H2O")
print(water.value)       # H2O
print(water.exact_mass)  # approximately 18.010565

compound = Compound.from_smiles("CCO")
print(compound.formula.value)  # C2H6O

adduct = Adduct.parse("[M+H]+")
ion = adduct.apply_to_formula(compound.formula)
print(ion.value)  # C2H7O+
print(adduct.apply_to_mz(compound.exact_mass))  # approximately 47.049690
```

Python exposes values such as `formula`, `charge`, and `exact_mass` as properties. Formula arithmetic supports `+`, `-`, and integer multiplication. Invalid formula, compound, or adduct input raises `ValueError`.

## Verify

After installation:

```sh
python tests/python/smoke.py
```

To build and test every binding together, install pybind11 in the selected interpreter and provide the Java and C# prerequisites listed in their guides:

```sh
PYTHON=python ./run_tests.sh --bindings
```

The integration test is [tests/python/smoke.py](../tests/python/smoke.py). Release packaging is described in [RELEASE.md](../RELEASE.md).
