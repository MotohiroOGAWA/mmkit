# Versioned releases

All languages are distributed as assets of the same GitHub Release:

`https://github.com/MotohiroOGAWA/mmkit/releases/tag/v<VERSION>`

Package versions use `MAJOR.MINOR.PATCH`; Git tags use `vMAJOR.MINOR.PATCH`. Consumers select a fixed tag and matching asset filename. Publishing to PyPI, Maven Central, or NuGet.org is not required. No release is published by the local build script.

## 1. Choose the release version

Run the script from the repository root and select `major`, `minor`, or `patch` at the prompt:

```sh
./release.sh
```

The script increments the selected component, updates package versions and versioned documentation examples, runs all tests, and prepares every language's release assets. From `0.2.0`, the choices are:

| Choice | New version |
| --- | --- |
| major | `1.0.0` |
| minor | `0.3.0` |
| patch | `0.2.1` |

You can also select the increment without a prompt:

```sh
./release.sh patch
./release.sh minor --dry-run
```

`--dry-run` previews the new version and output path without changing files or running builds. `./release.sh --current` builds the existing version without incrementing it, and refuses to overwrite an existing output directory.

`VERSION` is checked against CMake, Python metadata, and the Python runtime version. Java and C# versions come from CMake. To set an explicit version manually, use `python scripts/version.py 0.2.0`; to check consistency, use `python scripts/version.py`.

Increment PATCH for compatible fixes, MINOR for compatible features, and MAJOR for breaking changes. Never replace published assets; publish a new version instead. If a build fails after the version update, correct the failure and retry with `--current` after moving any incomplete output directory, rather than incrementing again. Review and commit the generated version changes before tagging and publishing in step 5.

## 2. Prepare the environment

The release script currently targets Linux. Build separately for each supported OS/architecture and Python ABI; do not claim support for untested platforms.

Required tools:

- C++17 compiler, CMake 3.20+, RDKit C++ headers and shared libraries.
- Python 3.10+, development headers, pybind11, scikit-build-core, build, RDKit, and bidict.
- JDK 11+, SWIG 4+, .NET SDK 8+, tar, and sha256sum.
- GitHub CLI for the optional publication commands below.

The Python and C++ RDKit versions must be compatible for molecule interchange. Native artifacts currently depend on external RDKit shared libraries; they are not self-contained runtime bundles. Record the exact build platform and RDKit version in the release notes.

## 3. Build and test every language

```sh
./release.sh patch
```

Set `PYTHON=/path/to/python` and `JOBS=2` as needed. `release.sh` checks version consistency, runs all six suites in Release mode, builds a Python source distribution and wheel, stages C++, and collects Java and C# packages. It stops on errors and refuses to overwrite an existing output directory.

For version `0.2.0` on Linux x86_64, output is `dist/v0.2.0/linux-x86_64/`:

| Language | Release asset |
| --- | --- |
| Python | `mmkit-0.2.0.tar.gz` and `mmkit-0.2.0-<python>-<abi>-<platform>.whl` |
| Java | `mmkit-0.2.0-java-linux-x86_64.tar.gz` containing JAR and JNI library |
| C# | `Mmkit.0.2.0.linux-x86_64.nupkg` containing managed and native libraries |
| C++ | `mmkit-0.2.0-cpp-linux-x86_64.tar.gz` containing headers, static library, and CMake package |
| Shared | `LICENSE`, `RELEASE.md`, and `SHA256SUMS` |

For multiple platforms, collect uniquely named native assets and wheels into the same release. Upload the source distribution once. Keep platform checksum manifests distinct when combining output directories.

## 4. Verify distribution files

The integration tests exercise bindings built from source. Also test the actual release packages in clean environments on the target platform before publication:

- Python: install the wheel by local path, run `tests/python/smoke.py`, then build and install from the `.tar.gz` in a separate clean environment. Ensure pip does not reuse an earlier package.
- Java: extract the release archive, compile `tests/java/Smoke.java` against its JAR, and run with `java.library.path` pointing to the extracted JNI library.
- C#: download or copy the package into a local package feed, reference it with `dotnet add package Mmkit --version 0.2.0 --source <feed-directory>`, and run the example from docs/csharp.md. This verifies package consumption rather than a project reference.
- C++: extract the archive and build a consumer with `find_package(mmkit CONFIG REQUIRED)` and `mmkit::mmkit`, setting `CMAKE_PREFIX_PATH` to the extraction directory.

Run `ldd` on the wheel extension and JNI/PInvoke libraries to identify runtime dependencies. Include compatible RDKit installation instructions and platform requirements in the release notes. Python source-archive consumers also need the native build tools. Wheels are platform-specific and are not automatically repaired or bundled for broad Linux compatibility.

Check the release files before uploading:

```sh
cd dist/v0.2.0/linux-x86_64
sha256sum -c SHA256SUMS
```

## 5. Publish to GitHub Releases

After committing and reviewing the exact release revision, create its version tag and draft release. Run these commands only when ready to publish; the build script does not execute them:

```sh
git tag -a v0.2.0 -m "Release 0.2.0"
git push origin v0.2.0
gh release create v0.2.0 dist/v0.2.0/linux-x86_64/* \
  --repo MotohiroOGAWA/mmkit --verify-tag --draft \
  --title "mmkit 0.2.0" --notes-file release-notes.md
```

Write `release-notes.md` before running the command. Include changes, supported platform and Python ABI, RDKit runtime requirements, and links to all language guides. Review the draft assets, then publish:

```sh
gh release edit v0.2.0 --repo MotohiroOGAWA/mmkit --draft=false
```

## 6. Verify the published URLs

The immutable versioned asset URL pattern is:

`https://github.com/MotohiroOGAWA/mmkit/releases/download/v<VERSION>/<ASSET>`

Use the commands in [README.md](README.md), [Python](docs/python.md), [Java](docs/java.md), and [C#](docs/csharp.md) against the public release from a clean environment. For Python, this includes:

```sh
python -m pip install --no-cache-dir "https://github.com/MotohiroOGAWA/mmkit/releases/download/v0.2.0/mmkit-0.2.0.tar.gz"
```

A compatible wheel can be installed by its exact asset URL instead. Package dependencies may still be resolved from the consumer's configured package index. Confirm the installed `mmkit.__version__` is the requested version and run the integration tests. A successful local build alone does not establish that public download URLs are available.
