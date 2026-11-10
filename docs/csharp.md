# C# usage

Source-build commands run from the repository root. Release installation does not require a repository checkout. The binding calls the shared C++ implementation of Formula, Compound, and Adduct. Commands and native filenames below target Linux.

## Install a specific release

From a .NET application directory:

```sh
mkdir -p vendor/packages
curl -fL "https://github.com/MotohiroOGAWA/mmkit/releases/download/v1.0.0/Mmkit.0.2.0.linux-x86_64.nupkg" -o vendor/packages/Mmkit.0.2.0.nupkg
dotnet add package Mmkit --version 1.0.0 --source "$PWD/vendor/packages"
dotnet run
```

The downloaded NuGet package includes the managed classes and native library for Linux x86_64. Consumers need .NET 8 and compatible RDKit shared libraries; a compiler, CMake, and SWIG are only needed to build from source. Replace the application's `Program.cs` with the example below before running it. If needed, add the RDKit library directory to `LD_LIBRARY_PATH`.

Change both the release tag and package version to select another release, and choose an asset matching your platform. The URL is available after the corresponding release is published.

## Source build requirements

.NET SDK 8+, SWIG 4+, a C++17 compiler, CMake 3.20+, and RDKit C++ development headers and libraries are required. The matching RDKit shared libraries must be available at runtime.

## Build

```sh
cmake -S . -B build/csharp -DCMAKE_BUILD_TYPE=Release -DMMKIT_BUILD_CSHARP=ON
cmake --build build/csharp --parallel 2
```

The build produces:

- `build/csharp/csharp/bin/Release/net8.0/Mmkit.dll`
- `build/csharp/libmmkit_csharp.so`
- `build/csharp/packages/Mmkit.0.2.0.nupkg`

The NuGet package includes the native library for the build platform. Substitute the selected release version when using another version.

## Use the classes

Create an application and reference the local package:

```sh
dotnet new console --framework net8.0 --output build/csharp-app
dotnet add build/csharp-app package Mmkit --version 1.0.0 --source "$PWD/build/csharp/packages"
```

Replace `build/csharp-app/Program.cs` with:

```csharp
using System;
using Mmkit;

using var water = Formula.parse("H2O");
using var compound = Compound.from_smiles("CCO");
using var formula = compound.formula();
using var adduct = Adduct.parse("[M+H]+");
using var ion = adduct.apply_to_formula(formula);

Console.WriteLine(water.to_string()); // H2O
Console.WriteLine(formula.to_string()); // C2H6O
Console.WriteLine(ion.to_string()); // C2H7O+
Console.WriteLine(adduct.apply_to_mz(compound.exact_mass()));
```

Run the application:

```sh
LD_LIBRARY_PATH="$PWD/build/csharp${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" dotnet run --project build/csharp-app
```

C# uses methods such as `formula()` and `exact_mass()`. Native objects implement `IDisposable`; use `using` declarations to release them. Invalid formula, compound, or adduct input raises `ArgumentException`.

If RDKit libraries are outside the system loader paths, add their directory to `LD_LIBRARY_PATH` as well.

## Verify

```sh
ctest --test-dir build/csharp -R '^csharp$' --output-on-failure --no-tests=error
```

The integration test is [tests/csharp/Program.cs](../tests/csharp/Program.cs), with its project in [Smoke.csproj](../tests/csharp/Smoke.csproj). Use `./run_tests.sh --bindings` to build and test all bindings together. Release packaging is described in [RELEASE.md](../RELEASE.md).
