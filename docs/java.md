# Java usage

Source-build commands run from the repository root. Release installation does not require a repository checkout. The binding calls the shared C++ implementation of Formula, Compound, and Adduct. Commands and native filenames below target Linux.

## Install a specific release

```sh
mkdir -p vendor/mmkit-1.0.0
curl -fL "https://github.com/MotohiroOGAWA/mmkit/releases/download/v1.0.0/mmkit-1.0.0-java-linux-x86_64.tar.gz" -o vendor/mmkit-java.tar.gz
tar -xzf vendor/mmkit-java.tar.gz -C vendor/mmkit-1.0.0
```

This downloads the prebuilt JAR and JNI library; consumers do not need a compiler, CMake, or SWIG. A Java runtime and compatible RDKit shared libraries are required. Compile and run the example below against the downloaded files:

```sh
javac -cp vendor/mmkit-1.0.0/mmkit-1.0.0.jar Main.java
java -Djava.library.path="$PWD/vendor/mmkit-1.0.0" -cp "vendor/mmkit-1.0.0/mmkit-1.0.0.jar:." Main
```

Change the tag, filenames, and paths to select another version. Choose an archive matching your platform. The URL is available after the corresponding release is published.

## Source build requirements

JDK 11+, SWIG 4+, a C++17 compiler, CMake 3.20+, and RDKit C++ development headers and libraries are required. The matching RDKit shared libraries must be available at runtime.

## Build

```sh
cmake -S . -B build/java -DCMAKE_BUILD_TYPE=Release -DMMKIT_BUILD_JAVA=ON
cmake --build build/java --parallel 2
```

The build produces `build/java/mmkit-1.0.0.jar` and `build/java/libmmkit_java.so`. Substitute the selected release version when using another version.

## Use the classes

Save the following as `Main.java`:

```java
import org.mmkit.*;

public class Main {
    public static void main(String[] args) {
        System.loadLibrary("mmkit_java");
        try (Formula water = Formula.parse("H2O");
             Compound compound = Compound.from_smiles("CCO");
             Formula formula = compound.formula();
             Adduct adduct = Adduct.parse("[M+H]+");
             Formula ion = adduct.apply_to_formula(formula)) {
            System.out.println(water.to_string()); // H2O
            System.out.println(formula.to_string()); // C2H6O
            System.out.println(ion.to_string()); // C2H7O+
            System.out.println(adduct.apply_to_mz(compound.exact_mass()));
        }
    }
}
```

Compile and run:

```sh
javac -cp build/java/mmkit-1.0.0.jar Main.java
java -Djava.library.path="$PWD/build/java" -cp "build/java/mmkit-1.0.0.jar:." Main
```

Java uses methods such as `formula()` and `exact_mass()`. Native objects implement `AutoCloseable`; use try-with-resources to release them. Invalid formula, compound, or adduct input raises `IllegalArgumentException`.

`java.library.path` locates the JNI library. If RDKit libraries are outside the system loader paths, also set `LD_LIBRARY_PATH` to their directory.

## Verify

```sh
ctest --test-dir build/java -R '^java$' --output-on-failure --no-tests=error
```

The integration test is [tests/java/Smoke.java](../tests/java/Smoke.java). Use `./run_tests.sh --bindings` to build and test all bindings together. Release packaging is described in [RELEASE.md](../RELEASE.md).
