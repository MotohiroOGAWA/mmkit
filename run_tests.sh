#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bindings=OFF
case "${1:-}" in
    "") ;;
    --bindings) bindings=ON ;;
    -h|--help)
        echo 'Usage: ./run_tests.sh [--bindings]'
        echo 'Runs Formula, Compound, and Adduct tests; --bindings also tests Python, Java, and C#.'
        echo 'Environment: BUILD_DIR, BUILD_TYPE (Debug), JOBS (2), PYTHON (python3).'
        exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
esac
if (( $# > 1 )); then
    echo 'Too many arguments. Use --help for usage.' >&2
    exit 2
fi
build_dir="${BUILD_DIR:-${project_dir}/build/tests-${bindings,,}}"
build_type="${BUILD_TYPE:-Debug}"
cmake_args=(
    -S "$project_dir" -B "$build_dir"
    "-DCMAKE_BUILD_TYPE=$build_type" -DBUILD_TESTING=ON
    "-DMMKIT_BUILD_PYTHON=$bindings"
    "-DMMKIT_BUILD_JAVA=$bindings"
    "-DMMKIT_BUILD_CSHARP=$bindings"
)
if [[ "$bindings" == ON ]]; then
    python_executable="$(command -v "${PYTHON:-python3}")"
    pybind11_dir="$("$python_executable" -m pybind11 --cmakedir)"
    cmake_args+=("-DPython_EXECUTABLE=$python_executable" "-Dpybind11_DIR=$pybind11_dir")
fi
cmake "${cmake_args[@]}"
cmake --build "$build_dir" --config "$build_type" --parallel "${JOBS:-2}"
ctest --test-dir "$build_dir" --build-config "$build_type" --output-on-failure --no-tests=error
printf '\nAll requested tests passed. Results: %s/Testing/Temporary/LastTest.log\n' "$build_dir"
