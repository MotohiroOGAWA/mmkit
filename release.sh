#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"
python_executable="${PYTHON:-python3}"
selection=""
dry_run=OFF
for argument in "$@"; do
    case "$argument" in
        major|minor|patch|--current)
            if [[ -n "$selection" ]]; then
                echo 'Choose only one of major, minor, patch, or --current.' >&2
                exit 2
            fi
            selection="$argument" ;;
        --dry-run) dry_run=ON ;;
        -h|--help)
            echo 'Usage: ./release.sh [major|minor|patch|--current] [--dry-run]'
            echo 'Without a version choice, prompts interactively. --current keeps the existing version.'
            echo '--dry-run previews the version and output path without modifying or building anything.'
            exit 0 ;;
        *) echo "Unknown argument: $argument" >&2; exit 2 ;;
    esac
done
current="$("$python_executable" scripts/version.py)"
if [[ -z "$selection" ]]; then
    printf 'Current version: %s\n' "$current"
    printf 'Select version increment:\n  1) major\n  2) minor\n  3) patch\n'
    if ! read -r -p 'Choice [1/2/3 or major/minor/patch]: ' selection; then
        echo 'A version choice is required. Use ./release.sh major, minor, or patch.' >&2
        exit 2
    fi
    case "$selection" in
        1) selection=major ;;
        2) selection=minor ;;
        3) selection=patch ;;
        major|minor|patch) ;;
        *) echo 'Invalid version choice.' >&2; exit 2 ;;
    esac
fi
version_args=()
if [[ "$selection" != --current ]]; then
    version_args=(--bump "$selection")
fi
version="$("$python_executable" scripts/version.py "${version_args[@]}" --dry-run)"
platform="$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)"
if [[ "$platform" != linux-* ]]; then
    echo 'The release packaging script currently supports Linux only.' >&2
    exit 1
fi
build_dir="$project_dir/build/release-$version-$platform"
output_dir="$project_dir/dist/v$version/$platform"
printf 'Release version: %s -> %s\nOutput: %s\n' "$current" "$version" "$output_dir"
if [[ "$dry_run" == ON ]]; then
    exit 0
fi
if [[ -e "$output_dir" ]]; then
    echo "Output already exists: $output_dir. Use a new version or move the previous output before rebuilding." >&2
    exit 1
fi
"$python_executable" scripts/version.py "${version_args[@]}" > /dev/null
BUILD_DIR="$build_dir" BUILD_TYPE=Release PYTHON="$python_executable" ./run_tests.sh --bindings
mkdir -p "$output_dir"
CMAKE_BUILD_PARALLEL_LEVEL="${JOBS:-2}" "$python_executable" -m build --outdir "$output_dir"
stage_dir="$build_dir/stage"
cmake --install "$build_dir" --config Release --prefix "$stage_dir" --component cpp
tar -czf "$output_dir/mmkit-$version-cpp-$platform.tar.gz" -C "$stage_dir" .
java_dir="$build_dir/java-release"
mkdir -p "$java_dir"
cp "$build_dir/mmkit-$version.jar" "$build_dir/libmmkit_java.so" LICENSE "$java_dir/"
tar -czf "$output_dir/mmkit-$version-java-$platform.tar.gz" -C "$java_dir" .
cp "$build_dir/packages/Mmkit.$version.nupkg" "$output_dir/Mmkit.$version.$platform.nupkg"
cp LICENSE RELEASE.md "$output_dir/"
(cd "$output_dir" && sha256sum ./*.whl ./*.tar.gz ./*.nupkg ./LICENSE ./RELEASE.md > SHA256SUMS)
printf '\nRelease files prepared in %s\nNo files have been uploaded.\n' "$output_dir"
