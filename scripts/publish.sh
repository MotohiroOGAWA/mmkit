#!/usr/bin/env bash
set -euo pipefail
export PATH="${HOME}/.local/bin:$PATH"
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$project_dir"
if ! command -v gh >/dev/null; then
    echo 'GitHub CLI (gh) is required. Install it and run: gh auth login' >&2
    exit 1
fi
gh auth status >/dev/null
if [[ -n "$(git status --porcelain)" ]]; then
    echo 'Commit all intended source and version changes before publishing. Run git status to review them.' >&2
    exit 1
fi
version="$("${PYTHON:-python3}" scripts/version.py)"
platform="$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)"
assets="$project_dir/dist/v$version/$platform"
if [[ ! -f "$assets/SHA256SUMS" ]]; then
    echo "Release assets are missing: $assets" >&2
    exit 1
fi
(cd "$assets" && sha256sum -c SHA256SUMS)
# Check that the packages belong to the version being tagged.
[[ -f "$assets/mmkit-$version.tar.gz" && -f "$assets/Mmkit.$version.$platform.nupkg" ]]
repo=MotohiroOGAWA/mmkit
tag="v$version"
revision="$(git rev-parse HEAD)"
if git rev-parse "$tag" >/dev/null 2>&1; then
    if [[ "$(git rev-list -n 1 "$tag")" != "$revision" ]]; then
        echo "$tag points to a different commit. Use a new version." >&2
        exit 1
    fi
else
    git tag -a "$tag" -m "Release $version"
fi
git push origin HEAD "refs/tags/$tag"
notes="$assets/release-notes.md"
printf '# mmkit %s\n\nPlatform: %s. Requires compatible RDKit shared libraries.\n\nInstallation: https://github.com/%s/blob/%s/README.md\n' "$version" "$platform" "$repo" "$tag" > "$notes"
if gh release view "$tag" --repo "$repo" >/dev/null 2>&1; then
    echo "Release $tag already exists. Assets will not be overwritten." >&2
    exit 1
fi
gh release create "$tag" "$assets"/*.whl "$assets"/*.tar.gz "$assets"/*.nupkg \
    "$assets/LICENSE" "$assets/RELEASE.md" "$assets/SHA256SUMS" \
    --repo "$repo" --verify-tag --draft --title "mmkit $version" --notes-file "$notes"
gh release edit "$tag" --repo "$repo" --draft=false
printf '\nPublished: https://github.com/%s/releases/tag/%s\n' "$repo" "$tag"
