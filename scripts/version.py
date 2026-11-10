#!/usr/bin/env python3
"""Check, set, or increment the shared MAJOR.MINOR.PATCH release version."""
import argparse
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
choice = parser.add_mutually_exclusive_group()
choice.add_argument('version', nargs='?', help='New MAJOR.MINOR.PATCH version')
choice.add_argument('--bump', choices=['major', 'minor', 'patch'])
parser.add_argument('--dry-run', action='store_true', help='Print the version without changing files')
args = parser.parse_args()
pattern = r'(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)'
current = (root / 'VERSION').read_text().strip()
if not re.fullmatch(pattern, current):
    parser.error('VERSION must contain MAJOR.MINOR.PATCH')
version = args.version or current
if args.bump:
    major, minor, patch = map(int, current.split('.'))
    version = {'major': f'{major + 1}.0.0',
               'minor': f'{major}.{minor + 1}.0',
               'patch': f'{major}.{minor}.{patch + 1}'}[args.bump]
if not re.fullmatch(pattern, version):
    parser.error('Expected MAJOR.MINOR.PATCH, for example 0.2.0')
patterns = {
    'CMakeLists.txt': r'(?<=project\(mmkit VERSION )[0-9]+\.[0-9]+\.[0-9]+',
    'pyproject.toml': r'(?m)(?<=^version = ")[0-9]+\.[0-9]+\.[0-9]+',
    'bindings/python/mmkit/__init__.py': r'(?m)(?<=^__version__ = ")[0-9]+\.[0-9]+\.[0-9]+',
}
updates = {}
for name, version_pattern in patterns.items():
    path = root / name
    text = path.read_text()
    match = re.search(version_pattern, text)
    if not match:
        raise SystemExit(f'Missing version in {name}')
    if not args.version and match.group() != current:
        raise SystemExit(f'Version mismatch in {name}; run python scripts/version.py {current}')
    updates[path] = re.sub(version_pattern, version, text)
if (args.version or args.bump) and not args.dry_run:
    # Keep release URLs and installation examples aligned with the package version.
    for path in [root / 'README.md', root / 'RELEASE.md', *(root / 'docs').glob('*.md')]:
        updates[path] = re.sub(r'(?<![0-9.])' + re.escape(current) + r'(?![0-9.])', version, path.read_text())
    for path, text in updates.items():
        path.write_text(text)
    (root / 'VERSION').write_text(version + '\n')
print(version)
