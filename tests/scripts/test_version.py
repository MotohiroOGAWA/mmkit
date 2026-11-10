"""Regression tests for release metadata and installation URL synchronization."""
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


class VersionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = pathlib.Path(self.temp.name)
        for name in ['VERSION', 'CMakeLists.txt', 'pyproject.toml',
                     'bindings/python/mmkit/__init__.py', 'scripts/version.py',
                     'README.md', 'RELEASE.md', 'docs/python.md', 'docs/java.md', 'docs/csharp.md']:
            target = self.root / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / name, target)

    def run_version(self, *args):
        return subprocess.check_output(
            [sys.executable, str(self.root / 'scripts/version.py'), *args], text=True
        ).strip()

    def test_stale_assets_extensions_and_unrelated_versions(self):
        readme = self.root / 'README.md'
        readme.write_text(readme.read_text() + '''
Python 3.11.2 is available.
https://github.com/MotohiroOGAWA/mmkit/releases/download/v0.2.0/mmkit-0.2.0.tar.gz
Mmkit.0.2.0.linux-x86_64.nupkg
mmkit-0.2.0.jar
''')
        self.assertEqual(self.run_version('2.4.6'), '2.4.6')
        text = readme.read_text()
        self.assertIn('download/v2.4.6/mmkit-2.4.6.tar.gz', text)
        self.assertIn('Mmkit.2.4.6.linux-x86_64.nupkg', text)
        self.assertIn('mmkit-2.4.6.jar', text)
        self.assertIn('Python 3.11.2', text)
        for name in ['README.md', 'RELEASE.md', 'docs/python.md', 'docs/java.md', 'docs/csharp.md']:
            text = (self.root / name).read_text()
            self.assertNotRegex(text, r'(?i)mmkit[-.]0\.2\.0')
        release = (self.root / 'RELEASE.md').read_text()
        for level, version in [('major', '3.0.0'), ('minor', '2.5.0'), ('patch', '2.4.7')]:
            self.assertIn(f'| {level} | `{version}` |', release)
        self.assertEqual(self.run_version(), '2.4.6')

    def test_bump_and_dry_run(self):
        self.run_version('2.4.6')
        before = {p: p.read_bytes() for p in self.root.rglob('*') if p.is_file()}
        self.assertEqual(self.run_version('--bump', 'patch', '--dry-run'), '2.4.7')
        self.assertEqual(before, {p: p.read_bytes() for p in before})
        self.assertEqual(self.run_version('--bump', 'patch'), '2.4.7')
        self.assertEqual(self.run_version('--bump', 'minor'), '2.5.0')
        self.assertEqual(self.run_version('--bump', 'major'), '3.0.0')
        self.assertIn('download/v3.0.0/mmkit-3.0.0.tar.gz', (self.root / 'README.md').read_text())


if __name__ == '__main__':
    unittest.main()
