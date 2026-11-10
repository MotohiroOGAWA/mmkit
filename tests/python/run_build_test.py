"""Run integration tests against the selected build, including with editable installs."""
import importlib.abc
import importlib.machinery
from pathlib import Path
import runpy
import sys


class BuildPackageFinder(importlib.abc.MetaPathFinder):
    def __init__(self, package_root):
        self.package_root = package_root

    def find_spec(self, fullname, path=None, target=None):
        if fullname == "mmkit" or fullname.startswith("mmkit."):
            spec = importlib.machinery.PathFinder.find_spec(
                fullname, [str(self.package_root)] if path is None else path
            )
            if spec is None:
                raise ModuleNotFoundError(f"{fullname} is missing from the selected build")
            return spec
        return None


if __name__ == "__main__":
    package_root = Path(sys.argv[1]).resolve()
    if not (package_root / "mmkit" / "__init__.py").is_file():
        raise RuntimeError(f"Built package not found: {package_root}")
    # Keep dependency imports intact while giving this build priority over installed mmkit.
    sys.meta_path.insert(0, BuildPackageFinder(package_root))
    runpy.run_path(str(Path(__file__).with_name("smoke.py")), run_name="__main__")
    import mmkit
    import mmkit._native

    for module in (mmkit, mmkit._native):
        origin = Path(module.__file__).resolve()
        if not origin.is_relative_to(package_root):
            raise RuntimeError(f"Loaded module outside the selected build: {origin}")
    print(f"Verified build package: {mmkit.__file__}")
