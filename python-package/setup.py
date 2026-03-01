"""Build helper: compile libpsm shared library before packaging."""

import platform
import shutil
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py

try:
    from wheel.bdist_wheel import bdist_wheel
except ImportError:
    bdist_wheel = None

# Paths relative to this file (python-package/)
_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parent  # primal-master/
_LIB_DIR = _HERE / "pyprimal" / "lib"


def _lib_name():
    """Return the platform-specific shared library name."""
    system = platform.system()
    if system == "Darwin":
        return "libpsm.dylib"
    elif system == "Windows":
        return "libpsm.dll"
    return "libpsm.so"


class BuildPyWithLib(build_py):
    """Custom build_py that compiles the C++ shared library first."""

    def run(self):
        self._ensure_lib()
        super().run()

    def _ensure_lib(self):
        lib_name = _lib_name()
        target = _LIB_DIR / lib_name

        # Already have a pre-built library — nothing to do
        if target.is_file():
            return

        makefile = _REPO_ROOT / "Makefile"
        if not makefile.is_file():
            raise RuntimeError(
                f"Cannot find pre-built {lib_name} in pyprimal/lib/ "
                f"and no Makefile found at {_REPO_ROOT} for compilation. "
                "Please build from source — see README for instructions."
            )

        # Compile the shared library
        subprocess.check_call(["make", "dylib"], cwd=str(_REPO_ROOT))

        # Copy to pyprimal/lib/
        built = _REPO_ROOT / "lib" / lib_name
        if not built.is_file():
            raise RuntimeError(
                f"make dylib succeeded but {built} was not created."
            )
        _LIB_DIR.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(built), str(target))


# Tag the wheel as platform-specific (it contains a compiled .so/.dylib)
cmdclass = {"build_py": BuildPyWithLib}

if bdist_wheel is not None:

    class PlatformWheel(bdist_wheel):
        """Mark the wheel as platform-specific but Python-version-agnostic.

        Since we use ctypes (not a compiled Python extension), the wheel
        works with any CPython/PyPy 3.x but is tied to the OS/arch.
        """

        def finalize_options(self):
            super().finalize_options()
            self.root_is_pure = False

        def get_tag(self):
            _, _, plat = super().get_tag()
            return "py3", "none", plat

    cmdclass["bdist_wheel"] = PlatformWheel

setup(cmdclass=cmdclass)
