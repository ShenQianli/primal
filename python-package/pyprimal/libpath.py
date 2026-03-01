"""Find the path to the libpsm shared library."""

import os
import sys
from typing import List


def find_lib_path() -> List[str]:
    """Find the path to the libpsm shared library.

    Returns
    -------
    list of str
        Paths to found library files.

    Raises
    ------
    RuntimeError
        If no library file is found.
    """
    curr_path = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
    lib_dir = os.path.join(curr_path, "lib")

    if sys.platform == "darwin":
        candidates = ["libpsm.dylib", "libpsm.so"]
    elif sys.platform == "win32":
        candidates = ["psm.dll", "libpsm.so"]
    else:
        candidates = ["libpsm.so"]

    lib_path = [
        os.path.join(lib_dir, name)
        for name in candidates
        if os.path.isfile(os.path.join(lib_dir, name))
    ]

    if not lib_path:
        raise RuntimeError(
            "Cannot find libpsm shared library. "
            "Please build from source: see README for instructions."
        )

    return lib_path
