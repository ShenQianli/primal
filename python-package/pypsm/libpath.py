# coding: utf-8
"""Find the path to psm dynamic library files."""

import os
import platform
import sys

def find_lib_path():
    """Find the path to psm dynamic library files.

    :return: List of all found library path to psm
    :rtype: list(string)
    """
    curr_path = os.path.dirname(os.path.abspath(os.path.expanduser(__file__)))
    dll_path = [os.path.join(curr_path, './lib/')]

    if sys.platform == 'win32':
        dll_path = [os.path.join(p, 'psm.dll') for p in dll_path] \
                    +[os.path.join(p, 'libpsm.so') for p in dll_path]
    elif sys.platform.startswith('linux'):
        dll_path = [os.path.join(p, 'libpsm.so') for p in dll_path]
    elif sys.platform == 'darwin':
        dll_path = [os.path.join(p, 'libpsm.so') for p in dll_path] \
                    +[os.path.join(p, 'libpsm.dylib') for p in dll_path]

    lib_path = [p for p in dll_path if os.path.exists(p) and os.path.isfile(p)]
    
    if not lib_path:
        print('Library file does not exist. Need to be updated!')
        return lib_path

    return lib_path
