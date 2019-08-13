# pylint: disable=invalid-name, exec-used
"""Setup psm package."""
from __future__ import absolute_import
import sys
import os
import shutil
from setuptools import setup, find_packages
# import subprocess
sys.path.insert(0, '.')

CURRENT_DIR = os.path.dirname(__file__)

#try to copy the complied lib files
libdir_candidate = [os.path.join(CURRENT_DIR, '../lib/')]

if sys.platform == 'win32':
    libcand_path = [os.path.join(p, 'psm.dll') for p in libdir_candidate]
    libcand_path = libcand_path + [os.path.join(p, 'libpsm.so') for p in libdir_candidate]
elif sys.platform.startswith('linux'):
    libcand_path = [os.path.join(p, 'libpsm.so') for p in libdir_candidate]
elif sys.platform == 'darwin':
    libcand_path = [os.path.join(p, 'libpsm.so') for p in libdir_candidate]
    libcand_path = libcand_path + [os.path.join(p, 'libpsm.dylib') for p in libdir_candidate]

lib_path = [p for p in libcand_path if os.path.exists(p) and os.path.isfile(p)]
if not os.path.exists('./pypsm/lib/'):
	os.mkdir('./pypsm/lib/')
for lib_file in lib_path:
    shutil.copy(lib_file,os.path.join(CURRENT_DIR, './pypsm/lib/'))


# We can not import `psm.libpath` in setup.py directly, since it will automatically import other package
# and case conflict to `install_requires`
libpath_py = os.path.join(CURRENT_DIR, 'pypsm/libpath.py')
libpath = {'__file__': libpath_py}
exec(compile(open(libpath_py, "rb").read(), libpath_py, 'exec'), libpath, libpath)
LIB_PATH = [os.path.relpath(libfile, CURRENT_DIR) for libfile in libpath['find_lib_path']()]
if not LIB_PATH:
    raise RuntimeError("libpsm does not exists")
else:
    print("libpsm already exists: %s" % LIB_PATH)


VERSION_PATH = os.path.join(CURRENT_DIR, 'pypsm/VERSION')

setup(name='pypsm',
	  author = 'Qianli Shen',
	  author_email = 'shenqianli@pku.edu.cn',
      version=open(VERSION_PATH).read().strip(),
      description="psm Python Package",
      install_requires=[
          'numpy',
      ],
      packages=find_packages(),
	  include_package_data=True,
      )
