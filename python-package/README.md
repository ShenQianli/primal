# PYPRIMAL


**PYPRIMAL**: **PY**thon package **P**a**R**ametric s**I**mplex **M**ethod for sp**A**rse **L**earning


Requirements
------------

- Linux or MacOS


Installation
------------

Install from source file (Github) with Makefile:

- Clone ``primal.git`` via ``git clone --recurse-submodules https://github.com/ShenQianli/primal.git``
- Make sure you have [setuptools](https://pypi.python.org/pypi/setuptools)  
- Run ``make Pyinstall`` command.


Install from source file (Github) with CMAKE:

- Clone ``primal.git`` via ``git clone --recurse-submodules https://github.com/ShenQianli/primal.git``
- Make sure you have [setuptools](https://pypi.python.org/pypi/setuptools) 
- Build the source file first via the ``cmake`` with ``CMakeLists.txt`` in the root directory. (You will see a ``.so`` or ``.dylib`` file under ``(root)/lib/`` )
- Run ``cd python-package; sudo python setup.py install`` command.


Install from PyPI:

- ``pip install pyprimal``
- **Note**: Owing to the setting on different OS, our distribution might not be working in your environment (especially in **Windows**). Thus please build from source.

You can test if the package has been successfully installed by:

```python
import pyprimal
pyprimal.test()

```

Usage
-----

```python
from pyprimal import SparseSVM
x = [[1,2,3], [4,5,6], [7,8,9]]
y = [-1, 1, 1]
solver = SparseSVM(x, y)
result = solver.coef()
solver.plot()
solver.plot('regpath')
```

See [tutorial](https://github.com/ShenQianli/primal/blob/master/tutorials/tutorial.ipynb)


Copy Right
----------

Author: Qianli Shen, Zichong Li  
Maintainer: Qianli Shen <shenqianli@pku.edu.cn>
