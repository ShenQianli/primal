# Installation

## From source (recommended)

```bash
git clone https://github.com/Gatech-Flash/primal.git
cd primal

# Build the C++ library
make clean && make dylib

# Install the Python package
cd python-package
pip install -e ".[viz,test]"
```

## Optional extras

```bash
pip install -e ".[viz]"   # matplotlib for plotting
pip install -e ".[test]"  # pytest for testing
```

## Verify install

```bash
python -c "import pyprimal; pyprimal.test()"
```

Expected output:

```
pyprimal 2.0.0 loaded successfully.
C library: libpsm loaded.
```

## PEP 668 environments (externally-managed)

If `pip install` shows `externally-managed-environment`, use a virtual env:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[viz,test]"
```

## Requirements

- Python >= 3.8
- NumPy >= 1.23
- Compiled `libpsm` shared library (Linux `.so` or macOS `.dylib`)
- Optional: matplotlib >= 3.5 for plotting
