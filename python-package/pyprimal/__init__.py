# coding: utf-8

import os

def test():
    """Show welcome information."""
    current_file = os.path.dirname(__file__)
    print("Python package psm has been successfully imported!")
    print("Version: "+open(os.path.join(current_file, './VERSION')).read().strip())

__all__ = ["model"]

from .model import QuantileRegression, SparseSVM, Dantzig, CompressedSensing
