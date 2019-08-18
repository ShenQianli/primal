# coding: utf-8
import time
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer

from .core import psm
from .libpath import find_lib_path

__all__ = ["QuantileRegression", "SparseSVM", "Dantzig", "CompressedSensing"]

def _load_lib():
	"""Load library."""
	lib_path = find_lib_path()
	lib = ctypes.cdll.LoadLibrary(lib_path[0])
	return lib

# _PSM_LIB = ctypes.cdll.LoadLibrary("/Users/shenqianli/Desktop/PSM/lib/libpsm.dylib")
_PSM_LIB = _load_lib()

class QuantileRegression(psm):
    """docstring for QuantileRegression"""
    def __init__(self, X, y, tau):
        super(QuantileRegression, self).__init__(X, y, 'QuantileRegression')
        self.tau = tau

    def _decor_cinterface(self, _function, max_it, lambda_threshold):
        FDoubleArray = ndpointer(ctypes.c_double, flags='F_CONTIGUOUS')
        CDoubleArray = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')
        CIntArray = ndpointer(ctypes.c_int, flags='C_CONTIGUOUS')
        _function.argtypes = [
            CIntArray, CIntArray, CDoubleArray, CDoubleArray, CDoubleArray, CIntArray, 
            CDoubleArray, CIntArray, CDoubleArray, CDoubleArray, CDoubleArray
        ]
        def wrapper():
            time_start = time.time()
            _function(np.array([self.n], dtype='int32'), np.array([self.d], dtype='int32'), 
                self.X, self.y, np.array([self.tau], dtype='double'), np.array([max_it], dtype='int32'),
                np.array([lambda_threshold], dtype='double'), self.result['T'], self.result['lambda_list'],
                self.result['theta_list'], self.result['target_list'])
            time_end = time.time()
            self.result['T'] = self.result['T'][0]
            for key in self.result.keys():
                if key.endswith('_list'):
                    self.result[key] = self.result[key][:self.result['T']]
            self.result['time'] = time_end - time_start

        return wrapper

    def train(self, max_it=100, lambda_threshold=1e-3):
        self.trainer = self._decor_cinterface(_PSM_LIB.QuantileRegression_api, max_it, lambda_threshold)
        super(QuantileRegression, self).train(max_it, lambda_threshold)

class SparseSVM(psm):
    """docstring for SparseSVM"""
    def __init__(self, X, y):
        super(SparseSVM, self).__init__(X, y, 'SparseSVM')
    
    def _decor_cinterface(self, _function, max_it, lambda_threshold):
        FDoubleArray = ndpointer(ctypes.c_double, flags='F_CONTIGUOUS')
        CDoubleArray = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')
        CIntArray = ndpointer(ctypes.c_int, flags='C_CONTIGUOUS')
        _function.argtypes = [
            CIntArray, CIntArray, CDoubleArray, CDoubleArray, CIntArray, 
            CDoubleArray, CIntArray, CDoubleArray, CDoubleArray, CDoubleArray, CDoubleArray
        ]
        def wrapper():
            time_start = time.time()
            _function(np.array([self.n], dtype='int32'), np.array([self.d], dtype='int32'), self.X, self.y,
                np.array([max_it], dtype='int32'), np.array([lambda_threshold], dtype='double'), 
                self.result['T'], self.result['lambda_list'],
                self.result['theta_list'], self.result['target_list'], self.result['theta0_list'])
            time_end = time.time()
            self.result['T'] = self.result['T'][0]
            for key in self.result.keys():
                if key.endswith('_list'):
                    self.result[key] = self.result[key][:self.result['T']]
            self.result['time'] = time_end - time_start

        return wrapper

    def train(self, max_it=100, lambda_threshold=1e-3):
        self.trainer = self._decor_cinterface(_PSM_LIB.SparseSVM_api, max_it, lambda_threshold)
        self.result['theta0_list'] = np.zeros(max_it, dtype='double')
        super(SparseSVM, self).train(max_it, lambda_threshold)


class Dantzig(psm):
    """docstring for Dantzig"""
    def __init__(self, X, y):
        super(Dantzig, self).__init__(X, y, 'Dantzig')
    def train(self, max_it=100, lambda_threshold=1e-3):
        self.trainer = self._decor_cinterface(_PSM_LIB.Dantzig_api, max_it, lambda_threshold)
        super(Dantzig, self).train(max_it, lambda_threshold)

class CompressedSensing(psm):
    """docstring for CompressedSensing"""
    def __init__(self, X, y):
        super(CompressedSensing, self).__init__(X, y, 'CompressedSensing')
    def train(self, max_it=100, lambda_threshold=1e-3):
        self.trainer = self._decor_cinterface(_PSM_LIB.CompressedSensing_api, max_it, lambda_threshold)
        super(CompressedSensing, self).train(max_it, lambda_threshold)
