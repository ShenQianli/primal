# coding: utf-8
"""
Main Interface of the package
"""

import time
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer

from .libpath import find_lib_path

def _load_lib():
    """Load library."""
    lib_path = find_lib_path()
    lib = ctypes.cdll.LoadLibrary(lib_path[0])
    return lib

class PSMsolver(object):
    """
    Parametric Simplex Method solver for solving Parametric Linear Programming

    """
    def __init__(self, A, b, b_bar, c, c_bar, B_init=None):
        self.A = np.array(A, dtype='double')
        self.M = self.A.shape[0]
        self.N = self.A.shape[1]
        self.b = np.array(b, dtype='double')
        self.b_bar = np.array(b_bar, dtype='double')
        self.c = np.array(c, dtype='double')
        self.c_bar = np.array(c_bar, dtype='double')
        self.B_init = np.array(B_init, dtype='int32') if B_init != None else np.array([-1], dtype='int32')
        self._PSM_LIB = _load_lib()
        self._function = self._PSM_LIB.ParametricSimplexMethod_api
        self.result = {'state': 'not trained'}

    def __str__(self):
        return "Parametric Simplex Method solver"

    def train(self, max_it=100, lambda_threshold=1e-3):
        """
        :param max_it(int): max iteration num
        :param lambda_threshold(double): iteration ends once lambda reach the threshold

        """
        self.result.update({
            'T' : np.zeros(1, dtype='int32'),
            'lambda_list' : np.zeros(max_it, dtype='double'),
            'theta_list': np.zeros((max_it, self.N), dtype='double'),
            'target_list': np.zeros(max_it, dtype='double'),
            'time': 0,
        })

        FDoubleArray = ndpointer(ctypes.c_double, flags='F_CONTIGUOUS')
        CDoubleArray = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')
        CIntArray = ndpointer(ctypes.c_int, flags='C_CONTIGUOUS')

        time_start = time.time()
        self._function.argtypes = [
            CIntArray, CIntArray, CDoubleArray, CDoubleArray, CDoubleArray, 
            CDoubleArray, CDoubleArray, CIntArray, CDoubleArray, CIntArray, 
            CDoubleArray, CDoubleArray, CDoubleArray, CIntArray,
        ]
        self._function(
            np.array([self.M], dtype='int32'), 
            np.array([self.N], dtype='int32'), 
            self.A,
            self.b,
            self.b_bar,
            self.c,
            self.c_bar,
            np.array([max_it], dtype='int32'), 
            np.array([lambda_threshold], dtype='double'), 
            self.result['T'], 
            self.result['lambda_list'],
            self.result['theta_list'], 
            self.result['target_list'],
            self.B_init,
            )
        time_end = time.time()
        self.result['T'] = self.result['T'][0]
        for key in self.result.keys():
            if key.endswith('_list'):
                self.result[key] = self.result[key][:self.result['T']] 
        self.result['time'] = time_end - time_start

        self.result['state'] = 'trained'

    def coef(self):
        '''
        Extract model coefficients
        '''
        if self.result['state'] == 'not trained':
            print("Warning: The model has not been trained yet!")
        return self.result

    def plot(self, mode='trajectory'):
        '''
        Visualize
        '''
        if self.result['state'] == 'not trained':
            print("The model has not been trained yet!")
        else:
            import matplotlib.pyplot as plt
            if mode == 'trajectory':
                plt.plot(self.result['lambda_list'], self.result['target_list'])
                plt.ylabel('Target')
                plt.xlabel('Regularization Parameter')
                plt.suptitle('Trajectory')
                plt.show()
            elif mode == 'regpath':
                plt.plot(self.result['lambda_list'], self.result['theta_list'])
                plt.ylabel('Coefficient')
                plt.xlabel('Regularization Parameter')
                plt.suptitle('Regularization Path')
                plt.show()
            else:
                raise NotImplementedError

    def predict(self):
        pass

class PSMmodel(object):
    """
    Parametric Simplex Method model basic class

    :param X(array): input data with shape [n, d]
    :param y(array): input label with shape [d, ]
    :param family(str): Options for model. 
        'QuantileRegression': Quantile regression
        'SparseSVM': Sparse Support Vector Machine
        'Dantzig': Dantzig selector
        'CompressedSensing': Compressed Sensing
    """
    def __init__(self, X, y, family):
        if family not in ['QuantileRegression', 'SparseSVM', 'Dantzig', 'CompressedSensing']:
            raise NotImplementedError
        # self.X = np.asfortranarray(X, dtype='double')
        # self.y = np.ascontiguousarray(y, dtype='double')
        self.X = np.array(X, dtype='double')
        self.y = np.array(y, dtype='double')
        self.n = self.X.shape[0]
        self.d = self.X.shape[1]
        self.family =family

        self.result = {'state': 'not trained'}

    def __str__(self):
        return "Parametric Simplex Method model of family %s" % (self.family)

    def _decor_cinterface(self, _function, max_it, lambda_threshold):
        '''
            Default C interface
        '''
        FDoubleArray = ndpointer(ctypes.c_double, flags='F_CONTIGUOUS')
        CDoubleArray = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')
        CIntArray = ndpointer(ctypes.c_int, flags='C_CONTIGUOUS')
        _function.argtypes = [
            CIntArray, CIntArray, CDoubleArray, CDoubleArray, CIntArray, 
            CDoubleArray, CIntArray, CDoubleArray, CDoubleArray, CDoubleArray
        ]
        def wrapper():
            time_start = time.time()
            _function(
                np.array([self.n], dtype='int32'), 
                np.array([self.d], dtype='int32'), 
                self.X, 
                self.y,
                np.array([max_it], dtype='int32'), 
                np.array([lambda_threshold], dtype='double'), 
                self.result['T'], 
                self.result['lambda_list'],
                self.result['theta_list'], 
                self.result['target_list'],
                )
            time_end = time.time()
            self.result['T'] = self.result['T'][0]
            for key in self.result.keys():
                if key.endswith('_list'):
                    self.result[key] = self.result[key][:self.result['T']] 
            self.result['time'] = time_end - time_start

        return wrapper

    def train(self, max_it=100, lambda_threshold=1e-3):
        """
        :param max_it(int): max iteration num
        :param lambda_threshold(double): iteration ends once lambda reach the threshold

        """
        self.result.update({
            'T' : np.zeros(1, dtype='int32'),
            'lambda_list' : np.zeros(max_it, dtype='double'),
            'theta_list': np.zeros((max_it, self.d), dtype='double'),
            'target_list': np.zeros(max_it, dtype='double'),
            'time': 0,
        })

        self.trainer()  #self.trainer constructed in child class
        self.result['state'] = 'trained'

    def coef(self):
        '''
        Extract model coefficients
        '''
        if self.result['state'] == 'not trained':
            print("Warning: The model has not been trained yet!")
        return self.result

    def plot(self, mode='trajectory'):
        '''
        Visualize
        '''
        if self.result['state'] == 'not trained':
            print("The model has not been trained yet!")
        else:
            import matplotlib.pyplot as plt
            if mode == 'trajectory':
                plt.plot(self.result['lambda_list'], self.result['target_list'])
                plt.ylabel('Target')
                plt.xlabel('Regularization Parameter')
                plt.suptitle('Trajectory')
                plt.show()
            elif mode == 'regpath':
                plt.plot(self.result['lambda_list'], self.result['theta_list'])
                plt.ylabel('Coefficient')
                plt.xlabel('Regularization Parameter')
                plt.suptitle('Regularization Path')
                plt.show()
            else:
                raise NotImplementedError

    def predict(self):
        pass
        



