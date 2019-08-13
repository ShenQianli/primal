# coding: utf-8
"""
Main Interface of the package
"""

import time
import numpy as np
import ctypes
from numpy.ctypeslib import ndpointer

class psm(object):
    """
    Parametric Simplex Method class for solving Parametric Linear Programming

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
        self.X = X
        self.y = y
        self.n = X.shape[0]
        self.d = X.shape[1]
        self.family =family

        self.result = {'state': 'not trained'}

    def __str__(self):
        return "Parametric Simplex Method class of family %s" % (self.family)

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
            _function(np.array([self.n], dtype='int32'), np.array([self.d], dtype='int32'), self.X, self.y,
                np.array([max_it], dtype='int32'), np.array([lambda_threshold], dtype='double'), 
                self.result['T'], self.result['lambda_list'],
                self.result['theta_list'], self.result['target_list'])
            time_end = time.time()
            self.result['T'] = self.result['T'][0]
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
        self.result['T'] = self.result['T'][0]
        for key in self.result.keys():
            if key.endswith('_list'):
                self.result[key] = self.result[key][:T] 

    def coef(self):
        '''
        Extract model coefficients
        '''
        if self.result['state'] == 'not trained':
            print("Warning: The model has not been trained yet!")
        return self.result

    def plot(self):
        '''
        Visualize
        '''
        if self.result['state'] == 'not trained':
            print("The model has not been trained yet!")
        else:
            import matplotlib.pyplot as plt
            plt.plot(self.result['lambda_list'], self.result['target_list'])
            plt.ylabel('Target')
            plt.xlabel('Regularization Parameter')
            plt.suptitle('Regularization Path')
            plt.show()

    def predict(self):
        pass
        



