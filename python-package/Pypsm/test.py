# coding: utf-8
import numpy as np
from model import Dantzig, SparseSVM

def main():
    n = 100
    d = 20
    X = np.random.rand(n, d).astype('double')
    y = np.random.rand(n).astype('double')
    solver = Dantzig(X, y)
    solver.train()
    print(solver.result['T'])
    print(solver.result['lambda_list'])
    print(solver.result['target_list'])
    print(solver.result['time'])

if __name__ == "__main__":
    main()