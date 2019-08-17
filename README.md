<h1 align="center">PRIMAL</h1>
<h4 align="center">High Performance R and Python Library for Sparse Learning</h4>

___PRIMAL___ (PaRametric sImplex Method for spArse Learning) implements a unified framework of parametric simplex method for a variety of sparse learning problems (e.g., Dantzig selector (for linear regression), sparse quantile regression, sparse support vector machines, and compressive sensing) combined with efficient hyper-parameter selection strategies. The core algorithm is implemented in C++ with Eigen3 support for portable high performance linear algebra. Runtime profiling is documented in the [__Performance__](#performance) section.

## Table of contents

- [Table of contents](#table-of-contents)
- [Introduction](#introduction)
- [Directory structure](#directory-structure)
- [Installation](#installation)
    - [Installing R package](#installing-r-package)
    - [Installing Python package](#installing-python-package)
- [Performance](#performance)
    - [R package](#r-package)
    - [Python package](#python-package)
- [References](#references)

## Introduction

Linear Programming (LP) based sparse learning methods, such as the Dantzig selector (for linear regression) [1], sparse quantile regression [2], sparse support vector machines [3], have been widely used in machine learning for high dimensional data analysis [2, 4, 5]. Despite of their popularity, their software implementations are quite limited. This library -- PaRametric sImplex Method for spArse Learning (PRIMAL) is proposed for the aforementioned LP-based sparse learning methods. It has the following two key features: 1) It provides a highly efficient optimization engine based on the parametric simplex method [6], which can efficiently solve large scale sparse learning problems; 2) Besides the estimation procedures, it provides additional functional modules such as data-dependent model selection and model visualization.

We also provide tutorials on the theoretical background and the code. Please see ``tutorial`` folder for tutorials.


## Directory structure
The directory is organized as follows:
* [__src__](src): C++ implementation of the PSM algorithm.
	* [__api.cpp__](api.cpp): C API as an interface for R and Python package.
	* [__PSM.cpp__](PSM.cpp): Core implemetation of the PSM solver. 
* [__include__](include)
	* [__PSM__](PSM): declarations of the C++ implementation
	* [__Eigen__](eigen3): Eigen3 header files for high performance linear algebra.
* [__R-package__](R-package): R wrapper for the source code.
* [__python-package__](python-package): Python wrapper for the source code.
* [__tutorials__](tutorials): tutorials for using the code in R and Python.
* [__profiling__](profiling): profiling the performance from R package.
* [__Makefile__](Makefile):Makefile local configurations.
* [__CmakeLists.txt__](CmakeLists.txt):Makefile local configurations.



## Installation
Third-party dependencies. The installation only depends on Eigen3's header files which are included in this github repo as a submodule. We use Eigen3 as an independent fully portable module so any existing Eigen3 installation will not have conflict with PSM installation.
### Installing R package
There are two ways to install the picasso R package.
- Installing from CRAN (recommended). The R package is hosted on CRAN. The easiest way to install R package is by running the following command in R
```R
install.packages("PSM")
```

- Installing from source code.
```bash
$ git clone --recurse-submodules https://github.com/ShenQianli/PSM.git
$ cd PSM; make Rinstall
```

### Installing Python package
There are two ways to install the PSM (pypsm) python package.
- Installing from PyPi (recommended). ``pip install pypsm --user``.
- Installing from source code.
 ```bash
 $git clone --recurse-submodules https://github.com/ShenQianli/PSM.git
 $cd PSM; make Pyinstall
 ```

You can test if the package has been successfully installed by ``python -c "import pypsm; psm.test()" ``

Details for installing python package can also be found in [document](https://hmjianggatech.github.io/picasso/) or [github](https://github.com/jasonge27/picasso/tree/master/python-package)


## Performance
```bash
$cd profiling
$Rscript benchmark.R
$python benchmark.py
```

### R package

We compare the timing performance of our package with R package "fastclime" and "flare". We fix the sample size n to be 200 and vary the data dimension d from 200 to 1000. Each entries of X is independent Gaussian and Gaussianized such that the column has uniform norm. We randomly select 2% entries from vector θ to be nonzero. Algorithm will stop when λ is less than $$2*\sigma*\sqrt{log(d)/n}$$. In "flare", we estimate θ when λ is equal to the smallest value in the Dantzig selector model. This means "flare" just give out one result with one specific λ, while our package give out a full solution path.
- Dantzig selector. PRIMAL achieves similar optimization performance to fastclime and flare. But PRIMAL is 2 times faster than fastclime and 3-30 times faster than flare.
- Compressed sensing. PRIMAL is 2 times faster than fastclime and achieves similar optimization.
![Performance_R](https://raw.githubusercontent.com/ShenQianli/primal/master/profiling/images/performance_R.png?token=AMWY2UDGNJRVOKHAM3NUUWS5K55JQ)

### Python package

Some description.

![A figure]()


## References

[1] Dantzig G. Linear programming and extensions, 2016.

[2] Belloni A, Chernozhukov V. ℓ1-penalized quantile regression in high-dimensional sparse models, 2011.

[2] Candes E, Tao T. The Dantzig selector: Statistical estimation when p is much larger than n, 2007.

[3] Wang L. The L1 penalized LAD estimator for high dimensional linear regression, 2013.

[4] Hudson N J, Reverter A, Dalrymple B P. A differential wiring analysis of expression data correctly identifies the gene containing the causal mutation, 2009.

[5] Bandyopadhyay S, Mehta M, Kuo D, et al. Rewiring of genetic networks in response to DNA damage, 2010.

[6] Pang H, Liu H, Vanderbei R, Zhao T. Parametric simplex method for sparse learning, 2017.




 
