# hIPPYlib - Inverse Problem PYthon library

[![Build Status](https://travis-ci.org/hippylib/hippylib.svg?branch=master)](https://travis-ci.org/hippylib/hippylib)
[![Doc Status](https://readthedocs.org/projects/hippylib/badge/?version=latest&style=flat)](https://hippylib.readthedocs.io/en/latest/)
[![status](http://joss.theoj.org/papers/053e0d08a5e9755e7b78898cff6f6208/status.svg)](http://joss.theoj.org/papers/053e0d08a5e9755e7b78898cff6f6208) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596931.svg)](https://doi.org/10.5281/zenodo.596931)

hIPPYlib implements state-of-the-art *scalable* *adjoint-based* algorithms for PDE-based *deterministic and Bayesian inverse problems*. It builds on <a href="http://www.fenicsproject.org" target="_blank">FEniCS</a> for the discretization of the PDE and on <a href="http://www.mcs.anl.gov/petsc/" target="_blank">PETSc</a> for scalable and efficient linear algebra operations and solvers.

## Features

- Friendly, compact, near-mathematical FEniCS notation to express the PDE and likelihood in weak form
- Automatic generation of efficient code for the discretization of weak forms using FEniCS
- Symbolic differentiation of weak forms to generate derivatives and adjoint information
- Globalized Inexact Newton-CG method to solve the inverse problem
- Low rank representation of the posterior covariace using randomized algorithms

See also our [tutorial](tutorial.md) and list of related [publications](research.md). For additional resources and tutorials please see the teaching material for the *2018 Gene Golub SIAM Summer School* on *Inverse Problems: Systematic Integration of Data with Models under Uncertainty* available [here](https://g2s3-2018.github.io/labs).

The complete API reference is available [here](http://hippylib.readthedocs.io/en/latest/index.html).

## Latest Release

- [Development version](https://github.com/hippylib/hippylib) (supports FEniCS-2019.1)
- Download [hippylib-2.2.1.zip](https://zenodo.org/record/2614052/files/hippylib/hippylib-2.2.1.zip?download=1)
- [Previous releases](download.md)

## News

- [Post-Doctoral Scholar position](postdoc_position.md) opening in Prof. Petra group at UC Merced.

## Contact

Developed by the [hIPPYlib team](about.md) at <a href="http://ices.utexas.edu" target="_blank">UT Austin</a> and <a href="http://naturalsciences.ucmerced.edu/" target="_blank">UC Merced</a>.

Please cite as 

```sh
@article{VillaPetraGhattas2016,
title = "{hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inverse Problems}",
author = {Villa, U. and Petra, N. and Ghattas, O.},
year = {2016},
url = {http://hippylib.github.io},
doi = {10.5281/zenodo.596931}
}

@article{VillaPetraGhattas2018,
title = "{hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inverse Problems}",
author = {Villa, U. and Petra, N. and Ghattas, O.},
journal = {Journal of Open Source Software},
volume = {3},
number = {30},
page = {940},
doi  = {10.21105/joss.00940},
year = {2018}
}
```


