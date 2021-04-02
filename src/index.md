# hIPPYlib - Inverse Problem PYthon library

[![Build Status](https://travis-ci.org/hippylib/hippylib.svg?branch=master)](https://travis-ci.org/hippylib/hippylib)
[![Doc Status](https://readthedocs.org/projects/hippylib/badge/?version=latest&style=flat)](https://hippylib.readthedocs.io/en/latest/)
[![status](http://joss.theoj.org/papers/053e0d08a5e9755e7b78898cff6f6208/status.svg)](http://joss.theoj.org/papers/053e0d08a5e9755e7b78898cff6f6208) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.596931.svg)](https://doi.org/10.5281/zenodo.596931)

hIPPYlib implements state-of-the-art *scalable* *adjoint-based* algorithms for PDE-based *deterministic and Bayesian inverse problems*. It builds on <a href="http://www.fenicsproject.org" target="_blank">FEniCS</a> for the discretization of the PDE and on <a href="http://www.mcs.anl.gov/petsc/" target="_blank">PETSc</a> for scalable and efficient linear algebra operations and solvers.

## Features

- Friendly, compact, near-mathematical FEniCS notation to
express, differentiate, and discretize the PDE forward model and
likelihood function

- Large scale optimization algorithms, such as globalized inexact
Newton-CG method, to solve the inverse problem

- Randomized algorithms for trace estimation, eigenvalues and singular values decomposition.

- Scalable sampling of Gaussian random fields

- Linearized Bayesian inversion with low-rank based
representation of the posterior covariance

- Hessian informed MCMC algorithms to explore the posterior
  distribution

- Forward propagation of uncertainty capabilities using Monte
  Carlo and Taylor expansion control variates


See also our [tutorial](tutorial.md) and list of related [publications](research.md). For additional resources and tutorials please see the teaching material for the *2018 Gene Golub SIAM Summer School* on *Inverse Problems: Systematic Integration of Data with Models under Uncertainty* available [here](https://g2s3-2018.github.io/labs).

The complete API reference is available [here](http://hippylib.readthedocs.io/en/latest/index.html).

## Latest Release

- [Development version](https://github.com/hippylib/hippylib)
- Download [hippylib-3.0.0.zip](https://zenodo.org/record/3634136/files/hippylib/hippylib-3.0.0.zip?download=1)
- [Previous releases](download.md)

## Contact

Developed by the [hIPPYlib team](about.md) at <a href="http://ices.utexas.edu" target="_blank">UT Austin</a>, <a href="http://naturalsciences.ucmerced.edu/" target="_blank">UC Merced</a>, and <a href="https://ese.wustl.edu/Pages/default.aspx" target="_blank">WUSTL</a>.

## Slack channel

The hIPPYlib slack channel is a good resource to request and receive help with using hIPPYlib. Everyone is invited to read and take part in discussions. Discussions about development of new features in hIPPYlib also take place here. You can join our Slack community by filling in [this form](https://forms.gle/w8B7uKSXxdVCmfZ99). 

Please cite as 

```sh
@article{10.1145/3428447,
    author = {Villa, Umberto and Petra, Noemi and Ghattas, Omar},
    title = "{HIPPYlib: An Extensible Software Framework for Large-Scale Inverse Problems Governed by PDEs: Part I: Deterministic Inversion and Linearized Bayesian Inference}",
    year = {2021},
    issue_date = {March 2021},
    publisher = {Association for Computing Machinery},
    address = {New York, NY, USA},
    volume = {47},
    number = {2},
    issn = {0098-3500},
    url = {https://doi.org/10.1145/3428447},
    doi = {10.1145/3428447},
    journal = {ACM Trans. Math. Softw.},
    month = apr,
    articleno = {16},
    numpages = {34}
    }

      


@article{VillaPetraGhattas19,
author = {{Villa}, Umberto and {Petra}, Noemi and {Ghattas}, Omar},
title = "{hIPPYlib: An Extensible Software Framework for Large-Scale Inverse Problems Governed by PDEs; Part I: Deterministic Inversion and Linearized Bayesian Inference}",
journal = {arXiv e-prints},
year = {2019},
archivePrefix = {arXiv},
eprint = {1909.03948}
}

@article{VillaPetraGhattas18,
title = "{hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inverse Problems}",
author = {Villa, U. and Petra, N. and Ghattas, O.},
journal = {Journal of Open Source Software},
volume = {3},
number = {30},
page = {940},
doi  = {10.21105/joss.00940},
year = {2018}
}

@article{VillaPetraGhattas16,
title = "{hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inverse Problems}",
author = {Villa, U. and Petra, N. and Ghattas, O.},
year = {2016},
url = {http://hippylib.github.io},
doi = {10.5281/zenodo.596931}
}
```
