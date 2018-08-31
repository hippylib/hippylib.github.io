# hIPPYlib - Inverse Problem PYthon library

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

- [Development version](https://github.com/hippylib/hippylib)
- Download [hippylib-2.1.0.tar.gz](https://goo.gl/N7g6wU)
- [Previous releases](download.md)

## Contact

Developed by the [hIPPYlib team](about.md) at <a href="http://ices.utexas.edu" target="_blank">UT Austin</a> and <a href="http://naturalsciences.ucmerced.edu/" target="_blank">UC Merced</a>.

Please cite as 

```sh
@article{VillaPetraGhattas2016,
title = "{hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inversion}",
author = {Villa, U. and Petra, N. and Ghattas, O.},
year = {2016},
url = {http://hippylib.github.io},
doi = {10.5281/zenodo.596931}
}

@article{VillaPetraGhattas2018,
title = "{hIPPYlib: an Extensible Software Framework for Large-scale Deterministic and Bayesian Inversion}",
author = {Villa, U. and Petra, N. and Ghattas, O.},
journal = {Journal of Open Source Software},
volume = {submitted},
year = {2018}
}
```


