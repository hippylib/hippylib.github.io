# hIPPYlib - Inverse Problem PYthon library

hIPPYlib implements state-of-the-art *scalable* *adjoint-based* algorithms for PDE-based *deterministic and Bayesian inverse problems*. It builds on [FEniCS](www.fenicsproject.org) for the discretization of the PDE and on [PETSc](http://www.mcs.anl.gov/petsc/) for scalable and efficient linear algebra operations and solvers.

## Features:

- Friendly, compact, near-mathematical FEniCS notation to express the PDE and likelihood in weak form
- Automatic generation of efficient code for the discretization of weak forms using FEniCS
- Symbolic differentiation of weak forms to generate derivatives and adjoint information
- Globalized Inexact Newton-CG method to solve the inverse problem
- Low rank representation of the posterior covariace using randomized algorithms.

See also our [tutorial](tutorial.md) and list of related [publications](publications.md).

## Latest Release:

- ** Initial release coming soon **
- [New features]()
- [Development version](https://github.com/hippylib/hippylib)
- [Download hippylib-1.0.0.tgz](https://github.com/hippylib/hippylib-releases)
- For older releases see the [download](download.md) section

## Contact

Developed by the [hIPPYlib team](about.md) at [UT Austin](http://ices.utexas.edu) and [UC Merced](http://naturalsciences.ucmerced.edu/).

To ask question and find answers see [here](https://groups.google.com/forum/#!forum/hippylib-support).

Please cite with 
```
@article{VillaPetraGhattas2016,
title = {"hIPPYlib: AN EXTENSIBLE SOFTWARE FRAMEWORK FOR LARGE-SCALE BAYESIAN INVERSION"},
author = {Villa U., and Petra, N., and Ghattas, O.},
year = {2016},
url = {http://hippylib.github.io}
}
```


