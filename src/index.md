# hIPPYlib - Inverse Problem PYthon library

hIPPYlib implements state-of-the-art *scalable* *adjoint-based* algorithms for PDE-based *deterministic and Bayesian inverse problems*. It builds on <a href="http://www.fenicsproject.org" target="_blank">FEniCS</a> for the discretization of the PDE and on <a href="http://www.mcs.anl.gov/petsc/" target="_blank">PETSc</a> for scalable and efficient linear algebra operations and solvers.

## Features

- Friendly, compact, near-mathematical FEniCS notation to express the PDE and likelihood in weak form
- Automatic generation of efficient code for the discretization of weak forms using FEniCS
- Symbolic differentiation of weak forms to generate derivatives and adjoint information
- Globalized Inexact Newton-CG method to solve the inverse problem
- Low rank representation of the posterior covariace using randomized algorithms.

See also our [tutorial](tutorial.md) and list of related [publications](outreach.md).

## Latest Release

- ** Initial release coming soon **
- [Development version](https://github.com/hippylib/hippylib)
- [Download hippylib-1.0.0.tgz](https://github.com/hippylib/hippylib-releases)

## Contact

Developed by the [hIPPYlib team](about.md) at <a href="http://ices.utexas.edu" target="_blank">UT Austin</a> and <a href="http://naturalsciences.ucmerced.edu/" target="_blank">UC Merced</a>.

To ask question and find answers see <a href="https://groups.google.com/forum/#!forum/hippylib-support" target="_blank">here</a>.

Please cite with 
```sh
@article{VillaPetraGhattas2016,
title = {"hIPPYlib: AN EXTENSIBLE SOFTWARE FRAMEWORK FOR LARGE-SCALE DETERMINISTIC AND LINEARIZED BAYESIAN INVERSION"},
author = {Villa U., and Petra, N., and Ghattas, O.},
year = {2016},
url = {http://hippylib.github.io}
}
```

