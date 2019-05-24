# Documentation

## Installation

`hIPPYlib` depends on [FEniCS](http://fenicsproject.org/) versions 1.6, 2016.1, 2016.2, 2017.1, 2017.2.

We recommend using `FEniCS 2017.2` with `hIPPYlib`.

`FEniCS` needs to be built with the following dependecies:

 - `numpy`, `scipy`, `matplotlib`, `mpi4py`
 - `PETSc` and `petsc4py` (version 3.7.0 or above)
 - `SLEPc` and `slepc4py` (version 3.7.0 or above)
 - PETSc dependencies: `parmetis`, `scotch`, `suitesparse`, `superlu_dist`, `ml`, `hypre`
 - (optional): `mshr`, `jupyter`

> For detailed installation instructions see [here](http://hippylib.readthedocs.io/en/latest/installation.html).

## FEniCS 2019.1.0

`hIPPYlib` now offers experimental support for `FEniCS 2019.1.0`, although it requires a customized [FEniCS docker image](https://hub.docker.com/r/hippylib/fenics). See branch [2019.1-dev2](https://github.com/hippylib/hippylib/tree/2019.1-dev2) for details.

## hiPPYlib Docker container

A Docker image `hIPPYlib`, `FEniCS` and their dependencies preinstalled is available [here](https://hub.docker.com/r/mparno/muq-hippylib/). The username is `user1` and password `Breckenridge1_g2s3`.

## Documentation

The complete API reference of `hIPPYlib` is available at [readthedocs](http://hippylib.readthedocs.io/en/latest/modules.html).

