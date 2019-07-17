# Documentation

## Installation of stable releases

The latest `hIPPYlib` release depends on [FEniCS](http://fenicsproject.org/) versions 2016.1, 2016.2, 2017.1, 2017.2.

We recommend using `FEniCS 2017.2` with the latest `hIPPYlib` release.

`FEniCS` needs to be built with the following dependecies:

 - `numpy`, `scipy`, `matplotlib`, `mpi4py`
 - `PETSc` and `petsc4py` (version 3.7.0 or above)
 - `SLEPc` and `slepc4py` (version 3.7.0 or above)
 - PETSc dependencies: `parmetis`, `scotch`, `suitesparse`, `superlu_dist`, `ml`, `hypre`
 - (optional): `mshr`, `jupyter`

> For detailed installation instructions of the latest stable release see [here](https://hippylib.readthedocs.io/en/2.2.0/installation.html).

## Installation of developement version (with FEniCS 2019.1.0 support)

The [development version](https://github.com/hippylib/hippylib) of `hIPPYlib` offers experimental support for `FEniCS 2019.1.0`. However, a customized [FEniCS docker image](https://hub.docker.com/r/hippylib/fenics) or customized [conda packages](https://anaconda.org/uvilla/fenics) is required at this time.

> For detailed installation instructions of the development version see [here](https://hippylib.readthedocs.io/en/latest/installation.html).

## hiPPYlib Docker container

A Docker image `hIPPYlib`, `FEniCS` and their dependencies preinstalled is available [here](https://hub.docker.com/r/mparno/muq-hippylib/). The username is `user1` and password `Breckenridge1_g2s3`.

## Documentation

The complete API reference of `hIPPYlib` is available at [readthedocs](http://hippylib.readthedocs.io/en/latest/modules.html).

