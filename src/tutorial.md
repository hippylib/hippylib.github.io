# Tutorial

These tutorials are the best place to learn about the basic features and the algorithms in `hIPPYlib`.

1. [FEniCS101](tutorials/1_FEniCS101/index.md) notebook illustrates the use of FEniCS for the solution of a linear boundary value problem.
2. [Poisson Deterministic](tutorials/2_PoissonDeterministic/index.md) notebook illustrates how to compute gradient/Hessian information and solve a non-linear parameter inversion for the Poisson equation in a deterministic setting.
3. [Subsurface Bayesian](tutorials/3_SubsurfaceBayesian/index.md) notebook illustrates how to solve a non-linear parameter inversion for the Poisson equation in a Bayesian setting.
4. [Advection-Diffusion Bayesian](tutorials/4_AdvectionDiffusionBayesian/index.md) notebook illustrates how to solve a time-dependent linear inverse problem in a Bayesian setting.
5. [Hessian Spectrum](tutorials/5_HessianSpectrum/index.md) notebook illustrates the spectral property of the Hessian operator for a linear source inversion problem.


The interactive ipython notebooks are located in the `tutorial` folder of the `hIPPYlib` release.

To run the notebooks follow these instructions.

1. Open a FEniCS terminal and type

```ssh
$ cd tutorial
$ jupyter notebook
```

2. A new tab will open in your web-brower showing the notebooks.
3. Click on the notebook you would like to use.
4. To run all the code in the notebook simply click on Cell --> Run All.

For more information on installing ipython and using notebooks see <a href="https://jupyter.readthedocs.io/en/latest/content-quickstart.html" target="_blank">here</a>.
