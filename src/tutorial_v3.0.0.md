# Tutorial - Development version

These tutorials are the best place to learn about the basic features and the algorithms in `hIPPYlib`.

> For the complete API reference click [here](http://hippylib.readthedocs.io/en/latest/index.html).

1. [FEniCS101](tutorials_v3.0.0/1_FEniCS101.md) notebook illustrates the use of FEniCS for the solution of a linear boundary value problem.
2. [Poisson Deterministic](tutorials_v3.0.0/2_PoissonDeterministic.md) notebook illustrates how to compute gradient/Hessian information and solve a non-linear parameter inversion for the Poisson equation in a deterministic setting.
3. [Subsurface Bayesian](tutorials_v3.0.0/3_SubsurfaceBayesian.md) notebook illustrates how to solve a non-linear parameter inversion for the Poisson equation in a Bayesian setting.
4. [Advection-Diffusion Bayesian](tutorials_v3.0.0/4_AdvectionDiffusionBayesian.md) notebook illustrates how to solve a time-dependent linear inverse problem in a Bayesian setting.
5. [Hessian Spectrum](tutorials_v3.0.0/5_HessianSpectrum.md) notebook illustrates the spectral property of the Hessian operator for a linear source inversion problem.

## Interactive tutorials

The interactive ipython notebooks are located in the `tutorial` folder of the `hIPPYlib` repository.

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

## Additional resources

For additional resources and tutorials please see the teaching material for the
*2018 Gene Golub SIAM Summer School* on *Inverse Problems: Systematic Integration of Data with Models under Uncertainty* available [here](https://g2s3-2018.github.io/labs).

> These tutorials require `hIPPYlib` versions 2.1.1 - 2.3.0.
