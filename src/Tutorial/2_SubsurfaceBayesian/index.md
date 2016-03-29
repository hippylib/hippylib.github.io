
$$\def\data{\bf d_\rm{obs}}
\def\vec{\bf}
\def\m{\bf m}
\def\map{\bf m_{\text{MAP}}}
\def\postcov{\bf \Gamma_{\text{post}}}
\def\prcov{\bf \Gamma_{\text{prior}}}
\def\matrix{\bf}
\def\Hmisfit{\bf H_{\text{misfit}}}
\def\HT{\tilde{\bf H}_{\text{misfit}}}
\def\diag{diag}
\def\Vr{\matrix V_r}
\def\Wr{\matrix W_r}
\def\Ir{\matrix I_r}
\def\Dr{\matrix D_r}
\def\H{\matrix H}
$$ 
# Example: Bayesian quantification of parameter uncertainty:
## Estimating the (Gaussian) posterior pdf of the coefficient parameter field in an elliptic PDE

In this example we tackle the problem of quantifying the
uncertainty in the solution of an inverse problem governed by an
elliptic PDE via the Bayesian inference framework. 
Hence, we state the inverse problem as a
problem of statistical inference over the space of uncertain
parameters, which are to be inferred from data and a physical
model.  The resulting solution to the statistical inverse problem
is a posterior distribution that assigns to any candidate set of
parameter fields our belief (expressed as a probability) that a
member of this candidate set is the ``true'' parameter field that
gave rise to the observed data.

For simplicity, in what follows we give finite-dimensional expressions (i.e., after
discretization of the parameter space) for the Bayesian
formulation of the inverse problem.

### Bayes' Theorem:

The posterior probability distribution combines the prior pdf
$\pi_{\text{prior}}(\m)$ over the parameter space, which encodes
any knowledge or assumptions about the parameter space that we may
wish to impose before the data are considered, with a likelihood pdf
$\pi_{\text{like}}(\vec{d}_{\text{obs}} \; | \; \m)$, which explicitly
represents the probability that a given set of parameters $\m$
might give rise to the observed data $\vec{d}_{\text{obs}} \in
\mathbb{R}^m$, namely:

$
\begin{align}
\pi_{\text{post}}(\m | \data) \propto
\pi_{\text{prior}}(\m) \pi_{\text{like}}(\data | \m).
\end{align}
$

Note that infinite-dimensional analog of Bayes' formula requires the use Radon-Nikodym derivatives instead of probability density functions.

### Gaussian prior and noise:

#### The prior:

We consider a Gaussian prior with mean $\vec m_{\text prior}$ and covariance $\bf \Gamma_{\text{prior}}$. The covariance is given by the discretization of the inverse of differential operator $\mathcal{A}^{-2} = (-\gamma \Delta + \delta I)^{-2}$, where $\gamma$, $\delta > 0$ control the correlation length and the variance of the prior operator. This choice of prior ensures that it is a trace-class operator, guaranteeing bounded pointwise variance and a well-posed infinite-dimensional Bayesian inverse problem

#### The likelihood:

$
\data =  \bf{f}(\m) + \bf{e }, \;\;\;  \bf{e} \sim \mathcal{N}(\bf{0}, \bf \Gamma_{\text{noise}} )
$

$
\pi_{\text like}(\data \; | \; \m)  = \exp \left( - \tfrac{1}{2} (\bf{f}(\m) - \data)^T \bf \Gamma_{\text{noise}}^{-1} (\bf{f}(\m) - \data)\right)
$

Here $\bf f$ is the parameter-to-observable map that takes a parameter vector $\m$ and maps
it to the space observation vector $\data$.

#### The posterior:

$
\pi_{\text{post}}(\m \; | \; \data)  \propto \exp \left( - \tfrac{1}{2} \parallel \bf{f}(\m) - \data \parallel^{2}_{\bf  \Gamma_{\text{noise}}^{-1}} \! - \tfrac{1}{2}\parallel \m - \m_{\text prior} \parallel^{2}_{\bf \Gamma_{\text{prior}}^{-1}} \right)
$

### The Gaussian approximation of the posterior: $\mathcal{N}(\vec{\map},\bf \Gamma_{\text{post}})$

The mean of this posterior distribution, $\vec{\map}$, is the
parameter vector maximizing the posterior, and
is known as the maximum a posteriori (MAP) point.  It can be found
by minimizing the negative log of the posterior, which amounts to
solving a deterministic inverse problem) with appropriately weighted norms,

$
\map := \underset{\m}{\arg \min} \; \mathcal{J}(\m) \;:=\;
\Big( 
-\frac{1}{2} \| \bf f(\m) - \data \|^2_{\bf \Gamma_{\text{noise}}^{-1}} 
-\frac{1}{2} \| \m -\m_{\text prior} \|^2_{\bf \Gamma_{\text{prior}}^{-1}} 
\Big).
$

The posterior covariance matrix is then given by the inverse of
the Hessian matrix of $\mathcal{J}$ at $\map$, namely

$
\bf \Gamma_{\text{post}} = \left(\Hmisfit(\map) + \bf \Gamma_{\text{prior}}^{-1} \right)^{-1}
$

#### The prior-preconditioned Hessian of the data misfit:

$
  \HT := \prcov^{1/2} \Hmisfit \prcov^{1/2}
$

#### The generalized eigenvalue problem:

$
 \Hmisfit \matrix{V} = \prcov^{-1} \matrix{V} \matrix{\Lambda},
$

where $\matrix{\Lambda} = diag(\lambda_i) \in \mathbb{R}^{n\times n}$
contains the generalized eigenvalues and the columns of $\matrix V\in
\mathbb R^{n\times n}$ the generalized eigenvectors such that 
$\matrix{V}^T \prcov^{-1} \matrix{V} = \matrix{I}$.

#### Randomized eigensolvers to construct the approximate spectral decomposition:  

When the generalized eigenvalues $\{\lambda_i\}$ decay rapidly, we can
extract a low-rank approximation of $\Hmisfit$ by retaining only the $r$
largest eigenvalues and corresponding eigenvectors,

$
 \HT = \prcov^{-1/2} \matrix{V}_r \matrix{\Lambda}_r \matrix{V}^T_r \prcov^{-1/2},
$

Here, $\matrix{V}_r \in \mathbb{R}^{n\times r}$ contains only the $r$
generalized eigenvectors of $\Hmisfit$ that correspond to the $r$ largest eigenvalues,
which are assembled into the diagonal matrix $\matrix{\Lambda}_r = \diag
(\lambda_i) \in \mathbb{R}^{r \times r}$.

#### The approximate posterior covariance::

Using the Sherman–Morrison–Woodbury formula, we write
$$
\begin{align}
  \notag \postcov = \left(\Hmisfit+ \prcov^{-1}\right)^{-1}
  = \prcov^{-1}-\matrix{V}_r \matrix{D}_r \matrix{V}_r^T +
  \mathcal{O}\left(\sum_{i=r+1}^{n} \frac{\lambda_i}{\lambda_i +
    1}\right),
\end{align}
$$

where $\matrix{D}_r :=\diag(\lambda_i/(\lambda_i+1)) \in
\mathbb{R}^{r\times r}$. The last term in this expression captures the
error due to truncation in terms of the discarded eigenvalues; this
provides a criterion for truncating the spectrum, namely that $r$ is
chosen such that $\lambda_r$ is small relative to 1. 

Therefore we can approximate the posterior covariance as

$$
\postcov \approx \prcov - \matrix{V}_r \matrix{D}_r
\matrix{V}_r^T
$$

#### Drawing samples from a Gaussian distribution with covariance $\H^{-1}$

Let $\bf x$ be a sample for the prior distribution, i.e. $\bf x \sim \mathcal{N}({\bf 0}, \prcov)$, then, using the low rank approximation of the posterior covariance, we compute a sample ${\bf v} \sim \mathcal{N}({\bf 0}, \H^{-1})$ as 
$$
  {\bf v} = \big\{ \Vr \big[ (\matrix{\Lambda}_r +
    \Ir)^{-1/2} - \Ir \big] \Vr^T\prcov^{-1}  + \bf I \big\} {\bf x} 
$$

## This tutorial shows:

- description of the inverse problem (the forward problem, the prior, and the misfit functional)
- convergence of the inexact Newton-CG algorithm
- low-rank-based approximation of the posterior covariance (built on a low-rank
approximation of the Hessian of the data misfit) 
- how to construct the low-rank approximation of the Hessian of the data misfit
- how to apply the inverse and square-root inverse Hessian to a vector efficiently
- samples from the Gaussian approximation of the posterior

## Goals:

By the end of this notebook, you should be able to:

- understand the Bayesian inverse framework
- visualise and understand the results
- modify the problem and code

## Mathematical tools used:

- Finite element method
- Derivation of gradiant and Hessian via the adjoint method
- inexact Newton-CG
- Armijo line search
- Bayes' formula
- randomized eigensolvers

## List of software used:

- <a href="http://fenicsproject.org/">FEniCS</a>, a parallel finite element element library for the discretization of partial differential equations
- <a href="http://www.mcs.anl.gov/petsc/">PETSc</a>, for scalable and efficient linear algebra operations and solvers
- <a href="http://matplotlib.org/">Matplotlib</a>, A great python package that I used for plotting many of the results
- <a href="http://www.numpy.org/">Numpy</a>, A python package for linear algebra.  While extensive, this is mostly used to compute means and sums in this notebook.

## 1. Load modules


```
import dolfin as dl
import math
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import sys
sys.path.append( "../" )
from hippylib import *

import nb

import logging
logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
dl.set_log_active(False)

np.random.seed(seed=1)
```

## 2. Generate the true parameter

This function generates a random field with a prescribed anysotropic covariance function.


```
def true_model(Vh, gamma, delta, anis_diff):
    prior = BiLaplacianPrior(Vh, gamma, delta, anis_diff )
    noise = dl.Vector()
    prior.init_vector(noise,"noise")
    noise_size = noise.array().shape[0]
    noise.set_local( np.random.randn( noise_size ) )
    atrue = dl.Vector()
    prior.init_vector(atrue, 0)
    prior.sample(noise,atrue)
    return atrue
```

## 3. Set up the mesh and finite element spaces

We compute a two dimensional mesh of a unit square with nx by ny elements.
We define a P2 finite element space for the *state* and *adjoint* variable and P1 for the *parameter*.


```
ndim = 2
nx = 64
ny = 64
mesh = dl.UnitSquareMesh(nx, ny)
Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
Vh = [Vh2, Vh1, Vh2]
print "Number of dofs: STATE={0}, PARAMETER={1}, ADJOINT={2}".format(Vh[STATE].dim(), Vh[PARAMETER].dim(), Vh[ADJOINT].dim())
```

    Number of dofs: STATE=16641, PARAMETER=4225, ADJOINT=16641


## 4. Set up the forward problem

To set up the forward problem we use the `PDEVariationalProblem` class, which requires the following inputs
- the finite element spaces for the state, parameter, and adjoint variables `Vh`
- the pde in weak form `pde_varf`
- the boundary conditions `bc` for the forward problem and `bc0` for the adjoint and incremental problems.

The `PDEVariationalProblem` class offer the following functionality:
- solving the forward/adjoint and incremental problems
- evaluate first and second partial derivative of the forward problem with respect to the state, parameter, and adojnt variables.


```
def u_boundary(x, on_boundary):
    return on_boundary and ( x[1] < dl.DOLFIN_EPS or x[1] > 1.0 - dl.DOLFIN_EPS)

u_bdr = dl.Expression("x[1]")
u_bdr0 = dl.Expression("0.0")
bc = dl.DirichletBC(Vh[STATE], u_bdr, u_boundary)
bc0 = dl.DirichletBC(Vh[STATE], u_bdr0, u_boundary)

f = dl.Expression("0.0")
    
def pde_varf(u,a,p):
    return dl.exp(a)*dl.inner(dl.nabla_grad(u), dl.nabla_grad(p))*dl.dx - f*p*dl.dx
    
pde = PDEVariationalProblem(Vh, pde_varf, bc, bc0)
```

## 4. Set up the prior

To obtain the synthetic true paramter $a_{\rm true}$ we generate a realization of a Gaussian random field with zero average and covariance matrix $\mathcal{C} = \widetilde{\mathcal{A}}^{-2}$, where $\widetilde{\mathcal{A}}$ is a differential operator of the form
$$ \widetilde{\mathcal{A}} = \gamma {\rm div}\, \Theta\, {\rm grad} + \delta I. $$
Here $\Theta$ is an s.p.d. anisotropic tensor of the form
$$ \Theta =
\begin{bmatrix}
\theta_1 \sin(\alpha)^2 & (\theta_1-\theta_2) \sin(\alpha) \cos{\alpha} \\
(\theta_1-\theta_2) \sin(\alpha) \cos{\alpha} & \theta_2 \cos(\alpha)^2.
\end{bmatrix} $$

For the prior model, we assume that we can measure the log-permeability coefficient at $N$ locations, and we denote with $a^1_{\rm true}$, $\ldots$, $a^N_{\rm true}$ such measures.
We also introduce the mollifier functions
$$ \delta_i(x) = \exp\left( -\frac{\gamma^2}{\delta^2} \| x - x_i \|^2_{\Theta^{-1}}\right), \quad i = 1, \ldots, N,$$
and we let
$$ \mathcal{A} = \widetilde{\mathcal{A}} + p \sum_{i=1}^N \delta_i I = \widetilde{\mathcal{A}} + p \mathcal{M},$$
where $p$ is a penalization costant (10 for this problem) and $ \mathcal{M} = \sum_{i=1}^N \delta_i I$.

We then compute $a_{\rm pr}$, the  mean  of  the  prior  measure,  as  a  regularized
least-squares fit of these point observations by solving
$$
a_{\rm pr} = arg\min_{m} \frac{1}{2}\langle a, \widetilde{\mathcal{A}} a\rangle + \frac{p}{2}\langle a_{\rm true} - a, \mathcal{M}(a_{\rm true}- a) \rangle.
$$

Finally the prior distribution is $\mathcal{N}(a_{\rm pr}, \mathcal{C}_{\rm prior})$, with $\mathcal{C}_{\rm prior} = \mathcal{A}^{-2}$.


```
gamma = .1
delta = .5
    
anis_diff = dl.Expression(code_AnisTensor2D)
anis_diff.theta0 = 2.
anis_diff.theta1 = .5
anis_diff.alpha = math.pi/4
atrue = true_model(Vh[PARAMETER], gamma, delta,anis_diff)
        
locations = np.array([[0.1, 0.1], [0.1, 0.9], [.5,.5], [.9, .1], [.9, .9]])
pen = 1e1
prior = MollifiedBiLaplacianPrior(Vh[PARAMETER], gamma, delta, locations, atrue, anis_diff, pen)
      
print "Prior regularization: (delta_x - gamma*Laplacian)^order: delta={0}, gamma={1}, order={2}".format(delta, gamma,2)    
            
objs = [dl.Function(Vh[PARAMETER],atrue), dl.Function(Vh[PARAMETER],prior.mean)]
mytitles = ["True Parameter", "Prior mean"]
nb.multi1_plot(objs, mytitles)
plt.show()

model = Model(pde,prior, misfit)
```

    Prior regularization: (delta_x - gamma*Laplacian)^order: delta=0.5, gamma=0.1, order=2



![png](output_10_1.png)


## 5. Set up the misfit functional and generate synthetic observations

To setup the observation operator, we generate *ntargets* random locations where to evaluate the value of the state.

To generate the synthetic observation, we first solve the forward problem using the true parameter $a_{\rm true}$. Synthetic observations are obtained by perturbing the state variable at the observation points with a random gaussian noise.
*rel_noise* is the signal to noise ratio.


```
ntargets = 300
rel_noise = 0.01


targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
print "Number of observation points: {0}".format(ntargets)
misfit = PointwiseStateObservation(Vh[STATE], targets)

utrue = pde.generate_state()
x = [utrue, atrue, None]
pde.solveFwd(x[STATE], x, 1e-9)
misfit.B.mult(x[STATE], misfit.d)
MAX = misfit.d.norm("linf")
noise_std_dev = rel_noise * MAX
randn_perturb(misfit.d, noise_std_dev)
misfit.noise_variance = noise_std_dev*noise_std_dev

vmax = max( utrue.max(), misfit.d.max() )
vmin = min( utrue.min(), misfit.d.min() )

plt.figure(figsize=(15,5))
nb.plot(dl.Function(Vh[STATE], utrue), mytitle="True State", subplot_loc=121, vmin=vmin, vmax=vmax)
nb.plot_pts(targets, misfit.d, mytitle="Observations", subplot_loc=122, vmin=vmin, vmax=vmax)
plt.show()
```

    Number of observation points: 300



![png](output_12_1.png)


## 6. Set up the model and test gradient and Hessian

The model is defined by three component:
- the `PDEVariationalProblem` `pde` which provides methods for the solution of the forward problem, adjoint problem, and incremental forward and adjoint problems.
- the `Prior` `prior` which provides methods to apply the regularization (*precision*) operator to a vector or to apply the prior covariance operator (i.e. to solve linear system with the regularization operator)
- the `Misfit` `misfit` which provides methods to compute the cost functional and its partial derivatives with respect to the state and parameter variables.

To test gradient and the Hessian of the model we use forward finite differences.


```
model = Model(pde, prior, misfit)

a0 = dl.interpolate(dl.Expression("sin(x[0])"), Vh[PARAMETER])
modelVerify(model, a0.vector(), 1e-12)
```

    (yy, H xx) - (xx, H yy) =  9.09491932242e-15



![png](output_14_1.png)


## 7. Compute the MAP point

We used the globalized Newtown-CG method to compute the MAP point.


```
a0 = prior.mean.copy()
solver = ReducedSpaceNewtonCG(model)
solver.parameters["rel_tolerance"] = 1e-9
solver.parameters["abs_tolerance"] = 1e-12
solver.parameters["max_iter"]      = 25
solver.parameters["inner_rel_tolerance"] = 1e-15
solver.parameters["c_armijo"] = 1e-4
solver.parameters["GN_iter"] = 5
    
x = solver.solve(a0)
    
if solver.converged:
    print "\nConverged in ", solver.it, " iterations."
else:
    print "\nNot Converged"

print "Termination reason: ", solver.termination_reasons[solver.reason]
print "Final gradient norm: ", solver.final_grad_norm
print "Final cost: ", solver.final_cost

plt.figure(figsize=(15,5))
nb.plot(dl.Function(Vh[STATE], x[STATE]), subplot_loc=121,mytitle="State")
nb.plot(dl.Function(Vh[PARAMETER], x[PARAMETER]), subplot_loc=122,mytitle="Parameter")
plt.show()
```

    
    It  cg_it cost            misfit          reg             (g,da)          ||g||L2        alpha          tolcg         
      1   1    1.205749e+03    1.205435e+03    3.147595e-01   -1.569088e+04   1.041993e+05   1.000000e+00   5.000000e-01
      2   3    3.456819e+02    3.444282e+02    1.253761e+00   -1.845351e+03   1.430874e+04   1.000000e+00   3.705684e-01
      3   1    2.745939e+02    2.732846e+02    1.309355e+00   -1.421202e+02   1.002730e+04   1.000000e+00   3.102127e-01
      4   7    1.691715e+02    1.647563e+02    4.415191e+00   -2.126866e+02   3.868977e+03   1.000000e+00   1.926929e-01
      5   6    1.573196e+02    1.522914e+02    5.028129e+00   -2.345175e+01   1.820008e+03   1.000000e+00   1.321613e-01
      6  14    1.424898e+02    1.297463e+02    1.274356e+01   -2.988290e+01   1.157435e+03   1.000000e+00   1.053940e-01
      7   2    1.421591e+02    1.294122e+02    1.274692e+01   -6.608825e-01   7.299626e+02   1.000000e+00   8.369856e-02
      8  22    1.407910e+02    1.253199e+02    1.547109e+01   -2.732846e+00   4.407936e+02   1.000000e+00   6.504072e-02
      9  16    1.407819e+02    1.253760e+02    1.540587e+01   -1.827524e-02   5.039967e+01   1.000000e+00   2.199285e-02
     10  29    1.407814e+02    1.253304e+02    1.545101e+01   -9.452146e-04   1.061811e+01   1.000000e+00   1.009465e-02
     11  36    1.407814e+02    1.253304e+02    1.545106e+01   -8.147549e-08   9.003262e-02   1.000000e+00   9.295390e-04
     12  62    1.407814e+02    1.253304e+02    1.545106e+01   -2.063794e-12   3.281322e-04   1.000000e+00   5.611670e-05
    
    Converged in  12  iterations.
    Termination reason:  Norm of the gradient less than tolerance
    Final gradient norm:  1.27740754831e-08
    Final cost:  140.781419038



![png](output_16_1.png)


## 8. Compute the low rank Gaussian approximation of the posterior
We used the *double pass* algorithm to compute a low-rank decomposition of the Hessian Misfit.
In particular, we solve

$$ \Hmisfit {\bf u} = \lambda \prcov^{-1} {\bf u}. $$

The Figure shows the largest *k* generalized eigenvectors of the Hessian misfit.
The effective rank of the Hessian misfit is the number of eigenvalues above the red line (y=1).
The effective rank is independent of the mesh size.


```
model.setPointForHessianEvaluations(x)
Hmisfit = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=True)
k = 50
p = 20
print "Single/Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)
Omega = np.random.randn(x[PARAMETER].array().shape[0], k+p)
d, U = doublePassG(Hmisfit, prior.R, prior.Rsolver, Omega, k)

posterior = GaussianLRPosterior(prior, d, U)
posterior.mean = x[PARAMETER]

plt.plot(range(0,k), d, 'b*', range(0,k+1), np.ones(k+1), '-r')
plt.yscale('log')
plt.xlabel('number')
plt.ylabel('eigenvalue')

nb.plot_eigenvectors(Vh[PARAMETER], U, mytitle="Eigenvector", which=[0,1,2,5,10,15])
```

    Single/Double Pass Algorithm. Requested eigenvectors: 50; Oversampling 20.



![png](output_18_1.png)



![png](output_18_2.png)


## 9. Prior and posterior pointwise variance fields


```
compute_trace = True
if compute_trace:
    post_tr, prior_tr, corr_tr = posterior.trace(method="Estimator", tol=5e-2, min_iter=20, max_iter=2000)
    print "Posterior trace {0:5e}; Prior trace {1:5e}; Correction trace {2:5e}".format(post_tr, prior_tr, corr_tr)
post_pw_variance, pr_pw_variance, corr_pw_variance = posterior.pointwise_variance("Exact")

objs = [dl.Function(Vh[PARAMETER], pr_pw_variance),
        dl.Function(Vh[PARAMETER], post_pw_variance)]
mytitles = ["Prior variance", "Posterior variance"]
nb.multi1_plot(objs, mytitles, logscale=True)
plt.show()
```

    Posterior trace 1.144005e-01; Prior trace 4.031887e-01; Correction trace 2.887882e-01



![png](output_20_1.png)


## 10. Generate samples from Prior and Posterior


```
nsamples = 5
noise = dl.Vector()
posterior.init_vector(noise,"noise")
noise_size = noise.array().shape[0]
s_prior = dl.Function(Vh[PARAMETER], name="sample_prior")
s_post = dl.Function(Vh[PARAMETER], name="sample_post")

range_pr = 2*math.sqrt( pr_pw_variance.max() )
ps_max   = 2*math.sqrt( post_pw_variance.max() ) + posterior.mean.max()
ps_min   = -2*math.sqrt( post_pw_variance.max() ) + posterior.mean.min()

for i in range(nsamples):
    noise.set_local( np.random.randn( noise_size ) )
    posterior.sample(noise, s_prior.vector(), s_post.vector())
    plt.figure(figsize=(15,5))
    nb.plot(s_prior, subplot_loc=121,mytitle="Prior sample", vmin=-range_pr, vmax=range_pr)
    nb.plot(s_post, subplot_loc=122,mytitle="Posterior sample", vmin=ps_min, vmax=ps_max)
    plt.show()
```


![png](output_22_0.png)



![png](output_22_1.png)



![png](output_22_2.png)



![png](output_22_3.png)



![png](output_22_4.png)


Copyright (c) 2016, The University of Texas at Austin & University of California, Merced.
All Rights reserved.
See file COPYRIGHT for details.

This file is part of the hIPPYlib library. For more information and source code
availability see https://hippylib.github.io.

hIPPYlib is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (as published by the Free Software Foundation) version 3.0 dated June 2007.
