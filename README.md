# Gradient-based marginal optimization

__gmo__ is an R package for fast optimization of marginal posterior
distributions. Using a stochastic gradient-based algorithm, gmo
estimates a set of parameters from a model while marginalizing out the
rest. This provides uncertainty over any nuisance parameters, and
generalizes parameter estimation using marginal densities. It acts as
a middleground between full Bayesian inference over all parameters and
point estimation over all parameters.

Here is an example for a mixture model from
[Rubin's 8 schools analysis (1981)](http://jeb.sagepub.com/content/6/4/377.short).
Each data point belongs to one of 8 groups, and gmo estimates the mean
and variance parameters of the prior on each group.
```R
library(gmo)
library(rstan)

data <- list(J = 8,
             K = 2,
             y = c(28,  8, -3,  7, -1,  1, 18, 12),
             sigma = c(15, 10, 16, 11,  9, 11, 10, 18))

fit.gmo <- gmo("models/8schools.stan", "models/8schools_local.stan", data=data)
```
The two Stan programs used above are found [here](demo/models/8schools.stan) and
[here](demo/models/8schools_local.stan). More examples are found in [`demo/`](demo/).

The gmo package uses a modified Stan program in order to distinguish
between parameters to estimate and parameters to marginalize out.
<!--[A tutorial for writing this Stan program is available-->
<!--here](https://github.com/gelman/gmo/wiki/Tutorial).   -->

## Features

The core feature of gmo is a fast way to maximize marginal posterior
densities, which includes, for example, maximum marginal likelihood, empirical
Bayes, and type II maximum likelihood. It is done in an iterative
scheme that is closely inspired by the EM algorithm. Here are
additional features it supports:

+ Support for all models written in [Stan](http://mc-stan.org)
<!--+ Specialized algorithms for mixed-effects models in-->
<!--  [lme4](https://github.com/lme4/lme4)              -->
+ Uncertainty using covariance estimates
+ Approximate maximum marginal likelihood
+ Penalized maximum marginal likelihood
+ Fully maximum marginal likelihood
+ Bayesian inference with data-dependent priors

## Installation
<!--(TODO not submitted to CRAN yet)-->
<!--To install the latest version from CRAN:                           -->
<!--```R                                                               -->
<!--install.packages("gmo")                                            -->
<!--```                                                                -->
GMO is experimental software and undergoing development. We plan to
submit to CRAN once it is ready.

To install the latest development version from Github:
```R
# install.packages("devtools")
devtools::install_github("stan-dev/gmo")
```

## Citation

We appreciate citations for GMO if you apply or build off it in your work.

+ Dustin Tran, Andrew Gelman, and Aki Vehtari. 2016. Gradient-based marginal optimization. In preparation.

```
@article{tran2016gmo,
  title = {Gradient-based marginal optimization},
  author = {Dustin Tran and Andrew Gelman and Aki Vehtari},
  journal = {In preparation},
  year = {2016}
}
```
