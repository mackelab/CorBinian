![Modified from Aups under GFDL , via Wikimedia Commons, https://commons.wikimedia.org/wiki/File:Heraldique_ours_passant2.svg](https://user-images.githubusercontent.com/8633855/29385509-462d1f06-82d8-11e7-84bf-1e97ffadd4f2.png "Bear, modified from Aups under GFDL , via Wikimedia Commons, https://commons.wikimedia.org/wiki/File:Heraldique_ours_passant2.svg")

# CorBinian: A toolbox for modelling and simulating high-dimensional binary data and count-data with correlations

CorBinian is a toolbox for statistical modelling of multivariate binary and count data with correlations. 

CorBinian supports multiple statistical models, including different variants of maximum entropy models (with fitting algorithms that scale to >100 dimensions), Dichotomized Gaussian models, as well as tools for estimating various statistical quantities of interest from these models. It has been developed and implemented with the goal of modelling spike-train recordings from neural populations, but many of the methods will be applicable more generally. 

This repository is based on --- and replaces --- the original toolbox [pop_spike](https://bitbucket.org/mackelab/pop_spike).


We tried to make conventions compatible across the different projects and to share utility functions, but there is still some mismatch in conventions and redundant functions.

See [github.com/mackelab](https://github.com/mackelab) or [mackelab.org/code](http://www.mackelab.org/code) for more repositories by the group.

## Documentation

### Maximum entropy modeling 

This repository contains efficient code for fitting maximum entropy models to multivariate binary data, including 
- second-order (a.k.a. Ising) models
- K-pairwise models (second order plus population-count terms)
- flat or K-synchronous models (only population-count terms, assuming homogeneous population)

Most of the code can be used with any other (binary) features computed from the data. 
Our implementations use MCMC and iterative scaling to scale to large (N > 100) populations. 

### Dichotomized Gaussian models

This repository contains also code for fitting dichotomized Gaussian models to multivariate binary data. 
Particularly fast code is available for 'homogeneous' models with shared mean- and interaction parameters across the population. 

### specific heat analysis

This repository also provides functions for computing/estimating the variance of log probabilities and related specific heat capacity from multivariate binary distributions.  

## Usage

To get started, change base_dir.m to the name of the directory that the code is sitting in, run startup.m to set the path, and run one of the demo-files in the demo-folder. 

**We are in the process of integrating new functions for 'K-pairwise' and 'flat'/'K-synchronous' maximum entropy models into the codebase. Pleas note that until this is completed, there may be some compatibility issues.** 
A stable version of many of the spike-train analyses featured here can also be found at [bitbucket.org/mackelab/pop_spike/](https://bitbucket.org/mackelab/pop_spike/src)

Much of this code was developed in collaboration with---or even by--- the co-authors on the various manuscript. Feel free to use the code, but please acknowledge the source and paper appropriately if you are using it for a publication. 

If you notice a bug, want to request a feature, or have a question or feedback, please make use of the issue-tracking capabilities of the repository. We love to hear from people using our code-- please send an email to info@mackelab.org.

The code is published under the GNU General Public License. The code is provided "as is" and has no warranty whatsoever. 

## Publications

The repository contains implementations of the methods presented in 

###  [M Nonnenmacher, C Behrens, P Berens, M Bethge, JH Macke: Signatures of criticality arise in simple neural population models with correlations, arXiv:1603.00097](http://www.mackelab.org/publications/#signatures-of-criticality-arise-in-simple-neural-population-models-with-correlations)

Currently implemented: Code for fast iterative scaling for K-pairwise maximum entropy models ([Tkacik et al. 2014](https://doi.org/10.1371/journal.pcbi.1003408)), based on pairwise Gibbs sampling with Rao-Blackwellized estimators for first and second moments of the distribution. 

* demo_iter_scaling.m: Fit a maximum entropy model (typically K-pairwise but code is flexible) using regularized fast iterative scaling with Rao-Blackwellization. 

* demo_MPF.m: Fit a K-pairwise maximum entropy model using minimum probability flow. 

###  [JH Macke*, P Berens*, AS Ecker, AS Tolias and M Bethge: Generating Spike Trains with Specified Correlation Coefficients. Neural Computation 21(2), 397-423, 02 2009](http://www.mackelab.org/publications/#signatures-of-criticality-arise-in-simple-neural-population-models-with-correlations)

Currently implemented: Functions for fitting and sampling from dichotomised Gaussian models both with binary and count observations. Some of this code is equivalent to a previous toolbox written primarily by P Berens and JH Macke in the lab of M Bethge, the original code-package can be found at http://bethgelab.org/software/mvd/.

Currently implemented: Functions for fitting and analysing dichotomised Gausian models to 'homogeneous' models in which all neurons are assumed to have the same firing rate and pairwise correlation, also some first functions for 'semi-homogeneous' models in which neurons are allowed to have different firing rates.

* demo_dich_gauss_01.m: Using dichotomised Gaussian on binary random variables
* demo_dich_gauss_counts.m: Using discretized Gaussian on correlated counts with arbitrary marginal distributions

### [JH Macke, M Opper, M Bethge: Common input explains higher-order correlations and entropy in a simple model of neural population activity. Physical Review Letters 106, 208102, 05 2011](http://www.mackelab.org/publications/#common-input-explains-higher-order-correlations-and-entropy-in-a-simple-model-of-neural-population-activity)

* demo_flat_models.m: Using the dichotomised Gaussian on homogeneous population models, i.e. models in which all neurons are assumed to to have the same mean firing rate and same pairwise correlation.

### [G Schwartz, JH Macke, D Amodei, H Tang, MJ Berry: Low error discrimination using a correlated population code. Journal of Neurophysiology, 108(4), 1069-1088, 04 2012](http://www.mackelab.org/publications/#low-error-discrimination-using-a-correlated-population-code)

Currently implemented: Code for fitting binary second order maximum entropy models (Ising models) which also works for large populations of neurons (N>100). Much of this code is based on a previous implemented by Tamara Broderick et al, http://arxiv.org/abs/0712.2437.

* demo_maxent_MCMC.m: Fit a second order maximum entropy model to a large population of neurons using MCMC and iterative scaling.

### [JH Macke, I Murray, P Latham: Estimation bias in maximum entropy models. Entropy 15:3109-3219, 08 2013](http://www.mackelab.org/publications/#estimation-bias-in-maximum-entropy-models)

Currently implemented: Functions for fitting maximum entropy models to small populations of neurons (N<15)

Missing: Methods for calculating bias 

* demo_maxent.m: Fit a maximum entropy model (typically second order but code is flexible) in a case which can be solved exactly, i.e. for which one does not need MCMC.



Additional code used specifically for [Nonnemacher et al. (2016)](https://arxiv.org/abs/1603.00097) can be found [here](https://github.com/mackelab/critical_retina). This includes scripts for fitting K-pairwise models to simulated population spike-data, and code for computing population measures such as the specific heat capacity from these models.
