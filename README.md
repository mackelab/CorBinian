# Welcome

This repository contains a number of different statistical methods for modelling multivariate binary and count data with correlations. 

It has been developed and implemented with the goal of modelling spike-train recordings from neural populations, but at least some of the methods will be applicable more generally. I have tried to make conventions compatible across the different projects and to share utility functions, but there is still some mismatch in conventions and redundant functions.

The current version is (and to some degree will always be) work in progress. 

See [https://bitbucket.org/mackelab/home](https://bitbucket.org/mackelab/home) for more repositories by the group.

## Usage

To get started, change base_dir.m to the name of the directory that the code is sitting in, run startup.m to set the path, and run one of the demo-files in the demo-folder. 


Much of this code was developed in collaboration with---or even by--- the co-authors on the various manuscript. Feel free to use the code, but please acknowledge the source and paper appropriately if you are using it for a publication. 

If you notice a bug, want to request a feature, or have a question or feedback, please make use of the issue-tracking capabilities of the repository. We love to hear from people using our code-- please send an email to info@mackelab.org.

The code is published under the GNU General Public License. The code is provided "as is" and has no warranty whatsoever. 

## Publications

(At least eventually), the repository will contain implementations of the methods presented in 

###  JH Macke*, P Berens*, AS Ecker, AS Tolias and M Bethge: Generating Spike Trains with Specified Correlation Coefficients. Neural Computation 21(2), 397-423, 02 2009

Currently implemented: Functions for fitting and sampling from dichotomised Gaussian models both with binary and count observations. Some of this code is equivalent to a previous toolbox written primarily by P Berens and JH Macke in the lab of M Bethge, the original code-package can be found at http://bethgelab.org/software/mvd/.

Currently implemented: Functions for fitting and analysing dichotomised Gausian models to 'homogeneous' models in which all neurons are assumed to have the same firing rate and pairwise correlation, also some first functions for 'semi-homogeneous' models in which neurons are allowed to have different firing rates.

Missing: Methods for asymptotic calculations, heat-capacity and Ising models. 


* demo_dich_gauss_01.m: Using dichotomised Gaussian on binary random variables
* demo_dich_gauss_counts.m: Using discretized Gaussian on correlated counts with arbitrary marginal distributions

### JH Macke, M Opper, M Bethge: Common input explains higher-order correlations and entropy in a simple model of neural population activity. Physical Review Letters 106, 208102, 05 2011

* demo_flat_models.m: Using the dichotomised Gaussian on homogeneous population models, i.e. models in which all neurons are assumed to to have the same mean firing rate and same pairwise correlation.


### G Schwartz, JH Macke, D Amodei, H Tang, MJ Berry: Low error discrimination using a correlated population code. Journal of Neurophysiology, 108(4), 1069-1088, 04 2012

Currently implemented: Code for fitting binary second order maximum entropy models (Ising models) which also works for large populations of neurons (N>100). Much of this code is based on a previous implemented by Tamara Broderick et al, http://arxiv.org/abs/0712.2437.

* demo_maxent_MCMC.m: Fit a second order maximum entropy model to a large population of neurons using MCMC and iterative scaling.

### JH Macke, I Murray, P Latham: Estimation bias in maximum entropy models. Entropy 15:3109-3219, 08 2013

Currently implemented: Functions for fitting maximum entropy models to small populations of neurons (N<15)

Missing: Methods for calculating bias 

* demo_maxent.m: Fit a maximum entropy model (typically second order but code is flexible) in a case which can be solved exactly, i.e. for which one does not need MCMC.