# Sieve
R package, Sieve. Perform nonparametric estimation by the method of sieves (generalized Fourier series).

The current version can solve regression and classification problems. I will make it able to handle time-to-event outcomes very soon.

I am planning to add some visulization functions for feature selection.

## Computationally tractable: 
The time and space expense both scale linearly in sample size and the number of basis functions specified by the users. Can directly handle 10k x 100 (sample size x dimension of features) data science problems.

## Theoretically guaranteed: 
Adaptive to the number of features/predictors truly associated with the outcome. Can achieve the information lower bounds (minimax rate) of estimation in many cases.

What we are doing: generating basis functions (something like multivariate Fourier basis), put everything in a LASSO solver (thank you glmnet!). That's it. 
