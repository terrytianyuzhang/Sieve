# Sieve

References:https://arxiv.org/abs/2206.02994

R package, Sieve. Perform nonparametric estimation by the method of sieves (estimation using multivariate orthogonal series). This type of estimators has been actively studied and applied in univariate feature settings, but in multivariate cases it hasn't received its deserved attention. 

Installing a package from GitHub can be tricky. But I found 80% of the errors can be solved by restarting RStudio.

The current version can solve regression and classification problems. The algorithm gives the estimated condition mean (regression) and estimated conditional probability functions (classification). I will make it able to handle time-to-event outcomes very soon.

## Computationally tractable: 
The time and space expense both scale linearly in sample size and the number of basis functions specified by the users. Can directly handle 10k x 100 (sample size x dimension of features) data science problems.

## Theoretically guaranteed: 
Adaptive to the number of features/predictors truly associated with the outcome. Can achieve the information lower bounds (minimax rate) of estimation in many cases. 

## What is penalized sieve estimation? 
Generating the proper basis functions (something like multivariate Fourier basis), put everything in a LASSO solver (thank you glmnet!). That's it. 

(Questions, suggestion, collaboration: shoot me an email: zty@uw.edu, Tianyu Zhang. Department of Biostatistics, University of Washington)
