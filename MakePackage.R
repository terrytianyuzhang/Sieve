library(Rcpp)
library(RcppArmadillo)
library(tictoc)

setwd('/Users/tianyuzhang/Documents/SievePackage/Sieve/')
# setwd('/Users/tianyuzhang/Documents/SievePackage/Sieve/src/')
# source('/Users/tianyuzhang/Documents/SievePackage/Sieve/R/SieveFittingModels.R')
# sourceCpp("PracticeRcppFunction.cpp")
# sourceCpp('/Users/tianyuzhang/Documents/SievePackage/Sieve/src/C_Functions.cpp')
# 
# # RcppArmadillo.package.skeleton('Sieve', cpp_files= '/Users/tianyu/Documents/SievePackage/', example_code = T) 
# RcppArmadillo.package.skeleton('Sieve', code_files = c('TensorProductFunction.R'))
# Rcpp::compileAttributes('Sieve')
# Now, move the cpp file into scr/

# To build a newer version of a package, I just need to run the following steps, nothing more
devtools::document()
devtools::build()
devtools::install()
library(Sieve)

# I DON'T NEED TO DO THE FOLLOWING AFTER RUNNING THE ABOVE devtools command
# Go to cmd, type 
# R CMD build Sieve/
# Then, 
# R CMD INSTALL Sieve/
