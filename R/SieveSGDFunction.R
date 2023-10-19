#' Preprocess the original data for sieve-SGD estimation.
#' 
#' @param X a data frame containing prediction features/ independent variables. The (i,j)-th element is the j-th dimension of the i-th sample's feature vector. 
#' So the number of rows equals to the sample size and the number of columns equals to the feature/covariate dimension. If the complete data set is large, this can be a representative subset of it (ideally have more than 1000 samples). 
#' @param s numerical array. Smoothness parameter, a smaller s corresponds to a more flexible model. Default is 2. The elements of this array should take values greater than 0.5. The larger s is, the smoother we are assuming the truth to be.
#' @param r0 numerical array. Initial learning rate/step size, don't set it too large. The step size at each iteration will be r0*(sample size)^(-1/(2s+1)), which is slowly decaying.
#' @param J numerical array. Initial number of basis functions, a larger J corresponds to a more flexible estimator The number of basis functions at each iteration will be J*(sample size)^(1/(2s+1)), which is slowly increasing. We recommend use J that is at least the dimension of predictor, i.e. the column number of the X matrix.
#' @param type a string. It specifies what kind of basis functions are used. The default is (aperiodic) cosine basis functions ('cosine'), which is enough for generic usage.
#' @param interaction_order a number. It also controls the model complexity. 1 means fitting an additive model, 2 means fitting a model allows, 3 means interaction terms between 3 dimensions of the feature, etc. The default is 3. 
#' For large sample size, lower dimension problems, try a larger value (but need to be smaller than the dimension of original features); for smaller sample size and higher dimensional problems, try set it to a smaller value (1 or 2).
#' @param omega the rate of dimension-reduction parameter. Default is 0.51, usually do not need to change.
#' @param norm_feature a logical variable. Default is TRUE. It means sieve_preprocess will rescale the each dimension of features to 0 and 1. Only set to FALSE when user already manually rescale them between 0 and 1.
#' @param norm_para a matrix. It specifies how the features are normalized. For training data, use the default value NULL.
#' @param lower_q lower quantile used in normalization. Default is 0.01 (1\% quantile).
#' @param upper_q upper quantile used in normalization. Default is 0.99 (99\% quantile).
#' 
#' @return A list containing the necessary information for next step model fitting. Typically, the list is used as the main input of sieve.sgd.solver.
#' \item{s.size.sofar}{a number. Number of samples has been processed so far.}
#' \item{type}{a string. The type of basis funtion.}
#' \item{hyper.para.list}{a list of hyperparameters.}
#' \item{index.matrix}{a matrix. Identifies the multivariate basis functions used in fitting.}
#' \item{index.row.prod}{the index product for each basis function. It is used in calculating basis function - specific learning rates.}
#' \item{inf.list}{a list storing the fitted results. It has a length of "number of unique combinations of the hyperparameters". The component of inf.list is itself a list, it has a hyper.para.index domain to specify its corresponding hyperparameters (need to be used together with hyper.para.list). Its rolling.cv domain is the progressive validation statistics for hyperparameter tuning; beta.f is the regression coefficients for the first length(beta.f) basis functions, the rest of the basis have 0 coefficients.}
#' \item{norm_para}{a matrix. It records how each dimension of the feature/predictor is rescaled, which is useful when rescaling the testing sample's predictors.}

#' @examples 
#' xdim <- 1 #1 dimensional feature
#' #generate 1000 training samples
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim)
#' sieve.model <- sieve.sgd.preprocess(X = TrainData[,2:(xdim+1)])
#' @export
#'
sieve.sgd.preprocess <- function(X, s = c(2), r0 = c(2), J = c(1), type = c('cosine'), 
                             interaction_order = c(3), omega = c(0.51),
                             norm_feature = TRUE, norm_para = NULL,
                             lower_q = 0.005, upper_q = 0.995){
  
  #s.size is sample size
  if(is(X,'numeric') | is(X, 'factor')){
    s.size <- length(X) #this is a univariate regression problem, sometimes the input data.frame is automatically conveted to an array
  }
  else{
    s.size <- dim(X)[1]
  }
  
  
  X <- sapply(X,unclass) #convert categorical variables into numeric variables
  X <- matrix(X, nrow = s.size)
  
  if(norm_feature){
    norm_list <- normalize_X(X, norm_para = norm_para, 
                                     lower_q = lower_q, upper_q = upper_q)
    # X <- norm_list$X
    norm_para <- norm_list$norm_para
  }
  
  #####set the number of observations counter back to 0
  s.size.sofar <- 0
  
  #####hyperparameter list
  hyper.para.list <- list(s = s, r0 = r0,
                          J = J, interaction_order = interaction_order,
                          omega = omega) #omega is related to a basis function-specific learning rate
  
  #####generate hyperparameter index matrix
  hyper.para.index.matrix <- expand.grid(1:length(s), 1:length(r0), 1:length(J),
                                         1:length(interaction_order), 1:length(omega))
  colnames(hyper.para.index.matrix) <- c('s', 'r0', 'J', 'interaction_order', 'omega')
  
  #####generate basis index matrix
  max.basisN <- ceiling( max(J) * s.size^(1/(2*min(s) + 1))
                        ) #the maximum number of basis functino needed to process all the data in X
  xdim <- dim(X)[2] #dimension of predictors
  index_matrix <- as.matrix(create_index_matrix(xdim, basisN = max.basisN, #the dimension of the index_matrix is specified by D (not xdim)
                                                interaction_order = max(interaction_order))[1:max.basisN,], 
                            ncol = xdim)
  index.row.prod <- index_matrix[, 1]#index product will be used when determining the basis function-specific learning rate
  index_matrix <- as.matrix(index_matrix[, -1], ncol = xdim)
  
  ####generate the list storing fitted beta and rolling cross-validation stat
  M <- dim(hyper.para.index.matrix)[1] #number of hyperparameter combinations
  inf.list <- vector("list", M) #for each combination of hyperparameter, create a list for storing the estiamted beta,etc.
  for(i in 1:M){
    tmp.inf <- list(hyper.para.index = hyper.para.index.matrix[i,],
                    rolling.cv = 0,
                    beta.f.int = 0, #beta of the internal fhat
                    beta.f = 0) #this is the estimator's beta  
    inf.list[[i]] <- tmp.inf
  }
  
  return(list(s.size.sofar = s.size.sofar,
              type = type,
              hyper.para.list = hyper.para.list,
              index.matrix = index_matrix,
              index.row.prod = index.row.prod,
              inf.list = inf.list,
              norm.para = norm_para))
}

#' Fit sieve-SGD estimators, using progressive validation for hyperparameter tuning.
#' 
#' @param sieve.model a list initiated using sieve.sgd.preprocess. Check the documentation of sieve.sgd.preprocess for more information.
#' @param X a data frame containing prediction features/ independent variables. 
#' @param Y training outcome.
#' @param cv_weight_rate this governs the divergence rate of rolling validation statistics. Default is set to be 1 and in general does not need to be changed.
#' 
#' @return A list. It contains the fitted regression coefficients and progressive validation statistics for each hyperparameter combination.
#' \item{s.size.sofar}{a number. Number of samples has been processed so far.}
#' \item{type}{a string. The type of basis funtion.}
#' \item{hyper.para.list}{a list of hyperparameters.}
#' \item{index.matrix}{a matrix. Identifies the multivariate basis functions used in fitting.}
#' \item{index.row.prod}{the index product for each basis function. It is used in calculating basis function - specific learning rates.}
#' \item{inf.list}{a list storing the fitted results. It has a length of "number of unique combinations of the hyperparameters". The component of inf.list is itself a list, it has a hyper.para.index domain to specify its corresponding hyperparameters (need to be used together with hyper.para.list). Its rolling.cv domain is the progressive validation statistics for hyperparameter tuning; beta.f is the regression coefficients for the first length(beta.f) basis functions, the rest of the basis have 0 coefficients.}
#' \item{norm_para}{a matrix. It records how each dimension of the feature/predictor is rescaled, which is useful when rescaling the testing sample's predictors.}

#' @examples 
#' frho.para <- xdim <- 1 ##predictor dimension
#' frho <- 'additive' ###truth is a sum of absolute functions 
#' type <- 'cosine' ###use cosine functions as the basis functions
#' #generate training data
#' TrainData <- GenSamples(s.size = 1e3, xdim = xdim, 
#'                                 frho.para = frho.para, 
#'                                 frho = frho, noise.para = 0.1)
#' #preprocess the model
#' sieve.model <- sieve.sgd.preprocess(X = TrainData[,2:(xdim+1)], 
#'                                     type = type,
#'                                     s = c(1,2),
#'                                     r0 = c(0.5, 2, 4),
#'                                     J = c(1, 4, 8))
#' 
#' ##train the model
#' sieve.model <- sieve.sgd.solver(sieve.model = sieve.model, 
#'                                 X = TrainData[,2:(xdim+1)], 
#'                                 Y  = TrainData[,1])
#' 
#' ##sieve-SGD can do multiple passes over the data, just like other SGD methods.
#' ##usually a second pass can still improve the prediction accuracy
#' ##watch out overfitting when performing multiple passes!
#' sieve.model <- sieve.sgd.solver(sieve.model = sieve.model, 
#'                               X = TrainData[,2:(xdim+1)], 
#'                               Y  = TrainData[,1])

#' @export
#' 
sieve.sgd.solver <- function(sieve.model, X, Y,
                             cv_weight_rate = 1){
  ####this function process the training X,Y using sieve-SGD
  ####also it uses rolling cross-validation to choose hyper-parameters
  
  #####normalize the X. Make sure the format of everything is correct.
  #s.size is sample size
  if(is(X,'numeric') | is(X, 'factor')){
    s.size <- length(X) #this is a univariate regression problem, sometimes the input data.frame is automatically conveted to an array
  }else{
    s.size <- dim(X)[1]
  }
  
  X <- sapply(X,unclass) #convert categorical variables into numeric variables
  X <- matrix(X, nrow = s.size)
  
  norm_para <- sieve.model$norm.para
  if(is.null(norm_para)){
    warning('there is no normalization parameter for X, 
            please use sieve.sgd.preprocess to preprocess the data.')  
  }
  
  norm_list <- normalize_X(X, norm_para = norm_para)
  X <- norm_list$X
  
  #####feed in the data one by one
  J <- sieve.model$hyper.para.list$J
  s <- sieve.model$hyper.para.list$s
  r0 <- sieve.model$hyper.para.list$r0
  omega <- sieve.model$hyper.para.list$omega
  interaction_order <- sieve.model$hyper.para.list$interaction_order
  type <- sieve.model$type
  M <- length(sieve.model$inf.list)
  s.size.sofar <- sieve.model$s.size.sofar
  
  if(s.size.sofar == 0){ #first time process the data
    index_matrix <- sieve.model$index.matrix
    index.row.prod <- sieve.model$index.row.prod
  }else{
    print('not training from the beginning')  
    
    #####regenerate basis index matrix
    max.basisN <- ceiling( max(J) * (s.size + s.size.sofar)^(1/(2*min(s) + 1))
    ) #the maximum number of basis functino needed to process all the data in X
    xdim <- dim(X)[2] #dimension of predictors
    index_matrix <- as.matrix(create_index_matrix(xdim, basisN = max.basisN, #the dimension of the index_matrix is specified by D (not xdim)
                                                          interaction_order = max(interaction_order))[1:max.basisN,], 
                              ncol = xdim)
    index.row.prod <- index_matrix[, 1]#index product will be used when determining the basis function-specific learning rate
    index_matrix <- as.matrix(index_matrix[, -1], ncol = xdim)
    
    
    sieve.model$index.matrix <- index_matrix
    sieve.model$index.row.prod <- index.row.prod
  }
  
  for(i in 1:s.size){
    
    i.sofar <- i + s.size.sofar #first pass, i.sofar = i
    
    newx <- matrix(X[i,], nrow = 1)
    newy <- Y[i]
    
    max.J <- ceiling(max(J) * (i.sofar)^(1/(2*min(s) + 1))) #this is the maximum number of basis functions need to estimated at this step
    
    Phi <- Design_M_C(newx, max.J, type, index_matrix) #one row design "matrix"
    
    for(m in 1:M){
      
      #####calculate rolling CV
      
      ###number of basis functions for this hyper-para combination
      J.m <- J[sieve.model$inf.list[[m]]$hyper.para.index$J] * (i.sofar)^(1/(2* s[sieve.model$inf.list[[m]]$hyper.para.index$s]+ 1))
      J.m <- ceiling(J.m)
      
      beta.f <- sieve.model$inf.list[[m]]$beta.f
      beta.f.int <- sieve.model$inf.list[[m]]$beta.f.int
      
      if(length(beta.f) < J.m){ #check if more basis are used
        beta.f <- c(beta.f, rep(0, J.m - length(beta.f)))
        beta.f.int <- c(beta.f.int, rep(0, J.m - length(beta.f.int)))
      }
      
      fnewx <- crossprod(beta.f, Phi[1:J.m]) ###f_i(X_{i+1})
      sieve.model$inf.list[[m]]$rolling.cv <- sieve.model$inf.list[[m]]$rolling.cv + i.sofar^(cv_weight_rate) * (newy - fnewx)^2
      
      ##########update beta's
      ###overall learning rate
      rn.m <- r0[sieve.model$inf.list[[m]]$hyper.para.index$r0] * (i.sofar)^(-1/(2* s[sieve.model$inf.list[[m]]$hyper.para.index$s]+ 1))
      
      #####CORE UPDATING STEPS
      fnewx.int <- crossprod(beta.f.int, Phi[1:J.m])
      res <- as.numeric(newy - fnewx.int)
      
      beta.f.int <- beta.f.int + rn.m * res * (index.row.prod[1:J.m])^(-2*omega)*Phi[1:J.m]
      
      beta.f <- (i.sofar-1)/i.sofar * beta.f + beta.f.int/i.sofar
      ################
      
      sieve.model$inf.list[[m]]$beta.f <- beta.f
      sieve.model$inf.list[[m]]$beta.f.int <- beta.f.int
      
    }
  }
  #update the current number of sample processed 
  sieve.model$s.size.sofar <- s.size.sofar + s.size
  
  # ######model comparison
  # rolling.cvs <- rep(1e10, M)
  # for(m in 1:M){
  #   rolling.cvs[m] <- sieve.model$inf.list[[m]]$rolling.cv
  # }
  # print(rolling.cvs/sieve.model$s.size.sofar)
  # print(which.min(rolling.cvs))
  sieve.model <- clean_up_result(sieve.model)
  return(sieve.model)

}

# xdim <- 2
# frho <- 'lineartensor'
# frho <- 'additive'
# TrainData <- GenSamples(s.size = 1e4, xdim = xdim, frho.para = 1, frho = frho, noise.para = 0.01)
# ###noise.para = 0.17 for lineartensor, snr 30, xdim = D = 4; =0.54, snr = 3
# X <- TrainData[,2:(xdim+1)]
# Y <- TrainData[,1]
# sieve.model <- sieve.sgd.preprocess(X = TrainData[,2:(xdim+1)], type = type,
#                                     s = c(1,2),
#                                     r0 = c(0.5, 2, 4),
#                                     J = c(0.5, 2, 4))
# sieve.model = sieve.sgd.solver(sieve.model = sieve.model, X = X, Y  = Y)
# sieve.model$inf.list[[13]]$hyper.para.index
# sieve.model$inf.list[[2]]
# sieve.model$inf.list[[4]]

# TestData <- GenSamples(s.size = 5e3, xdim = xdim, frho.para = 1, frho = frho, noise.para = 0)
# # X <- TestData[,2:(xdim+1)]
# # Y <- TestData[,1]
# sieve.model <- sieve.sgd.predict(sieve.model, TestData[,2:(xdim+1)])
# 
# mean((TestData[,1] - sieve.model$inf.list[[9]]$prdy)^2)
# mean((Y - sieve.model$inf.list[[15]]$prdy)^2)
# mean((Y - sieve.model$inf.list[[6]]$prdy)^2)
# plot(X, sieve.model$inf.list[[15]]$prdy)


#' Sieve-SGD makes prediction with new predictors.
#' 
#' @param sieve.model a list initiated using sieve.sgd.preprocess and sieve.sgd.solver. Check the documentation of sieve.sgd.preprocess for more information.
#' @param X a data frame containing prediction features/ independent variables. 

#' @return sieve.sgd.predict will update the given sieve.model input list.
#' \item{inf.list}{In each entry of the list inf.list, the array prdy is the predicted outcome under the given hyperparameter combination.}

#' @examples 
#' frho.para <- xdim <- 1 ##predictor dimension
#' frho <- 'additive' ###truth is a sum of absolute functions 
#' type <- 'cosine' ###use cosine functions as the basis functions
#' #generate training data
#' TrainData <- GenSamples(s.size = 1e3, xdim = xdim, 
#'                                 frho.para = frho.para, 
#'                                 frho = frho, noise.para = 0.1)
#' #preprocess the model
#' sieve.model <- sieve.sgd.preprocess(X = TrainData[,2:(xdim+1)], 
#'                                     type = type,
#'                                     s = c(1,2),
#'                                     r0 = c(0.5, 2, 4),
#'                                     J = c(1, 4, 8))
#' 
#' ##train the model
#' sieve.model <- sieve.sgd.solver(sieve.model = sieve.model, 
#'                                 X = TrainData[,2:(xdim+1)], 
#'                                 Y  = TrainData[,1])
#' ##generate new data
#' NewData <- GenSamples(s.size = 5e2, xdim = xdim, 
#'                       frho.para = frho.para, 
#'                       frho = frho, noise.para = 0.1)
#' sieve.model <- sieve.sgd.predict(sieve.model, X = NewData[, 2:(xdim+1)])
#' plot(NewData[, 2:(xdim+1)], sieve.model$best_model$prdy)
#' @export
#' 
sieve.sgd.predict <- function(sieve.model, X){
  
  #####normalize the X. Make sure the format of everything is correct.
  #s.size is sample size
  if(is(X,'numeric') | is(X, 'factor')){
    s.size <- length(X) #this is a univariate regression problem, sometimes the input data.frame is automatically conveted to an array
  }else{
    s.size <- dim(X)[1]
  }
  
  X <- sapply(X,unclass) #convert categorical variables into numeric variables
  X <- matrix(X, nrow = s.size)
  
  norm_para <- sieve.model$norm.para
  if(is.null(norm_para)){
    warning('there is no normalization parameter for X, 
            please use sieve.sgd.preprocess to preprocess the data.')  
  }
  
  norm_list <- normalize_X(X, norm_para = norm_para)
  X <- norm_list$X
  
  #####feed in the data one by one
  J <- sieve.model$hyper.para.list$J
  s <- sieve.model$hyper.para.list$s
  # r0 <- sieve.model$hyper.para.list$r0
  # omega <- sieve.model$hyper.para.list$omega
  index_matrix <- sieve.model$index.matrix
  # index.row.prod <- sieve.model$index.row.prod
  type <- sieve.model$type
  M <- length(sieve.model$inf.list)
  
  for(m in 1:M){
    #determine how many basis functions are needed
    J.m <- length(sieve.model$inf.list[[m]]$beta.f)
    
    Phi <- Design_M_C(X, J.m, type, index_matrix)
    
    tmp.prdy <- crossprod_C(Phi, matrix(sieve.model$inf.list[[m]]$beta.f))
    
    sieve.model$inf.list[[m]]$prdy <- tmp.prdy
  }
  
  sieve.model <- clean_up_result(sieve.model)
  print('find the best model prediction at')
  print('sieve.model$best_model$prdy')
  return(sieve.model)
  
}

#' Clean up the fitted model
#' 
#' @param sieve.model a sieve sgd model.
#' @return a processed sieve.model, adding function names and extract the best model
#' @export
clean_up_result <- function(sieve.model){
  index.matrix <- sieve.model$index.matrix
  index.row.prod <- sieve.model$index.row.prod
  basis_type <- sieve.model$type
  name_vector <- rep(NA, nrow(index.matrix))
  
  ##determine which model is the best according to rolling validation
  collect_cv_stat <- function(sieve.model){
    num_of_hyperpara <- length(sieve.model$inf.list)
    all_cv_stat <- NULL
    for(curr_hyperpara in 1:num_of_hyperpara){
      all_cv_stat <- c(all_cv_stat, sieve.model$inf.list[[curr_hyperpara]]$rolling.cv)
    }
    print('the best model index is')
    print(which.min(all_cv_stat))
    print('find it at sieve.model$best_model')
    return(which.min(all_cv_stat))
  }
  
  best_model_index <- collect_cv_stat(sieve.model)
  sieve.model$best_model <- sieve.model$inf.list[[best_model_index]]
  
  
  ####give the basis functions name
  spell_out_one_cosine_basis <- function(instruction_vector){
    basis_function_string <- NULL
    for(covariate_dimension in 1:length(instruction_vector)){
      
      if(instruction_vector[covariate_dimension] > 1){
        basis_function_string <- paste0(basis_function_string, '*cos(pi*',instruction_vector[covariate_dimension]- 1, '*x[', covariate_dimension,'])')
      }
    }
    if(is.null(basis_function_string)){
      basis_function_string <- '1'
    }else{
      basis_function_string <- gsub("^\\*(.+)", "\\1", basis_function_string)
    }
    
    return(basis_function_string)
  }
  
  if(basis_type == 'cosine'){
    for(basis_index in 1:length(index.row.prod)){
      instruction_vector <- index.matrix[basis_index,]
      name_vector[basis_index] <- spell_out_one_cosine_basis(instruction_vector)
    }
  }
  
  names(sieve.model$best_model$beta.f) <- name_vector[1:length(sieve.model$best_model$beta.f)]
  return(sieve.model)
}

######things i need to do
####prevent float number overflow
####able to process the data several times 
####other kinds of loss functions
####a better summary report. better documentation of the code.
####dynamically generating extra index matrix 
