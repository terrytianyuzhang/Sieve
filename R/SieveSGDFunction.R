#' Preprocess the original data for sieve estimation.
#'
#' Convert the covariates into numerical variables and normalize them.
#' 
#' @param X a data frame containing prediction features/ independent variables. The (i,j)-th element is the j-th dimension of the i-th sample's feature vector. 
#' So the number of rows equals to the sample size and the number of columns equals to the feature/covariate dimension. If the complete data set is large, this can be a representative subset of it (ideally have more than 1000 samples). 
#' @param s numerical array. Smoothness parameter. Default is 2. The elements of this array should take values greater than 0.5. The larger s is, the smoother we are assuming the truth to be.
#' @param r0 numerical array. Initial learning rate/step size. The step size at each iteration will be r0*(sample size)^(-1/(2s+1)), which is slowly decaying.
#' @param J numerical array. Initial number of basis functions. The number of basis functions at each iteration will be J*(sample size)^(1/(2s+1)), which is slowly increasing.
#' @param type a string. It specifies what kind of basis functions are used. The default is (aperiodic) cosine basis functions, which is suitable for most purpose.
#' @param interaction_order a number. It also controls the model complexity. 1 means fitting an additive model, 2 means fitting a model allows, 3 means interaction terms between 3 dimensions of the feature, etc. The default is 3. 
#' For large sample size, lower dimension problems, try a larger value (but need to be smaller than the dimension of original features); for smaller sample size and higher dimensional problems, try set it to a smaller value (1 or 2).
#' @param norm_feature a logical variable. Default is TRUE. It means sieve_preprocess will rescale the each dimension of features to 0 and 1. Only set to FALSE when user already manually rescale them between 0 and 1.
#' @param norm_para a matrix. It specifies how the features are normalized. For training data, use the default value NULL.
#' @param lower_q lower quantile used in normalization. Default is 0.01 (1% quantile).
#' @param upper_q upper quantile used in normalization. Default is 0.99 (99% quantile).
#'
#' @return A list containing the necessary information for next step model fitting. Typically, the list is used as the main input of Sieve::sieve_solver.
#' \item{type}{a string. The type of basis funtion.}
#' \item{interaction_order}{the highest order of interaction.}
#' \item{norm_para}{a matrix. It records how each dimension of the feature/predictor is rescaled, which is useful when rescaling the testing sample's predictors.}
#' \item{X}{a matrix. This is the rescaled original covariate matrix.} 
#' @examples 
#' xdim <- 1 #1 dimensional feature
#' #generate 1000 training samples
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim)
#' #use 50 cosine basis functions
#' type <- 'cosine'
#' sieve.model <- sieve.sgd.preprocess(X = TrainData[,2:(xdim+1)], type = type)
#' @export
#'
sieve.sgd.preprocess <- function(X, s = c(2), r0 = c(2), J = c(1), type = c('cosine'), 
                             interaction_order = c(3), 
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
    norm_list <- Sieve:::normalize_X(X, norm_para = norm_para, 
                                     lower_q = lower_q, upper_q = upper_q)
    # X <- norm_list$X
    norm_para <- norm_list$norm_para
  }
  
  #####set the number of observations counter back to 0
  s.size.sofar <- 0
  
  #####hyperparameter list
  hyper.para.list <- list(s = s, r0 = r0,
                          J = J, interaction_order = interaction_order,
                          omega = 0.51) #omega is related to a basis function-specific learning rate
  
  #####generate hyperparameter index matrix
  hyper.para.index.matrix <- expand.grid(1:length(s), 1:length(r0), 1:length(J),
                                         1:length(interaction_order), 1:length(omega))
  colnames(hyper.para.index.matrix) <- c('s', 'r0', 'J', 'interaction_order', 'omega')
  
  #####generate basis index matrix
  max.basisN <- ceiling( max(J) * s.size^(1/(2*min(s) + 1))
                        ) #the maximum number of basis functino needed to process all the data in X
  xdim <- dim(X)[2] #dimension of predictors
  index_matrix <- as.matrix(Sieve:::create_index_matrix(xdim, basisN = max.basisN, #the dimension of the index_matrix is specified by D (not xdim)
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

sieve.sgd.solver <- function(sieve.model, X, Y){
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
  
  norm_list <- Sieve:::normalize_X(X, norm_para = norm_para)
  X <- norm_list$X
  
  #####feed in the data one by one
  J <- sieve.model$hyper.para.list$J
  s <- sieve.model$hyper.para.list$s
  r0 <- sieve.model$hyper.para.list$r0
  omega <- sieve.model$hyper.para.list$omega
  index_matrix <- sieve.model$index.matrix
  index.row.prod <- sieve.model$index.row.prod
  type <- sieve.model$type
  M <- length(sieve.model$inf.list)
  
  for(i in 1:s.size){
    newx <- matrix(X[i,], nrow = 1)
    newy <- Y[i]
    
    max.J <- ceiling(max(J) * i^(1/(2*min(s) + 1))) #this is the maximum number of basis functions need to estimated at this step
    
    Phi <- Sieve:::Design_M_C(newx, max.J, type, index_matrix) #one row design "matrix"

    for(m in 1:M){
      
      #####calculate rolling CV
      
      ###number of basis functions for this hyper-para combination
      J.m <- J[sieve.model$inf.list[[m]]$hyper.para.index$J] * i^(1/(2* s[sieve.model$inf.list[[m]]$hyper.para.index$s]+ 1))
      J.m <- ceiling(J.m)
      
      beta.f <- sieve.model$inf.list[[m]]$beta.f
      beta.f.int <- sieve.model$inf.list[[m]]$beta.f.int
      
      if(length(beta.f) < J.m){ #check if more basis are used
        beta.f <- c(beta.f, rep(0, J.m - length(beta.f)))
        beta.f.int <- c(beta.f.int, rep(0, J.m - length(beta.f.int)))
      }
      
      fnewx <- crossprod(beta.f, Phi[1:J.m]) ###f_i(X_{i+1})
      sieve.model$inf.list[[m]]$rolling.cv <- sieve.model$inf.list[[m]]$rolling.cv + (newy - fnewx)^2
      
      ##########update beta's
      ###overall learning rate
      rn.m <- r0[sieve.model$inf.list[[m]]$hyper.para.index$r0] * i^(-1/(2* s[sieve.model$inf.list[[m]]$hyper.para.index$s]+ 1))
      
      #####CORE UPDATING STEP
      fnewx.int <- crossprod(beta.f.int, Phi[1:J.m])
      res <- as.numeric(newy - fnewx.int)
      
      beta.f.int <- beta.f.int + rn.m * res * (index.row.prod[1:J.m])^(-2*omega)*Phi[1:J.m]
      
      beta.f <- (i-1)/i * beta.f + beta.f.int/i
      ################
      
      sieve.model$inf.list[[m]]$beta.f <- beta.f
      sieve.model$inf.list[[m]]$beta.f.int <- beta.f.int
      
    }
  }
  
  ######model comparison
  rolling.cvs <- rep(1e10, M)
  for(m in 1:M){
    rolling.cvs[m] <- sieve.model$inf.list[[m]]$rolling.cv
  }
  print(rolling.cvs/s.size)
  print(which.min(rolling.cvs))
  return(sieve.model)
  
  
  
  
  
  
}
xdim <- 1
TrainData <- GenSamples(s.size = 1e4, xdim = xdim)
X <- TrainData[,2:(xdim+1)]
Y <- TrainData[,1]
sieve.model <- sieve.sgd.preprocess(X = TrainData[,2:(xdim+1)], type = type,
                                    s = c(1,2),
                                    r0 = c(1,2, 4),
                                    J = c(0.5,1,2))
sieve.model= sieve.sgd.solver(sieve.model = sieve.model, X = X, Y  = Y)
raw$inf.list[[12]]
raw$inf.list[[2]]
raw$inf.list[[4]]

TestData <- GenSamples(s.size = 1e4, xdim = xdim, noise.para = 0)
X <- TestData[,2:(xdim+1)]
Y <- TestData[,1]
sieve.model <- sieve.sgd.predict(sieve.model, TestData[,2:(xdim+1)])

mean((Y - sieve.model$inf.list[[13]]$prdy)^2)
mean((Y - sieve.model$inf.list[[14]]$prdy)^2)
mean((Y - sieve.model$inf.list[[6]]$prdy)^2)

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
  
  norm_list <- Sieve:::normalize_X(X, norm_para = norm_para)
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
    
    Phi <- Sieve:::Design_M_C(X, J.m, type, index_matrix)
    
    tmp.prdy <- Sieve:::crossprod_C(Phi, matrix(sieve.model$inf.list[[m]]$beta.f))
    
    sieve.model$inf.list[[m]]$prdy <- tmp.prdy
  }
  
  return(sieve.model)
  
}


