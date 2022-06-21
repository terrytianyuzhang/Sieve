#' Preprocess the original data for sieve estimation.
#'
#' Generate the design matrix for the downstream lasso-type penalized model fitting. 
#' 
#' @param X a data frame containing original features. The (i,j)-th element is the j-th dimension of the i-th sample's feature vector. 
#' So the number of rows equals to the sample size and the number of columns equals to the feature dimension.
#' @param basisN number of sieve basis function. It is in general larger than the dimension of the original feature. 
#' Default is 50*dimension of original feature. A larger value has a smaller approximation error but it is harder to estimate.
#' The computational time/memory requirement should scale linearly to \code{basisN}.
#' @param maxj a number. the maximum index product of the basis function. A larger value means more basisN. 
#' If basisN is already specified, do not need to provide value for this argument.
#' @param type a string. It specifies what kind of basis functions are used. The default is (aperiodic) cosine basis functions, which is suitable for most purpose.
#' @param interaction_order a number. It also controls the model complexity. 1 means fitting an additive model, 2 means fitting a model allows, 3 means interaction terms between 3 dimensions of the feature, etc. The default is 3. 
#' For large sample size, lower dimension problems, try a larger value (but need to be smaller than the dimension of original features); for smaller sample size and higher dimensional problems, try set it to a smaller value (1 or 2).
#' @param index_matrix a matrix. provide a pre-generated index matrix. The default is NULL, meaning sieve_preprocess will generate one for the user.  
#' @param norm_feature a logical variable. Default is TRUE. It means sieve_preprocess will rescale the each dimension of features to 0 and 1. Only set to FALSE when user already manually rescale them between 0 and 1.
#' @param norm_para a matrix. It specifies how the features are normalized. For training data, use the default value NULL.
#'
#' @return A list containing the necessary information for next step model fitting. Typically, the list is used as the main input of Sieve::sieve_solver.
#' \item{Phi}{a matrix. This is the design matrix directly used by the next step model fitting. The (i,j)-th element of this matrix is the evaluation of i-th sample's feature at the j-th basis function. The dimension of this matrix is sample size x basisN.} 
#' \item{X}{a matrix. This is the rescaled original feature/predictor matrix.} 
#' \item{type}{a string. The type of basis funtion.}
#' \item{index_matrix}{a matrix. It specifies what are the product basis functions used when constructing the design matrix Phi. It has a dimension basisN x dimension of original features. There are at most interaction_order many non-1 elements in each row.}
#' \item{basisN}{a number. Number of sieve basis functions.}
#' \item{norm_para}{a matrix. It records how each dimension of the feature/predictor is rescaled, which is useful when rescaling the testing sample's predictors.}
#' @examples 
#' xdim <- 1 #1 dimensional feature
#' #generate 1000 training samples
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim)
#' #use 50 cosine basis functions
#' type <- 'cosine'
#' basisN <- 50 
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type)
#' #sieve.model$Phi #Phi is the design matrix
#' 
#' xdim <- 5 #1 dimensional feature
#' #generate 1000 training samples
#' #only the first two dimensions are truly associated with the outcome
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim, 
#'                               frho = 'additive', frho.para = 2)
#'                               
#' #use 1000 basis functions
#' #each of them is a product of univariate cosine functions.
#' type <- 'cosine'
#' basisN <- 1000 
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type)
#' #sieve.model$Phi #Phi is the design matrix
#' 
#' #fit a nonaprametric additive model by setting interaction_order = 1
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type, 
#'                                 interaction_order = 1)
#' #sieve.model$index_matrix #for each row, there is at most one entry >= 2. 
#' #this means there are no basis functions varying in more than 2-dimensions 
#' #that is, we are fitting additive models without interaction between features.
#' @export
#'
sieve_preprocess <- function(X, basisN = NULL, maxj = NULL, 
                             type = 'cosine', 
                             interaction_order = 3, index_matrix = NULL,
                             norm_feature = TRUE, norm_para = NULL){
  
  if(is.null(basisN)){
    warning('user did not specify number of basis functions to use, default is xdim*50')
    warning('a theoretically good number of basis functions should be around (sample size)^(1/3)*(feature dimension)^3')
    basisN <- NCOL(X)*50
  } else if(basisN == 1){
    stop('need to provide a larger basisN, you can start with xdim*50')
  }
  
  if(is(X,'numeric') | is(X, 'factor')){
    s.size <- length(X) #this is a univariate regression problem, sometimes the input data.frame is automatically conveted to an array
  }
  else{
    s.size <- dim(X)[1]
  }
  
  
  X <- sapply(X,unclass) #convert categorical variables into numeric variables
  X <- matrix(X, nrow = s.size)
  
  if(norm_feature){
    norm_list <- normalize_X(X, norm_para)
    X <- norm_list$X
    norm_para <- norm_list$norm_para
  }
  xdim <- dim(X)[2]
  
  if(is.null(index_matrix)){
  index_matrix <- create_index_matrix(xdim = xdim,
                                          basisN = basisN,
                                          maxj = maxj,
                                          interaction_order = interaction_order)
  index_matrix <- as.matrix(index_matrix[,-1]) #the first column in the product of indices of the corresponding row.
  }
  
  Phi <- Design_M_C(X, basisN, type, index_matrix)
  
  return(list(Phi = Phi, X = X, type = type, index_matrix = index_matrix,
              basisN = basisN, norm_para = norm_para))
}

#' Calculate the coefficients for the basis functions
#'
#' This is the main function that performs sieve estimation. It calculate the coefficients by solving a penalized lasso type problem. 
#' 
#' @param model a list. Typically, it is the output of Sieve::sieve_preprocess.
#' @param Y a vector. The outcome variable. The length of Y equals to the training sample size, which should also match the row number of X in model.
#' @param l1 a logical variable. TRUE means calculating the coefficients by sovling a l1-penalized empirical risk minimization problem. FALSE means solving a least-square problem. Default is TRUE.
#' @param family a string. 'gaussian', mean-squared-error regression problem.
#' @param lambda same as the lambda of glmnet::glmnet.
#' @param nlambda a number. Number of penalization hyperparameter used when solving the lasso-type problem. Default is 100.
#'
#' @return a list. In addition to the preprocessing information, it also has the fitted value.
#' \item{Phi}{a matrix. This is the design matrix directly used by the next step model fitting. The (i,j)-th element of this matrix is the evaluation of i-th sample's feature at the j-th basis function. The dimension of this matrix is sample size x basisN.}
#' \item{X}{a matrix. This is the rescaled original feature/predictor matrix. }
#' \item{beta_hat}{a matrix. Dimension is basisN x nlambda. The j-th column corresponds to the fitted regression coeffcients using the j-th hyperparameter in lambda.}
#' \item{type}{a string. The type of basis funtion.}
#' \item{index_matrix}{a matrix. It specifies what are the product basis functions used when constructing the design matrix Phi. It has a dimension basisN x dimension of original features. There are at most interaction_order many non-1 elements in each row.}
#' \item{basisN}{a number. Number of sieve basis functions.}
#' \item{norm_para}{a matrix. It records how each dimension of the feature/predictor is rescaled, which is useful when rescaling the testing sample's predictors.}
#' \item{lambda}{a vector. It records the penalization hyperparameter used when solving the lasso problems. Default has a length of 100, meaning the algorithm tried 100 different penalization hyperparameters.}
#' \item{family}{a string. 'gaussian', continuous numerical outcome, regression probelm; 'binomial', binary outcome, classification problem.}
#' @examples 
#' xdim <- 1 #1 dimensional feature
#' #generate 1000 training samples
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim)
#' #use 50 cosine basis functions
#' type <- 'cosine'
#' basisN <- 50 
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type)
#' sieve.fit<- sieve_solver(model = sieve.model, Y = TrainData$Y)
#' 
#' ###if the outcome is binary, 
#' ###need to solve a nonparametric logistic regression problem
#' xdim <- 1
#' TrainData <- GenSamples(s.size = 1e3, xdim = xdim, y.type = 'binary', frho = 'nonlinear_binary')
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type)
#' sieve.fit<- sieve_solver(model = sieve.model, Y = TrainData$Y,
#'                          family = 'binomial')
#' @export
#'
sieve_solver <- function(model, Y, l1 = TRUE, family = "gaussian", 
                         lambda = NULL, nlambda = 100){
  if(l1 == FALSE){
    ## this is least square sieve estimator
    Phi <- model$Phi
    beta_hat <- least_square_C(Phi, Y)
    
    model$beta_hat <- beta_hat
    return(model)
  }else{
    #remove the intercept column in the design matrix
    #need to allow it
    mo <- glmnet::glmnet(model$Phi[,-1], Y, family = family, alpha = 1, 
                         nlambda = nlambda, intercept = TRUE, standardize = FALSE,
                         lambda = lambda)
    
    model$lambda <- mo$lambda
    model$beta_hat <- coef(mo)
    model$family <- family
    
    return(model)
  }
}


create_index_matrix <- function(xdim, basisN = NULL, maxj = NULL, interaction_order = 5){
  #xdim: the dimension of feature x
  #dimlimit: the working dimension, same as D' in the paper
  #basisN: default is NULL. We use this to specify the total number of rows in the index list.
  #maxj: default is NULL. We use this to specify the largest row product in the index list.
  
  index_matrix <- matrix(1, nrow = 1, ncol = xdim) #first row is all 1
  interaction_order <- min(interaction_order, xdim) #this determines the maximal number of non-1 elements each row
  
  if(!is.null(maxj)){#when we specify the max product value of the indices
    
    for(product_v in 2:maxj){#loop through all the product values 
      # product_v <- as.integer(2)
      
      # use C function to generate all possible factorization of product_v
      # each factorization uses at most interaction_order positive integers
      factors_v <- Generate_factors(product_v, interaction_order) 
      
      # this list store all greater than 1 parts of each factorization
      # it is order- sensitive (2,3) != (3,2)
      greater_than_1 <- list() 
      if(length(factors_v) > 0){
        for(i in 1:length(factors_v)){
          greater_than_1 <- c(greater_than_1, unique(combinat::permn(factors_v[[i]]))) #for each factorization, we find all its unique permutations
        }
      }
      
      greater_than_1 <- c(greater_than_1, c(product_v)) #add the trivival 1*1...1*product_v factorization
      
      for(i in 1:length(greater_than_1)){
        m <- length(greater_than_1[[i]]) #how many non-one position does this permutation need
        positions <- t(combn(1:xdim,m))
        new_index <- matrix(1, nrow = dim(positions)[1], ncol = xdim)
        
        for(j in 1:(dim(positions)[1])){
          new_index[j, positions[j,]] <- greater_than_1[[i]]
        }
        
        index_matrix <- rbind(index_matrix, new_index)
      }#loop in greater_than_1
    }#loop in product_v
    
  }else if(!is.null(basisN)){
    product_v <- 2
    while(dim(index_matrix)[1] < basisN){
      # product_v <- as.integer(2)
      factors_v <- Generate_factors(product_v, interaction_order)
      greater_than_1 <- list()
      if(length(factors_v) > 0){
        for(i in 1:length(factors_v)){
          greater_than_1 <- c(greater_than_1, unique(combinat::permn(factors_v[[i]])))
        }
      }
      greater_than_1 <- c(greater_than_1, c(product_v))
      greater_than_1
      
      for(i in 1:length(greater_than_1)){
        m <- length(greater_than_1[[i]]) #how many non-one position does this permutation need
        positions <- t(combn(1:xdim,m))
        new_index <- matrix(1, nrow = dim(positions)[1], ncol = xdim)
        for(j in 1:(dim(positions)[1])){
          new_index[j, positions[j,]] <- greater_than_1[[i]]
        }
        index_matrix <- rbind(index_matrix, new_index)
      }
      product_v <- product_v + 1
    }
  }
  index_matrix <- cbind(rep(0, dim(index_matrix)[1]),
                      index_matrix)
  
  for(j in 1:(dim(index_matrix)[1])){
    index_matrix[j,1] <- prod(index_matrix[j,-1]) #the first column is the row product
  }
  
  return(index_matrix)
}

normalize_X <- function(X, norm_para  = NULL, lower_q = 0.025, upper_q = 0.975){
  #normalize the support of covariate so that they are between 0,1.
  #
  
  if(is.null(norm_para)){
    norm_para <- matrix(0, nrow = 2, ncol = dim(X)[2])
    for(i in 1:(dim(X)[2])){
      lower_val <- quantile(X[,i],lower_q)
      upper_val <- quantile(X[,i],upper_q)
      X[,i] <- (X[,i] - lower_val)/(upper_val - lower_val)
      norm_para[1,i] <- lower_val
      norm_para[2,i] <- upper_val
    }
    X[X<=0] <- 0
    X[X>=1] <- 1
  }else{
    for(i in 1:(dim(X)[2])){
      lower_val <- norm_para[1,i]
      upper_val <- norm_para[2,i]
      X[,i] <- (X[,i] - lower_val)/(upper_val - lower_val)
    }
    X[X<=0] <- 0
    X[X>=1] <- 1
  }
  return(list(X = X, norm_para = norm_para))
}

#' Predict the outcome of interest for new samples
#'
#' Use the fitted sieve regression model from sieve_solver. It also returns the testing mean-squared errors.
#' 
#' @param model a list. Use the fitted model from sieve_solver.
#' @param testX a data frame. Dimension equals to test sample size x feature diemnsion. Should be of a similar format as the training feature provided to sieve_preprocess.
#' @param testY a vector. The outcome of testing samples (if known). Default is NULL. For regression problems, the algorithm also returns the testing mean-squared errors.
#'
#' @return a list. 
#' \item{predictY}{a matrix. Dimension is test sample size (# of rows) x number of penalty hyperparameter lambda (# of columns). 
#' For regression problem, that is, when family = "gaussian", each entry is the estimated conditional mean (or predictor of outcome Y). For classification problems (family = "binomial"), each entry is the predicted probability of having Y = 1 (which class is defined as "class 1" depends on the training data labeling). }
#' \item{MSE}{For regression problem, when testY is provided, the algorithm also calculates the mean-sqaured errors using testing data. Each entry of \code{MSE} correponds to one value of penalization hyperparameter \code{lambda}}
#' @examples 
#' xdim <- 1 #1 dimensional feature
#' #generate 1000 training samples
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim)
#' #use 50 cosine basis functions
#' type <- 'cosine'
#' basisN <- 50 
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type)
#' sieve.fit<- sieve_solver(model = sieve.model, Y = TrainData$Y)
#' #generate 1000 testing samples
#' TestData <- GenSamples(s.size = 1000, xdim = xdim)
#' sieve.prediction <- sieve_predict(model = sieve.fit, 
#'                                   testX = TestData[,2:(xdim+1)], 
#'                                   testY = TestData$Y)
#' ###if the outcome is binary, 
#' ###need to solve a nonparametric logistic regression problem
#' xdim <- 1
#' TrainData <- GenSamples(s.size = 1e3, xdim = xdim, y.type = 'binary', frho = 'nonlinear_binary')
#' sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
#'                                 basisN = basisN, type = type)
#' sieve.fit<- sieve_solver(model = sieve.model, Y = TrainData$Y,
#'                          family = 'binomial')
#'                          
#' ###the predicted value is conditional probability (of taking class 1).
#' TrainData <- GenSamples(s.size = 1e3, xdim = xdim, y.type = 'binary', frho = 'nonlinear_binary')
#' sieve.prediction <- sieve_predict(model = sieve.fit, 
#'                                   testX = TestData[,2:(xdim+1)])
#' @export
#'
sieve_predict <- function(model, testX, testY = NULL){
  
  if(!is.null(model$lambda)){
    lambda_num <- length(model$lambda)
    beta_hat <- matrix(model$beta_hat, ncol = lambda_num)
  }else{
    lambda_num <- 1
    beta_hat <- model$beta_hat
  }
  
  type <- model$type
  basisN <- model$basisN
  index_matrix <- model$index_matrix
  norm_para <- model$norm_para #specifies how to normalize the testing features/predictors.
  family <- model$family
  
  if(is(testX, 'numeric') | is(testX, 'factor')){
    test.size <- length(testX)
    xdim <- 1
  }else{
    test.size <- dim(testX)[1]
    xdim <- dim(testX)[2]
  }


  # testX <- as.matrix(testX, nrow = test.size)
  
  
  test_list <- sieve_preprocess(X = testX, basisN = basisN, 
                                       type = type, index_matrix = index_matrix, 
                                       norm_para = norm_para) #provide index_matrix

  predictY <- matrix(0, nrow = test.size, ncol = lambda_num)
  test_Phi <- test_list$Phi
  
  if(family == 'gaussian'){
    #predictY is the estimated conditional mean
    if(!is.null(testY)){
      MSE <- rep(0, lambda_num)
      for(i in 1:lambda_num){
        predictY[,i] <- crossprod_C(test_Phi, matrix(beta_hat[,i]))
        MSE[i] <- mean((predictY[,i] - testY)^2)
      }
      return(list(predictY = predictY, MSE = MSE))
    }else{
      for(i in 1:lambda_num){
        predictY[,i] <- crossprod_C(test_Phi, matrix(beta_hat[,i]))
      }
      return(list(predictY = predictY))
    }
  }else if(family == 'binomial'){
    #we first calculate the log odss function, then transform it back to probability
    for(i in 1:lambda_num){
      predictY[,i] <- crossprod_C(test_Phi, matrix(beta_hat[,i]))
    }
    predictY <- exp(predictY)/(1 + exp(predictY))
    return(list(predictY = predictY))
  }
  
}

#' Generate some simulation/testing samples with nonlinear truth.
#'
#' This function is used in several examples in the package.
#' 
#' @param s.size a number. Sample size.
#' @param xdim a number. Dimension of the feature vectors X.
#' @param x.dis a string. It specifies the distribution of feature X. The default is uniform distribution over \code{xdim}-dimensional unit cube.
#' @param x.para extra parameter to specify the feature distribution.
#' @param frho a string. It specifies the true regression/log odds functions used to generate the data set. The default is a linear function. 
#' @param frho.para extra parameter to specify the true underlying regression/log odds function.
#' @param y.type a string. Default is \code{y.type = 'continuous'}, meaning the outcome is numerical and the problem is regression. Set it to \code{y.type = 'binary'} for binary outcome.
#' @param noise.dis a string. For the distribution of the noise variable (under regression probelm settings). Default is Gaussian distribution.
#' @param noise.para a number. It specifies the magnitude of the noise in regression settings.
#' 
#' @return a \code{data.frame}. The variable \code{Y} is the outcome (either continuous or binary). Each of the rest of the variables corresponds to one dimension of the feature vector.
#' @examples 
#' xdim <- 1 #1 dimensional feature
#' #generate 1000 training samples
#' TrainData <- GenSamples(s.size = 1000, xdim = xdim)
#' #generate some noise-free testing samples
#' TestData <- GenSamples(s.size = 1000, xdim = xdim, noise.para = 0)
#' @export
#'

GenSamples <- function(s.size, xdim = 1, x.dis = "uniform",
                     x.para = NULL, frho = 'linear', frho.para = 1e2,
                     y.type = 'continuous', noise.dis = 'normal',
                     noise.para = 0.5){
  
  if(x.dis == 'uniform'){
    x <- matrix(runif(s.size * xdim, min = 0, max = 1), ncol = xdim)
  }else if(x.dis == 'dep12'){
    x <- matrix(runif(s.size * xdim, min = 0, max = 1), ncol = xdim)
    tempm <- matrix(0, nrow = s.size, ncol = 2)
    tempm[,1] <- x[,2]
    tempm[,2] <- x[,1]
    x[,1:2] <- (x[,1:1] + 0.5*tempm)/2
  }
  
  if(y.type == 'continuous'){
    if(noise.dis == 'normal'){
      if(frho == 'multiellipsoid'){
        basisN <- frho.para
        index_matrix <- as.matrix(create_index_matrix(xdim, basisN = basisN, 
                                                  interaction_order = xdim)[,-1], 
                                ncol = xdim)
        y <- apply(x, 1, truef, frho, index_matrix) + 
          rnorm(s.size, 0, noise.para)
      }else if(frho == 'STE'){
        basisN <- frho.para$basisN
        D <- frho.para$D
        index_matrix <- as.matrix(create_index_matrix(D, basisN = basisN, #the dimension of the index_matrix is specified by D (not xdim)
                                                  interaction_order = D)[1:basisN,-1], 
                                ncol = D)
        frho.para$index_matrix <- index_matrix
        y <- apply(x[,1:D], 1, truef, frho, frho.para) + #only the first D columns of x is relevant
          rnorm(s.size, 0, noise.para)
        }else{
        y <- apply(x, 1, truef, frho, frho.para) + 
          rnorm(s.size, 0, noise.para)
        }
     
    }
  }else if(y.type == 'binary'){
    y <- apply(x, 1, truef, frho)
  }
  
  traindata <- data.frame(cbind(y, x))
  
  #rename the columns
  names(traindata)[1] <- 'Y'
  for(i in 2:(xdim+1)){
    names(traindata)[i] <- paste0('X',i-1)
  }
  return(traindata)
}

truef <- function(x, FUN = 'linear', para = NULL){
  y <- 0
  xdim <- length(x)
  if(FUN == 'linear'){
    for(i in 1:xdim){
      y <- y + x[i]
    }
  }else if(FUN == 'additive'){
    D <- para
    for(i in 1:D){
      y <- y + as.numeric(i %% 2 == 1)*(0.5 - abs(x[i]-0.5)) + as.numeric(i %% 2 == 0)*exp(-x[i])
      # y <- y+ (x[i] * x[i])^2 + (x[i] * x[i+1])^3 + sin(3*x[i])*x[i+2] + x[i]^2*sin(x[i+2])
    }
  }else if(FUN == 'interaction'){
    D <- para
    for(i in 1:(D-1)){
      y <- y + psi(x[i], 2, 'legendre')*psi(x[i+1], 3, 'legendre')
      # y <- y+ (x[i] * x[i])^2 + (x[i] * x[i+1])^3 + sin(3*x[i])*x[i+2] + x[i]^2*sin(x[i+2])
    }
  }
  else if(FUN == 'sinsin'){ #this truth will ruin the additive model
      for(i in 1:(xdim-1)){
        y <- y + sin(2*pi*x[i]) * sin(2*pi*x[i+1])
        }
      }else if(FUN == 'sparse'){
        # y <- y+x[1]
        y <- y + 1
      # y <- y + cos(2*pi*x[1])*cos(2*pi*x[2])
        # y <- y + x[1]*x[2]
        # y <- y+ (6*x[1]-3)*sin(12*x[1]-6) *(6*x[2]-3)*sin(12*x[2]-6)
      # y <- y + (x[3]-0.5)*x[4]
      # this is the one in sparse additive model (down)
      y <- y  - sin(1.5*x[1])*(x[2])^3+ 1.5*(x[2]-0.5)^2
      y <- y+ sin(x[1])*cos(x[2])*((x[3])^3+ 1.5*(x[3]-0.5)^2) * (sin(exp(-0.5*x[4])))
      y <- y+ sin(x[2])*cos(x[3])*((x[4])^3+ 1.5*(x[4]-0.5)^2) * (sin(exp(-0.5*x[1])))
      }else if(FUN == 'lineartensor'){
        D <- para
    for(i in 1:(D-1)){
      # y <- y + psi(x[i], 3, 'cosine') + psi(x[i], 2, 'cosine')*psi(x[i+1], 2, 'cosine')
      y <- y + psi(x[i], 3, 'legendre') + psi(x[i], 2, 'legendre')*psi(x[i+1], 2, 'legendre')
    }
  }else if(FUN == 'turntensor'){
    for(i in 1:xdim){
      for(j in i:xdim){
        y <- y+ (0.5 - abs(x[i]-0.5)) * (0.5 - abs(x[j]-0.5))
      }
    }
  }else if(FUN == 'multiellipsoid'){
      y <- 0
      index_matrix <- para
      for(i in 1:(dim(index_matrix)[1])){ #index_matrix is created in the Gen_Train function
        # y <- y + (prod(index_matrix[i,]))^(-1.5) * multi_psi(x, index_matrix[i,], 'legendre')
        y <- y + (prod(index_matrix[i,]))^(-1.5) * multi_psi(x, index_matrix[i,], 'sobolev1')
        # y <- y + 3*(prod(index_matrix[i,]))^(-1.5) * multi_psi(max(x-0.5,0), index_matrix[i,], 'sobolev1')
        }
  }else if(FUN == 'STE'){
      y <- 0
      # y <- y+ (0.5 - abs(x[1]-0.5)) * (0.5 - abs(x[2]-0.5))
      index_matrix <- para$index_matrix
      for(i in 1:para$basisN){ #index_matrix is created in the Gen_Train function
        # y <- y + (prod(index_matrix[i,]))^(-1.6) * multi_psi(x, index_matrix[i,], 'legendre')
        if(prod(index_matrix[i,]) <= 8){
          # y <- y + multi_psi(x, index_matrix[i,], 'cosine')
          # print(index_matrix[i,])
          y <- y + 1* multi_psi(x, index_matrix[i,], 'cosine')
          }else{
            #######!!!!!!!!!!!!!!!!###########
          y <- y +0*para$coefs[i] * (prod(index_matrix[i,]))^(-1.6) * multi_psi(x, index_matrix[i,], 'cosine')
          # y <- y + para$coefs[i] * (prod(index_matrix[i,]))^(-1.6) * multi_psi(x, index_matrix[i,], 'cosine')
          }
        # y <- y + 3*(prod(index_matrix[i,]))^(-1.5) * multi_psi(max(x-0.5,0), index_matrix[i,], 'sobolev1')
      }
  }else if(FUN == 'sinproduct'){
    y1 <- y2 <- 1
    for(i in 1:xdim){
      y1 <- y1 * cos(3*x[i])
      y2 <- y2 * cos(6*x[i])
    }
    y <- y1 + y2
  }else if(FUN == 'linear_binary'){
    xdim <- length(x)
    y <- rbinom(1, 1, sum(x)/xdim)
  }else if(FUN == 'nonlinear_binary'){
    xdim <- length(x)
    y <- rbinom(1, 1, sum(abs(x-0.5))/xdim+0.2)
    # y <- rbinom(1,1, as.numeric(x[1] < 0.5))
  }
  
  
  return(y)
}

# my.legendre <- function(x, j){
#   x <- (x-0.5)*2
#   if(j == 1){
#     y <- 1
#   }
#   else if(j == 2){
#     y <- x
#   }
#   else if(j == 3){
#     y <- (3*x^2 - 1)/2
#   }
#   else if(j == 4){
#     y <- (5*x^3 - 3*x)/2
#   }
#   else if(j == 5){
#     y <- (35*x^4 - 30*x^2 + 3)/8
#   }
#   else if(j == 6){
#     y <- (63*x^5 - 70*x^3 + 15*x)/8
#   }
#   else if(j == 7){
#     y <- (231*x^6 - 315*x^4 + 105*x^2 - 5)/16
#   }
#   else if(j == 8){
#     y <- (429*x^7 - 693*x^5 + 315*x^3 - 35*x)/16
#   }
#   else if(j == 9){
#     y <- (6435*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35)/128
#   }
#   else{
#     print(paste0("trying to calculate ",j,"-th order Legendre polynomial"))
#     print(paste0("formula not provided"))
#     return(NA)
#   }
#   return(y * sqrt((2*j+1)/2)) #normalization step
# }

# psi <- function(x, j, type){
#   #this is univariate basis function
#   if(type == 'legendre'){
#     return(my.legendre(x, j))
#   }else if(type == 'sobolev1'){
#     if(j == 1){
#       return(1)
#     }else{
#       return(sin((2*(j-1) - 1)*pi*x/2))
#       }
#   }else if(type == 'polytri'){
#     if(j == 1){
#       return(1)
#     }
#     else if(j == 2){
#       return(x)
#     }
#     else if(j %% 2 == 1){
#       return(sin((j - 1)/2 *pi* x))
#     }
#     else{
#       return(cos((j - 2)/2 *pi* x))  
#       }
#     }
# }




# x <- c(-0.5,0.3,0.4)
# z <- c(-0.5,0.2,0.4)
# type <- 'sobolev1'
# tensor_kernel(x,z,type)


# x <- c(-0.5,0.2,0.4)
# type <- 'legendre'
# index <- c(1,1,5)
# multi_psi(x,index,type)

# multi_psi <- function(x, index, type){
#   xdim <- length(x)
#   psix <- 1
#   for(i in 1:xdim){
#     psix <- psix * psi(x[i], index[i], type)
#   }
#   return(psix)
# }

tensor_kernel <- function(x, z, type){
  my.kernel <- function(x, z, type){
    if(type == 'sobolev1'){
      return(1+min(x,z))
    }else if(type == 'cubicspline'){
      return(1 + x*z + as.numeric(x <= z)*(z*x^2/2 - x^3/6) + as.numeric(x > z)*(x*z^2/2 - z^3/6))
      }else if(type =='gaussian'){
      return(exp(-(x-z)^2))
      }
  }
  
  xdim <- length(x)
  pik <- 1
  for(i in 1:xdim){
    pik <- pik * my.kernel(x[i], z[i], type)
  }
  return(pik)
}

KRR_preprocess <- function(X, type, kernel.para = -1){
  KnSVD <- Kernel_M_C(X, type, kernel.para) #call the C function to do svd
  lambdas <- exp(seq(log(1e-4*min(KnSVD$s)),
                     log(max(KnSVD$s)),
                     length = 10))
  return(list(K = KnSVD$K, U = KnSVD$U, s = KnSVD$s,
              lambdas = lambdas, X = X, type = type,
              kernel.para = kernel.para))
}


# KRR.fit2(model_pre, Y = TrainData$Y, lambda)
# U <- model_pre$U
# s <- model_pre$s
# lambda <- 0.01
# Y <- TrainData$Y
# KRR_cal_beta_C(U, s, lambda, Y)
KRR.fit_C <- function(model_pre, Y, lambda){
  U <- model_pre$U
  s <- model_pre$s
  beta <- KRR_cal_beta_C(U, s, lambda, Y)
  return(list(beta.hat = beta, X = model_pre$X, type = model_pre$type,
              kernel.para = model_pre$kernel.para, lambda = lambda))
}

KRR.fit <- function(X, Y, lambda, type){
  s.size <- dim(X)[1]
  X <- as.matrix(X, nrow = s.size) #accelerate the matrix filling
  # K <- matrix(0, nrow = s.size, ncol = s.size)
  # for(i in 1:s.size){
  #   for(j in 1:s.size){
  #     K[i,j] <- tensor_kernel(X[i,],X[j,],type)
  #   }
  # }
  K <- outer( 
    1:s.size, 1:s.size, 
    Vectorize( function(i,j) tensor_kernel(X[i,], X[j,], type = type) ) 
  )
  beta.hat <- solve(K + lambda * diag(1,nrow = s.size)) %*%  Y 
  return( list(beta.hat = beta.hat, K = K, X = X, Y = Y, type = type))
}

# testX <- matrix(c(0.5,0.2,0.4, 0.5,0.3,0.4), nrow = 2)
# testY <- c(2,3)
# model <- KRR.fit(X,Y,lambda, type)
# KRR.predict <- function(testX, testY, model){
#   beta.hat <- model$beta.hat
#   trainX <- model$X
#   type <- model$type
#   kernel.para <- model$kernel.para
#   
#   trainX <- as.matrix(trainX, nrow = train.size)
#   testX <- as.matrix(testX, nrow = test.size)
#   
#   predictY <- KRR_predict_C(trainX, testX, type, beta.hat, kernel.para)
#   # train.size <- dim(trainX)[1]
#   # test.size <- dim(testX)[1]
#   # K <- matrix(0, nrow = test.size, ncol = train.size)
#   # for(i in 1:test.size){
#   #   for(j in 1:train.size){
#   #     K[i,j] <- tensor_kernel(testX[i,],trainX[j,],type)
#   #   }
#   # }
#   # predictY <- K %*% beta.hat
#   
#   MSE <- mean((predictY - testY)^2)
#   return(list(predictY = predictY, MSE = MSE))
# }
# KRR.predict(testX, testY, model)

# X <- matrix(1:9/12, nrow = 3)
# Y <- c(1,2,3)
# type <- 'legendre'
# basisN <- 2
# model <- sieve.fit(X,Y,basisN, type)
# 
# X = TrainData[,2:(xdim+1)]
# Y = TrainData$Y
# basisN = 3
# type = 'sobolev1'
# model <- sieve.fit(X,Y,basisN, type)



sieve.fit <- function(X, Y, basisN, type){
  s.size <- dim(X)[1]
  X <- as.matrix(X, nrow = s.size)
  xdim <- dim(X)[2]
  
  index_matrix <- as.matrix(create_index_matrix(xdim, basisN = basisN)[,-1], ncol = xdim)
  
  K <- matrix(0, nrow = s.size, ncol = basisN)
  for(i in 1:s.size){
    for(j in 1:basisN){
      K[i,j] <- multi_psi(X[i,], index_matrix[j,], type)
    }
  }
  beta.hat <- solve(crossprod(K)) %*% t(K) %*% Y 
  return( list(beta.hat = beta.hat, K = K, X = X, Y = Y,
               type = type, basisN = basisN, index_matrix = index_matrix))
}


# testX <- matrix(c(0.5,0.2,0.4, 0.5,0.3,0.4), nrow = 2)
# testY <- c(2,3)
# sieve.predict(testX, testY, model)
sieve.predict <- function(testX, testY, model){
  beta.hat <- model$beta
  trainX <- model$X
  type <- model$type
  basisN <- model$basisN

  train.size <- dim(trainX)[1]
  test.size <- dim(testX)[1]
  trainX <- as.matrix(trainX, nrow = train.size)
  testX <- as.matrix(testX, nrow = test.size)
  xdim <- dim(trainX)[2]

  index_matrix <- as.matrix(create_index_matrix(xdim, basisN = basisN)[,-1], ncol = xdim)
  
  K <- matrix(0, nrow = test.size, ncol = basisN)
  for(i in 1:test.size){
    for(j in 1:basisN){
      K[i,j] <- multi_psi(testX[i,], index_matrix[j,], type)
    }
  }
  predictY <- K %*% beta.hat
  MSE <- mean((predictY - testY)^2)
  return(list(predictY = predictY, MSE = MSE))
}

# create_index_matrix_iso <- function(basisN, xdim){
#   lista <- matrix(rep(1,xdim), ncol = xdim)
#   listb <- all_add_one(lista)
#   lista <- data.table(cbind(1,lista))
#   listb <- data.table(cbind(rep(2,xdim), listb))
#   while(dim(lista)[1] < basisN){
#     temp <- listb[V1 == min(V1),]
#     lista <- rbind(lista, temp)
#     for(i in 1:(dim(temp)[1])){
#       tempi <- cbind(rep(-1,xdim), all_add_one(temp[i,2:(xdim+1)]))
#       for(j in 1:(dim(tempi)[1])){
#         tempi[j,1] <- norm(as.array(tempi[j,-1]), type = '2')
#         # tempi[j,1] <- prod(tempi[j,-1])
#         }
#       listb <- rbind(listb, tempi)
#     }
#     listb <- unique(listb[V1 > min(V1),])
#   }
#   return(lista) #the first column indicates the ordering
# }

all_add_one <- function(index){
  index <- as.matrix(index)
  xdim <- length(index)
  added <- matrix(rep(index, xdim), ncol = xdim, byrow = TRUE)
  for(i in 1:xdim){
    added[i,i] <- added[i,i] + 1
  }
  return(added)
}
