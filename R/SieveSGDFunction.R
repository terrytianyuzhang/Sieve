#' Preprocess the original data for sieve estimation.
#'
#' Convert the covariates into numerical variables and normalize them.
#' 
#' @param X a data frame containing original features. The (i,j)-th element is the j-th dimension of the i-th sample's feature vector. 
#' So the number of rows equals to the sample size and the number of columns equals to the feature dimension.
#' @param type a string. It specifies what kind of basis functions are used. The default is (aperiodic) cosine basis functions, which is suitable for most purpose.
#' @param interaction_order a number. It also controls the model complexity. 1 means fitting an additive model, 2 means fitting a model allows, 3 means interaction terms between 3 dimensions of the feature, etc. The default is 3. 
#' For large sample size, lower dimension problems, try a larger value (but need to be smaller than the dimension of original features); for smaller sample size and higher dimensional problems, try set it to a smaller value (1 or 2).
#' @param norm_feature a logical variable. Default is TRUE. It means sieve_preprocess will rescale the each dimension of features to 0 and 1. Only set to FALSE when user already manually rescale them between 0 and 1.
#' @param norm_para a matrix. It specifies how the features are normalized. For training data, use the default value NULL.
#' @param lower_q lower quantile used in normalization. Default is 0.01 (1% quantile).
#' @param upper_q upper quantile used in normalization. Default is 0.99 (99% quantile).
#' @param returnX a logical vairable. returnX = TRUE means returning the normalized covariate matrix. Default is FALSE.
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
#' sieve.model <- sieve_sgd_preprocess(X = TrainData[,2:(xdim+1)], type = type)
#' @export
#'
sieve_sgd_preprocess <- function(X, type = 'cosine', 
                             interaction_order = 3, 
                             norm_feature = TRUE, norm_para = NULL,
                             lower_q = 0.01, upper_q = 0.99,
                             returnX = FALSE){
  
  
  if(is(X,'numeric') | is(X, 'factor')){
    s.size <- length(X) #this is a univariate regression problem, sometimes the input data.frame is automatically conveted to an array
  }
  else{
    s.size <- dim(X)[1]
  }
  
  
  X <- sapply(X,unclass) #convert categorical variables into numeric variables
  X <- matrix(X, nrow = s.size)
  
  if(norm_feature){
    norm_list <- normalize_X(X, norm_para, lower_q, upper_q)
    X <- norm_list$X
    norm_para <- norm_list$norm_para
  }
  
  if(returnX){
    return(list(X = X, type = type, interaction_order = interaction_order, norm_para = norm_para))
  }else{
    return(list(type = type, interaction_order = interaction_order, norm_para = norm_para))
  }
}
