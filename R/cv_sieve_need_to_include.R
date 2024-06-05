##cross-validated sieve estimators
library(data.table)
xdim <- 10
AllData <- GenSamples(s.size = 1000, xdim = xdim)
type <- 'cosine'
##

index_spliter <- function(array, n_folds = 2){
  
  # array <- 1:99
  
  # Calculate the length of each part
  part_length <- length(array) %/% n_folds
  
  # Create an empty list to store the parts
  parts <- vector("list", n_folds)
  
  # Randomly shuffle the array
  shuffled_array <- sample(array)
  
  # Split the shuffled array into parts
  for (fold_index in 1:n_folds) {
    start_index <- (fold_index - 1) * part_length + 1
    
    if(fold_index < n_folds){
      end_index <- fold_index * part_length
    }else{
      end_index <- length(array)
    }
    
    parts[[fold_index]] <- shuffled_array[start_index:end_index]
  }
  
  return(parts)
}
# basis_numbers <- ceiling(c(ncol(AllData) * c(5, nrow(AllData)^(1/5), nrow(AllData)^(1/3)),
#                          ncol(AllData)^2 * c(5, nrow(AllData)^(1/5), nrow(AllData)^(1/3))))

cross_validated_sieve <- function(AllData, 
                                  basis_numbers = NULL,
                                  n_folds = 5,
                                  type = 'cosine'){
  
  if(is.null(basis_numbers)){
    basis_numbers = ceiling(c(ncol(AllData) * c(5, nrow(AllData)^(1/5), nrow(AllData)^(1/3)),
                            ncol(AllData)^2 * c(5, nrow(AllData)^(1/5), nrow(AllData)^(1/3))))
  }
  
  validation_split_index <- index_spliter(1:nrow(AllData),
                                        n_folds = n_folds)
  parameter_tuning_reference <- data.frame()
  for(basisN in basis_numbers){
    print(paste0("trying basis number = ", basisN))
    sieve_fitting_lambda <- NULL 
    for(split_index in 1:n_folds){
      TrainData <- AllData[-validation_split_index[[split_index]], ]
      ValidationData <- AllData[validation_split_index[[split_index]], ]
    
  
      sieve.model <- sieve_preprocess(X = TrainData[,2:(xdim+1)], 
                                      basisN = basisN, type = type)
      sieve.fit<- sieve_solver(model = sieve.model, Y = TrainData$Y,
                               lambda = sieve_fitting_lambda)
      sieve_fitting_lambda <- sieve.fit$lambda #100 lambdas automatically determined by glmnet
      
  
      sieve.validation <- sieve_predict(model = sieve.fit, 
                                        testX = ValidationData[,2:(xdim+1)], 
                                        testY = ValidationData$Y)
  
      parameter_tuning_reference <- rbind(parameter_tuning_reference, 
                                          data.frame(l1_penalty = sieve.fit$lambda,
                                                     l1_penalty_index = 1:length(sieve.fit$lambda),
                                                     basisN = rep(basisN, length(sieve.fit$lambda)),
                                                     validation_MSE = sieve.validation$MSE,
                                                     split_index = rep(split_index, length(sieve.fit$lambda))))
      
    }
  }
  parameter_tuning_reference <- data.table(parameter_tuning_reference)
  average_data <- parameter_tuning_reference[, .(average_MSE = mean(validation_MSE),
                                                 l1_penalty = mean(l1_penalty)), #there are some round-off error when I combine between the folds
                                             by = .(l1_penalty_index, basisN)]
  
  best_combination_index <- which.min(average_data$average_MSE)
  best_lambda <- average_data$l1_penalty[best_combination_index]
  best_basis_number <- average_data$basisN[best_combination_index]
  return(list(best_lambda = best_lambda, 
              best_basis_number = best_basis_number))
}

basis_numbers <- ceiling(c(ncol(AllData) * c(5, nrow(AllData)^(1/5), nrow(AllData)^(1/3)),
                         ncol(AllData)^2 * c(5, nrow(AllData)^(1/5), nrow(AllData)^(1/3))))

best_hyperparameter <- cross_validated_sieve(AllData,
                                             basis_numbers = basis_numbers)
print(best_hyperparameter)
sieve.model <- sieve_preprocess(X = AllData[,2:(xdim+1)], 
                                basisN = best_hyperparameter$best_basis_number, 
                                type = type)
sieve.fit<- sieve_solver(model = sieve.model, Y = AllData$Y,
                         lambda = best_hyperparameter$best_lambda)

TestData <- GenSamples(s.size = 500, xdim = xdim)
sieve.test <- sieve_predict(model = sieve.fit, 
                                  testX = TestData[,2:(xdim+1)], 
                                  testY = TestData$Y)
plot(sieve.test$predictY, TestData$Y)
