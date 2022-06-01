#include <RcppArmadillo.h>
#include <string>
#include <cmath>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;


// double tryto(String s1){
//   if(s1 == "sobolev1"){
//     cout << 1;
//   }else{
//     cout << 0;
//
//   }
//
//   return -1;
// }

// part one: generate all the possible ways
// to factor an natural number into a product of dimlimit many numbers.

// defining vector of vector for
// storing factor combinations
vector<vector<int>>resultant;

// current_factor is current factor to be considered.
// current_product is current product of factors.
void factorsListFunc (int first, int each_prod, int n,
                      vector<int>single_result_list)
{
  // base case of this recursive function
  if (first > n || each_prod > n)
    return;

  // When current_product is equal to n,
  // then current contain list of factors
  // so it will be added to the vect
  if (each_prod == n)
  {
    resultant.push_back(single_result_list);
    return;
  }

  // In this loop we first calculate factors
  // of n and then it's combination so that
  // we get the value n in a recursive way .
  for (int i = first; i < n; i++)
  {
    if (i * each_prod > n)
      break;

    // if i divides n
    // properly then it
    // is factor of n
    if (n % i == 0)
    {

      // it is added to 'single_result_list' list
      single_result_list.push_back(i);

      // Here function is called recursively
      // and when (i*each_prod) is equal to n
      // we will store the 'single_result_list'
      // (it is basically the list containing all
      // combination of factors) into result_list.
      factorsListFunc(i, i * each_prod, n,
                      single_result_list);

      // here we will empty the 'single_result_list'
      // List so that new combination of factors
      // get stored in it.
      single_result_list.pop_back();
    }
  }
}

// Returns a list containing all ways
// to factorize a number n.
void factComb(int n)
{
  // making list of lists to store all
  // possible combinations of factors
  vector<int>single_result_list;

  // function to calculate all the combinations
  // of factors
  factorsListFunc(2, 1, n, single_result_list);
}

// [[Rcpp::export]]
List Generate_factors(int n, int dimlimit){

  // calling function for computing
  // vector of vector
  factComb(n);

  // printing all possible combination
  // of factors stored in vect

  List L = List::create();
  int k = 0;
  for (int i = 0; i < resultant.size(); i++)
  {
    vector<int> temp = resultant[i];
    if(temp.size() <= dimlimit){
      L.insert(k, temp);
      k++;
    }
    // for (int j = 0; j < resultant[i].size(); j++){
    //   cout << resultant[i][j] << " ";
    //  cout << endl;
    // }
  }
  resultant.clear();
  return L;
}

// This code is contributed by
// Atul kumar Shrivastava


///// part two: calculate basis functions and their tensor product.

// [[Rcpp::export]]
double psicos(double x, int j){
    return cos((j-1)*M_PI*x);
}

// [[Rcpp::export]]
double psisin(double x, int j){
  return sin((2*(j-1) - 1)*M_PI*x/2);
}

// [[Rcpp::export]]
double psipolytri(double x, int j){
  if(j ==2){
    return x;
  }else if(j % 2 == 1){
    return pow(2,0.5)*cos(2*M_PI*(j-1)/2*x);
  }else{
    return pow(2,0.5)*sin(2*M_PI*(j-2)/2*x);
  }
}

// [[Rcpp::export]]
double psipoly(double x, int j){
  return pow(x,j-1);
}

// [[Rcpp::export]]
double psi(double x, int j, String type){
  if(type == "sobolev1"){ // sobolev 1
    if(j == 1){
      return 1;
    }else{
      return sin((2*(j-1) - 1)*M_PI*x/2);
    }
  }else if(type == "cosine"){ // sobolev 1 using cosine function
    if(j == 1){
      return 1;
    }else{
      return cos((j-1)*M_PI*x);
    }
  }else if(type == "tri"){ // sobolev 1 using cosine function
    if(j == 1){
      return 1;
    }else if(j ==2){
      return x;
    }else if(j % 2 == 1){
      return pow(2,0.5)*cos(2*M_PI*(j-1)/2*x);
    }else{
      return pow(2,0.5)*sin(2*M_PI*(j-2)/2*x);
    }
  }else if(type == "legendre"){
    double y = 0;
    x = (x - 0.5)*2;
    if(j == 1){
      y = 1;
    }
    else if(j == 2){
      y = x;
    }
    else if(j == 3){
      y = (3*pow(x,2) - 1)/2;
    }
    else if(j == 4){
      y = (5*pow(x,3) - 3*x)/2;
    }
    else if(j == 5){
      y = (35*pow(x,4) - 30*pow(x,2) + 3)/8;
    }
    else if(j == 6){
      y = (63*pow(x,5) - 70*pow(x,3) + 15*x)/8;
    }
    else if(j == 7){
      y = (231*pow(x,6) - 315*pow(x,4) + 105*pow(x,2) - 5)/16;
    }
    else if(j == 8){
      y = (429*pow(x,7) - 693*pow(x,5) + 315*pow(x,3) - 35*x)/16;
    }
    else if(j == 9){
      y = (6435*pow(x,8) - 12012*pow(x,6) + 6930*pow(x,4) - 1260*pow(x,2) + 35)/128;
    }
    else{
      //std::cout << "trying to calculate " << j << "-th order Legendre polynomial" << "\n";
      //std::cout << "formula not provided" << "\n";
      y = -1;
    }
    return(y);
  }else{
    //std::cout << "basis function formula not provided" << "\n";
    return(-1);
  }
}

// [[Rcpp::export]]
double multi_psi(arma::vec x,arma::vec index,
                 String type){
  int xdim = x.n_rows;
  double psix = 1;
  for(int i = 0; i< xdim; i++){
    psix = psix * psi(x(i), index(i), type);
    // arma::cout << psix;
  }
  return psix;
}

double my_kernel(double x, double z, String type, double kernel_para){
  if(type == "sobolev1"){ //sobolev1
    return (1+std::min(x,z)) ;
  }
  else if(type == "gaussian"){
    return (exp(-kernel_para*pow(x-z, 2)));
  }
  else{
    //std::cout << "kernel type not specified" << "\n";
    return(-1) ;
  }
}

double tensor_kernel(arma::vec x, arma::vec z, String type, double kernel_para){
  int xdim = x.n_rows;
  double pik = 1;
  for (int i = 0; i < xdim; i++) {
    pik  = pik * my_kernel(x(i), z(i), type, kernel_para) ;
  }
  return pik;
}

// psi <- function(x, j, type){
// #this is univariate basis function
//   if(type == 'legendre'){
//     return(my.legendre(x, j))
//   }else if(type == 'sobolev1'){
//     if(j == 1){
//       return(1)
//     }else{
//       return(sin((2*(j-1) - 1)*pi*x/2))
//     }
//   }else if(type == 'polytri'){
//     if(j == 1){
//       return(1)
//     }
//     else if(j == 2){
//       return(x)
//     }
//     else if(j %% 2 == 1){
//       return(sin((j - 1)/2 *pi* x))
//     }
//     else{
//       return(cos((j - 2)/2 *pi* x))
//     }
//   }
// }



// [[Rcpp::export]]
List Kernel_M_C(arma::mat X, String type, double kernel_para){
  List ret ;
  int size = X.n_rows;
  arma::mat K = arma::mat(size, size);

  for (int i = 0; i < size; i++) {
    arma::vec Xrowi = arma::vectorise(X.row(i));
    for (int j = i; j < size; j++){
      K(i,j) = tensor_kernel(Xrowi,arma::vectorise(X.row(j)), type, kernel_para);
      K(j,i) = K(i,j);
    }
  }

  ret["K"] = K ;

  arma::mat U ;
  arma::vec s ;
  arma::mat V ;

  arma::svd(U, s, V, K) ;

  ret["U"] = U ;
  ret["s"] = s ;

  return(ret) ;
}

// [[Rcpp::export]]
arma::mat Design_M_C(arma::mat X, int basisN, String type, arma::mat index_matrix){
  int size = X.n_rows;
  int xdim = X.n_cols;
  
  arma::mat Phi = arma::mat(size, basisN);
  
  ////cosine basis
  if(type == "cosine"){
    for (int i = 0; i < size; i++) {
      // if (i % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
      arma::vec Xrowi = arma::vectorise(X.row(i));
      for (int j = 0; j < basisN; j++){ // loop for Jn many basis functions
        if (j % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
        // Phi(i,j) = multi_psi(Xrowi, arma::vectorise(index_matrix.row(j)),type);
        double psix = 1;
        for(int k = 0; k < xdim; k++){
          if(index_matrix(j,k) > 1){
            psix = psix * psicos(Xrowi(k), index_matrix(j,k));
          }
        }
        Phi(i,j) = psix;
      }
    }
  }else if(type == "sine"){
    for (int i = 0; i < size; i++) {
      // if (i % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
      arma::vec Xrowi = arma::vectorise(X.row(i));
      for (int j = 0; j < basisN; j++){ // loop for Jn many basis functions
        if (j % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
        // Phi(i,j) = multi_psi(Xrowi, arma::vectorise(index_matrix.row(j)),type);
        double psix = 1;
        for(int k = 0; k < xdim; k++){
          if(index_matrix(j,k) > 1){
            psix = psix * psisin(Xrowi(k), index_matrix(j,k));
          }
        }
        Phi(i,j) = psix;
      }
    }
  }else if(type == "polytri"){
    for (int i = 0; i < size; i++) {
      // if (i % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
      arma::vec Xrowi = arma::vectorise(X.row(i));
      for (int j = 0; j < basisN; j++){ // loop for Jn many basis functions
        if (j % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
        // Phi(i,j) = multi_psi(Xrowi, arma::vectorise(index_matrix.row(j)),type);
        double psix = 1;
        for(int k = 0; k < xdim; k++){
          if(index_matrix(j,k) > 1){
            psix = psix * psipolytri(Xrowi(k), index_matrix(j,k));
          }
        }
        Phi(i,j) = psix;
      }
    }
  }else if(type == "poly"){
    for (int i = 0; i < size; i++) {
      // if (i % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
      arma::vec Xrowi = arma::vectorise(X.row(i));
      for (int j = 0; j < basisN; j++){ // loop for Jn many basis functions
        if (j % 1000 == 0){Rcpp::checkUserInterrupt();} //interrupt when the data is too large
        // Phi(i,j) = multi_psi(Xrowi, arma::vectorise(index_matrix.row(j)),type);
        double psix = 1;
        for(int k = 0; k < xdim; k++){
          if(index_matrix(j,k) > 1){
            psix = psix * psipoly(Xrowi(k), index_matrix(j,k));
          }
        }
        Phi(i,j) = psix;
      }
    }
  }
  return Phi;
}

// [[Rcpp::export]]
arma::vec least_square_C(arma::mat Phi, arma::vec Y){
  arma::vec betahat = solve(Phi, Y);
  return betahat;

}

// [[Rcpp::export]]
arma::vec crossprod_C(arma::mat Phi, arma::vec betahat){
  arma::vec fhat = Phi * betahat;
  return fhat;

}






// [[Rcpp::export]]
arma::vec KRR_cal_beta_C(arma::mat U, arma::vec s,
                     double lambda, arma::vec Y){
  int size = U.n_rows;
  arma::mat D = arma::mat(size, size, arma::fill::zeros);

  for (int i = 0; i < size; i++){
    D(i,i) = 1/(s(i)+lambda);
  }

  arma::mat temp = U * D * U.t();

  arma::vec beta = temp * Y;
  return beta;
}

// [[Rcpp::export]]
arma::vec KRR_predict_C(arma::mat trainX, arma::mat testX,
                        String type, arma::vec beta_hat, double kernel_para){
  int trainsize = trainX.n_rows;
  int testsize = testX.n_rows;

  arma::mat Phi = arma::mat(testsize, trainsize);

  for (int i = 0; i < testsize; i++) {
    arma::vec Xrowi = arma::vectorise(testX.row(i));
    for (int j = 0; j < trainsize; j++){
      Phi(i,j) = tensor_kernel(Xrowi,arma::vectorise(trainX.row(j)), type, kernel_para);
    }
  }

  arma::vec predictY  = Phi * beta_hat;
  return predictY;
}

