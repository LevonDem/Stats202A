/*
#########################################################
## Stat 202A - Homework 3
## Author: 
## Date : 
## Description: This script implements QR decomposition,
## linear regression, and eigen decomposition / PCA 
## based on QR.
#########################################################
 
###########################################################
## INSTRUCTIONS: Please fill in the missing lines of code
## only where specified. Do not change function names, 
## function inputs or outputs. MAKE SURE TO COMMENT OUT ALL 
## OF YOUR EXAMPLES BEFORE SUBMITTING.
##
## Very important: Do not change your working directory
## anywhere inside of your code. If you do, I will be unable 
## to grade your work since R will attempt to change my 
## working directory to one that does not exist.
###########################################################
 
*/ 


# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


/* ~~~~~~~~~~~~~~~~~~~~~~~~~ 
 Sign function for later use 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

// [[Rcpp::export()]]
double signC(double d){
  return d<0?-1:d>0? 1:0;
}



/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 1: QR decomposition 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */  
  

// [[Rcpp::export()]]
List myQRC(const mat A){ 
  
  /*
  Perform QR decomposition on the matrix A
  Input: 
  A, an n x m matrix (mat)

  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
  */ 
  
  int n = A.n_rows;
  int m = A.n_cols;
  int s = 0;
  mat R = A;
  mat Q = eye(n, n);
  mat x(n, 1);
  mat v(n, 1);
  mat u(n, 1);
  List output;
  
  for(int k = 0; k < (m - 1); k++){
    
    x = 0 * x;
    for(int j = k; j < n; j ++){
      x(j, 0) = R(j, k);
    }
    s = -1 * signC(x(k, 0));
    v = x;
    v(k, 0) = x(k, 0) - s * norm(x);
    u = v / norm(v);
    
    R -= 2 * (u * (u.t() * R)); 
    Q -= 2 * (u * (u.t() * Q)); 
    
  }
  
  // Function should output a List 'output', with 
  // Q.transpose and R
  // Q is an orthogonal n x n matrix
  // R is an upper triangular n x m matrix
  // Q and R satisfy the equation: A = Q %*% R
  output["Q"] = Q.t();
  output["R"] = R;
  return(output);
  

}
  
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   Problem 2: Linear regression using QR 
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  
  
// [[Rcpp::export()]]
mat myLinearRegressionC(const mat X, const mat Y){
    
  /*  
  Perform the linear regression of Y on X
  Input: 
  X is an n x p matrix of explanatory variables
  Y is an n dimensional vector of responses
  Do NOT simulate data in this function. n and p
  should be determined by X.
  Use myQRC inside of this function
  
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
  
  */  
  
  // Let me start things off for you...
  int n = X.n_rows;
  int p = X.n_cols;
  
  mat z1(n, 1); 
  z1.fill(1);
  
  mat Z1 = join_rows(z1, X);
  mat Z  = join_rows(Z1, Y);
  mat R  = myQRC(Z)[1];
    
  mat R1 = R.submat(0, 0, p, p);
  mat Y1 = R.submat(0, p + 1, p, p + 1);

  // Function returns the 'p+1' by '1' matrix 
  // beta_ls of regression coefficient estimates
  mat beta_ls = inv(R1) * Y1;
  return(beta_ls.t());
  
}  

/* ~~~~~~~~~~~~~~~~~~~~~~~~ 
 Problem 3: PCA based on QR 
 ~~~~~~~~~~~~~~~~~~~~~~~~~~ */


// [[Rcpp::export()]]
List myEigen_QRC(const mat A, const int numIter = 1000){
  
  /*  
  
  Perform PCA on matrix A using your QR function, myQRC.
  Input:
  A: Square matrix
  numIter: Number of iterations
   
  #############################################
  ## FILL IN THE BODY OF THIS FUNCTION BELOW ##
  #############################################
   
   */  
  
  
  int r = A.n_rows;
  //mat V(r, r); //
  mat V = randu<mat>(r,r);
  mat Q(r, r);
  List QR_output;
  List output;

  for(int i = 0; i < numIter; i++){
    Q = as<mat>(myQRC(V)[0]);
    V = A * Q;
  }
  
  QR_output = myQRC(V);
  Q     = as<mat>(QR_output[0]);
  mat R = QR_output[1];
  vec D = R.diag();

  // Function should output a list with D and V
  // D is a vector of eigenvalues of A
  // V is the matrix of eigenvectors of A (in the 
  // same order as the eigenvalues in D.)
  output["D"] = D;
  output["V"] = Q;
  return(output);

}
  
