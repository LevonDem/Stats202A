"""

 Stat 202A - Homework 1
 Author: 
 Date : 
 Description: This script implements linear regression 
 using the sweep operator

 INSTRUCTIONS: Please fill in the missing lines of code
 only where specified. Do not change function names, 
 function inputs or outputs. You can add examples at the
 end of the script (in the "Optional examples" section) to 
 double-check your work, but MAKE SURE TO COMMENT OUT ALL 
 OF YOUR EXAMPLES BEFORE SUBMITTING.

 Do not use any of Python's built in functions for matrix 
 inversion or for linear modeling (except for debugging or 
 in the optional examples section).
 
"""

import numpy as np


################################
## Function 1: Sweep operator ##
################################

def mySweep(A, m):
    
    """
    Perform a SWEEP operation on A with the pivot element A[m,m].
    
    :param A: a square matrix (np.array).
    :param m: the pivot element is A[m, m].
    :returns a swept matrix (np.array). Original matrix is unchanged.
    """
    
    B = np.copy(A)   
    n = B.shape[0]
    for k in range(m):
        for i in range(n):
            for j in range(n):
                if i!=k and j!=k:
                    B[i,j] = B[i,j] - B[i,k]*B[k,j] / B[k,k]
        for i in range(n):
            if i!=k:
                 B[i,k] = B[i,k] / B[k,k]
        for j in range(n):
            if j!=k:
                B[k,j] = B[k,j] / B[k,k]
        B[k,k] = -1/B[k,k]
    
    ## The function outputs the matrix (np.array) B
    return(B)
  



########################################################
## Function 2: Linear regression using Sweep operator ##
########################################################

def myLinearRegression(X, Y):
  
  """
  Find the regression coefficient estimates beta_hat
  corresponding to the model Y = X * beta + epsilon
  Your code must use the sweep operator you coded above.
  Note: we do not know what beta is. We are only 
  given a matrix X and a vector Y and we must come 
  up with an estimate beta_hat.
  
  X: an 'n row' by 'p column' matrix (np.array) of input variables.
  Y: an n-dimensional vector (np.array) of responses

  FILL IN THE BODY OF THIS FUNCTION BELOW
  """
  
  ## Let me start things off for you...
  n = X.shape[0]
  p = X.shape[1]
  
  Z = np.hstack([np.ones([n, 1]), X, Y])
  A = np.dot(np.transpose(Z),  Z)
  S = mySweep(A, p+1)
  
  ## Function returns the (p+1)-dimensional vector (np.array) 
  ## beta_hat of regression coefficient estimates
  beta_hat = S[:(p+1), p+1]
  return beta_hat
  


########################################################
## Optional examples (comment out before submitting!) ##
########################################################

## def testing_Linear_Regression():
  
  ## This function is not graded; you can use it to 
  ## test out the 'myLinearRegression' function 

  ## You can set up a similar test function as was 
  ## provided to you in the R file.