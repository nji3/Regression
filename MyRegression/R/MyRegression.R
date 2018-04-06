#' QR Function
#'
#' This function allows you to do QR decomposition on the given matrix.
#' The output Q and R should satisfy: A=Q%*%R.
#' @param A The given matrix to do QR decomposition with
#' @export Q Orthogonal n by n matrix
#' @export R Upper triangular n by m matrix
#' myQR()

myQR <- function(A){
  
  n = nrow(A)
  m = ncol(A)
  R = A
  Q = diag(1,n,n)
  
  for (k in 1:(m-1)){
    x = as.matrix(rep(0,n))
    x[k:n,1] = R[k:n,k]
    v = x
    v[k] = x[k]+sign(x[k,1])*sqrt(sum(x^2))
    s = sqrt(sum(v^2))
    u = v/s
    R = R-2*u%*%t(u)%*%R
    Q = Q-2*u%*%t(u)%*%Q
  }
  
  return(list("Q" = t(Q), "R" = R))
  
}

#' Linear Regression Intercept Free Function
#'
#' This function allows you to do linear regression without intercept by QR decomposition.
#' Function returns the least squares solution vector.
#' @param X An n x p matrix of explanatory variables
#' @param Y An n dimensional vector of responses
#' @export beta_ls The least squares solution vector
#' myLM_0()

myLM_0 <- function(X, Y){
  
  n = dim(X)[1]
  p = dim(X)[2]
  Z = cbind(X,Y)
  R = myQR(Z)$R
  R1 = R[1:p,1:p]
  Y1 = R[1:p,(p+1)]
  beta_ls = solve(R1,Y1)
  
  return(beta_ls)
  
}

#' Linear Regression Function
#'
#' This function allows you to do linear regression with intercept by QR decomposition.
#' Function returns the least squares solution vector.
#' @param X An n x p matrix of explanatory variables
#' @param Y An n dimensional vector of responses
#' @export beta_ls The least squares solution vector
#' @export beta_ls_se The  standard errors of the least squares solutions
#' myLM()

myLM <- function(X, Y){
  
  n = dim(X)[1]
  p = dim(X)[2]
  Z = cbind(rep(1,n),X,Y)
  R = myQR(Z)$R
  R1 = R[1:(p+1),1:(p+1)]
  Y1 = R[1:(p+1),(p+2)]
  beta_ls = solve(R1,Y1)
  
  D = cbind(rep(1,n),X)
  Y_hat = D%*%beta_ls
  scale = solve(t(D)%*%D)
  beta_ls_se = sqrt(diag(sum((Y-Y_hat)^2)/(n-p)*scale))
  
  return(list(beta_ls = beta_ls, beta_ls_se = beta_ls_se))
  
}

#' Expit/sigmoid function
#'
#' This function used to compute the logit values for logistic regression function
#' @param x is a vector with probability values
#' expit()

expit <- function(x){
  1 / (1 + exp(-x))
}

#' Logistic Regression Function
#'
#' This function allows you to do Logistic regression with binomial family without intercept by solving iterated reweighed least squares.
#' Function returns the logistic regression solution vector.
#' @param X An n x p matrix of explanatory variables
#' @param Y An n dimensional vector of binary responses
#' @export beta The logistic regression solution vector
#' @export beta_se The  standard errors of the least squares solutions
#' myLogistic()

myLogistic <- function(X, Y){
  
  n = nrow(X)
  m = ncol(X)
  beta = rep(0,m)
  epsilon = 10^(-6)
  Continue = T
  
  while (Continue) {
    
    score = X%*%beta
    p = expit(score)
    w = p*(1-p)
    Y_hat = score+(Y-p)/w
    
    w_sqrt = unlist(lapply(w,sqrt))
    X_new = X*w_sqrt
    Y_new = Y_hat*w_sqrt
    beta_new = myLM_0(X_new,Y_new)
    if(sum(beta_new-beta)<=epsilon){
      Continue = F
    }
    beta = beta_new
  }
  
  p_hat = expit(X%*%beta)
  V = matrix(0,n,n)
  diag(V) = p_hat*(1-p_hat)
  beta_se = sqrt(diag(solve(t(X)%*%V%*%X)))
  
  return(list(beta = beta, beta_se = beta_se))
  
}

#' PCA Function
#'
#' This function allows you to perfrom PCA by using the QR function.
#' Function should output a list with D and V.
#' @param A A square matrix
#' @param numIter Number of iterations
#' @export D An vector of eigenvalues of A
#' @export V A matrix of eigenvectors of A (in the same order as the eigenvalues in D.)
#' myPCA()

myPCA <- function(A, numIter = 1000){
  
  T = numIter
  n = nrow(A)
  p = ncol(A)
  V = matrix(rnorm(n^2),n,n)
  for (i in 1:T) {
    Q = myQR(V)$Q
    R = myQR(V)$R
    V = A%*%Q
  }
  
  return(list("D" = diag(R), "V" = Q))
  
}

#' Sweep Function
#'
#' This function allows you to perfrom a Sweep operation on A with pivot element A[m,m].
#' Function should output a swept matrix.
#' @param A A square matrix
#' @param m The pivot element is A[m,m]
#' @export B A swept matrix
#' mySweep()

mySweep <- function(A, m){
  
  B <- A
  n <- nrow(B)
  
  for(k in 1:m){ 
    for(i in 1:n)     
      for(j in 1:n)   
        if(i != k  & j != k)     
          B[i,j] <- B[i,j] - B[i,k]*B[k,j]/B[k,k]    
        
        for(i in 1:n) 
          if(i != k) 
            B[i,k] <- B[i,k]/B[k,k]  
          
          for(j in 1:n) 
            if(j != k) 
              B[k,j] <- B[k,j]/B[k,k]
            
            B[k,k] <- - 1/B[k,k]
  }
  
  return(B)
  
}

#' Ridge Regression Function
#'
#' This function allows you to perfrom ridge regression of Y on X.
#' Function should return beta, the ridge regression solution.
#' @param X An n x p matrix of explanatory variables
#' @param Y An n vector of dependent variables or a matrix as long as the function works
#' @param lambda Regularization parameter (lambda >= 0)
#' @export beta_ridge A vector storing ridge regression solution
#' myRidge()

myRidge <- function(X, Y, lambda){
  
  n = dim(X)[1]
  p = dim(X)[2]
  Z = cbind(rep(1, n), X, Y) 
  A = t(Z) %*% Z
  D = diag(rep(lambda, p+2)) 
  D[1, 1] = 0
  D[p+2, p+2] = 0
  A=A+D
  S = mySweep(A, p+1)
  beta_ridge = S[1:(p+1), p+2]
  
  ## Function should output the vector beta_ridge, the 
  ## solution to the ridge regression problem. beta_ridge
  ## should have p + 1 elements.
  return(beta_ridge)
  
}

#' Spline Regression Function
#'
#' This function allows you to perfrom spline regression of Y on X.
#' Function should return a list containing the spline regression beta vector and the predicted Y values.
#' @param x An n x 1 vector or matrix of explanatory variables
#' @param Y An n x 1 vector of dependent variables or a matrix as long as the function works
#' @param lambda Regularization parameter (lambda >= 0)
#' @param p Number of cuts to make to the x-axis (default p = 100)
#' @export output A list containing the spline regression beta vector and the predicted Y values
#' mySpline()

mySpline <- function(x, Y, lambda, p = 100){
  
  n = length(x)
  x = sort(x)
  X = matrix(x, nrow=n)
  
  for (k in (1:(p-1))/p)
    X = cbind(X, (x>k)*(x-k))
  
  beta_spline = myRidge(X, Y, lambda)
  Yhat = cbind(rep(1, n), X)%*%beta_spline
  
  output <- list(beta_spline = beta_spline, predicted_y = Yhat)
  return(output)
  
}

#' Lasso solution path Function
#'
#' This function allows you to find the lasso solution path for various values of the regularization parameter lambda.
#' Function should return a matrix containing the lasso solution vector beta for each regularization parameter.
#' @param X An n x p matrix of explanatory variables
#' @param Y An n dimensional response vector
#' @param lambda_all A vector of regularization parameters (Make sure to sort lambda_all in decreasing order for efficiency)
#' @export beta_all The solution to the lasso regression problem for all the regularization parameters
#' myLasso()

myLasso <- function(X, Y, lambda_all){
  
  X = scale(X)
  X = cbind(rep(1,nrow(X)),X)
  Y = scale(Y)
  n = nrow(X)
  p = ncol(X)
  T = 500
  beta = matrix(rep(0, p), nrow = p)
  beta_all = matrix(rep(0, p*length(lambda_all)), nrow = p,ncol = length(lambda_all))
  
  L = length(lambda_all)
  R = Y
  ss = rep(0,p)
  for (j in 1:p) {
    ss[j] = sum(X[,j]^2)
  }
  for (l in 1:L) {
    lambda = lambda_all[l]
    for (t in 1:T) {
      db = sum(R*X[,1])/ss[1]
      b = beta[1]+db
      b = sign(b)*max(0,abs(b))
      db = b-beta[1]
      R = R-X[,1]*db
      beta[1] = b
      for (k in 2:p) {
        db = sum(R*X[,k])/ss[k]
        b = beta[k]+db
        b = sign(b)*max(0,abs(b)-lambda/ss[k])
        db = b-beta[k]
        R = R-X[,k]*db
        beta[k] = b
      }
    }
    beta_all[,l] = beta
  }
  
  par(mfrow = c(1, 2))
  matplot(t(matrix(rep(1, p), nrow = 1)%*%abs(beta_all)), t(beta_all), type = 'l',xlab='l1 norm',ylab='beta')
  
  Y_hat = X%*%beta_all
  e = c()
  for(i in 1:L) {
    k = sum(beta_all[,i]!=0)
    e[i] = sqrt(sum((Y_hat[,i]-Y)^2)/(n-k))
  }
  e = e[order(lambda_all)]
  plot(sort(lambda_all),e,xlab = 'Lambda',ylab = 'Estimation Error',main='Estimation Error against Lambda')
  lines(sort(lambda_all),e,col='red')
  
  return(beta_all)
  
}