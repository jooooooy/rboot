#' A Cat Function

boot_resid <- function(Y,X,x_new,R,seed){
  
  require(matlib)
  
  n = dim(X)[1]
  m = dim(X)[2]
  
  beta_hat <- inv(t(X)%*%X)%*%t(X)%*%Y  # Lease square estimate of beta - a numeric vector of size m
  Y_hat <- X%*%beta_hat # Predicted Value of Y - a numeric vector of size n
  mat_hat <- X%*%inv(t(X)%*%X)%*%t(X) # Hat Matrix - a numeric matrix of dimension n-by-n
  
  h <- diag(mat_hat) # Leverage - a numeric vector of size n
  
  raw_err <- Y - Y_hat # Raw residual - a numeric vector of size n
  mod_err <- raw_err/sqrt(1-h) # modified residual - a numeric vector of size n
  cen_err <- mod_err - mean(mod_err) # centered modified residual - a numeric vector of size n
  
  y_new_hat <- x_new%*%beta_hat # prediction for new observation - a float number
  
  #######   Resampling Residual   #########
  
  # Resampling Index
  set.seed(seed)
  err_rsp <- c()
  for(i in 1:R){
    cen_err_rsp <- sample(cen_err,size = n, replace = TRUE)
    Y_rsp <- Y_hat + cen_err_rsp
    beta_rsp <- inv(t(X)%*%X)%*%t(X)%*%Y_rsp
    err_rsp[i] <- x_new%*%(beta_rsp-beta_hat)-sample(cen_err,size = 1)
  }
  
  return(list(beta_hat = beta_hat,
              Y_hat = Y_hat,
              raw_err = raw_err,
              mod_err = mod_err,
              y_new = y_new_hat, 
              interval=rep(y_new_hat,2)+quantile(err_rsp,c(0.025,0.975))))
}
