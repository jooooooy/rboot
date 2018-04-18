#' Resampling residual via boostrap
#'
#' This function return LSE of linear model Y = X beta + epsilon, the fitted value of Y,
#' the raw residual, the modified residual and the predicted value for new observation x_new,
#' as well as its 95% confidence interval,
#' @author Jingyi Lin
#' @param Y the vector of response
#' @param X the matrix of covariates (including intercept)
#' @param x_new a new observation of covariates
#' @param R the totoal times of resampling
#' @param seed the random number generator (RNG) state
#' @keywords Boostrap
#' @export beta_hat LSE of linear model Y = X beta + epsilon
#' @export Y_hat the fitted value of Y
#' @export raw_err the raw residual
#' @export mod_err the modified residual
#' @export y_new the predicted for new observation x_new
#' @export interval the 95% confidence interval for y_new
#' @examples
#' cat_function()
#' library(boot)
#' nuke <- nuclear[, c(1, 2, 5, 7, 8, 10, 11)]
#' Y <- log(nuke$cost) # True value of Y
#' X <- data.frame(1,nuke$date,log(nuke$cap),nuke$ne,nuke$ct,log(nuke$cum.n),nuke$pt)
#' X <- data.matrix(X) # Desgin matrix X
#' colnames(X) <- NULL 
#' x_new <- c(1,73,log(886),0,0,log(11),1) # new observation
#' R <- 1000 # Resampling sample size
#' seed <- 921 # Random seed - a integer

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
