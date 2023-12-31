library(mvtnorm)
# library(tidyverse)

simulate_x <- function(n, p, k, sigmasq, lambda_generator, return_lambda = FALSE){
  lambda = lambda_generator(k, p)
  sigma = sigmasq*(diag(p) + tcrossprod(lambda))
  # sigma = cor(sigma) * sigmasq
  x = rmvnorm(n = n, mean = rep(0, p), sigma = sigma)
  if(return_lambda){
    return(list("x" = x, "lambda" = lambda, "sigma" = sigma))
  }
  x
}

simulate_lambda <- function(k, p, sigmasq){
  sapply(1:k, function(j) simulate_lambda_column(p, j))
}


simulate_lambda_sparse <- function(k, p, pr,val=2.5){
  sapply(1:k, function(j) simulate_lambda_column_sp(p, pr,val))
}

simulate_lambda_sparse2 <- function(k, p, pr,val=2.5){
  sapply(1:k, function(j) c(simulate_lambda_column_sp(p/2, pr/2,val),simulate_lambda_column_sp(p/2, pr/2,val) ) ) 
}

simulate_lambda_sparse3 <- function(k, p, val=2.5){
  props=seq(.65,.2,length.out=k)
  loading=sapply(props,function(prop,p,val) {
    nonzero=floor(p*prop)
    c(runif(nonzero,min=-val,max=val),rep(0,p-nonzero))
    } , p=p,val=val)
  loading[,sample.int(k,k)]
}


simulate_lambda_column_sp <- function(p, pr, val=2.5){
  k=floor(p*pr)
  x=numeric(p)
  j= sample.int(p-k,1)  #sample.int(floor(1/pr),1 )
  
  x[(j-1)+ (1:k)  ]= runif(k,-val,val)
  # x[  ((j-1)*k+1) : (j*k)]= runif(k,-val,val)
  
  return (x);
}


simulate_lambda_column <- function(p, j){
  value = runif(n = p, min = .5, max = 1) * sample(c(-1, 1), size = p, replace=TRUE)
  nonzero = rbinom(n = p,size=1, p = .4 + .2/j)
  # nonzero = rbernoulli(n = p, p = .4 + .2/j)
  value[!nonzero] = 0
  value
}

simulate_lambda_EFA_example <- function(k, p){
  sapply(1:k, function(j) simulate_lambda_EFA_example_column(k, p, j))
}

simulate_lambda_EFA_example_column <- function(k, p, j){
  map_dbl(1:p, function(d) assign_lambdadj(k, p, j, d))
}

rule1 <- function(k, p, j, d){
  ((1 + (j-1)*(p/k)) <= d) & (d <= j*(p/k))
}

rule2 <- function(k, p, j, d){
  ((1 + j*(p/k)) <= d) & (d <= (j+1)*(p/k)) & (j <= (k - 1))
}

rule3 <- function(k, p, j, d){
  (((j-1)*(p/k)) <= d) & (d <= (j*(p/k) - 1)) & (2 <= j)
}

assign_lambdadj <- function(k, p, j, d){
  if(rule1(k, p, j, d)){
    2*(6 - j)
  } else if(rule2(k, p, j, d) || rule3(k, p, j, d)){
    -2*(6 - j)
  } else{
    0
  }
}