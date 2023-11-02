## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----set_seed-----------------------------------------------------------------
set.seed(35)

## ----load_SUFA----------------------------------------------------------------
library(SUFA)

## ----load_R_funcs, message=FALSE, results=FALSE, warning=FALSE, comment=FALSE----
library(here)
source(here("vignettes","simulate_data_fxns.R")) #We load some R functions

## ----set_dimension------------------------------------------------------------
S=5 ## no. of studies
d= 200 ## observed dimension
q= 20 ## latent dimension of the shared subspace
n=d;ns=sapply(rpois(S,n/S),max,n/S) ## sample-sizes in each study

## ----set_dimension_latent-----------------------------------------------------
qs=sapply(rpois(S,floor(q/S))+1,min,floor(q/S))

## ----gen_lambda---------------------------------------------------------------
val=2
  lam=simulate_lambda_sparse (k = q,p=d,pr = .25,val = val)
  blank_rows=which(rowSums(lam)==0)
  for(bl in blank_rows){
    lam[bl,sample.int(q,5)]= runif(5,-val,val)
  }

## ----plot_lambda, fig.dim = c(3, 3.75)----------------------------------------
library(ggplot2)
heatplot(lam)

## ----idiosyncratic------------------------------------------------------------
diag_sd=sqrt(.5) ##sd

## ----null_lam-----------------------------------------------------------------
null_lam=MASS::Null(lam) ##null-space of lam

Y=replicate(S,list(1))
    lambda=replicate(S,list(1))
    covmat=replicate(S,list(1))
    for(s in 1:S){
      eta=matrix(rnorm(ns[s]*q), nrow=ns[s],ncol=q)
      
      lambda_s=matrix (rnorm(q*qs[s],sd=.25),  nrow=q,ncol=qs[s])
      lambda[[s]]=lam %*% lambda_s + matrix(rnorm(d*qs[s],sd=.1),nrow = d,ncol=qs[s])
      
      ls=matrix(rnorm(ns[s]*qs[s]), nrow=ns[s],ncol=qs[s])
      
      Y[[s]] =tcrossprod(eta,lam)+tcrossprod(ls,lambda[[s]]) + matrix(rnorm(d*ns[s],sd=diag_sd),nrow=ns[s],ncol=d)
    }
    
    covmat_shared=tcrossprod(lam) + diag_sd*diag(d) ##True shared covariance matrix

## ----fit_SUFA, message=FALSE, results=FALSE, warning=FALSE, comment=FALSE-----
time.taken=system.time(res<-fit_SUFA(Y,qmax=25,nthreads = 5,nrun = 7.5e3,
                                     nleapfrog = 4, leapmax = 9, del_range = c(0,.01)))

## -----------------------------------------------------------------------------
time.taken

res$acceptance_probability

## ----set_burnin---------------------------------------------------------------
burnin=500

## ----shared_sigma_mean--------------------------------------------------------
covmat_shared_est=SUFA_shared_covmat(res,burn = burnin)

## ----plot_sigmas, fig.dim = c(8, 5)-------------------------------------------
cormat_true=cov2cor(covmat_shared)
cormat_est=cov2cor(covmat_shared_est)
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----plot_sigmas_sparse, fig.dim = c(8, 5)------------------------------------
cormat_true[abs(cormat_true)<.1] =NA
cormat_est[abs(cormat_est)<.1] =NA
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----remove_ggplots, echo=FALSE-----------------------------------------------
rm(sig_true,sig_est,cormat_true,cormat_est)

## ----posterior_lambda---------------------------------------------------------
est_lam=lam.est(res$Lambda,burn = burnin)

## ----plot_lambdass, fig.dim = c(6, 5)-----------------------------------------
lam.sp=lam; lam.sp[lam==0]=NA
est_lam.sp=est_lam; est_lam.sp[est_lam==0]=NA

lam_true.gg=heatplot(lam.sp)+labs(title="True")
lam_est.gg=heatplot(est_lam.sp)+labs(title="Estimated")
ggpubr::ggarrange(lam_true.gg,lam_est.gg,ncol = 2)

## ----echo=F-------------------------------------------------------------------
rm(lam.sp,est_lam.sp, lam_true.gg, lam_est.gg)

## ----coeff_r2-----------------------------------------------------------------
fit_lm=summary(lm(lam~est_lam)) ## regress true lambda on estimated lambda
summary(sapply(fit_lm, function(x) x$r.squared)) ## report the r-squares r-squares

## ----wbic---------------------------------------------------------------------
WBIC(res ,Y, model="SUFA",burn=5e2,ncores=1)

## ----gen_lambda_FM2-----------------------------------------------------------
val=2
  lam=simulate_lambda_sparse2 (k = q,p=d,pr = .5,val = val)
  blank_rows=which(rowSums(lam)==0)
  for(bl in blank_rows){
    lam[bl,sample.int(q,5)]= runif(5,-val,val)
  }

## ----plot_lambda_FM2, fig.dim = c(3, 3.75)------------------------------------
library(ggplot2)
heatplot(lam)

## ----idiosyncratic_FM2--------------------------------------------------------
diag_sd=sqrt(.5) ##sd

## ----null_lam_FM2-------------------------------------------------------------
null_lam=MASS::Null(lam) ##null-space of lam

Y=replicate(S,list(1))
    lambda=replicate(S,list(1))
    covmat=replicate(S,list(1))
    for(s in 1:S){
      eta=matrix(rnorm(ns[s]*q), nrow=ns[s],ncol=q)
      
      lambda_s=matrix (rnorm(q*qs[s],sd=.25),  nrow=q,ncol=qs[s])
      lambda[[s]]=lam %*% lambda_s + matrix(rnorm(d*qs[s],sd=.1),nrow = d,ncol=qs[s])
      
      ls=matrix(rnorm(ns[s]*qs[s]), nrow=ns[s],ncol=qs[s])
      
      Y[[s]] =tcrossprod(eta,lam)+tcrossprod(ls,lambda[[s]]) + matrix(rnorm(d*ns[s],sd=diag_sd),nrow=ns[s],ncol=d)
    }
    
    covmat_shared=tcrossprod(lam) + diag_sd*diag(d) ##True shared covariance matrix

## ----fit_SUFA_FM2, message=FALSE, results=FALSE, warning=FALSE, comment=FALSE----
time.taken=system.time(res<-fit_SUFA(Y,qmax=25,nthreads = 5,nrun = 7.5e3,
                                     nleapfrog = 4, leapmax = 9, del_range = c(0,.01)))

## -----------------------------------------------------------------------------
time.taken

res$acceptance_probability

## ----set_burnin_FM2-----------------------------------------------------------
burnin=500

## ----shared_sigma_mean_FM2----------------------------------------------------
covmat_shared_est=SUFA_shared_covmat(res,burn = burnin)

## ----plot_sigmas_FM2, fig.dim = c(8, 5)---------------------------------------
cormat_true=cov2cor(covmat_shared)
cormat_est=cov2cor(covmat_shared_est)
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----plot_sigmas_sparse_FM2, fig.dim = c(8, 5)--------------------------------
cormat_true[abs(cormat_true)<.1] =NA
cormat_est[abs(cormat_est)<.1] =NA
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----remove_ggplots_FM2, echo=FALSE-------------------------------------------
rm(sig_true,sig_est,cormat_true,cormat_est)

## ----posterior_lambda_FM2-----------------------------------------------------
est_lam=lam.est(res$Lambda,burn = burnin)

## ----plot_lambdass_FM2, fig.dim = c(6, 5)-------------------------------------
lam.sp=lam; lam.sp[lam==0]=NA
est_lam.sp=est_lam; est_lam.sp[est_lam==0]=NA

lam_true.gg=heatplot(lam.sp)+labs(title="True")
lam_est.gg=heatplot(est_lam.sp)+labs(title="Estimated")
ggpubr::ggarrange(lam_true.gg,lam_est.gg,ncol = 2)

## ----echo=F-------------------------------------------------------------------
rm(lam.sp,est_lam.sp, lam_true.gg, lam_est.gg)

## ----coeff_r2_FM2-------------------------------------------------------------
fit_lm=summary(lm(lam~est_lam)) ## regress true lambda on estimated lambda
summary(sapply(fit_lm, function(x) x$r.squared)) ## report the r-squares r-squares

## ----wbic_FM2-----------------------------------------------------------------
WBIC(res ,Y, model="SUFA",burn=5e2,ncores=1)

## ----gen_lambda_FM3-----------------------------------------------------------
lam=matrix(0,nrow = d, ncol = q)
  nonzero=floor(d*.15); val=2
  for(i in 1:q){
    lam[sample.int(size=nonzero,n = d) ,i]=runif(nonzero,-val,val)
  }
  blank_rows=which(rowSums(lam)==0)
  for(bl in blank_rows){
    lam[bl,sample.int(q,5)]= runif(5,-val,val)
  }

## ----plot_lambda_FM3, fig.dim = c(3, 3.75)------------------------------------
library(ggplot2)
heatplot(lam)

## ----idiosyncratic_FM3--------------------------------------------------------
diag_sd=sqrt(.5) ##sd

## ----null_lam_FM3-------------------------------------------------------------
null_lam=MASS::Null(lam) ##null-space of lam

Y=replicate(S,list(1))
    lambda=replicate(S,list(1))
    covmat=replicate(S,list(1))
    for(s in 1:S){
      eta=matrix(rnorm(ns[s]*q), nrow=ns[s],ncol=q)
      
      lambda_s=matrix (rnorm(q*qs[s],sd=.25),  nrow=q,ncol=qs[s])
      lambda[[s]]=lam %*% lambda_s + matrix(rnorm(d*qs[s],sd=.1),nrow = d,ncol=qs[s])
      
      ls=matrix(rnorm(ns[s]*qs[s]), nrow=ns[s],ncol=qs[s])
      
      Y[[s]] =tcrossprod(eta,lam)+tcrossprod(ls,lambda[[s]]) + matrix(rnorm(d*ns[s],sd=diag_sd),nrow=ns[s],ncol=d)
    }
    
    covmat_shared=tcrossprod(lam) + diag_sd*diag(d) ##True shared covariance matrix

## ----fit_SUFA_FM3, message=FALSE, results=FALSE, warning=FALSE, comment=FALSE----
time.taken=system.time(res<-fit_SUFA(Y,qmax=25,nthreads = 5,nrun = 7.5e3,
                                     nleapfrog = 4, leapmax = 9, del_range = c(0,.01)))

## -----------------------------------------------------------------------------
time.taken

res$acceptance_probability

## ----set_burnin_FM3-----------------------------------------------------------
burnin=500

## ----shared_sigma_mean_FM3----------------------------------------------------
covmat_shared_est=SUFA_shared_covmat(res,burn = burnin)

## ----plot_sigmas_FM3, fig.dim = c(8, 5)---------------------------------------
cormat_true=cov2cor(covmat_shared)
cormat_est=cov2cor(covmat_shared_est)
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----plot_sigmas_sparse_FM3, fig.dim = c(8, 5)--------------------------------
cormat_true[abs(cormat_true)<.1] =NA
cormat_est[abs(cormat_est)<.1] =NA
sig_true=heatplot(cormat_true)+ theme(legend.position='bottom')+labs(title="True")
sig_est=heatplot(cormat_est)+ theme(legend.position='bottom')+labs(title="Estimated")
ggpubr::ggarrange(sig_true,sig_est,ncol = 2)

## ----remove_ggplots_FM3, echo=FALSE-------------------------------------------
rm(sig_true,sig_est,cormat_true,cormat_est)

## ----posterior_lambda_FM3-----------------------------------------------------
est_lam=lam.est(res$Lambda,burn = burnin)

## ----plot_lambdass_FM3, fig.dim = c(6, 5)-------------------------------------
lam.sp=lam; lam.sp[lam==0]=NA
est_lam.sp=est_lam; est_lam.sp[est_lam==0]=NA

lam_true.gg=heatplot(lam.sp)+labs(title="True")
lam_est.gg=heatplot(est_lam.sp)+labs(title="Estimated")
ggpubr::ggarrange(lam_true.gg,lam_est.gg,ncol = 2)

## ----echo=F-------------------------------------------------------------------
rm(lam.sp,est_lam.sp, lam_true.gg, lam_est.gg)

## ----coeff_r2_FM3-------------------------------------------------------------
fit_lm=summary(lm(lam~est_lam)) ## regress true lambda on estimated lambda
summary(sapply(fit_lm, function(x) x$r.squared)) ## report the r-squares r-squares

## ----wbic_FM3-----------------------------------------------------------------
WBIC(res ,Y, model="SUFA",burn=5e2,ncores=1)

