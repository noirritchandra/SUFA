spca=function(Y,nv=35,pr=.9){
  pca.results=irlba::irlba(Y,nv=min(nv,nrow(Y)-1, ncol(Y)-1))
  cum_eigs= cumsum(pca.results$d)/sum(pca.results$d)

  d=raster::clamp(min(which(cum_eigs>pr)),2,nv)

  eta= pca.results$u [,1:d]%*% diag(pca.results$d [1:d])
  lambda=pca.results$v[,1:d]
  list(lambda=lambda,eta=eta)
}

initiate_MSFA_genedata=function(...,nv=50,ds=25,do.scale=T,do.center=F){
  Y=list(...)
  # Y=lapply(Y,scale,scale=F)
  cell.names=lapply(Y,rownames)
  S=length(Y)
  Y_comb=scale(Reduce(rbind,Y),scale=do.scale,center = do.center)


  pc_Y=spca(Y_comb,nv)
  J=tcrossprod( pc_Y$eta,pc_Y$lambda )

  E=Y_comb-J
  ns=sapply(Y,nrow)
  cs=cumsum(ns)

  lambda_list=rep(list(1),length(Y))
  l_list=rep(list(1),length(Y))
  f_list=rep(list(1),length(Y))

  for(j in 1:S){
    f_list[[j]]=pc_Y$eta[(cs[j]-ns[j]+1):cs[j],]

    l=spca(E[(cs[j]-ns[j]+1):cs[j],],ds)
    lambda_list[[j]] =l$lambda
    l_list[[j]]=l$eta
  }

  Y=lapply(cell.names,function(x,dat) dat[x,], dat=Y_comb)
  list(Y=Y,Phi=pc_Y$lambda,f=f_list,Lambda=lambda_list,l=l_list)
}


initiate_MSFA=function(...,nv=50,ds=25,pr=.9,do.scale=T,do.center=F){
  Y=list(...)
  Y=lapply(Y,scale,scale=do.scale,center=do.center)
  S=length(Y)
  Y_comb=Reduce(rbind,Y)


  pc_Y=spca(Y_comb,nv,pr)
  J=tcrossprod( pc_Y$eta,pc_Y$lambda )

  E=Y_comb-J
  ns=sapply(Y,nrow)
  cs=cumsum(ns)

  lambda_list=rep(list(1),length(Y))
  l_list=rep(list(1),length(Y))
  f_list=rep(list(1),length(Y))

  for(j in 1:S){
    f_list[[j]]=pc_Y$eta[(cs[j]-ns[j]+1):cs[j],]

    l=spca(E[(cs[j]-ns[j]+1):cs[j],],ds,.9)
    lambda_list[[j]] =l$lambda
    l_list[[j]]=l$eta
  }

  list(Y=Y,Phi=pc_Y$lambda,f=f_list,Lambda=lambda_list,l=l_list)
}


# msfa_covmat_est=function(msfa_HMC,burn=1e3){
#   S=length(msfa_HMC$A)
#
#   # covmat=replicate(S,list(1))
#   p=nrow(msfa_HMC$Phi[,,1])
#   rep=nrow(msfa_HMC$residuals)
#   sigmat=array(0,c(p,p,S))
#
#   for(j in (burn+1):rep){
#     phi_cross=tcrossprod(msfa_HMC$Phi[,,j]) + diag(msfa_HMC$residuals[j,])
#     for(s in 1:S)
#       sigmat[,,s]=sigmat[,,s]+ phi_cross + tcrossprod(msfa_HMC$Phi[,,j] %*% msfa_HMC$A[[s]][,,j])
#   }
#
#   return(sigmat/(rep-burn))
# }


#' Computes posterior means of the marginal covariance matrix from the output of \code{\link{fit_SUFA}}
#'
#' @param res Output of \code{\link{fit_SUFA}}.
#' @param burn Number of burnin samples to discard.
#'
#' @return An array with the marginal covariance matrices in the studies concatenated across the third dimension.
#' @export
sufa_marginal_covs=function(res,burn=5e2){
  # imsfa_marginal_covs=function(res,burn=5e2){
  library(plyr); library(abind)

  S=length(res$A)
  p=nrow(res$Lambda); q=ncol(res$Lambda)
  covmats=array(dim=c(p,p,S))
  lam_list=apply(res$Lambda,3,identity, simplify = F)
  resids=apply(res$residuals,1,identity, simplify = F)

  for(s in 1:S){
    A_list=apply( res$A[[s]],3,identity, simplify = F)
    covs=mapply(function(phi,A, ps){
      chol_a= chol( diag(q) +tcrossprod(A))
      tcrossprod(tcrossprod(phi, chol_a)) + diag(ps)
    }, lam_list, A_list, resids ,SIMPLIFY = F)
    nrep=length(covs)
    covmats[,,s]= Reduce("+",covs[-(1:burn)])/(nrep-burn)
  }

  covmats
}


msfa_marginal_covs=function(res,burn=0){
  S=length(res$Lambda)
  p=dim(res$Phi)[1]
  sigmat=array(0,dim=c(p,p,S))
  # nrep=dim(res$Phi)[3]

  phi_cross=aaply(res$Phi[,,-(1:burn)],3,tcrossprod,.parallel = F) ##stitches along the first dim
  for(s in 1:S){
    lam_cross=phi_cross+aaply(res$Lambda[[s]] [,,-(1:burn)],3,tcrossprod,.parallel = F)   +aaply(res$psi[[s]][,,-(1:burn)],2,diag,.parallel = F)
    # lam_cross_list=apply( lam_cross,1,identity,simplify = F)
    # sigmat[,,s]= Reduce("+",lam_cross_list)/length(lam_cross_list)
    sigmat[,,s]= aaply(lam_cross,c(2,3),mean,.parallel = F)
    print(s)
  }
  return(sigmat)
}


#' Compute the widely applicable Bayesian information criterion (WBIC)
#' 
#' Computes the widely applicable Bayesian information criterion (WBIC) \insertCite{wbic2013}{SUFA} from the output of \code{\link{fit_SUFA}}
#' 
#' @param res Output of \code{\link{fit_SUFA}}.
#' @param y List of datasets which was fed into \code{\link{fit_SUFA}}.
#' @param model The default value should be used.
#' @param burn Number of burnin samples to discard.
#' @param ncores Should be set at 1.
#'
#' @return The WBIC value calculated from the MCMC samples.
#' @references
#' \insertAllCited{}
#' @export
WBIC=function(res,y, model="SUFA",burn=5e2,ncores=1){
  #res=output of SUFA/MSFA
  #y=list of datasets fed into the models.
  library(mvtnorm)
  if(model=="SUFA"){
    S=length(res$A)
    p=nrow(res$Lambda); q=ncol(res$Lambda)
    # covmats=array(dim=c(p,p,S))
    wbic_s=numeric(S)
    lam_list=apply(res$Lambda,3,identity, simplify = F)
    resids=apply(res$residuals,1,identity, simplify = F)


    for(s in 1:S){
      A_list=apply( res$A[[s]],3,identity, simplify = F)
      covs=parallel::mcmapply(function(phi,A, ps){
        chol_a= chol( diag(q) +tcrossprod(A))
        tcrossprod(tcrossprod(phi, chol_a)) + diag(ps)
      }, lam_list, A_list, resids ,SIMPLIFY = F,mc.cores=ncores)
      nrep=length(covs)
      # covmats[,,s]= Reduce("+",covs[-(1:burn)])/(nrep-burn)
      # wbic_s[s]= mean(sapply(covs, dmvnorm, x=y[[s]],log=T, mean=rep(0,p) ))
      tmp=parallel::mclapply(covs, mvtnorm::dmvnorm, x=y[[s]],log=T, mean=rep(0,p) , mc.cores = ncores)
      wbic_s[s]= sum(Reduce("+",tmp)/length(tmp))
    }
    wbic=sum(wbic_s)
  }

  if(model=="MSFA"){
    library(plyr)
    S=length(res$Lambda)
    p=dim(res$Phi)[1]
    # sigmat=array(0,dim=c(p,p,S))
    # nrep=dim(res$Phi)[3]
    wbic_s=numeric(S)

    phi_cross=aaply(res$Phi[,,-(1:burn)],3,tcrossprod,.parallel = F) ##stiches along the first dim
    for(s in 1:S){
      covs=phi_cross+aaply(res$Lambda[[s]] [,,-(1:burn)],3,tcrossprod,.parallel = F)   +aaply(res$psi[[s]][,,-(1:burn)],2,diag,.parallel = F)
      covs=apply( covs,1,identity,simplify = F)
      # sigmat[,,s]= Reduce("+",lam_cross_list)/length(lam_cross_list)
      tmp=parallel::mclapply(covs, dmvnorm, x=y[[s]],log=T, mean=rep(0,p) , mc.cores = ncores)
      wbic_s[s]= sum(Reduce("+",tmp)/length(tmp))
      # wbic_s[s]= mean( apply(lam_cross,3,dmvnorm, x=y[[s]],log=T, mean=rep(0,p)))
      # print(s)
    }
    wbic=sum(wbic_s)
  }

  if(model=="PFA"){
    library(plyr)
    S=ncol(res$Pertmat[[1]]) ##nstudy

    Qmats=lapply(1:S, function(s, pertmat){
      lapply(pertmat, function(each_pert_mc,s){
        matrix(each_pert_mc[,s],nrow=sqrt(length(each_pert_mc[,s])))
      } ,s=s)
    },res$Pertmat)


    chol_shared_cov=mapply(function(lam, sig, sig_eta ){
      # t(chol(tcrossprod(lam %*% diag(sig_eta))+ diag(sig^2)))
      t(chol(tcrossprod(lam %*% diag(sqrt(sig_eta)))+ diag(sig)))
    } , res$Loading, res$Errorsigma, res$Latentsigma,SIMPLIFY = F)

    marg_covs=lapply(Qmats, function(Qmat_list,chol_shared_cov_list){
      nmc=length(Qmat_list)
      marg_cov_list= mapply(function(Qmat, ch_shared_cov){
        tcrossprod(solve(Qmat,ch_shared_cov))
      } ,Qmat_list,chol_shared_cov_list, SIMPLIFY = F)

    } , chol_shared_cov_list=chol_shared_cov)


    wbic=sum(mapply(function(covs,yy){
      sum(rowMeans(sapply(covs, mvtnorm::dmvnorm,x=yy,log=T,mean=rep(0,ncol(yy)))))
    } , marg_covs, y))
  }


  if(!(model %in% c("SUFA","MSFA","PFA")))
    stop("model NOT in SUFA,MSFA,PFA!!")

  return(-wbic)
}


# get_covmat_msfa_adj=function(phi, eta, lambda_s,ls, ps){
#   p=nrow(phi[,,1])
#   rep=nrow(ps)
#   sigmat=matrix(0,p,p)
#   for(j in 1:rep){
#     phi_adj=adjust_loading(phi[,,j],eta[,,j])
#     lambda_adj=adjust_loading(lambda_s[,,j],ls[,,j])
#     sigmat=sigmat+ tcrossprod(phi_adj)+ tcrossprod(lambda_adj) +diag(ps[j,])
#   }
#
#   return(sigmat/rep)
# }


shared_covmat=function(phi, resids,burn=5e2){
  phis=apply(phi,3,identity, simplify = F)
  diags=apply(resids,1,identity, simplify = F)

  covmats=mapply(function(phi, ps)  tcrossprod(phi)+diag(ps) ,phis,diags, SIMPLIFY = F)

  Reduce("+", covmats[-(1:burn)])/ (length(covmats)-burn )
}

#' Computes posterior mean of the shared covariance matrix across studies from the output of \code{\link{fit_SUFA}}
#' @details
#' Computes posterior mean of the shared covariance matrix \eqn{\mathbf{\Sigma}=\mathbf{\Lambda}\mathbf{\Lambda}^{T}+\mathbf{\Delta}} across studies from the output of \code{\link{fit_SUFA}}.
#' 
#' @param res Output of \code{\link{fit_SUFA}}.
#' @param burn Number of burnin samples to discard.
#'
#' @return A matrix representing the posterior mean of \eqn{\mathbf{\Sigma}=\mathbf{\Lambda}\mathbf{\Lambda}^{T}+\mathbf{\Delta}}.
#' @export
SUFA_shared_covmat=function(res,burn=5e2){
  shared_covmat(res$Lambda,res$residuals,burn=burn)
}

###Estimate loading matrix using infinitefactor package
#' Post-process MCMC samples of loading matrices to resole rotational ambiguity
#'
#' The varimax-based post-processing method of \insertCite{poworoznek2021;textual}{SUFA}. This function can be generally applied to any MCMC-based factor analysis model.
#'
#' @param lam_array An array of MCMC samples of the loading matrix. The MCMC samples MUST be stored along the third dimension.
#' @param burn Number of burnin samples to discard.
#' @param alpha The elements of the loading matrices are coded as 0 if their respective 100(1-\code{alpha})\% posterior credible intervals (after fixing for rotational ambiguity) include 0.
#'
#' @return A point estimate matrix 
#' @export 
#' @references
#' \insertAllCited{}
#'
#' @examples
lam.est=function(lam_array, burn=5e2,alpha=.05){
  loads =  apply(lam_array[,,-(1:burn)],3, function(lam) (varimax(lam)[[1]]),
                 simplify = F) ##MCMC samples MUST be stored across the third dimension

  norms = sapply(loads, norm, "2")
  pivot = loads[order(norms)][[median(1:length(norms))]]

  aligned = lapply(loads, infinitefactor::msf, pivot)
  (lam.est<-infinitefactor::summat(aligned,alpha = alpha))
  lam.est[,order(colSums( lam.est^2),decreasing = T)]
}


#' Estimate shared and study-specific loading matrices from the output of \code{\link{fit_SUFA}}
#'
#' @param res Output of \code{\link{fit_SUFA}}.
#' @param burn Number of burnin samples to discard.
#' @param alpha The elements of the loading matrices are coded as 0 if their respective 100(1-\code{alpha})\% posterior credible intervals (after fixing for rotational ambiguity) include 0.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{\code{Shared}}{Point estimate of the shared loading matrix \eqn{\mathbf{\Lambda}}.}
#' \item{\code{Study_specific}}{A list comprising the point estimates of the study specific loading matrices \eqn{\mathbf{\Lambda A}_s}.}
#' }
#' @export 
lam.est.all=function(res, burn=5e2,alpha=.05){
  # n.cores=round(n.cores)
  # if(n.cores<1)
  #   stop("n.cores<1")
  lam_array=res$Lambda[,,-(1:burn)]
  A_lists=lapply(res$A,function(x,burn) x[,,-(1:burn)], burn=burn )

  var_loads =  apply(lam_array,3, varimax, simplify = F) ##MCMC samples MUST be stored across the third dimension

  loads=lapply(var_loads, '[[',1) #the loading matrices following varimax
  rots=lapply(var_loads, '[[',2) #the varimax rotation matrix

  ######obtaining summary of the shared loading
  norms = sapply(loads, norm, "2")
  pivot = loads[order(norms)][[median(1:length(norms))]] #we use the same pivot later
  align_mats = lapply(loads, matchalign_permmat, pivot) #permutation matrices for alignment
  aligned=mapply("%*%", loads,  align_mats,SIMPLIFY = F) #aligned loadings
  lam.est<-infinitefactor::summat(aligned,alpha=alpha)
  lam.est=lam.est[,order(colSums( lam.est^2),decreasing = T)]
  ############

  ######obtaining summary of the study-specific loadings
  rotmats=mapply("%*%",rots, align_mats,SIMPLIFY = F)

  library(doParallel)
  # cl <- makeCluster(min(n.cores,length(res$A)), type="FORK")
  # registerDoParallel(cl)
  A_lists_lam_aligned=foreach(A =A_lists,
                              .packages = c("infinitefactor", "OpenMPController","inline")) %do%{
                                a=apply(A,3,identity,simplify = F)
                                as_list=mapply(crossprod, rotmats, a,SIMPLIFY = F)
                                var_loads =  lapply(as_list, varimax)
                                loads=lapply(var_loads, '[[',1) #the loading matrices following varimax
                                # rots=lapply(var_loads, '[[',2) #the varimax rotation matrix

                                pivot = loads[order(norms)][[median(1:length(norms))]]
                                aligned_s = lapply(loads, infinitefactor::msf, pivot)
                                aligned_lam_as =mapply("%*%", aligned,  aligned_s,SIMPLIFY = F)
                                lam_as.est<- infinitefactor::summat(aligned_lam_as,alpha=alpha)
                                lam_as.est[,order(colSums( lam_as.est^2),decreasing = T)]
                              }
  # stopCluster(cl)
  list(Shared=lam.est, Study_specific=A_lists_lam_aligned)
}


#############Estimating latent factors
get_latent_factors=function(res, Y, burn=500){
  library(plyr); library(abind)

  S=length(res$A)
  p=nrow(res$Lambda); q=ncol(res$Lambda)
  covmats=array(dim=c(p,p,S))
  lam_list=apply(res$Lambda,3,identity, simplify = F)
  resids=apply(res$residuals,1,identity, simplify = F)
  lambda=lam.est(res$Lambda, burn=burn)

  eta_list=Y

  for(s in 1:S){
    A_list=apply( res$A[[s]],3,identity, simplify = F)
    covs=mapply(function(phi,A, ps){
      tcrossprod( phi %*% A ) + diag(ps)
    }, lam_list[-(1:burn)], A_list[-(1:burn)], resids[-(1:burn)] ,SIMPLIFY = F)
    nrep=length(covs)
    covmats[,,s]= Reduce("+",covs)/(nrep)
    # sigmas_inv=  chol2inv(chol(covmats[,,s]))

    Ws=  solve(covmats[,,s], lambda)

    eta_list[[s]]=t(solve(crossprod(lambda, Ws) +diag(q ) , t(Y[[s]] %*% Ws) ))
  }
  return(eta_list)
}


#' Fit the SUFA model
#'
#' @param Y List of datasets. In each dataset independent samples MUST be stored along the rows. All of the datasets MUST comprise the same features.
#' @param qmax Maximum allowed column dimension of the shared factor loading matrix.
#' @param a Dirichlet-Laplace shrinkage parameter.
#' @param bAs Prior variance hyperparameter of the elements in the study-specific \eqn{\mathbf{A}_{s}} matrices.
#' @param ms Prior expectation of the idiosyncratic variance parameters \eqn{\delta_{j}^{2}}'s for all \eqn{j=1,\dots,d}.
#' @param ss Prior variance of the idiosyncratic variance parameters \eqn{\delta_{j}^{2}}'s for all \eqn{j=1,\dots,d}.
#' @param nthreads Number of parallel threads to use.
#' @param nrun Number of MCMC iterations to perform.
#' @param thin Every \code{thin}-th MCMC sample is stored.
#' @param nleapfrog The number of leapfrog steps in each MCMC iteration is sampled from \eqn{\min\{1,\text{Poisson}(\texttt{nleapfrog}) \} }.
#' @param leapmax The number of leapfrog steps is bounded above by \code{leapmax}.
#' @param del_range The step size of \eqn{\delta t} in  each HMC step is simulated as \eqn{ \delta t \sim} Uniform\code{(del_range[1],del_range[2])}.
#' @param col_prop The proportion of randomly selected columns of \eqn{\mathbf{\Lambda}} to update in each HMC step. If the dimension of the data is high, we suggest lowering this value from 1.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{\code{Lambda}}{An array with the MCMC samples of \eqn{\mathbf{\Lambda}}. MCMC samples are stored along the third dimension.}
#' \item{\code{A}}{A list of the MCMC samples of the study-specific \eqn{\mathbf{A}_s} matrices. The \eqn{s}-th element of the list corresponds to \eqn{\mathbf{A}_s} with MCMC samples stored in the same manner as that of \eqn{\mathbf{\Lambda}}.}
#' \item{\code{residuals}}{A matrix with the MCMC samples of \eqn{\mathrm{diag}(\mathbf{\Delta})=(\delta_{1}^{2},\dots,\delta_{d}^{2})^{T}}. The MCMC samples are stored along the rows of \code{residuals}.}
#' \item{\code{acceptance_probability}}{The acceptance probability of the HMC-within-Gibbs sampler.}
#' }
#' @export
#' @examples
#' vignette(topic="Intro_simulation",package = "SUFA")
fit_SUFA=function(Y,qmax,a=0.5, bAs=1, ms=1,ss=7, nthreads=6,
                   nrun=7.5e3, thin=5,
                  nleapfrog=6,leapmax=10,del_range=c(.001,.0075),col_prop=1){

  Y_comb=Reduce(rbind,Y)
  p=ncol(Y_comb)
  init.vals=spca(Y_comb, pr=.95,nv=qmax)
  d_est=ncol(init.vals$lambda)
  nstudy=length(Y)
  ds_est= rep(floor(d_est/nstudy) ,nstudy)  #sapply(res$l,ncol)
  Phi.init=init.vals$lambda
  #################################################

  RcppArmadillo::armadillo_set_number_of_omp_threads(nthreads);
  SUFA_HMC(nrun=nrun, thin=thin, nleapfrog=nleapfrog,del_range=del_range,
                                      ps_hyper=c(ms,ss),  A_hyper=c(0,bAs), a=a,
                                      Y, ks = ds_est, Phi.init,leapmax=leapmax,leapmin=5,col_prob=col_prop ,nthreads = min(nstudy,nthreads) )
}

#' Fit a Bayesian sparse factor analysis model
#'
#' A fast Hamiltonian Monte Carlo-based sparse Bayesian factor model implementation
#'\deqn{Y_{i}\sim \mathrm{N} (0, \mathbf{\Lambda}\mathbf{\Lambda}^{T}+\mathbf{\Delta}) \text{ where } \mathrm{diag}(\mathbf{\Delta})=(\delta_{1}^{2},\dots,\delta_{d}^{2})^{T}}
#'\deqn{\mathrm{vectorise}(\mathbf{\Lambda})\sim \mathrm{Dirichlet-Laplace}(a), }
#'\deqn{\log\delta_{j}^{2} \sim \mathrm{N}(\mu_{\delta},\sigma_{\delta}^{2}), } such that \eqn{\mu_{\delta}} and \eqn{\sigma_{\delta}^{2}} are chosen such that a priori \eqn{E(\delta_{j}^{2})=}\code{ms} and \eqn{\mathrm{var}(\delta_{j}^{2})=}\code{ss} for all \eqn{j=1,\dots,d}.
#'
#' @param Y Data matrix
#' @param qmax Maximum allowed column dimension of the factor loading matrix.
#' @param a Dirichlet-Laplace shrinkage parameter.
#' @param ms Prior expectation of the idiosyncratic variance parameters \eqn{\delta_{j}^{2}}'s for all \eqn{j=1,\dots,d}.
#' @param ss Prior variance of the idiosyncratic variance parameters \eqn{\delta_{j}^{2}}'s for all \eqn{j=1,\dots,d}.
#' @param nrun Number of MCMC iterations to perform.
#' @param thin Every \code{thin}-th MCMC sample is stored.
#' @param nleapfrog The number of leapfrog steps in each MCMC iteration is sampled from \eqn{\min\{1,\text{Poisson}(\texttt{nleapfrog}) \} }.
#' @param leapmax The number of leapfrog steps is bounded above by \code{leapmax}.
#' @param del_range The step size of \eqn{\delta t} in  each HMC step is simulated as \eqn{ \delta t \sim} Uniform\code{(del_range[1],del_range[2])}.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{\code{Lambda}}{An array with the MCMC samples of \eqn{\mathbf{\Lambda}}. MCMC samples are stored along the third dimension.}
#' \item{\code{residuals}}{A matrix with the MCMC samples of \eqn{\mathrm{diag}(\mathbf{\Delta})=(\delta_{1}^{2},\dots,\delta_{d}^{2})^{T}}. The MCMC samples are stored along the rows of \code{residuals}.}
#' \item{\code{acceptance_probability}}{The acceptance probability of the HMC-within-Gibbs sampler.}
#' }
#' @export
#'
#' @examples
#' vignette(topic="sparse_BFA",package = "SUFA")
fit_FA=function(Y,qmax,a=0.5, ms=1,ss=7, 
                  nrun=7.5e3, thin=5,
                  nleapfrog=6,leapmax=10,del_range=c(.001,.0075)){
  p=ncol(Y)
  init.vals=spca(Y, pr=.95,nv=qmax)
  Phi.init=init.vals$lambda
  #################################################
  RcppArmadillo::armadillo_set_number_of_omp_threads(6);
  cov_est_HMC( a, c(ms, ss),
               nrun=nrun, thin=thin, nleapfrog=nleapfrog, del_range=del_range,
               phimat=Phi.init, Y=Y,leapmin=3,leapmax=leapmax)
}

#' Heatmap of a matrix
#'
#' @param x A matrix
#'
#' @return \code{ggplot} object representing the heatmap of \code{x}
#' @export
#'
#' @examples
#' x=matrix(rnorm(50*5),nrow=50)
#' x[sample.int(length(x), 100 )]=NA
#' heatplot(x)
heatplot=function(x){
  library(reshape2);library(ggplot2)
  melted_x <- reshape2::melt(t(x))
  
  ggplot(data = melted_x, aes(x=Var1, y=Var2, fill=value))+
    geom_tile(aes(fill = value)) +theme (panel.grid.major = element_blank(),panel.border = element_blank(),panel.background = element_blank())+
    scale_fill_viridis_c(option="B",alpha=1)
}

#' Gene expression datasets from the Immgen project
#' 
#'  
#' We integrate data from three studies analyzing gene expressions. 
#' The first is the GSE109125 bulkRNAseq dataset, collected from 103 highly purified immunocyte populations representing all lineages and several differentiation cascades and profiled using the ImmGen  pipeline \insertCite{bulk_train}{SUFA}. 
#' The second study is  a microarray dataset GSE15907 \insertCite{micro1_train1,micro1_train2}{SUFA}, measured on multiple \emph{ex-vivo} immune lineages, primarily from adult B6 male mice. 
#' Finally, we include the GSE37448 \insertCite{micro2_train}{SUFA} microarray dataset, also  part of the Immgen project. 
#' 
#' @docType data
#' @format A \code{list} with the following elements:
#' \describe{
#' \item{bulk}{The GSE109125 bulkRNASeq data in gene x cell format}
#' \item{array1}{First microarray dataset GSE15907 data in gene x cell format}
#' \item{array2}{Second microarray dataset GSE37448 data in gene x cell format}
#' \item{bulk.types}{A string vector indicating the cell-types in the bulkRNASeq data}
#' \item{array1.types}{A string vector indicating the cell-types in the first microarray data}
#' \item{array2.types}{A string vector indicating the cell-types in the second microarray data}
#' }
#' 
#' @source {GSE109125, GSE15907 and GSE37448}
#' @note In these datasets, we only include those cell-types with at least 3 observations. We have also \code{log2}-transformed the original expression-count data.
#' @references
#' \insertAllCited{}
#' @examples
#' download_genedata()

download_genedata=function(){
  library(httr)
  url <- "https://utdallas.box.com/shared/static/tuqwc8i0mzixs83wkvtla7qx0sg34365.rda"
  temp <- tempfile(fileext = ".rda")
  
  # Download the file to a temporary location
  httr::GET(url, write_disk(temp, overwrite = TRUE))
  
  # Load the file from the temporary location
  load(temp)
  
  # Clean up the temporary file
  unlink(temp)
  
  return(genedata)
}
