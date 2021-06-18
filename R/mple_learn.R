library(doParallel)
library(Matrix)

mple_learn <- function(
  config
){
  
}


#' Update theta using New's method
#' @param config a list, admm_alpha
#' @param theta
#' @return theta^{t+!} 
theta_update_hess <- function(
  config,
  theta, 
  X
){
  cl <- makeCluster(8)
  registerDoParallel(cl)
  # compute data
  H <- foreach(i=1:config$tslen, .packages = c('network', 'ergm')) %dopar%{
    if (i == 1 | i == config$tslen){
      gstat_delta(config, yt = X[[i]], yt0 = X[[i]])
    }else{
      gstat_delta(config, yt = X[[i+1]], yt0 = X[[i]])
    }
  }
  stopCluster(cl)
  # compute hessian
  # estimated logodds 
  Z = lapply(H, function(A, x) A%*%x, A = theta)
  mu = lapply(Z, sigmoid)
  W = lapply(mu, function(x) x*(1-x))
  hess0 = mapply(function(A, B) diag(as.numeric(A))%*%t(B), A=W, B=H, SIMPLIFY = F)
  hess1 = mapply(function(A, B) A%*%B, A=H, B=hess0, SIMPLIFY = F)
  hess = bdiag(hess1) + config$admm_alpha * diag(p*config$tslen)
  return(list(hess=hess, H=H, W=W, Z=Z, mu=mu))  
}



#' Compute the change stat in mple for all edges  (This function is not necessary. Use ergmMPLE)
#'
#' @param yt 
#' @param yt0 
#'
#' @return gdelta, a matrix of dim (p1+p2) x n^2
#' @export
#'
#' @examples
gstat_delta <- function(
  config,
  yt0,
  yt
){
  n = dim(yt)[1]
  p = length(config$myf_form) + length(config$myf_diss)
  gdelta <- matrix(0, nrow=p, ncol=n^2)
  for(j in 1:n){
    for(i in 1:n){
      if (i != j){
        y1 <- yt
        y1[i, j] <- 1
        y2 <- yt
        y2[i, j] <-  0  
        g1 = gstat(config, y1, yt0)
        g2 = gstat(config, y2, yt0)
        gdelta[, (j-1)*n+i] =  g1 - g2 
      }
    }
  }
  return(gdelta)
}



#' Title
#'
#' @param config 
#' @param X 
#' @param theta 
#' @param hess_res 
#'
#' @return
#' @export
#'
#' @examples
theta_update_grad <- function(
  config,
  X,
  theta,
  hess_res,
  z,
  u
){
  H = hess_res$H
  y
}

#' Compute the suff stat in stergm at time t
#'
#' @param formula, formula, e.g. nw ~ edge + mutual, indicate suff stat
#' @param yt matrix, adjacency matrix for network at time t, dimension NxN
#' @param yt0 matrix, adjacency matrix for network at time t-1
#'
#' @return a vector of length p1 + p2 containing g+ and g-
gstat <- function(
  config,
  yt0,
  yt
){
  ypos <- 1*(yt | yt0)
  yneg <- 1*(yt & yt0)
  nwpos <- network(ypos)
  nwneg <- network(yneg)
  fml_form = as.formula(paste('nwpos ~',paste(config$myf_form, collapse=" + ")))
  fml_diss = as.formula(paste('nwneg ~',paste(config$myf_diss, collapse=" + ")))
  gpos <- summary_formula(fml_form)
  gneg <- summary_formula(fml_diss)
  return(c(gpos, gneg))
}



theta_update <- function(
  config,
  X,
  theta,
  z,
  u
){
  
}



theta_update_f <- function(){
  
}

matmul <- function (A, B){
  return(A %*% B)
} 