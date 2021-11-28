rho_a = 0.9
num_time = 150
n = 200
library(dplyr)
library("Matrix")

data  = matrix(0,n^2, num_time)
v =  c( floor(num_time/3)+1,2*floor(num_time/3)+1 )

# t = 1
E =  list()
for(t in 1:num_time)
{
  ###
  if(t==1 ||  t== v[2]+1 )
  {
    P =  matrix(0.3,n,n)
    P[1:floor(n/4), 1:floor(n/4)] = 0.5
    P[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = 0.5
    P[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = 0.5
    P[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = 0.5
    diag(P) = 0

    # heatmap(P)
    # P    
    A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),P),n,n)
    aux = drop(A)
    aux[lower.tri(aux)] = aux[lower.tri(aux)]  
    diag(aux) = 0
    
    data[,t] = drop(matrix(A,n^2,1 ))
    
    # temp = svd(A)
    # xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) ) 
  }
  if( (t > 1      && t <=  v[1])  ||  ( t >v[2]+1 ) )
  {
    aux1 = P +  (1-P)*rho_a
    aux2 = P*(1-rho_a)
    
    aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
    aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
    A =  aux1*A + aux2*(1-A)
    
    aux = drop(A)
    aux[lower.tri(aux)] = aux[lower.tri(aux)]  
    diag(aux) = 0
    
    data[,t] = drop(matrix(A,n^2,1 ))
    
    # temp = svd(A)
    # xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
  }
  if(t ==  v[1]+1)
  {
    Q =  matrix(0.2,n,n)
    Q[1:floor(n/4), 1:floor(n/4)] = 0.45
    Q[(1+floor(n/4)):(2*floor(n/4)),(1+floor(n/4)):(2*floor(n/4)) ] = 0.45
    Q[(1+2*floor(n/4)):(3*floor(n/4)),(1+2*floor(n/4)):(3*floor(n/4)) ] = 0.45
    Q[(1+3*floor(n/4)):n,(1+3*floor(n/4)):n ] = 0.45
    diag(Q) = 0
    
    A = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),Q),n,n)
    aux = drop(A)
    aux[lower.tri(aux)] = aux[lower.tri(aux)]  
    diag(aux) = 0
    data[,t] = drop(matrix(A,n^2,1 ))
    
    # temp = svd(A)
    # xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) ) 
  }
  if(t > v[1]+1 && t <= v[2] )
  {
    aux1 = Q +  (1-Q)*rho_a
    aux2 = Q*(1-rho_a)
    
    aux1 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux1),n,n)
    aux2 = matrix(rbinom(matrix(1,n,n),matrix(1,n,n),aux2),n,n)
    A =  aux1*A + aux2*(1-A)
    
    aux = drop(A)
    aux[lower.tri(aux)] = aux[lower.tri(aux)]  
    diag(aux) = 0
    
    data[,t] = drop(matrix(A,n^2,1 ))
    
    # temp = svd(A)
    # xhat[t,,] =  temp$u[,1:d] %*%  diag( sqrt(temp$d[1:d]) )
  }
  ##################
  
  E[[t]] = A
}####for t
# E # not used

df = vector(mode = "list", length = num_time)
for(t in 1:num_time){
  df[[t]] = matrix(data[,t], n, n)
}
saveRDS(df, paste0('../data/scenario1_4_rho09n300.rds'))
# as(df[[29]] , "CsparseMatrix") %>% image()
