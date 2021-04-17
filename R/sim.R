# install.packages('statnet')
# install.packages('network')
# install.packages('ergm')
# install.packages('sna')
# library(devtools)
# devtools::install_github('jihuilee/VCERGM')

# VCERGM example ----------------------------------------------------------
library(ergm)
library(VCERGM)
library(network)
library(statnet)
library(tergm)
#' @importFrom network network
#' 
# Two network statistics (edges and reciprocity) in the model
# Use the names of statistics specified in R package 'ergm'
# Use the same strategy to specify the model 'object' as in R package 'ergm'
object = net ~ edges + mutual + transitiveties + cyclicalties
K = 120 # number of time points
n = 30 # number of nodes / network size
num.nodes.K = rep(n, K) # Number of nodes at each time point
directed = T # Directed network
nsim = 1 # Simulate one sequence of networks
seed = 123

# True phi(t) under VCERGM
# phi1 = sin(((1:K) + 20)/5) + 1 # phi(t) for edges
# phi2 = sin(((1:K) + 20)/3) * 0.6 + 0.4 # phi(t) for reciprocity
phi1 = rep(c(-0.4,0.2,0.3), each = K/3)
phi2 = rep(c(0.3,-0.1,0.2), each = K/3)
phi3 = rep(c(0.3,-0.1,0.2), each = K/3)
phi4 = rep(c(0.3,-0.1,0.2), each = K/3)
# phi3 = rep(c(0.3,-0.2,0.1),each=K/3)
phi = rbind(phi1, phi2, phi3, phi4)
plot(1:K, phi1, type = "l", ylim = range(phi), xlab = "Time", ylab = "Phi", main = "True phi(t)")
lines(1:K, phi2, col = 2)
# lines(1:K,phi3, col=3)

network = simulate_vcergm(object = object, num.nodes.K = num.nodes.K, phi = phi, 
                          nsim = nsim, seed = 1, directed = directed)
# for(i in seq(5,K,10)){
#   image(network$Networks[[1]][[i]])
#   cat(i,'\n')
# }
network$Statistics
networks = network$Networks

vcergmest = mple(object = object, 
                 network = networks,
                 directed = directed, 
                 constant = FALSE)

