# install.packages('statnet')
# install.packages('network')
# install.packages('ergm')
# install.packages('sna')

library(ergm)




# VCERGM example ----------------------------------------------------------
# library(devtools)
# devtools::install_github('jihuilee/VCERGM')

library(VCERGM)

object = net ~ edges + mutual
K = 50 # number of time points
n = 30 # number of nodes / network size
num.nodes.K = rep(n, K) # Number of nodes at each time point
degree.spline = 3
interior.knot = 20
directed = TRUE # Directed network
nsim = 1 # Simulate one sequence of networks
seed = 123

# True phi(t) under VCERGM
phi1 = sin(((1:K) + 20)/5) + 1 # phi(t) for edges
phi2 = sin(((1:K) + 20)/3) * 0.6 + 0.4 # phi(t) for reciprocity
phi = rbind(phi1, phi2)
plot(1:K, phi1, type = "l", ylim = range(phi), xlab = "Time", ylab = "Phi", main = "True phi(t)")
lines(1:K, phi2, col = 2)


network = simulate_vcergm(object = object, num.nodes.K = num.nodes.K, phi = phi, 
                          nsim = nsim, seed = seed, directed = directed)
network$Statistics
