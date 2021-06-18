rm(list=ls())
library(tergm)
library(network)
library(ndtv)
source('R/util.R')

#' Some candidate models
#' form_model = ~ edges + mutual + transitiveties + cyclicalties
#' diss_model = ~ offset(edges)
# To get an average duration of 5...
# coef.diss <- log(5 - 1)

# user defined parameters -------------------------------------------------
num_nodes = 50
num_timepts = 100
num_changepts = 5
form_model=~edges+mutual
diss_model=~edges+mutual

phi1 = c(-1,1,2,-3,1)
phi2 = c(-2,2,1,-1,3)
phi3 = c(log(8),log(3), -log(2), log(2), log(10))
phi4 = c(-2,2,-1,-1,3)
coefs_pos = rbind(phi1, phi2)
coefs_neg = rbind(phi3, phi4)
# plot_model_param(num_timepts, num_changepts, coefs_pos, coefs_neg)

# simulate data -----------------------------------------------------------
# generate initial network 
target.stats <- num_nodes
g0<-network.initialize(num_nodes, dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats, verbose=TRUE)
plot(g1)
ginit = g1
# sapply(res_adj_list, sum) # count number of edges for each graph
# plot_nw_seq(stergm.sim, end=50)
res <- simulate_nw_ts(
  form_model = form_model,
  diss_model = diss_model,
  coefs_pos = coefs_pos,
  coefs_neg = coefs_neg,
  num_timepts = num_timepts,
  num_changepts = num_changepts,
  num_nodes = num_nodes,
  desc = "form: edge+mutual. diss: edge+mutual",
  vers = 3, 
  ginit = g1
)

df = as.data.frame(stergm.sim)
df$duration %>% mean()


