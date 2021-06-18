X = readRDS('data/edge_mutual__mutual.rds')
library(tergm)
library(network)
library(ndtv)
source('R/util.R')

config = list()
config$admm_alpha = 100
config$tslen = length(X)
config$num_nodes = dim(X[[1]])[1]
config$E =config$num_nodes^2
config$myf_form = c('edges', 'mutual')
config$myf_diss = c('edges')

p_pos = 2
p_neg = 1
p = p_pos + p_neg



theta_pos <- 
theta_neg <- 
theta_agg <- numeric(length = )
# E =  choose(num_nodes, 2) # undirected
E =  config$num_nodes^2 # directed A_{ii} = 0 for all i and the corresponding change stats are zero

# test get_suff_stats -----------------------------------------------------



theta = rep(rnorm(1), times=p_pos+p_neg) # need to be feasible???


H = matrix(0, nrow = num_timepts * E, ncol = num_timepts * p)

delta = matrix(0, p, E) # p x E, filled column first
# delta = (d_11, d_21, d_31, ... d_1n, d_2n, ..d_nn)

delta_ijt = 

yt0 = X[[1]]
yt = X[[2]]


# calculate T-1 g ---------------------------------------------------------
#test gtstat
gg = gstat(config, yt0 = yt0, yt = yt1)
gg
#test gstat_delta
res = gstat_delta(config, yt0 = yt0, yt = yt1)
res[,1:4]
