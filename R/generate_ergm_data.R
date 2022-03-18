rm(list=ls())
library(tergm)
library(network)
library(igraph)
library(ndtv)
library(dplyr)
source('util.R')

#' Some candidate models
#' form_model = ~ edges + mutual + transitiveties + cyclicalties
#' diss_model = ~ offset(edges)
# To get an average duration of 5...
# coef.diss <- log(5 - 1)

# user defined parameters -------------------------------------------------
num_nodes = 30
num_timepts = 100
num_changepts = 4
# form_model= ~edges+triangle+transitiveties
# diss_model= ~edges+triangle+transitiveties
form_model= ~edges+triangle + nodematch("Sex")
diss_model= ~edges+triangle + nodematch("Sex")

phi1 = c(-1,1,2,-3,1)[1:num_changepts]
phi2 = c(-2,-1,1,-1,-2)[1:num_changepts]
phi3 = c(log(8),-log(3), -log(2), log(2), log(5))[1:num_changepts]
phi4 = c(-2, 2, -1, -1,1)[1:num_changepts]
phi5= c(1, 2, 3, -4, -2)[1:num_changepts]
phi6 = c(1, 2, 3, 3, -3)[1:num_changepts]

coefs_pos = rbind(phi1, phi2, phi5)
coefs_neg = rbind(phi3, phi4, phi6)
# coefs_pos = rbind(phi1, phi2, phi5)
# coefs_neg = rbind(phi3, phi4, phi6)
# plot_model_param(num_timepts, num_changepts, coefs_pos, coefs_neg)

# simulate data -----------------------------------------------------------
# generate initial network 
set.seed(1)
target.stats <- num_nodes
g0<-network.initialize(num_nodes, dir=FALSE)
g1<-san(g0~edges,target.stats=target.stats, verbose=TRUE)
plot(g1)
# ginit = g1
sex = c("Male", "Female")[rbinom(num_nodes, 1, 0.5) + 1]
network::set.vertex.attribute(g1,"Sex", sex)

# sapply(res_adj_list, sum) # count number of edges for each graph
# plot_nw_seq(stergm.sim, end=50)
sim <- simulate_nw_ts(
  form_model = form_model,
  diss_model = diss_model,
  coefs_pos = coefs_pos,
  coefs_neg = coefs_neg,
  num_timepts = num_timepts,
  num_changepts = num_changepts,
  num_nodes = num_nodes,
  desc = "form: edge+triangle+sex diss: edge+triangle+sex",
  # desc = "form: edge+mutual + nodematch_sex. diss: edge+mutual+nodematch_sex",
  output_path='../data/ergm/sim',
  vers = 2000, 
  exogenous=T,
  ginit = g1
)
# show heatmap of all networks. This takes a while to run
plot_adj_ts(res_adj_list = sim$nw)

# experiment --------------------------------------------------------------

samp.labs <- substr(get.node.attr(g1,"Sex"),1,1)
plot(g1, label=samp.labs, vertex.col='Sex')


plot(g1, label=samp.labs, vertex.col='Sex')
get.node.attr(g1,"Sex") %>% table()

for(i in 1:10){
  print(get.node.attr(res_adj_list[[i]], 'Sex') %>% table()  )
}


time_stable = num_timepts / num_changepts
cur_end = 0
res_adj_list = vector(mode = 'list', length = num_timepts)
# num_changepts
# time_stable
i = 1
ginit = g1
g1
for(i in 1:num_changepts){
  cat('[INFO] Simulate from ', cur_end, ' to ', cur_end + time_stable, '\n')
  stergm.sim <- simulate(
    # if the input graph already has 'net.obs.period', to continue the simulation,
    # we need to specify the time.start to be the end of 'net.obs.period'
    ginit, 
    formation=form_model,
    dissolution=diss_model,
    coef.form=coefs_pos[,i],
    coef.diss=coefs_neg[,i],
    nsim = 1, 
    time.slices = time_stable,
    time.start = cur_end # needs to be update everytime a new dynamic is generated
  )
  # newstart_nw%n%'net.obs.period' check https://github.com/statnet/tergm/blob/master/R/simulate.stergm.R for details
  
  for(t in (1 + cur_end) : (time_stable + cur_end)){
    tmpnet = network.extract(stergm.sim, at = t) %>% network() 
    res_adj_list[[t]] <- tmpnet
  }
  cur_end = cur_end + time_stable
  ginit = network.extract(stergm.sim, at = cur_end) %>% network() 
}

# experiment ends here ----------------------------------------------------



# network_list = lapply(sim$nw[1:10], network)
# res = tergm(
#   network_list~
#     Form(~edges + triangle)+
#     Diss(~edges + triangle),
#   estimate = "CMLE"
# )

res

# df = as.data.frame(sim)

# 
# sim <- simulate_nw_ts_random(
#   form_model = form_model,
#   diss_model = diss_model,
#   coefs_pos = coefs_pos,
#   coefs_neg = coefs_neg,
#   num_timepts = num_timepts,
#   changepts = c(0, 5, 15, 30, 41),
#   num_nodes = num_nodes,
#   desc = "form: edge+mutual. diss: edge+mutual",
#   vers = 101, 
#   ginit = g1
# )
# 
# df = as.data.frame(stergm.sim)
# df$duration %>% mean()
# 


# ergm with exogenous variables -------------------------------------------

data("faux.mesa.high")
library(dplyr)
faux.mesa.high %>% plot()
model4 <- ergm(
  faux.mesa.high ~ edges + nodematch("Grade") ++ gwesp(0.5, fixed = TRUE),
  verbose = TRUE
  # seed =789
)
