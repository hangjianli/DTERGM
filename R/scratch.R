
attr.s = vector("list", 1); attr.s[[1]] = attr[[s]]; names(attr.s) = "attr1"
nets = network(networks[[1]], vertex.attr = attr.s, directed = F)

temp = ergmMPLE(formula.s, output = "matrix")

temp$predictor %>% dim

nets
data(faux.mesa.high)

mynw = network(X[[1]])
mynw
# faux.mesa.high
formula <- mynw ~ edges + triangle

# + nodematch("Sex") + nodefactor("Grade")
mplesetup <- ergmMPLE(formula, output='matrix')
mplesetup$predictor[,,'triangle'] %>%  as.numeric() %>% unique()
mplesetup$weights
mplesetup$response
mplesetup$weights


a = summary(faux.mesa.high)
faux.mesa.high[, ] %>% sum()



# fused lasso -------------------------------------------------------------
library(genlasso)
set.seed(1)
n = 100
p = 10


set.seed(1)
n = 10
t = 5
nw_list = vector(mode = 'list', t)
for(i in 1:t){
  nw_list[[i]] = matrix(rbinom(n=n^2, size=1, prob = 0.5), n, n)
}
sim1nw = lapply(nw_list, network)
sim1NS = NetSeries(sim1nw)
res = ergm.bridge.dindstart.llk(
  object = sim1NS ~ Form(~edges+mutual)+Diss(~edges+mutual),
  coef = c(1,1,1,1),
  verbose = T,
  llkonly = F
)

res$llk.dind
res$llr
res$llk


res2 = ergm.bridge.dindstart.llk(
  object = sim1NS ~ Form(~edges)+Diss(~edges),
  coef = c(0,0),
  verbose = T,
  llkonly = T,
)

sim1nw = lapply(scenario1_2, network)
sim1NS = NetSeries(sim1nw)
res = ergm.bridge.dindstart.llk(
  object = sim1NS ~ Form(~edges+mutual)+Diss(~edges+mutual),
  coef = c(-1,-2,3,1),
  verbose = T,
  llkonly = F
)
res
