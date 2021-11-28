# install.packages("tergm")
# install.packages('ergm')
library(ergm)
library(tergm)
library(dplyr)

coef(samplk2)

set.seed(1)
n = 10
t = 5
nw_list = vector(mode = 'list', t)
for(i in 1:t){
  nw_list[[i]] = matrix(rbinom(n=n^2, size=1, prob = 0.5), n, n)
}
sim1nw = lapply(nw_list, network)
sim1NS = NetSeries(sim1nw)
res <- ergm.bridge.dindstart.llk(
  object = sim1NS ~ Form(~edges+mutual)+Diss(~edges+mutual),
  coef = c(-2,-1,1,1),
  verbose = T
)

remotes::install_github("statnet/tergm")
# ergm.bridge.dindstart.llk<-function(object, response=NULL, constraints=~., 
#                                     coef, obs.constraints=~.-observed, 
#                                     target.stats=NULL, dind=NULL, coef.dind=NULL, 
#                                     basis=ergm.getnetwork(object), ..., 
#                                     llkonly=TRUE, control=control.ergm.bridge(), verbose=FALSE){

# sim101_n20_t40_p2 is a list of adjacency matrices
sim1nw = lapply(sim101_n20_t40_p2, network)
sim1NS = NetSeries(sim1nw[1:5])
ergm.bridge.dindstart.llk(
  object = sim1NS ~ Form(~edges+mutual)+Diss(~edges+mutual),
  coef = c(1,1,1,1),
  verbose = T
)
object = sim1NS ~ Form(~edges+mutual)+Diss(~edges+mutual)
basis = ergm.getnetwork(object)

nw = ergm.getnetwork(sim1NS ~ edges)
ergm_preprocess_response(nw, NULL)

ergm.bridge.llr(sim1NS ~ Form(~edges+mutual)+Diss(~edges+mutual),
                from=c(1,1),
                to = c(1,1),
                verbose = T)



sim_settings <- simulate(object, coef=c(1,1), nsm=1, 
                         reference=~Bernoulli,
                         # constraints=list(constraints, obs.constraints),
                         observational=FALSE, output="ergm_state",
                         verbose=max(verbose-1,0),
                         ergm.getnetwork(object))


logit<-function(p)log(p/(1-p))
coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)
# Construct a network with 20 nodes and 20 edges
n<-20
target.stats<-edges<-20
g0<-network.initialize(n,dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

S<-10
# To get an average duration of 10...
duration<-10
coef.diss<-logit(1-1/duration)
# To get an average of 20 edges...
dyads<-network.dyadcount(g1)
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)
# ... coefficients.
print(coef.form)
print(coef.diss)
# Simulate a networkDynamic

dynsim<-simulate(g1,formation=~edges,dissolution=~edges,
                 coef.form=coef.form,coef.diss=coef.diss,
                 time.slices=S,verbose=TRUE)