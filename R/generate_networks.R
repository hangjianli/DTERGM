library(tergm)
library(network)
library(ndtv)

#' Some candidate models
#' form_model = ~ edges + mutual + transitiveties + cyclicalties
#' diss_model = ~ offset(edges)

num_nodes = 10
num_timepts = 50
num_changepts = 5

# Construct an initial network with n nodes and n edges
target.stats <- num_nodes
g0<-network.initialize(num_nodes, dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)
plot(g1)

# To get an average duration of 5...
coef.diss <- log(5 - 1)

form_model=~edges+mutual
diss_model=~edges

phi1 = c(-1,1,2,-3,1)
phi2 = c(-2,1,2,-1,1)
coefs = rbind(phi1, phi2)
phi1_ = rep(phi1, each = num_timepts / num_changepts)
phi2_ = rep(phi2, each = num_timepts / num_changepts)
coefs_ = rbind(phi1_, phi2_)
plot(1:num_timepts, phi1_, type = "l", 
     ylim = range(coefs_), xlab = "Time", ylab = "Phi", main = "True phi(t)")
lines(1:num_timepts, phi2_, col = 2)

time_stable = num_timepts / num_changepts

cur_end = 0
res_adj_list = vector(mode = 'list', length = num_timepts)
g0<-network.initialize(num_nodes, dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)
for(i in 1:num_changepts){
  cat('[INFO] Simulate from ', cur_end, ' to ', cur_end + time_stable, '\n')
  stergm.sim <- simulate(
    g1, 
    formation=form_model,
    dissolution=diss_model,
    coef.form=coefs[,i],
    coef.diss=coef.diss,
    nsim = 1, 
    time.slices = time_stable,
    time.start = cur_end # needs to be update everytime a new dynamic is generated
  )
  # newstart_nw%n%'net.obs.period' check https://github.com/statnet/tergm/blob/master/R/simulate.stergm.R for details
  
  for(t in (1 + cur_end) : (time_stable + cur_end)){
    tmpnet = network.extract(stergm.sim, at = t) %>% network() 
    res_adj_list[[t]] <- tmpnet[, ]
  }
  cur_end = cur_end + time_stable
  g1 = network.extract(stergm.sim, at = cur_end) %>% network() 
}

for(i in 1:num_timepts){
  image(res_adj_list[[i]],
        xlab = paste0('Time at ', i))
  Sys.sleep(1)
}

saveRDS(res_adj_list, 'data/simulated_nw_adj.rds')
sapply(res_adj_list, sum)

plot_nw_seq(stergm.sim, end=50)

plot_nw_seq <- function(nw_dynamic, start=0, end=10, nw_stats="~edges"){
  #' interactive plots of generated nw time series
  slice.par = list(start = start, 
                   end = end, 
                   interval = 1, 
                   aggregate.dur = 1, 
                   rule = "any")
  
  compute.animation(nw_dynamic, slice.par = slice.par)
  render.par = list(tween.frames = 10,
                    show.time = T,
                    show.stats = nw_stats)
  
  elabels<-lapply(get.edge.activity(stergm.sim),
                  function(spl){
                    paste("(",spl[,1],"-",spl[,2],")",sep='')
                  })
  
  plot.par = list(edge.col = "darkgray",
                  displaylabels = T,
                  label.cex = .8,
                  label.col = "blue",
                  edge.label = elabels
                  # vertex.cex = wealthsize)
                  )
  render.d3movie(nw_dynamic,
                 render.par = render.par,
                 plot.par = plot.par,
                 output.mode = 'htmlWidget')
  
}
