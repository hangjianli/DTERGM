sigmoid <- function(x){
  return(1 / (1 + exp(-x)))
}

plot_adj_ts <- function(
  res_adj_list
){
  for(i in 1:length(res_adj_list)){
    image(res_adj_list[[i]],
          xlab = paste0('Time at ', i))
    Sys.sleep(0.5)
  }
}


plot_model_param <- function(
  num_timepts,
  num_changepts,
  coefs_pos,
  coefs_neg
){
  phi_pos_1 = rep(coefs_pos[1, ], each = num_timepts / num_changepts)
  phi_pos_2 = rep(coefs_pos[2, ], each = num_timepts / num_changepts)
  phi_neg_1 = rep(coefs_neg[1, ], each = num_timepts / num_changepts)
  # phi_neg_2 = rep(coefs_neg[2, ], each = num_timepts / num_changepts)
  
  coefs_ = rbind(phi_pos_1, phi_pos_2)
  plot(
    1:num_timepts, phi_pos_1, type = "l", 
    ylim = range(coefs_), 
    xlab = "Time",
    ylab = "Phi",
    main = "Change points of parameters",
    col = 1)
  lines(1:num_timepts, phi_pos_2, col = 2)
  lines(1:num_timepts, phi_neg_1, col = 3)
  # lines(1:num_timepts, phi_neg_2, col = 4)
  
  legend(
    "bottomleft",
    legend = c("pos 1", "pos 2",
               "neg 1",
               # "neg 2"
               ),
    col = c(1,2,3,
            # 4
            ),
    lty = 1,
    cex = 0.7,
    bty = "n"
  )
}


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


#' Simulate network sequence with changing parameters
#'
#' @param form_model 
#' @param diss_model 
#' @param coefs_pos 
#' @param coefs_neg 
#' @param num_timepts 
#' @param num_changepts 
#' @param num_nodes 
#' @param g_init 
#'
#' @return
#' @export
#'
#' @examples
simulate_nw_ts <- function(
  form_model,
  diss_model,
  coefs_pos,
  coefs_neg,
  num_timepts,
  num_changepts,
  desc = 'form:[], diss: []',
  vers=1,
  num_nodes,
  ginit){
  time_stable = num_timepts / num_changepts
  cur_end = 0
  res_adj_list = vector(mode = 'list', length = num_timepts)
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
      res_adj_list[[t]] <- tmpnet[, ]
    }
    cur_end = cur_end + time_stable
    ginit = network.extract(stergm.sim, at = cur_end) %>% network() 
  }
  sim = list()
  sim$nw = res_adj_list
  sim$coefs_pos = coefs_pos
  sim$coefs_neg = coefs_neg
  sim$desc = desc
  saveRDS(sim, paste0('../data/sim', vers, '.rds'))
  return(sim)
}

simulate_nw_ts_random <- function(
  form_model,
  diss_model,
  coefs_pos,
  coefs_neg,
  num_timepts,
  changepts,
  desc = 'form:[], diss: []',
  vers=1,
  num_nodes,
  ginit){
  res_adj_list = vector(mode = 'list', length = num_timepts)
  for(i in 1:num_changepts){
    cat('[INFO] Simulate from ', changepts[i], ' to ', changepts[i+1]-1, '\n')
    stergm.sim <- simulate(
      # if the input graph already has 'net.obs.period', to continue the simulation,
      # we need to specify the time.start to be the end of 'net.obs.period'
      ginit, 
      formation=form_model,
      dissolution=diss_model,
      coef.form=coefs_pos[,i],
      coef.diss=coefs_neg[,i],
      nsim = 1, 
      time.slices = changepts[i+1] - changepts[i],
      time.start = changepts[i] # needs to be update everytime a new dynamic is generated
    )
    # newstart_nw%n%'net.obs.period' check https://github.com/statnet/tergm/blob/master/R/simulate.stergm.R for details
    
    for(t in  max(2, changepts[i]+1) : changepts[i+1]){
      tmpnet = network.extract(stergm.sim, at = t) %>% network() 
      res_adj_list[[t-1]] <- tmpnet[, ]
    }
    ginit = network.extract(stergm.sim, at = changepts[i+1]) %>% network() 
  }
  sim = list()
  sim$nw = res_adj_list
  sim$coefs_pos = coefs_pos
  sim$coefs_neg = coefs_neg
  sim$desc = desc
  sim$changepts = changepts
  saveRDS(sim, paste0('../data/sim', vers, '.rds'))
  return(sim)
}



simulate_random_graphs <- function(n, nsim=1000){
  nw_init = network(n)
  gsim = simulate(~edges + mutual, nsim = nsim, 
                  coef=c(0, 0), basis=nw_init)
  res_adj_list = vector(mode = 'list', length = nsim)
  for(i in 1:nsim){
    res_adj_list[[t-1]] <- gsim[[i]][,]
  }
  
}
