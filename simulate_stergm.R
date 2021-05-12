simulate_stergm <- function(
  object, 
  attr = NULL,
  num.nodes.K,
  order = 1,
  nsim = 100, 
  MCMC.burnin = 10000,
  MCMC.interval = 1000,
  seed = 123, 
  directed = c(TRUE, FALSE)
){
  set.seed(seed)
  #  Adjust the input according to DEFAULT if not specified
  directed = directed[1]
  # Check that input is OK
  # 1) Check that object is a formula object
  if (class(object) != "formula" & class(object) != "ergm") {
    stop("argument object must be either a formula or ergm object. Please respecify.") }
  
  K = length(num.nodes.K) # length of time series
  
  #Initialize the list of simulated networks and the statistics for each
  network.sims = vector("list", nsim)
  if (is.null(phi) == FALSE) {
    p = nrow(phi)
  } else{
    p = nrow(phicoef)
  }  # number of network statistics
  
  h.statistics = rep(list(matrix(NA, K, p)), nsim) #store sufficient statistics of all networks
  
  # Go through each time point and simulate nsims networks for each
  for (s in 1:K)
  {
    cat("Simulating networks for time", s, "\n")
    num.nodes = num.nodes.K[s]
    nets = network(num.nodes, directed = directed)
    
    # Attributes
    if (is.null(attr) == FALSE) {
      if (is.vector(attr[[s]])) {
        attr.s = vector("list", 1)
        attr.s[[1]] = attr[[s]]
        names(attr.s) = "attr1"
      } else {
        attr.s = vector("list", ncol(attr[[s]]))
        for (l in 1:ncol(attr[[s]])) {
          attr.s[[l]] = attr[[s]][, l]
        }
        names(attr.s) = paste("attr", 1:ncol(attr[[s]]), sep = "")
      }
      nets = network(nets, vertex.attr = attr.s, directed = directed)
    } else {
      nets = network(nets, directed = directed)
    }
    
    #replacing object with current network formula
    #    formula.s = ergm.update.formula(object, nets ~ ., from.new = TRUE)
    formula.s = nonsimp_update.formula(object, nets ~ ., from.new = TRUE)
    
    # Use an existing function in package 'ergm'
    if(is.null(attr)){
      sims = simulate(object = formula.s, coef = coef.s, nsim = nsim, seed = seed,
                      control = control.simulate(MCMC.burnin = MCMC.burnin,
                                                 MCMC.interval = MCMC.interval))
    } else{
      sims = simulate(object = formula.s, coef = coef.s, attr = attr.s, nsim = nsim, seed = seed,
                      control = control.simulate(MCMC.burnin = MCMC.burnin,
                                                 MCMC.interval = MCMC.interval))
    }
    if (nsim == 1)
    {
      network.sims[[1]][[s]] = as.matrix.network(sims, matrix.type = "adjacency")
      h.statistics[[1]][s,] = summary(as.formula(paste("sims ~ ", deparse(object[[3]]), sep = "")))
      if (is.null(rownames(phicoef)) == FALSE) {colnames(h.statistics[[1]]) = rownames(phicoef)}
      #      colnames(h.statistics[[1]]) = stat
    }
    
    if (nsim > 1)
    {
      for (i in 1:nsim)
      {
        network.sims[[i]][[s]] = as.matrix(sims[[i]], matrix.type = "adjacency")
        h.statistics[[i]][s, ] = as.matrix(attr(sims, "stats"))[i, ]
        if (is.null(names(coef.s)) == FALSE) {colnames(h.statistics[[i]]) = names(coef.s)}
        #        colnames(h.statistics[[i]]) = stat
      }
    }
    
  }
  return(list(Networks = network.sims, Statistics = h.statistics))
}
