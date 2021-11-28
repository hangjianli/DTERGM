#!/usr/bin/env Rscript
suppressMessages({
  library(ergm)
  library(network, verbose = F)
  library(optparse)
})


option_list = list(
  make_option(c("-s", "--nw_size"), type="integer", default=10, help="Size of the network."),
  make_option(c("-N", "--number_of_nw"), type="integer", default=1000, help="Number of MCMC networks to generate"),
  make_option(c("-o", "--outfile"), type="character", help="Output file name."),
  make_option(c("-e", "--experiment_dir"), type="character", help="experiment directory"),
  make_option(c("-s", "--model_statistics"), type="character", help="ergm terms")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



# Command:
# Rscript gen_mcmc_samples.R -o out -e "../output/exp2"

terms = stringr::str_split(s, " ", simplify = T)
nw_init = network(opt$nw_size, directed = F)
gsim = simulate(
  as.formula(paste('~ ',paste(terms, collapse = " + "))), 
  nsim = opt$N,
  coef=c(0, 0), 
  basis=nw_init)


adj = matrix(0, 10, 10)
for(i in 1:5){
  summary(as.formula(paste("gsim[[i]] ~ ",paste(terms, collapse = " + ") ))) %>% print()
}
