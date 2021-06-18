
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


genlasso()
