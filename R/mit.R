

rm(list = ls())


## Loading data
load("../data/phonedata.RData")

## data goes from 14 Sept 2004 to 5 May 2005
dim(mydata)
#data[,,1] is the matrix correpsonding to the first  four hour block on day 1
#data[,,2] is the matrix correpsonding to the second four hour block on day 1
# ...
#data[,,7] is the matrix correpsonding to the first  four hour block on day 2


library(R.matlab)

p=96 ## number of students/faculty

data=arra(0,nrow=96*96,ncol=232) ## 232 is number of days

data = array(0, dim = (c(232, 96, 96)))

for ( i in 1: 232){
  ##  all the data in the same are aggregated together
  data[i, , ]= mydata[,,((i-1)*6)+1] + 
    mydata[,,((i-1)*6)+2]+ 
    mydata[,,((i-1)*6)+3]+
    mydata[,,((i-1)*6)+4]+
    mydata[,,((i-1)*6)+5] +
    mydata[,,((i-1)*6)+6]
}
data[which(data>1)]=1
data %>% dim
saveRDS(data, '../data/mit.rds')
data[,1] %>% matrix(nrow = 96) %>% dim

data <- readRDS('../data/mit.rds')

# apply tergm -------------------------------------------------------------

network_list = apply(data, 1, network)
res = tergm(
  network_list~
    Form(~edges + mutual)+
    Diss(~edges + mutual),
  estimate = "CMLE"
)



yt0 =data[1,,] 
testnw = network(yt0, directed = F)
a = ergm(testnw ~ edges + triangles)

logLik(a)

a %>% summary()

data(florentine)
m1 = ergm(flomarriage~edges)
m1 %>% summary()

mynw <- network(matrix(0,7,7),dir=FALSE)
system.time(a <- ergm.allstats(mynw~edges+triangles))
cbind(a$weights, a$statmat)
