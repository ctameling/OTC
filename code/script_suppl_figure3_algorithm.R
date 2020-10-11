### compute transport plans with different algorithms
library(transport)
library(ggplot2)
library(collections)
RNGversion("3.5.3")  # change to old sampling default
set.seed(13)
dim <- 64
pxsize <- 15

# set up random images
a <- rbinom(dim^2, 1, p = 0.6)
a <- matrix(a, dim, dim)/ sum(a)

b <- rbinom(dim^2, 1, p = 0.4)
b <- matrix(b, dim, dim)/sum(b)

# algorithms
algos <- c('aha', 'shielding', 'revsimplex', 'primaldual', 'shortsimplex')
tplans <- collections::dict()

# calculate transportplans with differenrt algorithms
for (alg in algos){
  tplan <- transport::transport(transport::pgrid(a), transport::pgrid(b), p=2, method = alg)
  coord <- as.data.frame(arrayInd(tplan$from, c(dim, dim)) - arrayInd(tplan$to, c(dim,dim)))
  tplan$dist <- sqrt(pxsize^2*coord$V1^2 + pxsize^2*coord$V2^2)
  tplans$set(alg, tplan)
}

# evaluate OTC on transport plans from different algorithms
results <- data.frame(class = character(0), t = numeric(0), OTC = numeric(0))

for (alg in algos){
  print(alg)
  tplan <- tplans$get(alg)
  for (t in seq(0,pxsize*dim*sqrt(2)/10)){
      results <- rbind(results, data.frame(class = alg, t = t, OTC = 1- sum(tplan[tplan$dist > t,]$mass)))  
  }
}

# plot results
ggplot(data = results, aes(x = t, y = OTC, color = class)) + 
  geom_line() + 
  xlab('threshold in nm') + 
  ylab('Optimal Transport Colocalization') + 
  labs(color = 'Algorithm')