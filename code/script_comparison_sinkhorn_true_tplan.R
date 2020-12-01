setwd("/home/ctamelin/Dokumente/2-color-STED/Code")
library(Barycenter)
library(transport)

# set up a, b and c
n <- 10
a <- b <- c <- matrix(rep(0, n^2), n, n)
a[8, 2:8] <- 1
a[4:7,8] <- 1
a <- a/sum(a)

b[7, 1:7] <- 1
b[3:6,7] <- 1
b <- b/sum(b)

c[10, 2:8] <- 1
c[4:7,9] <- 1
c <- c/sum(c)

tplan_ab <- transport::transport(transport::pgrid(a), transport::pgrid(b), p = 2)
tplan_ab_mat <- matrix(rep(0, n^4), n^2, n^2)
for (i in 1:nrow(tplan_ab)){
  tplan_ab_mat[tplan_ab$from[i], tplan_ab$to[i]] <- tplan_ab$mass[i]
}
tplan_ac <- transport::transport(transport::pgrid(a), transport::pgrid(c), p = 2)
tplan_ac_mat <- matrix(rep(0, n^4), n^2, n^2)
for (i in 1:nrow(tplan_ac)){
  tplan_ac_mat[tplan_ac$from[i], tplan_ac$to[i]] <- tplan_ac$mass[i]
}
# grid 
grid <- expand.grid(y=seq(0.05, 0.95, 0.1), x=seq(0.05, 0.95, 0.1))
costm <- as.matrix(dist(grid, diag=TRUE, upper=TRUE))
results_ab <- Barycenter::Greenkhorn(a, b, costm)
results_ac <- Barycenter::Greenkhorn(a, c, costm)
image(results_ab$Transportplan, axes=FALSE, asp=1)
image(tplan_ab_mat, axes=FALSE, asp=1)
image(results_ac$Transportplan, axes=FALSE, asp=1)
image(tplan_ac_mat, axes=FALSE, asp=1)
