tryCatch(
  {
    current_path = rstudioapi::getActiveDocumentContext()$path
    setwd(dirname(current_path ))
  }, 
  error=function(cond){
    if (identical(cond, "RStudio not running")){
      this.dir <- dirname(parent.frame(2)$ofile)
      setwd(this.dir)
    }
  })

library(Barycenter)
library(transport)
library(plot.matrix)

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

# calculate tplans 
tplan_ab <- transport::transport(transport::pgrid(a), transport::pgrid(b), p = 2)
tplan_ac <- transport::transport(transport::pgrid(a), transport::pgrid(c), p = 2)

# convert tplans to full matrices
tplan_ab_mat <- matrix(rep(0, n^4), n^2, n^2)
for (i in 1:nrow(tplan_ab)){
  tplan_ab_mat[tplan_ab$from[i], tplan_ab$to[i]] <- tplan_ab$mass[i]
}
tplan_ac_mat <- matrix(rep(0, n^4), n^2, n^2)
for (i in 1:nrow(tplan_ac)){
  tplan_ac_mat[tplan_ac$from[i], tplan_ac$to[i]] <- tplan_ac$mass[i]
}

# calculate sinkhorn plans
# set up grid and cost matrix 
grid <- expand.grid(y=seq(0.05, 0.95, 0.1), x=seq(0.05, 0.95, 0.1))
costm <- as.matrix(dist(grid, diag=TRUE, upper=TRUE))
# do actual calculations
results_ab <- Barycenter::Greenkhorn(a, b, costm)
results_ac <- Barycenter::Greenkhorn(a, c, costm)


##### do the plotting #################################################
pdf("../results/sinhorn_image1_image2.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), pty="s") # adapt margins
plot(results_ab$Transportplan, main="", xlab = "from", ylab = "to")
dev.off()

pdf("../results/tplan_image1_image2.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), pty="s") # adapt margins
plot(tplan_ab_mat, main="", xlab = "from", ylab = "to")
dev.off()

pdf("../results/sinhorn_image1_image3.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), pty="s") # adapt margins
plot(results_ac$Transportplan, main="", xlab = "from", ylab = "to")
dev.off()

pdf("../results/tplan_image1_image3.pdf")
par(mar=c(5.1, 4.1, 4.1, 4.1), pty="s") # adapt margins
plot(tplan_ac_mat, main="", xlab = "from", ylab = "to")
dev.off()


