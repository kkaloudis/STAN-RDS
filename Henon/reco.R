library("data.table")
library("plot3D")
library("graphics")
library("ggplot2")
library("gridExtra")
library("RColorBrewer")
library("colorRamps")
library("tidyverse")
library("nleqslv")
library("dplyr")
library("utils")
library("zoo")
library("pracma")
library("rgl")
library("knitr")
library("rglwidget")
library("mvtnorm")
library("MASS")
##
# library("matrixStats")
##
library("rstan")
library("gdata")
library("bayesplot")
##
rm(list=ls(all=T))
setwd("/home/kkaloudis/Documents/Mathematics/PostDoc Napoli/Literature/Misc/STAN/my henon")
#############################################################################

henon <- function(theta,x){
  return(c(theta[1] + theta[2] * x[1] + theta[3] * x[1]^2 + x[2],theta[4] * x[1]))
}

lorbit <- function(n,thetatrue,x0,Sigma){
  orbit <- zeros(n,2)
  if (sum(diag(Sigma)) == 0){
    orbit[1,] <- henon(thetatrue, x0)
    for (i in 2:n){
      orbit[i,] <- henon(thetatrue, orbit[i-1,])
    }
  } else if (sum(diag(Sigma)) > 0){
    orbit[1,] <- henon(thetatrue, x0) + rmvnorm(1,zeros(2,1),Sigma)
    for (i in 2:n){
      orbit[i,] <- henon(thetatrue, orbit[i-1,]) + rmvnorm(1,zeros(2,1),Sigma)
    }
  }
  return(orbit)  
}
#############################################################################

x0 <- c(-0.7, 0.4)
n <- 350
thetatrue <- c(1.38, 0, -1, 0.27)
# thetatrue <- c(1., 0., -1.31, 0.23)
var1 <- 1e-04 # x-coordinate variance
var2 <- 1e-04 # y-coordinate variance
rho <- 0.25 # correlation
tau_scales <- c(sqrt(var1),sqrt(var2))
realOmega <- matrix(c(1,rho,rho,1),nrow=2)
## Construct covariance matrix
Sigmatrue <- diag(tau_scales) * realOmega * diag(tau_scales)
## Generate (deterministic or) random data
orbitdet <- lorbit(n,thetatrue,x0,zeros(2))
orbitnoi <- lorbit(n,thetatrue,x0,Sigmatrue)
##########################################################
par(mfrow=c(1,1),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
plot(orbitdet,type="p",xlab=expression(x[n]),ylab=expression(y[n]),pch=20)
points(orbitnoi,col="pink",pch=20)
legend("topright",c("deterministic","noisy"),col=c("black","pink"),pch=c(20,20))
#############################################################################
stanc("myhenon.stan")
stan_data <- list(N = n, x = orbitnoi, nu0=2, Sigma0 = matrix(c(1e03, -1, -1, 1e03),nrow=2))
###########################################
fit <- stan(file = "myhenon.stan", data = stan_data, warmup = 500, iter = 2500, chains = 2, cores = 2, thin = 2)
fit
####
posterior <- extract(fit)
str(posterior)
# Output : parameters
par(mfrow=c(2,2),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
hist(posterior$theta[,1],col = "beige",freq = F,xlab=expression(theta[1]),main="")
abline(v=thetatrue[1],col="red",lwd=2)
hist(posterior$theta[,2],col = "beige",freq = F,xlab=expression(theta[2]),main="")
abline(v=thetatrue[2],col="red",lwd=2)
hist(posterior$theta[,3],col = "beige",freq = F,xlab=expression(theta[3]),main="")
abline(v=thetatrue[3],col="red",lwd=2)
hist(posterior$theta[,4],col = "beige",freq = F,xlab=expression(theta[4]),main="")
abline(v=thetatrue[4],col="red",lwd=2)
#
# hist(posterior$x0[,1],col = "beige",freq = F,xlab=expression(x[0]),main="")
# abline(v=x0[1],col="red",lwd=2)
# hist(posterior$x0[,2],col = "beige",freq = F,xlab=expression(y[0]),main="")
# abline(v=x0[2],col="red",lwd=2)
################################################################################
# Output : Initial conditions
# Calculate kernel density estimate
bivn.kde <- kde2d(posterior$x0[,1], posterior$x0[,2], n = 150)   # from MASS package
# Contour plot overlayed on heat map image of results
par(mfrow=c(1,1))
image(bivn.kde,col = jet.col(n=150),xlab = expression(x[0]),ylab = expression(y[0]))       # from base graphics package
points(x0[1],x0[2],col="white",pch="*",cex=2)
################################################################################
# Output : Noise covariance matrix (check standard deviations for the 2 coordinates and correlation)
par(mfrow=c(1,3),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
hist(posterior$tau[,1],col = "beige",freq = F,xlab=expression(sigma[x]),main="")
abline(v=1/sqrt(var1),col="red",lwd=2)
hist(posterior$tau[,2],col = "beige",freq = F,xlab=expression(sigma[y]),main="")
abline(v=1/sqrt(var2),col="red",lwd=2)
hist(posterior$Omega[,1,2],col = "beige",freq = F,xlab=expression(rho),main="")
abline(v=rho,col="red",lwd=2)
################################################################################


