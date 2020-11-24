library("data.table")
library("graphics")
library("colorRamps")
##
library("dplyr")
library("readxl")
##
library("tidyverse")
library("nleqslv")
##
library("rstan")
library("gdata")
library("bayesplot")
##
rm(list=ls(all=T))
setwd("PATH")
#############################################################################

# logistic map
logis <- function(theta,x){
  return(1-theta*x^2)
}

# deterministic and random logistic map orbits
lorbit <- function(n,thetatrue,x0,nvar){
  orbit <- numeric(n)
  if (nvar == 0){
    orbit[1] <- logis(thetatrue, x0)
    for (i in 2:n){
      orbit[i] <- logis(thetatrue, orbit[i-1])
    }
  } else if (nvar > 0){
    orbit[1] <- logis(thetatrue, x0) + rnorm(1,0,sqrt(nvar))
    for (i in 2:n){
      orbit[i] <- logis(thetatrue, orbit[i-1]) + rnorm(1,0,sqrt(nvar))
    }
  }
  return(orbit)  
}
#############################################################################

x0 <- 0.5 # initial condition
n <- 250 # sample size
thetatrue <- 1.71 # control parameter
nvar <- 1e-03 # noise variance (if "0" --> deterministic)
orbitdet <- lorbit(n,thetatrue,x0,0) # deterministic orbit
orbitnoi <- lorbit(n,thetatrue,x0,nvar) # random orbit

par(mfrow=c(1,1),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
plot(orbitdet,type="l",xlab="n",ylab=expression(x[n]))
lines(orbitnoi,col="red")
#############################################################################
stanc("mylogis.stan") # RECONSTRUCTION
stan_data <- list(N = n, x = orbitnoi)
###########################################
fit <- stan(file = "mylogis.stan", data = stan_data, warmup = 500, iter = 5000, chains = 4, cores = 2, thin = 2)
fit
## PLOT RESULTS BASED ON THE POSTERIOR DISTRIBUTION
## PARAMETERS: CONTROL PARAMETER, INITIAL CONDITION, NOISE VARIANCE
posterior <- extract(fit)
str(posterior)
par(mfrow=c(1,3),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
hist(posterior$theta,col = "beige",freq = F,xlab="theta",main="")
abline(v=thetatrue,col="red",lwd=2)
hist(posterior$x0,col = "beige",freq = F,xlab=expression(x[0]),main="")
abline(v=x0,col="red",lwd=2)
hist(posterior$sigma,col = "beige",freq = F,xlab=expression(sigma),main="")
abline(v=sqrt(nvar),col="red",lwd=2)
##
traceplot(fit)
##
stan_dens(fit)
##
stan_hist(fit)
##
par(mfrow=c(1,1),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
plot(density(posterior$x0), main = "Initial condition posterior KDE",xlab=expression(x[0]))
points(x0,0, pch=17)

#############################################################################
stanc("mylogis_pred.stan") # RECONSTRUCTION
npred <- 8 # PREDICTION HORIZON
x0 <- 0.5
n <- 250
thetatrue <- 1.71
nvar <- 5e-04
orbitnoi <- lorbit(n+npred,thetatrue,x0,nvar)
par(mfrow=c(1,1),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
plot(orbitnoi,type="l",xlab="n",ylab=expression(x[n]))

stan_data <- list(N = n, x = orbitnoi[1:n], N_new=npred)
###########################################
fit <- stan(file = "mylogis_pred.stan", data = stan_data, warmup = 2000, iter = 5000, chains = 2, cores = 2, thin = 2)
fit
posterior <- extract(fit)
str(posterior)
## PLOT RESULTS BASED ON THE POSTERIOR DISTRIBUTION
## PARAMETERS: CONTROL PARAMETER, INITIAL CONDITION, NOISE VARIANCE &
## (npred) FUTURE UNOBSERVED OBSERVATIONS
par(mfrow=c(1,3),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
hist(posterior$theta,col = "beige",freq = F,xlab="theta",main="")
abline(v=thetatrue,col="red",lwd=2)
hist(posterior$x0,col = "beige",freq = F,xlab=expression(x[0]),main="")
abline(v=x0,col="red",lwd=2)
hist(posterior$sigma,col = "beige",freq = F,xlab=expression(sigma),main="")
abline(v=sqrt(nvar),col="red",lwd=2)
##
traceplot(fit)
##
par(mfrow=c(1,1),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
plot(density(posterior$x0), main = "Initial condition posterior KDE",xlab=expression(x[0]))
points(x0,0, pch=17)
##
x_pred <- as.matrix(fit, pars = "x_pred")
dim(x_pred)
par(mfrow=c(4,2),mar=c(3.5,3.5,1,1),mgp=c(2.4,0.8,0)) 
for (i in 1:npred){
  plot(density(x_pred[,i]), main = paste0("Prediction KDE: ",i),xlab=expression(x[i]))
  points(orbitnoi[n+i],0, pch=17)
}


################################################################
################################################################

