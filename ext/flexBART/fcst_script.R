library(invgamma)
library(stochvol)
library(factorstochvol)
library(MASS)
library(dbarts)
library(mvtnorm)

# Source several helper functions 
set.seed(123) #set a seed to allow reproducibility
source("aux_func.R")
source("estimation_functions/flexBART.R")
source("estimation_functions/fullBART.R")
source("estimation_functions/errorBART.R")

#Create a trivial dataset
Ytrn <- matrix(rnorm(250*5), 250, 5) #Y has T=250 and M=5



sl.mod <- c("mixBART") #allows for mixBART and BART. if you want to use fullBART you need to change the function call in line 34 accordingly
sl.sv.mod <- c("homo") #you can also put in "SV" or "heteroBART" for different volatility models
sl.fc.type <- c("exact")
fhorz <- 12

#Number of burned and saved draws. Set this much larger
nsave <- 500 
nburn <- 500

prior.mean <- matrix(0, ncol(Ytrn), ncol(Ytrn))
prior.mean[2,2] <- 1 #In case you would like to introduce a random prior if one of the series has a unit root. Otherwise, just comment this

quiet <- FALSE

fcst.unconditional <- flexBART(Ytrn, nburn=nburn, nsave = nsave, thinfac = 1, prior="HS", prior.sig = c(3, 0.9), model = sl.mod, sv = sl.sv.mod, fc.approx=sl.fc.type, fhorz = fhorz, quiet = quiet, pr.mean = prior.mean)

plot(density(fcst.unconditional$fcst[,"fhorz1",1])) #Plots the one-step-ahead predictive distribution for the firs series in Ytrn

