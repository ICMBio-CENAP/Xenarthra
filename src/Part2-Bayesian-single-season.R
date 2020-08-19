# Part2-Single-season Bayesian occupancy model
# Elildo Carvalho Jr @ ICMBio/CENAP 2020-08-17
# code based on templates from the book "Bayesian population analysis using WinBUGS - a hierarchical perspective" by Marc Kéry & Michael Schaub (2012, Academic Press)
# and on templates from the JAGS translation available at:
# https://www.vogelwarte.ch/de/projekte/publikationen/bpa/code-for-running-bpa-using-jags

#----- 1 - Load libraries-----
#library(dplyr)
#library(lubridate)
library(here)
library(R2jags)
library(rjags)
library(ggplot2)


#----- 2 - Source files-----
#source(here("bin", "ahumada_codes.R"))
#source(here("bin", "f-matrix-creator-experimental-probably-ok-but-need-check.R"))
#source(here("src", "fix_species_names.R")) # fix some names and remove "false" species


#----- 3 - Read and prepare data -----
Mtridactyla <- readRDS(here("data", "Mtridactyla.rds"))
y <- Mtridactyla[,2:16]
str(y)
SiteCovs <- Mtridactyla[,17:22]

names(SiteCovs)
distWater <- SiteCovs[,2]
slope <- SiteCovs[,3]
elevation <- SiteCovs[,4]
treeDensity <- SiteCovs[,6]

# Standardize covariates
mean.slope <- mean(slope, na.rm = TRUE)
sd.slope <- sd(slope[!is.na(slope)])
slope <- (slope-mean.slope)/sd.slope     # Standardise slope
slope[is.na(slope)] <- 0               # Impute zeroes (means)

mean.distWater <- mean(distWater, na.rm = TRUE)
sd.distWater <- sd(distWater[!is.na(distWater)])
distWater <- (distWater-mean.distWater)/sd.distWater     # Standardise distWater
distWater[is.na(distWater)] <- 0               # Impute zeroes (means)

mean.elevation <- mean(elevation, na.rm = TRUE)
sd.elevation <- sd(elevation[!is.na(elevation)])
elevation <- (elevation-mean.elevation)/sd.elevation     # Standardise elevation
elevation[is.na(elevation)] <- 0               # Impute zeroes (means)

mean.treeDensity <- mean(treeDensity, na.rm = TRUE)
sd.treeDensity <- sd(treeDensity[!is.na(treeDensity)])
treeDensity <- (treeDensity-mean.treeDensity)/sd.treeDensity     # Standardise treeDensity
treeDensity[is.na(treeDensity)] <- 0               # Impute zeroes (means)


#----- 4 - Single-season occupancy model -----

# Specify model in JAGS language
sink(here("bin", "model.jags"))
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)
beta1.psi ~ dnorm(0, 0.01) # dist.water
beta2.psi ~ dnorm(0, 0.01) # slope
beta3.psi ~ dnorm(0, 0.01) # elevation
beta4.psi ~ dnorm(0, 0.01) # tree.density
alpha.p ~ dnorm(0, 0.01)
#beta1.p ~ dnorm(0, 0.01)
#beta2.p ~ dnorm(0, 0.01)
#beta3.p ~ dnorm(0, 0.01)
#beta4.p ~ dnorm(0, 0.01)

# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) {
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
   lpsi[i] <- alpha.psi + beta1.psi * distWater[i] + beta2.psi * slope[i] + beta3.psi * elevation[i] + beta4.psi * treeDensity[i]

   # Observation model for the observations
   for (j in 1:T) {
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      #lp[i,j] <- alpha.p + beta1.p * DATES[i,j] + beta2.p * pow(DATES[i,j], 2) + beta3.p * HOURS[i,j] + beta4.p * pow(HOURS[i,j], 2)
      lp[i,j] <- alpha.p
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = y, R = nrow(y), T = ncol(y), distWater = distWater, slope = slope, elevation = elevation, treeDensity = treeDensity)

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
zst[zst == -Inf] <- 0
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.psi", "beta1.psi", "beta2.psi", "beta3.psi", "beta4.psi", "mean.p", "occ.fs", "alpha.p")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(win.data, inits, params, here("bin", "model.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

# Posterior distribution of the number of occupied sites in actual sample
hist(out$BUGSoutput$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied sites", )
#abline(v = 10, lwd = 2) # The observed number


# Examine how mean species-specific occupancy changes by elevation

x = seq(min(elevation), max(elevation), by=0.1) # use max(recovery) instead of 0.06 if all sites are included
y = (x*sd.elevation) + mean(elevation)
beta3.psi = out$BUGSoutput$sims.list$beta3.psi # beta3.psi is the slope for elevation effecta1 <- a1bySpecies

plot(y,y, type="l", ylim=c(0,1), xlim=c(min(y),max(y)), main="Myrmecophaga trydactyla relationship with elevation",
     col="white", ylab="Occupancy probability", xlab="Elevation (standardized)")

#for (i in 1:n) {
  elev <- (mean(out$BUGSoutput$sims.list$occ.fs)/61) + mean(beta3.psi)*x
  #lines(y,plogis(LOGGED), type="l", ylim=c(0,1), xlim=c(min(y),max(y)), main=i)
  lines(y,plogis(elev), type="l", ylim=c(0,1), xlim=c(min(y),max(y)), lwd=0.5, col="dark grey") # because max recovery in logged.occ is 5
#}
# add bold line for mean across species
loggedBold <- mean(elev) + mean(elev)*x
lines(y,plogis(loggedBold), type="l", ylim=c(0,1), xlim=c(min(y),max(y)), lwd=3) # lwd for line thickness



#Pstar <- array(NA, dim = c(out$BUGSoutput$n.sims, 10))
x <- cbind(rep(1, 3000), rep(2, 3000), rep(3, 3000), rep(4, 3000), rep(5, 3000), rep(6, 3000), rep(7, 3000), rep(8, 3000), rep(9, 3000), rep(10, 3000)) 
for (i in 1:out$BUGSoutput$n.sims) {
  for (j in 1:10){
    Pstar[i,j] <- 1 - (1 - out$BUGSoutput$sims.list$mean.p[i])^j
  } #j
} #i

boxplot(Pstar ~ x, col = "gray", las = 1, ylab = "Pstar", xlab = "Number of surveys", outline = FALSE)
abline(h = 0.95, lty = 2, lwd = 2)

par(mfrow = c(2, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.psi), nclass = 40, col = "gray", main = "Forest interior", xlab = "Occupancy probability", xlim = c(0, 1))
hist(plogis(out$BUGSoutput$sims.list$alpha.psi+ out$BUGSoutput$sims.list$beta.psi), nclass = 40, col = "gray", main = "Forest edge", xlab = "Occupancy probability", xlim = c(0, 1))

# Predict effect of time of day with uncertainty
mcmc.sample <- out$BUGSoutput$n.sims

original.date.pred <- seq(0, 60, length.out = 30)
original.hour.pred <- seq(180, 540, length.out = 30)
date.pred <- (original.date.pred - mean.date)/sd.date
hour.pred <- (original.hour.pred - mean.hour)/sd.hour
p.pred.date <- plogis(out$BUGSoutput$mean$alpha.p + out$BUGSoutput$mean$beta1.p * date.pred + out$BUGSoutput$mean$beta2.p * date.pred^2 )
p.pred.hour <- plogis(out$BUGSoutput$mean$alpha.p + out$BUGSoutput$mean$beta3.p * hour.pred + out$BUGSoutput$mean$beta4.p * hour.pred^2 )

array.p.pred.hour <- array.p.pred.date <- array(NA, dim = c(length(hour.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.p.pred.date[,i] <- plogis(out$BUGSoutput$sims.list$alpha.p[i] + out$BUGSoutput$sims.list$beta1.p[i] * date.pred + out$BUGSoutput$sims.list$beta2.p[i] * date.pred^2)
  array.p.pred.hour[,i] <- plogis(out$BUGSoutput$sims.list$alpha.p[i] + out$BUGSoutput$sims.list$beta3.p[i] * hour.pred + out$BUGSoutput$sims.list$beta4.p[i] * hour.pred^2)
}

# Plot for a subsample of MCMC draws
sub.set <- sort(sample(1:ni, size = 200))

par(mfrow = c(2, 1))
plot(original.date.pred, p.pred.date, main = "", ylab = "Detection probability", xlab = "Date (1 = 1 July)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
  lines(original.date.pred, array.p.pred.date[,i], type = "l", lwd = 1, col = "gray")
}
lines(original.date.pred, p.pred.date, type = "l", lwd = 3, col = "blue")

plot(original.hour.pred, p.pred.hour, main = "", ylab = "Detection probability", xlab = "Time of day (mins after noon)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
  lines(original.hour.pred, array.p.pred.hour[,i], type = "l", lwd = 1, col = "gray")
}
lines(original.hour.pred, p.pred.hour, type = "l", lwd = 3, col = "blue")
