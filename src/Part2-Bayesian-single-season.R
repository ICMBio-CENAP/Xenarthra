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
Pmaximus <- readRDS(here("data", "Pmaximus.rds"))
y <- Mtridactyla[,2:16]
#y <- Pmaximus[,2:16]
str(y)
SiteCovs <- Mtridactyla[,17:23]

# check corr in SiteCovs
cor(SiteCovs)
# elevation and slope 0.7 corr

names(SiteCovs)
landCover <- SiteCovs[,1]
distWater <- SiteCovs[,2]
slope <- SiteCovs[,3]
elevation <- SiteCovs[,4]
treeBurned <- SiteCovs[,5]
basalArea <- SiteCovs[,6]
treeDensity <- SiteCovs[,7]


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

mean.treeBurned <- mean(treeBurned, na.rm = TRUE)
sd.treeBurned <- sd(treeBurned[!is.na(treeBurned)])
treeBurned <- (treeBurned-mean.treeBurned)/sd.treeBurned     # Standardise treeBurned
treeBurned[is.na(treeBurned)] <- 0               # Impute zeroes (means)

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
beta2.psi ~ dnorm(0, 0.01) # elevation
beta3.psi ~ dnorm(0, 0.01) # tree.density
beta4.psi ~ dnorm(0, 0.01) # burned.trees
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
   lpsi[i] <- alpha.psi + beta1.psi * distWater[i] + beta2.psi * elevation[i] + beta3.psi * treeDensity[i] + beta4.psi * treeBurned[i]

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
dataJAGS <- list(y = y, R = nrow(y), T = ncol(y), distWater=distWater, elevation=elevation, treeDensity=treeDensity, treeBurned=treeBurned)

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
zst[zst == -Inf] <- 0
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3))}

# Parameters monitored
params <- c("alpha.psi", "beta1.psi", "beta2.psi", "beta3.psi", "beta4.psi", "mean.p", "occ.fs", "alpha.p", "z")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call JAGS from R (BRT < 1 min)
out <- jags(dataJAGS, inits, params, here("bin", "model.jags"), n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

# Posterior distribution of the number of occupied sites in actual sample
hist(out$BUGSoutput$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied sites", )
#abline(v = 10, lwd = 2) # The observed number


# A template for a figure with effect of a covariate on occupancy with uncertainty:
# Using elevation as a model
mcmc.sample <- out$BUGSoutput$n.sims
original.elev <- SiteCovs[,4]
original.elev.pred <- seq(min(original.elev), max(original.elev), length.out = 30)
elev.pred <- (original.elev.pred - mean.elevation)/sd.elevation
p.pred.elev <- rep(NA, length(elev.pred))
for(i in 1:length(p.pred.elev)) {
  p.pred.elev[i] <- plogis(out$BUGSoutput$mean$alpha.psi + out$BUGSoutput$mean$beta2.psi*elev.pred[i])
}
array.p.pred.elev <- array(NA, dim = c(length(elev.pred), mcmc.sample))
for (i in 1:mcmc.sample){
  array.p.pred.elev[,i] <- plogis(out$BUGSoutput$sims.list$alpha.psi[i] + out$BUGSoutput$sims.list$beta2.psi[i]*elev.pred)
}

# Plot for a subsample of MCMC draws
# write as a function:
plot.test <- function() {
  sub.set <- sort(sample(1:mcmc.sample, size = 200))

  plot(original.elev.pred, p.pred.elev, main = "", ylab = "Occupancy probability", xlab = "Elevation (m)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
  for (i in sub.set){
    lines(original.elev.pred, array.p.pred.elev[,i], type = "l", lwd = 1, col = "gray")
  }
  lines(original.elev.pred, p.pred.elev, type = "l", lwd = 3, col = "blue")
}

plot.test()

# save as jpeg
jpeg(here("results", "elevation_effect_test.jpg"), width = 800, height = 400) # Open jpeg file
plot.test()
dev.off()

