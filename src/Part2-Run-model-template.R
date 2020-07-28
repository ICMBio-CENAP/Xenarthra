# Occupancy model for Jaguars at Parque Estadual do Rio Doce
# Juliana Pasa, Fernando Azevedo, Elildo Carvalho Jr

# Part2-Occupancy models

#the data basis used here was organized manually, and contains the following columns:Cam.Trap.Name, Photo.Date, bin, n, Start.Date, End.Date, Sampling.Period, Sampling.Unit.Name.
#which consists in respectively:name of the point (including those without records of the species), record date, species name, number of individuals (optional), start date of the sampling, end date of the sampling, year of the sampling, name of the point sampled. 

# Read in TEAM data set and create and format so it is ready for wildlife community model
# Written by Jorge Ahumada @ Conservation International
# Adapted by Elildo Carvalho Jr @ ICMBio/CENAP, 2020-04-02

##----- 1 - Load libraries-----
library(dplyr)
library(lubridate)
library(reshape2)
library(ggplot2)
library(unmarked)
library(here)


##----- 2 - Source files-----
source(here("bin", "camera trap analysis functions-10-06-18.R"))
source(here("bin", "f-matrix-creator.R"))

#Load libraries and start from here (using the data already organized manually):

## ---- 3- Load data-------   
jaguars <- read.csv(here("data", "jaguars_data_systematic.csv"))

# fix dates and times
jaguars$Photo.Date <- as.Date(jaguars$Photo.Date)
jaguars$Start.Date <- as.Date(jaguars$Start.Date)
jaguars$End.Date <- as.Date(jaguars$End.Date)
#jaguars$Photo.Time <- as.POSIXct(paste(jaguars$Photo.Date, jaguars$Photo.Time))

# create a few additional columns
#jaguars$td.photo <- ymd_hms(paste(jaguars$Photo.Date, jaguars$Photo.Time, sep=" "))

# filter non-independent records (30 minute interval)
#jaguars <- f.separate.events(jaguars, 30)
#jaguars <- distinct(jaguars, Camera.Trap.Name, bin, grp, .keep_all=TRUE) # keep only independent records


## ---- 4- Create a matrix of presence/absence-------    
# Create a matrix of presence/absence of each species
# rows are camera trap points (sampling units) and columns are dates
# Use function f.matrix.creator2 to do this - cria matriz de deteccao, colapsando os dados sempre em 15 dias
paMatJAGUAR <- f.matrix.creator4(jaguars, jaguars$Sampling.Unit.Name, "Panthera onca", 15) # 15 is the desired number of occasions, you can chabge it at your will 
print(paMatJAGUAR)

# Possivel explicação para linha acima não funcionar: Photo.Date anterior a Start.Date
#Confira:
#jaguars$bla <- jaguars$Photo.Date-jaguars$Start.Date
#View(jaguars) # veja os valores de bla, alguns sao negativos


# This creates a list were each element of the list is a presence/absence matrix for a species
summary(paMatJAGUAR)
names(paMatJAGUAR)

# Look at the matrix for Panthera onca
View(paMatJAGUAR[["Panthera onca"]])

## ---- 5- Load covariate data-------    
# Still dont have the covariates, so create dummy variaables just to test the code
Camera.Trap.Name <- rownames(paMatJAGUAR[["Panthera onca"]])
#Camera.Trap.Name <- sort(unique(jaguars$Camera.Trap.Name))
altitude <- rnorm(length(Camera.Trap.Name)); slope <- rexp(length(Camera.Trap.Name),1/2500) # dummy variables
covsPERD <- data.frame(Camera.Trap.Name, altitude, slope)

# We just need the cameras operating in 2017. Which ones are they?
suPERD <- unique(jaguars$Camera.Trap.Name)

# filter our covariate file - only keep the lines that matches suPERD
covsPERD <- filter(covsPERD, Camera.Trap.Name %in% suPERD)

# Normalize the two covariates - altitude and slope
covsPERD <- mutate(covsPERD, norm.altitude = (altitude - mean(altitude))/sd(altitude), norm.slope = (slope-mean(slope))/sd(slope))
head(covsPERD)

#Sort by camera trap unit name for analysis
covsPERD <- arrange(covsPERD, Camera.Trap.Name)
head(covsPERD)

## ---- 6- Run occupancy model -------    

# Convert these matrices into a special format to analyze them in unmarked
#neste comando que está dando o erro:
umPonca <- unmarkedFrameOccu(y = paMatJAGUAR[["Panthera onca"]], siteCovs = covsPERD)
summary(umPonca)

## fit a single season occupancy model with no covariates
occMod0_Ponca <- occu(~1 ~1, umPonca)
# Look at the model results
occMod0_Ponca
# Transform the estimates from log to linear
backTransform(occMod0_Ponca, "state")
backTransform(occMod0_Ponca, "det")

## fit model with covariate for occupancy
occMod1_Ponca <- occu(~norm.altitude ~1, umPonca)
occMod1_Ponca
backTransform(occMod1_Ponca, "state")

# We cannot use backTransform for p because it is dependent on elevation
# First create a dataframe with the average value of elevation (0)
newdata <- data.frame(norm.altitude = 0)
# Use predict to get the values
predict(occMod1_Ponca, type="det", newdata=newdata)

# How is detection varying with elevation
newdata <- data.frame(norm.altitude = seq(-2.8, 2, 0.05))
mod_pred <- predict(occMod1_Ponca, type="det", newdata=newdata)
mod_pred <- data.frame(mod_pred, elevation = newdata)

# Transform elevation back to original scale
mod_pred <- mutate(mod_pred, elev = sd(covsPERD$altitude)*norm.altitude + mean(covsPERD$altitude))
mod_pred

# plot it
ggplot(mod_pred, aes(x = elev, y = Predicted)) + geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "blue") + xlab("Elevation") + ylab("Detection probability")

# save as jpeg
jpeg(here("results", "dummy_elevation_psi.jpg"), width = 800, height = 600) # Open jpeg file
ggplot(mod_pred, aes(x = elev, y = Predicted)) + geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "blue") + xlab("Elevation") + ylab("Detection probability")
dev.off()

# Do another model where occupancy is a function of slope and detection is a function of altitude
occMod2_Ponca <- occu(~norm.slope ~norm.altitude, umPonca)
occMod2_Ponca

# Compare these models
modelList <- fitList('psi()p()' = occMod0_Ponca, 'psi(altitude)p()' = occMod1_Ponca, 'psi(slope)p(altitude)' = occMod2_Ponca)
modSel(modelList)

