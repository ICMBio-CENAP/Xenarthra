# Occupancy model for Xenarthrans at Gurupi Biological Reserve
# Elildo Carvalho Jr, Alexandre Martins, Eloisa Mendonca

# Part1-Prepare data

# Read in TEAM data set and create and format so it is ready for wildlife community model
# Written by Jorge Ahumada @ Conservation International
# Adapted by Elildo Carvalho Jr @ ICMBio/CENAP, 2020-04-02

##----- 1 - Load libraries-----
library(dplyr)
library(lubridate)
library(reshape2)
library(activity)
library(ggplot2)
library(ggmap)
library(here)


##----- 2 - Source files-----
source(here("bin", "camera trap analysis functions-10-06-18.R")) # using package here to build a path to the subdirectory "bin"

## ----3 - Load data-------
gurupi2016 <- read.csv(here("data", "gurupi2016.csv"))
#gurupi2016 <- f.readin.fix.data(here("data", "gurupi2016.csv"))

# fix dates and times
gurupi2016$Photo.Date <- as.Date(gurupi2016$Photo.Date)
gurupi2016$Camera.Start.Date <- as.Date(gurupi2016$Camera.Start.Date)
gurupi2016$Camera.End.Date <- as.Date(gurupi2016$Camera.End.Date)
#gurupi2016$Photo.Time <- as.POSIXct(paste(gurupi2016$Photo.Date, gurupi2016$Photo.Time))

# create a few additional columns
gurupi2016$bin <- factor(paste(gurupi2016$Genus,gurupi2016$Species, sep=" "))
gurupi2016$td.photo <- ymd_hms(paste(gurupi2016$Photo.Date, gurupi2016$Photo.Time, sep=" "))
gurupi2016$Start.Date <- gurupi2016$Camera.Start.Date
gurupi2016$End.Date <- gurupi2016$Camera.End.Date
gurupi2016$Sampling.Period <- gurupi2016$Sampling.Event
gurupi2016$Sampling.Unit.Name <- gurupi2016$Camera.Trap.Name

# How many rows and columns?
dim(gurupi2016)

# Look at the first 6 rows of data
head(gurupi2016)

# filter non-independent records (30 minute interval)
gurupi2016 <- f.separate.events(gurupi2016, 30)
gurupi2016 <- distinct(gurupi2016, Camera.Trap.Name, bin, grp, .keep_all=TRUE) # keep only independent records

## ----Number of images per camera trap-----------
imgsPerCT <- gurupi2016 %>% group_by(Camera.Trap.Name, Latitude, Longitude) %>% summarize(n = n()) %>% arrange(desc(n))
imgsPerCT

# plot it
p <- ggplot(imgsPerCT, aes(x = reorder(Camera.Trap.Name, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Camera Trap Name") + ylab("Number of images")
p

## ----plot deployment in time---------------------------------------------

# Look at how camera traps were deployed
# First create a table
tabCTN_PD <- with(gurupi2016, table(Camera.Trap.Name, Photo.Date))
head(tabCTN_PD)
# Get it ready for ggplot
tabCTN_PD <- melt(tabCTN_PD)
# Plot it
p <- ggplot(tabCTN_PD, aes(y = Camera.Trap.Name, x = Photo.Date)) + geom_raster(aes(fill=value), alpha = 0.8) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p


##---- Species richness ---------------------------
unique(gurupi2016$bin)

#Remove the blank species
# Let's do this in a copy of gurupi2016
gurupi2016.copy <- gurupi2016

gurupi2016.copy$bin <- droplevels(gurupi2016.copy$bin, exclude = " ")
unique(gurupi2016.copy$bin)

# How often they show up in the camera traps
imgsPerSp <- gurupi2016.copy %>% group_by(bin) %>% summarize(n = n()) %>% arrange(desc(n))
print(imgsPerSp, n=Inf)

#plot
ggplot(imgsPerSp, aes(x = reorder(bin, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Species") + ylab("Number of images")



##---- Activity patterns ---------------------------

# Write a couple of functions for activity patterns
calcActivity <- function(dataset, speciesName, threshold){
  # Extract the data from species
  spdata <- filter(dataset, bin == speciesName)
  
  # Group by events that are two minutes apart
  spdata <- f.separate.events(spdata, threshold)
  
  # Order the data sequentially
  spdata <- f.order.data(spdata)
  
  # Calculate for each event the time
  # First get the median time for each event
  acTime <- spdata %>% group_by(grp) %>% summarize(median = median(td.photo))
  # Then extract just the hh:mm:ss information
  acTime <- f.extracthms(acTime$median)
  # Extract the hour information for the circular plot
  acTime_hour <- data.frame(seq = 1:length(acTime), hour = hour(acTime))
  # Convert to radian time
  acTime <- f.timeformatfunc(acTime)
  # fit an activity kernel
  actModel <- fitact(acTime, sample = "model", reps=100)
  # Plot it in a couple ways
  
  plot(actModel, main = speciesName)
  
  list(actModel, acTime_hour)
}

plot.circular.activity <- function(actdata, speciesName) {
  ggplot(actdata, aes(x = hour)) + geom_histogram(breaks = seq(0, 24), colour = "grey") + coord_polar(start = 0)  + theme_minimal() + ggtitle(speciesName) + scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24))
}

# Run it for Myrmecophaga tridactyla
actMod_Myrtri <- calcActivity(gurupi2016, "Myrmecophaga tridactyla", 2)

# Do a ciruclar plot of activity
plot.circular.activity(actMod_Myrtri[[2]], "Myrmecophaga tridactyla")

# Run it for Tamandua tetradactyla
actMod_Tamtet <- calcActivity(gurupi2016, "Tamandua tetradactyla", 2)

# Compare the activity of two species
compareAct(list(actMod_Myrtri[[1]],actMod_Tamtet[[1]]))



##---- Generate spatial distributions ---------------------------

# Start with provide the lon/lat range of the data
lon <- range(gurupi2016$Longitude)
lat <- range(gurupi2016$Latitude)

# Extract the unique lat/lons and put them on a data frame
locationsGurupi <- unique(cbind(as.character(gurupi2016$Camera.Trap.Name), gurupi2016$Latitude,gurupi2016$Longitude))
locationsGurupi <- data.frame(Camera.Trap.Name = locationsGurupi[,1], Latitude = as.numeric(locationsGurupi[,2]), Longitude = as.numeric(locationsGurupi[,3]))
locationsGurupi <- dplyr::arrange(locationsGurupi, Camera.Trap.Name)

# If you have internet: Download the map from google
map <- get_map(location = c(c(lon[1],lat[1]),c(lon[2],lat[2])), zoom = 10, source = "google", maptype = "terrain")

# Plot the locations of the camera traps
ggmap(map, extent = "normal", maprange = T) + geom_point(data=locationsGurupi, aes(x = Longitude, y = Latitude), colour="red", size = 0.1)

# Plot the number of images per camera trap point
ggmap(map, extent = "normal", maprange = T) + geom_point(data = imgsPerCT, aes(x = Longitude, y = Latitude, color = n), size = 0.5)

# Plot as a surface
ggmap(map, extent = "device", legend = "topleft")  + stat_density2d(aes(x = Longitude, y = Latitude, fill = ..level..), data = gurupi2016, geom = "polygon", size = 2, bins = 10, alpha = 0.5)

# Do it for Myrmecophaga
# First extract the number of photographic events by camera trap
sevPerCT_Myrtri <- Myrtri %>% group_by(Camera.Trap.Name, Latitude, Longitude) %>% summarize(n = length(unique(grp)))

# Then put them in a map
ggmap(map) + geom_point(data = sevPerCT_Myrtri, aes(x = Longitude, y = Latitude, color = n), size = 0.5) + ggtitle("Dasyprocta prymnolopha")

