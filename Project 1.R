#HW 1
dir = ("C:/Users/Qudus Uthman/Desktop/Site B/scantic river")
setwd(dir)

#install.packages("ncdf4")
#install.packages("raster")
#install.packages("rgdal")
#install.packages("RNetCDF")

library(ncdf4)
library(raster)
library(rgdal)
library(RNetCDF)

#### GATHER DATA ####
#load shapefiles
new = readOGR(dsn=dir,layer="scant") #First site
new.extent = extent(new)
vaupes = readOGR(dsn=dir,layer="uaupes") #Second site
vaupes.extent = extent(vaupes)

## Manually download files, we want aet, pet, precipitation, max temp, min temp 
#Enter lat and lon ranges
lat.range=c(new.extent@ymin,new.extent@ymax)        
lon.range=c(new.extent@xmin,new.extent@xmax)

# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
vars = c("aet","pet","ppt","tmax","tmin")
var=vars[1]
#PART 1: run for each variable##############################################################################################
baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")
nc <- open.nc(baseurlagg)
lon <- var.get.nc(nc, "lon")
lat <- var.get.nc(nc, "lat")
lat.range <- sort(lat.range)                              #!sort user input values from low to high
lon.range <- sort(lon.range)
lat.index <- which(lat>=lat.range[1]&lat<=lat.range[2])    #! index values within specified range
lon.index <- which(lon>=lon.range[1]&lon<=lon.range[2])    
lat.n1 <- length(lat.index)                                #!value for count
lon.n1 <- length(lon.index)
start <- c(lon.index[1], lat.index[1], 1)
count <- c(lon.n1, lat.n1, NA)                            #! parameter change: 'NA' instead of '-1' to signify entire dimension
data1 <-var.get.nc(nc, variable = var,start = start, count,unpack=TRUE)   
#PART 2: Save variable after running #############################################################################################################
new.aet = apply(data1,MARGIN=3,FUN=mean) #average over each pixel value (not masking)
new.pet = apply(data1,MARGIN=3,FUN=mean)
new.ppt = apply(data1,MARGIN=3,FUN=mean)
new.tmax = apply(data1,MARGIN=3,FUN=mean)
new.tmin = apply(data1,MARGIN=3,FUN=mean)
new.t <- cbind(new.tmax,new.tmin)
new.tavg = rowMeans(new.t)
#PART 3: Combine all data, add dates##############################################################################################################
new.data = as.data.frame(cbind(new.aet,new.pet,new.ppt,new.tavg))
#new.data = new.data[c(265:624),] #extract years 1980-2010, (22*12+1):(52*12)
#assign dates
new.data$MY <- seq(as.Date("1958/1/1"), by = "month", length.out = length(new.data[,1]))
write.csv(new.data, file="site_1.csv")
