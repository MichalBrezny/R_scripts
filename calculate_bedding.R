# Script for calculation of interpolated dip direction and dip
# Based on Santangelo, Michele et al. 2015. „A Method for the Assessment of the 
# Influence of Bedding on Landslide Abundance and Types". Landslides 12(2): 295–309.


#Input rasters have to be prepared in GIS



#load packages
library(raster)
library(rgdal)


#input raster data from GIS
#modify path

#as vector
nx <- as.vector(raster("rasters/rx.tif"))
ny <- as.vector(raster("rasters/ry.tif"))
nz <- as.vector(raster("rasters/rz.tif"))
r <- as.vector(raster("rasters/r.tif"))
m <- as.vector(raster("rasters/m.tif"))

#as raster
rx <- raster("rasters/rxtif")

########################

# calculate dip direction function

dipdir_calc <- function(nx,ny, nz ,m){
        z <- rep(NA, length(nx))
        i <- which(nz >= 0 & nx>=0 & ny >= 0) #for normal bedding 
        z[i] <- asin(nx[i]/m[i])
        i <- which(nz >= 0 & nx<0 & ny >=0)
        z[i] <- 2*pi-asin(-nx[i]/m[i])
        i <- which(nz >= 0 & nx<0 & ny<0)
        z[i] <- pi+asin(-nx[i]/m[i])
        i <- which(nz >= 0 & nx>=0 & ny<0)
        z[i] <- pi-asin(nx[i]/m[i])
        i <- which(nz < 0 & nx>=0 & ny >= 0)  #overturned bedding
        z[i] <- pi + asin(nx[i]/m[i])
        i <- which(nz < 0 & nx<0 & ny >=0)
        z[i] <- pi-asin(-nx[i]/m[i])
        i <- which(nz < 0 & nx<0 & ny<0)
        z[i] <- asin(-nx[i]/m[i])
        i <- which(nz < 0 & nx>=0 & ny<0)
        z[i] <- asin(nx[i]/m[i])
        return(z/pi*180)
}

# calculate dip 

dip_calc <- function(nz,r){
        z <- rep(NA, length(nz))
        ia <- which(nz>=0) 
        z[ia] <- acos(nz[ia]/r[ia])
        ib <- which(nz<0)
        z[ib] <- pi-acos(-nz[ib]/r[ib])
        return(z/pi*180)
}

#####
#Calculate dip direction and dip
dipdir <- dipdir_calc(nx,ny,nz,m)
dip <- dip_calc(nz,r)

#####################
#convert to raster and save 
rxm <- as.matrix(rx)

#dip direction
ddm <- as.matrix(dipdir, nrow=nrow(rxm), ncol= ncol(rxm), byrow=T)
dipdir_raster <- raster(ddm, template = rx)
plot(azimuth_raster)
writeRaster(azimuth_raster, "azimuth_rast.tif", overwrite =T)

#dip
dm <- as.matrix(dip, nrow=nrow(rxm), ncol= ncol(rxm), byrow=T)
dip_raster <- raster(dm, template = rx)
plot(dip_raster)
writeRaster(dip_raster, "dip_rast.tif",overwrite =T)

#end of the script





