rm(list=ls())

# # Check package availability
# packages.needed<-c("rgdal","plotKML","chron","spatstat","maptools","raster")
# if (length(setdiff(packages.needed, rownames(installed.packages()))) > 0) {
#   install.packages(setdiff(packages.needed, rownames(installed.packages())))  
# }

# Load packages
library(rgdal)
# library(plotKML)
library(spatstat)
library(maptools)
library(chron)
library(raster)
library(ggplot2)
library(rlist)
library(plotrix)
library(zoo)
library(dplyr)
library(rgeos)
library(sf)

# README ------------------------------------------------------------------

# Environment settings ----------------------------------------------------

# Working directory
# wd<-"/home/elic/fclapuyt/tracking-gps/"
# wd<-"C:/Users/fclapuyt/Nextcloud/Students/Therese/tracking-gps"
wd<-"/home/franz/Nextcloud/Students/Therese/tracking-gps"
# wd<-"~/Nextcloud/Students/Therese/tracking-gps"

# Subfolders
dir.area<-"data/malia"
dir.inputs<-file.path(dir.area,"tracks-utm")
dir.outputs<-file.path(dir.area,"outputs")
dir.src<-"src"
dir.functions<-"functions"
dir.shp<-"shp"
dir.aoi<-"aoi"
dir.graphs<-"graphs"
dir.lines<-"lines"
dir.points<-"points"
dir.thresholded<-"thresholded"
dir.tables<-"tables"
dir.museum<-"museum"

# Variables ---------------------------------------------------------------

# Define variables

# Name of the shapefile containing the limits of the area
limits.area<-"limits-area"
# EPSG number of local projection coordinate system (UTM) 
id.proj<-"+init=epsg:32635"


# Working directory and setup ---------------------------------------------

# Define working directory
setwd(wd)

# Create missing directories
if(!dir.exists(file.path(wd,dir.outputs))){
  dir.create(file.path(wd, dir.outputs))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.shp))){
  dir.create(file.path(wd,dir.outputs,dir.shp))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.graphs))){
  dir.create(file.path(wd,dir.outputs,dir.graphs))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.shp,dir.aoi))){
  dir.create(file.path(wd,dir.outputs,dir.shp,dir.aoi))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.shp,dir.lines))){
  dir.create(file.path(wd,dir.outputs,dir.shp,dir.lines))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.shp,dir.points))){
  dir.create(file.path(wd,dir.outputs,dir.shp,dir.points))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.shp,dir.thresholded))){
  dir.create(file.path(wd,dir.outputs,dir.shp,dir.thresholded))  
}
if(!dir.exists(file.path(wd,dir.outputs,dir.shp,dir.museum))){
  dir.create(file.path(wd,dir.outputs,dir.shp,dir.museum))  
}

# Load custom functions
# source(file.path(wd,dir.src,"distLatLon.R"))
source(file.path(wd,dir.src,dir.functions,"shiftVector.R"))


# Script GPS Tracking -----------------------------------------------------

# Read museum shapefile
museum.spdf<-readOGR(dsn=file.path(wd,dir.area,"museum","museum.shp"),stringsAsFactors=FALSE)
museum.spdf<-spTransform(museum.spdf,CRS(id.proj))

# Read entrance shapefile
ins.spdf<-readOGR(dsn=file.path(wd,dir.area,"museum","museum-doors.shp"),stringsAsFactors=FALSE)
ins.spdf<-spTransform(ins.spdf,CRS(id.proj))

#list of gpx files in the working directory defined as input
lines.files<-list.files(file.path(wd,dir.outputs,dir.shp,dir.lines),pattern=".shp")
ids.files<-unname(sapply(lines.files, function(x) as.numeric(unlist(strsplit(x, split="-"))[2])))
# points.files<-list.files(file.path(wd,dir.outputs,dir.shp,dir.points),pattern=".shp")

id<-16
# for(id in ids.files){
  
  # Read lines and points for current track 
  gpx.lines<-spTransform(readOGR(file.path(wd,dir.outputs,dir.shp,dir.lines,paste0("Track-",id,"-lines.shp"))), id.proj)
  gpx.points<-spTransform(readOGR(file.path(wd,dir.outputs,dir.shp,dir.points,paste0("Track-",id,"-points.shp"))), id.proj)

  # Create intersection between track and museum doors (as lines)  
  ints.points<-gIntersection(gpx.lines,ins.spdf, byid=TRUE, checkValidity=FALSE)

  # Create attribute table and create spatial point data frame of intersections between tracks and entrances
  ints.attr<-data.frame(id=seq(1,length(ints.points),1))
  ints.spdf<-SpatialPointsDataFrame(ints.points,ints.attr)
  writeOGR(ints.spdf, dsn = file.path(wd,dir.outputs,dir.shp,dir.museum,paste0("intersections-",id,"-points.shp")), layer="test", driver = "ESRI Shapefile")
  
  # int.id<-1
  # # for(int.id in ints.spdf@data$id){
  #   
  #   int<-ints.spdf[ints.spdf@data$id == int.id,]
  #   
  # # }
  
  # Get gpx tracks segments which intersect entrances
  ints.lines<-gIntersects(ins.spdf,gpx.lines,byid=TRUE)
  ints.lines<-gpx.lines[as.numeric(rownames(data.frame(ints.lines)%>%filter_all(any_vars(. %in% TRUE))))+1,]
  # plot(ints.lines)
  # plot(ints.spdf,add=TRUE)

  # For each segment crossing an entrance, get points of gpx tracks that are just before and 
  # just after an intersection with a crossing (in order to compute time and stuff...).
  vertices.list<-list()
  line.id<-ints.lines@data$id[1]
  for(line.id in ints.lines@data$id){
    
    vertices<-st_intersects(st_as_sf(ints.lines[ints.lines@data$id == line.id,]),st_as_sf(gpx.points))
    vertices<-gpx.points[gpx.points@data$id == vertices[[1]],]
    
    df<-data.frame(vertices@data)
    df$time<-as.character(df$time)
    df$time_p2<-as.character(df$time_p2)
    
    vertices.list<-list.append(vertices.list,df)
    
  }
  
    
  
# }
  

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# .########..########.##.....##
# .##.....##.##.......##.....##
# .##.....##.##.......##.....##
# .##.....##.######...##.....##
# .##.....##.##........##...##.
# .##.....##.##.........##.##..
# .########..########....###...
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Point pattern analysis --------------------------------------------------
# Useful documentation: https://mgimond.github.io/Spatial/

# Read shapefile containing limits of the study area
# area.spdf<-readOGR(dsn=file.path(wd, dir.inputs, list.files(path=file.path(wd, dir.inputs, dir.shp), pattern=".shp")))

# # Divide study area into quadrats
# area.rstr<-raster(area.spdf)
# res(area.rstr)<-100
# area.quads <- as(rasterize(area.spdf, area.rstr), 'SpatialPolygons')
# 
# # Display quads and gpx tracks
# plot(quads)
# points(track.spdf.full, col="red", cex=.5)

# # Convert limits to owin type and get coordinates of points
# area.owin<-as.owin(area.spdf)
# gpx.pts<-coordinates(track.spdf.full)
# 
# # Point pattern analysis
# gpx.ppp<-ppp(gpx.pts[, 1], gpx.pts[, 2], window=area.owin)
# plot(gpx.ppp)
# gpx.density<-density(gpx.ppp)
# plot(gpx.density)
# 
# # Count points within quadrats
# pts.count<-quadratcount(gpx.ppp, nx=50, ny=50)
# # Plot count
# plot(gpx.ppp, pch=".", cols="grey70")
# plot(pts.count, add=TRUE)
# 
# # Compute density of points within quadrats
# pts.intensity<-intensity(pts.count)
# plot(intensity(pts.count, image=TRUE))
# plot(gpx.ppp, pch=".", add=TRUE, col="green", main="Density of points")
# 
# pts.tess<-as.tess(pts.count)
# pts.spdf<-as(pts.tess, "SpatialPolygons")
# 
# pts.count$tile
# 
# plot(pts.spdf)

# owin2Polygons <- function(x, id="1") {
#  stopifnot(is.owin(x))
#  x <- as.polygonal(x)
#  closering <- function(df) { df[c(seq(nrow(df)), 1),] }
#  pieces <- lapply(x$bdry, 
#          function(p) {
#           Polygon(coords=closering(cbind(p$x, p$y)), 
#               hole=is.hole.xypolygon(p)) })
#  z <- Polygons(pieces, id)
#  return(z)
# }
# 
# tess2SP <- function(x) {
#  stopifnot(is.tess(x))
#  y <- tiles(x)
#  nam <- names(y)
#  z <- list()
#  for(i in seq(y))
#   z[[i]] <- owin2Polygons(y[[i]], nam[i])
#  return(SpatialPolygons(z))
# }
# 
# owin2SP <- function(x) {
#  stopifnot(is.owin(x))
#  y <- owin2Polygons(x)
#  z <- SpatialPolygons(list(y))
#  return(z)
# }
# 
# test<-owin2SP(test)
# is(pts.count)