# Reset working directory
rm(list=ls())

# README ------------------------------------------------------------------

# The script "count_visitors" computes the number of visitors passing by 
# an area of interest, along with the time spent in the AOI.

# File structure 
# -- Working directory (wd)
#   -- src (everything but data)
#       - tracking-visitors.R
#     -- functions 
#       - shiftVector.R
#   -- inputs
#     -- tracks
#        - track-1.shp
#        - track-2.shp
#        - track-3.shp
#        - ...
#     -- aoi
#        - aoi_for_mean_dir.shp
#    -- outputs (will be created by the script)

# Outputs and required subfolders will automatically be created by the script. 
# Each time you run the script, these outputs folders and data will be overwritten.
# Therefore change dir.outputs name every time you want to save (not overwrite) previous runs... 

# Environment settings ----------------------------------------------------

# Working directory
# wd<-"/home/elic/fclapuyt/tracking-gps/"
wd<-"/home/franz/Documents/work/repos/tracking-visitors"

# Subfolders of working directory
dir.outputs<-"outputs"
dir.inputs<-"inputs"
# Subfolders of inputs directory
dir.tracks<-"tracks" 
dir.aoi<-"aoi"

# Variables ---------------------------------------------------------------

# Name of study site
area.name<-"malia"

# Delete outputs directory before running analysis
delete.outputs<-FALSE

# Extension of tracks in the input folder (.gpx or .shp)
track.ext<-".shp"
track.seg.id<-"track_se_1"

# Name of the shapefile containing the limits of the area
limits.area<-"limits-area"
# EPSG number of local projection coordinate system (UTM) 
id.proj<-"+init=epsg:32635"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Normally nothing to change from here... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# # Check package availability
# packages.needed<-c("rgdal", "plotKML", "chron", "spatstat", "maptools", "raster")
# if (length(setdiff(packages.needed, rownames(installed.packages()))) > 0) {
#  install.packages(setdiff(packages.needed, rownames(installed.packages())))
# }

# Load packages -----------------------------------------------------------
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
library(tidyr)

# Working directory and setup ---------------------------------------------

# Define working directory
setwd(wd)

# Define new directories
dir.src<-"src"
dir.functions<-"functions"
dir.graphs<-"graphs" # subfolder of outputs containing graphs of speed and elevation against time
dir.shp<-"shp" # subfolder of outputs containing all spatial data
dir.lines<-"lines" # tracks as line shp
dir.points<-"points" # tracks as point shp
dir.movement<-"movement" # will contain tracks with points which have a speed higher than threshold (lower threshold)
dir.meandir<-"mean-dir"
dir.tables<-"tables"
dir.stops<-"stops"
dir.centroids<-"centroids"
dir.full<-"full"

dir.outputs<-paste0(area.name, '-', dir.outputs)

# Delete outputs directory
if(delete.outputs == TRUE){
  
  unlink(file.path(wd, dir.outputs), recursive = TRUE)
  
}

# Create missing directories
# if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir))){
#   dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir))
# }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script Count visitors -------------------------------------------------

#list of gpx files in the working directory defined as input
# tracks.files<-list.files(file.path(wd, dir.inputs, dir.tracks), pattern=track.ext)

# Read AOI file for mean direction computation
aoi.spdf<-readOGR(dsn = file.path(wd, dir.inputs, dir.aoi, "aoi_for_mean_dir.shp"), stringsAsFactors = FALSE)
aoi.spdf<-spTransform(aoi.spdf, id.proj)

# Read shapefile containing points data for all tracks
track.spdf.full<-readOGR(dsn = file.path(wd, dir.outputs, dir.shp, dir.full, "tracks-full.shp"), stringsAsFactors = FALSE)

# List of gpx files in the working directory defined as input
lines.files<-list.files(file.path(wd,dir.outputs,dir.shp,dir.lines),pattern=".shp")
ids.tracks<-unique(track.spdf.full@data$track)
# ids.tracks<-unname(sapply(lines.files, function(x) as.numeric(unlist(strsplit(x, split="-"))[2])))
# points.files<-list.files(file.path(wd,dir.outputs,dir.shp,dir.points),pattern=".shp")

# Loop within aoi files
aoi.id<-1
# for(aoi.id in aoi.spdf@data$Id){
  
  # # Current AOI
  # aoi.id<-as.numeric(aoi.id)
  # aoi.cur<-aoi.spdf[aoi.id,]
  # 
  # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # # Select points within current aoi
  # pts.cur<-track.spdf.full[aoi.cur,]

  track.id<-16
  # for(track.id in ids.tracks){
  
    # Read lines and points for current track 
    gpx.lines<-spTransform(readOGR(file.path(wd,dir.outputs,dir.shp,dir.lines,paste0("lines-track-",track.id,".shp"))), id.proj)
    gpx.points<-spTransform(readOGR(file.path(wd,dir.outputs,dir.shp,dir.points,paste0("pts-track-",track.id,".shp"))), id.proj)
    
    # Create intersection between track and museum doors (as lines)  
    ints.points<-gIntersection(gpx.lines,ins.spdf, byid=TRUE, checkValidity=FALSE)
    
    # Create attribute table and create spatial point data frame of intersections between tracks and entrances
    ints.attr<-data.frame(id=seq(1,length(ints.points),1))
    ints.spdf<-SpatialPointsDataFrame(ints.points,ints.attr)
    writeOGR(ints.spdf, dsn = file.path(wd,dir.outputs,dir.shp,dir.museum,paste0("intersections-",track.id,"-points.shp")), layer="intersections", driver = "ESRI Shapefile")
    
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
   
# }
  





