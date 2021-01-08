# Reset working directory
rm(list=ls())

# README ------------------------------------------------------------------

# The script "process_tracks" reads raw data from visitors tracks. It computes
# for each point of the track : speed, direction, time, smoothed speed,... 
# If a lower speed threshold is set, each track is divided in "movement" and
# "non-movement" segments, to allow further analyses. 
# Main outputs are: 
#   - a table containing the summary for each track (summary-tracks.csv)
#   - tracks shapefiles as points and lines.
#   - point shapefiles of movement and non-movement segments.
#   - centroids for each non-movement segment for each track. 
#   - graphs of visitor's speed along the track. 

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

# Smooth elevation, speed and direction data
smooth.ele<-0.05
smooth.speed<-0.02
smooth.direction<-20

# Filter speed greather than a threshold ? 
# Points which have a greater speed will be discarded from the track and further analyses 
# because considered to be GPS errors as a pedestrian can not walk faster than the upper threshold
speed.threshold.upper.apply<-FALSE
speed.threshold.upper<-10
# Filter speed lower than a threshold ? 
# Points which have a lower speed than the threshold will be considered as "no movement" and 
# handled separately in the summary
speed.threshold.lower.apply<-TRUE
speed.threshold.lower<-1
# Filter data based on a time threshold ? 
# Apply time threshold to remove wrong start of GPS
time.threshold.apply<-TRUE
time.threshold.begin<-100

# Compute some facultative outputs?
create.graphs<-TRUE
create.lines<-TRUE


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
dir.stops<-"stops"
dir.centroids<-"centroids"
dir.full<-"full"
dir.speed<-"speed"
dir.outputs<-paste0(dir.outputs, '-', area.name)

# Delete outputs directory
if(delete.outputs == TRUE){
  
  unlink(file.path(wd, dir.outputs), recursive = TRUE)
  
}

# Create missing directories
if(!dir.exists(file.path(wd, dir.outputs))){
  dir.create(file.path(wd, dir.outputs)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp))){
  dir.create(file.path(wd, dir.outputs, dir.shp)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.graphs))){
  dir.create(file.path(wd, dir.outputs, dir.graphs)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.graphs, dir.speed))){
  dir.create(file.path(wd, dir.outputs, dir.graphs, dir.speed)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.lines))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.lines)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.points))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.points)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.movement))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.movement)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.stops))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.stops)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.centroids))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.centroids)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.full))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.full)) 
}

# Load custom functions
# source(file.path(wd, dir.src, "distLatLon.R"))
source(file.path(wd, dir.src, dir.functions, "shiftVector.R"))

# Script GPS Tracking -----------------------------------------------------

#list of gpx files in the working directory defined as input
tracks.files<-list.files(file.path(wd, dir.inputs, dir.tracks), pattern=track.ext)

# Create empty objects to collect results
tracks.spdf.list<-list()
summary.list<-list()
tracks.df.list<-list()
centroids.stops.all<-data.frame()

# # Read AOI shapefile
# aoi.spdf<-readOGR(dsn=file.path(wd, dir.inputs, dir.aoi, "AOI_for_LMD_poly-utm.shp"), stringsAsFactors=FALSE)
# aoi.spdf<-spTransform(aoi.spdf, CRS(id.proj))

# tracks.files<-tracks.files[1:2]
# track.file<-tracks.files[1]
for(track.file in tracks.files){
  
  # track name
  track.name<-unlist(strsplit(track.file, split="\\."))[1]
  # track.name<-gsub(pattern="Trac?", replacement="Track", track.name)
  # track.name<-gsub(pattern=" ", replacement="-", track.name)
  
  # Get id of track
  track.id<-unlist(strsplit(track.name, split="_"))[4]
  # track.id<-as.numeric(unlist(strsplit(track.name, split="_"))[4])
  
  # Read GPX file and project in UTM coordinates
  track.spdf<-readOGR(file.path(wd, dir.inputs, dir.tracks, track.file), stringsAsFactors = FALSE)
  track.spdf<-spTransform(track.spdf, CRS(id.proj))
  
  # Clean attribute table and join coordinates of points
  track.spdf@data<-cbind(track.spdf@coords, track.spdf@data[, c(track.seg.id, "ele", "time")])
  colnames(track.spdf@data)[1:2]<-c("lon", "lat")
  track.spdf@data$track_se_1<-as.numeric(track.spdf@data$track_se_1)
  # Add filename to points
  track.spdf@data$track.name<-track.name
  track.spdf@data$track.id<-track.id
  
  # Shift vectors for lat and lon so that each row also contains the next position.
  track.spdf@data$lat.p0<-shiftVector(track.spdf@data$lat, 1)
  track.spdf@data$lon.p0<-shiftVector(track.spdf@data$lon, 1)
  # Compute distances (in metres) between each point and the next one, and cumulated distance
  track.spdf@data$dist<-pointDistance(p1=track.spdf@data[, c("lon", "lat")], p2=track.spdf@data[, c("lon.p0", "lat.p0")], lonlat=FALSE)
  
  # Convert time in to POSIXct format
  track.spdf@data$time<-as.POSIXct(strptime(track.spdf@data$time, format="%Y/%m/%d %H:%M:%S"), format="%Y/%m/%d %H:%M:%S")
  
  # Shift time so that each row contains the timestamp of the next point
  track.spdf@data$time.p0<-shiftVector(track.spdf@data$time, 1)
  # Convert to readable time
  track.spdf@data$time.p0<-as.POSIXct(track.spdf@data$time.p0, origin = "1970-01-01")
  # Compute the number of seconds between each point and the next one.
  track.spdf@data$time.diff<-as.numeric(difftime(track.spdf@data$time, track.spdf@data$time.p0))
  
  # Remove first point
  track.spdf<-track.spdf[c(2:length(track.spdf)),]
  
  # Compute cumulated distance
  track.spdf@data$dist.cum<-cumsum(track.spdf@data$dist)

  # Compute speed and smoothed speed
  track.spdf@data$speed.ms<-track.spdf@data$dist/track.spdf@data$time.diff
  track.spdf@data$speed.kmh<-track.spdf@data$speed.ms*3.6
  track.spdf@data$speed.lowess<-lowess(track.spdf@data$speed.kmh, f=smooth.speed)$y
  
  # Threshold speed for high values that are not realistic
  if(speed.threshold.upper.apply == TRUE){
  
    track.spdf<-track.spdf[which(track.spdf@data$speed.kmh < speed.threshold.upper),]
    
    # Recompute speed and distance --------------------------------------------
    
    # Shift vectors for lat and lon so that each row also contains the next position.
    track.spdf@data$lat.p0<-shiftVector(track.spdf@data$lat, 1)
    track.spdf@data$lon.p0<-shiftVector(track.spdf@data$lon, 1)
    # Compute distances (in metres) between each point and the next one, and cumulated distance
    track.spdf@data$dist<-pointDistance(p1=track.spdf@data[, c("lon", "lat")], p2=track.spdf@data[, c("lon.p0", "lat.p0")], lonlat=FALSE)

    # Shift time so that each row contains the timestamp of the next point
    track.spdf@data$time.p0<-shiftVector(track.spdf@data$time, 1)
    # Convert to readable time
    track.spdf@data$time.p0<-as.POSIXct(track.spdf@data$time.p0, origin = "1970-01-01")
    # Compute the number of seconds between each point and the next one.
    track.spdf@data$time.diff<-as.numeric(difftime(track.spdf@data$time, track.spdf@data$time.p0))
    
    # Remove first point
    track.spdf<-track.spdf[c(2:length(track.spdf)),]
    
    # Compute cumulated distance
    track.spdf@data$dist.cum<-cumsum(track.spdf@data$dist)
    
    # Compute speed and smoothed speed
    track.spdf@data$speed.ms<-track.spdf@data$dist/track.spdf@data$time.diff
    track.spdf@data$speed.kmh<-track.spdf@data$speed.ms*3.6
    track.spdf@data$speed.lowess<-lowess(track.spdf@data$speed.kmh, f=smooth.speed)$y
  
  }
  
  # Compute smoothed elevation and total climb
  track.spdf@data$ele.lowess<-lowess(track.spdf@data$ele, f=smooth.ele)$y
  # Compute elevation difference between each point
  # track.spdf@data$ele.p0<-shiftVector(track.spdf@data$ele, 1)
  # track.spdf@data$climb<-track.spdf@data$ele-track.spdf@data$ele.p0
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Displacement vectors ----------------------------------------------------
  
  # Compute vector components of displacements between points
  track.spdf@data$vec.x<-track.spdf@data$lon-track.spdf@data$lon.p0
  track.spdf@data$vec.y<-track.spdf@data$lat-track.spdf@data$lat.p0
  track.spdf@data$azimut.rad<-atan2(track.spdf@data$vec.y, track.spdf@data$vec.x)
  track.spdf@data$azimut.deg<-atan2(track.spdf@data$vec.y, track.spdf@data$vec.x)*(180/pi)
  
  track.spdf@data$azimut.smooth<-rollmean(x=track.spdf@data$azimut.deg, k=smooth.direction, fill=NA)
  track.spdf@data$id.sort<-seq(1, nrow(track.spdf@data), 1)
  
  
  # Track "cleaning" --------------------------------------------------------
  
  # Remove points with NA values (missing time.diff)
  track.spdf<-track.spdf[!is.na(track.spdf@data$time.diff),]
  # track.spdf<-track.spdf[c(1:(nrow(track.spdf)-2)),]
  
  # Apply time threshold to remove wrong start of GPS
  if(time.threshold.apply == TRUE){
  
    # Find points which have a time difference greater than threshold
    time.id<-which(track.spdf@data$time.diff > time.threshold.begin)
    
    # If points are found, process them
    if(length(time.id) > 0){
    
      if(time.id[1] < 100){
      
      track.spdf<-track.spdf[c((time.id[1]+1):nrow(track.spdf)),]
      
      }
    
    }
  
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Movement analysis with lower speed threshold ----------------------------
  
  if(speed.threshold.lower.apply == TRUE){
  
    # Extract points which are characterised by a speed greater than the lower speed threshold
    track.spdf.movement<-track.spdf[which(track.spdf@data$speed.lowess > speed.threshold.lower),]
    # Save current track to POINT shapefile
    writeOGR(obj=track.spdf.movement, dsn=file.path(wd, dir.outputs, dir.shp, dir.movement), layer=paste0("pts-mvts-", track.id), driver="ESRI Shapefile", overwrite_layer=TRUE)
   
    # Extract points which are characterised by a speed lower than the lower speed threshold
    track.spdf.stops<-track.spdf[which(track.spdf@data$speed.lowess <= speed.threshold.lower),]
    # Save current track to POINT shapefile
    writeOGR(obj=track.spdf.stops, dsn=file.path(wd, dir.outputs, dir.shp, dir.stops, paste0("stops-", track.id)), layer=paste0("pts-stops-", track.id, "-full"), driver="ESRI Shapefile", overwrite_layer=TRUE)
    
    # In the main track file, add a "movement" column (movement == 1) for further analysis
    track.spdf@data$movement<-0
    track.spdf@data[which(track.spdf@data$speed.lowess > speed.threshold.lower), "movement"]<-1

    # Iterate along track to detect moments of movement and non-movement based on lower.speed threshold
    time<-0
    dist<-0
    moments<-list()
    stops.ids<-list()
    stop.ids<-c()
    speed.lowess<-c()

    # for(pt.id in track.spdf@data$id.sort[-length(track.spdf)]){
    for(pt.id in track.spdf@data$id.sort){

      # Get current and next points in the track
      pt.cur<-track.spdf@data[which(track.spdf@data$id.sort == pt.id),]
      pt.next<-track.spdf@data[which(track.spdf@data$id.sort == pt.id + 1),]
      
      # Compute cumulative time and distance for segment
      time<-time + pt.cur$time.diff
      dist<-dist + pt.cur$dist
      
      # Track speed
      speed.lowess<-c(speed.lowess, pt.cur$speed.lowess)
      
      if(pt.cur$movement == 0){

        stop.ids<-c(stop.ids, pt.cur$id.sort)

      }
      
      # Check if pt.cur is the end of a segment, i.e. from non-movement to movement and vice-versa.
      if(nrow(pt.next) != 0){
        
        if(pt.cur$movement != pt.next$movement){
        
          # Update moments list with data of the segment
          moments<-list.append(moments, c(time = time, dist = dist, 
                                          movement = pt.cur$movement,
                                          speed.lowess = mean(speed.lowess)))
          if(!is.null(stop.ids)){
            stops.ids<-list.append(stops.ids, stop.ids)
          }
          # Reset time and distance for next segment
          time<-0
          dist<-0
          stop.ids<-c()
          speed.lowess<-c()
        
        } 
        
      }else{
        
        # Update moments list with data of the segment
        moments<-list.append(moments, c(time = time, dist = dist, 
                                        movement = pt.cur$movement,
                                        speed.lowess = mean(speed.lowess)))
        
      }
      
    }
  
    if(length(moments) > 0){
    
      # Handle moments list and convert as dataframe. 
      moments<-t(data.frame(moments))
      rownames(moments)<-seq(1, nrow(moments), 1)
      moments<-as.data.frame(moments)
      
    }else{
      moments<-data.frame()
    }
  
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Write points output -----------------------------------------------------
  
  # Plot vector field for visitor's displacements
  # track.vec<-track.spdf[seq(from=(smooth.direction/2), to=nrow(track.spdf)-(smooth.direction/2), by=8),]
  # vectorField(xpos=track.vec@data$lon, ypos=track.vec@data$lat, u=track.vec@data$azimut.rad, v=1, vecspec="rad")
  
  # Save current gpx to POINT shapefile
  writeOGR(obj=track.spdf, dsn=file.path(wd, dir.outputs, dir.shp, dir.points, paste0("pts-track-",track.id,".shp")), layer=track.name, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  # Summary metrics ----------------------------------------------------------
  
  # Compute summary statistics and store in dataframe
  time.start<-track.spdf@data[1, "time.p0"]
  # duration.total<-as.numeric(difftime(track.spdf@data[nrow(track.spdf@data), "time"], time.start, units="mins"))
  duration.total<-sum(track.spdf@data$time.diff)/60
  dist.total<-sum(track.spdf@data$dist)
  speed.kmh.mean<-mean(track.spdf@data$speed.kmh)
  speed.lowess.mean<-mean(track.spdf@data$speed.lowess, na.rm=TRUE)
  # ele.mean<-mean(track.spdf@data$ele)
  # ele.lowess.mean<-mean(track.spdf@data$ele.lowess)
  # climb.total<-sum(track.spdf@data$climb[track.spdf@data$climb>=0])
  
  # Compute summary statistics for movements sections
  if(speed.threshold.lower.apply == TRUE & length(moments) > 0){
  
    # Compute movements stats
    movements<-moments[which(moments$movement == 1),]
    speed.mean.movement<-mean(movements$dist/movements$time*3.6)
    duration.movement<-sum(movements$time/60)
    dist.movement<-sum(movements$dist)
    speed.mean.movement.lowess<-mean(movements$speed.lowess)
    # Compute non-movements stats
    stops<-moments[which(moments$movement == 0),]
    speed.mean.stops<-mean(stops$dist/stops$time*3.6)
    duration.stops<-sum(stops$time/60)
    dist.stops<-sum(stops$dist)
    speed.mean.stops.lowess<-mean(stops$speed.lowess)
    stops.time.summary<-summary(moments[which(moments$movement == 0), "time"])
    stops.dist.summary<-summary(moments[which(moments$movement == 0), "dist"])
    stops.count<-nrow(stops)  
  
  }else{
  
    speed.lowess.mean.movement<-0
    duration.movement<-0
    dist.movement<-0
    stops.time.summary<-summary(0)
    stops.dist.summary<-summary(0)
  
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute centroids for non-movement sections -----------------------------
  
  if(!is.empty(stops.ids)){
  
    # Compute shapefile for each stop and centroid. 
    centroids<-data.frame()
    # Loop through stops
    # stop.id<-1
    for(stop.id in 1:length(stops.ids)){
    
      # Extract track part for current stop
      stop<-track.spdf[which(track.spdf@data$id.sort %in% stops.ids[[stop.id]]),]
      # Save stop as shapefile
      writeOGR(obj = stop, dsn = file.path(wd, dir.outputs, dir.shp, dir.stops, paste0("stops-", track.id), paste0("pts-stops-", track.id, "-stop-", stop.id, ".shp")), 
      layer=paste0(track.name, "-stop-", stop.id), overwrite_layer = TRUE, driver="ESRI Shapefile")
      
      # Compute centroid of current stop and collect in dataframe
      coords<-cbind(lon=mean(stop@coords[, 1]), lat=mean(stop@coords[, 2]), 
                    id=stop.id, duration = sum(stop@data$time.diff), 
                    var.x = sd(stop@coords[, 1]), var.y = sd(stop@coords[, 2]))
      centroids<-rbind(centroids, coords)
  
    }
    
    # Add track ID to centroids dataframe
    centroids<-cbind(centroids, track=track.id, stringsAsFactors=FALSE)
    
    # Gather centroids of stops for current track in a dataframe
    centroids.stops.all<-rbind(centroids.stops.all, centroids)
    # Convert centroids of stops for current track to spatial data frame
    coordinates(centroids) = ~ lon+lat
    projection(centroids) = CRS(id.proj)
    # Save centroids of track as shapefile
    writeOGR(obj = centroids, dsn = file.path(wd, dir.outputs, dir.shp, dir.centroids, paste0("stops-", track.id, "-centroids.shp")), 
    layer = paste0(track.name, "-stops-centroids"), overwrite_layer = TRUE, driver = "ESRI Shapefile")
  
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Save summary of current track in a list
  summary.list<-list.append(summary.list, list(name=track.name, 
                                               duration=duration.total, 
                                               distance=dist.total, 
                                               speed.kmh.mean=speed.kmh.mean, 
                                               speed.lowess.mean=speed.lowess.mean, 
                                               duration.movement=duration.movement, 
                                               dist.movement=dist.movement, 
                                               speed.mean.movement=speed.mean.movement,
                                               speed.mean.movement.lowess=speed.mean.movement.lowess,
                                               duration.stops=duration.stops,
                                               dist.stops=dist.stops,
                                               speed.mean.stops=speed.mean.stops,
                                               speed.mean.stops.lowess=speed.mean.stops.lowess,
                                               stops.count=stops.count,
                                               stops.dist.summary=stops.dist.summary, 
                                               stops.time.summary=stops.time.summary
                                               # ele.mean=ele.mean, ele.lowess.mean=ele.lowess.mean, climb.total=climb.total
                                               ))
  
  # Collect track.spdf in a list for further processing
  tracks.spdf.list<-list.append(tracks.spdf.list, track.spdf)
  tracks.df.list<-list.append(tracks.df.list, as.data.frame(track.spdf))
  
  # Graphs ------------------------------------------------------------------
  
  if(create.graphs){
  
    # Plot elevation against time and save in output folder
    plot.elev.lowess<-ggplot(data=track.spdf@data, aes(x=time, y=ele))+
    geom_point(size=1, color="darkgrey")+
    geom_point(aes(y=ele.lowess), size=1, color="red")+
    labs(x="Time", y="Elevation (m)", 
    title=paste0("Start time: ", time.start, " / Duration: ", round(duration.total, 2), " min"))
    # print(plot.elev.lowess)
    # ggsave(plot.elev.lowess, filename=file.path(wd, dir.outputs, dir.graphs, paste0("elevation-vs-time-", track.name, ".png")), device="png", dpi=300, units="cm", width=20, height=15)
    
    # Plot speed against time and save in output folder
    plot.speed<-ggplot(data=track.spdf@data, aes(x=time, y=speed.kmh))+
    geom_point(size=1, color="darkgrey")+
    geom_point(aes(y=speed.lowess), size=1, color="red")+
    labs(x="Time", y="Speed (km/h)", 
    title=paste0("Start time: ", time.start, " / Duration: ", round(duration.total, 2), " min"))
    # print(plot.speed)
    ggsave(plot.speed, filename=file.path(wd, dir.outputs, dir.graphs, dir.speed, paste0("speed-vs-time-track-", track.id, ".png")), device="png", dpi=300, units="cm", width=20, height=15)
  
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Line datasets -----------------------------------------------------------
  
  if(create.lines){
  
    # Create line shapefile from track
    
    # How many points between line segments (smoothing factor)
    lines.interval<-1
    
    # Loop within gpx to create vertices of lines
    lines.data<-data.frame()
    for(point in seq(1, nrow(track.spdf@data), lines.interval)){
      
      line.temp<-rbind(track.spdf@data[point,], track.spdf@data[point+lines.interval,])
      line.temp$id<-point
      
      lines.data<-rbind(lines.data, line.temp)
    
    }
    # Remove NA lines
    lines.data<-lines.data[which(lines.data$id != point),]
    
    # Convert lines.data into SpatialLinesDataframes
    d<-lines.data
    coordinates(d)<-~lon+lat
    x <- lapply(split(d, d$id), function(x) Lines(list(Line(coordinates(x))), x$id[1L]))
    
    lines <- SpatialLines(x)
    data <- data.frame(id = unique(d$id))
    rownames(data) <- data$id
    track.spdf.line <- SpatialLinesDataFrame(lines, data)
    proj4string(track.spdf.line)<-CRS(id.proj)
    
    # track.spdf.line<-SpatialLinesDataFrame(SpatialLines(list(Lines(list(Line(coordinates(track.spdf))), "id"))), 
    #                   data=data.frame(x=1, y=1, row.names = "id"))
    
    # Save current gpx to LINE shapefile
    writeOGR(obj=track.spdf.line, dsn=file.path(wd, dir.outputs, dir.shp, dir.lines, paste0("lines-track-", track.id, ".shp")), layer=paste0(track.name, "-lines"), driver="ESRI Shapefile", overwrite_layer=TRUE)
    
  }

}

# Process summary of tracks -----------------------------------------------
# Convert summary list to dataframe
summary<-data.frame(track.name=sapply(summary.list, "[[", "name"), 
                    duration=sapply(summary.list, "[[", "duration"), 
                    distance=sapply(summary.list, "[[", "distance"), 
                    speed.kmh.mean=sapply(summary.list, "[[", "speed.kmh.mean"), 
                    speed.lowess.mean=sapply(summary.list, "[[", "speed.lowess.mean"), 
                    duration.movement=sapply(summary.list, "[[", "duration.movement"), 
                    dist.movement=sapply(summary.list, "[[", "dist.movement"), 
                    speed.mean.movement=sapply(summary.list, "[[", "speed.mean.movement"),
                    speed.mean.movement.lowess=sapply(summary.list, "[[", "speed.mean.movement.lowess"),
                    duration.stops=sapply(summary.list, "[[", "duration.stops"), 
                    dist.stops=sapply(summary.list, "[[", "dist.stops"), 
                    speed.mean.stops=sapply(summary.list, "[[", "speed.mean.stops"),
                    speed.mean.stops.lowess=sapply(summary.list, "[[", "speed.mean.stops.lowess")
                    # ele.mean=sapply(summary.list, "[[", "ele.mean"), 
                    # ele.lowess.mean=sapply(summary.list, "[[", "ele.lowess.mean"), 
                    # climb.total=sapply(summary.list, "[[", "climb.total")
                    )

# Write summary in csv
write.csv2(x=summary, file=file.path(wd, dir.outputs, "summary-tracks.csv"))

# Process merged tracks ---------------------------------------------------
# Merge all gpx files in a data.frame and convert to spdf
gpx.df.full<-do.call("rbind", tracks.df.list)
track.spdf.full<-SpatialPointsDataFrame(gpx.df.full[, c("lon", "lat")], gpx.df.full)
proj4string(track.spdf.full)<-CRS(id.proj)
# Save all gpx to shapefile
writeOGR(obj=track.spdf.full, dsn=file.path(wd, dir.outputs, dir.shp, dir.full, "tracks-full.shp"), layer="tracks-full", driver="ESRI Shapefile", overwrite_layer=TRUE)

# Convert centroids of stops for all tracks to spatial data frame
coordinates(centroids.stops.all) = ~ lon+lat
projection(centroids.stops.all) = CRS(id.proj)
# Save centroids of track as shapefile
writeOGR(obj = centroids.stops.all, dsn = file.path(wd, dir.outputs, dir.shp, dir.centroids, paste0("stops-centroids-full.shp")), 
     layer = "stops-centroids-full", overwrite_layer = TRUE, driver = "ESRI Shapefile")

# Speed threshold ---------------------------------------------------------
# Extract points which are characterised by a speed greater than the speed threshold
track.spdf.full.stops<-track.spdf.full[which(track.spdf.full@data$speed.kmh > speed.threshold.lower),]
# Save speed thresholded gpx file to shapefile
writeOGR(obj=track.spdf.full.stops, dsn=file.path(wd, dir.outputs, dir.shp, dir.full, "tracks-full-movement.shp"), layer="tracks-full-movement", driver="ESRI Shapefile", overwrite_layer=TRUE)
# Extract points which are characterised by a speed greater than the speed threshold
track.spdf.full.stops<-track.spdf.full[which(track.spdf.full@data$speed.kmh <= speed.threshold.lower),]
# Save speed thresholded gpx file to shapefile
writeOGR(obj=track.spdf.full.stops, dsn=file.path(wd, dir.outputs, dir.shp, dir.full, "tracks-full-stops.shp"), layer="tracks-full-stops", driver="ESRI Shapefile", overwrite_layer=TRUE)


# Histogram all speeds ----------------------------------------------------
# Gather data of speed from all tracks
speed.full<-track.spdf.full@data[, c("speed.kmh", "speed.lowess")]
speed.full<-stack(speed.full)
# Create histogram
hist.speed<-ggplot(data = speed.full, aes(x = values, color = ind, fill = ind))+
  geom_histogram(position = "identity", binwidth = 0.5, alpha=0.5)+
  labs(title = "Distribution of visitors' speed of displacement", x = "Speed (km/h)", y = "Count")
hist.speed
# Save to png file in the output directory
ggsave(filename = file.path(wd, dir.outputs, "histogram-speed-full.png"))


# Correlation between distance and speed per track ------------------------

# Duration | Total distance in movement
cor.duration.distmvt<-ggplot(data = summary, aes(x = duration, y = dist.movement))+
  geom_point()+
  theme_bw()+
  labs(title = "Correlation between total duration of visit and distance of movement",
       x = "Total duration (min)",
       y = "Distance during movement (m)")
cor.duration.distmvt

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, "corr-duration-distmvt.png"), 
       plot = cor.duration.distmvt, device = "png")

# Duration | Total distance
cor.duration.disttot<-ggplot(data = summary, aes(x = duration, y = distance))+
  geom_point()+
  theme_bw()+
  labs(title = "Correlation between total duration of visit and total distance",
       x = "Total duration (minutes)",
       y = "Total distance (m)")
cor.duration.disttot

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, "corr-duration-disttot.png"), 
       plot = cor.duration.disttot, device = "png")

# Duration | Smoothed speed during movement phases
cor.duration.speedmvt<-ggplot(data = summary, aes(x = duration, y = speed.mean.movement.lowess))+
  geom_point()+
  theme_bw()+
  labs(title = "Correlation between total duration of visit and smoothed displacement speed",
       x = "Total duration (minutes)",
       y = "Smoothed displacement speed (km/h)")
cor.duration.speedmvt

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, "corr-duration-speedmvt.png"), 
       plot = cor.duration.speedmvt, device = "png")

# Duration | Smoothed speed during non-movement phases
cor.duration.speedstops<-ggplot(data = summary, aes(x = duration, y = speed.mean.stops.lowess))+
  geom_point()+
  theme_bw()+
  labs(title = "Correlation between total duration of visit and smoothed speed during \"stops\"",
       x = "Total duration (minutes)",
       y = "Smoothed speed during \"stops\" (km/h)")
cor.duration.speedstops

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, "corr-duration-speedstops.png"), 
       plot = cor.duration.speedstops, device = "png")


