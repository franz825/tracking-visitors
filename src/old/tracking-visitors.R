# Reset working directory
rm(list=ls())

# README ------------------------------------------------------------------

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

# Delete outputs directory before running analysis
delete.outputs<-TRUE

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
dir.meandir<-"aoi" 
dir.tables<-"tables"
dir.museum<-"museum"
dir.stops<-"stops"
dir.centroids<-"centroids"

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
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir)) 
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
# if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.museum))){
#   dir.create(file.path(wd, dir.outputs, dir.shp, dir.museum)) 
# }
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.stops))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.stops)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.centroids))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.centroids)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables)) 
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.shp))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.shp)) 
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

# track.file<-tracks.files[4]
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
    writeOGR(obj=track.spdf.movement, dsn=file.path(wd, dir.outputs, dir.shp, dir.movement, paste0(track.name, "-movement")), layer=track.name, driver="ESRI Shapefile", overwrite_layer=TRUE)
   
    # Extract points which are characterised by a speed lower than the lower speed threshold
    track.spdf.stops<-track.spdf[which(track.spdf@data$speed.lowess <= speed.threshold.lower),]
    # Save current track to POINT shapefile
    writeOGR(obj=track.spdf.stops, dsn=file.path(wd, dir.outputs, dir.shp, dir.stops, paste0(track.name, "-stops")), layer=track.name, driver="ESRI Shapefile", overwrite_layer=TRUE)
    
    # In the main track file, add a "movement" column (movement == 1) for further analysis
    track.spdf@data$movement<-0
    track.spdf@data[which(track.spdf@data$speed.lowess > speed.threshold.lower), "movement"]<-1

    # Iterate along track to detect moments of movement and non-movement based on lower.speed threshold
    time<-0
    dist<-0
    moments<-list()
    stops.ids<-list()
    stop.ids<-c()

    # for(pt.id in track.spdf@data$id.sort[-length(track.spdf)]){
    for(pt.id in track.spdf@data$id.sort){

      # Get current and next points in the track
      pt.cur<-track.spdf@data[which(track.spdf@data$id.sort == pt.id),]
      pt.next<-track.spdf@data[which(track.spdf@data$id.sort == pt.id + 1),]
      
      # Compute cumulative time and distance for segment
      time<-time + pt.cur$time.diff
      dist<-dist + pt.cur$dist
      
      if(pt.cur$movement == 0){

        stop.ids<-c(stop.ids, pt.cur$id.sort)

      }
      
      # Check if pt.cur is the end of a segment, i.e. from non-movement to movement and vice-versa.
      if(nrow(pt.next) != 0){
        
        if(pt.cur$movement != pt.next$movement){
        
          # Update moments list with data of the segment
          moments<-list.append(moments, c(time = time, dist = dist, movement = pt.cur$movement))
          if(!is.null(stop.ids)){
            stops.ids<-list.append(stops.ids, stop.ids)
          }
          # Reset time and distance for next segment
          time<-0
          dist<-0
          stop.ids<-c()
        
        } 
        
      }else{
        
        # Update moments list with data of the segment
        moments<-list.append(moments, c(time = time, dist = dist, movement = pt.cur$movement))
        
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
  writeOGR(obj=track.spdf, dsn=file.path(wd, dir.outputs, dir.shp, dir.points, paste0(track.name, "-points.shp")), layer=track.name, driver="ESRI Shapefile", overwrite_layer=TRUE)
  
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
    # Compute non-movements stats
    stops<-moments[which(moments$movement == 0),]
    speed.mean.stops<-mean(stops$dist/stops$time*3.6)
    duration.stops<-sum(stops$time/60)
    dist.stops<-sum(stops$dist)
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
      writeOGR(obj = stop, dsn = file.path(wd, dir.outputs, dir.shp, dir.stops, paste0(track.name, "-stop-", stop.id, ".shp")), 
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
    writeOGR(obj = centroids, dsn = file.path(wd, dir.outputs, dir.shp, dir.centroids, paste0(track.name, "-stops-centroids.shp")), 
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
                                               duration.stops=duration.stops,
                                               dist.stops=dist.stops,
                                               speed.mean.stops=speed.mean.stops,
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
    ggsave(plot.speed, filename=file.path(wd, dir.outputs, dir.graphs, paste0("speed-vs-time-", track.name, ".png")), device="png", dpi=300, units="cm", width=20, height=15)
  
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
    writeOGR(obj=track.spdf.line, dsn=file.path(wd, dir.outputs, dir.shp, dir.lines, paste0(track.name, "-lines.shp")), layer=paste0(track.name, "-lines"), driver="ESRI Shapefile", overwrite_layer=TRUE)
    
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
                    duration.stops=sapply(summary.list, "[[", "duration.stops"), 
                    dist.stops=sapply(summary.list, "[[", "dist.stops"), 
                    speed.mean.stops=sapply(summary.list, "[[", "speed.mean.stops")
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
writeOGR(obj=track.spdf.full, dsn=file.path(wd, dir.outputs, dir.shp, "gpx-full.shp"), layer="gpx-full", driver="ESRI Shapefile", overwrite_layer=TRUE)

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
writeOGR(obj=track.spdf.full.stops, dsn=file.path(wd, dir.outputs, dir.shp, "gpx-full-movement.shp"), layer="gpx-full-movement", driver="ESRI Shapefile", overwrite_layer=TRUE)
# Extract points which are characterised by a speed greater than the speed threshold
track.spdf.full.stops<-track.spdf.full[which(track.spdf.full@data$speed.kmh <= speed.threshold.lower),]
# Save speed thresholded gpx file to shapefile
writeOGR(obj=track.spdf.full.stops, dsn=file.path(wd, dir.outputs, dir.shp, "gpx-full-stops.shp"), layer="gpx-full-stops", driver="ESRI Shapefile", overwrite_layer=TRUE)


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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read AOI file for mean direction computation
aoi.spdf<-readOGR(dsn = file.path(wd, dir.inputs, dir.aoi, "aoi_for_mean_dir.shp"), stringsAsFactors = FALSE)
aoi.spdf<-spTransform(aoi.spdf, id.proj)

# Compute data within aois for full data

# Azimut classes to categorize directions
azimut.classes<-seq(-180, 180, 45)

# List to store results 
mean.dir.aois<-list()

# Two columns to store results in aoi shp
# aoi.spdf@data$tan<-0
aoi.spdf@data$tan2<-0
# aoi.spdf@data$tan.th<-0
aoi.spdf@data$tan2.th<-0

aoi.summary<-data.frame()

# Iterate within aois
# aoi.id<-1
for(aoi.id in aoi.spdf@data$Id){
 
 # Current AOI
 aoi.id<-as.numeric(aoi.id)
 aoi.cur<-aoi.spdf[aoi.id,]
 
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # Select points within current aoi
 pts.cur<-track.spdf.full[aoi.cur,]
 pts.cur.th<-track.spdf.movement[aoi.cur,]
 
 # Compute total x and y components of movement
 sum.x<-sum(pts.cur@data$vec.x)
 sum.y<-sum(pts.cur@data$vec.y)
 sum.x.th<-sum(pts.cur.th@data$vec.x)
 sum.y.th<-sum(pts.cur.th@data$vec.y)

 # Compute mean total movement within aoi
 mean.dir.tan<-atan(sum.y/sum.x)*(180/pi)
 mean.dir.tan2<-atan2(sum.y, sum.x)*(180/pi)
 mean.dir.tan.th<-atan(sum.y.th/sum.x.th)*(180/pi)
 mean.dir.tan2.th<-atan2(sum.y.th, sum.x.th)*(180/pi)
 
 # Store in aoi spdf
 # aoi.spdf@data[aoi.id, "mean.tan"]<-mean.dir.tan
 aoi.spdf@data[aoi.id, "tan2"]<-mean.dir.tan2
 # aoi.spdf@data[aoi.id, "tan.th"]<-mean.dir.tan.th
 aoi.spdf@data[aoi.id, "tan2.th"]<-mean.dir.tan2.th
 
 # Store in list
 mean.dir.aois[[aoi.id]]<-list(mean.dir.tan=mean.dir.tan, 
                                mean.dir.tan2=mean.dir.tan2, 
                                mean.dir.tan.th=mean.dir.tan.th, 
                                mean.dir.tan2.th=mean.dir.tan2.th
                                )
 
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 # Create vectors per AOI
 
 # Create Line dataframe for AOI
aoi.lines = SpatialLinesDataFrame(
  sl=SpatialLines(
    LinesList=list(
        Lines(
          list(
            Line(
            coords=matrix(c(0, 1, 0, 1), 2, 2)
            )
          ), ID=1)
        )
      ), 
    data=data.frame(a=1, b=2, c=3))[-1,]
# Project empty dataframe in utm
proj4string(aoi.lines)<-CRS(id.proj)

# Get unique values from tracks within AOI
tracks.cur<-unique(pts.cur@data$track.name)
id.line<-1
# track.name<-tracks.cur[2]
for(track.name in tracks.cur){

  # Selects points from current track 
  track.cur<-pts.cur[which(pts.cur@data$track.name == track.name),]
  track.cur<-track.cur[c(1, length(track.cur)),]
  
  # Compute direction from first to last point of the track within AOI
  vec.x<-track.cur@data$lon[2]-track.cur@data$lon[1]
  vec.y<-track.cur@data$lat[2]-track.cur@data$lat[1]
  azimut.deg<-atan2(vec.y, vec.x)*(180/pi)
  
  # Loop within gpx to create vertices of lines
  lines.data<-data.frame()
  for(point in seq(1, nrow(track.cur@data), 1)){
  
    # line.id<-line.id+1
    
    line.temp<-rbind(track.cur@data[point,], track.cur@data[point+1,])
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
  data <- data.frame(id = unique(d$id), id2=id.line, track=track.name, azimut=azimut.deg)
  rownames(data) <- data$id
  track.cur.line <- SpatialLinesDataFrame(lines, data)
  proj4string(track.cur.line)<-CRS(id.proj)
  
  # Collect line to aoi.lines
  aoi.lines<-rbind(aoi.lines, track.cur.line)
  id.line<-id.line+1

}

# Extract aoi lines data frame 
# aoi.lines.df<-aoi.lines@data
# aoi.lines.df$id<-NULL
# colnames(aoi.lines.df)[1]<-"id"  

# Names azimuts
azimut.names<-data.frame(class=seq(1, 8, 1), azimut.name=c("OSO", "SSO", "SSE", "ESE", "ENE", "NNE", "NNO", "ONO"))

# Create linear scale bins
aoi.lines@data$azimut.class<-as.vector(cbind(bin=cut(aoi.lines@data$azimut, azimut.classes)))

# Compute azimut frequency per binned direction
aoi.lines.freq<-aoi.lines@data %>%
count(azimut.class)
colnames(aoi.lines.freq)[2]<-"count"
aoi.lines.freq$freq.rel<-aoi.lines.freq$count/sum(aoi.lines.freq$count)
aoi.lines.freq$aoi.id<-aoi.id
# Add azimut name to table
aoi.lines.freq<-merge(x=aoi.lines.freq, y=azimut.names, by.x="azimut.class", by.y="class")

aoi.summary<-rbind(aoi.summary, aoi.lines.freq)

# Save tables of frequencies and aoi.lines
write.csv2(aoi.lines@data, file=file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables, paste0("aoi-", aoi.id, "-vectors.csv")))
write.csv2(aoi.lines.freq, file=file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables, paste0("aoi-", aoi.id, "-azimut-freq.csv")))

# Save current aoi vectors to list
writeOGR(obj=aoi.lines, dsn=file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.shp, paste0("aoi-", aoi.id, "-vectors.shp")), 
          layer=paste0("aoi-", aoi.id, "-vectors"), driver="ESRI Shapefile", overwrite_layer=TRUE)

# Save current aoi vectors to list
# aoi.lines.df$ws<-1
# aoi.lines.df$azimut.north<-0
# 
# aoi.lines.df[which(aoi.lines.df$azimut > 0)]
# pollutionRose(aoi.lines.df, ws="ws", wd="azimut")
# 
# view(pollutionRose)

}

# Save aoi spdf to shapefile
writeOGR(obj=aoi.spdf, dsn=file.path(wd, dir.outputs, dir.shp, "aoi-mean-directions.shp"), layer="aoi-mean-directions", 
          driver="ESRI Shapefile", overwrite_layer=TRUE)

# Save aoi summary to csv
write.csv2(aoi.summary, file=file.path(wd, dir.outputs, "summary-aoi.csv"))


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
area.spdf<-readOGR(dsn=file.path(wd, dir.inputs, list.files(path=file.path(wd, dir.inputs, dir.shp), pattern=".shp")))

# # Divide study area into quadrats
# area.rstr<-raster(area.spdf)
# res(area.rstr)<-100
# area.quads <- as(rasterize(area.spdf, area.rstr), 'SpatialPolygons')
# 
# # Display quads and gpx tracks
# plot(quads)
# points(track.spdf.full, col="red", cex=.5)

# Convert limits to owin type and get coordinates of points
area.owin<-as.owin(area.spdf)
gpx.pts<-coordinates(track.spdf.full)

# Point pattern analysis
gpx.ppp<-ppp(gpx.pts[, 1], gpx.pts[, 2], window=area.owin)
plot(gpx.ppp)
gpx.density<-density(gpx.ppp)
plot(gpx.density)

# Count points within quadrats
pts.count<-quadratcount(gpx.ppp, nx=50, ny=50)
# Plot count
plot(gpx.ppp, pch=".", cols="grey70")
plot(pts.count, add=TRUE)

# Compute density of points within quadrats
pts.intensity<-intensity(pts.count)
plot(intensity(pts.count, image=TRUE))
plot(gpx.ppp, pch=".", add=TRUE, col="green", main="Density of points")

pts.tess<-as.tess(pts.count)
pts.spdf<-as(pts.tess, "SpatialPolygons")

pts.count$tile

plot(pts.spdf)

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

