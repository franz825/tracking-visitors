# Reset working directory
# rm(list=setdiff(ls(), "track.spdf.full"))
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
dir.aoi.time<-"aoi-time"

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

# Filter speed lower than a threshold ? 
speed.threshold.lower.apply<-TRUE

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
# dir.meandir<-"mean-dir"
dir.tables<-"tables"
dir.stops<-"stops"
dir.centroids<-"centroids"
dir.full<-"full"
dir.temp<-"temp"
dir.rdata<-"rdata"
dir.meantime<-"mean-time"

dir.outputs<-paste0(dir.outputs, '-', area.name)

# Delete outputs directory
if(delete.outputs == TRUE){
  
  unlink(file.path(wd, dir.outputs), recursive = TRUE)
  
}

# Create missing directories
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.temp))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.temp))
}
if(!dir.exists(file.path(wd, dir.outputs, dir.rdata))){
  dir.create(file.path(wd, dir.outputs, dir.rdata))
}
if(!dir.exists(file.path(wd, dir.outputs, dir.tables, dir.meantime))){
  dir.create(file.path(wd, dir.outputs, dir.tables, dir.meantime))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Script Mean time in AOI ------------------------------------------------------

# Read AOI files in aoi-time directory
aoi.spdf<-st_as_sf(spTransform(readOGR(dsn = file.path(wd, dir.inputs, dir.aoi.time, "aoi-time-all.shp"), stringsAsFactors = FALSE), id.proj))

# Read shapefile containing points data for all tracks
track.spdf.full<-st_as_sf(readOGR(dsn = file.path(wd, dir.outputs, dir.shp, dir.full, "tracks-full.shp"), stringsAsFactors = FALSE))

# List of gpx files in the working directory defined as input
lines.files<-list.files(file.path(wd,dir.outputs,dir.shp,dir.lines),pattern=".shp")
ids.tracks<-unique(track.spdf.full$track)

# Container for results
summary.aois<-list()

# Loop within aois to compute time spent by each visitor in each aoi
# aoi.id<-aoi.spdf$aoi_id[1]
for(aoi.id in aoi.spdf$aoi_id){
  
  # Current AOI
  aoi.cur<-aoi.spdf[which(aoi.spdf$aoi_id == aoi.id),]
  # Container for results
  summary.aoi<-list()

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Loop within tracks
  # track.id<- 13
  for(track.id in ids.tracks){
  
    # Read lines and points for current track 
    # track.lines<-st_as_sf(spTransform(readOGR(file.path(wd,dir.outputs,dir.shp,dir.lines,paste0("lines-track-",track.id,".shp"))), id.proj))
    track.pts<-st_as_sf(spTransform(readOGR(file.path(wd,dir.outputs,dir.shp,dir.points,paste0("pts-track-",track.id,".shp"))), id.proj))
    
    # int.line<-st_intersection(aoi.cur, track.lines)
    int.pts<-st_intersection(aoi.cur, track.pts)
    
    # st_write(int.line, file.path(wd, dir.outputs, dir.shp, dir.temp, "int-line.shp"), append = FALSE)
    # st_write(int.pts, file.path(wd, dir.outputs, dir.shp, dir.temp, "int-pts.shp"), append = FALSE)
    
    # Iterate along track to detect moments of movement and non-movement based on lower.speed threshold
    time<-0
    dist<-0
    moments<-list()
    stops.ids<-list()
    stop.ids<-c()
    speed.lowess<-c()
    segment.id<-1
    speed.mean.movement<-0
    speed.mean.movement.lowess<-0
    movements<-NULL
    stops<-NULL
    id.start<-as.integer(st_drop_geometry(int.pts[1,"id_sort"]))
    get.id.start = FALSE
    
    # Loop along the track
    for(pt.id in int.pts$id_sort){
      
      if(get.id.start == TRUE){
        id.start<-pt.id
        get.id.start<-FALSE
      }
      
      # Get current and next points in the track
      pt.cur<-int.pts[which(int.pts$id_sort == pt.id),]
      pt.next<-int.pts[which(int.pts$id_sort == pt.id + 1),]
      
      # Compute cumulative time and distance for segment
      time<-time + pt.cur$tim_dff
      dist<-dist + pt.cur$dist
      
      # Track speed
      speed.lowess<-c(speed.lowess, pt.cur$spd_lws)
      
      if(pt.cur$movemnt == 0){
        
        stop.ids<-c(stop.ids, pt.cur$id.sort)
        
      }
      
      # Check if pt.cur is the end of a segment, i.e. from non-movement to movement and vice-versa.
      if(nrow(pt.next) != 0){
        
        if(pt.cur$movemnt != pt.next$movemnt){

          # Compute beginning and ending ids for moment
          id.start<-id.start
          id.end<-pt.cur$id_sort
          
          # Update moments list with data of the segment
          moments<-list.append(moments, c(time = time, dist = dist, 
                                          movement = pt.cur$movemnt,
                                          speed.lowess = mean(speed.lowess),
                                          segment.id = segment.id,
                                          id.start = id.start,
                                          id.end = id.end))
          
          if(!is.null(stop.ids)){
            stops.ids<-list.append(stops.ids, stop.ids)
          }
          # Reset time and distance for next segment
          time<-0
          dist<-0
          stop.ids<-c()
          speed.lowess<-c()
          id.start<-pt.next$id_sort
          
        }
        
      }else{
        
        # Compute beginning and ending ids for moment
        id.start<-id.start
        id.end<-pt.cur$id_sort
        
        # Update moments list with data of the segment
        moments<-list.append(moments, c(time = time, dist = dist, 
                                        movement = pt.cur$movemnt,
                                        speed.lowess = mean(speed.lowess),
                                        segment.id = segment.id,
                                        id.start = id.start,
                                        id.end = id.end))
        
        # Reset time and distance for next segment
        time<-0
        dist<-0
        stop.ids<-c()
        speed.lowess<-c()
        get.id.start<-TRUE
        
        # Update segment value
        segment.id<-segment.id + 1

      }
      
    }
    
    if(length(moments) > 0){
      
      # Handle moments list and convert as dataframe. 
      moments<-t(data.frame(moments))
      rownames(moments)<-seq(1, nrow(moments), 1)
      moments<-as.data.frame(moments)
      moments$time<-moments$time/60
      
    }else{
      moments<-data.frame()
    }
    
    # Compute summary statistics and store in dataframe
    time.start<-int.pts[1, "time_p0"]
    # duration.total<-as.numeric(difftime(int.pts[nrow(int.pts), "time"], time.start, units="mins"))
    duration.total<-sum(int.pts$tim_dff)/60
    dist.total<-sum(int.pts$dist)
    speed.kmh.mean<-mean(int.pts$spd_kmh)
    speed.lowess.mean<-mean(int.pts$spd_lws, na.rm=TRUE)

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
      # Compute number of segments
      segments.nb<-max(movements$segment.id)
      
      # Compute segments characteristics
      segments<-moments %>% group_by(segment.id) %>% summarise(aoi.id=aoi.id,
                                                               aoi.name=aoi.cur$layer,
                                                               track.id=track.id,
                                                               track.name=as.character(track.pts$trck_nm[1]),
                                                               segment.time = sum(time), 
                                                               segment.dist = sum(dist),
                                                               stops.nb = sum(movement == 0),
                                                               id.start = min(id.start),
                                                               id.end = max(id.end))
      
      # Compute timestamps for segments
      int.pts.mutated <- st_drop_geometry(int.pts) %>% mutate_if(is.factor, as.character)
      segments$time.start<-int.pts.mutated[which(int.pts.mutated$id_sort %in% segments$id.start),"time"]
      segments$time.end<-int.pts.mutated[which(int.pts.mutated$id_sort %in% segments$id.end),"time"]
      
      
    
    }else{
      
      speed.lowess.mean.movement<-0
      duration.movement<-0
      dist.movement<-0
      stops.time.summary<-summary(0)
      stops.dist.summary<-summary(0)
      segments<-0
      segments.nb<-0
      
    }
    
    # timestamp<-as.POSIXct(strptime(as.character(int.pts[1,]$time), format="%Y-%m-%d %H:%M:%S"))
    
    # Gather all data for the current track in a data container
    summary.aoi<-list.append(summary.aoi,
                             list(aoi.id=aoi.id,
                                  aoi.name=aoi.cur$layer,
                                  track.id=track.id,
                                  track.name=as.character(track.pts$trck_nm[1]),
                                  duration.total=duration.total,
                                  distance.total=dist.total,
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
                                  stops.time.summary=stops.time.summary,
                                  stops.count=stops.count,
                                  segments.nb=segments.nb,
                                  moments=moments,
                                  segments=segments
                                  )
                             )
    
    
    
  }
  
  # Convert summary list to dataframe
  summary<-data.frame(aoi.name=sapply(summary.aoi,"[[", "aoi.name"),
                      track.name=sapply(summary.aoi, "[[", "track.name"), 
                      duration=sapply(summary.aoi, "[[", "duration.total"), 
                      distance=sapply(summary.aoi, "[[", "distance.total"), 
                      speed.kmh.mean=sapply(summary.aoi, "[[", "speed.kmh.mean"), 
                      speed.lowess.mean=sapply(summary.aoi, "[[", "speed.lowess.mean"), 
                      duration.movement=sapply(summary.aoi, "[[", "duration.movement"), 
                      dist.movement=sapply(summary.aoi, "[[", "dist.movement"), 
                      speed.mean.movement=sapply(summary.aoi, "[[", "speed.mean.movement"),
                      speed.mean.movement.lowess=sapply(summary.aoi, "[[", "speed.mean.movement.lowess"),
                      duration.stops=sapply(summary.aoi, "[[", "duration.stops"), 
                      dist.stops=sapply(summary.aoi, "[[", "dist.stops"), 
                      speed.mean.stops=sapply(summary.aoi, "[[", "speed.mean.stops"),
                      speed.mean.stops.lowess=sapply(summary.aoi, "[[", "speed.mean.stops.lowess"),
                      stops.nb=sapply(summary.aoi, "[[", "stops.count"),
                      segments.nb=sapply(summary.aoi, "[[", "segments.nb"))
  
  # Convert factors to characters
  summary<-summary %>% mutate_if(is.factor, as.character)
  # Save summary as csv file
  write.csv2(x = summary, file = file.path(wd, dir.outputs, dir.tables, dir.meantime, paste0("summary-time-in-aoi-",summary.aoi[[1]]$aoi.id,".csv")))
  
  # Convert segments to dataframe
  segments<-do.call("rbind", lapply(summary.aoi, "[[", "segments"))
  segments<-segments[which(segments$segment.id != 0),]
  # Save segments as csv file
  write.csv2(x = segments, file = file.path(wd, dir.outputs, dir.tables, dir.meantime, paste0("segments-time-in-aoi-",summary.aoi[[1]]$aoi.id,".csv")))
  
  # Save summary as Rdata file
  saveRDS(object = summary.aoi, file = file.path(wd, dir.outputs, dir.rdata, paste0("summary-time-in-aoi-", aoi.id, ".rds")))
  
  # Gather all data for the current aoi in a data container
  summary.aois<-list.append(summary.aois, summary.aoi)

}

# Save outputs ------------------------------------------------------------
# Save summary as Rdata file
saveRDS(object = summary.aois, file = file.path(wd, dir.outputs, dir.rdata, "summary-time-in-aois.rds"))
# summary.aois<-readRDS(file = file.path(wd, dir.outputs, dir.rdata, "summary-time-in-aois.rds"))





