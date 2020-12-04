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
dir.meandir<-"aoi"
dir.tables<-"tables"
dir.stops<-"stops"
dir.centroids<-"centroids"
dir.full<-"full"

# Delete outputs directory
if(delete.outputs == TRUE){
  
  unlink(file.path(wd, dir.outputs), recursive = TRUE)
  
}

# Create missing directories
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir))
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables))
}
if(!dir.exists(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.shp))){
  dir.create(file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.shp))
}
if(!dir.exists(file.path(wd, dir.outputs, dir.graphs, dir.aoi))){
  dir.create(file.path(wd, dir.outputs, dir.graphs, dir.aoi))
}
if(!dir.exists(file.path(wd, dir.outputs, dir.tables))){
  dir.create(file.path(wd, dir.outputs, dir.tables))
}

# Script AOI computation --------------------------------------------------

#list of gpx files in the working directory defined as input
# tracks.files<-list.files(file.path(wd, dir.inputs, dir.tracks), pattern=track.ext)

# Read AOI file for mean direction computation
aoi.spdf<-readOGR(dsn = file.path(wd, dir.inputs, dir.aoi, "aoi_for_mean_dir.shp"), stringsAsFactors = FALSE)
aoi.spdf<-spTransform(aoi.spdf, id.proj)

# Read shapefile containing points data for all tracks
track.spdf.full<-readOGR(dsn = file.path(wd, dir.outputs, dir.shp, dir.full, "gpx-full.shp"), stringsAsFactors = FALSE)

# Read Visitors classes
visitors<-read.csv(file = file.path(wd, dir.inputs, dir.tables, "visitors-categories.csv"), stringsAsFactors = FALSE)
colnames(visitors)<-c('date', 'class', 'id')
visitors<-visitors[,c(2,3)]

# Merge visitors category to spdf full
track.spdf.full@data<-merge(x = track.spdf.full@data, y = visitors, by.x = "track_d", by.y = "id")

# Movements points
track.spdf.mvt<-track.spdf.full[which(track.spdf.full@data$movemnt == 1),]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute data within aois for full data ----------------------------------

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
aoi.summary.prior<-data.frame()

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
  pts.cur.th<-track.spdf.mvt[aoi.cur,]
  
  # Compute total x and y components of movement
  sum.x<-sum(pts.cur@data$vec_x)
  sum.y<-sum(pts.cur@data$vec_y)
  sum.x.th<-sum(pts.cur.th@data$vec_x)
  sum.y.th<-sum(pts.cur.th@data$vec_y)
  
  # Compute mean total movement within aoi
  # mean.dir.tan<-atan(sum.y/sum.x)*(180/pi)
  mean.dir.tan2<-atan2(sum.y, sum.x)*(180/pi)
  # mean.dir.tan.th<-atan(sum.y.th/sum.x.th)*(180/pi)
  mean.dir.tan2.th<-atan2(sum.y.th, sum.x.th)*(180/pi)
  
  # Store in aoi spdf
  # aoi.spdf@data[aoi.id, "mean.tan"]<-mean.dir.tan
  aoi.spdf@data[aoi.id, "tan2"]<-mean.dir.tan2
  # aoi.spdf@data[aoi.id, "tan.th"]<-mean.dir.tan.th
  aoi.spdf@data[aoi.id, "tan2.th"]<-mean.dir.tan2.th
  
  # Store in list
  mean.dir.aois[[aoi.id]]<-list(mean.dir.tan2=mean.dir.tan2, 
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
  tracks.cur<-unique(pts.cur@data$track_d)
  id.line<-1
  # track.name<-tracks.cur[2]
  for(track.name in tracks.cur){
        
    # Selects points from current track 
    track.cur<-pts.cur[which(pts.cur@data$track_d == track.name),]
    track.cur<-track.cur[c(1, length(track.cur)),]
    
    # Get prior status for visitor
    prior.status<-track.cur@data$class[1]
    
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
    data <- data.frame(id = unique(d$id), id2=id.line, track=track.name, azimut=azimut.deg, prior=prior.status)
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
  
  # Compute azimut frequency per binned direction grouped by prior visit status
  aoi.lines.freq.prior<-aoi.lines@data %>% 
    group_by(prior) %>%
    count(azimut.class)
  colnames(aoi.lines.freq.prior)[3]<-"count"
  
  aoi.lines.freq.prior$freq.rel<-aoi.lines.freq.prior$count/sum(aoi.lines.freq.prior$count)
  aoi.lines.freq.prior$aoi.id<-aoi.id
  
  # Add azimut name to table
  aoi.lines.freq.prior<-merge(x=aoi.lines.freq.prior, y=azimut.names, by.x="azimut.class", by.y="class")
  
  aoi.summary.prior<-rbind(aoi.summary.prior, aoi.lines.freq.prior)
  
  # Save tables of frequencies and aoi.lines
  write.csv2(aoi.lines@data, file=file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables, paste0("aoi-", aoi.id, "-vectors.csv")))
  write.csv2(aoi.lines.freq, file=file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables, paste0("aoi-", aoi.id, "-azimut-freq.csv")))
  write.csv2(aoi.lines.freq.prior, file=file.path(wd, dir.outputs, dir.shp, dir.meandir, dir.tables, paste0("aoi-", aoi.id, "-azimut-freq-prior.csv")))
  
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
writeOGR(obj=aoi.spdf, dsn=file.path(wd, dir.outputs, dir.shp, dir.aoi, "aoi-mean-directions.shp"), layer="aoi-mean-directions", 
         driver="ESRI Shapefile", overwrite_layer=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Graphs - all visitors ---------------------------------------------------

# Arrange summary by azimut class for absolute number of people
aoi.summary.count <- aoi.summary %>% select(-freq.rel) %>% spread(aoi.id, count)
# Arrange summary by azimut class for relative number of people
aoi.summary.freq<- aoi.summary %>% select(-count) %>% spread(aoi.id, freq.rel)

# Histogram for directions within AOIS - ABSOLUTE VALUES
aoi.azimuts.count<-ggplot(data = aoi.summary, aes(x = azimut.name, y = count, fill = factor(aoi.id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  labs(title = "Analysis of visitors direction within AOI - Relative values",
       x = "Direction",
       y = "Nombre relatif de visiteurs",
       fill = "AOI")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#4B5F07", "#DA5F07"))
aoi.azimuts.count

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, dir.aoi, "aoi-directions-count.png"), plot = aoi.azimuts.count, device = "png")

# Histogram for directions within AOIS - RELATIVE VALUES
aoi.azimuts.freq<-ggplot(data = aoi.summary, aes(x = azimut.name, y = freq.rel, fill = factor(aoi.id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  labs(title = "Analysis of visitors direction within AOI - Relative values",
       x = "Direction",
       y = "Nombre  relatif de visiteurs",
       fill = "AOI")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#4B5F07", "#DA5F07"))
aoi.azimuts.freq

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, dir.aoi, "aoi-directions-freq.png"), plot = aoi.azimuts.freq, device = "png")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Graphs - Prior visits ---------------------------------------------------

# Arrange summary by azimut class for absolute number of people
aoi.summary.prior.count <- aoi.summary.prior %>% select(-freq.rel) %>% spread(aoi.id, count)
# Arrange summary by azimut class for relative number of people
aoi.summary.prior.freq<- aoi.summary.prior %>% select(-count) %>% spread(aoi.id, freq.rel)

# Histogram for directions within AOIS - ABSOLUTE VALUES
aoi.azimuts.prior.count.yes<-ggplot(data = aoi.summary.prior %>% filter(prior == "yes"), aes(x = azimut.name, y = count, fill = factor(aoi.id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  labs(title = "Analysis of visitors direction within AOI - Absolute values - PRIOR = YES",
       x = "Direction",
       y = "Nombre absolu de visiteurs (Prior = yes)",
       fill = "AOI")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#4B5F07", "#DA5F07"))
aoi.azimuts.prior.count.yes

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, dir.aoi, "aoi-directions-count-prior-yes.png"), plot = aoi.azimuts.count.yes, device = "png")

# Histogram for directions within AOIS - ABSOLUTE VALUES
aoi.azimuts.prior.count.no<-ggplot(data = aoi.summary.prior %>% filter(prior == "no"), aes(x = azimut.name, y = count, fill = factor(aoi.id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  labs(title = "Analysis of visitors direction within AOI - Absolute values - PRIOR = NO",
       x = "Direction",
       y = "Nombre absolu de visiteurs (Prior = no)",
       fill = "AOI")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#4B5F07", "#DA5F07"))
aoi.azimuts.prior.count.no

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, dir.aoi, "aoi-directions-count-prior-no.png"), plot = aoi.azimuts.count.no, device = "png")

# Histogram for directions within AOIS - RELATIVE VALUES
aoi.azimuts.prior.freq.yes<-ggplot(data = aoi.summary.prior %>% filter(prior == "yes"), aes(x = azimut.name, y = freq.rel, fill = factor(aoi.id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  labs(title = "Analysis of visitors direction within AOI - Relative values - PRIOR = YES",
       x = "Direction",
       y = "Nombre relatif de visiteurs (Prior = yes)",
       fill = "AOI")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#4B5F07", "#DA5F07"))
aoi.azimuts.prior.freq.yes

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, dir.aoi, "aoi-directions-freq-prior-yes.png"), plot = aoi.azimuts.freq.yes, device = "png")

# Histogram for directions within AOIS - RELATIVE VALUES
aoi.azimuts.prior.freq.no<-ggplot(data = aoi.summary.prior %>% filter(prior == "no"), aes(x = azimut.name, y = freq.rel, fill = factor(aoi.id)))+
  geom_bar(stat="identity", position = "dodge")+
  theme_bw()+
  labs(title = "Analysis of visitors direction within AOI - Relative values - PRIOR = NO",
       x = "Direction",
       y = "Nombre relatif de visiteurs (Prior = no)",
       fill = "AOI")+
  scale_fill_manual(values=c("#00AFBB", "#E7B800", "#4B5F07", "#DA5F07"))
aoi.azimuts.prior.freq.no

ggsave(filename = file.path(wd, dir.outputs, dir.graphs, dir.aoi, "aoi-directions-freq-prior-no.png"), plot = aoi.azimuts.freq.no, device = "png")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save tables -------------------------------------------------------------

# Save aoi summaries to csv
write.csv2(aoi.summary, file=file.path(wd, dir.outputs, dir.tables, "aoi-directions-summary.csv"))
write.csv2(aoi.summary.prior, file=file.path(wd, dir.outputs, dir.tables, "aoi-directions-summary-prior.csv"))

write.csv2(aoi.summary.count, file=file.path(wd, dir.outputs, dir.tables, "aoi-directions-summary-count.csv"))
write.csv2(aoi.summary.freq, file=file.path(wd, dir.outputs, dir.tables, "aoi-directions-summary-freq.csv"))

write.csv2(aoi.summary.prior.count, file=file.path(wd, dir.outputs, dir.tables, "aoi-directions-summary-prior-count.csv"))
write.csv2(aoi.summary.prior.freq, file=file.path(wd, dir.outputs, dir.tables, "aoi-directions-summary-prior-freq.csv"))


