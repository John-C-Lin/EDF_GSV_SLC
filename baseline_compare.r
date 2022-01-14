# Compare calculated baseline against baseline sampled by GSV vehicles in locations outside of SLV core (e.g., Emigration Canyon)
# January 7th, 2022 by John C. Lin (John.Lin@utah.edu) 

require(ggplot2);require(reshape2)
require(fields); require(ggmap); require(arrow)
require(raster); require(tidyverse)

#################
GSVdat <- read_feather("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/gsv-data/by_receptor/by_receptor.feather")
#################

# add NOx column
nox_ppb <- GSVdat$no_ppb + GSVdat$no2_ppb
nox_ppb_base <- GSVdat$no_ppb_base + GSVdat$no2_ppb_base
nox_ppb_ex <- GSVdat$no_ppb_ex + GSVdat$no2_ppb_ex
GSVdat <- data.frame(GSVdat,nox_ppb,nox_ppb_base,nox_ppb_ex)
GSVdat$time <- as.POSIXct(GSVdat$time)
observations_bbox <- extent(-111.79, -111.67, 40.75, 40.79)   # Emigration Canyon
#observations <- read_feather("gsv-data/by_receptor/by_receptor.feather") %>%
baseline.all <- GSVdat %>%
    filter(
        longitude >= observations_bbox@xmin, longitude <= observations_bbox@xmax,
        latitude >= observations_bbox@ymin, latitude <= observations_bbox@ymax,)
emigrationTF <- GSVdat$longitude >= observations_bbox@xmin & GSVdat$longitude <= observations_bbox@xmax
emigrationTF <- emigrationTF & (GSVdat$latitude >= observations_bbox@ymin & GSVdat$latitude <= observations_bbox@ymax)


observations_bbox <- extent(-112.19, -112.148, 40.735, 40.774)   # small section along shore of Great Salt Lake
#observations <- read_feather("gsv-data/by_receptor/by_receptor.feather") %>%
baseline.all2 <- GSVdat %>%
    filter(
        longitude >= observations_bbox@xmin, longitude <= observations_bbox@xmax,
        latitude >= observations_bbox@ymin, latitude <= observations_bbox@ymax,)
# plot(observations$longitude,observations$latitude,pch=16,cex=0.5)
lakeshoreTF <- GSVdat$longitude >= observations_bbox@xmin & GSVdat$longitude <= observations_bbox@xmax
lakeshoreTF <- lakeshoreTF & (GSVdat$latitude >= observations_bbox@ymin & GSVdat$latitude <= observations_bbox@ymax)

GSVdat <- data.frame(GSVdat,emigrationTF,lakeshoreTF)

# check which days have lots of observations 
YYYYMMDD <- format(baseline.all$time,"%Y-%m-%d")
N.obs <- tapply(YYYYMMDD,YYYYMMDD,length)
N.obs <- N.obs[order(names(N.obs))]
print(N.obs)

# plot a particular day's time series
YYYYMMDD.sel <- "2019-10-24"
base.obs <- baseline.all[YYYYMMDD%in%YYYYMMDD.sel,]
base.obs2 <- baseline.all2[format(baseline.all2$time,"%Y=%m-%d")%in%YYYYMMDD.sel,]
GSVdat.sub <- GSVdat[format(GSVdat$time,"%Y-%m-%d")%in%YYYYMMDD.sel,]
tracers <- c("co2d_ppm", "nox_ppb", "bc_ngm3","pm25_ugm3")[-3]
for(i in 1:length(tracers)){
  tracer <- tracers[i]
  pngfile <- paste0("baseline_",tracer,"_",YYYYMMDD.sel,".png")
  png(filename=pngfile)
  plot(GSVdat.sub$time,GSVdat.sub[,paste0(tracer,"_base")],xlab="Time [UTC]",ylab=tracer,
       pch=16,cex=0.7,main=YYYYMMDD.sel)
  points(GSVdat.sub$time,GSVdat.sub[,tracer],pch=16,cex=0.7,col="lightgray")
  points(GSVdat.sub$time,GSVdat.sub[,paste0(tracer,"_base")],pch=16,cex=0.7)
  points(base.obs$time,base.obs[,paste0(tracer,"_base")],pch=16,cex=0.9,col="red")
  if(nrow(base.obs>1))points(base.obs2$time,base.obs2[,paste0(tracer,"_base")],pch=16,cex=0.9,col="blue")
  dev.off()
  print(paste(pngfile,"generated"))
} # for(i in 1:length(tracers)){



# plot MULTIPLE days' time series
# Time.start <- "2019-05-01 00:00:00"; Time.end <-   "2020-03-31 23:59:59"
# Time.start <- "2019-08-01 00:00:00"; Time.end <-   "2019-10-31 23:59:59"
# Time.start <- "2019-08-01 00:00:00"; Time.end <-   "2019-08-31 23:59:59"
Time.start <- "2019-08-15 00:00:00"; Time.end <-   "2019-08-31 23:59:59"
Time.start <- as.POSIXct(Time.start); Time.end   <- as.POSIXct(Time.end)
GSVdat.sub <- GSVdat[(GSVdat$time >= Time.start) & (GSVdat$time <= Time.end),]
tracers <- c("co2d_ppm", "nox_ppb", "bc_ngm3","pm25_ugm3")[-3]
for(i in 1:length(tracers)){
  tracer <- tracers[i]
  pngfile <- paste0("baseline_",tracer,"_",substring(as.character(Time.start),1,10),"_to_",
                    substring(as.character(Time.end),1,10),".png")
  png(filename=pngfile)
  SEL <- !GSVdat.sub$emigrationTF & !GSVdat.sub$lakeshoreTF
  plot(GSVdat.sub$time[SEL],GSVdat.sub[SEL,paste0(tracer,"_base")],xlab="Time [UTC]",ylab=tracer,
       pch=16,cex=0.7,main=YYYYMMDD.sel)
  points(GSVdat.sub$time[SEL],GSVdat.sub[SEL,tracer],pch=16,cex=0.7,col="lightgray")
  points(GSVdat.sub$time[SEL],GSVdat.sub[SEL,paste0(tracer,"_base")],pch=16,cex=0.7)
  points(GSVdat.sub$time[!SEL],GSVdat.sub[!SEL,paste0(tracer,"_base")],pch=16,cex=0.9,col="red")
  dev.off()
  print(paste(pngfile,"generated"))
} # for(i in 1:length(tracers)){





