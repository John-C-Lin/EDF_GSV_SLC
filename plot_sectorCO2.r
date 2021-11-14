# Plots sectoral contributions in terms of CO2 from VULCAN's different sectors
# October 16th, 2021 by John C. Lin (John.Lin@utah.edu)

require(ggplot2); require(reshape2)
require(fields); require(ggmap)

#################
GSVdat <- readRDS("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/data/receptors/receptors_with_sectors.rds")
#YYYY.MM.DDsel <- "2019-07-10"
#YYYY.MM.DDsel <- "2019-08-08"
YYYY.MM.DDsel <- "2020-03-03"
#################

YYYY.MM.DD <- substring(as.character(GSVdat$time),1,10)
tapply(YYYY.MM.DD,YYYY.MM.DD,length) # list dates and how much data in each day
sel <- YYYY.MM.DD == YYYY.MM.DDsel
# dev.new()
# plot(GSVdat$longitude[sel],GSVdat$latitude[sel],pch=16,cex=0.5)
dat <- GSVdat[sel,]

register_google(key = "AIzaSyC8i2epZtRWGisCxJvOxKaimUf6s8GJctY")  #JCL's API key
#  entire SLC map
#map <- get_map(location = c(left=-112.1, bottom=40.54, right=-111.75, top=40.9), maptype = 'satellite', zoom = 10, color = 'bw')
#  dynamic map 
map <- get_map(location = c(left=min(dat$longitude,na.rm=T), bottom=min(dat$latitude,na.rm=T), 
                            right=max(dat$longitude,na.rm=T), top=max(dat$latitude,na.rm=T)), maptype = 'satellite', color = 'bw')


# (1) On-road
dev.new(); theme_set(theme_bw())
#g <- ggplot(dat,aes(x=longitude,y=latitude)) + geom_point(aes(col=onroad),size=1.0)
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=onroad),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

# (2) Residential
dev.new(); theme_set(theme_bw())
#g <- ggplot(dat,aes(x=longitude,y=latitude)) + geom_point(aes(col=residential),size=1.0)
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=residential),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

# (3) Industrial
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=industrial),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

# (4) Airport
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=airport),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

# (5) Commercial
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=commercial),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

# Calculate percentages from different sectors
percOnroad <- dat$onroad*100/dat$total
percIndustrial <- dat$industrial*100/dat$total
percResidential <- dat$residential*100/dat$total
dat <- data.frame(dat,percOnroad,percIndustrial,percResidential)

# (4) % Onroad
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=percOnroad),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

# (5) % Residential
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_point(data=dat,aes(x=longitude,y=latitude,col=percResidential),size=0.8)
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_color_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)


if(FALSE){
#generate gridded values 
dx<-0.002;dy<-0.002   #grid resolution [deg]
#dx<-0.005;dy<-0.005   #grid resolution [deg]
#NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
ymin<-40.49; ymax<-40.86    #grid limits [deg]
xmin<--112.19; xmax<--111.7
LATS<-seq(ymin,ymax,dy); LONS<-seq(xmin,xmax,dx)
#grid<-matrix(NA,nrow=length(LONS),ncol=length(LATS));dimnames(grid)<-list(LONS,LATS)
#  convert from lat/lon to grid coordinates
IIs<-floor((dat$longitude - xmin)/dx) + 1
JJs<-floor((dat$latitude  - ymin)/dy) + 1
# move LATS/LONS to represent CENTER of each gridcell, instead of SOUTHWEST (bottom-left) corner
LONS <- LONS + dx/2;  LATS <- LATS + dy/2

#(0) On-road CO2 
isNA<-is.na(dat$onroad)
N<-tapply(dat$onroad[!isNA],list(IIs[!isNA],JJs[!isNA]),length)   #sample size
dimnames(N)<-list(LONS[as.numeric(rownames(N))],LATS[as.numeric(colnames(N))])
med<-tapply(dat$onroad,list(IIs,JJs),median,na.rm=T)
dimnames(med)<-list(LONS[as.numeric(rownames(med))],LATS[as.numeric(colnames(med))])

# plot with image.plot
zlims<-NULL
#tmp<-med;tmp[tmp>zlims[2]]<-zlims[2];tmp[tmp<zlims[1]]<-zlims[1]
#X11();par(cex.lab=1.1,cex.axis=1.1,cex.main=1.3)
i#image.plot(as.numeric(rownames(med)),as.numeric(colnames(med)),tmp,xlab="Lon",ylab="Lat",zlim=zlims,main="Onroad CO2 from STILT-VULCAN[ppm]\n(median)")
#title(sub=paste(range(dat$time[!isNA]),collapse="   to   "))

# plot with geom_tile
mdat <- melt(med)
mdat <- mdat[!is.na(mdat$value),]
dev.new(); theme_set(theme_bw())
#g <- ggmap(map) + geom_tile(data=mdat,aes(x=Var1,y=Var2,fill=value))
g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=value)) + coord_cartesian()
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))
g <- g + scale_fill_gradientn(colours=tim.colors(7)) + labs(title=YYYY.MM.DDsel)
plot(g)

} # if(FALSE){

