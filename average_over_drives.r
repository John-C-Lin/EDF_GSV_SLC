# Average data over many drives/days
# October 16th, 2021 by John C. Lin (John.Lin@utah.edu)

require(ggplot2); require(reshape2)
require(fields); require(ggmap)
require(arrow)

#################
#GSVdat <- readRDS("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/data/receptors/receptors_with_sectors.rds")
#GSVdat <- readRDS("./receptors_with_sectors.rds")   # can't find this data on Ben's directory, so use backed up version from Google Drive
GSVdat <- read_feather("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/gsv-data/by_receptor/by_receptor.feather")
Nmin <- 10  # minimum number of receptors to retain 

# grid configuration
dx<-0.002;dy<-0.002   #grid resolution [deg]
#dx<-0.005;dy<-0.005    #grid resolution [deg]
#dx<-0.01;dy<-0.01    #grid resolution [deg]
#ymin<-40.49; ymax<-40.87    #grid limits [deg]
ymin<-40.57; ymax<-40.87    #grid limits [deg]
#xmin<--112.19; xmax<--111.7 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
xmin<--112.10; xmax<--111.76 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
LATS<-seq(ymin,ymax,dy); LONS<-seq(xmin,xmax,dx)
# move LATS/LONS to represent CENTER of each gridcell, instead of SOUTHWEST (bottom-left) corner
LONS <- LONS + dx/2;  LATS <- LATS + dy/2

# restrict date/times when want to focus analyses [GMT]
Time.start <- "2019-12-01 00:00:00"
Time.end <-   "2020-02-29 23:59:59"

outputdir <- "./out"   # where to store output
#################

dat <- GSVdat
# remove data outside of specified grid
sel <- dat$longitude < xmin | dat$longitude > xmax | dat$latitude < ymin | dat$latitude > ymax
dat <- dat[!sel,]

register_google(key = "AIzaSyC8i2epZtRWGisCxJvOxKaimUf6s8GJctY")  #JCL's API key
#  dynamic map 
#map <- get_map(location = c(left=min(dat$longitude,na.rm=T), bottom=min(dat$latitude,na.rm=T),
#                            right=max(dat$longitude,na.rm=T), top=max(dat$latitude,na.rm=T)), maptype = 'satellite', color = 'bw')
map <- get_map(location = c(left=xmin, bottom=ymin, right=xmax, top=ymax), maptype = 'satellite', color = 'bw')

# add NOx column
nox_ppb_ex <- dat$no_ppb_ex + dat$no2_ppb_ex
dat <- data.frame(dat,nox_ppb_ex)
dat$time <- as.POSIXct(dat$time)

# filter data to focus on time range
Time.start <- as.POSIXct(Time.start)
Time.end   <- as.POSIXct(Time.end)
sel <- (dat$time >= Time.start) & (dat$time <= Time.end)
dat <- dat[sel,]

#  convert from lat/lon to grid coordinates
IIs<-floor((dat$longitude - xmin)/dx) + 1
JJs<-floor((dat$latitude  - ymin)/dy) + 1

VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
          paste0(c("cmv","onroad","industrial","residential","commercial","elec_prod","nonroad","rail","airport","total"),"_ppm_ex"))

for(vv in 1:length(VARs)){

VAR <- VARs[vv]
print(paste("...... Processing",VAR," ......"))
isNA <- is.na(dat[,VAR])
N <- tapply(dat[!isNA,VAR],list(IIs[!isNA],JJs[!isNA]),length)   #sample size
dimnames(N) <- list(LONS[as.numeric(rownames(N))],LATS[as.numeric(colnames(N))])
med <- tapply(dat[,VAR],list(IIs,JJs),median,na.rm=T)
dimnames(med) <- list(LONS[as.numeric(rownames(med))],LATS[as.numeric(colnames(med))])
med[N < Nmin] <- NA   # not retain gridcells with too few datapoints
datfilenm <- paste0(VAR,".rds")
saveRDS(med,file=datfilenm);print(paste(datfilenm,"generated"))

print(quantile(N,na.rm=T))

ZLIMS <- NULL
if(VAR%in%c("onroad","industrial"))ZLIMS <- c(0,10)
if(VAR%in%c("residential","commercial","nonroad","rail"))ZLIMS <- c(0,2)
if(VAR=="pm25_ugm3_ex")ZLIMS <- c(0,10)
if(VAR=="co_ppb_ex")ZLIMS <- c(0,300)
if(VAR=="co2d_ppm_ex")ZLIMS <- c(0,60)
if(VAR=="ch4d_ppm_ex")ZLIMS <- c(0,0.1)
if(VAR=="total")ZLIMS <- c(0,60)
if(VAR%in%c("no_ppb_ex","no2_ppb_ex"))ZLIMS <- c(0,40)
if(VAR%in%c("nox_ppb_ex"))ZLIMS <- c(0,80)
if(VAR=="bc_ngm3_ex")ZLIMS <- c(0,2000)

# plot with image.plot
#zlims<-c(0,10)
#tmp<-med;tmp[tmp>zlims[2]]<-zlims[2];tmp[tmp<zlims[1]]<-zlims[1]
#X11();par(cex.lab=1.1,cex.axis=1.1,cex.main=1.3)
#image.plot(as.numeric(rownames(med)),as.numeric(colnames(med)),tmp,xlab="Lon",ylab="Lat",zlim=zlims,main="Onroad CO2 from STILT-VULCAN[ppm]\n(median)")
#title(sub=paste(range(dat$time[!isNA]),collapse="   to   "))

#!!!! test how good BC is in picking out major roadways !!!!
#if(VAR=="bc_ngm3_ex"){
#  sel <- med> 1000
#  med[sel] <- ZLIMS[2]
#  med[!sel] <- ZLIMS[1]
#} # if(VAR=="bc_ngm3_ex"){

# plot with geom_tile
mdat <- melt(med)
mdat <- mdat[!is.na(mdat$value),]
dev.new(); theme_set(theme_bw())
#g <- ggmap(map) + geom_tile(data=mdat,aes(x=Var1,y=Var2,fill=value))
g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=value)) + coord_cartesian()
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),plot.caption=element_text(size=12),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),plot.title=element_text(size=20))
g <- g + scale_fill_gradientn(colours=tim.colors(7), limits=ZLIMS)
g <- g + labs(title=VAR,caption=paste(range(dat$time[!isNA]),collapse="   to   "),tag=paste0("Nmin=",Nmin))
plot(g)
figfilenm <- paste0(VAR,".png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))
datfilenm <- paste0(VAR,".csv")
datout <- as.matrix(mdat)
dimnames(datout) <- list(NULL,c("longitude","latitude",VAR))
write.csv(datout,file=datfilenm);print(paste(datfilenm,"generated"))

gc()

} # for(vv in 1:length(VARs)){

file.copy(from=paste0(VARs,".png"),to=outputdir,overwrite=TRUE)
file.remove(paste0(VARs,".png"))
file.copy(from=paste0(VARs,".csv"),to=outputdir,overwrite=TRUE)
file.remove(paste0(VARs,".csv"))
file.copy(from=paste0(VARs,".rds"),to=outputdir,overwrite=TRUE)







