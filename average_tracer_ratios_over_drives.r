# Average tracer ratio data over many drives/days
# October 22nd, 2021 by John C. Lin (John.Lin@utah.edu)
# NOTE:  based on "average_over_drives.r"

require(ggplot2); require(reshape2)
require(fields); require(ggmap)

#################
# calculate ratios between VARs.2/VARs.1
#VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
#          "onroad","industrial","residential","commercial","total")
VARs.2 <- c("nox_ppb_ex","bc_ngm3_ex","pm25_ugm3_ex","ch4d_ppm_ex")
VARs.1 <- c("co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex")


GSVdat <- readRDS("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/data/receptors/receptors_with_sectors.rds")
Nmin <- 20  # minimum number of receptors to retain 

# grid configuration
dx<-0.002;dy<-0.002   #grid resolution [deg]
#ymin<-40.49; ymax<-40.87    #grid limits [deg]
ymin<-40.57; ymax<-40.87    #grid limits [deg]
#xmin<--112.19; xmax<--111.7 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
xmin<--112.10; xmax<--111.76 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
LATS<-seq(ymin,ymax,dy); LONS<-seq(xmin,xmax,dx)
# move LATS/LONS to represent CENTER of each gridcell, instead of SOUTHWEST (bottom-left) corner
LONS <- LONS + dx/2;  LATS <- LATS + dy/2
#################

if(length(VARs.1)!=length(VARs.2))stop("length of VARs.1 and VARs.2 need to be the same")


dat <- GSVdat
# remove data outside of specified grid
sel <- dat$longitude < xmin | dat$longitude > xmax | dat$latitude < ymin | dat$latitude > ymax
dat <- dat[!sel,]

register_google(key = "AIzaSyC8i2epZtRWGisCxJvOxKaimUf6s8GJctY")  #JCL's API key
map <- get_map(location = c(left=xmin, bottom=ymin, right=xmax, top=ymax), maptype = 'satellite', color = 'bw')

#  convert from lat/lon to grid coordinates
IIs<-floor((dat$longitude - xmin)/dx) + 1
JJs<-floor((dat$latitude  - ymin)/dy) + 1

nox_ppb_ex <- dat$no_ppb_ex + dat$no2_ppb_ex
dat <- data.frame(dat,nox_ppb_ex)

VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
          "onroad","industrial","residential","commercial","total")[6]

for(vv in 1:length(VARs.1)){

VAR1 <- VARs.1[vv]; VAR2 <- VARs.2[vv]
print(paste0("...... Processing ",VAR2,"/",VAR1," ......"))
isNA<-is.na(dat[,VAR1])|is.na(dat[,VAR2])
N<-tapply(dat[!isNA,VAR1],list(IIs[!isNA],JJs[!isNA]),length)   #sample size
dimnames(N)<-list(LONS[as.numeric(rownames(N))],LATS[as.numeric(colnames(N))])
med<-tapply(dat[,VAR2]/dat[,VAR1],list(IIs,JJs),median,na.rm=T)
dimnames(med)<-list(LONS[as.numeric(rownames(med))],LATS[as.numeric(colnames(med))])
med[N < Nmin] <- NA   # not retain gridcells with too few datapoints
LAB <- paste0("median(",VAR2,":",VAR1,")")

datfilenm <- paste0(LAB,".rds")
saveRDS(med,file=datfilenm);print(paste(datfilenm,"generated"))

print(quantile(N,na.rm=T))

ZLIMS <- NULL
if(VAR1=="co2d_ppm_ex"){
  if(VAR2%in%c("nox_ppb_ex"))ZLIMS <- c(0,3)
  if(VAR2=="bc_ngm3_ex")ZLIMS <- c(0,500)
  if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(0,1.0)
  #if(VAR2=="co_ppb_ex")ZLIMS <- c(0,300)
  if(VAR2=="ch4d_ppm_ex")ZLIMS <- c(0,0.004)
} # if(VAR1=="co2d_ppm_ex"){

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
g <- g + labs(title=LAB,caption=paste(range(dat$time[!isNA]),collapse="   to   "),tag=paste0("Nmin=",Nmin))
plot(g)
figfilenm <- paste0(LAB,".png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))
datfilenm <- paste0(LAB,".csv")
datout <- as.matrix(mdat)
dimnames(datout) <- list(NULL,c("longitude","latitude",LAB))
write.csv(datout,file=datfilenm);print(paste(datfilenm,"generated"))

gc()

} # for(vv in 1:length(VARs)){






