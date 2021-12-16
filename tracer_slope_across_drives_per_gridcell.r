# Calculate tracer slopes, correlation, and p-values across different drives within each gridcell
# October 24th, 2021 by John C. Lin (John.Lin@utah.edu)
# NOTE:  based on "average_tracer_ratios_over_drivesV1.r"

require(ggplot2); require(reshape2)
require(fields); require(ggmap)
require(arrow)

#################
# calculate ratios between VARs.2/VARs.1
#VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
#          paste0(c("cmv","onroad","industrial","residential","commercial","elec_prod","nonroad","rail","airport","total"),"_ppm_ex"))
VARs.2 <- c("nox_ppb_ex","bc_ngm3_ex","pm25_ugm3_ex","ch4d_ppm_ex")
VARs.1 <- c("co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex")

#GSVdat <- readRDS("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/data/receptors/receptors_with_sectors.rds")
#GSVdat <- readRDS("./receptors_with_sectors.rds")   # can't find this data on Ben's directory, so use backed up version from Google Drive
GSVdat <- read_feather("/uufs/chpc.utah.edu/common/home/lin-group8/btf/google-street-view/gsv-data/by_receptor/by_receptor.feather")
Nmin <- 20  # minimum number of receptors to retain 

# grid configuration
#dx<-0.002;dy<-0.002   #grid resolution [deg]
#dx<-0.005;dy<-0.005   #grid resolution [deg]
dx<-0.01;dy<-0.01   #grid resolution [deg]
#ymin<-40.49; ymax<-40.87    #grid limits [deg]
ymin<-40.57; ymax<-40.87    #grid limits [deg]
#xmin<--112.19; xmax<--111.7 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
xmin<--112.10; xmax<--111.76 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
LATS<-seq(ymin,ymax,dy); LONS<-seq(xmin,xmax,dx)
# move LATS/LONS to represent CENTER of each gridcell, instead of SOUTHWEST (bottom-left) corner
LONS <- LONS + dx/2;  LATS <- LATS + dy/2

p.value.threshold <- 0.05   # threshold for p-value for correlation strength used for filtering gridcells

# restrict date/times when want to focus analyses [GMT]
#Time.start <- "2019-12-01 00:00:00"
#Time.end <-   "2020-02-29 23:59:59"
Time.start <- "2019-05-01 00:00:00"
Time.end <-   "2020-03-31 23:59:59"

regressTF <- TRUE  # whether or not to carry out regression (time-consuming step); if set to FALSE, then will generate map of the tracer slopes
outputdir <- "./out"   # where to store output
#################

if(length(VARs.1)!=length(VARs.2))stop("length of VARs.1 and VARs.2 need to be the same")

regress.SMA.slope <- function(x,y){
  require(lmodel2)
  #print(paste("in regress.SMA.slope:",length(x),length(y)))
  try(xfit<-lmodel2(formula=y~x),silent=TRUE)
  if(exists("xfit")){
    sel<-xfit$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
    result <- xfit$regression.results[sel,"Slope"]
  } else {result <- NA} # if(exists(xfit)){
  return(result)
} # regress.SMA.slope <- function(x,y){

dat <- GSVdat
# remove data outside of specified grid
sel <- dat$longitude < xmin | dat$longitude > xmax | dat$latitude < ymin | dat$latitude > ymax
dat <- dat[!sel,]

register_google(key = "AIzaSyC8i2epZtRWGisCxJvOxKaimUf6s8GJctY")  #JCL's API key
map <- get_map(location = c(left=xmin, bottom=ymin, right=xmax, top=ymax), maptype = 'satellite', color = 'bw')

nox_ppb_ex <- dat$no_ppb_ex + dat$no2_ppb_ex
dat <- data.frame(dat,nox_ppb_ex)

# filter data to focus on time range
Time.start <- as.POSIXct(Time.start)
Time.end   <- as.POSIXct(Time.end)
sel <- (dat$time >= Time.start) & (dat$time <= Time.end)
dat <- dat[sel,]

#  convert from lat/lon to grid coordinates
IIs<-floor((dat$longitude - xmin)/dx) + 1
JJs<-floor((dat$latitude  - ymin)/dy) + 1


for(vv in 1:length(VARs.1)){

VAR1 <- VARs.1[vv]; VAR2 <- VARs.2[vv]
print(paste0("...... Processing ",VAR2,"/",VAR1," ......"))
isNA <- is.na(dat[,VAR1])|is.na(dat[,VAR2])
N <- tapply(dat[!isNA,VAR1],list(IIs[!isNA],JJs[!isNA]),length)   #sample size
dimnames(N) <- list(LONS[as.numeric(rownames(N))],LATS[as.numeric(colnames(N))])

LAB <- paste0(VAR2,"_SLOPE_",VAR1)
datfilenm <- paste0(LAB,".rds")

if(!regressTF){print("...............Skipping regression step!!...............")}
if(regressTF){
  slope <- N; slope[1:length(slope)] <- NA
  corr <- slope; p.value <- slope
  INDs <- unique(paste0(IIs[!isNA],",",JJs[!isNA]))
  for(k in 1:length(INDs)){
    print(paste(k,"out of",length(INDs)))
    i <- as.numeric(strsplit(INDs[k],split=",")[[1]][1])
    j <- as.numeric(strsplit(INDs[k],split=",")[[1]][2])
    if(i>nrow(slope) | j>ncol(slope))next
    sel <- IIs[!isNA]==i&JJs[!isNA]==j
    # print(paste(i,j,sum(sel)))
    if(sum(sel)<Nmin)next
    slope[i,j] <- regress.SMA.slope(x=dat[!isNA,VAR1][sel],y=dat[!isNA,VAR2][sel])      # regression slope
    # corr[i,j] <- cor(dat[!isNA,VAR1][sel],dat[!isNA,VAR2][sel], use="na.or.complete")   # correlation coefficient
    xcor <- cor.test(dat[!isNA,VAR1][sel],dat[!isNA,VAR2][sel], method="pearson")       # test for correlation, and calculate correlation coefficient
    corr[i,j] <- xcor$estimate
    p.value[i,j] <- xcor$p.value
  } # for(k in 1:length(INDs)){

  result <- list(slope=slope,corr=corr,N=N,p.value=p.value)
  saveRDS(result,file=datfilenm);print(paste(datfilenm,"generated"))
} # if(regressTF){

result <- readRDS(datfilenm)

ZLIMS <- NULL
if(VAR1=="co2d_ppm_ex"){
  if(VAR2%in%c("nox_ppb_ex"))ZLIMS <- c(0,3)
  if(VAR2=="bc_ngm3_ex")ZLIMS <- c(0,500)
  #if(VAR2=="bc_ngm3_ex")ZLIMS <- c(-1000,0)
  if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(0,1.0)
  #if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(-1.0,0)
  #if(VAR2=="co_ppb_ex")ZLIMS <- c(0,300)
  if(VAR2=="ch4d_ppm_ex")ZLIMS <- c(0,0.004)
} # if(VAR1=="co2d_ppm_ex"){

# 1) Plot slope with geom_tile
mdat <- melt(result$slope); mdat <- mdat[!is.na(result$slope),]
colnames(mdat)[colnames(mdat)=="value"] <- "Slope"
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=Slope)) + coord_cartesian()
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),plot.caption=element_text(size=12),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),plot.title=element_text(size=20))
g <- g + scale_fill_gradientn(colours=tim.colors(7), limits=ZLIMS)
g <- g + labs(title=LAB,tag=paste0("Nmin=",Nmin))
g <- g + labs(caption=paste(paste(range(dat$time[!isNA]),collapse=" to "),paste0(";  dx=",dx,"; dy=",dy)))
plot(g)
figfilenm <- paste0(LAB,".png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))
datfilenm <- paste0(LAB,".csv")
datout <- as.matrix(mdat)
dimnames(datout) <- list(NULL,c("longitude","latitude",LAB))
write.csv(datout,file=datfilenm);print(paste(datfilenm,"generated"))

# 1.5) Plot slope with geom_tile, after filtering for p.value
mdat <- melt(result$slope)
sel <- result$p.value < p.value.threshold & !is.na(result$slope)
mdat <- mdat[sel,]
colnames(mdat)[colnames(mdat)=="value"] <- "Slope"
xsub <- paste0("All Data; p<",p.value.threshold)
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=Slope)) + coord_cartesian()
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),plot.caption=element_text(size=12),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),plot.title=element_text(size=20))
g <- g + scale_fill_gradientn(colours=tim.colors(7), limits=ZLIMS)
g <- g + labs(title=LAB,tag=paste0("Nmin=",Nmin))
g <- g + labs(subtitle=xsub,caption=paste(paste(range(dat$time[!isNA]),collapse=" to "),paste0(";  dx=",dx,"; dy=",dy)))
plot(g)
figfilenm <- paste0(LAB,"p-value_filtered.png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))


# 2) Plot correlation with geom_tile
#ZLIMS.R <- NULL
ZLIMS.R <- c(-1.0,1.0)
mdat <- melt(result$corr); mdat <- mdat[!is.na(result$corr),]
colnames(mdat)[colnames(mdat)=="value"] <- "Correlation"
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=Correlation)) + coord_cartesian()
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),plot.caption=element_text(size=12),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),plot.title=element_text(size=20))
g <- g + scale_fill_gradientn(colours=two.colors(7,start="darkblue"), limits=ZLIMS.R)
g <- g + labs(title=LAB,tag=paste0("Nmin=",Nmin))
g <- g + labs(caption=paste(paste(range(dat$time[!isNA]),collapse=" to "),paste0(";  dx=",dx,"; dy=",dy)))
plot(g)
figfilenm <- paste0(LAB,"_R.png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))

# 3) Plot p.value with geom_tile
ZLIMS.R <- NULL
#ZLIMS.R <- c(-1.0,1.0)
mdat <- melt(result$p.value); mdat <- mdat[!is.na(result$p.value),]
colnames(mdat)[colnames(mdat)=="value"] <- "p.value"
dev.new(); theme_set(theme_bw())
g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=p.value)) + coord_cartesian()
g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),plot.caption=element_text(size=12),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),plot.title=element_text(size=20))
g <- g + scale_fill_gradientn(colours=two.colors(7,start="darkblue"), limits=ZLIMS.R)
g <- g + labs(title=LAB,tag=paste0("Nmin=",Nmin))
g <- g + labs(caption=paste(paste(range(dat$time[!isNA]),collapse=" to "),paste0(";  dx=",dx,"; dy=",dy)))
plot(g)
figfilenm <- paste0(LAB,"_p-value.png")
ggsave(figfilenm);print(paste(figfilenm,"generated"))

# 4) Plot p-value versus correlation coefficient 
dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
plot(result$corr,result$p.value,pch=16)
abline(h=p.value.threshold,col="red")
xsub <- paste(sum(result$p.value < p.value.threshold,na.rm=T),"out of total of",sum(!is.na(result$p.value)),"under p.value=",p.value.threshold)
title(main=LAB,sub=xsub)
figfilenm <- paste0(LAB,"_p-value_vs_R.png")
dev.copy(png,figfilenm);print(paste(figfilenm,"generated"));dev.off()

gc()

} # for(vv in 1:length(VARs)){
graphics.off()

xfiles <- list.files(pattern="_SLOPE_")
file.copy(from=xfiles,to=outputdir,overwrite=TRUE)




