# Calculate tracer:tracer RATIOS over different gridcells and correlates with STILT-VULCAN sectors [ppm of CO2]
# V2(211023):  add multi-panel plot correlating tracer ratio with STILT-VULCAN sector tracer using 'facet_wrap'
# Reads in output from "average_drive_drives.r"
# October 22nd, 2021 by John C. Lin (John.Lin@utah.edu)

require(ggplot2); require(reshape2)
require(fields); require(ggmap)

#####################
# calculate ratios between VARs.2/VARs.1
#VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
#          "onroad","industrial","residential","commercial","total")
VARs.2 <- c("nox_ppb_ex","bc_ngm3_ex","pm25_ugm3_ex","ch4d_ppm_ex")[1]
VARs.1 <- c("co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex")[1]

# grid configuration
ymin<-40.57; ymax<-40.87    #grid limits [deg]
xmin<--112.10; xmax<--111.76 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell
#################

register_google(key = "AIzaSyC8i2epZtRWGisCxJvOxKaimUf6s8GJctY")  #JCL's API key
map <- get_map(location = c(left=xmin, bottom=ymin, right=xmax, top=ymax), maptype = 'satellite', color = 'bw')


if(length(VARs.1)!=length(VARs.2))stop("length of VARs.1 and VARs.2 need to be the same")

regress.SMA <- function(x,y){
  require(lmodel2)
  xfit<-lmodel2(formula=y~x)
  sel<-xfit$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  result <- xfit$regression.results[sel,]
  return(result)
} #f<-(x,y)


##############################################################################################################
#       Tracer ratio as:   median(tracer2)/median(tracer1)
#------------------#
#  Create the objects storing tracer ratio as median(tracer2)/median(tracer1)
xsub <- "All Data"
for(i in 1:length(VARs.1)){
  VAR1 <- VARs.1[i]; VAR2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(VAR1,".rds"))
  vdat2 <- readRDS(paste0(VAR2,".rds"))

  ratio <- vdat2/vdat1
  LAB <- paste0("median(",VAR2,"):median(",VAR1,")")
  datfilenm <- paste0(LAB,".rds")
  saveRDS(ratio,file=datfilenm);print(paste(datfilenm,"generated"))

  #ZLIMS <- NULL
  ZLIMS <- quantile(ratio,na.rm=T,probs=c(0,0.95))
if(VAR1=="co2d_ppm_ex"){
  if(VAR2%in%c("nox_ppb_ex"))ZLIMS <- c(0,3)
  if(VAR2=="bc_ngm3_ex")ZLIMS <- c(0,500)
  if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(0,1.0)
  #if(VAR2=="co_ppb_ex")ZLIMS <- c(0,300)
  if(VAR2=="ch4d_ppm_ex")ZLIMS <- c(0,0.004)
} # if(VAR1=="co2d_ppm_ex"){

  # plot with geom_tile
  mdat <- melt(ratio)
  mdat <- mdat[!is.na(mdat$value),]
  dev.new(); theme_set(theme_bw())
  g <- ggmap(map) + geom_raster(data=mdat,aes(x=Var1,y=Var2,fill=value)) + coord_cartesian()
  g <- g + theme(axis.title.y=element_text(size=14), axis.title.x=element_text(size=14),plot.caption=element_text(size=12),
               axis.text.x=element_text(size=14), axis.text.y=element_text(size=14),plot.title=element_text(size=20))
  g <- g + scale_fill_gradientn(colours=tim.colors(7), limits=ZLIMS)
  g <- g + labs(title=LAB)
  plot(g)

  figfilenm <- paste0(LAB,".png")
  ggsave(figfilenm);print(paste(figfilenm,"generated"))
  datfilenm <- paste0(LAB,".csv")
  datout <- as.matrix(mdat)
  dimnames(datout) <- list(NULL,c("longitude","latitude",LAB))
  write.csv(datout,file=datfilenm);print(paste(datfilenm,"generated"))

  gc()
} # for(i in 1:length(VARs.1)){

# create data object that can be used in 'facet_grid'
VULCANvars <- c("onroad","industrial","residential","commercial") # contributions from STILT-VULCAN sectors [ppm] that want to correlate with tracer ratio
for(i in 1:length(VARs.1)){
  VAR1 <- VARs.1[i]; VAR2 <- VARs.2[i]
  LAB <- paste0("median(",VAR2,"):median(",VAR1,")")
  ratio <- readRDS(file=paste0(LAB,".rds"))
  DX <- unique(signif(diff(as.numeric(rownames(ratio))),4))[1]
  DY <- unique(signif(diff(as.numeric(colnames(ratio))),4))[1]

  ZLIMS <- NULL
  if(VAR1=="co2d_ppm_ex"){
    if(VAR2%in%c("nox_ppb_ex"))ZLIMS <- c(0,3)
    if(VAR2=="bc_ngm3_ex")ZLIMS <- c(0,500)
    if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(0,1.0)
    #if(VAR2=="co_ppb_ex")ZLIMS <- c(0,300)
    if(VAR2=="ch4d_ppm_ex")ZLIMS <- c(0,0.004)
  } # if(VAR1=="co2d_ppm_ex"){

  ###############
  #(1) ALL DATA
  DAT <- NULL
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(ratio)[1] | dim(vulcan)[2]!=dim(ratio)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer ratio grid")
  tmp <- data.frame(as.vector(ratio),VULCANvar,as.vector(vulcan))
  DAT <- rbind(DAT,tmp)
} # for(j in 1:length(VULCANvars))
  colnames(DAT) <- c("ratio","VULCANsector","CO2_ppm")
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors

  isNA <- is.na(dat$CO2_ppm)|is.na(dat[,1])
  dat <- dat[!isNA,]
  xsub <- "All Data"
  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=ratio)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(c(0,20))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
   ggsave("tracer_ratio_1.png")


  ###############
  #  (2) Use BC to select for MAJOR ROADWAYS (e.g., highways)
  BCthresh <- 1000
  BCdat <- readRDS("bc_ngm3_ex.rds"); SEL.BC <- BCdat > BCthresh
  SEL <- SEL.BC
  ratio <- readRDS(file=paste0(LAB,".rds"))
  if(length(ratio)!=length(BCdat))stop("dimensions of BC data and ratio data are NOT the same")
  xsub <- paste0("bc_ngm3_ex > ",BCthresh)
  DAT <- NULL
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(ratio)[1] | dim(vulcan)[2]!=dim(ratio)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer ratio grid")
  tmp <- data.frame(as.vector(ratio[SEL]),VULCANvar,as.vector(vulcan[SEL]))
  DAT <- rbind(DAT,tmp)
} # for(j in 1:length(VULCANvars))
  colnames(DAT) <- c("ratio","VULCANsector","CO2_ppm")
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors
  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=ratio)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(c(0,20))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
  ggsave("tracer_ratio_2.png")


  ###############
  #  (3) Use BC to filter OUT MAJOR ROADWAYS (e.g., highways)
  BCthresh <- 1000; BCdat <- readRDS("bc_ngm3_ex.rds")
  SEL.BC <- BCdat > BCthresh
  SEL <- !SEL.BC
  ratio <- readRDS(file=paste0(LAB,".rds"))
  if(length(ratio)!=length(BCdat))stop("dimensions of BC data and ratio data are NOT the same")
  xsub <- paste0("bc_ngm3_ex <= ",BCthresh)
  DAT <- NULL
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(ratio)[1] | dim(vulcan)[2]!=dim(ratio)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer ratio grid")
  tmp <- data.frame(as.vector(ratio[SEL]),VULCANvar,as.vector(vulcan[SEL]))
  DAT <- rbind(DAT,tmp)
} # for(j in 1:length(VULCANvars))
  colnames(DAT) <- c("ratio","VULCANsector","CO2_ppm")
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors
  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=ratio)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(c(0,20))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
  ggsave("tracer_ratio_3.png")
} # for(i in 1:length(VARs.1)){



##############################################################################################################
#       Tracer ratio as:   median(tracer2/tracer1)
#------------------#
for(i in 1:length(VARs.1)){
  VAR1 <- VARs.1[i]; VAR2 <- VARs.2[i]
  LAB <- paste0("median(",VAR2,":",VAR1,")")
  ratio <- readRDS(paste0(LAB,".rds"))
  DX <- unique(signif(diff(as.numeric(rownames(ratio))),4))[1]
  DY <- unique(signif(diff(as.numeric(colnames(ratio))),4))[1]

  # create data object that can be used in 'facet_grid'
VULCANvars <- c("onroad","industrial","residential","commercial") # contributions from STILT-VULCAN sectors [ppm] that want to correlate with tracer ratio

  ZLIMS <- NULL
  if(VAR1=="co2d_ppm_ex"){
    if(VAR2%in%c("nox_ppb_ex"))ZLIMS <- c(0,3)
    if(VAR2=="bc_ngm3_ex")ZLIMS <- c(0,500)
    if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(0,1.0)
    #if(VAR2=="co_ppb_ex")ZLIMS <- c(0,300)
    if(VAR2=="ch4d_ppm_ex")ZLIMS <- c(0,0.004)
  } # if(VAR1=="co2d_ppm_ex"){

  ###############
  #(1) ALL DATA
  xsub <- "All Data"
  DAT <- NULL
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(ratio)[1] | dim(vulcan)[2]!=dim(ratio)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer ratio grid")
  tmp <- data.frame(as.vector(ratio),VULCANvar,as.vector(vulcan))
  DAT <- rbind(DAT,tmp)
} # for(j in 1:length(VULCANvars))
  colnames(DAT) <- c("ratio","VULCANsector","CO2_ppm")
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors

  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=ratio)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(c(0,20))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
  ggsave("tracer_ratio_4.png")


  ###############
  #  (3) Use BC to filter OUT MAJOR ROADWAYS (e.g., highways)
  BCthresh <- 1000; BCdat <- readRDS("bc_ngm3_ex.rds")
  SEL.BC <- BCdat > BCthresh
  SEL <- !SEL.BC
  ratio <- readRDS(file=paste0(LAB,".rds"))
  if(length(ratio)!=length(BCdat))stop("dimensions of BC data and ratio data are NOT the same")
  xsub <- paste0("bc_ngm3_ex <= ",BCthresh)
  DAT <- NULL
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(ratio)[1] | dim(vulcan)[2]!=dim(ratio)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer ratio grid")
  tmp <- data.frame(as.vector(ratio[SEL]),VULCANvar,as.vector(vulcan[SEL]))
  DAT <- rbind(DAT,tmp)
} # for(j in 1:length(VULCANvars))
  colnames(DAT) <- c("ratio","VULCANsector","CO2_ppm")
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors
  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=ratio)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(c(0,20))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
  ggsave("tracer_ratio_6.png")

} # for(i in 1:length(VARs.1)){



