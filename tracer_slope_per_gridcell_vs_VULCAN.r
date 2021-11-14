# Calculate tracer:tracer SLOPES over different gridcells and correlates with STILT-VULCAN sectors [ppm of CO2]
# Reads in output from "tracer_slope_across_drives_per_gridcell.r"
# V2(211102): Add statistics (R, slope, N) to multi-panel ggplot (facet_wrap)
# October 25th, 2021 by John C. Lin (John.Lin@utah.edu)

require(ggplot2); require(reshape2)
require(fields); require(ggmap)
require(lmodel2)

#####################
# calculate ratios between VARs.2/VARs.1
#VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
#          "onroad","industrial","residential","commercial","total")
VARs.2 <- c("nox_ppb_ex","bc_ngm3_ex","pm25_ugm3_ex","ch4d_ppm_ex")[1:3]
VARs.1 <- c("co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex")[1:3]

# contributions from STILT-VULCAN sectors [ppm] that want to correlate with tracer slope
VULCANvars <- paste0(c("onroad","nonroad","industrial","rail","residential","commercial"),"_ppm_ex")

# grid configuration
ymin<-40.57; ymax<-40.87    #grid limits [deg]
xmin<--112.10; xmax<--111.76 #NOTE:  currently using BOTTOM-LEFT (southwest) corner's lat/lon coordinates to label gridcell

# threshold for p-value (calculated in "tracer_slope_across_drives_per_gridcell.r"); when p>threshold, then THROW AWAY the calculated regression slope
p.value.threshold <- 0.05

XLIMS <- c(0,20)   # range of VULCAN sector contributions [ppm]
outputdir <- "./out"   # where to store output
#################

register_google(key = "AIzaSyC8i2epZtRWGisCxJvOxKaimUf6s8GJctY")  #JCL's API key
map <- get_map(location = c(left=xmin, bottom=ymin, right=xmax, top=ymax), maptype = 'satellite', color = 'bw')


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

regress.SMA <- function(x,y){
  require(lmodel2)
  xfit<-lmodel2(formula=y~x)
  sel<-xfit$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  result <- xfit$regression.results[sel,]
  return(result)
} #f<-(x,y)


bootlmodel2 <- function(x,y,N=1000,method="SMA") {
  if(length(x)!=length(y))stop(paste("x and y need to be the same length!"))
  # bootstrap to calculate errors in regression slope and intercept, based on Model-II regression
  result <- matrix(NA,nrow=N,ncol=2); colnames(result) <- c("m","interc")
  result <- data.frame(result)
  for(i in 1:N){ 
    ind.s <- sample(1:length(x),replace=TRUE)
    x2 <- x[ind.s]
    y2 <- y[ind.s]
    #xlm <- lmodel2(formula=y ~ x)
    xlm <- lmodel2(formula=y2 ~ x2)
    sel <- xlm$regression.results[,"Method"]==method #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
    result[i,"m"] <- xlm$regression.results[sel,"Slope"]
    result[i,"interc"] <- xlm$regression.results[sel,"Intercept"]
  } # for(i in 1:N){ 
  return(result)
} # bootlmodel2 <- (x,y,N=1000) {

##############################################################################################################
#       Read in tracer2:tracer1 regression SLOPES
#------------------#
for(i in 1:length(VARs.1)){
  VAR1 <- VARs.1[i]; VAR2 <- VARs.2[i]
  LAB <- paste0(VAR2,"_SLOPE_",VAR1)
  print(paste("--------- Processing:",LAB,"------------"))
  datfilenm <- paste0(LAB,".rds")
  RESULT <- readRDS(paste0(LAB,".rds"))   # read in output from  "tracer_slope_across_drives_per_gridcell.r"
  slope <- RESULT$slope
  p.value <- RESULT$p.value
  DX <- unique(signif(diff(as.numeric(rownames(slope))),4))[1]
  DY <- unique(signif(diff(as.numeric(colnames(slope))),4))[1]

  # create data object that can be used in 'facet_grid'

  ZLIMS <- NULL
  if(VAR1=="co2d_ppm_ex"){
    if(VAR2%in%c("nox_ppb_ex"))ZLIMS <- c(0,3)
    if(VAR2=="bc_ngm3_ex")ZLIMS <- c(0,500)
    if(VAR2=="pm25_ugm3_ex")ZLIMS <- c(0,1.0)
    #if(VAR2=="co_ppb_ex")ZLIMS <- c(0,300)
    if(VAR2=="ch4d_ppm_ex")ZLIMS <- c(0,0.004)
  } # if(VAR1=="co2d_ppm_ex"){

  ###############
  # (1) ALL DATA
  xsub <- paste0("All Data; p<",p.value.threshold)
  SEL0 <- p.value < p.value.threshold
  SEL0 <- SEL0&(slope>ZLIMS[1]&slope<ZLIMS[2])
  DAT <- NULL
  data.table <- NULL  # statistics table
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(slope)[1] | dim(vulcan)[2]!=dim(slope)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer slope grid")
  SEL <- SEL0&(vulcan > XLIMS[1] & vulcan < XLIMS[2])
  if(sum(SEL,na.rm=T)<3){print(paste(VULCANvar,"not enough data"));next}
  #  create factor that separates high vs low contributions from VULCAN sector
  hi.perctile <- 0.90
  xquant <- quantile(vulcan[SEL],probs=c(0,hi.perctile,1),na.rm=T)
  hi.low <- cut(vulcan[SEL],breaks=xquant,include.lowest=TRUE)
  tmp <- data.frame(as.vector(slope[SEL]),VULCANvar,as.vector(vulcan[SEL]),hi.low)
  DAT <- rbind(DAT,tmp)
  #  calculate statistics table
  xcor <- cor.test(slope[SEL],vulcan[SEL], method="pearson")       # test for correlation, and calculate correlation coefficient
  R <- xcor$estimate
  p.value <- xcor$p.value
  xlm <- lmodel2(formula=slope[SEL] ~ vulcan[SEL])
  sel <- xlm$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  m <- xlm$regression.results[sel,"Slope"]
  interc <- xlm$regression.results[sel,"Intercept"]
  N <- xlm$n
  #tmp <- bootlmodel2(x=vulcan[SEL],y=slope[SEL])
  #m.sd <- sqrt(var(tmp$m))   # error in regression slope from bootstrap
  m.sd <- NA
  data.table <- rbind(data.table,data.frame(VULCANvar,R,p.value,m,m.sd,interc,N))
} # for(j in 1:length(VULCANvars))
  if(is.null(DAT)){print("not enough data; skip plotting");next}
  colnames(DAT) <- c("slope","VULCANsector","CO2_ppm","hi.low")
  colnames(data.table)[1] <- "VULCANsector"
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  #levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors
  data.table$VULCANsector <- as.factor(data.table$VULCANsector)
  #levels(data.table$VULCANsector) <- VULCANvars

  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=slope)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(XLIMS)
  g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = 0.5*XLIMS[2], y = 0.8*ZLIMS[2], label = paste0('R=', round(R, 2), '\n',
                      'p.value=',signif(p.value,3),'\n','Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
  figfilenm <- paste0(LAB,"_vs_VULCAN.png")
  ggsave(figfilenm);print(paste(figfilenm,"generated"))

  # add box-and-whiskers plot
  theme_set(theme_classic())
  g <- ggplot(dat, aes(VULCANsector, slope))
  g + geom_boxplot(aes(fill=as.factor(hi.low))) + 
      theme(axis.text.x = element_text(angle=65,vjust=0.6))
  figfilenm <- paste0(LAB,"_vs_VULCAN_box.png")
  ggsave(figfilenm);print(paste(figfilenm,"generated"))


  ###############
  # (2) Use BC to filter OUT MAJOR ROADWAYS (e.g., highways)
  BCthresh <- 1000; BCdat <- readRDS("bc_ngm3_ex.rds")
  SEL.BC <- BCdat > BCthresh
  SEL0 <- !SEL.BC & (p.value < p.value.threshold)
  SEL0 <- SEL0&(slope>ZLIMS[1]&slope<ZLIMS[2])
  xsub <- paste0("bc_ngm3_ex <= ",BCthresh,"; p<",p.value.threshold)
  DAT <- NULL
  data.table <- NULL  # statistics table
for(j in 1:length(VULCANvars)){
  VULCANvar <- VULCANvars[j]
  print(VULCANvar)
  vulcan <- readRDS(file=paste0(VULCANvar,".rds"))
  if(dim(vulcan)[1]!=dim(slope)[1] | dim(vulcan)[2]!=dim(slope)[2])stop("dimensions not the same between STILT-VULCAN grid and tracer slope grid")
  SEL <- SEL0&(vulcan > XLIMS[1] & vulcan < XLIMS[2])
  if(sum(SEL,na.rm=T)<3){print(paste(VULCANvar,"not enough data"));next}
  #  create factor that separates high vs low contributions from VULCAN sector
  hi.perctile <- 0.90
  xquant <- quantile(vulcan[SEL],probs=c(0,hi.perctile,1),na.rm=T)
  hi.low <- cut(vulcan[SEL],breaks=xquant,include.lowest=TRUE)
  tmp <- data.frame(as.vector(slope[SEL]),VULCANvar,as.vector(vulcan[SEL]),hi.low)
  DAT <- rbind(DAT,tmp)
  #  calculate statistics table
  xcor <- cor.test(slope[SEL],vulcan[SEL], method="pearson")       # test for correlation, and calculate correlation coefficient
  R <- xcor$estimate
  p.value <- xcor$p.value
  xlm <- lmodel2(formula=slope[SEL] ~ vulcan[SEL])
  sel <- xlm$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  m <- xlm$regression.results[sel,"Slope"]
  interc <- xlm$regression.results[sel,"Intercept"]
  N <- xlm$n
  #tmp <- bootlmodel2(x=vulcan[SEL],y=slope[SEL])
  #m.sd <- sqrt(var(tmp$m))   # error in regression slope from bootstrap
  m.sd <- NA
  data.table <- rbind(data.table,data.frame(VULCANvar,R,p.value,m,m.sd,interc,N))
} # for(j in 1:length(VULCANvars))
  if(is.null(DAT)){print("not enough data; skip plotting");next}
  colnames(DAT) <- c("slope","VULCANsector","CO2_ppm","hi.low")
  colnames(data.table)[1] <- "VULCANsector"
  dat <- DAT
  dat$VULCANsector <- as.factor(dat$VULCANsector)
  #levels(dat$VULCANsector) <- VULCANvars  # arrange order of sectors
  data.table$VULCANsector <- as.factor(data.table$VULCANsector)
  #levels(data.table$VULCANsector) <- VULCANvars

  dev.new()
  theme_set(theme_bw())
  g <- ggplot(dat, aes(x=CO2_ppm,y=slope)) + geom_point(,size=1.5) +
       stat_smooth(method = "lm") # + facet_grid(. ~ VULCANsector)
  g <- g + labs(title=paste0(LAB),caption=paste0("dx=",DX,"; dy=",DY))
  g <- g + labs(subtitle=xsub)
  g <- g + theme(strip.text.x = element_text(size = 16, colour = "black"), axis.title.y=element_text(size=16), axis.title.x=element_text(size=16),
               axis.text.x=element_text(size=10), axis.text.y=element_text(size=14),plot.caption=element_text(size=12),title=element_text(size=16))
  g <- g + ylim(ZLIMS) + xlim(XLIMS)
  g <- g + geom_text(data = data.table, hjust = 0, col='blue',
                  aes(x = 0.5*XLIMS[2], y = 0.8*ZLIMS[2], label = paste0('R=', round(R, 2), '\n',
                      'p.value=',signif(p.value,3),'\n','Slope=', signif(m, 3),'+/-',signif(m.sd,3), '\n','N=',N)))
  g + facet_wrap(.~VULCANsector, nrow = 2, scales='free')
  figfilenm <- paste0(LAB,"_vs_VULCAN_BCfilter.png")
  ggsave(figfilenm);print(paste(figfilenm,"generated"))

  # add box-and-whiskers plot
  theme_set(theme_classic())
  g <- ggplot(dat, aes(VULCANsector, slope))
  g + geom_boxplot(aes(fill=as.factor(hi.low))) + 
      theme(axis.text.x = element_text(angle=65,vjust=0.6))
  figfilenm <- paste0(LAB,"_vs_VULCAN_BCfilter_box.png")
  ggsave(figfilenm);print(paste(figfilenm,"generated"))

} # 

xfiles <- list.files(pattern="_vs_VULCAN.png")
file.copy(from=xfiles,to=outputdir,overwrite=TRUE)



