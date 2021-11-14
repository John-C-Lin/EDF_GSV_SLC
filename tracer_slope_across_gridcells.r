# Calculate tracer:tracer CORRELATIONS and REGRESSION SLOPES (not RATIOS) over different gridcells 
# October 18th, 2021 by John C. Lin (John.Lin@utah.edu)

#####################
# calculate correlations between VARS.2 (y) versus VARS.1 (x)
#VARs <- c("pm25_ugm3_ex","co_ppb_ex","bc_ngm3_ex","no_ppb_ex","no2_ppb_ex","nox_ppb_ex","co2d_ppm_ex","ch4d_ppm_ex",
#          "onroad","industrial","residential","commercial","total")
VARs.2 <- c("nox_ppb_ex","bc_ngm3_ex","pm25_ugm3_ex","ch4d_ppm_ex")[1]
VARs.1 <- c("co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex","co2d_ppm_ex")[1]
#####################

if(length(VARs.1)!=length(VARs.2))stop("length of VARs.1 and VARs.2 need to be the same")

regress.SMA <- function(x,y){
  require(lmodel2)
  xfit<-lmodel2(formula=y~x)
  sel<-xfit$regression.results[,"Method"]=="SMA" #from Nick Murdoch's test code, appears that the SMA method (standard major axis regression)
  result <- xfit$regression.results[sel,]
  return(result)
} #f<-(x,y)


####################
#  ALL the data
xsub <- "All Data"
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1; y <- vdat2
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_1.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()


####################
#  Use BC to select for MAJOR ROADWAYS (e.g., highways)
BCthresh <- 1000
BCdat <- readRDS("bc_ngm3_ex.rds")
SEL.BC <- BCdat > BCthresh
SEL <- SEL.BC
xsub <- paste0("bc_ngm3_ex > ",BCthresh)
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1[SEL]; y <- vdat2[SEL]
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_2.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()


####################
#  Use BC to filter OUT major roadways (e.g., highways)
SEL <- !SEL.BC
xsub <- paste0("bc_ngm3_ex <= ",BCthresh)
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1[SEL]; y <- vdat2[SEL]
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_3.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()


####################
#  Use BC to filter OUT major roadways (e.g., highways) AND Onroad CO2 % is > threshold  of total modeled CO2
PercOnroad.min <- 80
xsub <- paste0("bc_ngm3_ex <= ",BCthresh," AND PercOnroad >= ",PercOnroad.min,"%")
PercOnroad <- 100*readRDS("onroad_ppm_ex.rds")/readRDS("total_ppm_ex.rds")
SEL <- !SEL.BC & (PercOnroad >= PercOnroad.min)
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1[SEL]; y <- vdat2[SEL]
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_4.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()


####################
#  Use BC to filter OUT major roadways (e.g., highways) AND Industrial CO2 % is > threshold  of total modeled CO2
PercIndustrial.min <- 50
xsub <- paste0("bc_ngm3_ex <= ",BCthresh," AND PercIndustrial >= ",PercIndustrial.min,"%")
PercIndustrial <- 100*readRDS("industrial_ppm_ex.rds")/readRDS("total_ppm_ex.rds")
SEL <- !SEL.BC & (PercIndustrial >= PercIndustrial.min)
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1[SEL]; y <- vdat2[SEL]
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_5.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()

####################
#  Use BC to filter OUT major roadways (e.g., highways) AND Resdiential CO2 % is > threshold  of total modeled CO2
PercResidential.min <- 10
xsub <- paste0("bc_ngm3_ex <= ",BCthresh," AND PercResidential >= ",PercResidential.min,"%")
PercResidential <- 100*readRDS("residential_ppm_ex.rds")/readRDS("total_ppm_ex.rds")
SEL <- !SEL.BC & (PercResidential >= PercResidential.min)
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1[SEL]; y <- vdat2[SEL]
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_6.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()

####################
#  Use BC to filter OUT major roadways (e.g., highways) AND Commercial CO2 % is > threshold  of total modeled CO2
PercCommercial.min <- 10
xsub <- paste0("bc_ngm3_ex <= ",BCthresh," AND PercCommercial >= ",PercCommercial.min,"%")
PercCommercial <- 100*readRDS("commercial_ppm_ex.rds")/readRDS("total_ppm_ex.rds")
SEL <- !SEL.BC & (PercCommercial >= PercCommercial.min)
for(i in 1:length(VARs.1)){
  var1 <- VARs.1[i]; var2 <- VARs.2[i]
  # read in output from "average_over_drives.r"
  vdat1 <- readRDS(paste0(var1,".rds"))
  vdat2 <- readRDS(paste0(var2,".rds"))
  x <- vdat1[SEL]; y <- vdat2[SEL]
  isNA <- is.na(x)|is.na(y)
  x <- x[!isNA];y<-y[!isNA]
  xcor <- cor(x,y,use="pairwise.complete.obs",method="pearson")
  xfit <- regress.SMA(x=x,y=y)
  slope <- as.numeric(xfit["Slope"]); intercept <- as.numeric(xfit["Intercept"])

  dev.new();par(cex.axis=1.3,cex.main=1.3,cex.lab=1.3,cex.sub=1.3)
  plot(x,y,pch=16,xlab=var1,ylab=var2)
  lines(x,slope*x + intercept, col="red")
  #xmain <- paste0("R=",signif(xcor,3),"; Slope=",signif(slope,4),"\n Intercept=",signif(intercept,4))
  xmain <- paste0("R=",signif(xcor,3),"; y=",signif(slope,4),"*x + ",signif(intercept,4))
  title(main=xmain,sub=xsub)
  dev.copy(png,"tracer_slope_across_gridcells_7.png");dev.off()
} # for(i in 1:length(VARs.1)){
gc()


