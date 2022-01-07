# What would error statistics mean in terms of OFFSET in terms of Monte Carlo source determination?  
# Jan. 6th, 2022 by John C. Lin (John.Lin@utah.edu)

#####################
starts <- c("201908010000","201910010000","201912010000")
ends <-   c("201909302300","201911302300","202002292300")
Output <- "/uufs/chpc.utah.edu/common/home/u0791084/PROJECTS/Lin_Uintah_CH4_trends_2021/Transporterr_stiltread/Output/"
L <- seq(0,1000)   # lengthscale [m]
#####################

if(length(starts)!=length(ends))stop("starts and ends not the same length")

# Key MesoWest sites in SLV recommended by Alex Jacques
stids <- c("FPN","FPS","HERUT","MTMET","NAA","NHMU","SUNUT","TRJO","UFD09","WBB","FARM","TPC","K36U","KSLC",
           "KTVY","KU42","QBV","QED","QH3","QHW","QLN","QMG","QNP","QRP","QSA","BAC","CEN","KIJ","UT11","UT12",
           "UT20","UT201","UT23","UT248","UT3","UT5","UT7","UT9","UTALP","UTBIG","UTCDF","UTCHL","UTCOL","UTDAN",
           "UTDCD","UTDWY","UTHEB","UTJUN","UTLGP","UTLPC","UTMFS","UTORM","UTPCR","UTPKL","UTPLC","UTPR4","UTQRY",
           "UTSTR","UTSVC","UTTPD","UTWAN","UTWLC","UCC14")
# Additional MesoWest sites in northern SLC
stids <- c(stids,c("USDR2"))  #UofU MiniSodar2
mettype <- "HRRR"

colnms <- expand.grid(c("Uobs.","Vobs."),stids)
colnms <- paste0(colnms[,1],colnms[,2])

Wspd.all <- NULL; err.all <- NULL
for(i in 1:length(starts)){
  start <- starts[i]; end <- ends[i]
  # observed and simulated meteorology 
  resultname <- paste0(Output,"/",mettype,"_obs_",start,"to",end,".RDS")
  obs <- readRDS(resultname)$dat
  obs <- obs[,colnames(obs)%in%colnms]
  stations <- unique(substring(colnames(obs),6,nchar(colnames(obs))))
  Wspd.sub <- matrix(NA,nrow=nrow(obs),ncol=length(stations))
  colnames(Wspd.sub) <- stations
  for(j in 1:length(stations)){
    Uobs <- obs[,paste0("Uobs.",stations[j])]
    Vobs <- obs[,paste0("Vobs.",stations[j])]
    Wspd <- sqrt(Uobs^2+Vobs^2)
    Wspd.sub[,j] <- Wspd
  } # for(j in 1:length(stations)){
  Wspd.ave <- mean(Wspd.sub,na.rm=T)
  Wspd.all <- c(Wspd.all,Wspd.ave)

  # transport error statistics
  resultname <- paste0(Output,"/",mettype,"_Errstats_",start,"to",end,".RDS")
  tmp <- readRDS(resultname)$dat
  dat <- tmp[rownames(tmp)%in%stids,]
  err.ave <- apply(dat,2,mean,na.rm=T)
  err.all <- rbind(err.all,err.ave)
} # for(i in 1:length(starts)){
rownames(err.all) <- paste0(starts,"_to_",ends)
names(Wspd.all) <- paste0(starts,"_to_",ends)
print(err.all)
print(Wspd.all)
print(mean(err.all[,"U.rmse"]))
print(mean(Wspd.all))

# error in source location due to transport error
err.L <- mean(err.all[,"U.rmse"])*(L/mean(Wspd.all))
plot(L,err.L,xlab="Distance to Source [m]",ylab="Error in Source Location [m]")