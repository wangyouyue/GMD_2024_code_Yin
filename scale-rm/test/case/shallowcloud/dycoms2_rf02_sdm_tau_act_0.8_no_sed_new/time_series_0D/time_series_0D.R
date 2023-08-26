############################################################################
###
### Program to calculate the time series of domain avaraged variables (0D)
###
############################################################################
###
### History: 171218, S.Shima,     initital version
###          180330, S.Takahashi  add following parameters
###                                 TKE
###                                 Max W Variance
###                                 cfrac
###                                 Entrainment Rate
###
############################################################################


############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters
TIME_SKIP <- 300.0

############################################################################
########## Read Data from Files
##### Make a list of files,  mpiranks
allfiles = dir("../",pattern="^history.")
tmp = strsplit(allfiles,"\\history.pe|\\.nc")
allmpiranks = unique(matrix(unlist(tmp),nrow=2)[2,])
names(allfiles) = allmpiranks
MPINUM <- length(allmpiranks)

##### Open the first file
ncin <- nc_open(paste("../",allfiles[1],sep=""))

##### Times
alltimes <- ncvar_get(ncin,"time")
time_units <- ncatt_get(ncin,"time","units")
alltimes_h <- alltimes / 3600.0

##### Grids
IMAX <- length(ncvar_get(ncin,"CX"))-4
JMAX <- length(ncvar_get(ncin,"CY"))-4
KMAX <- length(ncvar_get(ncin,"CZ"))-4
CDX <- ncvar_get(ncin,"CDX")
CDY <- ncvar_get(ncin,"CDY")
CDZ <- ncvar_get(ncin,"CDZ")
CZ <- ncvar_get(ncin,"CZ")
AREA <- sum(CDX[3:(IMAX+2)])*sum(CDY[3:(JMAX+2)])*MPINUM

##### Close
nc_close(ncin)

############################################################################
########## Read History Files
cat(sprintf("start reading history files\n"))
UALL <- NULL
VALL <- NULL
WALL <- NULL
QCALL <- NULL
QRALL <- NULL
QVALL <- NULL
DENSALL <- NULL
RAINALL <- NULL

##### loop of MPI rank
for(mpirank in allmpiranks){
    cat(sprintf("processing the rank = %s \n",mpirank))

    ##### Open
    ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

    ##### Read DENS
    DENS <- ncvar_get(ncin,"DENS")
    DENS_units <- ncatt_get(ncin,"DENS","units")
    dimnames(DENS)[[4]] <- alltimes
    if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
        cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
        return()
    }
    DENSALL <- append(DENSALL,DENS)

    ##### Read RAIN_sd
    RAIN <- ncvar_get(ncin,"RAIN_ACC_sd")
    RAIN_units <- ncatt_get(ncin,"RAIN_ACC_sd","units")
    dimnames(RAIN)[[3]] <- alltimes
    RAINALL <- append(RAINALL,RAIN)

    ##### Read QC
    QC <- ncvar_get(ncin,"QC_sd")
    QC_units <- ncatt_get(ncin,"QC_sd","units")
    dimnames(QC)[[4]] <- alltimes
    QCALL <- append(QCALL,QC)

    ##### Read QR
    QR <- ncvar_get(ncin,"QR_sd")
    QR_units <- ncatt_get(ncin,"QR_sd","units")
    dimnames(QR)[[4]] <- alltimes
    QRALL <- append(QRALL,QR)

    ##### Read QV
    QV <- ncvar_get(ncin,"QV")
    QV_units <- ncatt_get(ncin,"QV","units")
    dimnames(QV)[[4]] <- alltimes
    QVALL <- append(QVALL,QV)

    ##### Read U
    U <- ncvar_get(ncin,"U")
    U_units <- ncatt_get(ncin,"U","units")
    dimnames(U)[[4]] <- alltimes
    UALL <- append(UALL,U)

    ##### Read V
    V <- ncvar_get(ncin,"V")
    V_units <- ncatt_get(ncin,"V","units")
    dimnames(V)[[4]] <- alltimes
    VALL <- append(VALL,V)

    ##### Read W
    W <- ncvar_get(ncin,"W")
    W_units <- ncatt_get(ncin,"W","units")
    dimnames(W)[[4]] <- alltimes
    WALL <- append(WALL,W)

    ##### Close
    nc_close(ncin)
}
cat(sprintf("end reading history files\n\n"))

############################################################################
########## TKE
cat(sprintf("start plotting tke\n"))

##### calculate parameters
UALL_ARY <- array(UALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
VALL_ARY <- array(VALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
WALL_ARY <- array(WALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
DENSALL_ARY <- array(DENSALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(UALL)
rm(VALL)
rm(WALL)
rm(DENSALL)

U_AVE <- matrix(0.0,KMAX,length(alltimes))
U_VAL <- matrix(0.0,KMAX,length(alltimes))
V_AVE <- matrix(0.0,KMAX,length(alltimes))
V_VAL <- matrix(0.0,KMAX,length(alltimes))
W_AVE <- matrix(0.0,KMAX,length(alltimes))
W_VAL <- matrix(0.0,KMAX,length(alltimes))
DENS_AVE <- matrix(0.0,KMAX,length(alltimes))

for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        U_AVE[k,tm] <- sum(UALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
        V_AVE[k,tm] <- sum(VALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
        W_AVE[k,tm] <- sum(WALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
        DENS_AVE[k,tm] <- sum(DENSALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
    }
}
U_AVE <- U_AVE / (IMAX * MPINUM * JMAX)
V_AVE <- V_AVE / (IMAX * MPINUM * JMAX)
W_AVE <- W_AVE / (IMAX * MPINUM * JMAX)
DENS_AVE <- DENS_AVE / (IMAX * MPINUM * JMAX)

for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        for(rk in 1:MPINUM){
            for(j in 1:JMAX){
                for(i in 1:IMAX){
                    U_VAL[k,tm] <- U_VAL[k,tm] + (U_AVE[k,tm] - UALL_ARY[i,j,k,tm,rk])^2
                    V_VAL[k,tm] <- V_VAL[k,tm] + (V_AVE[k,tm] - VALL_ARY[i,j,k,tm,rk])^2
                    W_VAL[k,tm] <- W_VAL[k,tm] + (W_AVE[k,tm] - WALL_ARY[i,j,k,tm,rk])^2
                }
            }
        }
    }
}
U_VAL <- U_VAL / (IMAX * MPINUM * JMAX)
V_VAL <- V_VAL / (IMAX * MPINUM * JMAX)
W_VAL <- W_VAL / (IMAX * MPINUM * JMAX)

TKE <- rep(0.0,length(alltimes))
for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        TKE[tm] <- TKE[tm] +
               (U_VAL[k,tm] + V_VAL[k,tm] + W_VAL[k,tm]) / 2 * DENS_AVE[k,tm] * CDZ[2+k]
    }
}

##### plot parameters
pdf("tke_tmsrs.pdf")
par(mgp=c(2.0,1.0,0))
plot(
    alltimes_h,
    TKE,
    type="l",
    main="Domain Averaged TKE",
    xlab="Time [h]",
    ylab=expression(paste("TKE [kg/",s^2,"]"))
)

dev.off()

cat(sprintf("end plotting tke\n\n"))

##### remove used objects (except W_VAL, DENSALL_ARY)
rm(UALL_ARY)
rm(VALL_ARY)
rm(WALL_ARY)
rm(U_AVE)
rm(V_AVE)
rm(W_AVE)
rm(DENS_AVE)
rm(U_VAL)
rm(V_VAL)
rm(TKE)

############################################################################
########## Max. W Variance
cat(sprintf("start plotting Max. W valiance\n"))

##### calculate parameters
W_VAL_MAX <- rep(0.0,length(alltimes))
for(tm in 1:length(alltimes)){
    W_VAL_MAX[tm] = max(W_VAL[1:KMAX,tm])
}

##### plot parameters
pdf("w_max_val_tmsrs.pdf")
par(mgp=c(2.0,1.0,0))
plot(
    alltimes_h,
    W_VAL_MAX,
    type="l",
    main="Domain Averaged Max. W Variance",
    xlab="Time [h]",
    ylab=expression(paste("Max. W Variance [",m^2,"/",s^2,"]"))
)

dev.off()

cat(sprintf("end plotting Max. W valiance\n\n"))

##### remove used objects
rm(W_VAL)
rm(W_VAL_MAX)

############################################################################
########## Entrainment Rate
cat(sprintf("start plotting entrainment rate\n"))

##### calculate parameters
QCALL_ARY <- array(QCALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
QRALL_ARY <- array(QRALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
QVALL_ARY <- array(QVALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(QCALL)
rm(QRALL)
rm(QVALL)

QTALL_ARY <- QCALL_ARY + QRALL_ARY + QVALL_ARY
INVZ <- array(0.0,dim=c(IMAX,JMAX,length(alltimes),MPINUM))
for(tm in 1:length(alltimes)){
    for(rk in 1:MPINUM){
        for(i in 1:IMAX){
            for(j in 1:JMAX){
                for(k in 1:KMAX){
                    if(QTALL_ARY[i,j,k,tm,rk] < 8.0e-3){
                        INVZ[i,j,tm,rk] <- CZ[k+2]
                        break
                    }
                }
            }
        }
    }
}

INVZ_AVE <- rep(0.0,length(alltimes))
for(tm in 1:length(alltimes)){
    INVZ_AVE[tm] <- sum(INVZ[1:IMAX,1:JMAX,tm,1:MPINUM])
}
INVZ_AVE <- INVZ_AVE / (IMAX * JMAX * MPINUM)

ENT_RATE <- rep(0.0,length(alltimes))
# Note : ENT_RATE on t=0 is regarded as 0.0.
for(tm in 2:length(alltimes)){
    ENT_RATE[tm] <- (INVZ_AVE[tm] - INVZ_AVE[tm-1]) / TIME_SKIP + 3.75e-6 * INVZ_AVE[tm]
}
ENT_RATE <- ENT_RATE * 100

##### plot parameters
pdf("entrainment_rate_tmsrs.pdf")
par(mgp=c(2.0,1.0,0))
plot(
    alltimes_h,ENT_RATE,
    type="l",
    main="Domain Averaged Entrainment Rate",
    xlab="Time [h]",
    ylab="Entrainment Rate [cm/s]"
)

dev.off()

cat(sprintf("end plotting entrainment rate\n\n"))

##### remove used objects
rm(QTALL_ARY)
rm(QVALL_ARY)

############################################################################
########## LWP,CWP,RWP
cat(sprintf("start plotting water path\n"))

##### calculate parameters
TOT_CLOUD_MASS <- rep(0.0,length(alltimes))
TOT_RAIN_MASS  <- rep(0.0,length(alltimes))

for(rk in 1:MPINUM){
    for(tm in 1:length(alltimes)){

        CLOUD_MASS <- QCALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                      DENSALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])
        TOT_CLOUD_MASS[tm] <- TOT_CLOUD_MASS[tm] + sum(CLOUD_MASS[1:IMAX,1:JMAX,1:KMAX])

        RAIN_MASS <- QRALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                     DENSALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                     (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])
        TOT_RAIN_MASS[tm] <- TOT_RAIN_MASS[tm] + sum(RAIN_MASS[1:IMAX,1:JMAX,1:KMAX])
    }
}
CWP_AVE <- TOT_CLOUD_MASS/AREA
RWP_AVE <- TOT_RAIN_MASS/AREA
LWP_AVE <- CWP_AVE + RWP_AVE

##### plot parameters
pdf("water_path_tmsrs.pdf")

cols   <- c("black","blue", "red")
ltys   <- c(1, 2, 4)
labels <- c("LWP","CWP","RWP")
ymax <- 0.2

par(mgp=c(2.0,1.0,0))
plot(
    alltimes_h,
    LWP_AVE,
    type="l",lty = ltys[1], col = cols[1],
    ylim=c(0.0,ymax),
    main="Domain Averaged Water Path",
    xlab="Time [h]",
    ylab=expression(paste("Water Path [kg/",m^2,"]")),
    ann=F)

par(new=T)
plot(
    alltimes_h,
    CWP_AVE,
    type="l",lty = ltys[2], col = cols[2],
    ylim=c(0.0,ymax),
    main="Domain Averaged Water Path",
    xlab="Time [h]",
    ylab=expression(paste("Water Path [kg/",m^2,"]")),
    ann=F)

par(new=T)
plot(
    alltimes_h,
    RWP_AVE,
    type="l",lty = ltys[3], col = cols[3],
    ylim=c(0.0,ymax),
    main="Domain Averaged Water Path",
    xlab="Time [h]",
    ylab=expression(paste("Water Path [kg/",m^2,"]"))
)

legend("topleft", legend = labels, col = cols, lty = ltys)

dev.off()

cat(sprintf("end plotting water path\n\n"))

##### remove used objects
rm(CWP_AVE)
rm(RWP_AVE)
rm(LWP_AVE)
rm(TOT_CLOUD_MASS)
rm(TOT_RAIN_MASS)
rm(QRALL_ARY)

############################################################################
########## cfrac
cat(sprintf("start plotting cfrac\n"))

##### calculate parameters
TOT_CFRAC  <- rep(0.0,length(alltimes))

for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        for(j in 1:JMAX){
            for(i in 1:IMAX){
                for(rk in 1:MPINUM){
                    if(QCALL_ARY[i,j,k,tm,rk] * 1000 > 0.01){
                        TOT_CFRAC[tm] <- TOT_CFRAC[tm] + 1;
                    }
                }
            }
        }
    }
}
TOT_CFRAC <- TOT_CFRAC / (IMAX * JMAX * KMAX * MPINUM)

##### plot parameters
pdf("cfrac_tmsrs.pdf")
par(mgp=c(2.0,1.0,0))
plot(
    alltimes_h,
    TOT_CFRAC,
    type="l",
    main="Domain Averaged cfrac",
    xlab="Time [h]",
    ylab="cfrac"
)

dev.off()

cat(sprintf("end plotting cfrac\n\n"))

##### remove used objects
rm(TOT_CFRAC)
rm(QCALL_ARY)

############################################################################
########## Plot Domain Averaged Precipitation Rate
cat(sprintf("start plotting precipitation rate\n"))

##### calculate parameters
RAINALL_ARY <- array(RAINALL,dim=c(IMAX,JMAX,length(alltimes),MPINUM))
rm(RAINALL)

TOT_ACC_RAIN <- rep(0.0,length(alltimes))
for(rk in 1:MPINUM){
    for(tm in 1:length(alltimes)){
        ACC_RAIN <-  RAINALL_ARY[1:IMAX,1:JMAX,tm,rk] * (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)])
        TOT_ACC_RAIN[tm] <- TOT_ACC_RAIN[tm] + sum(ACC_RAIN[1:IMAX,1:JMAX])
    }
}
ACC_RAIN_AVE <- TOT_ACC_RAIN/AREA

##### CONVERSION to Rain rate [mm/day]
RATE_RAIN_AVE <- rep(0.0,length(alltimes))
RATE_RAIN_AVE[1] <- 0.0 
for(tm in 2:length(alltimes)){
        RATE_RAIN_AVE[tm] <- (ACC_RAIN_AVE[tm]-ACC_RAIN_AVE[tm-1])/(alltimes[tm]-alltimes[tm-1])*3600.0*24.0*1000.0
}

##### plot parameters
pdf("precipitation_tmsrs.pdf")
par(mgp=c(2.0,1.0,0))
plot(
    alltimes_h,
    RATE_RAIN_AVE,
    type="o",
    main="Domain Averaged Precipitation Rate",
    xlab="Time [h]",
    ylab="Precipitation Rate [mm/day]"
)

dev.off()

cat(sprintf("end plotting precipitation rate\n\n"))

