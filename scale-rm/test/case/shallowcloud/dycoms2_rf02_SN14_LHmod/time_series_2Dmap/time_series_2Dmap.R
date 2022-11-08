############################################################################
###
### Program to calculate the time series of 2D map
### (for single run)
###
############################################################################
###
### History: 180318, S.Takahashi,   initital version
###
############################################################################


############################################################################
########## Loading Libraries
library(ncdf4)
library(fields)

############################################################################
########## Parameters

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

##### Grids
IMAX <- length(ncvar_get(ncin,"CX"))-4
JMAX <- length(ncvar_get(ncin,"CY"))-4
KMAX <- length(ncvar_get(ncin,"CZ"))-4
CDX <- ncvar_get(ncin,"CDX")
CDY <- ncvar_get(ncin,"CDY")
CDZ <- ncvar_get(ncin,"CDZ")
CAREA <- CDX[3] * CDY[3]

##### Area informaion
CXG <- ncvar_get(ncin,"CXG")
CYG <- ncvar_get(ncin,"CYG")
XMAX <- CXG[3] + CXG[length(CXG)-2]
YMAX <- CYG[3] + CYG[length(CYG)-2]
XNUM <- length(CXG)-4
YNUM <- length(CYG)-4

##### Number of processors(Use X-direction only)
XPRC <- XNUM %/% IMAX

##### Close
nc_close(ncin)

############################################################################
########## Read History Files
cat(sprintf("start reading history files\n"))
DENSALL <- NULL
QCALL <- NULL
QRALL <- NULL

##### loop of MPI rank
for(mpirank in allmpiranks){

    cat(sprintf("processing the rank = %s \n",mpirank))

    ##### Open history file
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

    ##### Read QC
    QC <- ncvar_get(ncin,"QC")
    QC_units <- ncatt_get(ncin,"QC","units")
    dimnames(QC)[[4]] <- alltimes
    QCALL <- append(QCALL,QC)

    ##### Read QR
    QR <- ncvar_get(ncin,"QR")
    QR_units <- ncatt_get(ncin,"QR","units")
    dimnames(QR)[[4]] <- alltimes
    QRALL <- append(QRALL,QR)

    ##### Close
    nc_close(ncin)
}

############################################################################
########## LWP,CWP,RWP
cat(sprintf("start plotting water path\n"))

##### calculate parameters
DENSALL_ARY <- array(DENSALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(DENSALL)

QCALL_ARY <- array(QCALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(QCALL)

QRALL_ARY <- array(QRALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(QRALL)

TOT_CLOUD_MASS <- array(0, dim=c(XNUM, YNUM, length(alltimes)))
TOT_RAIN_MASS  <- array(0, dim=c(XNUM, YNUM, length(alltimes)))

for(rk in 1:MPINUM){
    pno <- rk - 1
    for(tm in 1:length(alltimes)){

        CLOUD_MASS <- QCALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                      DENSALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                      (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])

        RAIN_MASS <- QRALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                     DENSALL_ARY[1:IMAX,1:JMAX,1:KMAX,tm,rk] *
                     (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)])

        for(i in 1:IMAX){
            x <- (pno %% XPRC)*IMAX + i
            for(j in 1:JMAX){
                y <- (pno %/% XPRC)*JMAX + j
                TOT_CLOUD_MASS[x,y,tm] <- sum(CLOUD_MASS[i,j,1:KMAX])
                TOT_RAIN_MASS[x,y,tm] <- sum(RAIN_MASS[i,j,1:KMAX])
            }
        }
    }
}
CWP <- TOT_CLOUD_MASS/CAREA
RWP <- TOT_RAIN_MASS/CAREA
LWP <- CWP + RWP
zmax <- max(LWP)

##### plot parameters
pdf("lwp_2Dmap_tmsrs.pdf")

for(tm in 1:length(alltimes)){
    image.plot(
        CXG[3:(XNUM+2)],CYG[3:(YNUM+2)],LWP[,,tm],
        bty="o",
        xlim=c(0, XMAX),
        ylim=c(0, YMAX),
        zlim=c(0, zmax),
        xlab="X[m]",
        ylab="Y[m]",
        main=expression(paste("Horizontal Distribution of LWP[kg/",m^2,"]")),
        sub=sprintf("Time = %d [s]",(tm-1)*60),
        col=gray(seq(0, 1, 1/128))
    )
}

dev.off()

cat(sprintf("end plotting water path\n\n"))
