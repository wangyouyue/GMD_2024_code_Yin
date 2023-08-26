###########################################################################
###
### Program to calculate the vertical profiles of variables (0D)
###
############################################################################
###
### History: 180310, S.Takahashi, initital version
###
############################################################################


############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters
G <- 9.81 #gravitational acceleration[m/s^2]

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
CZ <- ncvar_get(ncin,"CZ")

##### Close
nc_close(ncin)

############################################################################
########## Read History Files
cat(sprintf("start reading history files\n"))

QCALL <- NULL
QRALL <- NULL
QVALL <- NULL
UALL <- NULL
VALL <- NULL
WALL <- NULL
PTALL <- NULL

##### loop of MPI rank
for(mpirank in allmpiranks){
    cat(sprintf("processing the rank = %s \n",mpirank))

    ##### Open
    ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

    ##### Read QC
    QC <- ncvar_get(ncin,"QC_sd")
    QC_units <- ncatt_get(ncin,"QC_sd","units")
    dimnames(QC)[[4]] <- alltimes
    if(!setequal(dim(QC[,,,1]),c(IMAX,JMAX,KMAX))){
        cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
        return()
    }
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

    ##### Read PT
    PT <- ncvar_get(ncin,"PT")
    PT_units <- ncatt_get(ncin,"PT","units")
    dimnames(PT)[[4]] <- alltimes
    PTALL <- append(PTALL,PT)

    ##### Close
    nc_close(ncin)
}

############################################################################
########## QL
cat(sprintf("start plotting QL\n"))

##### calculate parameters
QCALL_ARY <- array(QCALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
QRALL_ARY <- array(QRALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(QCALL)
rm(QRALL)

QLALL_ARY <- QCALL_ARY + QRALL_ARY
QL_AVE <- matrix(0.0,nrow=KMAX, ncol=length(alltimes))
for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        QL_AVE[k,tm] <- mean(QLALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
    }
}

# convert kg/kg to g/kg
QL_AVE <- QL_AVE * 1000

##### plot parameters
pdf("QL_vprof.pdf")

for(tm in 1:length(alltimes)){
    plot(
        QL_AVE[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main="Vertical Profile of QL",
        xlim=c(0, max(QL_AVE)),
        xlab="QL[g/kg]",
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)
}

dev.off()

cat(sprintf("end plotting QL\n\n"))

##### remove used objects (QCALL_ARY,QLALL_ARY will use later.)
rm(QRALL_ARY)
rm(QL_AVE)

############################################################################
########## QT
cat(sprintf("start plotting QT\n"))

##### calculate parameters
QVALL_ARY <- array(QVALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(QVALL)

QTALL_ARY <- QLALL_ARY + QVALL_ARY
QT_AVE <- matrix(0.0,nrow=KMAX, ncol=length(alltimes))
for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        QT_AVE[k,tm] <- mean(QTALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
    }
}

# convert kg/kg to g/kg
QT_AVE <- QT_AVE * 1000

##### plot parameters
pdf("QT_vprof.pdf")

for(tm in 1:length(alltimes)){
    plot(
        QT_AVE[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main="Vertical Profile of QT",
        xlim=c(0, max(QT_AVE)),
        xlab="QT[g/kg]",
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)
}

dev.off()

cat(sprintf("end plotting QT\n\n"))

##### remove used objects (QCALL_ARY will use later.)
rm(QVALL_ARY)
rm(QLALL_ARY)
rm(QTALL_ARY)
rm(QT_AVE)

############################################################################
########## cfrac
cat(sprintf("start plotting cfrac\n"))

##### calculate parameters
CFRAC  <- matrix(0.0,nrow=KMAX, ncol=length(alltimes))
for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        for(i in 1:IMAX){
            for(j in 1:JMAX){
                for(rk in 1:MPINUM){
                    if(QCALL_ARY[i,j,k,tm,rk] * 1000 > 0.01){
                        CFRAC[k,tm] <- CFRAC[k,tm] + 1
                    }
                }
            }
        }
    }
}
CFRAC <- CFRAC / (IMAX*JMAX*MPINUM)

##### plot parameters
pdf("cfrac_vprof.pdf")

for(tm in 1:length(alltimes)){
    plot(
        CFRAC[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main="Vertical Profile of cfrac",
        xlim=c(0, max(CFRAC)),
        xlab="cfrac",
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)
}

dev.off()
cat(sprintf("end plotting cfrac\n\n"))

##### remove used objects
rm(QCALL_ARY)
rm(CFRAC)

############################################################################
########## TKE
cat(sprintf("start plotting TKE\n"))

##### calculate parameters
UALL_ARY <- array(UALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
VALL_ARY <- array(VALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
WALL_ARY <- array(WALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),MPINUM))
rm(UALL)
rm(VALL)
rm(WALL)

U_AVE <- matrix(0.0,KMAX,length(alltimes))
U_VAL <- matrix(0.0,KMAX,length(alltimes))
V_AVE <- matrix(0.0,KMAX,length(alltimes))
V_VAL <- matrix(0.0,KMAX,length(alltimes))
W_AVE <- matrix(0.0,KMAX,length(alltimes))
W_VAL <- matrix(0.0,KMAX,length(alltimes))

for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        U_AVE[k,tm] <- mean(UALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
        V_AVE[k,tm] <- mean(VALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
        W_AVE[k,tm] <- mean(WALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
    }
}

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
U_VAL <- U_VAL / (IMAX * JMAX * MPINUM)
V_VAL <- V_VAL / (IMAX * JMAX * MPINUM)
W_VAL <- W_VAL / (IMAX * JMAX * MPINUM)

TKE <- (U_VAL + V_VAL + W_VAL) / 2.0

##### plot parameters
pdf("tke_vprof.pdf")

for(tm in 1:length(alltimes)){
    par(mgp=c(2.0,1.0,0))

    plot(
        TKE[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main="Vertical Profile of TKE",
        xlim=c(0,max(TKE)),
        xlab=expression(paste("TKE [",m^2,"/",s^2,"]")),
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)

}
dev.off()

cat(sprintf("end plotting TKE\n\n"))

##### remove used objects (W_ALL_ARY and W_AVE will use later.)
rm(UALL_ARY)
rm(VALL_ARY)
rm(U_AVE)
rm(V_AVE)
rm(U_VAL)
rm(V_VAL)
rm(W_VAL)
rm(TKE)

############################################################################
########## W
cat(sprintf("start plotting W\n"))

##### plot parameters
pdf("W_vprof.pdf")

for(tm in 1:length(alltimes)){
    plot(
        W_AVE[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main="Vertical Profile of W",
        xlim=c(min(W_AVE), max(W_AVE)),
        xlab="W [m/s]",
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)
}
dev.off()

cat(sprintf("end plotting W\n\n"))

# W_ALL_ARY and W_AVE will use later.

############################################################################
########## Potential Temperature
cat(sprintf("start plotting Potential Temperature\n"))

##### calculate parameters
PTALL_ARY <- array(PTALL,dim=c(IMAX,JMAX,KMAX,length(alltimes),length(allmpiranks)))
rm(PTALL)

PT_AVE <- matrix(0.0,nrow=KMAX, ncol=length(alltimes))
for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        PT_AVE[k,tm] <- mean(PTALL_ARY[1:IMAX,1:JMAX,k,tm,1:MPINUM])
    }
}

##### plot parameters
pdf("PT_vprof.pdf")

for(tm in 1:length(alltimes)){
    plot(
        PT_AVE[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main=expression(paste("Vertical Profile of ",theta)),
        xlim=c(min(PT_AVE), max(PT_AVE)),
        xlab=expression(paste(theta," [K]")),
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)
}
dev.off()
cat(sprintf("end plotting Potential Temperature\n\n"))

# PTALL_ARY and PT_AVE will use later.

############################################################################
########## tot_boy
cat(sprintf("start plotting tot_boy\n"))

##### calculate parameters
TOT_BOY <- matrix(0.0,nrow=KMAX, ncol=length(alltimes))
for(tm in 1:length(alltimes)){
    for(k in 1:KMAX){
        for(rk in 1:MPINUM){
            for(j in 1:JMAX){
                for(i in 1:IMAX){
                    TOT_BOY[k,tm] <- TOT_BOY[k,tm] +
                        (WALL_ARY[i,j,k,tm,rk] - W_AVE[k,tm]) *
                        (PTALL_ARY[i,j,k,tm,rk] - PT_AVE[k,tm]) / PT_AVE[k,tm]
                }
            }
        }
    }
}
TOT_BOY <- TOT_BOY * G * 10000.0 / (IMAX*JMAX*MPINUM)

##### plot parameters
pdf("tot_boy_vprof.pdf")

for(tm in 1:length(alltimes)){
    par(mgp=c(2.0,1.0,0))
    tm_real <- (tm - 1)*60
    plot(
        TOT_BOY[1:KMAX,tm],CZ[3:(KMAX+2)],
        type="l",
        main="Vertical Profile of tot_boy",
        xlim=c(min(TOT_BOY), max(TOT_BOY)),
        xlab=expression(paste("tot_boy[",cm^3,"/",s^2,"]")),
        ylab="Altitude[m]"
    )
    mtext(sprintf("Time = %d [s]",(tm-1)*60),side=3,line=0.25)
}
dev.off()

cat(sprintf("end plotting tot_boy\n\n"))
