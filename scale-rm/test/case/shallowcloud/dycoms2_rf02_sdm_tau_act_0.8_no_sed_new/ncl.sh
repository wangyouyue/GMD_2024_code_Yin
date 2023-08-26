#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=sa
#PJM -L node=8
#PJM --mpi proc=320
#PJM -L elapse=24:00:00
#PJM -g r22775
#PJM -j
#------- Program execution -------#
module load intel impi hdf5 netcdf netcdf-fortran
source ~/.bashrc
conda activate ncl_stable
# run
rm -rf SDM_history_*.ncl
MPI=7680
core_num=320
let p=MPI/core_num
i=0
while(( $i<$(($core_num-1)) ))
do
    cp ./SDM_history.ncl ./SDM_history_$(($i*$p)).ncl
    sed -i "s/0,MPI/$(($i*$p)),$((($i+1)*$p))/" ./SDM_history_$(($i*$p)).ncl
    ncl SDM_history_$(($i*$p)).ncl &
    let "i++"
done
cp ./SDM_history.ncl ./SDM_history_$(($i*$p)).ncl
sed -i "s/0,MPI/$(($i*$p)),$((($i+1)*$p))/" ./SDM_history_$(($i*$p)).ncl
ncl SDM_history_$(($i*$p)).ncl
wait
