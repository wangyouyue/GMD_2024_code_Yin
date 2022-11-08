#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=ha
#PJM -L node=192
#PJM --mpi proc=7680
#PJM -L elapse=240:00:00
#PJM -g r22775
#PJM -j
#------- Program execution -------#
module load intel impi hdf5 netcdf netcdf-fortran
source ~/.bashrc
# run
mpiexec.hydra -n ${PJM_MPI_PROC} ./scale-rm_init init.conf || exit
mpiexec.hydra -n ${PJM_MPI_PROC} ./scale-rm run.conf || exit
