#!/bin/bash
#PBS -q S
#PBS -l select=2:ncpus=20:mpiprocs=20
#PBS -l walltime=24:00:00
#PBS -N scale-sdm

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}
module load intel/17.0.0 mpt hdf5/1.8.12 netcdf/4.4.1

# run
mpiexec_mpt dplace -s1 ./scale-rm_init init.conf || exit
mpiexec_mpt dplace -s1 ./scale-rm  run.conf || exit
