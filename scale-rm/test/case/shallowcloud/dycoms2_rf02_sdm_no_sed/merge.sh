 #!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=sa
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=10:00:00
#PJM -g r22775
#PJM -j
#------- Program execution -------#
module load intel impi hdf5 netcdf netcdf-fortran
source ~/.bashrc
conda activate ncl_stable
# run
ncl merge.ncl
