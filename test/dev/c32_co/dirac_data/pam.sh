#!/bin/sh

## module load dirac/19.0
## Replacement of 'module load dirac/19.0' is given below sentences.
#DIRAC=/home/noda/DIRAC-21.1-Source
#PAM=$DIRAC/bin/pam-dirac
module load DIRAC/22.0
# openmpi (i8)
#export PATH="$DIRAC/openmpi310_i8/bin:$PATH"
#export LD_LIBRARY_PATH="$DIRAC/openmpi310_i8/lib:$LD_LIBRARY_PATH"
export DIRAC_TMPDIR=$HOME/dirac_scr
export OMPI_MCA_hwloc_base_binding_policy=nonei


NPROCS=4

# PAM=pam-dirac
INPFILE=CO.inp
MOLFILE=CO.mol
# parallel
pam --mpi=$NPROCS --get="MRCONEE MDCIN*" '--keep_scratch' --mol=${MOLFILE} --inp=${INPFILE} &> CO.log
# serial
# $PAM --get="MRCONEE MDCINT" '--keep_scratch' --mol=${MOLFILE} --inp=${INPFILE} >& H2.log
#$PAM --mol=${MOLFILE} --inp=${INPFILE} --noarch
