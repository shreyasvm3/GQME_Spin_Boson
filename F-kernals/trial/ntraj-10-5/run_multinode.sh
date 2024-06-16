#!/bin/bash

#SBATCH -J SB_G2_GQME 
#SBATCH -t 168:00:00
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -p astra
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH --exclusive
#SBATCH -w c0055
# Enter the working directory
cd ${SLURM_SUBMIT_DIR}

echo "starting $SLURM_JOB_ID at `date` on `hostname`"

echo "tmp directory is ${TMPDIR}"

echo "$USER"

MYTMP=/tmp/${USER}/${SLURM_JOB_ID}

mpiexec -pernode /usr/bin/mkdir -vp $MYTMP #|| exit $?

#echo "Copying data over... "

mpiexec -pernode scp -v ${SLURM_SUBMIT_DIR}/dyn.x  $MYTMP/. #|| exit $?
mpiexec -pernode scp -v ${SLURM_SUBMIT_DIR}/input  $MYTMP/. #|| exit $?

echo "$(pwd)"

cd $MYTMP

echo "$(pwd)"

module swap gnu8 intel
#run fortran code
time mpiexec -np 12 ./dyn.x

# Copy output files back from the temp directory to working directory
mpiexec -pernode rsync -r $MYTMP/ $SLURM_SUBMIT_DIR/ #|| exit $?

rm -rf $MYTMP

exit 0 


 
