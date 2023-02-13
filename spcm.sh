#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --exclusive
#SBATCH --export="NONE"

set -x

ulimit -s unlimited
export OMP_STACK_SIZE=4G

cd $SLURM_SUBMIT_DIR

\rm *.grb

export SLURM_EXPORT_ENV=ALL
export MPIAUTOCONFIG=mpiauto.PGI.conf
#xport MPIAUTOCONFIG=mpiauto.DDT.conf

 /opt/softs/mpiauto/mpiauto --nouse-slurm-mpi --verbose -np 8 --wrap --wrap-stdeo -- ./spcm.x --case t0031l015-008mpi --write-grib-1 --write-grib-2 --stat-gp
#/opt/softs/mpiauto/mpiauto --nouse-slurm-mpi --verbose -np 8 --wrap --wrap-stdeo -- ./spcm.x --case t0107l070-008mpi --stat-gp
#/opt/softs/mpiauto/mpiauto --nouse-slurm-mpi --verbose -np 1 --wrap --wrap-stdeo -- ./spcm.x --case t0031l015-001mpi --stat-gp


#/opt/softs/mpiauto/mpiauto --nouse-slurm-mpi --verbose -np 8 -openmp 1 --wrap --wrap-stdeo -- ./spcm.x --case t0149l105-008mpi --write-grib-1 --write-grib-2 --stat-gp --stat-sp

