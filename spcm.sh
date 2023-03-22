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

#xport SLURM_EXPORT_ENV=ALL
#xport MPIAUTOCONFIG=mpiauto.PGI.conf
#xport MPIAUTOCONFIG=mpiauto.DDT.conf

#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015c2.4-008mpi --stat-gp
 /opt/softs/mpiauto/mpiauto --verbose -np 1 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015c1.0-001mpi+cor --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015c1.0-008mpi+cor --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015c1.0-008mpi+cor --write-grib-1 --write-grib-2 --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015c1.0-008mpi --write-grib-1 --write-grib-2 --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0107l070c1.0-008mpi --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 1 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015c1.0-001mpi --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 8 -openmp 1 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0149l105c1.0-008mpi --write-grib-1 --write-grib-2 --stat-gp --stat-sp
