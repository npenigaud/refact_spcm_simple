#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --exclusive
#SBATCH --export="NONE"

set -x

ulimit -s unlimited
export OMP_STACK_SIZE=4G

cd /home/gmap/mrpm/marguina/pack/48t3_sidyn-spcm.05.IMPIIFC2018.x/spcm_simple


\rm *.grb

#xport MPIAUTOCONFIG=mpiauto.DDT.conf
 /opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015-008mpi --write-grib-1 --write-grib-2 --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0107l070-008mpi --stat-gp
#/opt/softs/mpiauto/mpiauto --verbose -np 1 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0031l015-001mpi --stat-gp


#/opt/softs/mpiauto/mpiauto --verbose -np 8 -openmp 1 --wrap --wrap-stdeo -- ./intel/spcm.x --case t0149l105-008mpi --write-grib-1 --write-grib-2 --stat-gp --stat-sp

