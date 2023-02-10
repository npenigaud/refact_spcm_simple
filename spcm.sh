#!/bin/bash
#SBATCH -N1
#SBATCH --time 00:05:00
#SBATCH --exclusive
#SBATCH --export="NONE"

set -x

cd /home/gmap/mrpm/marguina/pack/48t3_sidyn-spcm.05.IMPIIFC2018.x/spcm_simple


\rm *.grb

#xport MPIAUTOCONFIG=mpiauto.DDT.conf
#/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./spcm.x --case t0031l015-8mpi --write-grib-1 --write-grib-2 --stat
/opt/softs/mpiauto/mpiauto --verbose -np 8 --wrap --wrap-stdeo -- ./spcm.x --case t0107l070-8mpi --stat



