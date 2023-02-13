
# Intel

FIAT=/home/gmap/mrpm/marguina/pack/48t3_mrg48R1.05.IMPIIFC2018.x/hub/local/install/Fiat
ECTRANS=/home/gmap/mrpm/marguina/pack/48t3_mrg48R1.05.IMPIIFC2018.x/hub/local/install/Ectrans
ECCODES=/home/gmap/mrpm/khatib/opt/i-2018.5.274/eccodes-2.27.0

FFLAGS= \
   -align array64byte,all -traceback -fpic -qopenmp -qopenmp-threadprivate compat \
   -fp-model source -ftz -diag-disable=remark,cpu-dispatch -qopt-report=5 -qopt-report-phase=vec -g -O2 -march=core-avx2 \
   -finline-functions -finline-limit=500 -Winline -qopt-prefetch=4 -fast-transcendentals -fimf-use-svml -no-fma 

#FFLAGS=-g -O0 -qopenmp -qopenmp-threadprivate compat

FC=/home/gmap/mrpm/khatib/public/bin/mimpifc-18.0.5.274 \
   -I$(FIAT)/include \
   -I$(FIAT)/module/parkind_dp \
   -I$(FIAT)/module/fiat \
   -I$(ECTRANS)/include/ectrans \
   -I$(ECTRANS)/module/trans_dp \
   -I$(ECCODES)/include \
   -convert big_endian -assume byterecl \
   $(FFLAGS)

LIBS=-L$(FIAT)/lib64 -lparkind_dp -lfiat -Wl,-rpath,$(FIAT)/lib64 -L$(ECTRANS)/lib64 -ltrans_dp -Wl,-rpath,$(ECTRANS)/lib64 -L$(ECCODES)/lib64 -leccodes_f90 -leccodes -Wl,-rpath,$(ECCODES)/lib64

# PGI

FIAT=/home/gmap/mrpm/marguina/hub/install
ETRANS=/home/gmap/mrpm/marguina/hub/install
ECCODES=/home/gmap/mrpm/marguina/install/PGI2211/eccodes-2.14.0

FFLAGS= -mp -fPIC -Minfo -gopt -Mlarge_arrays -Mlist -traceback -Mnofma -Mbyteswapio -Mbackslash

FC=/home/gmap/mrpm/marguina/install/gmkpack_support/wrapper/PGI2211/mpif90 \
   -I$(FIAT)/include \
   -I$(FIAT)/module/parkind_dp \
   -I$(FIAT)/module/fiat \
   -I$(ECTRANS)/include/ectrans \
   -I$(ECTRANS)/module/trans_dp \
   -I$(ECCODES)/include \
   $(FFLAGS)

LIBS=-L$(FIAT)/lib64 -lparkind_dp -lfiat -Wl,-rpath,$(FIAT)/lib64 -L$(ECTRANS)/lib64 -ltrans_dp -Wl,-rpath,$(ECTRANS)/lib64 -L$(ECCODES)/lib -leccodes_f90 -leccodes -Wl,-rpath,$(ECCODES)/lib \
        -L/opt/softs/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/lib -llapack_lp64 -lblas_lp64 -Wl,-rpath,/opt/softs/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/lib

