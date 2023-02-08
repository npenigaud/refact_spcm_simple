
FIAT=/home/gmap/mrpm/marguina/pack/48t3_mrg48R1.05.IMPIIFC2018.x/hub/local/install/Fiat

FC=/home/gmap/mrpm/khatib/public/bin/mimpifc-18.0.5.274 \
   -I$(FIAT)/include \
   -I$(FIAT)/module/parkind_dp \
   -I$(FIAT)/module/fiat \
   -convert big_endian -assume byterecl -align array64byte,all -traceback -fpic -qopenmp -qopenmp-threadprivate compat \
   -fp-model source -ftz -diag-disable=remark,cpu-dispatch -qopt-report=5 -qopt-report-phase=vec -g -O2 -march=core-avx2 \
   -finline-functions -finline-limit=500 -Winline -qopt-prefetch=4 -fast-transcendentals -fimf-use-svml -no-fma 
LIBS=-L$(FIAT)/lib64 -lparkind_dp -lfiat -Wl,-rpath,$(FIAT)/lib64
