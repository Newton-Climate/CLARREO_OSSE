#!/bin/sh
make -f Makefile_modtran5.intel
./makelib5.sh
make -f Makefile.intel clean
make -f Makefile.intel
rm b30.*.nc 
rm blah.txt
rm modtran_spectrum_out.nc
rm modtran_out_b30*.tp5
mpiexec -np 1 ./radiation >blah.txt
