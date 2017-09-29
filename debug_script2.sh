#!/bin/sh
make -f Makefile.intel
rm b30.*.nc 
rm blah.txt
rm modtran_spectrum_out.nc
rm modtran_out_b30*.tp5 *.flx
mpiexec -np 1 ./radiation >blah.txt
