#!/bin/sh


make -f Makefile_modtran4.intel
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../

date >blah.txt
make -f Makefile.intel
rm out.nc *.tp5 *.tp6
nice ./radiation /home1/dfeldman/esg/b30.042a.cam2.h0.2000-01.nc /home1/dfeldman/esg/b30.042a.cam2.h0.2000-01.nc out.nc < /global/scratch/drfeldma/qsub/settings_forcing.inp

