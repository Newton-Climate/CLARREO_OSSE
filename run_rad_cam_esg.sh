#!/bin/sh

rm *.o radiation
rm ./mod4_obj/*.o
make -f Makefile.intel clean

./mkSrcfiles
./mkDepends ./ Srcfiles >Depends

make -f Makefile_modtran4.intel
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../

date >blah.txt
make -f Makefile.intel
rm out.nc *.tp5 *.tp6
nice ./radiation /global/scratch/drfeldma/esg/b30.042a.cam2.h0.2050-10.nc /global/scratch/drfeldma/esg/b30.042a.cam2.h0.2050-10.nc out.nc < /global/scratch/drfeldma/qsub/settings_forcing.inp

