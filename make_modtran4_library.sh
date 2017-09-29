#!/bin/sh
rm *.o radiation
rm  Mod4v3r1_F90.exe
make -f Makefile_modtran4
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../
cp params_rtmip.h params.h
make
date
./radiation /global/home/users/drfeldma/input_files/exp4.nc /global/home/users/drfeldma/input_files/exp4.nc ./OUTPUT/exp4.nc <./INPUT/settings_forcing.inp >blah.txt
date
#rm diff_file.nc
#ncdiff ./OUTPUT/exp4.nc ./OUTPUT/standard.nc diff_file.nc
#ncdump diff_file.nc | more
