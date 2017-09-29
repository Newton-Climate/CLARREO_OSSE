#!/bin/sh

#modtran

make -f Makefile_modtran4
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../
cp params_rtmip.h params.h

#exp1
make
./radiation ../../../input_files/exp1.nc ../../../input_files/exp1.nc ./OUTPUT/exp1.nc <./INPUT/settings_forcing.inp >blah.txt

#exp2a
./radiation ../../../input_files/exp2a.nc ../../../input_files/exp2a.nc ./OUTPUT/exp2a.nc <./INPUT/settings_forcing.inp >blah.txt

#exp2b
./radiation ../../../input_files/exp2b.nc ../../../input_files/exp2b.nc ./OUTPUT/exp2b.nc <./INPUT/settings_forcing.inp >blah.txt

#exp3a
./radiation ../../../input_files/exp3a.nc ../../../input_files/exp3a.nc ./OUTPUT/exp3a.nc <./INPUT/settings_forcing.inp >blah.txt

#exp3b
./radiation ../../../input_files/exp3b.nc ../../../input_files/exp3b.nc ./OUTPUT/exp3b.nc <./INPUT/settings_forcing.inp >blah.txt

#exp3c
./radiation ../../../input_files/exp3c.nc ../../../input_files/exp3c.nc ./OUTPUT/exp3c.nc <./INPUT/settings_forcing.inp >blah.txt

#exp3d
./radiation ../../../input_files/exp3d.nc ../../../input_files/exp3d.nc ./OUTPUT/exp3d.nc <./INPUT/settings_forcing.inp >blah.txt

#exp4
./radiation ../../../input_files/exp4.nc ../../../input_files/exp4.nc ./OUTPUT/exp4.nc <./INPUT/settings_forcing.inp >blah.txt
