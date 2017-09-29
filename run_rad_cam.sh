#!/bin/sh
rm *.o radiation
make -f Makefile_modtran4
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../

date >blah.txt
cp params_cam.h params.h
make -f Makefile
nice ./radiation /global/home/users/drfeldma/input_files/cam_output.nc /global/home/users/drfeldma/input_files/cam_output.nc ./OUTPUT/out_clear.nc < ./INPUT/settings_forcing.inp >>blah.txt
date >>blah.txt

