#!/bin/sh
#rm *.o radiation
make -f Makefile_modtran4
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../

date >blah.txt
#cp params_cam_single_column.h params.h
make
nice ./radiation /global/home/users/drfeldma/input_files/cam_output_subset_aerosol.nc /global/home/users/drfeldma/input_files/cam_output_subset_aerosol.nc out_cam_single_column_aerosol.nc < ./INPUT/settings_forcing.inp  >>blah.txt
date >>blah.txt

