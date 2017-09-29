#!/bin/sh
rm radiation
make -f Makefile_modtran4
cd ./mod4_obj 
ar rc libmodtran4.a *.o
mv libmodtran4.a ../modtran_lib
cd ../

make -f Makefile

date >blah.txt
cp params_cam_single_column.h params.h
make
nice ./radiation ../../../input_files/cam_output_single_column_hlat.nc ../../../input_files/cam_output_single_column_hlat.nc out_clear_single_column_hlat.nc < ./INPUT/settings_forcing.inp >>blah.txt
date >>blah.txt

