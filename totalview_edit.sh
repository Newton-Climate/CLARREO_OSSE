#!/bin/sh

date >blah.txt
rm *.tp5 *.nc *.tp6
/global/scratch/kmuriki/TotalView/workbench/toolworks/totalview.8.6.2-1/bin/totalview radiation -a /global/scratch/drfeldma/esg/b30.042a.cam2.h0.2050-10.nc /global/scratch/drfeldma/esg/b30.042a.cam2.h0.2050-10.nc out.nc < /global/scratch/drfeldma/qsub/settings_forcing.inp &
date >>blah.txt

