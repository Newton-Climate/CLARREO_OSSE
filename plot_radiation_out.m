

nc = netcdf('./OUTPUT/exp4.nc','read');

wvl = nc{'WAVELENGTH_HRES'}(:);
rad = nc{'RADIANCE_HRES'}(:);

figure(1)
%clf

plot(wvl,rad)

close(nc);
clear nc
