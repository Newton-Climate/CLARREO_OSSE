
clear

nc = netcdf('../input_files/modis_brdf_2007.nc','read');

brdf = nc{'brdf'}(:,:,:,:,:);
count = nc{'count'}(:,:,:,:,:);

figure(1)
clf

mask = ones(128,256);
for i=1:128
  for j=1:256
    if brdf(1,1,1,i,j)>482 & brdf(1,1,1,i,j)<483
       [i j]
       brdf(1,1,1,i,j)
    end
  end
end


