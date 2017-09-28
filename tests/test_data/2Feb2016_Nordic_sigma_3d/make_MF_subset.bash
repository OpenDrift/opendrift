e0=0
e1=20
x0=40
x1=70
ncks -d eta_rho,$e0,$e1 -d eta_psi,$e0,$e1 -d eta_v,$e0,$e1 -d eta_u,$e0,$e1 -d xi_rho,$x0,$x1 -d xi_psi,$x0,$x1 -d xi_v,$x0,$x1 -d xi_u,$x0,$x1 Nordic-4km_SLEVELS_avg_00_subset2Feb2016.nc Nordic_subset.nc -O --fl_fmt=netcdf4_classic

ncks -O -d ocean_time,0 Nordic_subset.nc Nordic_subset_day1.nc
ncks -O -d ocean_time,1 Nordic_subset.nc Nordic_subset_day2.nc
ncks -O -d ocean_time,2 Nordic_subset.nc Nordic_subset_day3.nc

readerinfo.py Nordic_subset.nc -p x_sea_water_velocity
