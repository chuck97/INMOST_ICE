import numpy as np
from netCDF4 import Dataset

for k in range(2010, 2020):
	nc_f = "./concentration_data/seaice_conc_monthly_nh_f17_"+ str(k) + "03_v03r01.nc"
	nc_fid = Dataset(nc_f, 'r')
	concs = nc_fid.variables['seaice_conc_monthly_cdr'][:]
	size_x = nc_fid.dimensions['xgrid']
	size_y = nc_fid.dimensions['ygrid']
	proj = nc_fid.variables['projection']
	step_in_km = 25.0

	xmin = proj.grid_boundary_left_projected_x
	xmax = proj.grid_boundary_right_projected_x
	ymin = proj.grid_boundary_bottom_projected_y
	ymax = proj.grid_boundary_top_projected_y

	hd =                               \
	'ncols=' + str(size_x.size) +      \
	',nrows=' + str(size_y.size) +     \
	',xmin=' + str(xmin) +             \
	',xmax=' + str(xmax) +             \
	',ymin=' + str(ymin) +             \
	',ymax=' + str(ymax) +             \
	',step=' + str(step_in_km*1000)    \

	np.savetxt('./concentration_data/concentrations_' + str(k) + '_3.txt', concs[0], header = hd)

	nc_f = "./concentration_data/seaice_conc_monthly_nh_f17_"+ str(k) + "04_v03r01.nc"
	nc_fid = Dataset(nc_f, 'r')
	concs = nc_fid.variables['seaice_conc_monthly_cdr'][:]
	size_x = nc_fid.dimensions['xgrid']
	size_y = nc_fid.dimensions['ygrid']
	proj = nc_fid.variables['projection']
	step_in_km = 25.0

	xmin = proj.grid_boundary_left_projected_x
	xmax = proj.grid_boundary_right_projected_x
	ymin = proj.grid_boundary_bottom_projected_y
	ymax = proj.grid_boundary_top_projected_y

	hd =                               \
	'ncols=' + str(size_x.size) +      \
	',nrows=' + str(size_y.size) +     \
	',xmin=' + str(xmin) +             \
	',xmax=' + str(xmax) +             \
	',ymin=' + str(ymin) +             \
	',ymax=' + str(ymax) +             \
	',step=' + str(step_in_km*1000)    \

	np.savetxt('./concentration_data/concentrations_' + str(k) + '_4.txt', concs[0], header = hd)
