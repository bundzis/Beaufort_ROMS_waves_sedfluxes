############## Find Ubar ##################
# The purpose of this script is to use the output 
# from HYCOM2ROMS_u_currents_bryclm.py to find
# ubar and create the ubar clm and bry netcdf files
# using xesmf and the wonderful barotropic.f90 from
# model2roms Github repo


##########################################


# Load in the packages
import numpy as np
import xarray as xr
import xesmf as xe
import ESMF
from datetime import datetime, timedelta
from cftime import num2date, date2num
from netCDF4 import Dataset


# Load in the data
# interpolated u currents 
u_currents = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/Attempt001/u_currents_clm_001.nc') #UPDATE PATH

# vertical grid stuff
grid_vertical = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/ROMS_grid_depth_hpluszeta_2020_003.nc') #, drop_variables='z_w')  #UPDATE PATH 

# Need z_w to be positive so multiply by -1
z_w_pos = grid_vertical.z_w.values * (-1)
# check that this worked
print('z_w: ', grid_vertical.z_w[0,:,200,200].values, flush=True)
print('z_w_pos: ', z_w_pos[0,:,200,200], flush=True)

# Read in the dimensions
#time_len = len(grid_vertical.time)
time_len = len(u_currents.v3d_time)
eta_rho_len = len(grid_vertical.eta_rho) # 206
xi_rho_len = len(grid_vertical.xi_rho) # 608
eta_u_len = len(grid_vertical.eta_u) # 206
xi_u_len = len(grid_vertical.xi_u) # 607
eta_v_len = len(grid_vertical.eta_v) # 205
xi_v_len = len(grid_vertical.xi_v) # 608

# Define other dimension lengths
# eta rho
Mp = len(grid_vertical.eta_rho)

# xi rho
Lp = len(grid_vertical.xi_rho)

# for u/v points 
Lm, Mm = (Lp-2), (Mp-2) #number/dimension of cells
L,  M  = Lm+1, Mm+1 #number/dimension of psi points

# srho
N = len(grid_vertical.s_rho)

# time
# this is done in u_currents so copy from there
time_tmp_len = np.copy(time_len)

# copy the actual time values - seconds since 2000-01-01
#time_tmp = grid_vertical.time.values
time_tmp = u_currents.v3d_time.values


# Set up the netcdf for the climatology 
# ------------------------------- Create the netCDF files ---------------------------

#name of file I am writing to
vert_ubar_currents_clm = '/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/Attempt001/ubar_currents_clm_001.nc' #UPDATE PATH

#create file to write to
nc1 = Dataset(vert_ubar_currents_clm, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006.nc',
                      type = 'ROMS grid vertically interpolated HYCOM ubar currents climatology',
                      history = 'Created by Brianna Undzis',
                      Conventions = 'CF',
                      Institution = 'University of Colorado Boulder',
                      date = str(datetime.today()))

#create dictionary for model
d = {}
d = global_defaults

for att, value in d.items():
    setattr(nc1, att, value)
    
# Create dimensions
nc1.createDimension('xi_rho',  Lp)   # RHO
nc1.createDimension('eta_rho', Mp)
nc1.createDimension('s_rho', N)
nc1.createDimension('s_w',   (N+1))
nc1.createDimension('xi_u',    L)    # u
nc1.createDimension('eta_u',   Mp)
nc1.createDimension('v2d_time', None)
nc1.createDimension('one',     1)

# Create variables # this took several minutes (~9)
# --------------------
# Coordinate Variables
# --------------------
# xi rho
xi_rho = nc1.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
xi_rho.long_name = 'xi coordinate of RHO-points'
xi_rho.standard_name = 'projection_xi_coordinate'
xi_rho.units = 'meter'
xi_rho_tmp = np.arange(0, Lp)
xi_rho[:] = xi_rho_tmp[:]

# eta rho
eta_rho = nc1.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
eta_rho.long_name = 'eta coordinate of RHO-points'
eta_rho.standard_name = 'projection_eta_coordinate'
eta_rho.units = 'meter'
eta_rho_tmp = np.arange(0, Mp)
eta_rho[:] = eta_rho_tmp[:]

# xi u
xi_u = nc1.createVariable('xi_u', 'd', ('xi_u',), zlib = True)
xi_u.long_name = 'xi coordinate of U-points'
xi_u.standard_name = 'projection_xi_coordinate'
xi_u.units = 'meter'
xi_u[:]=np.arange(0,L)

# eta u
eta_u = nc1.createVariable('eta_u', 'd', ('eta_u',), zlib = True)
eta_u.long_name = 'eta coordinate of U-points'
eta_u.standard_name = 'projection_eta_coordinate'
eta_u.units = 'meter'
eta_u[:]=np.arange(0,Mp)

# s rho
s_rho = nc1.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)
s_rho[:] = s_rho_tmp[:]

# v2d_time (in seconds)
v2d_time_g = nc1.createVariable('v2d_time', None, ('v2d_time'), zlib=True)
v2d_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
v2d_time_g.units = 'second'
v2d_time_g.field = 'time, scalar, series'
v2d_time_g[:] = time_tmp[:]

# --------------------
# Boundary Currents
# --------------------

# ubar current
ubar_interp_g = nc1.createVariable('ubar', 'f8', ('v2d_time', 'eta_u', 'xi_u'), zlib=True)
ubar_interp_g.long_name = 'depth-averaged water u momentum'
ubar_interp_g.units = 'meter per second'

# Set up the netcdf for the boundary 
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
vert_ubar_currents_bry = '/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/Attempt001/ubar_currents_bry_001.nc' #UPDATE PATH

#create file to write to
nc2 = Dataset(vert_ubar_currents_bry, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006.nc',
                      type = 'ROMS grid vertically interpolated HYCOM ubar currents boundary',
                      history = 'Created by Brianna Undzis',
                      Conventions = 'CF',
                      Institution = 'University of Colorado Boulder',
                      date = str(datetime.today()))

#create dictionary for model
d = {}
d = global_defaults

for att, value in d.items():
    setattr(nc2, att, value)
    
# Create dimensions
nc2.createDimension('xi_rho',  Lp)   # RHO
nc2.createDimension('eta_rho', Mp)
nc2.createDimension('s_rho', N)
nc2.createDimension('s_w',   (N+1))
nc2.createDimension('xi_u',    L)    # u
nc2.createDimension('eta_u',   Mp)
nc2.createDimension('v2d_time', None)
nc2.createDimension('one',     1)

# Create variables
# --------------------
# Coordinate Variables
# --------------------
# xi rho
xi_rho = nc2.createVariable('xi_rho', 'd', ('xi_rho',), zlib=True)
xi_rho.long_name = 'xi coordinate of RHO-points'
xi_rho.standard_name = 'projection_xi_coordinate'
xi_rho.units = 'meter'
xi_rho_tmp = np.arange(0, Lp)
xi_rho[:] = xi_rho_tmp[:]

# eta rho
eta_rho = nc2.createVariable('eta_rho', 'd', ('eta_rho',), zlib=True)
eta_rho.long_name = 'eta coordinate of RHO-points'
eta_rho.standard_name = 'projection_eta_coordinate'
eta_rho.units = 'meter'
eta_rho_tmp = np.arange(0, Mp)
eta_rho[:] = eta_rho_tmp[:]

# xi u
xi_u = nc2.createVariable('xi_u', 'd', ('xi_u',), zlib = True)
xi_u.long_name = 'xi coordinate of U-points'
xi_u.standard_name = 'projection_xi_coordinate'
xi_u.units = 'meter'
xi_u[:]=np.arange(0,L)

# eta u
eta_u = nc2.createVariable('eta_u', 'd', ('eta_u',), zlib = True)
eta_u.long_name = 'eta coordinate of U-points'
eta_u.standard_name = 'projection_eta_coordinate'
eta_u.units = 'meter'
eta_u[:]=np.arange(0,Mp)

# s rho
s_rho = nc2.createVariable('s_rho', 'd', ('s_rho',), zlib=True)
s_rho.long_name = 's coordinate of RHO-points'
s_rho.standard_name = 'projection_s_coordinate'
s_rho.units = 'meter'
s_rho_tmp = np.arange(0, N)

# v2d_time (in seconds)
v2d_time_g = nc2.createVariable('v2d_time', None, ('v2d_time'), zlib=True)
v2d_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
v2d_time_g.units = 'second'
v2d_time_g.field = 'ubar_time, scalar, series'
v2d_time_g[:] = time_tmp[:]

# --------------------
# Boundary Currents
# --------------------

# ubar current
# ubar_west
ubar_west_g = nc2.createVariable('ubar_west', 'f8', ('v2d_time', 'eta_u'), zlib=True)
ubar_west_g.long_name = 'depth-averaged water u momentum western boundary'
ubar_west_g.units = 'meter per second'

# ubar_north
ubar_north_g = nc2.createVariable('ubar_north', 'f8', ('v2d_time', 'xi_u'), zlib=True)
ubar_north_g.long_name = 'depth-averaged water u momentum northern boundary'
ubar_north_g.units = 'meter per second'

# ubar_east
ubar_east_g = nc2.createVariable('ubar_east', 'f8', ('v2d_time', 'eta_u'), zlib=True)
ubar_east_g.long_name = 'depth-averaged water u momentum eastern boundary'
ubar_east_g.units = 'meter per second'

# --------------------------------- End Netcdf set up -----------------------------------


# We need to interpolate z_w from rho points onto u points 
# Use the weights we already have from  HYCOM2ROMS_u_currents_bryclm.py
# import the weights 
rho2u_weights = '/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Bryclm_conds/regrid_rho2u_weights.nc'

# Make the rho grid 
ds_in_rho = grid_vertical.copy()
ds_in_rho['lat'] = (('eta_rho', 'xi_rho'), ds_in_rho.lat_rho.values)
ds_in_rho['lon'] = (('eta_rho', 'xi_rho'), ds_in_rho.lon_rho.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (ROMS rho)
ds_in_rho['mask'] = (('eta_rho', 'xi_rho'), ds_in_rho.mask_rho.values)

# Make the u grid 
# Output grid (ROMS, but keeps HYCOM vertical levels)
#ds_out_u = grid_vertical
ds_out_u = grid_vertical.copy()
ds_out_u['lat'] = (('eta_u', 'xi_u'), ds_out_u.lat_u.values)
ds_out_u['lon'] = (('eta_u', 'xi_u'), ds_out_u.lon_u.values)

# Add masks
ds_out_u['mask'] = (('eta_u', 'xi_u'), ds_out_u.mask_u.values)

# Make the rho2u regridder
regridder_rho2u = xe.Regridder(ds_in_rho, ds_out_u, "bilinear", extrap_method='nearest_s2d', weights=rho2u_weights)

# Regrid z_w_pos 
dr_rho2u_zwpos = np.copy(z_w_pos)
dr_out_rho2u_zwpos = regridder_rho2u(dr_rho2u_zwpos)
#dr_out_rho2u_zwpos

# Now all our inputs are ready
# Compile and import model2roms barotropic.f90
# Use barotropic.f90 to calculate ubar
from numpy import f2py
with open('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Bryclm_conds/barotropic.f90') as sourcefile1:
    sourcecode1 = sourcefile1.read()
f2py.compile(sourcecode1, modulename='barotropic', extension='.f90')
import barotropic

# Call the function to calculate ubar
# Make an array to fill with the new values 
ubar_interp1 = np.empty((time_len, eta_u_len, xi_u_len))

# Loop through time
for t in range(time_len):
    # Print the time to see where we are
    print('t: ', t)
    
    # Call the function 
    ubar_interp1[t,:,:] = barotropic.velocity.ubar(u_currents.u[t,:,:,:].values, ubar_interp1[t,:,:], 
                                               z_w_pos[t,:,:,:], dr_out_rho2u_zwpos[t,:,:,:])
    
    # Save this to the netcdf
    # climatology 
    ubar_interp_g[t,:,:] = ubar_interp1[t,:,:]
    
    # boundary
    ubar_west_g[t,:] = ubar_interp1[t,:,0]
    ubar_north_g[t,:] = ubar_interp1[t,-1,:]
    ubar_east_g[t,:] = ubar_interp1[t,:,-1]
    
    # force save to netcdf
    nc1.sync()
    nc2.sync()

#print(ubar_interp1[0,200,200])

# Close the netcdfs
nc1.close()
nc2.close()

# Should probably manually check this calculation before writing to netcdf and assuming it is correct
# print the depths 
zws = dr_out_rho2u_zwpos[100,:,180,480]
print('zws: ', zws, flush=True)

# print the distances between the depths (cell thicknesses)
diffs = np.diff(zws)
print('\ndiffs: ', diffs, flush=True)

# Print the currents at these depths 
curr = u_currents.u[100,:,180,480].values
print('\ncurrents: ', curr, flush=True)

ubrr = np.sum(diffs*curr)/np.sum(diffs)
print('\nubrr: ', ubrr, flush=True)

print('\ncalculated ubrr: ', ubar_interp1[100,180,480], flush=True)

# It looks like it is working as expected! Woohoo!

