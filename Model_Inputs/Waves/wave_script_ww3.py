"""
Created on Wed Sep 28 14:38:49 2022

@author: brun1463
"""

################ Wave Forcing - ERA5 #####################
# The purpose of this script is to look 
# at the ERA5 wave data for 2019
# to see if it would be good to use for the wave 
# forcing. If it is, this script will then make 
# the wave forcing file for ROMS.
#
#
###################################################


# Load in the packages 
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd
import ESMF
import math
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from cftime import num2date, date2num

# Load in the HYCOM wave data
# OG
#wave_data = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_scratch/External_data/ERA5_2019_wavedata_hsperdir01.nc')  # UPDATE PATH
wave_data = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Waves/ww3_global_WaveWatch_III_Global_Wave_Model_best.nc')  # UPDATE PATH
# TEMP
#wave_data = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_scratch/External_data/ERA5_2020_wavedata_hsperdir01.nc')

# Look at the data
#print(wave_data)

# Cut out just the time period we care about for the model - the open water season
# waves 
#era5_waves = wave_data.sel(time=slice('2019-07-01','2019-11-01 21:00:00'))
# OG
ww3_waves = wave_data.sel(time=slice('2020-07-01 09:00:00', '2020-11-03 09:00:00')) # this is 2019-07-01 01:00:00 - 2019-11-03 01:00:00 in AKDT
# TEMP
#era5_waves = wave_data.sel(time=slice('2020-07-01 09:00:00', '2020-11-02 09:00:00')) # this is 2020-07-01 01:00:00 - 2020-11-02 01:00:00 in AKDT

# Make a datetime array of the corresponding datetime values in AKDT to use for the netcdf
# OG
#time_akdt = np.arange(datetime(2020,7,1,hour=1,minute=0,second=0), datetime(2020,11,2,hour=4,minute=0, second=0),timedelta(hours=3)) # OG
time_akdt = np.arange(datetime(2020,7,1,hour=1,minute=0,second=0), datetime(2020,11,3,hour=2,minute=0, second=0),timedelta(hours=1)) # try hourly
# TEMP
#time_akdt = np.arange(datetime(2020,7,1,hour=1,minute=0,second=0), datetime(2020,11,2,hour=4,minute=0, second=0),timedelta(hours=3))
time_akdt_dt = pd.to_datetime(time_akdt)
#print(time_akdt_dt[0:7]) # local time
#print(time_akdt_dt[-9:-1]) # local time
#print(time_akdt_dt[-1]) # local time
#print(len(time_akdt_dt))
#print('era5 time length: ', len(era5_waves.time.values)) # GMT
#print('start time: ', era5_waves.time[0].values) # GMT
#print('end time: ', era5_waves.time[-1].values) # GMT

# Check the lat/lon convention of the era5 data
#print('era5 lat: ', era5_waves.latitude.values)
#print('era5 long: ', era5_waves.longitude.values)

# Look at the era5 data
# Plot significant wave height 
#fig1, ax1 = plt.subplots()
#cs1 = ax1.contourf(era5_waves.longitude.values, era5_waves.latitude.values, era5_waves.swh[0,:,:].values)
#cbar1 = plt.colorbar(cs1, label='Significant Wave Height (m?)')
#ax1.set_title('Significant Wahe Height (swh) (m?)')
#ax1.set_xlabel('Longitude')
#ax1.set_ylabel('Latitude')

# Load in the ROMS grid 
# List variables to drop 
#drop_vars = ['z_w', 'z_rho']
grid = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Include/KakAKgrd_shelf_big010_smooth006.nc')  # UPDATE PATH
#'/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/ROMS_grid_depth_hpluszeta_011.nc', drop_variables=drop_vars)

# Check the lat.lon convention of the grid 
#print('grid lat: ', grid.lat_rho.values)
#print('grid long: ', grid.lon_rho.values)
# Looks like both grids are in the same convention so we can leave it

# Pull out the angle to rotate the currents to match the grid's u,v
phi = grid.angle[0,0].values # radians 

# Read in the dimensions
#time_len = len(grid_vertical.time)
eta_rho_len = len(grid.eta_rho) # 206
xi_rho_len = len(grid.xi_rho) # 608
eta_u_len = len(grid.eta_u) # 206
xi_u_len = len(grid.xi_u) # 607
eta_v_len = len(grid.eta_v) # 205
xi_v_len = len(grid.xi_v) # 608

# Define other dimension lengths
# eta rho
Mp = len(grid.eta_rho)

# xi rho
Lp = len(grid.xi_rho)

# for u/v points 
Lm, Mm = (Lp-2), (Mp-2) #number/dimension of cells
L,  M  = Lm+1, Mm+1 #number/dimension of psi points

# latitude
lat_len = len(grid.lat_rho)

# longitude
lon_len = len(grid.lon_rho)

# time
# Use the time from the era5 data
#time_len = len(era5_waves.time) # OG
#datetime1 = pd.to_datetime(era5_waves.time.values) # OG
time_len = len(time_akdt_dt)
datetime1 = time_akdt_dt

# Convert all time to seconds since 2000-01-01 (really 1999-12-31 so it starts at beginning
# of year=hour 0)
time_tmp = ((datetime1[:] - datetime(1999,12,31)).total_seconds() - 86400)
# Make sure this worked
#print(time_tmp[0:5])
#print(time_tmp[-5:-1])

# name the variables to fill the netcdf 
# latitude
lat_tmp = grid.lat_rho.values

# longitude
lon_tmp = grid.lon_rho.values

# ww3 number of lats
ww3_lat_len = len(ww3_waves.lat.values)
#print('era5_lat_len', era5_lat_len)

# era5 number of lons
ww3_lon_len = len(ww3_waves.lon.values)
#print('era5_lon_len', era5_lon_len)


# Set up the netcdf for the forcing
# ------------------------------- Create the netCDF file ---------------------------

#name of file I am writing to
# OG
wave_frc = '/projects/brun1463/ROMS/Kakak3_Alpine_2020/Include/wave_forcing_file_kaktovik_shelf_ww3_2020_double_waves_data002.nc'  # UPDATE PATH
# TEMP
#wave_frc = '/projects/brun1463/ROMS/Kakak3_Alpine/Include/wave_forcing_file_kaktovik_shelf_era5_2020_data001.nc'  # UPDATE PATH
#create file to write to
nc1 = Dataset(wave_frc, 'w', format='NETCDF4')

#Global attributes
global_defaults = dict(gridname = 'KakAKgrd_shelf_big010_smooth006.nc',
                      type = 'ROMS grid wave forcing file',
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
nc1.createDimension('wave_time', None)
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

# wave_time (in seconds)
wave_time_g = nc1.createVariable('wave_time', None, ('wave_time'), zlib=True)
wave_time_g.long_name = 'seconds since 2000-01-01 00:00:00' #with initialization of 2000-01-01 00:00:00
wave_time_g.units = 'second'
wave_time_g.field = 'time, scalar, series'
wave_time_g[:] = time_tmp[:]
    
# =============================================================================
# # time (HYCOM version)
# time_g2 = nc1.createVariable('time_HYCOM', None, ('time_HYCOM'), zlib=True)
# time_g2.long_name = '3-hour time steps' 
# time_g2.units = 'datetime'
# time_g2.field = 'time_HYCOM, scalar, series'
# time_g2[:] = time_tmp2[:]
# =============================================================================

# --------------------
# Waves
# --------------------
# ************* copying and editing HYCOM2ROMS_salt_bryclm.py *****************8
# Significant wave height (Hwave)
swh_interp_g = nc1.createVariable('Hwave', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
swh_interp_g.long_name = 'wind-induced significant wave height'
swh_interp_g.units = 'meter' 

# REAL NAME OF WAVE DATA FROM ERA5
# Mean wind-induced wave direction (Dwave)
dir_interp_g = nc1.createVariable('Dwave', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
dir_interp_g.long_name = 'wind-induced wave direction - mean'
dir_interp_g.units = 'degrees' # **CHECK THE CONVENTION

# *** FAKE NAME OF WAVE DATA FROM ERA5 ***
# Mean wind-induced wave direction (Dwave)
# Renaming as peak to test something in ROMS
pdir_interp_g = nc1.createVariable('Dwavep', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
pdir_interp_g.long_name = 'wind-induced wave direction - peak'
pdir_interp_g.units = 'degrees' # **CHECK THE CONVENTION - I think this is right
# ****************************************

# Peak wind-induced surface wave period (Pwave_top)
pwavet_interp_g = nc1.createVariable('Pwave_top', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
pwavet_interp_g.long_name = 'wind-induced peak surface wave period'
pwavet_interp_g.units = 'second' 

# Bottom orbital velocity 
ubr_g = nc1.createVariable('Uwave_rms', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
ubr_g.long_name = 'wind-induced bottom orbital velocity'
ubr_g.units = 'meter second-1'

# Bottom wave period 
pwaveb_g = nc1.createVariable('Pwave_bot', 'f8', ('wave_time', 'eta_rho', 'xi_rho'), zlib=True)
pwaveb_g.long_name = 'wind-induced bottom wave period'
pwaveb_g.units = 'second'

# ------------------------------- End netCDF file setup ---------------------------

# Set the input and output grids, and sepcify the lat/lon
# Since we are looking at waves, we will use lon_rho and lat_rho as the primary lat/lon for the grid 
# Input grid (era5)
ds_in_ww3 = ww3_waves.copy() 
#ds_in_era5['lon_360'] = ds_in_hycom.lon.values
ds_in_ww3['lon'] = ds_in_ww3.lon.values
ds_in_ww3['lat'] = ds_in_ww3.lat.values

# Output grid (ROMS rho grid)
#ds_out_rho = grid_vertical
ds_out_rho = grid.copy()
ds_out_rho['lat'] = (('eta_rho', 'xi_rho'), ds_out_rho.lat_rho.values)
ds_out_rho['lon'] = (('eta_rho', 'xi_rho'), ds_out_rho.lon_rho.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (era5)
# this is only a surface mask - which is what we want 
ds_in_ww3_mask = xr.where(~np.isnan(ds_in_ww3['Thgt'][0,0,:,:].values), 1, 0) 
print('got through nan mask', flush=True)
ds_in_ww3['mask'] = (('lat', 'lon'), ds_in_ww3_mask)

# Output grid (ROMS rho grid)
ds_out_rho['mask'] = (('eta_rho', 'xi_rho'), ds_out_rho.mask_rho.values)

# Regrid from era5 grid to rho grid with the masks included and extrapolation used 
regridder_ww32rho = xe.Regridder(ds_in_ww3, ds_out_rho, method="bilinear", extrap_method='inverse_dist') #extrap_method="nearest_s2d"
regridder_ww32rho

# Save the weights - only need to do this once
fn_ww32rho = regridder_ww32rho.to_netcdf('regrid_ww32rho_weights.nc')
#print(fn_era52rho)

# Now use the regridder/weights to regrid the significant wave height
dr_ww32rho_Thgt = ww3_waves['Thgt'][:,0,:,:].copy() # significant hieght of combined wind waves and swell
dr_out_ww32rho_Thgt = regridder_ww32rho(dr_ww32rho_Thgt) 
dr_out_ww32rho_Thgt

# Now use the regridder/weights to regrid the mean wave direction
dr_ww32rho_Tdir = ww3_waves['Tdir'][:,0,:,:].copy() # peak wave direction 
dr_out_ww32rho_Tdir = regridder_ww32rho(dr_ww32rho_Tdir) 
dr_out_ww32rho_Tdir

# Now use the regridder/weights to regrid the peak wave period
dr_ww32rho_Tper = ww3_waves['Tper'][:,0,:,:].copy() # peak wave period 
dr_out_ww32rho_Tper = regridder_ww32rho(dr_ww32rho_Tper) 
dr_out_ww32rho_Tper 

# Use fill.f90 to fill the nans in the array
# Import fill.f90 from model2roms to see how to 
# use this/if it can be used to get rid of nans 
from numpy import f2py
with open('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/fill.f90') as sourcefile2:
    sourcecode2 = sourcefile2.read()
f2py.compile(sourcecode2, modulename='fill', extension='.f90')
import fill

# Define a function to call to do the filling, taken from model2roms
def laplacefilter(field, threshold, toxi, toeta):
    undef = 2.0e+35 
    tx = 0.9 * undef
    critx = 0.01
    cor = 1.6
    mxs = 10

    field = np.where(abs(field) > threshold, undef, field)

    field = fill.extrapolate.fill(int(1), int(toxi),
                                int(1), int(toeta),
                                float(tx), float(critx), float(cor), float(mxs),
                                np.asarray(field, order='F'),
                                int(toxi),
                                int(toeta))
    return field

# Loop through time to fill in the nans
# Make some variables first
toxi = xi_rho_len
toeta = eta_rho_len

# Make an array to hold the new data without nans
#print('got here 1')
#dr_out_era52rho_zeta_nonan = np.empty((time_len, eta_rho_len, xi_rho_len))
# Significant wave height
dr_out_ww32rho_Thgt_nonan = np.empty((time_len, eta_rho_len, xi_rho_len))

# Mean wave direction
dr_out_ww32rho_Tdir_nonan = np.empty((time_len, eta_rho_len, xi_rho_len)) 

# Peak wave period 
dr_out_ww32rho_Tper_nonan = np.empty((time_len, eta_rho_len, xi_rho_len))
#print('got here 2')

# Make a copy of the OG array to work with
#print('got here 3')
# Significant wave height
dr_out_ww32rho_Thgt_cp1 = dr_out_ww32rho_Thgt.copy() 

# Mean wave direction
dr_out_ww32rho_Tdir_cp1 = dr_out_ww32rho_Tdir.copy() 

# Peak wave period 
dr_out_ww32rho_Tper_cp1 = dr_out_ww32rho_Tper.copy() 
#print('got here 4')

# Loop through depth to replace all the nans with real values 
# Loop through time
for t in range(time_len):
    # Print the time 
    print('t: ', t, flush=True)

    # Pull out the horizontal 'field' for that time
    # Significant wave height
    #field = dr_out_hycom2rho_zeta_cp1[t,:,:]
    field1 = dr_out_ww32rho_Thgt_cp1[t,:,:]
    
    # Mean wave direction 
    field2 = dr_out_ww32rho_Tdir_cp1[t,:,:]
    
    # Peak wave period
    field3 = dr_out_ww32rho_Tper_cp1[t,:,:]


    # Use the Laplace Filter to get rid of nans
    #field = laplacefilter(field, 1000, toxi, toeta)
    field1 = laplacefilter(field1, 1000, toxi, toeta) # Thgt
    field2 = laplacefilter(field2, 1000, toxi, toeta) # Tdir
    field3 = laplacefilter(field3, 1000, toxi, toeta) # Tper
    

    # Multiply by the rho mask 
    #field = field * ds_out_rho.mask.values
    field1 = field1 * ds_out_rho.mask.values # Thgt
    field2 = field2 * ds_out_rho.mask.values # Tdir
    field3 = field3 * ds_out_rho.mask.values # Tper

    # Check to see if there are any nans
    #print('nans: ', np.where(np.isnan(field)))
    #print('nanmin: ', np.nanmin(field))
    #print('nanmax: ', np.nanmax(field))
    #input('press enter to continue...')

    # Save this field to a new array
    #dr_out_hycom2rho_zeta_nonan[t,:,:] = field
    dr_out_ww32rho_Thgt_nonan[t,:,:] = field1 # Thgt
    dr_out_ww32rho_Tdir_nonan[t,:,:] = field2 # Tdir
    dr_out_ww32rho_Tper_nonan[t,:,:] = field3 # Tper
    

# # TEMP FIX UNTIL THINGS ARE UPDATED/I find better wave data/SWAN is run
# # Try resmapling over time here
# # shww
# # Make a dataset
# dr_out_ww32rho_Thgt_nonan_xr = xr.Dataset(
#     data_vars=dict(Thgt_nonan=(['time','eta_rho', 'xi_rho'], dr_out_ww32rho_Thgt_nonan)
#                    ), 
#     coords=dict(
#         time=('time', time_akdt), 
#         eta_rho=('eta_rho', grid.eta_rho.values), 
#         xi_rho=('xi_rho', grid.xi_rho.values))
#     )

# # Inteprolate over the nans
# #dr_out_era52rho_shww_nonan2 = dr_out_era52rho_shww_nonan_xr.shww_nonan.resample(time='1H').interpolate('linear')
# dr_out_ww32rho_Thgt_nonan2 = dr_out_ww32rho_Thgt_nonan_xr.Thgt_nonan.interpolate_na(dim='time', method='linear', max_gap=None)
# dr_out_ww32rho_Thgt_nonan3 = dr_out_ww32rho_Thgt_nonan2.fillna(0.0) 

# # dir
# # Make a dataset
# dr_out_ww32rho_Tdir_nonan_xr = xr.Dataset(
#     data_vars=dict(Tdir_nonan=(['time','eta_rho', 'xi_rho'], dr_out_ww32rho_Tdir_nonan)
#                    ), 
#     coords=dict(
#         time=('time', time_akdt), 
#         eta_rho=('eta_rho', grid.eta_rho.values), 
#         xi_rho=('xi_rho', grid.xi_rho.values))
#     )

# # Inteprolate over the nans
# dr_out_www32rho_Tdir_nonan2 = dr_out_ww32rho_Tdir_nonan_xr.Tdir_nonan.interpolate_na(dim='time', method='linear', max_gap=None)
# dr_out_ww32rho_Tdir_nonan3 = dr_out_ww32rho_Tdir_nonan2.fillna(0.0)

# # pwavet
# # Make a dataset
# dr_out_ww32rho_Tper_nonan_xr = xr.Dataset(
#     data_vars=dict(Tper_nonan=(['time','eta_rho', 'xi_rho'], dr_out_ww32rho_Tper_nonan)
#                    ), 
#     coords=dict(
#         time=('time', time_akdt), 
#         eta_rho=('eta_rho', grid.eta_rho.values), 
#         xi_rho=('xi_rho', grid.xi_rho.values))
#     )

# # Inteprolate over the nans
# dr_out_ww32rho_Tper_nonan2 = dr_out_ww32rho_Tper_nonan_xr.Tper_nonan.interpolate_na(dim='time', method='linear', max_gap=None)
# dr_out_ww32rho_Tper_nonan3 = dr_out_ww32rho_Tper_nonan2.fillna(0.0)
# # END TEMP 

    
# Save the interpolated/regridded data to the netcdf files
# to the climatology file
#zeta_interp_g[:,:,:] = dr_out_hycom2rho_zeta_nonan[:,:,:]
#swh_interp_g[:,:,:] = dr_out_ww32rho_Thgt_nonan[:,:,:] # Thgt OG
swh_interp_g[:,:,:] = dr_out_ww32rho_Thgt_nonan[:,:,:]*2.0 # Thgt, double waves version
dir_interp_g[:,:,:] = dr_out_ww32rho_Tdir_nonan[:,:,:] # Tdir
pdir_interp_g[:,:,:] = dr_out_ww32rho_Tdir_nonan[:,:,:] # Tdir
pwavet_interp_g[:,:,:] = dr_out_ww32rho_Tper_nonan[:,:,:] # Tper
# TEMP FIX
# swh_interp_g[:,:,:] = dr_out_ww32rho_Thgt_nonan3[:,:,:] # Thgt
# dir_interp_g[:,:,:] = dr_out_ww32rho_Tdir_nonan3[:,:,:] # Tdir
# pdir_interp_g[:,:,:] = dr_out_ww32rho_Tdir_nonan3[:,:,:] # Tdir
# pwavet_interp_g[:,:,:] = dr_out_ww32rho_Tper_nonan3[:,:,:] # Tper


# NEW ATTEMPT AT PWAVE_BOT
# Now that we have the regirdded wave data, use this along with equations
# from Wiberg & Sherwood 2008 Appendices D & E to calculate the bottom
# wave period and bottom orbital velocity
# Load in the function from the python script
from Wiberg_Sherwood_ubspecpar import qkhfs, ubspecpar
#from ubspecpar import ubspecpar
#from qkhfs import qkhfs

# Create empty arrays to hold the values
# Bottom orbital velocity
ubr_tmp = np.empty((time_len, eta_rho_len, xi_rho_len)) 
ubr_tmp_mask = np.empty((time_len, eta_rho_len, xi_rho_len)) 
# Bottom peak wave period
tbr_tmp = np.empty((time_len, eta_rho_len, xi_rho_len))
tbr_tmp_mask = np.empty((time_len, eta_rho_len, xi_rho_len))
# Iterations for calculation
iter_tmp = np.empty((time_len, eta_rho_len, xi_rho_len))
iter_tmp_mask = np.empty((time_len, eta_rho_len, xi_rho_len))

# Loop through time, eta, and xi to calculate the bottom orbital
# velocity and wave period
# Loop through time
for tt in range(time_len):
    # Loop through eta
    for ee in range(eta_rho_len):
        ## Loop through xi
        #for xx in range(xi_rho_len):
            
        # Pull out the inputs for this grid cell
        # Significant wave height 
        #hs_tmp = dr_out_ww32rho_Thgt_nonan[tt,ee,:] # OG
        hs_tmp = dr_out_ww32rho_Thgt_nonan[tt,ee,:]*2.0 # double wave height version 
        #print('hs_tmp: ', hs_tmp)
        # Peak wave period 
        tp_tmp = dr_out_ww32rho_Tper_nonan[tt,ee,:]
        #print('tp_tmp: ', tp_tmp)
        # Water depth
        depth_tmp = grid.h[ee,:].values
        #print('depth_tmp: ', depth_tmp)
        
        
        # Call the function and save the output
        ubr_tmp[tt,ee,:], tbr_tmp[tt,ee,:], iter_tmp[tt,ee,:] = ubspecpar(hs_tmp, tp_tmp, depth_tmp, 'D')
            
    # Multiply the results by the rho mask
    ubr_tmp_mask[tt,:,:] = ubr_tmp[tt,:,:]*ds_out_rho.mask.values
    tbr_tmp_mask[tt,:,:] = tbr_tmp[tt,:,:]*ds_out_rho.mask.values
    iter_tmp_mask[tt,:,:] = iter_tmp[tt,:,:]*ds_out_rho.mask.values
    
    print('tt: ', tt, flush=True)
            
# Save these values to the netcdf
ubr_g[:,:,:] = ubr_tmp_mask[:,:,:] # bottom orbital velocity 
pwaveb_g[:,:,:] = tbr_tmp_mask[:,:,:]# bottom wave period

    

# Close the netcdfs
nc1.close()









