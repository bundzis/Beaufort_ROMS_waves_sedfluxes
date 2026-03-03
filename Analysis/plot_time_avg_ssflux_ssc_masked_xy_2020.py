#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 10:36:14 2023

@author: brun1463
"""

########## SSC Flux & SSC For Surf, 1m Above Seafloor, Depth-Integrated #########
# The purpose of this script is to plot the time-averaed suspended sediment flux
# as quivers on top of the time-averaged ssc concentration for the surface,
# 1 m above seafloor, and depth-integrated. 
#
# Notes:
# - This script needs to be run in xroms env since it uses xroms to interpolate
# in depth
##################################################################################



# Load in the packages
import numpy as np
import xarray as xr
import xroms
import xesmf as xe
from glob import glob
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import cmocean 


# Set a universal fontsize
fontsize = 25

# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize) 
matplotlib.rc('ytick', labelsize=fontsize)

# Prevent tick labels from overlapping
matplotlib.rcParams['xtick.major.pad'] = 12
matplotlib.rcParams['ytick.major.pad'] = 12

# Load in the grid
#grid = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Include/KakAKgrd_shelf_big010_smooth006.nc')
grid = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Grids/KakAKgrd_shelf_big010_smooth006.nc') # UPDATE PATH


# Pull out some dimensions
eta_rho_len = len(grid.eta_rho)
xi_rho_len = len(grid.xi_rho)
s_rho_len = int(20)

# Load in the rho masks 
#mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
#mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')


# Need: time-averaged ss flux in uv at surface, 1 m above seafloor, depth-integrated; 
# time-averaegd ssc at surface, 1 m above seafloor, depth-integrated (or averaged)
# Surface and depth-integrated/depth-averaged should be fairly straight forward form 
# previous ss flux and ssc plots. For 1 m above seafloor, need to first interpoalte 
# redgrd(u), regrid(v), and ssc to  1 m above seafloor, then use functions to calculate flux
# in u and v directions at this level/of the output, then can plot


# Before defining functions, need to find regridding weights to go from 
# u and v points to rho points 
# Set the input and output grids, and sepcify the lat/lon
# Since we are looking at u and v, we will use lon_u and lat_u as the primary lat/lon for the grid 
# and make another with lat_v and lon_v
# U input grid 
ds_in_u = grid.copy() # need to use lon_180 for this grid 
ds_in_u['lon'] = (('eta_u', 'xi_u'), ds_in_u.lon_u.values)
ds_in_u['lat'] = (('eta_u', 'xi_u'), ds_in_u.lat_u.values)
# V input grid 
ds_in_v = grid.copy() # need to use lon_180 for this grid 
ds_in_v['lon'] = (('eta_v', 'xi_v'), ds_in_v.lon_v.values)
ds_in_v['lat'] = (('eta_v', 'xi_v'), ds_in_v.lat_v.values)

# Output grid (ROMS rho grid)
ds_out_rho = grid.copy()
ds_out_rho['lat'] = (('eta_rho', 'xi_rho'), ds_out_rho.lat_rho.values)
ds_out_rho['lon'] = (('eta_rho', 'xi_rho'), ds_out_rho.lon_rho.values)

# Add masks 
# ex: ds["mask"] = xr.where(~np.isnan(ds["zeta"].isel(ocean_time=0)), 1, 0)
# Input grid (HYCOM)
# this is only a surface mask - which is what we want 
# U
ds_in_u['mask'] = (('eta_u', 'xi_u'), ds_in_u.mask_u.values)
# V
ds_in_v['mask'] = (('eta_v', 'xi_v'), ds_in_v.mask_v.values)
# Output grid (ROMS rho grid)
ds_out_rho['mask'] = (('eta_rho', 'xi_rho'), ds_out_rho.mask_rho.values)

# Regrid from u grid to rho grid with the masks included and extrapolation used 
regridder_u2rho = xe.Regridder(ds_in_u, ds_out_rho, method="bilinear", extrap_method='inverse_dist') #extrap_method="nearest_s2d"
regridder_u2rho
# Regrid from v grid to rho grid with the masks included and extrapolation used 
regridder_v2rho = xe.Regridder(ds_in_v, ds_out_rho, method="bilinear", extrap_method='inverse_dist') #extrap_method="nearest_s2d"
regridder_v2rho



# Make a bunch of functions 
# Surface
# Make a function to calculate the sediment flux in the 
# u direction for all sediment classes combined
def calc_u_surf_ssc_flux_allsed(filename, regridder_u2rho):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the surface suspended 
    sediment flux in the u direction for all sediment classes added
    together.

    Parameters
    ----------
    filename : The name/path of the model output file.

    Returns
    -------
    ssc_flux_u_surf_allsed: Surface SSC flux in u direction
    for all sediment classes combined 
    ssc_allsed_surf_tmp: Time series of ssc at the surface 

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
    
    # To collapse to horizontal, pull out hte surface values
    ssc_allsed_surf_tmp = ssc_allsed_tmp[:,-1,:,:]
    
    # Pull out the surface u velocities at all times, spaces
    u_surf_tmp = model_output.u[:,-1,:,:]
    
    # Interpolate them onto rho points 
    u_surf_tmp_rho = regridder_u2rho(u_surf_tmp)
    
    # Pull out the thickness of the cell in the y direction 
    #dy = 1.0/model_output.pn
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_u_surf_allsed = ((ssc_allsed_surf_tmp*u_surf_tmp_rho))
    
    # Return the surface u flux for all sediment classes
    return(ssc_flux_u_surf_allsed, ssc_allsed_surf_tmp)


# Make a function to calculate the sediment flux in the 
# v direction for all sediment classes combined
def calc_v_surf_ssc_flux_allsed(filename, regridder_v2rho):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the time-averaged, depth-integrated suspended 
    sediment flux in the v direction for all sediment classes added
    together.

    Parameters
    ----------
    filename : The name/path of the model output file.

    Returns
    -------
    ssc_flux_v_surf_allsed: Surface SSC flux in v direction
    for all sediment classes combined 

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
    
    # To collapse to horizontal, pull out hte surface values
    ssc_allsed_surf_tmp = ssc_allsed_tmp[:,-1,:,:]
    
    # Pull out the surface v velocities at all times, spaces
    v_surf_tmp = model_output.v[:,-1,:,:]
    
    # Interpolate them onto rho points 
    v_surf_tmp_rho = regridder_v2rho(v_surf_tmp)
    
    # Pull out the thickness of the cell in the x direction 
    #dx = 1.0/model_output.pm
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_v_surf_allsed = ((ssc_allsed_surf_tmp*v_surf_tmp_rho))
    
    
    # Return the surface v flux for all sediment classes
    return(ssc_flux_v_surf_allsed)


# 1 meter from Surface
# Make a function to calculate the sediment flux in the 
# u direction for all sediment classes combined
def calc_u_1m_fromsurf_ssc_flux_allsed(filename, regridder_u2rho, i):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the surface suspended 
    sediment flux in the u direction for all sediment classes added
    together 1 m from the surface.

    Parameters
    ----------
    filename : The name/path of the model output file.

    Returns
    -------
    ssc_flux_u_surf_allsed: Surface SSC flux in u direction
    for all sediment classes combined 
    ssc_allsed_surf_tmp: Time series of ssc at the surface 

    """
    
    # Load in the ROMS output 
    ds = xr.open_dataset(filename)
    ds['h'] = ds.bath
    ds, xgrid = xroms.roms_dataset(ds)
    
    # Set two depths, 1 and 2 m above seafloor
    #depths_surf = np.asarray([-1.0])
    #depths_surf = np.asarray([-2.0])
    depths_surf = np.asarray([-0.5, -1.5])
    
    # Add all the sediment classes together
    ssc_allsed_tmp = ds.mud_01 + ds.mud_02 + ds.sand_01 + ds.sand_02 + ds.sand_03
   
    # Interopolate onto 1 m from the surface
    # SSC
    ssc_allsed_1m_fromsurf_tmp = xroms.isoslice(ssc_allsed_tmp, depths_surf, xgrid, axis='Z')
    # U 
    u_interp_1m_fromsurf = xroms.isoslice(ds.u, depths_surf, xgrid, axis='Z')
    
    # Interpolate them onto rho points 
    u_surf_tmp_rho = regridder_u2rho(u_interp_1m_fromsurf)
    
    # If taking the value over two depths, take the average 
    if i == 6:
        u_surf_tmp_rho = np.nanmean(u_surf_tmp_rho, axis=0)
        ssc_allsed_1m_fromsurf_tmp = np.nanmean(ssc_allsed_1m_fromsurf_tmp, axis=0)
    else:
        u_surf_tmp_rho = np.nanmean(u_surf_tmp_rho, axis=1)
        ssc_allsed_1m_fromsurf_tmp = np.nanmean(ssc_allsed_1m_fromsurf_tmp, axis=1)
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_u_surf_allsed = ((ssc_allsed_1m_fromsurf_tmp*u_surf_tmp_rho))
    
    # Return the surface u flux for all sediment classes
    return(ssc_flux_u_surf_allsed, ssc_allsed_1m_fromsurf_tmp)


# Make a function to calculate the sediment flux in the 
# v direction for all sediment classes combined
def calc_v_1m_fromsurf_ssc_flux_allsed(filename, regridder_v2rho, i):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the time-averaged, depth-integrated suspended 
    sediment flux in the v direction for all sediment classes added
    together 1 m from the surface.

    Parameters
    ----------
    filename : The name/path of the model output file.

    Returns
    -------
    ssc_flux_v_surf_allsed: Surface SSC flux in v direction
    for all sediment classes combined 

    """
    
    # Load in the ROMS output 
    ds = xr.open_dataset(filename)
    ds['h'] = ds.bath
    ds, xgrid = xroms.roms_dataset(ds)
    
    # Set two depths, 1 and 2 m above seafloor
    #depths_surf = np.asarray([-1.0])
    #depths_surf = np.asarray([-2.0])
    depths_surf = np.asarray([-0.5, -1.5])
    
    # Add all the sediment classes together
    ssc_allsed_tmp = ds.mud_01 + ds.mud_02 + ds.sand_01 + ds.sand_02 + ds.sand_03

    # Interopolate onto 1 m from the surface
    # SSC
    ssc_allsed_1m_fromsurf_tmp = xroms.isoslice(ssc_allsed_tmp, depths_surf, xgrid, axis='Z')
    # V
    v_interp_1m_fromsurf = xroms.isoslice(ds.v, depths_surf, xgrid, axis='Z')
    
    # Interpolate them onto rho points 
    v_surf_tmp_rho = regridder_v2rho(v_interp_1m_fromsurf)
    
    # If taking the value over two depths, take the average 
    if i == 6:
        v_surf_tmp_rho = np.nanmean(v_surf_tmp_rho, axis=0)
        ssc_allsed_1m_fromsurf_tmp = np.nanmean(ssc_allsed_1m_fromsurf_tmp, axis=0)
    else:
        v_surf_tmp_rho = np.nanmean(v_surf_tmp_rho, axis=1)
        ssc_allsed_1m_fromsurf_tmp = np.nanmean(ssc_allsed_1m_fromsurf_tmp, axis=1)
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_v_surf_allsed = ((ssc_allsed_1m_fromsurf_tmp*v_surf_tmp_rho))
    
    # Return the surface u flux for all sediment classes
    return(ssc_flux_v_surf_allsed)


# Depth-integrated 
# Make a function to calculate the sediment flux in the 
# u direction for all sediment classes combined
def calc_u_depth_int_ssc_flux_allsed(filename, regridder_u2rho):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the time-averaged, depth-integrated suspended 
    sediment flux in the eastward direction for all sediment classes added
    together.

    Parameters
    ----------
    filename : The name/path of the model output file.

    Returns
    -------
    depth_int_ssc_flux_east_allsed: Depth-integrated SSC flux in eastward direction
    for all sediment classes combined 
    depth_int_ssc_allsed: Time series of depth-integrated ssc 

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
    
    # To collapse to horizontal, multiply each layer by its
    # thickness
    # Calculate the time-varying thickness of the cells
    dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
    
    # Pull out the u velocities at all times, depths, spaces
    u_tmp = model_output.u
    
    # Interpolate them onto rho points 
    u_tmp_rho = regridder_u2rho(u_tmp)
    
    # Pull out the thickness of the cell in the y direction 
    #dy = 1.0/model_output.pn
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_allsed = ((ssc_allsed_tmp*u_tmp_rho)*(dz))
    
    # Then depth-integrated by summing over depth and dividing by dy
    depth_int_ssc_flux_u_allsed = (ssc_flux_allsed.sum(dim='s_rho'))
    
    # Calculate depth-integrated ssc
    depth_int_ssc_allsed = (((ssc_allsed_tmp*dz)).sum(dim='s_rho'))
    
    # Divide by bathymetry to get depth-averaged SSC (kg/m3)
    depth_avg_ssc_allsed = depth_int_ssc_allsed/model_output.bath[:,:,:].values
    
    # Return the depth-integrated u flux for all sediment classes
    return(depth_int_ssc_flux_u_allsed, depth_avg_ssc_allsed)


# Make a function to calculate the sediment flux in the 
# v direction for all sediment classes combined
def calc_v_depth_int_ssc_flux_allsed(filename, regridder_v2rho):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the time-averaged, depth-integrated suspended 
    sediment flux in the v direction for all sediment classes added
    together.

    Parameters
    ----------
    filename : The name/path of the model output file.

    Returns
    -------
    depth_int_ssc_flux_north_allsed: Depth-integrated SSC flux in northward direction
    for all sediment classes combined 

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
    
    # To collapse to horizontal, multiply each layer by its
    # thickness
    # Calculate the time-varying thickness of the cells
    dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
    
    # Pull out the v velocities at all times, depths, spaces
    v_tmp = model_output.v
    
    # Interpolate them onto rho points 
    v_tmp_rho = regridder_v2rho(v_tmp)
    
    # Pull out the thickness of the cell in the x direction 
    #dx = 1.0/model_output.pm
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_allsed = ((ssc_allsed_tmp*v_tmp_rho)*(dz))
    
    # Then depth-integrated by summing over depth and dividing by dx
    depth_int_ssc_flux_v_allsed = (ssc_flux_allsed.sum(dim='s_rho'))
    
    # Return the depth-integrated v flux for all sediment classes
    return(depth_int_ssc_flux_v_allsed)



# 1 m above seafloor
# Make a function to interpolate ROMS u and v currents to 1 m above seafloor
def interp_roms_uv_to1m_rho(filename, regridder_u2rho, regridder_v2rho):
    """
    This function takes a given ROMS ocean_his file, opens it and 
    interpolates it onto the given depths, then returns the interpolated
    data. Right now, it is set up to interpolate u and v currents.
    

    Returns
    -------
    None.

    """
    
    # Load in the ROMS output 
    ds = xr.open_dataset(filename)
    ds['h'] = ds.bath
    ds, xgrid = xroms.roms_dataset(ds)
    
    # Get the depths as height from seabed
    height_from_seabed = ds.z_rho + ds.bath
    height_from_seabed.name = 'z_rho'
    
    # Set two depths, 1 and 2 m above seafloor
    depths = np.asarray([1.0])

    # Interopolate onto the given CODA depths
    u_rho_interp_1m = xroms.isoslice(regridder_u2rho(ds.u), depths, xgrid, 
                                         iso_array=height_from_seabed, axis='Z')
    v_rho_interp_1m = xroms.isoslice(regridder_v2rho(ds.v), depths, xgrid, 
                                         iso_array=height_from_seabed, axis='Z')
     
    # Return these currents 
    return(u_rho_interp_1m, v_rho_interp_1m)


# Make a function to interpolate ROMS u and v currents to 1 m above seafloor
def interp_roms_ssc_to1m_rho(filename):
    """
    This function takes a given ROMS ocean_his file, opens it and 
    interpolates it onto the given depths, then returns the interpolated
    data. Right now, it is set up to interpolate u and v currents.
    

    Returns
    -------
    None.

    """
    
    # Load in the ROMS output 
    ds = xr.open_dataset(filename)
    ds['h'] = ds.bath
    ds, xgrid = xroms.roms_dataset(ds)
    
    # Get the depths as height from seabed
    height_from_seabed = ds.z_rho + ds.bath
    height_from_seabed.name = 'z_rho'
    
    # Set two depths, 1 and 2 m above seafloor
    depths = np.asarray([1.0])

    # Interopolate onto the given CODA depths
    ssc_allsed_interp_1m = xroms.isoslice((ds.mud_01+ds.mud_02+ds.sand_01+ds.sand_02+ds.sand_03), depths, xgrid, 
                                         iso_array=height_from_seabed, axis='Z')
     
    # Return the ssc
    return(ssc_allsed_interp_1m)



# Define a function to pull out the length of time in the model run
# And the time steps
def get_model_time(filenames, num_files):
    """
    This function loops though model output and pulls
    out the entire length of the run, as well as the 
    individual time steps of the run.
    
    Inputs:
    - filenames: path and name of model output
    - num_files: the number of model output files
    
    Outputs:
    - time_len: length of time of model run (integer)
    - time_steps_list: list of time steps of full run (datetimes64)
    - time_lengths: array holding the lenght of time of each output file
    """

    # Create an array to hold the length of time in each output file
    time_lengths = np.empty((num_files))

    # Loop through output to pull out lenth of time
    for k in range(num_files):
        # Open the output file
        model_output = xr.open_dataset(filenames[k])

        # Pull out the length of time 
        time_lengths[k] = len(model_output.ocean_time)

    # Now sum over all the lengths to get total time
    time_len = np.sum(time_lengths, axis=0)

    # Convert from float to int
    time_len = int(time_len)

    # Loop back through the output to pull out the time step and save it
    # Make a list to hold the time steps 
    time_steps_list = []
    # Loop through output
    for h in range(num_files):
        # Open the output file
        model_output = xr.open_dataset(filenames[h])

        # Get the length of the run
        output_len = len(model_output.ocean_time)

        # Loop through each time step and append it to the list
        for g in range(output_len):
            time_steps_list.append(model_output.ocean_time[g].values)

    # Return this time length and time steps
    return(time_len, time_steps_list, time_lengths)




# Loop through model output and call the function
# First, get all the file names 
# ROMS 2020 output
# dbsed0003
#file_names = glob('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Output/dbsed0003/ocean_his_beaufort_shelf_2020_dbsed0003_*.nc')
file_names = glob('/Volumes/Model_outpu/Paper1/2020_dbsed0003/ocean_his_beaufort_shelf_2020_dbsed0003_*.nc')


# Sort them to be in order
file_names2 = sorted(file_names)

# Check to see if this worked
print(file_names2[0], flush=True)
print(file_names2[1], flush=True)
print(file_names2[2], flush=True)
print(file_names2[-1], flush=True)

# Pull out the number of files
num_files = len(file_names2)

# Pull out the length of time of the full run, the time steps, 
# and the length of time of each output file
full_time_len, time_steps, time_lengths = get_model_time(file_names2, num_files)


# Make some arrays to hold output
# Surface 
ssflux_u_surf_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ssc_surf_allsed_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ssflux_v_surf_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# Depth-integrated
ssflux_depth_int_u_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ssc_depth_avg_allsed_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ssflux_depth_int_v_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# 1 am above seafloor
u_1m_rho_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_1m_rho_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ssc_1m_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))

# Make an xarray dataset to merge dep_net timeseries into 
#dep_nets_xr = xr.Dataset()

# Set a time step to track which time step the loop is on
time_step = 0

# Loop through the model output
for j in range(num_files):
#for j in range(1):

    print('j: ', j)
    
    # Call the function to process the output
    # Get surface ss flux and ssc
    ssflux_u_surf_tmp, ssc_surf_allsed_tmp = calc_u_1m_fromsurf_ssc_flux_allsed(file_names2[j], regridder_u2rho, j)
    ssflux_v_surf_tmp = calc_v_1m_fromsurf_ssc_flux_allsed(file_names2[j], regridder_v2rho, j)
    
    # Get depth-integrated ss flux and ssc 
    ssflux_u_depth_int_tmp, ssc_depth_avg_allsed_tmp = calc_u_depth_int_ssc_flux_allsed(file_names2[j], regridder_u2rho)
    ssflux_v_depth_int_tmp = calc_v_depth_int_ssc_flux_allsed(file_names2[j], regridder_v2rho)
    
    # Get 1 m above seafloor currents and ssc
    u_1m_rho_tmp, v_1m_rho_tmp = interp_roms_uv_to1m_rho(file_names2[j], regridder_u2rho, regridder_v2rho)
    ssc_1m_tmp = interp_roms_ssc_to1m_rho(file_names2[j])
    
    
    # Save these to the arrays 
    #print('time_step: ', time_step)
    #print('time_step + time_lengths[j]: ', time_step+time_lengths[j])
    start = int(time_step)
    end = int(time_step+time_lengths[j])
    # Surface 
    ssflux_u_surf_full[start:end,:,:] = ssflux_u_surf_tmp
    ssflux_v_surf_full[start:end,:,:] = ssflux_v_surf_tmp
    ssc_surf_allsed_full[start:end,:,:] = ssc_surf_allsed_tmp
    # Depth-integrated
    ssflux_depth_int_u_full[start:end,:,:] = ssflux_u_depth_int_tmp
    ssflux_depth_int_v_full[start:end,:,:] = ssflux_v_depth_int_tmp
    ssc_depth_avg_allsed_full[start:end,:,:] = ssc_depth_avg_allsed_tmp
    # 1 m above seafloor 
    u_1m_rho_full[start:end,:,:] = u_1m_rho_tmp
    v_1m_rho_full[start:end,:,:] = v_1m_rho_tmp
    ssc_1m_full[start:end,:,:] = ssc_1m_tmp
    
    # Update the base time_step
    time_step = time_step + time_lengths[j]
    
    

# Now for 1 m above seafloor, get ss flux from time series
ssflux_u_1m_full = (u_1m_rho_full*ssc_1m_full)
ssflux_v_1m_full = (v_1m_rho_full*ssc_1m_full)

# Take the averages over time for all of these
# Surface 
ssflux_u_surf_avg = np.mean(ssflux_u_surf_full, axis=0)
ssflux_v_surf_avg = np.mean(ssflux_v_surf_full, axis=0)
ssc_surf_allsed_avg = np.mean(ssc_surf_allsed_full, axis=0)
# Depth-integrated
ssflux_depth_int_u_avg = np.mean(ssflux_depth_int_u_full, axis=0)
ssflux_depth_int_v_avg = np.mean(ssflux_depth_int_v_full, axis=0)
ssc_depth_avg_allsed_avg = np.mean(ssc_depth_avg_allsed_full, axis=0)
# Bottom 
ssflux_u_1m_avg = np.mean(ssflux_u_1m_full, axis=0)
ssflux_v_1m_avg = np.mean(ssflux_v_1m_full, axis=0)
ssc_1m_avg = np.mean(ssc_1m_full, axis=0)


# Now things should be ready to plot soooooooo try plotting...



# --------------------------------------------------------------------------------------
# ----- Plot 1: Time-Averaged Surface, 1m above seafloor, and depth-averaged ------------
# ------------------- currents and salinity, Masked XY ----------------------------------
# --------------------------------------------------------------------------------------
# Plot time-averaged surface, 1m above seafloor, and depth-averaged currents as 
# quivers and salinity as contours in background 

# Set the number of cells in the sponge on each open boundary
c_west = 36
c_north = 45
c_east = 36

# Make it so land will appear
temp_mask = grid.mask_rho.copy()
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)

# Prep the data by  ultiplying by the mask and trimming
# Multiply by mask
h_masked = grid.h.values*grid.mask_rho.values*mask_rho_nan.nudge_mask_rho_nan
# Trim 
lon_rho_trimmed = grid.lon_rho[:,c_west:-c_west].values
lat_rho_trimmed = grid.lat_rho[:,c_west:-c_west].values
h_masked_trimmed = h_masked[:,c_west:-c_west]

# Mask, trim, slice
# Mask
# Surface 
ssflux_u_surf_avg_wland_masked = ssflux_u_surf_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_v_surf_avg_wland_masked = ssflux_v_surf_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssc_surf_allsed_avg_wland_masked = ssc_surf_allsed_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 1 m 
ssflux_u_1m_avg_wland_masked = ssflux_u_1m_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_v_1m_avg_wland_masked = ssflux_v_1m_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssc_1m_avg_wland_masked = ssc_1m_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Depth-integrated 
ssflux_depth_int_u_avg_wland_masked = ssflux_depth_int_u_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_depth_int_v_avg_wland_masked = ssflux_depth_int_v_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssc_depth_avg_allsed_avg_masked = ssc_depth_avg_allsed_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask

# Trim
# Surface 
ssflux_u_surf_avg_wland_masked_trimmed = ssflux_u_surf_avg_wland_masked[:,c_west:-c_west]
ssflux_v_surf_avg_wland_masked_trimmed = ssflux_v_surf_avg_wland_masked[:,c_west:-c_west]
ssc_surf_allsed_avg_wland_masked_trimmed = ssc_surf_allsed_avg_wland_masked[:,c_west:-c_west]
# 1 m 
ssflux_u_1m_avg_wland_masked_trimmed = ssflux_u_1m_avg_wland_masked[:,c_west:-c_west]
ssflux_v_1m_avg_wland_masked_trimmed = ssflux_v_1m_avg_wland_masked[:,c_west:-c_west]
ssc_1m_avg_wland_masked_trimmed = ssc_1m_avg_wland_masked[:,c_west:-c_west]
# Depth-integrated
ssflux_depth_int_u_avg_wland_masked_trimmed = ssflux_depth_int_u_avg_wland_masked[:,c_west:-c_west]
ssflux_depth_int_v_avg_wland_masked_trimmed = ssflux_depth_int_v_avg_wland_masked[:,c_west:-c_west]
ssc_depth_avg_allsed_avg_trimmed = ssc_depth_avg_allsed_avg_masked[:,c_west:-c_west]

# Slice
nth_slice = 20 #15
# Surface 
ssflux_u_surf_avg_wland_masked_trimmed_slice = ssflux_u_surf_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
ssflux_v_surf_avg_wland_masked_trimmed_slice = ssflux_v_surf_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
# 1 m
ssflux_u_1m_avg_wland_masked_trimmed_slice = ssflux_u_1m_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
ssflux_v_1m_avg_wland_masked_trimmed_slice = ssflux_v_1m_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
# Depth-integrated
ssflux_depth_int_u_avg_wland_masked_trimmed_slice = ssflux_depth_int_u_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
ssflux_depth_int_v_avg_wland_masked_trimmed_slice = ssflux_depth_int_v_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
# Since fluxes will be a quiver plot, slice the data 
lon_rho_trimmed_slice = lon_rho_trimmed[::nth_slice,::nth_slice]
lat_rho_trimmed_slice = lat_rho_trimmed[::nth_slice,::nth_slice]

# Make a fake xy with the right resolution to be able to plot without the angle
x_rho_flat = np.arange(0,750*len(grid.x_rho[0,:]),750)
y_rho_flat = np.arange(0,600*len(grid.y_rho[:,0]),600)
# Prep the data by  ultiplying by the mask and trimming
# Trim 
x_rho_flat_trimmed = x_rho_flat[c_west:-c_west]
# Slice
x_rho_flat_trimmed_slice = x_rho_flat_trimmed[::nth_slice]
y_rho_flat_slice = y_rho_flat[::nth_slice]




# Set the colormap
cmap1=cmocean.cm.turbid
#cmap1.set_under('darkgray')
# Make the figure
fig1, ax1 = plt.subplots(3, figsize=(21, 18)) # (21,23)

# Set colorbar levels for all plots currents
lev1 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots salinity 
#lev2 = np.arange(0,0.06,0.001)
lev2 = np.arange(0,0.05,0.001)

# Surface 
# In the grid's u and v directions
# Plot bathymetry
# Plot surface ssc and flux
ax1[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax1[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs1 = ax1[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_surf_allsed_avg_wland_masked_trimmed, lev2, cmap=cmap1, extend='max')
# Plot bathymetry contours
ax1[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev1, colors='bisque')
# Plot currents
q1 = ax1[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, ssflux_u_surf_avg_wland_masked_trimmed_slice, 
                   ssflux_v_surf_avg_wland_masked_trimmed_slice, color='teal', width=2, 
                   angles='xy', scale_units='xy', units='xy')
ax1[0].quiverkey(q1, 0.65, 0.85, U=0.005, label='0.005 kg/m\u00b2s', fontproperties={'size':fontsize-2})

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax1[0].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m flux and ssc
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax1[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax1[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs2 = ax1[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_1m_avg_wland_masked_trimmed, lev2, cmap=cmap1, extend='max')
# Plot bathymetry contours 
ax1[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev1, colors='bisque')
# Plot currents
q2 = ax1[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_u_1m_avg_wland_masked_trimmed_slice, ssflux_v_1m_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax1[1].quiverkey(q2, 0.65, 0.85, U=0.001, label='0.001 kg/m\u00b2s', fontproperties={'size':fontsize-2})

# Label the plot
plt.setp(ax1[1].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[1].set_ylabel('Y (km)', fontsize=fontsize)

# Depth-integrated
# In the grid's u and v directions
# Plot depth-avg salinity 
ax1[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax1[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs3 = ax1[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_depth_avg_allsed_avg_trimmed, lev2, cmap=cmap1, extend='max')
# Plot bathymetry contours 
ax1[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev1, colors='bisque')
# Plot currents
q3 = ax1[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_depth_int_u_avg_wland_masked_trimmed_slice, ssflux_depth_int_v_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax1[2].quiverkey(q3, 0.65, 0.85, U=0.02, label='0.02 kg/ms', fontproperties={'size':fontsize-2})
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax1[2].get_xticklabels(), visible=True)
ax1[2].set_xlabel('X (km)', fontsize=fontsize)
ax1[2].set_ylabel('Y (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig1.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)

axes = ax1.ravel().tolist()
cbar1_ax = fig1.add_axes([0.85, bottom, 0.05, top-bottom])
cbar1 = plt.colorbar(cs1, ax=axes, cax=cbar1_ax, orientation='vertical').set_label(label='SSC (kg/m\u00b3)', size=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)



# =============================================================================
# # ------------------------------- Ratios -------------------------------------
# # Do some more analysis to see what percent of the SSC (and maybe sediment fluxes)
# # at each of these levels is mud01, mud02, sand01, and sand02
# 
# # Do the same interpolation and time-averaging but for each sediment class individually 
# # Make a function to interpolate ROMS ssc for given sediment class
# #  currents to 1 m above seafloor
# def interp_roms_ssc_to1m_onesed(filename, sedclass):
#     """
#     This function takes a given ROMS ocean_his file, opens it and 
#     interpolates it onto the given depths, then returns the interpolated
#     data. Right now, it is set up to interpolate u and v currents.
#     
# 
#     Returns
#     -------
#     None.
# 
#     """
#     
#     # Load in the ROMS output 
#     ds = xr.open_dataset(filename)
#     ds['h'] = ds.bath
#     ds, xgrid = xroms.roms_dataset(ds)
#     
#     # Get the depths as height from seabed
#     height_from_seabed = ds.z_rho + ds.bath
#     height_from_seabed.name = 'z_rho'
#     
#     # Set two depths, 1 and 2 m above seafloor
#     depths = np.asarray([1.0])
# 
#     # Interopolate onto the given CODA depths
#     ssc_allsed_interp_1m = xroms.isoslice(ds[sedclass], depths, xgrid, 
#                                          iso_array=height_from_seabed, axis='Z')
#      
#     # Return the ssc
#     return(ssc_allsed_interp_1m)
# 
# 
# # Make a function to pull out the surface values for SSC for each class
# def get_onesed_surf_timeseries(filename, sedclass):
#     """
#     
# 
#     Parameters
#     ----------
#     filename : TYPE
#         DESCRIPTION.
#     sedcalss : TYPE
#         DESCRIPTION.
# 
#     Returns
#     -------
#     None.
# 
#     """
#     
#     # Load in the ROMS output 
#     ds = xr.open_dataset(filename)
#     
#     # Get the sediment time series
#     sed_tmp = ds[sedclass][:,-1,:,:].values
#     
#     # Return this 
#     return(sed_tmp)
#     
#     
#     
# # Make a function to calculate the depth-averaged SSC for a given sediment class
# def calc_depth_avg_ssc_onesed(filename, sedclass):
#     """
#     The purpose of this function is to take a given model output file, load 
#     in the output, and caluclate the depth-averaged SSC for a given sediment class.
# 
#     Parameters
#     ----------
#     filename : The name/path of the model output file.
# 
#     Returns
#     -------
#     depth_avg_ssc_onesed: Time series of depth-integrated ssc 
# 
#     """
#     
#     # Load in the model output
#     model_output = xr.open_dataset(filename)
#     
#     # Add all the sediment classes together
#     ssc_onesed_tmp = model_output[sedclass] 
#     
#     # To collapse to horizontal, multiply each layer by its
#     # thickness
#     # Calculate the time-varying thickness of the cells
#     dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
#     
#     # Multiply by thickness
#     ssc_onesed_thick_tmp = ssc_onesed_tmp*dz
#     
#     # Sum over depth to get depth-integrated 
#     ssc_onesed_depth_int_tmp = ssc_onesed_thick_tmp.sum(dim='s_rho')
#     
#     # Divide by depth to get depth-averaged 
#     ssc_onesed_depth_avg_tmp = ssc_onesed_depth_int_tmp/model_output.bath[:,:,:].values
#     
#     # Return the depth-averaged ssc for the given class
#     return(ssc_onesed_depth_avg_tmp)
# 
# 
#     
# 
# # Delete some variables for memory
# # Now for 1 m above seafloor, get ss flux from time series
# del(ssflux_u_1m_full)
# del(ssflux_v_1m_full)
# 
# # Take the averages over time for all of these
# # Surface 
# del(ssflux_u_surf_full)
# del(ssflux_v_surf_full)
# del(ssc_surf_allsed_full)
# # Depth-integrated
# del(ssflux_depth_int_u_full)
# del(ssflux_depth_int_v_full)
# del(ssc_depth_avg_allsed_full)
# # Bottom 
# del(ssc_1m_full)
#     
# 
# 
# # Call the function for the classes
# # Make some arrays to hold output
# # Surface 
# mud01_surf_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# mud02_surf_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# sand01_surf_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# sand02_surf_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# # Depth-integrated
# mud01_davg_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# mud02_davg_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# sand01_davg_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# sand02_davg_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# # 1 am above seafloor
# mud01_1m_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# mud02_1m_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# sand01_1m_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# sand02_1m_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# 
# # Make an xarray dataset to merge dep_net timeseries into 
# #dep_nets_xr = xr.Dataset()
# 
# # Set a time step to track which time step the loop is on
# time_step2 = 0
# 
# # Loop through the model output
# for j in range(num_files):
# #for j in range(1):
# 
#     print('j: ', j)
#     
#     # Call the function to process the output
#     # Get surface sediments 
#     mud01_surf_tmp = get_onesed_surf_timeseries(file_names2[j], 'mud_01')
#     mud02_surf_tmp = get_onesed_surf_timeseries(file_names2[j], 'mud_02')
#     sand01_surf_tmp = get_onesed_surf_timeseries(file_names2[j], 'sand_01')
#     sand02_surf_tmp = get_onesed_surf_timeseries(file_names2[j], 'sand_02')
#     
#     # Get depth-averaged for each class
#     mud01_davg_tmp = calc_depth_avg_ssc_onesed(file_names2[j], 'mud_01')
#     mud02_davg_tmp = calc_depth_avg_ssc_onesed(file_names2[j], 'mud_02')
#     sand01_davg_tmp = calc_depth_avg_ssc_onesed(file_names2[j], 'sand_01')
#     sand02_davg_tmp = calc_depth_avg_ssc_onesed(file_names2[j], 'sand_02')
#     
#     # Get 1 m above seafloor for each class
#     mud01_1m_tmp = interp_roms_ssc_to1m_onesed(file_names2[j], 'mud_01')
#     mud02_1m_tmp = interp_roms_ssc_to1m_onesed(file_names2[j], 'mud_02')
#     sand01_1m_tmp = interp_roms_ssc_to1m_onesed(file_names2[j], 'sand_01')
#     sand02_1m_tmp = interp_roms_ssc_to1m_onesed(file_names2[j], 'sand_02')
#     
#     
#     # Save these to the arrays 
#     #print('time_step: ', time_step)
#     #print('time_step + time_lengths[j]: ', time_step+time_lengths[j])
#     start2 = int(time_step2)
#     end2 = int(time_step2+time_lengths[j])
#     # Surface 
#     mud01_surf_full[start2:end2,:,:] = mud01_surf_tmp
#     mud02_surf_full[start2:end2,:,:] = mud02_surf_tmp
#     sand01_surf_full[start2:end2,:,:] = sand01_surf_tmp
#     sand02_surf_full[start2:end2,:,:] = sand02_surf_tmp
#     # Depth-averaged
#     mud01_davg_full[start2:end2,:,:] = mud01_davg_tmp
#     mud02_davg_full[start2:end2,:,:] = mud02_davg_tmp
#     sand01_davg_full[start2:end2,:,:] = sand01_davg_tmp
#     sand02_davg_full[start2:end2,:,:] = sand02_davg_tmp
#     # 1 m above seafloor 
#     mud01_1m_full[start2:end2,:,:] = mud01_1m_tmp
#     mud02_1m_full[start2:end2,:,:] = mud02_1m_tmp
#     sand01_1m_full[start2:end2,:,:] = sand01_1m_tmp
#     sand02_1m_full[start2:end2,:,:] = sand02_1m_tmp
#     
#     # Update the base time_step
#     time_step2 = time_step2 + time_lengths[j]
# 
# 
# 
# # Take the average over time
# # Surface
# mud01_surf_avg = np.mean(mud01_surf_tmp, axis=0)
# mud02_surf_avg = np.mean(mud02_surf_tmp, axis=0)
# sand01_surf_avg = np.mean(sand01_surf_tmp, axis=0)
# sand02_surf_avg = np.mean(sand02_surf_tmp, axis=0)
# # Depth-averaged 
# mud01_davg_tavg = np.mean(mud01_davg_full, axis=0)
# mud02_davg_tavg = np.mean(mud02_davg_full, axis=0)
# sand01_davg_tavg = np.mean(sand01_davg_full, axis=0)
# sand02_davg_tavg = np.mean(sand02_davg_full, axis=0)
# # 1 m above seafloor
# mud01_1m_avg = np.mean(mud01_1m_full, axis=0)
# mud02_1m_avg = np.mean(mud02_1m_full, axis=0)
# sand01_1m_avg = np.mean(sand01_1m_full, axis=0)
# sand02_1m_avg = np.mean(sand02_1m_full, axis=0)
# 
# # Take the ratio of this over the ones from above 
# # -- Surface --
# frac_mud01_surf = mud01_surf_avg/ssc_surf_allsed_avg
# frac_mud02_surf = mud02_surf_avg/ssc_surf_allsed_avg
# frac_sand01_surf = sand01_surf_avg/ssc_surf_allsed_avg
# frac_sand02_surf = sand02_surf_avg/ssc_surf_allsed_avg
# #print('frac_mud01_surf: ', frac_mud01_surf)
# #print('frac_mud02_surf: ', frac_mud02_surf)
# #print('frac_sand01_surf: ', frac_sand01_surf)
# #print('frac_sand02_surf: ', frac_sand02_surf)
# # Mean
# print('frac_mud01_surf mean: ', np.nanmean(frac_mud01_surf))
# print('frac_mud02_surf mean: ', np.nanmean(frac_mud02_surf))
# print('frac_sand01_surf mean: ', np.nanmean(frac_sand01_surf))
# print('frac_sand02_surf mean: ', np.nanmean(frac_sand02_surf))
# # Std dev
# print('frac_mud01_surf stddev: ', np.nanstd(frac_mud01_surf))
# print('frac_mud02_surf stddev: ', np.nanstd(frac_mud02_surf))
# print('frac_sand01_surf stddev: ', np.nanstd(frac_sand01_surf))
# print('frac_sand02_surf stddev: ', np.nanstd(frac_sand02_surf))
# # Mins 
# print('frac_mud01_surf min: ', np.nanmin(frac_mud01_surf))
# print('frac_mud02_surf min: ', np.nanmin(frac_mud02_surf))
# print('frac_sand01_surf min: ', np.nanmin(frac_sand01_surf))
# print('frac_sand02_surf min: ', np.nanmin(frac_sand02_surf))
# # Maxes
# print('frac_mud01_surf max: ', np.nanmax(frac_mud01_surf))
# print('frac_mud02_surf max: ', np.nanmax(frac_mud02_surf))
# print('frac_sand01_surf max: ', np.nanmax(frac_sand01_surf))
# print('frac_sand02_surf max: ', np.nanmax(frac_sand02_surf))
# # -- Depth-averaged --
# frac_mud01_depth_avg = mud01_davg_tavg/ssc_depth_avg_allsed_avg
# frac_mud02_depth_avg = mud02_davg_tavg/ssc_depth_avg_allsed_avg
# frac_sand01_depth_avg = sand01_davg_tavg/ssc_depth_avg_allsed_avg
# frac_sand02_depth_avg = sand02_davg_tavg/ssc_depth_avg_allsed_avg
# #print('frac_mud01_depth_avg: ', frac_mud01_depth_avg)
# #print('frac_mud02_depth_avg: ', frac_mud02_depth_avg)
# #print('frac_sand01_depth_avg: ', frac_sand01_depth_avg)
# #print('frac_sand02_depth_avg: ', frac_sand02_depth_avg)
# # Mean
# print('frac_mud01_depth_avg mean: ', np.nanmean(frac_mud01_depth_avg))
# print('frac_mud02_depth_avg mean: ', np.nanmean(frac_mud02_depth_avg))
# print('frac_sand01_depth_avg mean: ', np.nanmean(frac_sand01_depth_avg))
# print('frac_sand02_depth_avg mean: ', np.nanmean(frac_sand02_depth_avg))
# # Std dev
# print('frac_mud01_depth_avg stddev: ', np.nanstd(frac_mud01_depth_avg))
# print('frac_mud02_depth_avg stddev: ', np.nanstd(frac_mud02_depth_avg))
# print('frac_sand01_depth_avg stddev: ', np.nanstd(frac_sand01_depth_avg))
# print('frac_sand02_depth_avg stddev: ', np.nanstd(frac_sand02_depth_avg))
# # Mins 
# print('frac_mud01_depth_avg min: ', np.nanmin(frac_mud01_depth_avg))
# print('frac_mud02_depth_avg min: ', np.nanmin(frac_mud02_depth_avg))
# print('frac_sand01_depth_avg min: ', np.nanmin(frac_sand01_depth_avg))
# print('frac_sand02_depth_avg min: ', np.nanmin(frac_sand02_depth_avg))
# # Maxes
# print('frac_mud01_depth_avg max: ', np.nanmax(frac_mud01_depth_avg))
# print('frac_mud02_depth_avg max: ', np.nanmax(frac_mud02_depth_avg))
# print('frac_sand01_depth_avg max: ', np.nanmax(frac_sand01_depth_avg))
# print('frac_sand02_depth_avg max: ', np.nanmax(frac_sand02_depth_avg))
# # -- 1 m above seafloor --
# frac_mud01_1m = mud01_1m_avg/ssc_1m_avg
# frac_mud02_1m = mud02_1m_avg/ssc_1m_avg
# frac_sand01_1m = sand01_1m_avg/ssc_1m_avg
# frac_sand02_1m = sand02_1m_avg/ssc_1m_avg
# #print('frac_mud01_1m: ', frac_mud01_1m)
# #print('frac_mud02_1m: ', frac_mud02_1m)
# #print('frac_sand01_1m: ', frac_sand01_1m)
# #print('frac_sand02_1m: ', frac_sand02_1m)
# # Mean
# print('frac_mud01_1m mean: ', np.nanmean(frac_mud01_1m))
# print('frac_mud02_1m mean: ', np.nanmean(frac_mud02_1m))
# print('frac_sand01_1m mean: ', np.nanmean(frac_sand01_1m))
# print('frac_sand02_1m mean: ', np.nanmean(frac_sand02_1m))
# # Std dev
# print('frac_mud01_1m stddev: ', np.nanstd(frac_mud01_1m))
# print('frac_mud02_1m stddev: ', np.nanstd(frac_mud02_1m))
# print('frac_sand01_1m stddev: ', np.nanstd(frac_sand01_1m))
# print('frac_sand02_1m stddev: ', np.nanstd(frac_sand02_1m))
# # Mins 
# print('frac_mud01_1m min: ', np.nanmin(frac_mud01_1m))
# print('frac_mud02_1m min: ', np.nanmin(frac_mud02_1m))
# print('frac_sand01_1m min: ', np.nanmin(frac_sand01_1m))
# print('frac_sand02_surf min: ', np.nanmin(frac_sand02_1m))
# # Maxes
# print('frac_mud01_1m max: ', np.nanmax(frac_mud01_1m))
# print('frac_mud02_1m max: ', np.nanmax(frac_mud02_1m))
# print('frac_sand01_1m max: ', np.nanmax(frac_sand01_1m))
# print('frac_sand02_1m max: ', np.nanmax(frac_sand02_1m))
# 
# 
# # Plot these
# # --------------------------------------------------------------------------------------
# # ----- Plot 2: Time-Averaged Surface Fractions of SSC ------------
# # --------------------------------------------------------------------------------------
# # Make a subplot of the fraction of the SSC in the surface that is mud01, mud02, sand01,
# # and sand02
# 
# # Set the colormap
# cmap2=cmocean.cm.tempo
# #cmap1.set_under('darkgray')
# # Make the figure
# fig2, ax2 = plt.subplots(2, 2, figsize=(28, 10)) # (21,23)
# fig2.suptitle('Surface', fontsize=fontsize+2)
# 
# # Set colorbar levels for all plots currents
# lev3 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# # Set colorbar levels for all plots salinity 
# lev4 = np.arange(0, 1.0, 0.01)
# 
# # Mud01
# cs4 = ax2[0,0].contourf(frac_mud01_surf, levels=lev4, cmap=cmap2, extend='both')
# ax2[0,0].contour(grid.h.values, lev3, colors='deeppink')
# ax2[0,0].set_title('Mud01', fontsize=fontsize)
# #ax2[0,0].set_xlabel('xi', fontsize=fontsize)
# ax2[0,0].set_ylabel('eta', fontsize=fontsize)
# ax2[0,0].set_facecolor('darkgray')
# 
# # Mud02
# cs5 = ax2[0,1].contourf(frac_mud02_surf, levels=lev4, cmap=cmap2, extend='both')
# ax2[0,1].contour(grid.h.values, lev3, colors='deeppink')
# ax2[0,1].set_title('Mud02', fontsize=fontsize)
# #ax2[0,1].set_xlabel('xi', fontsize=fontsize)
# #ax2[0,1].set_ylabel('eta', fontsize=fontsize)
# ax2[0,1].set_facecolor('darkgray')
# 
# # Sand01
# cs6 = ax2[1,0].contourf(frac_sand01_surf, levels=lev4, cmap=cmap2, extend='both')
# ax2[1,0].contour(grid.h.values, lev3, colors='deeppink')
# ax2[1,0].set_title('Sand01', fontsize=fontsize)
# ax2[1,0].set_xlabel('xi', fontsize=fontsize)
# ax2[1,0].set_ylabel('eta', fontsize=fontsize)
# ax2[1,0].set_facecolor('darkgray')
# 
# # Sand02
# cs7 = ax2[1,1].contourf(frac_sand02_surf, levels=lev4, cmap=cmap2, extend='both')
# ax2[1,1].contour(grid.h.values, lev3, colors='deeppink')
# ax2[1,1].set_title('Sand02', fontsize=fontsize)
# ax2[1,1].set_xlabel('xi', fontsize=fontsize)
# #ax2[1,1].set_ylabel('eta', fontsize=fontsize)
# ax2[1,1].set_facecolor('darkgray')
# 
# # Adjust colorbar placement
# bottom, top = 0.1, 0.95
# left, right = 0.1, 0.8
# fig2.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.16, hspace=0.3)
# fig2.subplots_adjust(top=0.88)
# 
# axes = ax2.ravel().tolist()
# cbar2_ax = fig2.add_axes([0.82, bottom, 0.02, top-bottom])
# cbar2 = plt.colorbar(cs4, ax=axes, cax=cbar2_ax, orientation='vertical').set_label(label='Ratio Over Total', size=fontsize)
# 
# 
# 
# # --------------------------------------------------------------------------------------
# # ----- Plot 3: Time-Averaged Depth-Avg Fractions of SSC ------------
# # --------------------------------------------------------------------------------------
# # Make a subplot of the fraction of the SSC depth-averaged that is mud01, mud02, sand01,
# # and sand02
# 
# # Set the colormap
# cmap3=cmocean.cm.tempo
# #cmap1.set_under('darkgray')
# # Make the figure
# fig3, ax3 = plt.subplots(2, 2, figsize=(28, 10)) # (21,23)
# fig3.suptitle('Depth-Averaged', fontsize=fontsize+2)
# 
# # Set colorbar levels for all plots
# lev5 = np.arange(0, 1.0, 0.01)
# 
# # Mud01
# cs8 = ax3[0,0].contourf(frac_mud01_depth_avg, levels=lev5, cmap=cmap3, extend='both')
# ax3[0,0].contour(grid.h.values, lev3, colors='deeppink')
# ax3[0,0].set_title('Mud01', fontsize=fontsize)
# #ax2[0,0].set_xlabel('xi', fontsize=fontsize)
# ax3[0,0].set_ylabel('eta', fontsize=fontsize)
# ax3[0,0].set_facecolor('darkgray')
# 
# # Mud02
# cs9 = ax3[0,1].contourf(frac_mud02_depth_avg, levels=lev5, cmap=cmap3, extend='both')
# ax3[0,1].contour(grid.h.values, lev3, colors='deeppink')
# ax3[0,1].set_title('Mud02', fontsize=fontsize)
# #ax2[0,1].set_xlabel('xi', fontsize=fontsize)
# #ax2[0,1].set_ylabel('eta', fontsize=fontsize)
# ax3[0,1].set_facecolor('darkgray')
# 
# # Sand01
# cs10 = ax3[1,0].contourf(frac_sand01_depth_avg, levels=lev5, cmap=cmap3, extend='both')
# ax3[1,0].contour(grid.h.values, lev3, colors='deeppink')
# ax3[1,0].set_title('Sand01', fontsize=fontsize)
# ax3[1,0].set_xlabel('xi', fontsize=fontsize)
# ax3[1,0].set_ylabel('eta', fontsize=fontsize)
# ax3[1,0].set_facecolor('darkgray')
# 
# # Sand02
# cs11 = ax3[1,1].contourf(frac_sand02_depth_avg, levels=lev5, cmap=cmap3, extend='both')
# ax3[1,1].contour(grid.h.values, lev3, colors='deeppink')
# ax3[1,1].set_title('Sand02', fontsize=fontsize)
# ax3[1,1].set_xlabel('xi', fontsize=fontsize)
# #ax2[1,1].set_ylabel('eta', fontsize=fontsize)
# ax3[1,1].set_facecolor('darkgray')
# 
# # Adjust colorbar placement
# bottom, top = 0.1, 0.95
# left, right = 0.1, 0.8
# fig3.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.16, hspace=0.3)
# fig3.subplots_adjust(top=0.88)
# 
# axes = ax3.ravel().tolist()
# cbar3_ax = fig3.add_axes([0.82, bottom, 0.02, top-bottom])
# cbar3 = plt.colorbar(cs8, ax=axes, cax=cbar3_ax, orientation='vertical').set_label(label='Ratio Over Total', size=fontsize)
# 
# 
# 
# 
# # --------------------------------------------------------------------------------------
# # ----- Plot 4: Time-Averaged 1 m Above Seafloor Fractions of SSC ------------
# # --------------------------------------------------------------------------------------
# # Make a subplot of the fraction of the SSC 1 m above seafloor that is mud01, mud02, sand01,
# # and sand02
# 
# # Set the colormap
# cmap4=cmocean.cm.tempo
# #cmap1.set_under('darkgray')
# # Make the figure
# fig4, ax4 = plt.subplots(2, 2, figsize=(28, 10)) # (21,23)
# fig4.suptitle('1 m Above Seafloor', fontsize=fontsize+2)
# 
# # Set colorbar levels for all plots
# lev6 = np.arange(0, 1.0, 0.01)
# 
# # Mud01
# cs12 = ax4[0,0].contourf(frac_mud01_1m, levels=lev6, cmap=cmap4, extend='both')
# ax4[0,0].contour(grid.h.values, lev3, colors='deeppink')
# ax4[0,0].set_title('Mud01', fontsize=fontsize)
# #ax2[0,0].set_xlabel('xi', fontsize=fontsize)
# ax4[0,0].set_ylabel('eta', fontsize=fontsize)
# ax4[0,0].set_facecolor('darkgray')
# 
# # Mud02
# cs13 = ax4[0,1].contourf(frac_mud02_1m, levels=lev6, cmap=cmap4, extend='both')
# ax4[0,1].contour(grid.h.values, lev3, colors='deeppink')
# ax4[0,1].set_title('Mud02', fontsize=fontsize)
# #ax2[0,1].set_xlabel('xi', fontsize=fontsize)
# #ax2[0,1].set_ylabel('eta', fontsize=fontsize)
# ax4[0,1].set_facecolor('darkgray')
# 
# # Sand01
# cs14 = ax4[1,0].contourf(frac_sand01_1m, levels=lev6, cmap=cmap4, extend='both')
# ax4[1,0].contour(grid.h.values, lev3, colors='deeppink')
# ax4[1,0].set_title('Sand01', fontsize=fontsize)
# ax4[1,0].set_xlabel('xi', fontsize=fontsize)
# ax4[1,0].set_ylabel('eta', fontsize=fontsize)
# ax4[1,0].set_facecolor('darkgray')
# 
# # Sand02
# cs15 = ax4[1,1].contourf(frac_sand02_1m, levels=lev6, cmap=cmap4, extend='both')
# ax4[1,1].contour(grid.h.values, lev3, colors='deeppink')
# ax4[1,1].set_title('Sand02', fontsize=fontsize)
# ax4[1,1].set_xlabel('xi', fontsize=fontsize)
# #ax2[1,1].set_ylabel('eta', fontsize=fontsize)
# ax4[1,1].set_facecolor('darkgray')
# 
# # Adjust colorbar placement
# bottom, top = 0.1, 0.95
# left, right = 0.1, 0.8
# fig4.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.16, hspace=0.3)
# fig4.subplots_adjust(top=0.88)
# 
# axes = ax4.ravel().tolist()
# cbar4_ax = fig4.add_axes([0.82, bottom, 0.02, top-bottom])
# cbar4 = plt.colorbar(cs12, ax=axes, cax=cbar4_ax, orientation='vertical').set_label(label='Ratio Over Total', size=fontsize)
# 
# 
# =============================================================================

################### Add On Net Erosion & Deposition Analyses then #############
############################# Plot Together ###################################

# +++++++++++++++++++ Method B: Use ROMS Bed_Thickness +++++++++++++++++++++++++++++++++++
# Find net erosion/deposition by using bed_thickness from ROMS since dep_net is elusive

# Define a function to loop through the model output and calculate the change 
# in kg of sediment in the seabed (all sediment classes added together) and 
# use this mass to back out the depth of seabed change from that mass
def get_mass_and_thickness_change_allsed2(filenames, num_files, time_len, eta_rho_len, xi_rho_len):
    """
    This function takes all of the output files from ROMS, loops through them to pull out
    the change in the mass of sediment in the seabed, both compared to the start of the 
    run and in between time steps. It also pulls out the total mass of sediment in the 
    seabed at each time step and the erosion/deposition in m of the seabed between
    time steps (like dep_net). This function does this for all sediment classes
    combined.
    
    Right now, this function does not use porosity in any way but the output will be 
    cmompared with dep_net to try to figure out where porosity should come into play.

    Parameters
    ----------
    filenames : Names/paths of all of the model output files 
    num_files : Number of model output files
    time_len : length of time of the full model run (# of time steps)
    eta_rho_len : Length of eta_rho points in the grid
    xi_rho_len : Length of xi_rho points in the grid 

    Returns
    -------
    seabed_mass_allsed : Total mass of sediment in seabed at each time step in space
                         (all sediment classes added together, [ocean_time, eta_rho, xi_rho])
    seabed_mass_allsed_diff_from_start : Difference in total mass of sediment in seabed
                                         since start of the model run (all sediment classes
                                         added together, [ocean_time, eta_rho, xi_rho])
    seabed_mass_allsed_diff_between_steps : Difference in total mass of sediment between time
                                            steps (all sediment classes added together,
                                            [ocean_time, eta_rho, xi_rho])
    seabed_erodepo_allsed_m : Erosion/Deposition of seabed in m between time steps
                              (all sediment classes added together, [ocean_time, eta_rho, xi_rho])

    """
    
    # Before looping through all files to get the mass of sediment in the seabed
    # at each time step, make empy arrays to hold this value for the full model run
    # Time series of depth of the seabed from all seds
    depth_seabed_allsed = np.empty((time_len, eta_rho_len, xi_rho_len))
    # Time series of bed thickness
    bed_thickness = np.empty((time_len, eta_rho_len, xi_rho_len))
    
    # Set a time step to track which time step the loop is on
    time_step = 0
    
    # Loop through the model output
    for i in range(num_files):
        # Open the model output 
        model_output = xr.open_dataset(filenames[i])
        
        # Pull out the length of time in this output file
        time_len_tmp = len(model_output.ocean_time)
        
        # Set densities and porosities for each class
        rho_mud01 = 1750 # kg/m3
        rho_mud02 = 1750 # kg/m3
        rho_sand01 = 2650 # kg/m3
        rho_sand02 = 2650 # kg/m3
        rho_sand03 = 2650 # kg/m3
        
        por_mud01 = 0.85
        por_mud02 = 0.85
        por_sand01 = 0.5
        por_sand02 = 0.5
        por_sand03 = 0.45
        
        # For each class, divide the mass in each layer by the density times
        # (1- porosity)
        depth_mud01 = ((model_output.mudmass_01)/(rho_mud01*(1-por_mud01))) # m
        depth_mud02 = ((model_output.mudmass_02)/(rho_mud02*(1-por_mud02))) # m
        depth_sand01 = ((model_output.sandmass_01)/(rho_sand01*(1-por_sand01))) # m
        depth_sand02 = ((model_output.sandmass_02)/(rho_sand02*(1-por_sand02))) # m
        depth_sand03 = ((model_output.sandmass_03)/(rho_sand03*(1-por_sand03))) # m
        
        # Add all of the classes together
        depth_all = depth_mud01 + depth_mud02 + depth_sand01 + depth_sand02 + depth_sand03
        
        # Sum over all layers to get a net depth 
        depth_all_tot = depth_all.sum(dim='Nbed')
        
        # Pull out bed thickness
        bed_thick = model_output.bed_thickness
        
        # Sum over all layers to get total thickness
        bed_thick_tot = bed_thick.sum(dim='Nbed')
        
        # Now we have the depth of the seabed from all of the sediment classes
        # combined 
        # Save this to the array
        start = int(time_step)
        end = int(time_step+time_len_tmp)
        depth_seabed_allsed[start:end,:,:] = depth_all_tot
        bed_thickness[start:end,:,:] = bed_thick_tot
                 
        # Update the time step
        time_step = time_step + time_len_tmp
        
    # To compare this with dep net, take the difference between neighbors 
    depth_tot_all_diff = depth_seabed_allsed[1:,:,:] - depth_seabed_allsed[:-1,:,:]
    
    # Return the time series of the total depth of the seabed at each point 
    # and the time series of the differences in depth between time neighbors
    # and return time series of bed thickness to compare
    return(depth_seabed_allsed, depth_tot_all_diff, bed_thickness)

# Do improved calculation of dep_net/seabed change since dep_net is elusive
# Get the depth of the seabed at each time step and the difference between time steps 
seabed_depth_allsed, seabed_depth_diff_allsed, roms_bed_thick = get_mass_and_thickness_change_allsed2(file_names2, num_files, full_time_len, eta_rho_len, xi_rho_len)
# Find the net erosion/deposition by taking the difference in bed thickness between
# the end and beginning of the run 
roms_net_erodepo_bedthick = roms_bed_thick[-1,:,:] - roms_bed_thick[0,:,:]


# ----------------------------------------------------------------------------------------------
# -------------------------- Plot 5: Sediment Flux, SSC, and -----------------------------------
# ----------------------- Net Erosion & Deposition Masked XY -----------------------------------
# ----------------------------------------------------------------------------------------------
# Plot the ssc flux, ssc, and net erosiion deposition but masking the 
# regions we don't trust and in xy
# First prep the ero depo data
roms_net_erodepo_bedthick_wland_nan = roms_net_erodepo_bedthick*temp_mask
roms_net_erodepo_bedthick_wland_nan = np.where(np.isnan(roms_net_erodepo_bedthick_wland_nan), -100, roms_net_erodepo_bedthick_wland_nan)

# Make it so land appears 
# Mask, trim
roms_net_erodepo_bedthick_wland_masked = roms_net_erodepo_bedthick*temp_mask*mask_rho_nan.nudge_mask_rho_nan
roms_net_erodepo_bedthick_wland_nan_masked = roms_net_erodepo_bedthick_wland_nan*mask_rho_nan.nudge_mask_rho_nan
roms_net_erodepo_bedthick_wland_masked_trimmed = roms_net_erodepo_bedthick_wland_masked[:,c_west:-c_west]
roms_net_erodepo_bedthick_wland_nan_masked_trimmed = roms_net_erodepo_bedthick_wland_nan_masked[:,c_west:-c_west]


# Make the figure
fig5, ax5 = plt.subplots(4, figsize=(22,21)) # (18,8) (26,12) (26,8) (26,10)

# Set the colormaps
cmap5 = cmocean.cm.turbid
cmap6 = cmocean.cm.delta # ero depo

# Set colorbar levels for all plots currents
lev5 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots salinity 
#lev2 = np.arange(0,0.06,0.001)
lev6 = np.arange(0,0.05,0.001) # ssc
lev7 = np.arange(-1, 1, 0.01) # ero depo

# Surface ssflux and ssc
# In the grid's u and v directions
# Plot bathymetry
# Plot surface ssc and flux
ax5[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs16 = ax5[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_surf_allsed_avg_wland_masked_trimmed, lev6, cmap=cmap5, extend='max')
# Plot bathymetry contours
ax5[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='bisque')
# Plot currents
q4 = ax5[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, ssflux_u_surf_avg_wland_masked_trimmed_slice, 
                   ssflux_v_surf_avg_wland_masked_trimmed_slice, color='teal', width=2, 
                   angles='xy', scale_units='xy', units='xy')
ax5[0].quiverkey(q4, 0.65, 0.85, U=0.005, label='0.005 kg/m\u00b2s', fontproperties={'size':fontsize-2})

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax5[0].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax5[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m flux and ssc
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax5[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs17 = ax5[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_1m_avg_wland_masked_trimmed, lev6, cmap=cmap5, extend='max')
# Plot bathymetry contours 
ax5[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='bisque')
# Plot currents
q5 = ax5[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_u_1m_avg_wland_masked_trimmed_slice, ssflux_v_1m_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax5[1].quiverkey(q5, 0.65, 0.85, U=0.001, label='0.001 kg/m\u00b2s', fontproperties={'size':fontsize-2})

# Label the plot
plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax5[1].set_ylabel('Y (km)', fontsize=fontsize)

# Depth-integrated
# In the grid's u and v directions
# Plot depth-avg salinity 
ax5[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs18 = ax5[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_depth_avg_allsed_avg_trimmed, lev6, cmap=cmap5, extend='max')
# Plot bathymetry contours 
ax5[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='bisque')
# Plot currents
q6 = ax5[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_depth_int_u_avg_wland_masked_trimmed_slice, ssflux_depth_int_v_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax5[2].quiverkey(q6, 0.65, 0.85, U=0.02, label='0.02 kg/ms', fontproperties={'size':fontsize-2})
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax5[2].get_xticklabels(), visible=False)
#ax5[2].set_xlabel('X (km)', fontsize=fontsize)
ax5[2].set_ylabel('Y (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig5.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)

#axes = ax1.ravel().tolist()
cbar5_ax = fig5.add_axes([0.82, 0.32, 0.02, 0.63]) #[left, bottom, width, height]
cbar5 = plt.colorbar(cs18, ax=[ax5[0], ax5[1], ax5[2]], cax=cbar5_ax, orientation='vertical', pad=0.015).set_label(label='SSC (kg/m\u00b3)', size=fontsize, labelpad=15)

# Plot net erosion & deposition
ax5[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[3].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs19 = ax5[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  roms_net_erodepo_bedthick_wland_masked_trimmed*100, lev7, cmap=cmap6, extend='max')
ax5[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgray')
# Colville River
eta_col_idx = 39
xi_col_idx = 166-36
s1 = ax5[3].scatter(x_rho_flat_trimmed[xi_col_idx-2]/1000, y_rho_flat[eta_col_idx-9]/1000, 
            marker='d',  s=900, linewidth=4, color='brown', label='Colville')
cbar6_ax = fig5.add_axes([0.82, 0.1, 0.02, 0.2]) #[left, bottom, width, height]
cbar6 = plt.colorbar(cs19, cax=cbar6_ax, ax=ax5[3], orientation='vertical', pad=0.015).set_label('Net Erosion \n& Deposition (cm)', fontsize=fontsize, labelpad=15)

# Label the plot
#ax13.set_title('Time-Averaged Depth-Integrated SSC Flux for All Seds Masked', fontsize=fontsize)
ax5[3].set_xlabel('X (km)', fontsize=fontsize)
ax5[3].set_ylabel('Y (km)', fontsize=fontsize)

ax5[3].legend(fontsize=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
# =============================================================================
# plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.246, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# =============================================================================
# Add labels
# Bottom right 
plt.text(0.773, 0.760, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.542, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.329, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.113, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# ----------------------------------------------------------------------------------------------
# -------------------------- Plot 6: Sediment Flux, SSC, and -----------------------------------
# ----------------------- Net Erosion & Deposition Masked XY (slightly different) --------------
# ----------------------------------------------------------------------------------------------
# Plot the ssc flux, ssc, and net erosiion deposition but masking the 
# regions we don't trust and in xy

# Make the figure
fig6, ax6 = plt.subplots(4, figsize=(22,21)) # (18,8) (26,12) (26,8) (26,10)

# Set the colormaps
cmap5 = cmocean.cm.turbid
cmap6 = cmocean.cm.delta # ero depo

# Set colorbar levels for all plots currents
lev5 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots salinity 
#lev2 = np.arange(0,0.06,0.001)
lev6 = np.arange(0,0.05,0.001) # ssc
lev7 = np.arange(-1.6, 1.6, 0.01) # np.arange(-1, 1, 0.01) # ero depo

# Depth-integrated
# In the grid's u and v directions
# Plot depth-avg salinity 
ax6[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs18 = ax6[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_depth_avg_allsed_avg_trimmed, lev6, cmap=cmap5, extend='max')
# Plot bathymetry contours 
ax6[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgrey') # bisque
# Plot currents
q6 = ax6[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_depth_int_u_avg_wland_masked_trimmed_slice, ssflux_depth_int_v_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax6[0].quiverkey(q6, 0.65, 0.85, U=0.02, label='0.02 kg $m^{-1}$ $s^{-1}$', fontproperties={'size':fontsize-2})
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax6[0].get_xticklabels(), visible=False)
#ax5[2].set_xlabel('X (km)', fontsize=fontsize)
ax6[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)


# Surface ssflux and ssc
# In the grid's u and v directions
# Plot bathymetry
# Plot surface ssc and flux
ax6[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs16 = ax6[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_surf_allsed_avg_wland_masked_trimmed, lev6, cmap=cmap5, extend='max')
# Plot bathymetry contours
ax6[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgrey') # bisque 
# Plot currents
q4 = ax6[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, ssflux_u_surf_avg_wland_masked_trimmed_slice, 
                   ssflux_v_surf_avg_wland_masked_trimmed_slice, color='teal', width=2, 
                   angles='xy', scale_units='xy', units='xy')
ax6[1].quiverkey(q4, 0.65, 0.85, U=0.005, label='0.005 kg $m^{-2}$ $s^{-1}$', fontproperties={'size':fontsize-2})

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax6[1].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax6[1].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m flux and ssc
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax6[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs17 = ax6[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_1m_avg_wland_masked_trimmed, lev6, cmap=cmap5, extend='max')
# Plot bathymetry contours 
ax6[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgrey') # bisque
# Plot currents
q5 = ax6[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_u_1m_avg_wland_masked_trimmed_slice, ssflux_v_1m_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax6[2].quiverkey(q5, 0.65, 0.85, U=0.001, label='0.001 kg $m^{-2}$ $s^{-1}$', fontproperties={'size':fontsize-2})

# Label the plot
plt.setp(ax6[2].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax6[2].set_ylabel('Y (km)', fontsize=fontsize)

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig6.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
#axes = ax1.ravel().tolist()
cbar7_ax = fig6.add_axes([0.82, 0.32, 0.015, 0.63]) #[left, bottom, width, height]
cbar7 = plt.colorbar(cs18, ax=[ax6[0], ax6[1], ax6[2]], cax=cbar7_ax, orientation='vertical', pad=0.015).set_label(label='SSC \n(kg $m^{-3}$)', size=fontsize, labelpad=55, rotation='horizontal')

# Plot net erosion & deposition
ax6[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[3].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs19 = ax6[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  roms_net_erodepo_bedthick_wland_masked_trimmed*100, lev7, cmap=cmap6, extend='max')
ax6[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgray')
# Colville River
eta_col_idx = 39
xi_col_idx = 166-36
s1 = ax6[3].scatter(x_rho_flat_trimmed[xi_col_idx-2]/1000, y_rho_flat[eta_col_idx-9]/1000, 
            marker='d',  s=900, linewidth=4, color='brown', label='Colville')
cbar8_ax = fig6.add_axes([0.82, 0.1, 0.015, 0.2]) #[left, bottom, width, height]
cbar8 = plt.colorbar(cs19, cax=cbar8_ax, ax=ax6[3], orientation='vertical', pad=0.015).set_label('Net \nErosion & \nDeposition \n(cm)', fontsize=fontsize, labelpad=85, rotation='horizontal')

# Label the plot
#ax13.set_title('Time-Averaged Depth-Integrated SSC Flux for All Seds Masked', fontsize=fontsize)
ax6[3].set_xlabel('X (km)', fontsize=fontsize)
ax6[3].set_ylabel('Y (km)', fontsize=fontsize)

ax6[3].legend(fontsize=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
# =============================================================================
# plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.246, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# =============================================================================
# Add labels
# Bottom right 
plt.text(0.773, 0.760, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.542, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.329, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.113, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# ----------------------------------------------------------------------------------------------
# -------------------------- Plot 7: Sediment Flux, SSC, and -----------------------------------
# ----------------------- Net Erosion & Deposition Masked XY (slightly different) --------------
# ----------------------------------------------------------------------------------------------
# Make the same plot as above but with some slight edits
# (mg/L SSC, "net deposition" on colorbar, and add in river locations and legend 
# to all panels)

# Set the river indices
# Go from West to East
# Kalikpik River
eta_kal_idx = 23 #22
xi_kal_idx = 87 - c_west

# Fish Creek
eta_fis_idx = 20
xi_fis_idx = 117 - c_west #116

# Colville River
eta_col_idx = 39
xi_col_idx = 166 - c_west #166

# Sakonowyak River
eta_sak_idx = 46 #45
xi_sak_idx = 234 - c_west 

# Kuparik
# Kukpuk
eta_kuk_idx = 41 #40
xi_kuk_idx = 239 - c_west 

# Kuparuk
eta_kup_idx = 41 #40
xi_kup_idx = 242 - c_west

# Fawn Creek
#eta_faw_idx = 44 #43
#xi_faw_idx = 249

# Putuligayuk River
eta_put_idx = 28 #27
xi_put_idx = 264 - c_west 

# Sagavanirktok River
eta_sag_idx = 37 #36
xi_sag_idx = 279 - c_west 

# Canning River
# Staines River
eta_sta_idx = 27 #26
xi_sta_idx = 393 - c_west 

# Canning River
eta_can_idx = 20 #19
xi_can_idx = 416 - c_west 

# Katakturuk River
eta_kat_idx = 9 #8
xi_kat_idx = 447 - c_west

# Hulahula River
eta_hul_idx = 40
xi_hul_idx = 489 - c_west 

# Jago River
eta_jag_idx = 62 #61
xi_jag_idx = 528 - c_west 

# Siksik River
eta_sik_idx = 46
xi_sik_idx = 574 - c_west  #573


# Make the figure
fig7, ax7 = plt.subplots(4, figsize=(22,21)) # (18,8) (26,12) (26,8) (26,10)

# Depth-integrated
# In the grid's u and v directions
# Plot depth-avg salinity 
ax7[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
# kg/m3
cs20 = ax7[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_depth_avg_allsed_avg_trimmed, lev6, cmap=cmap5, extend='max')
# mg/L
#cs20 = ax7[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
 #                 ssc_depth_avg_allsed_avg_trimmed*1000, lev6*1000, cmap=cmap5, extend='max') 
# Plot bathymetry contours 
ax7[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgrey') # bisque
# Plot currents
q7 = ax7[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_depth_int_u_avg_wland_masked_trimmed_slice, ssflux_depth_int_v_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax7[0].quiverkey(q7, 0.65, 0.85, U=0.02, label='0.02 kg $m^{-1}$ $s^{-1}$', fontproperties={'size':fontsize-2})
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax7[0].get_xticklabels(), visible=False)
#ax5[2].set_xlabel('X (km)', fontsize=fontsize)
ax7[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Plot the river mouths
# Kalikpik River
s1 = ax7[0].scatter(x_rho_flat_trimmed[xi_kal_idx]/1000, y_rho_flat[eta_kal_idx]/1000, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
s2 = ax7[0].scatter(x_rho_flat_trimmed[xi_fis_idx]/1000, y_rho_flat[eta_fis_idx]/1000, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
s3 = ax7[0].scatter(x_rho_flat_trimmed[xi_col_idx]/1000, y_rho_flat[eta_col_idx]/1000, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
s4 = ax7[0].scatter(x_rho_flat_trimmed[xi_sak_idx]/1000, y_rho_flat[eta_sak_idx]/1000, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s5 = ax7[0].scatter(x_rho_flat_trimmed[xi_kuk_idx]/1000, y_rho_flat[eta_kuk_idx]/1000, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax7[0].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
s8 = ax7[0].scatter(x_rho_flat_trimmed[xi_put_idx]/1000, y_rho_flat[eta_put_idx]/1000, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
s9 = ax7[0].scatter(x_rho_flat_trimmed[xi_sag_idx]/1000, y_rho_flat[eta_sag_idx]/1000, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
s10 = ax7[0].scatter(x_rho_flat_trimmed[xi_sta_idx]/1000, y_rho_flat[eta_sta_idx]/1000, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
s11 = ax7[0].scatter(x_rho_flat_trimmed[xi_can_idx]/1000, y_rho_flat[eta_can_idx]/1000, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
s12 = ax7[0].scatter(x_rho_flat_trimmed[xi_kat_idx]/1000, y_rho_flat[eta_kat_idx]/1000, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
s13 = ax7[0].scatter(x_rho_flat_trimmed[xi_hul_idx]/1000, y_rho_flat[eta_hul_idx]/1000, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
s14 = ax7[0].scatter(x_rho_flat_trimmed[xi_jag_idx]/1000, y_rho_flat[eta_jag_idx]/1000, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
#s15 = ax7[0].scatter(x_rho_flat_trimmed[xi_sik_idx]/1000, y_rho_flat[eta_sik_idx]/1000, 
#            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')



# Surface ssflux and ssc
# In the grid's u and v directions
# Plot bathymetry
# Plot surface ssc and flux
ax7[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
# kg/m3
cs21 = ax7[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_surf_allsed_avg_wland_masked_trimmed, lev6, cmap=cmap5, extend='max')
# mg/L
#cs21 = ax7[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
 #                 ssc_surf_allsed_avg_wland_masked_trimmed*1000, lev6*1000, cmap=cmap5, extend='max')
# Plot bathymetry contours
ax7[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgrey') # bisque 
# Plot currents
q8 = ax7[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, ssflux_u_surf_avg_wland_masked_trimmed_slice, 
                   ssflux_v_surf_avg_wland_masked_trimmed_slice, color='teal', width=2, 
                   angles='xy', scale_units='xy', units='xy')
ax7[1].quiverkey(q8, 0.65, 0.85, U=0.005, label='0.005 kg $m^{-2}$ $s^{-1}$', fontproperties={'size':fontsize-2})

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax7[1].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax7[1].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Plot the river mouths
# Kalikpik River
s1 = ax7[1].scatter(x_rho_flat_trimmed[xi_kal_idx]/1000, y_rho_flat[eta_kal_idx]/1000, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
s2 = ax7[1].scatter(x_rho_flat_trimmed[xi_fis_idx]/1000, y_rho_flat[eta_fis_idx]/1000, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
s3 = ax7[1].scatter(x_rho_flat_trimmed[xi_col_idx]/1000, y_rho_flat[eta_col_idx]/1000, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
s4 = ax7[1].scatter(x_rho_flat_trimmed[xi_sak_idx]/1000, y_rho_flat[eta_sak_idx]/1000, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s5 = ax7[1].scatter(x_rho_flat_trimmed[xi_kuk_idx]/1000, y_rho_flat[eta_kuk_idx]/1000, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax7[1].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
s8 = ax7[1].scatter(x_rho_flat_trimmed[xi_put_idx]/1000, y_rho_flat[eta_put_idx]/1000, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
s9 = ax7[1].scatter(x_rho_flat_trimmed[xi_sag_idx]/1000, y_rho_flat[eta_sag_idx]/1000, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
s10 = ax7[1].scatter(x_rho_flat_trimmed[xi_sta_idx]/1000, y_rho_flat[eta_sta_idx]/1000, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
s11 = ax7[1].scatter(x_rho_flat_trimmed[xi_can_idx]/1000, y_rho_flat[eta_can_idx]/1000, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
s12 = ax7[1].scatter(x_rho_flat_trimmed[xi_kat_idx]/1000, y_rho_flat[eta_kat_idx]/1000, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
s13 = ax7[1].scatter(x_rho_flat_trimmed[xi_hul_idx]/1000, y_rho_flat[eta_hul_idx]/1000, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
s14 = ax7[1].scatter(x_rho_flat_trimmed[xi_jag_idx]/1000, y_rho_flat[eta_jag_idx]/1000, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
#s15 = ax7[1].scatter(x_rho_flat_trimmed[xi_sik_idx]/1000, y_rho_flat[eta_sik_idx]/1000, 
#            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')


# Bottom  1m flux and ssc
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax7[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
# kg/m3
cs22 = ax7[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  ssc_1m_avg_wland_masked_trimmed, lev6, cmap=cmap5, extend='max')
# mg/L
#s22 = ax7[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
 #                 ssc_1m_avg_wland_masked_trimmed*1000, lev6*1000, cmap=cmap5, extend='max')
# Plot bathymetry contours 
ax7[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgrey') # bisque
# Plot currents
q9 = ax7[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ssflux_u_1m_avg_wland_masked_trimmed_slice, ssflux_v_1m_avg_wland_masked_trimmed_slice, 
                   color='teal', width=2,
                   angles='xy', scale_units='xy', units='xy')
ax7[2].quiverkey(q9, 0.65, 0.85, U=0.001, label='0.001 kg $m^{-2}$ $s^{-1}$', fontproperties={'size':fontsize-2})

# Label the plot
plt.setp(ax7[2].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax7[2].set_ylabel('Y (km)', fontsize=fontsize)

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig7.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
#axes = ax1.ravel().tolist()
cbar7_ax = fig7.add_axes([0.82, 0.32, 0.015, 0.63]) #[left, bottom, width, height]
# kg/m3
cbar7 = plt.colorbar(cs20, ax=[ax7[0], ax7[1], ax7[2]], cax=cbar7_ax, orientation='vertical', pad=0.015).set_label(label='SSC \n(kg $m^{-3}$)', size=fontsize, labelpad=55, rotation='horizontal')
# mg/L
#cbar7 = plt.colorbar(cs20, ax=[ax7[0], ax7[1], ax7[2]], cax=cbar7_ax, orientation='vertical', pad=0.015).set_label(label='SSC \n(mg $L^{-1}$)', size=fontsize, labelpad=55, rotation='horizontal')

# Plot the river mouths
# Kalikpik River
s1 = ax7[2].scatter(x_rho_flat_trimmed[xi_kal_idx]/1000, y_rho_flat[eta_kal_idx]/1000, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
s2 = ax7[2].scatter(x_rho_flat_trimmed[xi_fis_idx]/1000, y_rho_flat[eta_fis_idx]/1000, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
s3 = ax7[2].scatter(x_rho_flat_trimmed[xi_col_idx]/1000, y_rho_flat[eta_col_idx]/1000, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
s4 = ax7[2].scatter(x_rho_flat_trimmed[xi_sak_idx]/1000, y_rho_flat[eta_sak_idx]/1000, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s5 = ax7[2].scatter(x_rho_flat_trimmed[xi_kuk_idx]/1000, y_rho_flat[eta_kuk_idx]/1000, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax7[2].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
s8 = ax7[2].scatter(x_rho_flat_trimmed[xi_put_idx]/1000, y_rho_flat[eta_put_idx]/1000, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
s9 = ax7[2].scatter(x_rho_flat_trimmed[xi_sag_idx]/1000, y_rho_flat[eta_sag_idx]/1000, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
s10 = ax7[2].scatter(x_rho_flat_trimmed[xi_sta_idx]/1000, y_rho_flat[eta_sta_idx]/1000, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
s11 = ax7[2].scatter(x_rho_flat_trimmed[xi_can_idx]/1000, y_rho_flat[eta_can_idx]/1000, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
s12 = ax7[2].scatter(x_rho_flat_trimmed[xi_kat_idx]/1000, y_rho_flat[eta_kat_idx]/1000, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
s13 = ax7[2].scatter(x_rho_flat_trimmed[xi_hul_idx]/1000, y_rho_flat[eta_hul_idx]/1000, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
s14 = ax7[2].scatter(x_rho_flat_trimmed[xi_jag_idx]/1000, y_rho_flat[eta_jag_idx]/1000, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
#s15 = ax7[2].scatter(x_rho_flat_trimmed[xi_sik_idx]/1000, y_rho_flat[eta_sik_idx]/1000, 
#            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')


# Plot net erosion & deposition
ax7[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[3].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs23 = ax7[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  roms_net_erodepo_bedthick_wland_masked_trimmed*100, lev7, cmap=cmap6, extend='both')
ax7[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='dimgray')
# Colorbar
cbar9_ax = fig7.add_axes([0.82, 0.1, 0.015, 0.2]) #[left, bottom, width, height]
cbar9 = plt.colorbar(cs23, cax=cbar9_ax, ax=ax7[3], orientation='vertical', pad=0.015).set_label('Net \nDeposition \n(cm)', fontsize=fontsize, labelpad=85, rotation='horizontal')
#cbar9 = plt.colorbar(cs22, cax=cbar9_ax, ax=ax6[3], orientation='vertical', pad=0.015).set_label('Net \nErosion & \nDeposition \n(cm)', fontsize=fontsize, labelpad=85, rotation='horizontal')

# Label the plot
#ax13.set_title('Time-Averaged Depth-Integrated SSC Flux for All Seds Masked', fontsize=fontsize)
ax7[3].set_xlabel('X (km)', fontsize=fontsize)
ax7[3].set_ylabel('Y (km)', fontsize=fontsize)

# Plot the river mouths
# Kalikpik River
s1 = ax7[3].scatter(x_rho_flat_trimmed[xi_kal_idx]/1000, y_rho_flat[eta_kal_idx]/1000, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
s2 = ax7[3].scatter(x_rho_flat_trimmed[xi_fis_idx]/1000, y_rho_flat[eta_fis_idx]/1000, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
s3 = ax7[3].scatter(x_rho_flat_trimmed[xi_col_idx]/1000, y_rho_flat[eta_col_idx]/1000, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
s4 = ax7[3].scatter(x_rho_flat_trimmed[xi_sak_idx]/1000, y_rho_flat[eta_sak_idx]/1000, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s5 = ax7[3].scatter(x_rho_flat_trimmed[xi_kuk_idx]/1000, y_rho_flat[eta_kuk_idx]/1000, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax7[3].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
s8 = ax7[3].scatter(x_rho_flat_trimmed[xi_put_idx]/1000, y_rho_flat[eta_put_idx]/1000, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
s9 = ax7[3].scatter(x_rho_flat_trimmed[xi_sag_idx]/1000, y_rho_flat[eta_sag_idx]/1000, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
s10 = ax7[3].scatter(x_rho_flat_trimmed[xi_sta_idx]/1000, y_rho_flat[eta_sta_idx]/1000, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
s11 = ax7[3].scatter(x_rho_flat_trimmed[xi_can_idx]/1000, y_rho_flat[eta_can_idx]/1000, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
s12 = ax7[3].scatter(x_rho_flat_trimmed[xi_kat_idx]/1000, y_rho_flat[eta_kat_idx]/1000, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
s13 = ax7[3].scatter(x_rho_flat_trimmed[xi_hul_idx]/1000, y_rho_flat[eta_hul_idx]/1000, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
s14 = ax7[3].scatter(x_rho_flat_trimmed[xi_jag_idx]/1000, y_rho_flat[eta_jag_idx]/1000, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
#s15 = ax7[3].scatter(x_rho_flat_trimmed[xi_sik_idx]/1000, y_rho_flat[eta_sik_idx]/1000, 
#            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Put a legend for the rivers
ax7[3].legend(fontsize=fontsize-2, loc='center left', ncol=4, columnspacing=0.1, 
              labelspacing=0.1,  bbox_to_anchor=(0.1, -0.45))

# Add subplot labels
# =============================================================================
# plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# plt.text(0.775, 0.246, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# =============================================================================
# Add labels
# Bottom right 
plt.text(0.773, 0.760, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.542, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.329, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.113, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)












# Calculate and print a bunch of stats 
# Make a function to mask the data
def masked_array(data, threshold):
    """
    This function takes an array and masks all values that are less
    than a certain given threshold. The functions returns 1 for areas that meet 
    the condition and 0 for areas that don't. So areas where the array is less
    than the threshold get returned as 1 and areas greater than the threshold
    are returned as 0. This function maintains the shape of the array.
    
    """
    return (data <= threshold).astype(int)

# Make a function to mask the data but that takes two thresholds
def masked_array_lowhigh_2dloop(data, lower, upper):
    """
    This function takes an array and masks all values that are less
    than a certain given threshold. The functions returns 1 for areas that meet 
    the condition and 0 for areas that don't. So areas where the array is less
    than the threshold get returned as 1 and areas greater than the threshold
    are returned as 0. This function maintains the shape of the array.
    
    """
    mask_tmp = np.empty_like((data))
    
    # Loop through dimension 1
    for i in range(len(data[:,0])):
        # Loop through dimension 2
        for j in range(len(data[0,:])):
            # Compare against threshold 
            value = data[i,j]
            if lower < value <= upper:
                mask_tmp[i,j] = 1
            else:
                mask_tmp[i,j] = 0
    
    
    return (mask_tmp).astype(int)
 
# Call the function to make the mask
# First make an array of the bathymetry of the grid pre-masked
h_masked = grid.h * mask_rho_nan.nudge_mask_rho_nan
# Inner shelf
h_masked1 = h_masked.copy()
#inner_shelf_mask_rho = masked_array(h_masked1, 20)
inner_shelf_mask_rho = masked_array_lowhigh_2dloop(h_masked1, 2, 20)
# Mid shelf
h_masked2 = h_masked.copy()
mid_shelf_mask_rho = masked_array_lowhigh_2dloop(h_masked2, 20, 40)
# Outer shelf 
h_masked3 = h_masked.copy()
outer_shelf_mask_rho = masked_array_lowhigh_2dloop(h_masked3, 40, 60)
# 0 - 10 m depth
h_masked4 = h_masked.copy()
#inner_10m_mask_rho = masked_array(h_masked4, 10)
inner_10m_mask_rho = masked_array_lowhigh_2dloop(h_masked4, 2, 10)
# 10 - 60 m depth
h_masked5 = h_masked.copy()
outer_10_60m_mask_rho = masked_array_lowhigh_2dloop(h_masked5, 10, 60)

# Make the masks nan where they are 0 so that these out of bounds areas are 
# nanned out 
# Inner Shelf
inner_shelf_mask_rho_nan_idx = np.where(inner_shelf_mask_rho == 0.0)
inner_shelf_mask_rho_nan = inner_shelf_mask_rho.copy()
inner_shelf_mask_rho_nan = inner_shelf_mask_rho_nan.astype('float')
inner_shelf_mask_rho_nan[inner_shelf_mask_rho_nan_idx] = np.nan
# Mid Shelf
mid_shelf_mask_rho_nan_idx = np.where(mid_shelf_mask_rho == 0.0)
mid_shelf_mask_rho_nan = mid_shelf_mask_rho.copy()
mid_shelf_mask_rho_nan = mid_shelf_mask_rho_nan.astype('float')
mid_shelf_mask_rho_nan[mid_shelf_mask_rho_nan_idx] = np.nan
# Outer Shelf 
outer_shelf_mask_rho_nan_idx = np.where(outer_shelf_mask_rho == 0.0)
outer_shelf_mask_rho_nan = outer_shelf_mask_rho.copy()
outer_shelf_mask_rho_nan = outer_shelf_mask_rho_nan.astype('float')
outer_shelf_mask_rho_nan[outer_shelf_mask_rho_nan_idx] = np.nan
# 0 - 10 m depth
inner_10m_mask_rho_nan_idx = np.where(inner_10m_mask_rho == 0.0)
inner_10m_mask_rho_nan = inner_10m_mask_rho.copy()
inner_10m_mask_rho_nan = inner_10m_mask_rho_nan.astype('float')
inner_10m_mask_rho_nan[inner_10m_mask_rho_nan_idx] = np.nan
# 10 - 60 m depth
outer_10_60m_mask_rho_nan_idx = np.where(outer_10_60m_mask_rho == 0.0)
outer_10_60m_mask_rho_nan = outer_10_60m_mask_rho.copy()
outer_10_60m_mask_rho_nan = outer_10_60m_mask_rho_nan.astype('float')
outer_10_60m_mask_rho_nan[outer_10_60m_mask_rho_nan_idx] = np.nan


# Now multiply by the mask to get the different regions 
# For depth-averaged ssc
# Inner 
ssc_depth_avg_allsed_avg_masked_inner = ssc_depth_avg_allsed_avg * inner_shelf_mask_rho_nan
# Mid
ssc_depth_avg_allsed_avg_masked_mid = ssc_depth_avg_allsed_avg * mid_shelf_mask_rho_nan
# Outer
ssc_depth_avg_allsed_avg_masked_outer = ssc_depth_avg_allsed_avg * outer_shelf_mask_rho_nan
# 0 - 10 m
ssc_depth_avg_allsed_avg_masked_10m = ssc_depth_avg_allsed_avg * inner_10m_mask_rho_nan
# 10 - 60 m
ssc_depth_avg_allsed_avg_masked_10_60m = ssc_depth_avg_allsed_avg * outer_10_60m_mask_rho_nan

# For surface SSC 
# Inner 
ssc_surf_allsed_avg_masked_inner = ssc_surf_allsed_avg * inner_shelf_mask_rho_nan
# Mid
ssc_surf_allsed_avg_masked_mid = ssc_surf_allsed_avg * mid_shelf_mask_rho_nan
# Outer
ssc_surf_allsed_avg_masked_outer = ssc_surf_allsed_avg * outer_shelf_mask_rho_nan
# 0 - 10 m
ssc_surf_allsed_avg_masked_10m = ssc_surf_allsed_avg * inner_10m_mask_rho_nan
# 10 - 60 m
ssc_surf_allsed_avg_masked_10_60m = ssc_surf_allsed_avg * outer_10_60m_mask_rho_nan

# For depth-integrated ssflux
# 0 - 10 m eastern
# u
ssflux_depth_int_u_avg_masked_10m = ssflux_depth_int_u_avg * inner_10m_mask_rho
ssflux_depth_int_u_avg_masked_10m_east = ssflux_depth_int_u_avg_masked_10m[:,304:]
# v 
ssflux_depth_int_v_avg_masked_10m = ssflux_depth_int_v_avg * inner_10m_mask_rho
ssflux_depth_int_v_avg_masked_10m_east = ssflux_depth_int_v_avg_masked_10m[:,304:]
# 0 - 10 m western
# u
ssflux_depth_int_u_avg_masked_10m_west = ssflux_depth_int_u_avg_masked_10m[:,:304]
# v 
ssflux_depth_int_v_avg_masked_10m_west = ssflux_depth_int_v_avg_masked_10m[:,:304]
# Inner 
ssflux_depth_int_u_avg_masked_inner = ssflux_depth_int_u_avg * inner_shelf_mask_rho_nan
ssflux_depth_int_v_avg_masked_inner = ssflux_depth_int_v_avg * inner_shelf_mask_rho_nan
# Mid
ssflux_depth_int_u_avg_masked_mid = ssflux_depth_int_u_avg * mid_shelf_mask_rho_nan
ssflux_depth_int_v_avg_masked_mid = ssflux_depth_int_v_avg * mid_shelf_mask_rho_nan
# Outer
ssflux_depth_int_u_avg_masked_outer = ssflux_depth_int_u_avg * outer_shelf_mask_rho_nan
ssflux_depth_int_v_avg_masked_outer = ssflux_depth_int_v_avg * outer_shelf_mask_rho_nan
# 0 - 10 m
ssflux_depth_int_u_avg_masked_10m = ssflux_depth_int_u_avg * inner_10m_mask_rho_nan
ssflux_depth_int_v_avg_masked_10m = ssflux_depth_int_v_avg * inner_10m_mask_rho_nan
# 10 - 60 m
ssflux_depth_int_u_avg_masked_10_60m = ssflux_depth_int_u_avg * outer_10_60m_mask_rho_nan
ssflux_depth_int_v_avg_masked_10_60m = ssflux_depth_int_v_avg * outer_10_60m_mask_rho_nan
# Surface ssflux
# Inner 
ssflux_u_surf_avg_masked_inner = ssflux_u_surf_avg * inner_shelf_mask_rho_nan
ssflux_v_surf_avg_masked_inner = ssflux_v_surf_avg * inner_shelf_mask_rho_nan

# Mask and trim these so that we are just looking at the regions in the plot
# Mask, trim, slice
# Mask
# SSC
# Depth-averaged SSC
# Inner
ssc_depth_avg_allsed_avg_masked_inner_masked = ssc_depth_avg_allsed_avg_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Mid
ssc_depth_avg_allsed_avg_masked_mid_masked = ssc_depth_avg_allsed_avg_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Outer
ssc_depth_avg_allsed_avg_masked_outer_masked = ssc_depth_avg_allsed_avg_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 0 - 10 m
ssc_depth_avg_allsed_avg_masked_10m_masked = ssc_depth_avg_allsed_avg_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 60 m
ssc_depth_avg_allsed_avg_masked_10_60m_masked = ssc_depth_avg_allsed_avg_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Surface SSC
# Inner
ssc_surf_allsed_avg_masked_inner_masked = ssc_surf_allsed_avg_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Mid
ssc_surf_allsed_avg_masked_mid_masked = ssc_surf_allsed_avg_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Outer
ssc_surf_allsed_avg_masked_outer_masked = ssc_surf_allsed_avg_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 0 - 10 m
ssc_surf_allsed_avg_masked_10m_masked = ssc_surf_allsed_avg_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 60 m
ssc_surf_allsed_avg_masked_10_60m_masked = ssc_surf_allsed_avg_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask

# SSflux - depth-integrated 
# Inner
ssflux_depth_int_u_avg_masked_inner_masked = ssflux_depth_int_u_avg_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_depth_int_v_avg_masked_inner_masked = ssflux_depth_int_v_avg_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Mid
ssflux_depth_int_u_avg_masked_mid_masked = ssflux_depth_int_u_avg_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_depth_int_v_avg_masked_mid_masked = ssflux_depth_int_v_avg_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Outer
ssflux_depth_int_u_avg_masked_outer_masked = ssflux_depth_int_u_avg_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_depth_int_v_avg_masked_outer_masked = ssflux_depth_int_v_avg_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 0 - 10 m
ssflux_depth_int_u_avg_masked_10m_masked = ssflux_depth_int_u_avg_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_depth_int_v_avg_masked_10m_masked = ssflux_depth_int_v_avg_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 60 m
ssflux_depth_int_u_avg_masked_10_60m_masked = ssflux_depth_int_u_avg_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_depth_int_v_avg_masked_10_60m_masked = ssflux_depth_int_v_avg_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Surface ssflux
# Inner
ssflux_u_surf_avg_masked_inner_masked = ssflux_u_surf_avg_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
ssflux_v_surf_avg_masked_inner_masked = ssflux_v_surf_avg_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask


# Trim
# SSC
# Depth-averaged SSC
# Inner
ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed = ssc_depth_avg_allsed_avg_masked_inner_masked[:,c_west:-c_west]
# Mid
ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed = ssc_depth_avg_allsed_avg_masked_mid_masked[:,c_west:-c_west]
# Outer
ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed = ssc_depth_avg_allsed_avg_masked_outer_masked[:,c_west:-c_west]
# 0 - 10 m
ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed = ssc_depth_avg_allsed_avg_masked_10m_masked[:,c_west:-c_west]
# 10 - 60 m
ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed = ssc_depth_avg_allsed_avg_masked_10_60m_masked[:,c_west:-c_west]
# Surface SSC
# Inner
ssc_surf_allsed_avg_masked_inner_masked_trimmed = ssc_surf_allsed_avg_masked_inner_masked[:,c_west:-c_west]
# Mid
ssc_surf_allsed_avg_masked_mid_masked_trimmed = ssc_surf_allsed_avg_masked_mid_masked[:,c_west:-c_west]
# Outer 
ssc_surf_allsed_avg_masked_outer_masked_trimmed = ssc_surf_allsed_avg_masked_outer_masked[:,c_west:-c_west]
# 0 - 10 m
ssc_surf_allsed_avg_masked_10m_masked_trimmed = ssc_surf_allsed_avg_masked_10m_masked[:,c_west:-c_west]
# 10 - 60 m
ssc_surf_allsed_avg_masked_10_60m_masked_trimmed = ssc_surf_allsed_avg_masked_10_60m_masked[:,c_west:-c_west]

# SSflux
# Inner
ssflux_depth_int_u_avg_masked_inner_masked_trimmed = ssflux_depth_int_u_avg_masked_inner_masked[:,c_west:-c_west]
ssflux_depth_int_v_avg_masked_inner_masked_trimmed = ssflux_depth_int_v_avg_masked_inner_masked[:,c_west:-c_west]
# Mid
ssflux_depth_int_u_avg_masked_mid_masked_trimmed = ssflux_depth_int_u_avg_masked_mid_masked[:,c_west:-c_west]
ssflux_depth_int_v_avg_masked_mid_masked_trimmed = ssflux_depth_int_v_avg_masked_mid_masked[:,c_west:-c_west]
# Outer 
ssflux_depth_int_u_avg_masked_outer_masked_trimmed = ssflux_depth_int_u_avg_masked_outer_masked[:,c_west:-c_west]
ssflux_depth_int_v_avg_masked_outer_masked_trimmed = ssflux_depth_int_v_avg_masked_outer_masked[:,c_west:-c_west]
# 0 - 10 m
ssflux_depth_int_u_avg_masked_10m_masked_trimmed = ssflux_depth_int_u_avg_masked_10m_masked[:,c_west:-c_west]
ssflux_depth_int_v_avg_masked_10m_masked_trimmed = ssflux_depth_int_v_avg_masked_10m_masked[:,c_west:-c_west]
# 10 - 60 m
ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed = ssflux_depth_int_u_avg_masked_10_60m_masked[:,c_west:-c_west]
ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed = ssflux_depth_int_v_avg_masked_10_60m_masked[:,c_west:-c_west]
# Surface ssflux
# Inner
ssflux_u_surf_avg_masked_inner_masked_trimmed = ssflux_u_surf_avg_masked_inner_masked[:,c_west:-c_west]
ssflux_v_surf_avg_masked_inner_masked_trimmed = ssflux_v_surf_avg_masked_inner_masked[:,c_west:-c_west]



# Print some statistics 
# SSC
# Depth-averaged SSC
# Mean 
ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_avg = np.nanmean(ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_avg)
ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_avg = np.nanmean(ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_avg)
ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_avg = np.nanmean(ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_avg)
ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_avg = np.nanmean(ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_avg)
ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_avg = np.nanmean(ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_avg)

# Standard deviation 
ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_std = np.nanstd(ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) std ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_std)
ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_std = np.nanstd(ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) std ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_std)
ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_std = np.nanstd(ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) std ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_std)
ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_std = np.nanstd(ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m std total ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_std)
ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_std = np.nanstd(ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m std ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_std)

# Min
ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_min = np.nanmin(ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) min ssc_depth_avg_allsed_avg(kg/m3): ', ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_min)
ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_min = np.nanmin(ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) min ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_min)
ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_min = np.nanmin(ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) min ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_min)
ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_min = np.nanmin(ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m min ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_min)
ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_min = np.nanmin(ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m min ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_min)

# Max
ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_max = np.nanmax(ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) max ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_inner_masked_trimmed_max)
ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_max = np.nanmax(ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) max ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_mid_masked_trimmed_max)
ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_max = np.nanmax(ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) max ssc_depth_avg_allsed_avg (kg/m3: ', ssc_depth_avg_allsed_avg_masked_outer_masked_trimmed_max)
ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_max = np.nanmax(ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m max tssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10m_masked_trimmed_max)
ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_max = np.nanmax(ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m max ssc_depth_avg_allsed_avg (kg/m3): ', ssc_depth_avg_allsed_avg_masked_10_60m_masked_trimmed_max)


# Surface SSC
# Mean 
ssc_surf_allsed_avg_masked_inner_masked_trimmed_avg = np.nanmean(ssc_surf_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_inner_masked_trimmed_avg)
ssc_surf_allsed_avg_masked_mid_masked_trimmed_avg = np.nanmean(ssc_surf_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_mid_masked_trimmed_avg)
ssc_surf_allsed_avg_masked_outer_masked_trimmed_avg = np.nanmean(ssc_surf_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_outer_masked_trimmed_avg)
ssc_surf_allsed_avg_masked_10m_masked_trimmed_avg = np.nanmean(ssc_surf_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10m_masked_trimmed_avg)
ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_avg = np.nanmean(ssc_surf_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_avg)

# Standard deviation 
ssc_surf_allsed_avg_masked_inner_masked_trimmed_std = np.nanstd(ssc_surf_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) std ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_inner_masked_trimmed_std)
ssc_surf_allsed_avg_masked_mid_masked_trimmed_std = np.nanstd(ssc_surf_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) std ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_mid_masked_trimmed_std)
ssc_surf_allsed_avg_masked_outer_masked_trimmed_std = np.nanstd(ssc_surf_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) std ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_outer_masked_trimmed_std)
ssc_surf_allsed_avg_masked_10m_masked_trimmed_std = np.nanstd(ssc_surf_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m std total ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10m_masked_trimmed_std)
ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_std = np.nanstd(ssc_surf_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m std ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_std)

# Min
ssc_surf_allsed_avg_masked_inner_masked_trimmed_min = np.nanmin(ssc_surf_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) min ssc_surf_allsed_avg(kg/m3): ', ssc_surf_allsed_avg_masked_inner_masked_trimmed_min)
ssc_surf_allsed_avg_masked_mid_masked_trimmed_min = np.nanmin(ssc_surf_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) min ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_mid_masked_trimmed_min)
ssc_surf_allsed_avg_masked_outer_masked_trimmed_min = np.nanmin(ssc_surf_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) min ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_outer_masked_trimmed_min)
ssc_surf_allsed_avg_masked_10m_masked_trimmed_min = np.nanmin(ssc_surf_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m min ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10m_masked_trimmed_min)
ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_min = np.nanmin(ssc_surf_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m min ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_min)

# Max
ssc_surf_allsed_avg_masked_inner_masked_trimmed_max = np.nanmax(ssc_surf_allsed_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) max ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_inner_masked_trimmed_max)
ssc_surf_allsed_avg_masked_mid_masked_trimmed_max = np.nanmax(ssc_surf_allsed_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) max ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_mid_masked_trimmed_max)
ssc_surf_allsed_avg_masked_outer_masked_trimmed_max = np.nanmax(ssc_surf_allsed_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) max ssc_surf_allsed_avg (kg/m3: ', ssc_surf_allsed_avg_masked_outer_masked_trimmed_max)
ssc_surf_allsed_avg_masked_10m_masked_trimmed_max = np.nanmax(ssc_surf_allsed_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m max tssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10m_masked_trimmed_max)
ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_max = np.nanmax(ssc_surf_allsed_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m max ssc_surf_allsed_avg (kg/m3): ', ssc_surf_allsed_avg_masked_10_60m_masked_trimmed_max)



# Ssflux 
# Mean in 0 - 10 m depth in east
# u
ssflux_depth_int_u_avg_masked_10m_east_avg = np.nanmean(ssflux_depth_int_u_avg_masked_10m_east, axis=(0,1))
print('0 - 10 m mean total ssflux_depth_int_u_avg_masked_10m_east (kg/ms):', ssflux_depth_int_u_avg_masked_10m_east_avg)
# v
ssflux_depth_int_v_avg_masked_10m_east_avg = np.nanmean(ssflux_depth_int_v_avg_masked_10m_east, axis=(0,1))
print('0 - 10 m mean total ssflux_depth_int_v_avg_masked_10m_east (kg/ms):', ssflux_depth_int_v_avg_masked_10m_east_avg)
# Mean in 0 - 10 m depth in west
# u
ssflux_depth_int_u_avg_masked_10m_west_avg = np.nanmean(ssflux_depth_int_u_avg_masked_10m_west, axis=(0,1))
print('0 - 10 m mean total ssflux_depth_int_u_avg_masked_10m_west (kg/ms):', ssflux_depth_int_u_avg_masked_10m_west_avg)
# v
ssflux_depth_int_v_avg_masked_10m_west_avg = np.nanmean(ssflux_depth_int_v_avg_masked_10m_west, axis=(0,1))
print('0 - 10 m mean total ssflux_depth_int_v_avg_masked_10m_west (kg/ms):', ssflux_depth_int_v_avg_masked_10m_west_avg)
# -- Mean --
# Depth-integrated
# Inner
ssflux_depth_int_u_avg_masked_inner_masked_trimmed_avg = np.nanmean(ssflux_depth_int_u_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_u_avg (kg/ms): ', ssflux_depth_int_u_avg_masked_inner_masked_trimmed_avg)
ssflux_depth_int_v_avg_masked_inner_masked_trimmed_avg = np.nanmean(ssflux_depth_int_v_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_v_avg (kg/ms): ', ssflux_depth_int_v_avg_masked_inner_masked_trimmed_avg)
# Mid
ssflux_depth_int_u_avg_masked_mid_masked_trimmed_avg = np.nanmean(ssflux_depth_int_u_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_u_avg (kg/ms): ', ssflux_depth_int_u_avg_masked_mid_masked_trimmed_avg)
ssflux_depth_int_v_avg_masked_mid_masked_trimmed_avg = np.nanmean(ssflux_depth_int_v_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_v_avg (kg/ms): ', ssflux_depth_int_v_avg_masked_mid_masked_trimmed_avg)
# Outer
ssflux_depth_int_u_avg_masked_outer_masked_trimmed_avg = np.nanmean(ssflux_depth_int_u_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_u_avg (kg/ms): ', ssflux_depth_int_u_avg_masked_outer_masked_trimmed_avg)
ssflux_depth_int_v_avg_masked_outer_masked_trimmed_avg = np.nanmean(ssflux_depth_int_v_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_v_avg (kg/ms): ', ssflux_depth_int_v_avg_masked_outer_masked_trimmed_avg)
# 0 - 10
ssflux_depth_int_u_avg_masked_10m_masked_trimmed_avg = np.nanmean(ssflux_depth_int_u_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_u_avg (kg/ms): ', ssflux_depth_int_u_avg_masked_10m_masked_trimmed_avg)
ssflux_depth_int_v_avg_masked_10m_masked_trimmed_avg = np.nanmean(ssflux_depth_int_v_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_v_avg (kg/ms): ', ssflux_depth_int_v_avg_masked_10m_masked_trimmed_avg)
# 10 - 60
ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_avg = np.nanmean(ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_u_avg (kg/ms): ', ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_avg)
ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_avg = np.nanmean(ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_v_avg (kg/ms): ', ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_avg)
# Surface 
# Inner
ssflux_u_surf_avg_masked_inner_masked_trimmed_avg = np.nanmean(ssflux_u_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_u_surf_avg (kg/ms): ', ssflux_u_surf_avg_masked_inner_masked_trimmed_avg)
ssflux_v_surf_avg_masked_inner_masked_trimmed_avg = np.nanmean(ssflux_v_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_v_surf_avg (kg/ms): ', ssflux_v_surf_avg_masked_inner_masked_trimmed_avg)


# -- Standard deviation--
# Depth-integrated
# Inner
ssflux_depth_int_u_avg_masked_inner_masked_trimmed_std = np.nanstd(ssflux_depth_int_u_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_u_std (kg/ms): ', ssflux_depth_int_u_avg_masked_inner_masked_trimmed_std)
ssflux_depth_int_v_avg_masked_inner_masked_trimmed_std = np.nanstd(ssflux_depth_int_v_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_v_std (kg/ms): ', ssflux_depth_int_v_avg_masked_inner_masked_trimmed_std)
# Mid
ssflux_depth_int_u_avg_masked_mid_masked_trimmed_std = np.nanstd(ssflux_depth_int_u_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_u_std (kg/ms): ', ssflux_depth_int_u_avg_masked_mid_masked_trimmed_std)
ssflux_depth_int_v_avg_masked_mid_masked_trimmed_std = np.nanstd(ssflux_depth_int_v_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_v_std (kg/ms): ', ssflux_depth_int_v_avg_masked_mid_masked_trimmed_std)
# Outer
ssflux_depth_int_u_avg_masked_outer_masked_trimmed_std = np.nanstd(ssflux_depth_int_u_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_u_std (kg/ms): ', ssflux_depth_int_u_avg_masked_outer_masked_trimmed_std)
ssflux_depth_int_v_avg_masked_outer_masked_trimmed_std = np.nanstd(ssflux_depth_int_v_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_v_std (kg/ms): ', ssflux_depth_int_v_avg_masked_outer_masked_trimmed_std)
# 0 - 10
ssflux_depth_int_u_avg_masked_10m_masked_trimmed_std = np.nanstd(ssflux_depth_int_u_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_u_std (kg/ms): ', ssflux_depth_int_u_avg_masked_10m_masked_trimmed_std)
ssflux_depth_int_v_avg_masked_10m_masked_trimmed_std = np.nanstd(ssflux_depth_int_v_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_v_std (kg/ms): ', ssflux_depth_int_v_avg_masked_10m_masked_trimmed_std)
# 10 - 60
ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_std = np.nanstd(ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_u_std (kg/ms): ', ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_std)
ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_std = np.nanstd(ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_v_std (kg/ms): ', ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_std)
# Surface 
# Inner
ssflux_u_surf_avg_masked_inner_masked_trimmed_std = np.nanstd(ssflux_u_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_u_surf_std (kg/ms): ', ssflux_u_surf_avg_masked_inner_masked_trimmed_std)
ssflux_v_surf_avg_masked_inner_masked_trimmed_std = np.nanstd(ssflux_v_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_v_surf_std (kg/ms): ', ssflux_v_surf_avg_masked_inner_masked_trimmed_std)


# -- Min --
# Depth-integrated
# Inner
ssflux_depth_int_u_avg_masked_inner_masked_trimmed_min = np.nanmin(ssflux_depth_int_u_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_u_min (kg/ms): ', ssflux_depth_int_u_avg_masked_inner_masked_trimmed_min)
ssflux_depth_int_v_avg_masked_inner_masked_trimmed_min = np.nanmin(ssflux_depth_int_v_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_v_min (kg/ms): ', ssflux_depth_int_v_avg_masked_inner_masked_trimmed_min)
# Mid
ssflux_depth_int_u_avg_masked_mid_masked_trimmed_min = np.nanmin(ssflux_depth_int_u_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_u_min (kg/ms): ', ssflux_depth_int_u_avg_masked_mid_masked_trimmed_min)
ssflux_depth_int_v_avg_masked_mid_masked_trimmed_min = np.nanmin(ssflux_depth_int_v_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_v_min (kg/ms): ', ssflux_depth_int_v_avg_masked_mid_masked_trimmed_min)
# Outer
ssflux_depth_int_u_avg_masked_outer_masked_trimmed_min = np.nanmin(ssflux_depth_int_u_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_u_min (kg/ms): ', ssflux_depth_int_u_avg_masked_outer_masked_trimmed_min)
ssflux_depth_int_v_avg_masked_outer_masked_trimmed_min = np.nanmin(ssflux_depth_int_v_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_v_min (kg/ms): ', ssflux_depth_int_v_avg_masked_outer_masked_trimmed_min)
# 0 - 10
ssflux_depth_int_u_avg_masked_10m_masked_trimmed_min = np.nanmin(ssflux_depth_int_u_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_u_min (kg/ms): ', ssflux_depth_int_u_avg_masked_10m_masked_trimmed_min)
ssflux_depth_int_v_avg_masked_10m_masked_trimmed_min = np.nanmin(ssflux_depth_int_v_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_v_min (kg/ms): ', ssflux_depth_int_v_avg_masked_10m_masked_trimmed_min)
# 10 - 60
ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_min = np.nanmin(ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_u_min (kg/ms): ', ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_min)
ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_min = np.nanmin(ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_v_min (kg/ms): ', ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_min)
# Surface 
# Inner
ssflux_u_surf_avg_masked_inner_masked_trimmed_min = np.nanmin(ssflux_u_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_u_surf_min (kg/ms): ', ssflux_u_surf_avg_masked_inner_masked_trimmed_min)
ssflux_v_surf_avg_masked_inner_masked_trimmed_min = np.nanmin(ssflux_v_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_v_surf_min (kg/ms): ', ssflux_v_surf_avg_masked_inner_masked_trimmed_min)

# -- Max -- 
# Depth-integrated
# Inner
ssflux_depth_int_u_avg_masked_inner_masked_trimmed_max = np.nanmax(ssflux_depth_int_u_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_u_max (kg/ms): ', ssflux_depth_int_u_avg_masked_inner_masked_trimmed_max)
ssflux_depth_int_v_avg_masked_inner_masked_trimmed_max = np.nanmax(ssflux_depth_int_v_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_depth_int_v_max (kg/ms): ', ssflux_depth_int_v_avg_masked_inner_masked_trimmed_max)
# Mid
ssflux_depth_int_u_avg_masked_mid_masked_trimmed_max = np.nanmax(ssflux_depth_int_u_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_u_max (kg/ms): ', ssflux_depth_int_u_avg_masked_mid_masked_trimmed_max)
ssflux_depth_int_v_avg_masked_mid_masked_trimmed_max = np.nanmax(ssflux_depth_int_v_avg_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged ssflux_depth_int_v_max (kg/ms): ', ssflux_depth_int_v_avg_masked_mid_masked_trimmed_max)
# Outer
ssflux_depth_int_u_avg_masked_outer_masked_trimmed_max = np.nanmax(ssflux_depth_int_u_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_u_max (kg/ms): ', ssflux_depth_int_u_avg_masked_outer_masked_trimmed_max)
ssflux_depth_int_v_avg_masked_outer_masked_trimmed_max = np.nanmax(ssflux_depth_int_v_avg_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged ssflux_depth_int_v_max (kg/ms): ', ssflux_depth_int_v_avg_masked_outer_masked_trimmed_max)
# 0 - 10
ssflux_depth_int_u_avg_masked_10m_masked_trimmed_max = np.nanmax(ssflux_depth_int_u_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_u_max (kg/ms): ', ssflux_depth_int_u_avg_masked_10m_masked_trimmed_max)
ssflux_depth_int_v_avg_masked_10m_masked_trimmed_max = np.nanmax(ssflux_depth_int_v_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean ssflux_depth_int_v_max (kg/ms): ', ssflux_depth_int_v_avg_masked_10m_masked_trimmed_max)
# 10 - 60
ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_max = np.nanmax(ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_u_max (kg/ms): ', ssflux_depth_int_u_avg_masked_10_60m_masked_trimmed_max)
ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_max = np.nanmax(ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean ssflux_depth_int_v_max (kg/ms): ', ssflux_depth_int_v_avg_masked_10_60m_masked_trimmed_max)
# Surface 
# Inner
ssflux_u_surf_avg_masked_inner_masked_trimmed_max = np.nanmax(ssflux_u_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_u_surf_max (kg/ms): ', ssflux_u_surf_avg_masked_inner_masked_trimmed_max)
ssflux_v_surf_avg_masked_inner_masked_trimmed_max = np.nanmax(ssflux_v_surf_avg_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged ssflux_v_surf_max (kg/ms): ', ssflux_v_surf_avg_masked_inner_masked_trimmed_max)



# Erosoion & Deposition
roms_net_erodepo_wland_masked_trimmed_meantot = np.mean(roms_net_erodepo_bedthick_wland_masked_trimmed)
print('mean net erosion/deposition over full domain (m): ', roms_net_erodepo_wland_masked_trimmed_meantot)
# Make functions to find positive and negative
# Mean erosion (negative valules)
def mean_negative(L):
    # Get all negative numbers into another list
    neg_only = [x for x in L if x < 0]
    if neg_only:
        return sum(neg_only) /  len(neg_only)
    raise ValueError('No negative numbers in input')
    
# Call the function
roms_net_erodepo_bedthick_wland_masked_trimmed_meanneg = mean_negative(roms_net_erodepo_bedthick_wland_masked_trimmed.values.ravel())
print('roms_net_erodepo_bedthick_wland_masked_trimmed_mean neg/erosion (m): ', roms_net_erodepo_bedthick_wland_masked_trimmed_meanneg)

# Mean deposition (positive values)
def mean_positive(L):
    # Get all positive numbers into another list
    pos_only = [x for x in L if x > 0]
    if pos_only:
        return sum(pos_only) /  len(pos_only)
    raise ValueError('No postive numbers in input')
    
# Call the function
roms_net_erodepo_bedthick_wland_masked_trimmed_meanpos = mean_positive(roms_net_erodepo_bedthick_wland_masked_trimmed.values.ravel())
print('roms_net_erodepo_bedthick_wland_masked_trimmed_mean pos/deposition (m): ', roms_net_erodepo_bedthick_wland_masked_trimmed_meanpos)

# Standard deviation 
# Standard deviation erosion (negative valules)
def std_negative(L):
    # Get all negative numbers into another list
    neg_only = [x for x in L if x < 0]
    if neg_only:
        return np.std(neg_only)
    raise ValueError('No negative numbers in input')
    
# Call the function
roms_net_erodepo_bedthick_wland_masked_trimmed_stdneg = std_negative(roms_net_erodepo_bedthick_wland_masked_trimmed.values.ravel())
print('roms_net_erodepo_bedthick_wland_masked_trimmed std neg/erosion (m): ', roms_net_erodepo_bedthick_wland_masked_trimmed_stdneg)

# Standard deviation deposition (positive values)
def std_positive(L):
    # Get all positive numbers into another list
    pos_only = [x for x in L if x > 0]
    if pos_only:
        return np.std(pos_only)
    raise ValueError('No postive numbers in input')
    
# Call the function
roms_net_erodepo_bedthick_wland_masked_trimmed_stdpos = std_positive(roms_net_erodepo_bedthick_wland_masked_trimmed.values.ravel())
print('roms_net_erodepo_bedthick_wland_masked_trimmed std pos/deposition (m): ', roms_net_erodepo_bedthick_wland_masked_trimmed_stdpos)


# Calculate percent of space that time-averaged, depth-integrated ss fluxes were
# westward 
# Percent negative flux
def pcnt_negative(L):
    # Get all negative numbers into another list
    neg_only = [x for x in L if x < 0]
    print('len neg only: ', len(neg_only))
    print('len L: ', len(L))
    # Get length of L without nans
    nans = [y for y in L if np.isnan(y)]
    print('len nan: ', len(nans))
    if neg_only:
        return len(neg_only)/(len(L)-len(nans))
    raise ValueError('No negative numbers in input')

# Call the function
ssflux_depth_int_u_avg_wland_masked_trimmed_slice_pcnt_west = pcnt_negative(ssflux_depth_int_u_avg_wland_masked_trimmed_slice.values.ravel())
print('Percent of space that depth-integrated ssc flux is westward: ', ssflux_depth_int_u_avg_wland_masked_trimmed_slice_pcnt_west)

# Percent positive flux
def pcnt_positive(L):
    # Get all negative numbers into another list
    pos_only = [x for x in L if x > 0]
    print('len pos only: ', len(pos_only))
    print('len L: ', len(L))
    # Get length of L without nans
    nans = [y for y in L if np.isnan(y)]
    print('len nan: ', len(nans))
    if pos_only:
        return len(pos_only)/(len(L)-len(nans))
    raise ValueError('No positive numbers in input')

# Call the function
ssflux_depth_int_u_avg_wland_masked_trimmed_slice_pcnt_east = pcnt_positive(ssflux_depth_int_u_avg_wland_masked_trimmed_slice.values.ravel())
print('Percent of space that depth-integrated ssc flux is eastward: ', ssflux_depth_int_u_avg_wland_masked_trimmed_slice_pcnt_east)




# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Set up the data
roms_time_avg_ssc_depth_avg_surf_bot_erodepo = xr.Dataset(
    data_vars=dict(
        ssc_depth_avg_time_avg=(['y','x'], ssc_depth_avg_allsed_avg_trimmed.values),
        ssflux_u_depth_int_time_avg=(['y_slice','x_slice'], ssflux_depth_int_u_avg_wland_masked_trimmed_slice.values),
        ssflux_v_depth_int_time_avg=(['y_slice','x_slice'], ssflux_depth_int_v_avg_wland_masked_trimmed_slice.values),
        ssc_surf_time_avg=(['y','x'], ssc_surf_allsed_avg_wland_masked_trimmed.values),
        ssflux_u_surf_time_avg=(['y_slice','x_slice'], ssflux_u_surf_avg_wland_masked_trimmed_slice.values),
        ssflux_v_surf_time_avg=(['y_slice','x_slice'], ssflux_v_surf_avg_wland_masked_trimmed_slice.values),
        ssc_1m_above_seafloor_time_avg=(['y','x'], ssc_1m_avg_wland_masked_trimmed.values),
        ssflux_u_1m_above_seafloor_time_avg=(['y_slice','x_slice'], ssflux_u_1m_avg_wland_masked_trimmed_slice.values),
        ssflux_v_1m_above_seafloor_time_avg=(['y_slice','x_slice'], ssflux_v_1m_avg_wland_masked_trimmed_slice.values),
        net_ero_depo=(['y','x'], roms_net_erodepo_bedthick_wland_masked_trimmed.values*100),
        ),
    coords=dict(
        x_full=('x', x_rho_flat_trimmed),
        x_slice=('x_slice', x_rho_flat_trimmed_slice),
        y_full=('y', y_rho_flat), 
        y_slice=('y_slice', y_rho_flat_slice)
        ),
    attrs=dict(description='Time-averaged ROMS output including depth-averaged, surface, and 1 meter above seafloor suspended sediment concentrations (mg/m3) and fluxes, and net deposition (cm)')) 
# Add more metadata?
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssc_depth_avg_time_avg.name='depth-averaged, time-averaged suspended sediment concentration (kg/m3)'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssflux_u_depth_int_time_avg.name='depth-integrated, time-averaged suspended sediment flux in u direction (kg/m/s)'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssflux_v_depth_int_time_avg.name='depth-integrated, time-averaged suspended sediment flux in v direction (kg/m/s)'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssc_surf_time_avg.name='time-averaged surface suspended sediment concentration (kg/m3)'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssflux_u_surf_time_avg.name='time-averaged surface suspended sediment flux in u direction (kg/m2/s)'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssflux_v_surf_time_avg.name='time-averaged surface suspended sediment flux in v direction (kg/m2/s)'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssc_1m_above_seafloor_time_avg.name='time-averaged suspended sediment concentration (kg/m3) 1m above seafloor'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssflux_u_1m_above_seafloor_time_avg.name='time-averaged suspended sediment flux in u direction (kg/m2/s) 1m above seafloor'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.ssflux_v_1m_above_seafloor_time_avg.name='time-averaged suspended sediment flux in v direction (kg/m2/s) 1m above seafloor'
roms_time_avg_ssc_depth_avg_surf_bot_erodepo.net_ero_depo.name='net deposition (cm)'

# Save to a netcdf
#roms_time_avg_ssc_depth_avg_surf_bot_erodepo.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Paper1_Take2/fig8_roms_depth_avg_surf_bot_ssc_ssflux_erodepo.nc')


