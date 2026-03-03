#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 15:08:59 2023

@author: brun1463
"""

################## Scatter Plot: Bed Stress & Sediment Flux ###################
# The purpose of this script is to make a scatter plot of the bed shear stress
# on the x-axis and depth-integrated SSC on the y-axis. This will be done for 
# a spatially-averaged (masked) domain, just the inner shelf, just the middle
# shelf, and just the outer shelf (which will all be defined based on water 
# depth).
#
# Notes:
# - The script will first be set up to do the depth-integrated SSC flux 
#   for all sediment classes added together 
# - Also look into doing Hsig on the x axis, u currents on y axis, and dot 
#   color being the magnitude of the sediment flux
#   - Do this for the spatially-averaged (masked) shelf, then also for east 
#     (masked) and west (masked)
# - Do east/west NOT u/v since this is a scatter plot
# - This script has been updated to use the 2020 ROMS output
###############################################################################


# Load in the packages 
import numpy as np
import xarray as xr
import pandas as pd
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec
#from matplotlib import cm, ticker
#from matplotlib import transforms 
import cmocean
#import matplotlib.ticker as tick
#import matplotlib.patches as patches
#from matplotlib.colors import LinearSegmentedColormap
#import ESMF
#import xesmf as xe
from glob import glob


# Set a universal fontsize
fontsize = 20

# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)

# Prevent tick labels from overlapping
matplotlib.rcParams['xtick.major.pad'] = 12
matplotlib.rcParams['ytick.major.pad'] = 12


# Load in the model grid
#grid = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Include/KakAKgrd_shelf_big010_smooth006.nc') # UPDATE PATH
grid = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Grids/KakAKgrd_shelf_big010_smooth006.nc') # CHANGE PATH

# Pull out some dimensions
eta_rho_len = len(grid.eta_rho)
xi_rho_len = len(grid.xi_rho)

# Load in ROMS wave input
#roms_wave = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Include/wave_forcing_file_kaktovik_shelf_era5_data005.nc')
wave_frc = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/Forcing_files/wave_forcing_file_kaktovik_shelf_ww3_2020_data002.nc')
# Convert time
wave_frc['wave_time'] = pd.to_datetime(wave_frc.wave_time.values+86400, origin=pd.datetime(1999,12,31), unit='s')

# Load in the rho masks 
mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')



# So we want a time series of the bed shear stress and the depth-integrated SSC flux
# for all space and time in the domain
# Also want time series of significant wave height and eastward currents 
# If we get it for the whole domain first, then we can sort it into the different
# sections of the domain 

# Make a bunch of functions that will open a given model output and pull out
# the time series that we want 

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
    #print(time_lengths)

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


# Make a function to pull out the spatial time series of bed stress from both
# wave and currents 
def get_bstrcwmax_ubar_vbar_time_series(filename):
    """
    This function opens a given model output and pulls out the time series
    for the wave-current bed shear stress and the eastward depth-averaged
    currents.

    Parameters
    ----------
    filename : Path to model output file (string)

    Returns
    -------
    bstrcwmax_tmp: Wave-Current bed shear stress time series on rho points
    ubar_eastward_tmp: Depth-averaged currents in east-west direction on rho points

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Pull out the bstrcwmax time series
    bstrcwmax_tmp = model_output.bstrcwmax[:,:,:]
    #print('bastrcwmax_tmp shape: ', bstrcwmax_tmp.shape)
    #print('model_output.bstrcwmax shape: ', model_output.bstrcwmax.shape)
    
    # Pull out the depth-averaged eastward currents 
    ubar_eastward_tmp = model_output.ubar_eastward[:,:,:]
    #print('ubar_eastward_tmp shape: ', ubar_eastward_tmp.shape)
    #print('model_output.ubar_eastward shape: ', model_output.ubar_eastward.shape)
    
    # Pull out the depth-averaged northward currents 
    vbar_northward_tmp = model_output.vbar_northward[:,:,:]
    
    # Return these arrays
    return(bstrcwmax_tmp, ubar_eastward_tmp, vbar_northward_tmp)
    
    
# Make a function to get time series of bottom u_eastward and v_northward currents
def get_bottom_ueast_vnorth_time_series(filename):
    """
    

    Parameters
    ----------
    filename : Path to model output file (string)

    Returns
    -------
    None.
    u_eastward_bottom_tmp: Time series of eastward bottom currents 
    v_northward_bottom_tmp: Time series of northward bottom currents

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Pull out the bottom u_eastward currents
    u_eastward_bottom_tmp = model_output.u_eastward[:,0,:,:]
    
    # Pull out the bottom v_northward currents
    v_northward_bottom_tmp = model_output.v_northward[:,0,:,:]
    
    # Return these arrays
    return(u_eastward_bottom_tmp, v_northward_bottom_tmp)


# Make a function to calculate the sediment flux in the 
# eastward direction for all sediment classes combined
def calc_eastward_ssc_flux_allsed(filename):
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

    """

    # Load in the model output
    model_output = xr.open_dataset(filename)

    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03

    # To collapse to horizontal, multiply each layer by its
    # thickness
    # Calculate the time-varying thickness of the cells
    dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)

    # Pull out the eastward velocities at all times, depths, spaces
    u_eastward_tmp = model_output.u_eastward

    # Pull out the thickness of the cell in the y direction 
    dy = 1.0/model_output.pn

    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_allsed = ((ssc_allsed_tmp*u_eastward_tmp)*(dz))
    
    # Then depth-integrated by summing over depth and dividing by dy
    depth_int_ssc_flux_east_allsed = (ssc_flux_allsed.sum(dim='s_rho'))

    # Return the depth-integrated eastward flux for all sediment classes
    return(depth_int_ssc_flux_east_allsed)
    

# Make a function to calculate the sediment flux in the 
# northward direction for all sediment classes combined
def calc_northward_ssc_flux_allsed(filename):
    """
    The purpose of this function is to take a given model output file, load 
    in the output, and caluclate the time-averaged, depth-integrated suspended 
    sediment flux in the northward direction for all sediment classes added
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

    # Pull out the eastward velocities at all times, depths, spaces
    v_northward_tmp = model_output.v_northward

    # Pull out the thickness of the cell in the y direction 
    dx = 1.0/model_output.pm

    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_allsed = ((ssc_allsed_tmp*v_northward_tmp)*(dz))

    # Then depth-integrated by summing over depth and dividing by dx
    depth_int_ssc_flux_north_allsed = (ssc_flux_allsed.sum(dim='s_rho'))

    # Return the depth-integrated eastward flux for all sediment classes
    return(depth_int_ssc_flux_north_allsed)



# Loop through model output to get the different
# time series and merge them each into an array
# First, get all the output file names
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

# Before looping through to make the movie, get the length of the model time
full_time_len, time_steps, time_lengths = get_model_time(file_names2, num_files)

# Make empty arrays to hold the full time series 
bstrcwmax_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ubar_eastward_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))
vbar_northward_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_east_allsed_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_north_allsed_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))
u_eastward_bottom_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_northward_bottom_fullrun = np.empty((full_time_len, eta_rho_len, xi_rho_len))

# Loop through the model output files to call the functions that pull out the 
# time series and save them into the arrays
# Initiate time
time_step = 0

# Loop through the model output files
for j in range(num_files):
    print('j: ', j, flush=True)
    
    # Bed stress and eastward currents 
    # Call the function
    bstrcwmax_tmp, ubar_eastward_tmp, vbar_northward_tmp = get_bstrcwmax_ubar_vbar_time_series(file_names2[j])
    
    # Bottom currents 
    u_eastward_bottom_tmp, v_northward_bottom_tmp = get_bottom_ueast_vnorth_time_series(file_names2[j])
    
    # Depth-integrated sediment flux 
    # Call the function to get ssc flux
    # Note: This returns the flux for all time steps in this output file 
    depth_int_ssc_flux_east_allsed_tmp2 = calc_eastward_ssc_flux_allsed(file_names2[j])
    depth_int_ssc_flux_north_allsed_tmp2 = calc_northward_ssc_flux_allsed(file_names2[j])
    
    # Save the output into arrays
    start = int(time_step)
    end = int(time_step+time_lengths[j])
    bstrcwmax_fullrun[start:end,:,:] = bstrcwmax_tmp
    ubar_eastward_fullrun[start:end,:,:] = ubar_eastward_tmp
    vbar_northward_fullrun[start:end,:,:] = vbar_northward_tmp
    depth_int_ssc_flux_east_allsed_fullrun[start:end,:,:] = depth_int_ssc_flux_east_allsed_tmp2
    depth_int_ssc_flux_north_allsed_fullrun[start:end,:,:] = depth_int_ssc_flux_north_allsed_tmp2
    u_eastward_bottom_fullrun[start:end,:,:] = u_eastward_bottom_tmp
    v_northward_bottom_fullrun[start:end,:,:] = v_northward_bottom_tmp
    
    # Update the base time_step
    time_step = time_step + time_lengths[j]

    # Check time
    print('j in dep net and ssc flux loop: ', j, flush=True)
    print('start: ', start, flush=True)
    print('end: ', end, flush=True)
    print('time_step: ', time_step)


# Pull out the time series of significant wave height
Hsig_fullrun = wave_frc.Hwave
# Pull out the time series of bottom wave orbital velocity 
Uwave_rms_fullrun = wave_frc.Uwave_rms

# Now that we have all of the time series, multiply them all by the masks 
# that cut out the nudged areas on east/west and depths greater than 60 m
# Let's try the nan mask since that is good for math
# Make empty arrays to hold new values
bstrcwmax_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ubar_eastward_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))
vbar_northward_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_east_allsed_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_north_allsed_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))
Hsig_fullrun_masked = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
Uwave_rms_fullrun_masked = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
u_eastward_bottom_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_northward_bottom_fullrun_masked = np.empty((full_time_len, eta_rho_len, xi_rho_len))

# Loop through time to multiply by mask
# Time series from model output
for k in range(full_time_len):
    bstrcwmax_fullrun_masked[k,:,:] = bstrcwmax_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    ubar_eastward_fullrun_masked[k,:,:] = ubar_eastward_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    vbar_northward_fullrun_masked[k,:,:] = vbar_northward_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    depth_int_ssc_flux_east_allsed_fullrun_masked[k,:,:] = depth_int_ssc_flux_east_allsed_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    depth_int_ssc_flux_north_allsed_fullrun_masked[k,:,:] = depth_int_ssc_flux_north_allsed_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    u_eastward_bottom_fullrun_masked[k,:,:] = u_eastward_bottom_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    v_northward_bottom_fullrun_masked[k,:,:] = v_northward_bottom_fullrun[k,:,:] * mask_rho_nan.nudge_mask_rho_nan
    
# Time series for wave forcing 
for j in range(len(wave_frc.wave_time)):
    Hsig_fullrun_masked[j,:,:] = Hsig_fullrun[j,:,:] * mask_rho_nan.nudge_mask_rho_nan
    Uwave_rms_fullrun_masked[j,:,:] = Uwave_rms_fullrun[j,:,:] * mask_rho_nan.nudge_mask_rho_nan

# Calculate the magnitude of the depth-integrated sediment flux
depth_int_sscflux_allsed_mag_fullrun_masked = np.sqrt(((depth_int_ssc_flux_east_allsed_fullrun_masked)**2)+((depth_int_ssc_flux_north_allsed_fullrun_masked)**2))

# Calculate the magnitude of the depth-averaged currents 
depth_avg_cur_mag_fullrun_masked = np.sqrt(((ubar_eastward_fullrun_masked)**2)+((vbar_northward_fullrun_masked)**2))

# Calculate the magnitude of the bottom currents 
bottom_cur_mag_fullrun_masked = np.sqrt(((u_eastward_bottom_fullrun_masked)**2)+((v_northward_bottom_fullrun_masked)**2))

# Spatially-average the results 
bstrcwmax_fullrun_masked_avg = np.nanmean(bstrcwmax_fullrun_masked, axis=(1,2))
ubar_eastward_fullrun_masked_avg = np.nanmean(ubar_eastward_fullrun_masked, axis=(1,2))
vbar_northward_fullrun_masked_avg = np.nanmean(vbar_northward_fullrun_masked, axis=(1,2))
depth_int_sscflux_allsed_mag_fullrun_masked_avg = np.nanmean(depth_int_sscflux_allsed_mag_fullrun_masked, axis=(1,2))
Hsig_fullrun_masked_avg = np.nanmean(Hsig_fullrun_masked, axis=(1,2))
Uwave_rms_fullrun_masked_avg = np.nanmean(Uwave_rms_fullrun_masked, axis=(1,2))
depth_avg_cur_mag_fullrun_masked_avg = np.nanmean(depth_avg_cur_mag_fullrun_masked, axis=(1,2))
u_eastward_bottom_fullrun_masked_avg = np.nanmean(u_eastward_bottom_fullrun_masked, axis=(1,2))
v_northward_bottom_fullrun_masked_avg = np.nanmean(v_northward_bottom_fullrun_masked, axis=(1,2))
bottom_cur_mag_fullrun_masked_avg = np.nanmean(bottom_cur_mag_fullrun_masked, axis=(1,2))

# Also split results into different sections of the shelf, then take average o
# over just those sections 
# Slice to those regions
# East
Hsig_fullrun_masked_east = Hsig_fullrun_masked[:,:,304:]
Uwave_rms_fullrun_masked_east = Uwave_rms_fullrun_masked[:,:,304:]
ubar_eastward_fullrun_masked_east = ubar_eastward_fullrun_masked[:,:,304:]
vbar_northward_fullrun_masked_east = vbar_northward_fullrun_masked[:,:,304:]
depth_int_sscflux_allsed_mag_fullrun_masked_east = depth_int_sscflux_allsed_mag_fullrun_masked[:,:,304:]
bstrcwmax_fullrun_masked_east = bstrcwmax_fullrun_masked[:,:,304:]
u_eastward_bottom_fullrun_masked_east = u_eastward_bottom_fullrun_masked[:,:,304:]
v_northward_bottom_fullrun_masked_east = v_northward_bottom_fullrun_masked[:,:,304:]
# West
Hsig_fullrun_masked_west = Hsig_fullrun_masked[:,:,:304]
Uwave_rms_fullrun_masked_west = Uwave_rms_fullrun_masked[:,:,:304]
ubar_eastward_fullrun_masked_west = ubar_eastward_fullrun_masked[:,:,:304]
vbar_northward_fullrun_masked_west = vbar_northward_fullrun_masked[:,:,:304]
depth_int_sscflux_allsed_mag_fullrun_masked_west = depth_int_sscflux_allsed_mag_fullrun_masked[:,:,:304]
bstrcwmax_fullrun_masked_west = bstrcwmax_fullrun_masked[:,:,:304]
u_eastward_bottom_fullrun_masked_west = u_eastward_bottom_fullrun_masked[:,:,:304]
v_northward_bottom_fullrun_masked_west = v_northward_bottom_fullrun_masked[:,:,:304]
# Inner? (these are more complicated so maybe skip for now and do this for cumsum
# plots instead in a different script)
# Mid?
# Outer?
# Take the avgerage in these regions 
# East
Hsig_fullrun_masked_east_avg = np.nanmean(Hsig_fullrun_masked_east, axis=(1,2))
Uwave_rms_fullrun_masked_east_avg = np.nanmean(Uwave_rms_fullrun_masked_east, axis=(1,2))
ubar_eastward_fullrun_masked_east_avg = np.nanmean(ubar_eastward_fullrun_masked_east, axis=(1,2))
vbar_northward_fullrun_masked_east_avg = np.nanmean(vbar_northward_fullrun_masked_east, axis=(1,2))
depth_int_sscflux_allsed_mag_fullrun_masked_east_avg = np.nanmean(depth_int_sscflux_allsed_mag_fullrun_masked_east, axis=(1,2))
bstrcwmax_fullrun_masked_east_avg = np.nanmean(bstrcwmax_fullrun_masked_east, axis=(1,2))
u_eastward_bottom_fullrun_masked_east_avg = np.nanmean(u_eastward_bottom_fullrun_masked_east, axis=(1,2))
v_northward_bottom_fullrun_masked_east_avg = np.nanmean(v_northward_bottom_fullrun_masked_east, axis=(1,2))
# West 
Hsig_fullrun_masked_west_avg = np.nanmean(Hsig_fullrun_masked_west, axis=(1,2))
Uwave_rms_fullrun_masked_west_avg = np.nanmean(Uwave_rms_fullrun_masked_west, axis=(1,2))
ubar_eastward_fullrun_masked_west_avg = np.nanmean(ubar_eastward_fullrun_masked_west, axis=(1,2))
vbar_northward_fullrun_masked_west_avg = np.nanmean(vbar_northward_fullrun_masked_west, axis=(1,2))
depth_int_sscflux_allsed_mag_fullrun_masked_west_avg = np.nanmean(depth_int_sscflux_allsed_mag_fullrun_masked_west, axis=(1,2))
bstrcwmax_fullrun_masked_west_avg = np.nanmean(bstrcwmax_fullrun_masked_west, axis=(1,2))
u_eastward_bottom_fullrun_masked_west_avg = np.nanmean(u_eastward_bottom_fullrun_masked_west, axis=(1,2))
v_northward_bottom_fullrun_masked_west_avg = np.nanmean(v_northward_bottom_fullrun_masked_west, axis=(1,2))


# --- Split into inner, mid, and outer shelf by depth ---
# Make masks for each of these regions 
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

# Now multiply by the mask to get the different regions 
# Make empty arrays to hold new values
# Inner
Hsig_fullrun_masked_inner = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
Uwave_rms_fullrun_masked_inner = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
ubar_eastward_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
vbar_northward_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_east_allsed_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_north_allsed_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
bstrcwmax_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
u_eastward_bottom_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_northward_bottom_fullrun_masked_inner = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# Mid 
Hsig_fullrun_masked_mid = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
Uwave_rms_fullrun_masked_mid = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
ubar_eastward_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
vbar_northward_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_east_allsed_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_north_allsed_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
bstrcwmax_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
u_eastward_bottom_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_northward_bottom_fullrun_masked_mid = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# Outer
Hsig_fullrun_masked_outer = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
Uwave_rms_fullrun_masked_outer = np.empty((len(wave_frc.wave_time), eta_rho_len, xi_rho_len))
ubar_eastward_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
vbar_northward_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_east_allsed_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
depth_int_ssc_flux_north_allsed_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
bstrcwmax_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
u_eastward_bottom_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_northward_bottom_fullrun_masked_outer = np.empty((full_time_len, eta_rho_len, xi_rho_len))
# Loop through time to multiply by mask
# Time series from model output
for k in range(full_time_len):
    # Inner
    #Hsig_fullrun_masked_inner[k,:,:] = Hsig_fullrun[k,:,:] * inner_shelf_mask_rho_nan
    #Uwave_rms_fullrun_masked_inner[k,:,:] = Uwave_rms_fullrun[k,:,:] * inner_shelf_mask_rho_nan
    ubar_eastward_fullrun_masked_inner[k,:,:] = ubar_eastward_fullrun[k,:,:] * inner_shelf_mask_rho_nan
    vbar_northward_fullrun_masked_inner[k,:,:] = vbar_northward_fullrun_masked_inner[k,:,:] * inner_shelf_mask_rho_nan
    depth_int_ssc_flux_east_allsed_fullrun_masked_inner[k,:,:] = depth_int_ssc_flux_east_allsed_fullrun[k,:,:] * inner_shelf_mask_rho_nan
    depth_int_ssc_flux_north_allsed_fullrun_masked_inner[k,:,:] = depth_int_ssc_flux_north_allsed_fullrun[k,:,:] * inner_shelf_mask_rho_nan
    bstrcwmax_fullrun_masked_inner[k,:,:] = bstrcwmax_fullrun_masked_inner[k,:,:] * inner_shelf_mask_rho_nan
    u_eastward_bottom_fullrun_masked_inner[k,:,:] = u_eastward_bottom_fullrun_masked_inner[k,:,:] * inner_shelf_mask_rho_nan
    v_northward_bottom_fullrun_masked_inner[k,:,:] = v_northward_bottom_fullrun_masked_inner[k,:,:]  * inner_shelf_mask_rho_nan
    # Mid
    #Hsig_fullrun_masked_mid[k,:,:] = Hsig_fullrun[k,:,:] * mid_shelf_mask_rho_nan
    #Uwave_rms_fullrun_masked_mid[k,:,:] = Uwave_rms_fullrun[k,:,:] * mid_shelf_mask_rho_nan
    ubar_eastward_fullrun_masked_mid[k,:,:] = ubar_eastward_fullrun[k,:,:] * mid_shelf_mask_rho_nan
    vbar_northward_fullrun_masked_mid[k,:,:] = vbar_northward_fullrun_masked_mid[k,:,:] * mid_shelf_mask_rho_nan
    depth_int_ssc_flux_east_allsed_fullrun_masked_mid[k,:,:] = depth_int_ssc_flux_east_allsed_fullrun[k,:,:] * mid_shelf_mask_rho_nan
    depth_int_ssc_flux_north_allsed_fullrun_masked_mid[k,:,:] = depth_int_ssc_flux_north_allsed_fullrun[k,:,:] * mid_shelf_mask_rho_nan
    bstrcwmax_fullrun_masked_mid[k,:,:] = bstrcwmax_fullrun_masked_mid[k,:,:] * mid_shelf_mask_rho_nan
    u_eastward_bottom_fullrun_masked_mid[k,:,:] = u_eastward_bottom_fullrun_masked_mid[k,:,:] * mid_shelf_mask_rho_nan
    v_northward_bottom_fullrun_masked_mid[k,:,:] = v_northward_bottom_fullrun_masked_mid[k,:,:]  * mid_shelf_mask_rho_nan
    # Inner
    #Hsig_fullrun_masked_outer[k,:,:] = Hsig_fullrun[k,:,:] * outer_shelf_mask_rho_nan
    #Uwave_rms_fullrun_masked_outer[k,:,:] = Uwave_rms_fullrun[k,:,:] * outer_shelf_mask_rho_nan
    ubar_eastward_fullrun_masked_outer[k,:,:] = ubar_eastward_fullrun[k,:,:] * outer_shelf_mask_rho_nan
    vbar_northward_fullrun_masked_outer[k,:,:] = vbar_northward_fullrun_masked_outer[k,:,:] * outer_shelf_mask_rho_nan
    depth_int_ssc_flux_east_allsed_fullrun_masked_outer[k,:,:] = depth_int_ssc_flux_east_allsed_fullrun[k,:,:] * outer_shelf_mask_rho_nan
    depth_int_ssc_flux_north_allsed_fullrun_masked_outer[k,:,:] = depth_int_ssc_flux_north_allsed_fullrun[k,:,:] * outer_shelf_mask_rho_nan
    bstrcwmax_fullrun_masked_outer[k,:,:] = bstrcwmax_fullrun_masked_outer[k,:,:] * outer_shelf_mask_rho_nan
    u_eastward_bottom_fullrun_masked_outer[k,:,:] = u_eastward_bottom_fullrun_masked_outer[k,:,:] * outer_shelf_mask_rho_nan
    v_northward_bottom_fullrun_masked_outer[k,:,:] = v_northward_bottom_fullrun_masked_outer[k,:,:]  * outer_shelf_mask_rho_nan

# Separate loop for waves 
for kk in range(len(wave_frc.wave_time)):
    # Inner
    Hsig_fullrun_masked_inner[kk,:,:] = Hsig_fullrun[kk,:,:] * inner_shelf_mask_rho_nan
    Uwave_rms_fullrun_masked_inner[kk,:,:] = Uwave_rms_fullrun[kk,:,:] * inner_shelf_mask_rho_nan
    # Mid
    Hsig_fullrun_masked_mid[kk,:,:] = Hsig_fullrun[kk,:,:] * mid_shelf_mask_rho_nan
    Uwave_rms_fullrun_masked_mid[kk,:,:] = Uwave_rms_fullrun[kk,:,:] * mid_shelf_mask_rho_nan
    # Outer
    Hsig_fullrun_masked_outer[kk,:,:] = Hsig_fullrun[kk,:,:] * outer_shelf_mask_rho_nan
    Uwave_rms_fullrun_masked_outer[kk,:,:] = Uwave_rms_fullrun[kk,:,:] * outer_shelf_mask_rho_nan

# Calculate the magnitude of the depth-integrated sediment flux
depth_int_sscflux_allsed_mag_fullrun_masked_inner = np.sqrt(((depth_int_ssc_flux_east_allsed_fullrun_masked_inner)**2)+((depth_int_ssc_flux_north_allsed_fullrun_masked_inner)**2))
depth_int_sscflux_allsed_mag_fullrun_masked_mid = np.sqrt(((depth_int_ssc_flux_east_allsed_fullrun_masked_mid)**2)+((depth_int_ssc_flux_north_allsed_fullrun_masked_mid)**2))
depth_int_sscflux_allsed_mag_fullrun_masked_outer = np.sqrt(((depth_int_ssc_flux_east_allsed_fullrun_masked_outer)**2)+((depth_int_ssc_flux_north_allsed_fullrun_masked_outer)**2))

# Calculate the magnitude of the bottom currents
bottom_cur_mag_fullrun_masked_inner = np.sqrt(((u_eastward_bottom_fullrun_masked_inner)**2)+((v_northward_bottom_fullrun_masked_inner)**2))
bottom_cur_mag_fullrun_masked_mid = np.sqrt(((u_eastward_bottom_fullrun_masked_mid)**2)+((v_northward_bottom_fullrun_masked_mid)**2))
bottom_cur_mag_fullrun_masked_outer = np.sqrt(((u_eastward_bottom_fullrun_masked_outer)**2)+((v_northward_bottom_fullrun_masked_outer)**2))

# Take the spatial average over each region
# Inner
Hsig_fullrun_masked_inner_avg = np.nanmean(Hsig_fullrun_masked_inner, axis=(1,2))
Uwave_rms_fullrun_masked_inner_avg = np.nanmean(Uwave_rms_fullrun_masked_inner, axis=(1,2))
ubar_eastward_fullrun_masked_inner_avg = np.nanmean(ubar_eastward_fullrun_masked_inner, axis=(1,2))
vbar_northward_fullrun_masked_inner_avg = np.nanmean(vbar_northward_fullrun_masked_inner, axis=(1,2))
depth_int_sscflux_allsed_mag_fullrun_masked_inner_avg = np.nanmean(depth_int_sscflux_allsed_mag_fullrun_masked_inner, axis=(1,2))
bstrcwmax_fullrun_masked_inner_avg = np.nanmean(bstrcwmax_fullrun_masked_inner, axis=(1,2))
u_eastward_bottom_fullrun_masked_inner_avg = np.nanmean(u_eastward_bottom_fullrun_masked_inner, axis=(1,2))
v_northward_bottom_fullrun_masked_inner_avg = np.nanmean(v_northward_bottom_fullrun_masked_inner, axis=(1,2))
bottom_cur_mag_fullrun_masked_inner_avg = np.nanmean(bottom_cur_mag_fullrun_masked_inner, axis=(1,2))
# Mid
Hsig_fullrun_masked_mid_avg = np.nanmean(Hsig_fullrun_masked_mid, axis=(1,2))
Uwave_rms_fullrun_masked_mid_avg = np.nanmean(Uwave_rms_fullrun_masked_mid, axis=(1,2))
ubar_eastward_fullrun_masked_mid_avg = np.nanmean(ubar_eastward_fullrun_masked_mid, axis=(1,2))
vbar_northward_fullrun_masked_mid_avg = np.nanmean(vbar_northward_fullrun_masked_mid, axis=(1,2))
depth_int_sscflux_allsed_mag_fullrun_masked_mid_avg = np.nanmean(depth_int_sscflux_allsed_mag_fullrun_masked_mid, axis=(1,2))
bstrcwmax_fullrun_masked_mid_avg = np.nanmean(bstrcwmax_fullrun_masked_mid, axis=(1,2))
u_eastward_bottom_fullrun_masked_mid_avg = np.nanmean(u_eastward_bottom_fullrun_masked_mid, axis=(1,2))
v_northward_bottom_fullrun_masked_mid_avg = np.nanmean(v_northward_bottom_fullrun_masked_mid, axis=(1,2))
bottom_cur_mag_fullrun_masked_mid_avg = np.nanmean(bottom_cur_mag_fullrun_masked_mid, axis=(1,2))
# Outer
Hsig_fullrun_masked_outer_avg = np.nanmean(Hsig_fullrun_masked_outer, axis=(1,2))
Uwave_rms_fullrun_masked_outer_avg = np.nanmean(Uwave_rms_fullrun_masked_outer, axis=(1,2))
ubar_eastward_fullrun_masked_outer_avg = np.nanmean(ubar_eastward_fullrun_masked_outer, axis=(1,2))
vbar_northward_fullrun_masked_outer_avg = np.nanmean(vbar_northward_fullrun_masked_outer, axis=(1,2))
depth_int_sscflux_allsed_mag_fullrun_masked_outer_avg = np.nanmean(depth_int_sscflux_allsed_mag_fullrun_masked_outer, axis=(1,2))
bstrcwmax_fullrun_masked_outer_avg = np.nanmean(bstrcwmax_fullrun_masked_outer, axis=(1,2))
u_eastward_bottom_fullrun_masked_outer_avg = np.nanmean(u_eastward_bottom_fullrun_masked_outer, axis=(1,2))
v_northward_bottom_fullrun_masked_outer_avg = np.nanmean(v_northward_bottom_fullrun_masked_outer, axis=(1,2))
bottom_cur_mag_fullrun_masked_outer_avg = np.nanmean(bottom_cur_mag_fullrun_masked_outer, axis=(1,2))



# Play with plotting! 


# ------------------------------------------------------------------------------
# ----------- Plot 1: Scatter Plot Bed Stress, Eastward Currents ---------------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) bed shear stress on
# the x-axis and the spatially-averaged (masked) eastward currents on 
# the y-axis 

# Make the figure 
fig1, ax1 = plt.subplots(figsize=(10,10))
ax1.scatter(bstrcwmax_fullrun_masked_avg, ubar_eastward_fullrun_masked_avg, color='darkturquoise',
            edgecolors='black')
ax1.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax1.set_xlabel('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax1.set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)
ax1.set_title('Spatially-Averaged Bed Shear Stress \nand Depth-Averaged Currents', fontsize=fontsize)



# ------------------------------------------------------------------------------
# ----------- Plot 2: Scatter Plot Hsig, East Current, SSC Flux ----------------
# ------------------------------------------------------------------------------
# Make a scatter plot with Hsig on the x-axis, depth-averaged eastward currents
# on the y-axis, and the colors indicating the magnitude of SSC flux 

# Make the figure 
# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly
# First trim to have same ending 
Hsig_fullrun_masked_avg_same_times = Hsig_fullrun_masked_avg[:-117]
# Take every 3rd time
Hsig_fullrun_masked_avg_same_times = Hsig_fullrun_masked_avg_same_times[::3]

cmap2 = cmocean.cm.turbid
fig2, ax2 = plt.subplots(figsize=(12,10))
sc1 = ax2.scatter(Hsig_fullrun_masked_avg_same_times, ubar_eastward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  cmap=cmap2, edgecolors='black')
cbar1 = plt.colorbar(sc1, extend='max').set_label(label='SSC Flux Magnitude (kg/ms)', fontsize=fontsize)
ax2.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax2.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax2.set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)
ax2.set_title('Spatially-Averaged', fontsize=fontsize)


# ------------------------------------------------------------------------------
# ----------- Plot 3: Scatter Plot Hsig, North Current, SSC Flux ----------------
# ------------------------------------------------------------------------------
# Make a scatter plot with Hsig on the x-axis, depth-averaged northward currents
# on the y-axis, and the colors indicating the magnitude of SSC flux 

# Make the figure 
# Note: Had to trim wave time since the model did not run until November 2 but 
# rather October 29 hour 4 
# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly

cmap3 = cmocean.cm.turbid
fig3, ax3 = plt.subplots(figsize=(12,10))
sc2 = ax3.scatter(Hsig_fullrun_masked_avg_same_times, vbar_northward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  cmap=cmap2, edgecolors='black')
cbar2 = plt.colorbar(sc2, extend='max').set_label(label='SSC Flux Magnitude (kg/ms)', fontsize=fontsize)
ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax3.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax3.set_ylabel('Depth-Averaged \nNorthward Currents (m/s)', fontsize=fontsize)
ax3.set_title('Spatially-Averaged', fontsize=fontsize)


# ------------------------------------------------------------------------------
# --------------- Plot 4: Scatter Plot Bed Stress, SSC Flux -------------------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) bed shear stress on
# the x-axis and the spatially-averaged (masked) sediment flux magntiude on 
# the y-axis 

# Make the figure 
fig4, ax4 = plt.subplots(figsize=(8,8))
ax4.scatter(bstrcwmax_fullrun_masked_avg, depth_int_sscflux_allsed_mag_fullrun_masked_avg, color='darkturquoise',
            edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax4.set_xlabel('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax4.set_ylabel('Depth-Integrated \nSSC Flux (kg/ms)', fontsize=fontsize)
ax4.set_title('Spatially-Averaged', fontsize=fontsize)


# ------------------------------------------------------------------------------
# --------------- Plot 5: Scatter Plot Bed Stress, SSC Flux -------------------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) bed shear stress on
# the x-axis and the spatially-averaged (masked) sediment flux magntiude on 
# the y-axis and the colors representing the current magnitude 

# Make the figure 
cmap5 = cmocean.cm.tempo
fig5, ax5 = plt.subplots(figsize=(8,8))
sc5 = ax5.scatter(bstrcwmax_fullrun_masked_avg, depth_int_sscflux_allsed_mag_fullrun_masked_avg,
            c=depth_avg_cur_mag_fullrun_masked_avg, cmap=cmap5, edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
cbar5 = plt.colorbar(sc5, extend='max').set_label(label='Depth-Avg Current Magnitude (m/s)', fontsize=fontsize)
ax5.set_xlabel('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax5.set_ylabel('Depth-Integrated \nSSC Flux (kg/ms)', fontsize=fontsize)
ax5.set_title('Spatially-Averaged', fontsize=fontsize)


# ----------------------------------------------------------------------------------
# -- Plot 6: Scatter Plot Hsig, East Current, SSC Flux for Whole, East, and West --
# ----------------------------------------------------------------------------------
# Make a scatter plot with Hsig on the x-axis, depth-averaged eastward currents
# on the y-axis, and the colors indicating the magnitude of SSC flux 

# Make the figure 
# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly
# First trim to have same ending 
Hsig_fullrun_masked_east_avg_same_times = Hsig_fullrun_masked_east_avg[:-117]
Hsig_fullrun_masked_west_avg_same_times = Hsig_fullrun_masked_west_avg[:-117]
# Take every 3rd time
Hsig_fullrun_masked_east_avg_same_times = Hsig_fullrun_masked_east_avg_same_times[::3]
Hsig_fullrun_masked_west_avg_same_times = Hsig_fullrun_masked_west_avg_same_times[::3]

cmap6 = cmocean.cm.turbid
fig6, ax6 = plt.subplots(figsize=(12,10))
sc6 = ax6.scatter(Hsig_fullrun_masked_avg_same_times, ubar_eastward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  cmap=cmap2, marker='.', label='Whole Grid')
sc7 = ax6.scatter(Hsig_fullrun_masked_east_avg_same_times, ubar_eastward_fullrun_masked_east_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_east_avg,
                  cmap=cmap2, marker='x', label='East')
sc8 = ax6.scatter(Hsig_fullrun_masked_west_avg_same_times, ubar_eastward_fullrun_masked_west_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_west_avg,
                  cmap=cmap2, marker='+', label='West')
cbar6 = plt.colorbar(sc6, extend='max').set_label(label='SSC Flux Magnitude (kg/ms)', fontsize=fontsize)
ax6.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax6.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax6.set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)
ax6.set_title('Spatially-Averaged', fontsize=fontsize)
ax6.legend(fontsize=fontsize)


# ----------------------------------------------------------------------------------
# -- Plot 7: Scatter Plot Hsig, East Current, SSC Flux for Whole, East, and West  --
# ------------------------------------- Subplots -----------------------------------
# ----------------------------------------------------------------------------------
# Same as above but split into different subplots

# Make the figure 
# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly 
cmap7 = cmocean.cm.turbid
fig7, ax7 = plt.subplots(1, 3, figsize=(20,8))
gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
# Whole Grid
ax7_0 = plt.subplot(gs[0])
sc9 = ax7_0.scatter(Hsig_fullrun_masked_avg_same_times, ubar_eastward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  cmap=cmap7, marker='o', s=50, edgecolors='black', label='Whole Grid',
                  vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
ax7_0.set_ylim(-0.65, 0.4)
ax7_0.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax7_0.set_title('Whole Grid', fontsize=fontsize)
ax7_0.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax7_0.set_ylabel('Depth-Averaged Eastward \nCurrents (m/s)', fontsize=fontsize)

# East
ax7_1 = plt.subplot(gs[1], sharey=ax7_0)
sc10 = ax7_1.scatter(Hsig_fullrun_masked_east_avg_same_times, ubar_eastward_fullrun_masked_east_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_east_avg,
                  cmap=cmap7, marker='D', s=50, edgecolors='black', label='East',
                  vmin=0.0, vmax=0.4)
ax7_1.set_ylim(-0.65, 0.4)
ax7_1.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax7_1.set_title('East Half', fontsize=fontsize)
ax7_1.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
plt.setp(ax7_1.get_yticklabels(), visible=False)
#ax7[1].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# West 
ax7_2 = plt.subplot(gs[2], sharey=ax7_0)
sc11 = ax7_2.scatter(Hsig_fullrun_masked_west_avg_same_times, ubar_eastward_fullrun_masked_west_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_west_avg,
                  cmap=cmap7, marker='s', s=50, edgecolors='black', label='West',
                  vmin=0.0, vmax=0.4)
ax7_2.set_ylim(-0.65, 0.4)
ax7_2.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax7_2.set_title('West Half', fontsize=fontsize)
ax7_2.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
plt.setp(ax7_2.get_yticklabels(), visible=False)
#ax7[2].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# Set the colorbar
fig7.subplots_adjust(right=0.83)
cbar_ax7 = fig7.add_axes([0.85, 0.15, 0.02, 0.7])
cbar7 = plt.colorbar(sc11, extend='max', ax=[ax7_0, ax7_1, ax7_2], cax=cbar_ax7, pad=0.025).set_label(label='SSC Flux Magnitude (kg/ms)', fontsize=fontsize)
#ax7.legend(fontsize=fontsize)


plt.subplots_adjust(wspace=0.05) # 0.1



# ------------------------------------------------------------------------------
# ----------- Plot 8: Scatter Plot Eastward & Nortward  Currents ---------------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) eastward currents on
# the x-axis and the spatially-averaged (masked) northward currents on 
# the y-axis 

# Make the figure 
fig8, ax8 = plt.subplots(figsize=(10,10))
ax8.scatter(ubar_eastward_fullrun_masked_avg, vbar_northward_fullrun_masked_avg, color='darkturquoise',
            edgecolors='black')
ax8.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax8.axvline(x=0.0, linestyle='--', color='k', linewidth=2)
ax8.set_xlabel('Depth-Averaged Eastward Currents (m/s)', fontsize=fontsize)
ax8.set_ylabel('Depth-Averaged Northward Currents (m/s)', fontsize=fontsize)
ax8.set_title('Spatially-Averaged Depth-Averaged Currents', fontsize=fontsize)


# ------------------------------------------------------------------------------
# ----------- Plot 9: Scatter Plot Eastward & Nortward  Currents ---------------
# ------------------- Whole, East, West Subplots -------------------------------
# ------------------------------------------------------------------------------
# Make the same plot as above with ubar on x-axis and vbar on y-axis but for 
# just certain regions of the grid 

# Make the figure 
cmap9 = cmocean.cm.turbid
fig9, ax9 = plt.subplots(1, 3, figsize=(20,8))
gs2 = gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
# Whole Grid
ax9_0 = plt.subplot(gs2[0])
sc12 = ax9_0.scatter(ubar_eastward_fullrun_masked_avg, vbar_northward_fullrun_masked_avg,
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  cmap=cmap9, marker='o', s=50, edgecolors='black', label='Whole Grid',
                  vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
ax9_0.set_ylim(-0.2, 0.3)
ax9_0.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax9_0.axvline(x=0.0, linestyle='--', color='k', linewidth=2)
ax9_0.set_title('Whole Grid', fontsize=fontsize)
ax9_0.set_xlabel('Depth-Averaged Eastward \nCurrents (m/s)', fontsize=fontsize)
ax9_0.set_ylabel('Depth-Averaged Northward \nCurrents (m/s)', fontsize=fontsize)

# East
ax9_1 = plt.subplot(gs2[1], sharey=ax9_0)
sc13 = ax9_1.scatter(ubar_eastward_fullrun_masked_east_avg, vbar_northward_fullrun_masked_east_avg,
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_east_avg,
                  cmap=cmap9, marker='D', s=50, edgecolors='black', label='East',
                  vmin=0.0, vmax=0.4)
ax9_1.set_ylim(-0.2, 0.3)
ax9_1.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax9_1.axvline(x=0.0, linestyle='--', color='k', linewidth=2)
ax9_1.set_title('East Half', fontsize=fontsize)
ax9_1.set_xlabel('Depth-Averaged Eastward \nCurrents (m/s)', fontsize=fontsize)
plt.setp(ax9_1.get_yticklabels(), visible=False)
#ax7[1].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# West 
ax9_2 = plt.subplot(gs2[2], sharey=ax9_0)
sc14 = ax9_2.scatter(ubar_eastward_fullrun_masked_west_avg, vbar_northward_fullrun_masked_west_avg,
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_west_avg,
                  cmap=cmap9, marker='s', s=50, edgecolors='black', label='West',
                  vmin=0.0, vmax=0.4)
ax9_2.set_ylim(-0.2, 0.3)
ax9_2.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax9_2.axvline(x=0.0, linestyle='--', color='k', linewidth=2)
ax9_2.set_title('West Half', fontsize=fontsize)
ax9_2.set_xlabel('Depth-Averaged Eastward \nCurrents (m/s)', fontsize=fontsize)
plt.setp(ax9_2.get_yticklabels(), visible=False)
#ax7[2].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# Set the colorbar
fig9.subplots_adjust(right=0.83)
cbar_ax9 = fig9.add_axes([0.85, 0.15, 0.02, 0.7])
cbar9 = plt.colorbar(sc14, extend='max', ax=[ax9_0, ax9_1, ax9_2], cax=cbar_ax9, pad=0.025).set_label(label='SSC Flux Magnitude (kg/ms)', fontsize=fontsize)
#ax7.legend(fontsize=fontsize)

plt.subplots_adjust(wspace=0.05) # 0.1


# ------------------------------------------------------------------------------
# ------------ Plot 10: Scatter Plot Hsig, Cur Mag, Bed Stress -----------------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) significant wave height on
# the x-axis and the spatially-averaged (masked) depth-averaged current magntiude on 
# the y-axis and the colors representing the bed shear stress
# (It would make more sense to do bottom currents instead of depth-averaged but
# start here for now?)

# Make the figure 
cmap10 = cmocean.cm.tempo
fig10, ax10 = plt.subplots(figsize=(8,8))
sc15 = ax10.scatter(Hsig_fullrun_masked_avg_same_times, depth_avg_cur_mag_fullrun_masked_avg,
            c=bstrcwmax_fullrun_masked_avg, cmap=cmap10, edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
cbar15 = plt.colorbar(sc15, extend='max').set_label(label='Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax10.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax10.set_ylabel('Depth-Averaged Current \nMagnitude (m/s)', fontsize=fontsize)
ax10.set_title('Spatially-Averaged Whole Grid', fontsize=fontsize)


# ------------------------------------------------------------------------------
# -------- Plot 11: Scatter Plot Bottom Wave Orbital Velocity, Bottom ----------
# ----------------------- Current Magnitude, Bed Stress ------------------------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) bottom wave orbital
# velocity on the x-axis and the spatially-averaged (masked) bottom current 
# magntiude on the y-axis and the colors representing the bed shear stress

# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly
# First trim to have same ending 
Uwave_rms_fullrun_masked_avg_same_times = Uwave_rms_fullrun_masked_avg[:-117]
# Take every 3rd time
Uwave_rms_fullrun_masked_avg_same_times = Uwave_rms_fullrun_masked_avg_same_times[::3]

# Make the figure 
cmap11 = cmocean.cm.tempo
fig11, ax11 = plt.subplots(figsize=(8,8))
sc16 = ax11.scatter(Uwave_rms_fullrun_masked_avg_same_times, bottom_cur_mag_fullrun_masked_avg,
            c=bstrcwmax_fullrun_masked_avg, cmap=cmap11, edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
cbar16 = plt.colorbar(sc16, extend='max').set_label(label='Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax11.set_xlabel('Bottom Wave Orbital \nVelocity (m/s)', fontsize=fontsize)
ax11.set_ylabel('Bottom Current \nMagnitude (m/s)', fontsize=fontsize)
ax11.set_title('Spatially-Averaged Whole Grid', fontsize=fontsize)


# ----------------------------------------------------------------------------------
# -- Plot 12: Scatter Plot Hsig, East Current, SSC Flux for Inner, Mid, Outer  --
# ------------------------------------- Subplots -----------------------------------
# ----------------------------------------------------------------------------------
# Scatter plot of wave height, eastward currents, and SSC flux for the inner,
# mid, and outer shelf 

# Make the figure 
# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly
# First trim to have same ending 
Hsig_fullrun_masked_inner_avg_same_times = Hsig_fullrun_masked_inner_avg[:-117]
Hsig_fullrun_masked_mid_avg_same_times = Hsig_fullrun_masked_mid_avg[:-117]
Hsig_fullrun_masked_outer_avg_same_times = Hsig_fullrun_masked_outer_avg[:-117]
# Take every 3rd time
Hsig_fullrun_masked_inner_avg_same_times = Hsig_fullrun_masked_inner_avg_same_times[::3]
Hsig_fullrun_masked_mid_avg_same_times = Hsig_fullrun_masked_mid_avg_same_times[::3]
Hsig_fullrun_masked_outer_avg_same_times = Hsig_fullrun_masked_outer_avg_same_times[::3]


cmap12 = cmocean.cm.turbid #turbid
fig12, ax12 = plt.subplots(1, 3, figsize=(20,8)) #(20,8)
gs3 = gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
# Whole Grid
ax12_0 = plt.subplot(gs3[0])
sc17 = ax12_0.scatter(Hsig_fullrun_masked_inner_avg_same_times, ubar_eastward_fullrun_masked_inner_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_inner_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e0),
                  cmap=cmap12, marker='o', s=50, edgecolors='black', label='Whole Grid') # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
ax12_0.set_ylim(-0.4, 0.4)
ax12_0.set_xlim(0, 1.0)
ax12_0.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax12_0.set_title('Inner: 5 - 20 m', fontsize=fontsize)
ax12_0.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax12_0.set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# East
ax12_1 = plt.subplot(gs3[1], sharey=ax12_0)
sc18 = ax12_1.scatter(Hsig_fullrun_masked_mid_avg_same_times, ubar_eastward_fullrun_masked_mid_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_mid_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e0),
                  cmap=cmap12, marker='D', s=50, edgecolors='black', label='East')
ax12_1.set_ylim(-0.4, 0.4)
ax12_1.set_xlim(0, 1.0)
ax12_1.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax12_1.set_title('Mid: 20 - 40 m', fontsize=fontsize)
ax12_1.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
plt.setp(ax12_1.get_yticklabels(), visible=False)
#ax7[1].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# West 
ax12_2 = plt.subplot(gs3[2], sharey=ax12_0)
sc19 = ax12_2.scatter(Hsig_fullrun_masked_outer_avg_same_times, ubar_eastward_fullrun_masked_outer_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_outer_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-3, vmax=1e0),
                  cmap=cmap12, marker='s', s=50, edgecolors='black', label='West')
ax12_2.set_ylim(-0.4, 0.4)
ax12_2.set_xlim(0, 1.0)
ax12_2.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax12_2.set_title('Outer: 40 - 60 m', fontsize=fontsize)
ax12_2.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
plt.setp(ax12_2.get_yticklabels(), visible=False)
#ax7[2].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# Set the colorbar
fig12.subplots_adjust(right=0.83)
cbar_ax12 = fig12.add_axes([0.85, 0.15, 0.02, 0.7])
cbar12 = plt.colorbar(sc19, extend='max', ax=[ax12_0, ax12_1, ax12_2], cax=cbar_ax12, pad=0.025).set_label(label='Suspended Sediment \nFlux Magnitude (kg/ms)', fontsize=fontsize)
#ax7.legend(fontsize=fontsize)

# Label subplots
plt.text(0.327, 0.840, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.565, 0.840, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.807, 0.840, 'c)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)


plt.subplots_adjust(wspace=0.05) # 0.1, 0.05


# Do some quick caluclations of percent of time with eastward/westward currents 
# in mid shelf and the associated percent of sediment flux 
# Get the total flux
total_ssflux_mid = depth_int_sscflux_allsed_mag_fullrun_masked_mid_avg.sum()

# Get the eastward and westward indices 
eastward_idx = np.where(ubar_eastward_fullrun_masked_mid_avg>0)
westward_idx = np.where(ubar_eastward_fullrun_masked_mid_avg<0)

# Get the number of events 
east_cnt = len(eastward_idx[0][:])
west_cnt = len(westward_idx[0][:])

# Get eastward and westward total fluxes
eastward_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_mid_avg[eastward_idx[0][:]].sum()
westward_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_mid_avg[westward_idx[0][:]].sum()

# Find the fraction of time of each 
frac_time_east = east_cnt/(east_cnt+west_cnt)
frac_time_west = west_cnt/(east_cnt+west_cnt)

# Find the fraction of ssflux from each
frac_ssflux_east = eastward_ssflux/total_ssflux_mid
frac_ssflux_west = westward_ssflux/total_ssflux_mid

# Print results:
print('frac_time_east: ', frac_time_east)
print('frac_time_west: ', frac_time_west)
print('frac_ssflux_east: ', frac_ssflux_east)
print('frac_ssflux_west: ', frac_ssflux_west)



# ------------------------------------------------------------------------------
# -------- Plot 13: Scatter Plot Bottom Wave Orbital Velocity, Bottom ----------
# ----- Current Magnitude, Bed Stress for Innter, Mid, and Outer Shelf ---------
# ------------------------------------------------------------------------------
# Make a scatter plot with the spatially-averaged (masked) bottom wave orbital
# velocity on the x-axis and the spatially-averaged (masked) bottom current 
# magntiude on the y-axis and the colors representing the bed shear stress

# STILL FIXING THIS PLOT AND NEED TO ADD CODE ABOVE 

# New note for 2020 version - need to take the wave variables on the same time 
# as the model output since the wave data is hourly and the output is 3 hourly
# First trim to have same ending 
Uwave_rms_fullrun_masked_inner_avg_same_times = Uwave_rms_fullrun_masked_inner_avg[:-117]
Uwave_rms_fullrun_masked_mid_avg_same_times = Uwave_rms_fullrun_masked_mid_avg[:-117]
Uwave_rms_fullrun_masked_outer_avg_same_times = Uwave_rms_fullrun_masked_outer_avg[:-117]
# Take every 3rd time
Uwave_rms_fullrun_masked_inner_avg_same_times = Uwave_rms_fullrun_masked_inner_avg_same_times[::3]
Uwave_rms_fullrun_masked_mid_avg_same_times = Uwave_rms_fullrun_masked_mid_avg_same_times[::3]
Uwave_rms_fullrun_masked_outer_avg_same_times = Uwave_rms_fullrun_masked_outer_avg_same_times[::3]

# Make the figure 
cmap13 = cmocean.cm.tempo
fig13, ax13 = plt.subplots(figsize=(8,8))
sc20 = ax13.scatter(Uwave_rms_fullrun_masked_inner_avg_same_times, bottom_cur_mag_fullrun_masked_inner_avg,
            c=bstrcwmax_fullrun_masked_inner_avg, cmap=cmap11, edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
cbar16 = plt.colorbar(sc16, extend='max').set_label(label='Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax11.set_xlabel('Bottom Wave Orbital \nVelocity (m/s)', fontsize=fontsize)
ax11.set_ylabel('Bottom Current \nMagnitude (m/s)', fontsize=fontsize)

fig13, ax13 = plt.subplots(1, 3, figsize=(20,8))
gs4 = gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
# Whole Grid
ax13_0 = plt.subplot(gs4[0])
sc20 = ax13_0.scatter(Uwave_rms_fullrun_masked_inner_avg_same_times, bottom_cur_mag_fullrun_masked_inner_avg,
                      c=bstrcwmax_fullrun_masked_inner_avg,
                  cmap=cmap13, marker='o', s=50, edgecolors='black', label='Whole Grid',
                  vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
#ax13_0.set_ylim(-0.65, 0.4)
ax13_0.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax13_0.set_title('Inner: 5 - 20 m', fontsize=fontsize)
ax13_0.set_xlabel('Bottom Wave Orbital \nVelocity (m/s)', fontsize=fontsize)
ax13_0.set_ylabel('Bottom Current \nMagnitude (m/s)', fontsize=fontsize)

# East
ax13_1 = plt.subplot(gs4[1], sharey=ax13_0)
sc21 = ax13_1.scatter(Uwave_rms_fullrun_masked_mid_avg_same_times, bottom_cur_mag_fullrun_masked_mid_avg,
                      c=bstrcwmax_fullrun_masked_mid_avg,
                  cmap=cmap13, marker='D', s=50, edgecolors='black', label='East',
                  vmin=0.0, vmax=0.4)
#ax13_1.set_ylim(-0.65, 0.4)
ax13_1.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax13_1.set_title('Mid: 20 - 40 m', fontsize=fontsize)
ax13_1.set_xlabel('Bottom Wave Orbital \nVelocity (m/s)', fontsize=fontsize)
plt.setp(ax13_1.get_yticklabels(), visible=False)
#ax7[1].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# West 
ax13_2 = plt.subplot(gs4[2], sharey=ax13_0)
sc22 = ax13_2.scatter(Uwave_rms_fullrun_masked_outer_avg_same_times, bottom_cur_mag_fullrun_masked_outer_avg,
                      c=bstrcwmax_fullrun_masked_outer_avg,
                  cmap=cmap13, marker='s', s=50, edgecolors='black', label='West',
                  vmin=0.0, vmax=0.4)
#ax13_2.set_ylim(-0.65, 0.4)
ax13_2.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax13_2.set_title('Outer: 40 - 60 m', fontsize=fontsize)
ax13_2.set_xlabel('Bottom Wave Orbital \nVelocity (m/s)', fontsize=fontsize)
plt.setp(ax13_2.get_yticklabels(), visible=False)
#ax7[2].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# Set the colorbar
fig13.subplots_adjust(right=0.83)
cbar_ax13 = fig13.add_axes([0.85, 0.15, 0.02, 0.7])
cbar13 = plt.colorbar(sc22, extend='max', ax=[ax13_0, ax13_1, ax13_2], cax=cbar_ax13, pad=0.025).set_label(label='Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
#ax7.legend(fontsize=fontsize)


plt.subplots_adjust(wspace=0.05) # 0.1



# ----------------------------------------------------------------------------------
# ------ Plot 14: Scatter Plot Hsig, East Current, SSC Flux for Whole Grid   -------
# ------------------------------------- Log Scale -----------------------------------
# ----------------------------------------------------------------------------------
# Same as above but split into different subplots

# Make the figure 
# Note: Had to trim wave time since the model did not run until November 2 but 
# rather October 29 hour 4 
cmap14 = cmocean.cm.turbid
fig14, ax14 = plt.subplots(figsize=(12,10))
# Whole Grid
sc23 = ax14.scatter(Hsig_fullrun_masked_avg_same_times, ubar_eastward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-2, vmax=1e0), cmap=cmap14, marker='o', s=50, edgecolors='black', label='Whole Grid')
                  #vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
ax14.set_ylim(-0.65, 0.4)
ax14.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
#ax14.set_title('Whole Grid', fontsize=fontsize)
ax14.set_xlabel('Significant Wave Height (m)', fontsize=fontsize)
ax14.set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# Set the colorbar
# (10e-2, 10e0)
cbar14 = plt.colorbar(sc23, extend='max').set_label(label='Suspended Sediment \nFlux Magnitude (kg/ms)', fontsize=fontsize)



# ----------------------------------------------------------------------------------
# ------ Plot 15: Scatter Plot Bstrcwmax, East Current, SSC Flux for Whole Grid   -------
# ------------------------------------- Log Scale -----------------------------------
# ----------------------------------------------------------------------------------

# Make the figure 
# Note: Had to trim wave time since the model did not run until November 2 but 
# rather October 29 hour 4 
cmap15 = cmocean.cm.turbid
fig15, ax15 = plt.subplots(figsize=(12,10))
# Whole Grid
sc24 = ax15.scatter(bstrcwmax_fullrun_masked_avg, ubar_eastward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-2, vmax=1e0), cmap=cmap14, marker='o', s=50, edgecolors='black', label='Whole Grid')
                  #vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
#ax15.set_ylim(-0.65, 0.4)
ax15.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax15.set_title('Whole Grid', fontsize=fontsize)
ax15.set_xlabel('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax15.set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)

# Set the colorbar
# (10e-2, 10e0)
cbar15 = plt.colorbar(sc24, extend='max').set_label(label='Suspended Sediment \nFlux Magnitude (kg/ms)', fontsize=fontsize)



# ----------------------------------------------------------------------------------
# ------ Plot 16: Subplot of two plots together: Scatter Plot Bstrcwmax, -----------
# ----------- East Current, SSC Flux for Whole Grid and wave orbitals, -------------
# ------------------------------ current magnitude, and bed stress   ---------------
# ----------------------------------------------------------------------------------

# Make the figure 
# Note: Had to trim wave time since the model did not run until November 2 but 
# rather October 29 hour 4 
cmap16 = cmocean.cm.tempo
cmap17 = cmocean.cm.turbid
fig16, ax16 = plt.subplots(1, 2, figsize=(16,9))
# Whole Grid
# Bottom wave orbital velocity and bottom current magnitude with bed shear 
# stress as shading 
sc25 = ax16[0].scatter(Uwave_rms_fullrun_masked_avg_same_times, bottom_cur_mag_fullrun_masked_avg,
            c=bstrcwmax_fullrun_masked_avg, cmap=cmap16, edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
cbar16 = plt.colorbar(sc25, extend='max', ax=ax16[0], orientation='horizontal', pad=0.18).set_label(label='Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax16[0].set_xlabel('Bottom Wave Orbital \nVelocity (m/s)', fontsize=fontsize)
ax16[0].set_ylabel('Bottom Current Magnitude \n(m/s)', fontsize=fontsize)


# Bed shear stress, depth-avg eastward current, ss flux as shading 
sc26 = ax16[1].scatter(bstrcwmax_fullrun_masked_avg, ubar_eastward_fullrun_masked_avg, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-2, vmax=1e0), cmap=cmap17, marker='o', s=50, edgecolors='black', label='Whole Grid')
                  #vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
#ax15.set_ylim(-0.65, 0.4)
ax16[1].axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax16[1].set_xlabel('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax16[1].set_ylabel('Depth-Averaged \nEastward Currents (m/s)', fontsize=fontsize)
# Set the colorbar
# (10e-2, 10e0)
cbar17 = plt.colorbar(sc26, extend='max', ax=ax16[1], orientation='horizontal', pad=0.18).set_label(label='Suspended Sediment \nFlux Magnitude (kg/ms)', fontsize=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.1, wspace=0.35) #0.08

# Add subplot labels
plt.text(0.429, 0.846, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.873, 0.846, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# Same as above but slightly different format
# ----------------------------------------------------------------------------------
# ------ Plot 17: Subplot of two plots together: Scatter Plot Bstrcwmax, -----------
# ----------- East Current, SSC Flux for Whole Grid and wave orbitals, -------------
# ------------------------------ current magnitude, and bed stress   ---------------
# ----------------------------------------------------------------------------------

# Make the figure 
# Note: Had to trim wave time since the model did not run until November 2 but 
# rather October 29 hour 4 
cmap16 = cmocean.cm.tempo
cmap17 = cmocean.cm.turbid
cmap17.set_under(cmap17(0))
fig17, ax17 = plt.subplots(1, 2, figsize=(16,6)) # (16,9)
# Whole Grid
# Bottom wave orbital velocity and bottom current magnitude with bed shear 
# stress as shading 
sc25 = ax17[0].scatter(Uwave_rms_fullrun_masked_avg_same_times*100, bottom_cur_mag_fullrun_masked_avg*100,
            c=bstrcwmax_fullrun_masked_avg, cmap=cmap16, edgecolors='black')
#ax3.axhline(y=0.0, linestyle='--', color='k', linewidth=2)
#cax17 = fig17.add_axes([0.13, .80, 0.3, 0.03])
#fig17.colorbar(surf, orientation='horizontal', cax=cax)
cax17 = fig17.add_axes([0.13, .96, 0.29, 0.03])
cbar17 = plt.colorbar(sc25, extend='max', cax=cax17, orientation='horizontal', pad=0.02) #.set_label(label='Bed Shear Stress (N/m\u00b2)', fontsize=fontsize, labelpad=8)
ax17[0].set_xlabel('Bottom Wave Orbital \nVelocity (cm/s)', fontsize=fontsize)
ax17[0].set_ylabel('Bottom \nCurrent \nMagnitude \n(cm/s)', fontsize=fontsize, rotation=0, labelpad=60)
ax17[0].set_title('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize, pad=60)


# Bed shear stress, depth-avg eastward current, ss flux as shading 
#depth_int_sscflux_allsed_mag_fullrun_masked_avg_min_adjust = np.where(depth_int_sscflux_allsed_mag_fullrun_masked_avg == 0, 1e-2, depth_int_sscflux_allsed_mag_fullrun_masked_avg)
sc26 = ax17[1].scatter(bstrcwmax_fullrun_masked_avg, ubar_eastward_fullrun_masked_avg*100, 
                  c=depth_int_sscflux_allsed_mag_fullrun_masked_avg,
                  norm=matplotlib.colors.LogNorm(vmin=1e-2, vmax=1e0), cmap=cmap17, marker='o', s=50, edgecolors='black', label='Whole Grid')
                  #vmin=0.0, vmax=0.4) # s=depth_int_sscflux_allsed_mag_fullrun_masked_avg*500
#ax15.set_ylim(-0.65, 0.4)
ax17[1].axhline(y=0.0, linestyle='--', color='k', linewidth=2)
ax17[1].set_xlabel('Bed Shear Stress (N/m\u00b2)', fontsize=fontsize)
ax17[1].set_ylabel('Depth \nAveraged \nEastward \nCurrents \n(cm/s)', fontsize=fontsize, rotation=0, labelpad=50)
# Set the colorbar
# (10e-2, 10e0)
cax18 = fig17.add_axes([0.605, .96, 0.29, 0.03])
cbar18 = plt.colorbar(sc26, extend='both', cax=cax18, orientation='horizontal', pad=0.03, ticklocation='bottom') #.set_label(label='Suspended Sediment Flux Magnitude (kg/ms)', fontsize=fontsize, labelpad=8)
ax17[1].set_title('Suspended Sediment Flux Magnitude (kg/ms)', fontsize=fontsize, pad=60)
# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.1, wspace=0.55) #0.08

# Add subplot labels
plt.text(0.402, 0.817, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.870, 0.817, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)







# Do some quick caluclations of percent of time with eastward/westward currents 
# in the whole shelf 
# Get the total flux
total_ssflux_whole = depth_int_sscflux_allsed_mag_fullrun_masked_avg.sum()

# Get the percent of time that the bottom current mangitude is greater than 0.1 m/s
# Find the indices
bottom_cur_mag_idx_morethan10 = np.where(bottom_cur_mag_fullrun_masked_avg>0.1)
# Get the number of events
bottom_cur_mag_morethan10_cnt = len(bottom_cur_mag_idx_morethan10[0][:])
# Get the total fluxes during these times
bottom_cur_mag_morethan10_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_avg[bottom_cur_mag_idx_morethan10[0][:]].sum()
# Find the fraction of ssflux during these events compared to total
frac_ssflux_bottom_cur_mag_morethan10 = bottom_cur_mag_morethan10_ssflux/total_ssflux_whole
print('fraction of ssflux during bottom current magnitudes exceeding 10 cm/s: ', frac_ssflux_bottom_cur_mag_morethan10 )
# Find the fraction of time
frac_time_bottom_cur_mag_morethan10 = bottom_cur_mag_morethan10_cnt/full_time_len
print('fraction of time bottom current magnitude exceeds 10 cm/s: ', frac_time_bottom_cur_mag_morethan10)

# Get the eastward and westward indices 
eastward_idx = np.where(ubar_eastward_fullrun_masked_avg>0)
westward_idx = np.where(ubar_eastward_fullrun_masked_avg<0)
westward_idx_morethan10 = np.where(ubar_eastward_fullrun_masked_avg<-0.1)
westward_idx_morethan20 = np.where(ubar_eastward_fullrun_masked_avg<-0.2)

# Get the number of events 
east_cnt = len(eastward_idx[0][:])
west_cnt = len(westward_idx[0][:])
west_morethan10_cnt = len(westward_idx_morethan10[0][:])
west_morethan20_cnt = len(westward_idx_morethan20[0][:])

# Get eastward and westward total fluxes
eastward_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_avg[eastward_idx[0][:]].sum()
westward_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_avg[westward_idx[0][:]].sum()
westward_morethan10_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_avg[westward_idx_morethan10[0][:]].sum()
westward_morethan20_ssflux = depth_int_sscflux_allsed_mag_fullrun_masked_avg[westward_idx_morethan20[0][:]].sum()

# Find the fraction of time of each 
frac_time_east = east_cnt/(east_cnt+west_cnt)
frac_time_west = west_cnt/(east_cnt+west_cnt)
frac_time_west_morethan10 = west_morethan10_cnt/(east_cnt+west_cnt)
frac_time_west_morethan20 = west_morethan20_cnt/(east_cnt+west_cnt)

# Find the fraction of ssflux from each
frac_ssflux_east = eastward_ssflux/total_ssflux_whole
frac_ssflux_west = westward_ssflux/total_ssflux_whole
frac_ssflux_west_morethan10 = westward_morethan10_ssflux/total_ssflux_whole
frac_ssflux_west_morethan20 = westward_morethan20_ssflux/total_ssflux_whole

# Print results:
print('Whole Grid:')
print('frac_time_east: ', frac_time_east)
print('frac_time_west: ', frac_time_west)
print('frac_time_west_morethan10: ', frac_time_west_morethan10)
print('frac_time_west_morethan20: ', frac_time_west_morethan20)
print('frac_ssflux_east: ', frac_ssflux_east)
print('frac_ssflux_west: ', frac_ssflux_west)
print('frac_ssflux_west_morethan10: ', frac_ssflux_west_morethan10)
print('frac_ssflux_west_morethan20: ', frac_ssflux_west_morethan20)



# Get the fraction of time/indices when significant wave height exceeded 0.5 m on inner shelf
# Get the indices for greater than 0.5 m
# Slice first to fit 
Hsig_fullrun_masked_inner_avg_slice = Hsig_fullrun_masked_inner_avg_same_times
Hsig_fullrun_masked_inner_avg_slice_morethanp5_idx = np.where(Hsig_fullrun_masked_inner_avg_slice>0.1)

# Get the number of events 
Hsig_fullrun_masked_inner_avg_slice_morethanp5_cnt = len(Hsig_fullrun_masked_inner_avg_slice_morethanp5_idx[0][:])

# Find cnt of eastward and westward during these times
#Hsig_fullrun_masked_inner_avg_slice_morethanp5 = Hsig_fullrun_masked_inner_avg_slice[Hsig_fullrun_masked_inner_avg_slice_morethanp5_idx[0][:]]
ubar_eastward_fullrun_masked_inner_avg_morethanp5 = ubar_eastward_fullrun_masked_inner_avg[Hsig_fullrun_masked_inner_avg_slice_morethanp5_idx[0][:]]
eastward_depthavg_cur_morethanp5_idx = np.where(ubar_eastward_fullrun_masked_inner_avg_morethanp5>0)
westward_depthavg_cur_morethanp5_idx = np.where(ubar_eastward_fullrun_masked_inner_avg_morethanp5<0)

# Find the fraction of these times that currents were eastward or westward 
frac_time_eastward_morethanp5 = len(eastward_depthavg_cur_morethanp5_idx[0][:])/Hsig_fullrun_masked_inner_avg_slice_morethanp5_cnt
frac_time_westward_morethanp5 = len(westward_depthavg_cur_morethanp5_idx[0][:])/Hsig_fullrun_masked_inner_avg_slice_morethanp5_cnt

# Print results
print('fraction of time eastward depth-avg currents when swh > 0.1 m: ', frac_time_eastward_morethanp5)
print('fraction of time westward depth-avg currents when swh > 0.1 m: ', frac_time_westward_morethanp5)



# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Set up the data
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux = xr.Dataset(
    data_vars=dict(
        uwave_rms_spatial_avg=(['ocean_time'], Uwave_rms_fullrun_masked_avg_same_times),
        bottom_cur_mag_spatial_avg=(['ocean_time'], bottom_cur_mag_fullrun_masked_avg),
        bed_stress_spatial_avg=(['ocean_time'], bstrcwmax_fullrun_masked_avg),
        Hsig_spatial_avg=(['ocean_time'], Hsig_fullrun_masked_avg_same_times),
        depth_avg_u_current_spatial_avg=(['ocean_time'], ubar_eastward_fullrun_masked_avg),
        depth_int_ssflux_spatial_avg=(['ocean_time'], depth_int_sscflux_allsed_mag_fullrun_masked_avg),
        ),
    coords=dict(
        ocean_time = time_steps
        ),
    attrs=dict(description='Time-averaged ROMS output including depth-averaged, surface, and 1 meter above seafloor suspended sediment concentrations (mg/m3) and fluxes, and net deposition (cm)')) 
# Add more metadata?
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux.uwave_rms_spatial_avg.name='spatially-averaged bottom wave orbital velocity (m/s)'
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux.bottom_cur_mag_spatial_avg.name='spatially-averaged bottom current magnitude (m/s)'
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux.bed_stress_spatial_avg.name='spatially-averaged bed shear stress (N/m2)'
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux.Hsig_spatial_avg.name='spatially-averaged significant wave height (m)'
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux.depth_avg_u_current_spatial_avg.name='spatially-averaged, depth-averaged u current (m/s)'
roms_spatial_avg_bedstress_uwave_bot_cur_ssflux.depth_int_ssflux_spatial_avg.name='spatially-averaged, depth-integrated suspended sediment flux magnitude (kg/m/s)'

# Save to a netcdf
#roms_spatial_avg_bedstress_uwave_bot_cur_ssflux .to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Paper1_Take2/Data/fig9_roms_spatial_avg_uwave_cur_mag_bed_stress_hsig_ssflux.nc')










