############################# Plot Flow Integrated Transects ###########################
# The purpose of this script is to plot flow-integrated transects for different regions
# of the domain (cross-shore) for sediment fluxes. This script will then find the same
# plots for the sensitivity tests and plot the differences.
#
# Notes:
# - 
#########################################################################################


# Load in the packages
import numpy as np
import xarray as xr
import xesmf as xe
#import ESMF
from glob import glob
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator, LogLocator, LogFormatterExponent, LogFormatterSciNotation)
import cmocean
from matplotlib import colors

# Set a universal fontsize
fontsize = 25

# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize)
matplotlib.rc('ytick', labelsize=fontsize)

# Prevent tick labels from overlapping
matplotlib.rcParams['xtick.major.pad'] = 12
matplotlib.rcParams['ytick.major.pad'] = 12

# Load in the grid
grid = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Include/KakAKgrd_shelf_big010_smooth006.nc')
#grid = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Grids/KakAKgrd_shelf_big010_smooth006.nc')

# Pull out some dimensions
eta_rho_len = len(grid.eta_rho)
xi_rho_len = len(grid.xi_rho)

# Set number of vertical dimensions
s_rho_len = int(20) # Set manually 

# Load in the rho masks
mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
#mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
#mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')


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


# So we want to have flow/time-integrated sediment fluxes for different transects
# Should I pick the transects first, then do the math? or do it for all space
# then pull out the transects? I think picking transects first would run much 
# faster so maybe start there...

# ------------------- Define a Bunch of Functions -----------------------
# Define a function to pull out the lenght of time in the model run
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


# Make a function to loop through the output, pull out the time-integrated
# sediemnt flux for all sediment classes at 3 given transects 
def get_time_int_depth_int_ss_flux_allsed_3transects_cross_shore(filenames, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho):
    """
    The purpose of this function is to loop through model output and pull out the 
    total sediment flux at given transects, integrated it over time, then return it.
    
    Inputs:
    - filenames: Names/Paths of the model output
    - eta_rho_len: Length of eta_rho coordinates
    - s_rho_len: Length of s_rho coordinates
    - xi1: xi_rho coordinate of transect 1
    - xi2: xi_coordinate of transect 2
    - xi3: xi_coordinate of transect 3
    - regridder_u2rho: Weights to regrid from u points to rho points 
    Outputs:
    - time_int_ssflux_tran1: Time-integrated ss flux at transect 1 over all depths
    - time_int_ssflux_tran2: Time-integrated ss flux at transect 2 over all depths
    - time_int_ssflux_tran3: Time-integrated ss flux at transect 3 over all depths
    - time_int_depth_int_ssflux_tran1: Time-integrated, depth-integrated ss flux at transect 1 over all depths
    - time_int_depth_int_ssflux_tran2: Time-integrated, depth-integrated ss flux at transect 2 over all depths
    - time_int_depth_int_ssflux_tran3: Time-integrated, depth-integrated ss flux at transect 3 over all depths
    - time_avg_dz_tran1: Time-averaged depths for transect 1
    - time_avg_dz_tran2: Time-averaged depths for transect 2
    - time_avg_dz_tran3: Time-averaged depths for transect 3
    """
    
    # Find the number of files to loop through 
    num_files = len(filenames)
    
    # Pull out the length of time of the full run, the time steps, 
    # and the length of time of each output file
    full_time_len, time_steps, time_lengths = get_model_time(filenames, num_files)
    
    # Make empty arrays to hold the output
    # SS flux at transects for all time and depths at a given xi
    ssflux_tran1 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    ssflux_tran2 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    ssflux_tran3 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    # Depth
    dz_full = np.empty((full_time_len, s_rho_len, eta_rho_len, xi_rho_len))
    
    # Set the time step
    dt = 10800 # seconds (3 hours)
    
    # Set the time step
    time_step = 0 
    
    # Loop through the model output
    for j in range(num_files):
        print('j: ', j)
        
        # Open the model output
        model_output = xr.open_dataset(filenames[j])
        
        # Pull out u and interpolate it onto rho points
        u_tmp = model_output.u
        u_rho = regridder_u2rho(u_tmp)
        
        # Calculate the time-varying thickness of the cells
        dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
        
        # Pull out the total SSC
        tot_ssc = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
        
        # Get the values for each transect
        ssflux_tran1_tmp = tot_ssc[:,:,:,xi1] * u_rho[:,:,:,xi1]
        ssflux_tran2_tmp = tot_ssc[:,:,:,xi2] * u_rho[:,:,:,xi2]
        ssflux_tran3_tmp = tot_ssc[:,:,:,xi3] * u_rho[:,:,:,xi3]
        print('ssflux_tran1 shape: ', np.shape(ssflux_tran1))
        
        # Save these to arrays (need whole time serie before intergating)
        start = int(time_step)
        end = int(time_step+time_lengths[j])
        ssflux_tran1[start:end,:,:] = ssflux_tran1_tmp
        ssflux_tran2[start:end,:,:] = ssflux_tran2_tmp
        ssflux_tran3[start:end,:,:] = ssflux_tran3_tmp
        dz_full[start:end,:,:,:] = dz[:,:,:,:]
        
        # Update the time step
        time_step = time_step + time_lengths[j]
        
    # Once we have the full time series, we can integrated over time
    # Integrated over time
    time_int_ssflux_tran1 = np.sum(ssflux_tran1*dt, axis=0)
    time_int_ssflux_tran2 = np.sum(ssflux_tran2*dt, axis=0)
    time_int_ssflux_tran3 = np.sum(ssflux_tran3*dt, axis=0)
    print('time_int_ssflux_tran1_tmp shape: ', np.shape(time_int_ssflux_tran1))

    # Integrated over depth
    time_int_depth_int_ssflux_tran1 = np.sum(np.sum(ssflux_tran1*dt*dz[:,:,:,xi1], axis=0), axis=0)
    time_int_depth_int_ssflux_tran2 = np.sum(np.sum(ssflux_tran2*dt*dz[:,:,:,xi2], axis=0), axis=0)
    time_int_depth_int_ssflux_tran3 = np.sum(np.sum(ssflux_tran3*dt*dz[:,:,:,xi3], axis=0), axis=0)
    print('time_int_depth_int_ssflux_tran1_tmp shape: ', np.shape(time_int_depth_int_ssflux_tran1))
    
    # Find the time-averaged depth at each spot to use for plotting
    #time_avg_dz_tran1 = np.nanmean(dz[:,:,:,xi1], axis=0)
    #time_avg_dz_tran2 = np.nanmean(dz[:,:,:,xi2], axis=0)
    #time_avg_dz_tran3 = np.nanmean(dz[:,:,:,xi3], axis=0)
    
    # Return these values
    return(time_int_ssflux_tran1, time_int_ssflux_tran2, time_int_ssflux_tran3, time_int_depth_int_ssflux_tran1, time_int_depth_int_ssflux_tran2, time_int_depth_int_ssflux_tran3)


# Make a function that is similar to above but gets the time-integrated, depth-integrated
# flux for each sediment class 
def get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(filenames, sed_class, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho):
    """
    The purpose of this function is to get the time-integrated, depth-integrated
    suspended sediment (SS) flux for a given sediment class at three given cross-shore
    transects.
    
    Inputs:
    - filenames: Names/Paths of the model output
    - sed_class: String name of sediment class of interest (ex: 'mud_01')
    - eta_rho_len: Length of eta_rho coordinates
    - s_rho_len: Length of s_rho coordinates
    - xi1: xi_rho coordinate of transect 1
    - xi2: xi_coordinate of transect 2
    - xi3: xi_coordinate of transect 3
    - regridder_u2rho: Weights to regrid from u points to rho points 
    Outputs:
    - time_int_ssflux_1sed_tran1: Time-integrated ss flux at transect 1 over all depths
    - time_int_ssflux_1sed_tran2: Time-integrated ss flux at transect 2 over all depths
    - time_int_ssflux_1sed_tran3: Time-integrated ss flux at transect 3 over all depths
    - time_int_depth_int_ssflux_1sed_tran1: Time-integrated, depth-integrated ss flux at transect 1 over all depths
    - time_int_depth_int_ssflux_1sed_tran2: Time-integrated, depth-integrated ss flux at transect 2 over all depths
    - time_int_depth_int_ssflux_1sed_tran3: Time-integrated, depth-integrated ss flux at transect 3 over all depths
    """
    
    # Find the number of files to loop through 
    num_files = len(filenames)
    
    # Pull out the length of time of the full run, the time steps, 
    # and the length of time of each output file
    full_time_len, time_steps, time_lengths = get_model_time(filenames, num_files)
    
    # Make empty arrays to hold the output
    # SS flux at transects for all time and depths at a given xi
    ssflux_1sed_tran1 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    ssflux_1sed_tran2 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    ssflux_1sed_tran3 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    # Depth
    dz_full = np.empty((full_time_len, s_rho_len, eta_rho_len, xi_rho_len))
    
    # Set the time step
    dt = 10800 # seconds (3 hours)
    
    # Set the time step
    time_step = 0 
    
    # Loop through the model output
    for j in range(num_files):
        print('j: ', j)
        
        # Open the model output
        model_output = xr.open_dataset(filenames[j])
        
        # Pull out u and interpolate it onto rho points
        u_tmp = model_output.u
        u_rho = regridder_u2rho(u_tmp)
        
        # Calculate the time-varying thickness of the cells
        dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
        
        # Pull out the total SSC
        ssc_1sed = model_output[sed_class]
        
        # Get the values for each transect
        ssflux_1sed_tran1_tmp = ssc_1sed[:,:,:,xi1] * u_rho[:,:,:,xi1]
        ssflux_1sed_tran2_tmp = ssc_1sed[:,:,:,xi2] * u_rho[:,:,:,xi2]
        ssflux_1sed_tran3_tmp = ssc_1sed[:,:,:,xi3] * u_rho[:,:,:,xi3]
        print('ssflux_1sed_tran1 shape: ', np.shape(ssflux_1sed_tran1))
        
        # Save these to arrays (need whole time serie before intergating)
        start = int(time_step)
        end = int(time_step+time_lengths[j])
        ssflux_1sed_tran1[start:end,:,:] = ssflux_1sed_tran1_tmp
        ssflux_1sed_tran2[start:end,:,:] = ssflux_1sed_tran2_tmp
        ssflux_1sed_tran3[start:end,:,:] = ssflux_1sed_tran3_tmp
        dz_full[start:end,:,:,:] = dz[:,:,:,:]
        
        # Update the time step
        time_step = time_step + time_lengths[j]
        
    # Once we have the full time series, we can integrated over time
    # Integrated over time
    time_int_ssflux_1sed_tran1 = np.sum(ssflux_1sed_tran1*dt, axis=0)
    time_int_ssflux_1sed_tran2 = np.sum(ssflux_1sed_tran2*dt, axis=0)
    time_int_ssflux_1sed_tran3 = np.sum(ssflux_1sed_tran3*dt, axis=0)
    print('time_int_ssflux_1sed_tran1_tmp shape: ', np.shape(time_int_ssflux_1sed_tran1))

    # Integrated over depth
    time_int_depth_int_ssflux_1sed_tran1 = np.sum(np.sum(ssflux_1sed_tran1*dt*dz[:,:,:,xi1], axis=0), axis=0)
    time_int_depth_int_ssflux_1sed_tran2 = np.sum(np.sum(ssflux_1sed_tran2*dt*dz[:,:,:,xi2], axis=0), axis=0)
    time_int_depth_int_ssflux_1sed_tran3 = np.sum(np.sum(ssflux_1sed_tran3*dt*dz[:,:,:,xi3], axis=0), axis=0)
    print('time_int_depth_int_ssflux_1sed_tran1_tmp shape: ', np.shape(time_int_depth_int_ssflux_1sed_tran1))
    
    # Return these values
    return(time_int_ssflux_1sed_tran1, time_int_ssflux_1sed_tran2, time_int_ssflux_1sed_tran3, time_int_depth_int_ssflux_1sed_tran1, time_int_depth_int_ssflux_1sed_tran2, time_int_depth_int_ssflux_1sed_tran3)
        

    
# Make a function to get the time-averaged u currents at each transect 
# for a given run
def get_time_avg_u_currents_at_transects(filenames, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho):
    """
    The purpose of this function is to get the time-averaged u currents at 
    each transect for a given model run.
    
    Inputs:
    - filenames: Names/Paths of the model output
    - eta_rho_len: Length of eta_rho coordinates
    - s_rho_len: Length of s_rho coordinates
    - xi1: xi_rho coordinate of transect 1
    - xi2: xi_coordinate of transect 2
    - xi3: xi_coordinate of transect 3
    - regridder_u2rho: Weights to regrid from u points to rho points 
    Outputs:
    - time_avg_u_tran1: Time-averaged u velocity at transect 1
    - time_avg_u_tran2: Time-averaged u velocity at transect 2
    - time_avg_u_tran3: Time-averaged u velocity at transect 3
    
    """
    # Find the number of files to loop through 
    num_files = len(filenames)
    
    # Pull out the length of time of the full run, the time steps, 
    # and the length of time of each output file
    full_time_len, time_steps, time_lengths = get_model_time(filenames, num_files)
    
    # Make empty arrays to hold the output
    # SS flux at transects for all time and depths at a given xi
    u_tran1 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    u_tran2 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    u_tran3 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    # Depth
    dz_full = np.empty((full_time_len, s_rho_len, eta_rho_len, xi_rho_len))
    
    # Set the time step
    time_step = 0 
    
    # Loop through the model output
    for j in range(num_files):
        print('j: ', j)
        
        # Open the model output
        model_output = xr.open_dataset(filenames[j])
        
        # Pull out u and interpolate it onto rho points
        u_tmp = model_output.u
        u_rho = regridder_u2rho(u_tmp)
        
        # Calculate the time-varying thickness of the cells
        #dz = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
        
        # Get the values for each transect
        u_tran1_tmp = u_rho[:,:,:,xi1]
        u_tran2_tmp = u_rho[:,:,:,xi2]
        u_tran3_tmp = u_rho[:,:,:,xi3]
        print('u_tran1 shape: ', np.shape(u_tran1_tmp))
        
        # Save these to arrays (need whole time serie before intergating)
        start = int(time_step)
        end = int(time_step+time_lengths[j])
        u_tran1[start:end,:,:] = u_tran1_tmp
        u_tran2[start:end,:,:] = u_tran2_tmp
        u_tran3[start:end,:,:] = u_tran3_tmp
        #dz_full[start:end,:,:,:] = dz[:,:,:,:]
        
        # Update the time step
        time_step = time_step + time_lengths[j]
        
    # Once we have the full time series, we can get the averaged over time
    # Average over time
    time_avg_u_tran1 = np.nanmean(u_tran1, axis=0)
    time_avg_u_tran2 = np.nanmean(u_tran2, axis=0)
    time_avg_u_tran3 = np.nanmean(u_tran3, axis=0)
    print('time_avg_u_tran1_tmp shape: ', np.shape(time_avg_u_tran1))
    
    # Return these values
    return(time_avg_u_tran1, time_avg_u_tran2, time_avg_u_tran3)


# Make a function to get time series of z_rho at the transects
def get_z_rho_at_transects(filenames, eta_rho_len, s_rho_len, xi1, xi2, xi3):
    """
    The purpose of this function is to get a time series of the SSC at 
    each transect for a given model run.
    
    Inputs:
    - filenames: Names/Paths of the model output
    - eta_rho_len: Length of eta_rho coordinates
    - s_rho_len: Length of s_rho coordinates
    - xi1: xi_rho coordinate of transect 1
    - xi2: xi_coordinate of transect 2
    - xi3: xi_coordinate of transect 3 
    Outputs:
    - tot_ssc_tran1: Time series of ssc at transect 1
    - tot_ssc_tran2: Time series of ssc at transect 2
    - tot_ssc_tran3: Time series of ssc at transect 3
    
    """
    # Find the number of files to loop through 
    num_files = len(filenames)
    
    # Pull out the length of time of the full run, the time steps, 
    # and the length of time of each output file
    full_time_len, time_steps, time_lengths = get_model_time(filenames, num_files)
    
    # Make empty arrays to hold the output
    # Depth
    z_rho_tran1 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    z_rho_tran2 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    z_rho_tran3 = np.empty((full_time_len, s_rho_len, eta_rho_len))
    
    # Set the time step
    time_step = 0 
    
    # Loop through the model output
    for j in range(num_files):
        print('j: ', j)
        
        # Open the model output
        model_output = xr.open_dataset(filenames[j])
        
        # Get the z_rho
        z_rho_tmp = model_output.z_rho
        
        # Get the values for each transect
        z_rho_tran1_tmp = z_rho_tmp[:,:,:,xi1]
        z_rho_tran2_tmp = z_rho_tmp[:,:,:,xi2]
        z_rho_tran3_tmp = z_rho_tmp[:,:,:,xi3]
        print('z_rho_tran1_tmp shape: ', np.shape(z_rho_tran1_tmp))
        
        # Save these to arrays (need whole time serie before intergating)
        start = int(time_step)
        end = int(time_step+time_lengths[j])
        z_rho_tran1[start:end,:,:] = z_rho_tran1_tmp
        z_rho_tran2[start:end,:,:] = z_rho_tran2_tmp
        z_rho_tran3[start:end,:,:] = z_rho_tran3_tmp
        
        # Update the time step
        time_step = time_step + time_lengths[j]
    
    # Return these values
    return(z_rho_tran1, z_rho_tran2, z_rho_tran3)
    
    
    
# ------------------------ End Function Definitions --------------------------

# Get all the file names 
# ROMS 2020 output
# dbsed0003 (standard)
file_names = glob('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0003/ocean_his_beaufort_shelf_2020_dbsed0003_*.nc')
# Sort them to be in order
file_names2 = sorted(file_names)

# dbsed0005 (no waves)
file_names_nowaves = glob('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0005_nowaves/ocean_his_beaufort_shelf_2020_dbsed0005_11234017_*.nc')
# Sort them to be in order
file_names3 = sorted(file_names_nowaves)

# dbsed0006 (double waves)
file_names_double_waves = glob('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0006_double_waves/ocean_his_beaufort_shelf_2020_dbsed0006_*.nc')
# Sort them to be in order
file_names4 = sorted(file_names_double_waves)

# Set the transects
xi1 = 152
xi2 = 304
xi3 = 456

# Call the functions 
# --- Standard Run ---
# All sed 
time_int_ssflux_tran1_std, time_int_ssflux_tran2_std, time_int_ssflux_tran3_std, time_int_depth_int_ssflux_tran1_std, time_int_depth_int_ssflux_tran2_std, time_int_depth_int_ssflux_tran3_std = get_time_int_depth_int_ss_flux_allsed_3transects_cross_shore(file_names2, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Mud01
time_int_ssflux_mud01_tran1_std, time_int_ssflux_mud01_tran2_std, time_int_ssflux_mud01_tran3_std, time_int_depth_int_ssflux_mud01_tran1_std, time_int_depth_int_ssflux_mud01_tran2_std, time_int_depth_int_ssflux_mud01_tran3_std = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names2, 'mud_01', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Mud02
time_int_ssflux_mud02_tran1_std, time_int_ssflux_mud02_tran2_std, time_int_ssflux_mud02_tran3_std, time_int_depth_int_ssflux_mud02_tran1_std, time_int_depth_int_ssflux_mud02_tran2_std, time_int_depth_int_ssflux_mud02_tran3_std = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names2, 'mud_02', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand01
time_int_ssflux_sand01_tran1_std, time_int_ssflux_sand01_tran2_std, time_int_ssflux_sand01_tran3_std, time_int_depth_int_ssflux_sand01_tran1_std, time_int_depth_int_ssflux_sand01_tran2_std, time_int_depth_int_ssflux_sand01_tran3_std = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names2, 'sand_01', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand02
time_int_ssflux_sand02_tran1_std, time_int_ssflux_sand02_tran2_std, time_int_ssflux_sand02_tran3_std, time_int_depth_int_ssflux_sand02_tran1_std, time_int_depth_int_ssflux_sand02_tran2_std, time_int_depth_int_ssflux_sand02_tran3_std = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names2, 'sand_02', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand03
time_int_ssflux_sand03_tran1_std, time_int_ssflux_sand03_tran2_std, time_int_ssflux_sand03_tran3_std, time_int_depth_int_ssflux_sand03_tran1_std, time_int_depth_int_ssflux_sand03_tran2_std, time_int_depth_int_ssflux_sand03_tran3_std = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names2, 'sand_03', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Get the depths
# Standard Run (use this for all since it is similar)
z_rho_tran1_std, z_rho_tran2_std, z_rho_tran3_std = get_z_rho_at_transects(file_names2, eta_rho_len, s_rho_len, xi1, xi2, xi3)

# Take the time-averaged z_rho to get average depth for all times
time_avg_z_rho_tran1_std = np.mean(z_rho_tran1_std, axis=0)
time_avg_z_rho_tran2_std = np.mean(z_rho_tran2_std, axis=0)
time_avg_z_rho_tran3_std = np.mean(z_rho_tran3_std, axis=0)

#time_avg_dz_tran1_std_cumsum = np.cumsum(time_avg_dz_tran1_std, axis=0)
#time_avg_dz_tran2_std_cumsum = np.cumsum(time_avg_dz_tran2_std, axis=0)
#time_avg_dz_tran3_std_cumsum = np.cumsum(time_avg_dz_tran3_std, axis=0)


# --- No Waves Run ---
# Call the functions for the no waves run 
# All sed 
time_int_ssflux_tran1_nowaves, time_int_ssflux_tran2_nowaves, time_int_ssflux_tran3_nowaves, time_int_depth_int_ssflux_tran1_nowaves, time_int_depth_int_ssflux_tran2_nowaves, time_int_depth_int_ssflux_tran3_nowaves = get_time_int_depth_int_ss_flux_allsed_3transects_cross_shore(file_names3, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Mud01
time_int_ssflux_mud01_tran1_nowaves, time_int_ssflux_mud01_tran2_nowaves, time_int_ssflux_mud01_tran3_nowaves, time_int_depth_int_ssflux_mud01_tran1_nowaves, time_int_depth_int_ssflux_mud01_tran2_nowaves, time_int_depth_int_ssflux_mud01_tran3_nowaves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names3, 'mud_01', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Mud02
time_int_ssflux_mud02_tran1_nowaves, time_int_ssflux_mud02_tran2_nowaves, time_int_ssflux_mud02_tran3_nowaves, time_int_depth_int_ssflux_mud02_tran1_nowaves, time_int_depth_int_ssflux_mud02_tran2_nowaves, time_int_depth_int_ssflux_mud02_tran3_nowaves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names3, 'mud_02', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand01
time_int_ssflux_sand01_tran1_nowaves, time_int_ssflux_sand01_tran2_nowaves, time_int_ssflux_sand01_tran3_nowaves, time_int_depth_int_ssflux_sand01_tran1_nowaves, time_int_depth_int_ssflux_sand01_tran2_nowaves, time_int_depth_int_ssflux_sand01_tran3_nowaves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names3, 'sand_01', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand02
time_int_ssflux_sand02_tran1_nowaves, time_int_ssflux_sand02_tran2_nowaves, time_int_ssflux_sand02_tran3_nowaves, time_int_depth_int_ssflux_sand02_tran1_nowaves, time_int_depth_int_ssflux_sand02_tran2_nowaves, time_int_depth_int_ssflux_sand02_tran3_nowaves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names3, 'sand_02', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand03
time_int_ssflux_sand03_tran1_nowaves, time_int_ssflux_sand03_tran2_nowaves, time_int_ssflux_sand03_tran3_nowaves, time_int_depth_int_ssflux_sand03_tran1_nowaves, time_int_depth_int_ssflux_sand03_tran2_nowaves, time_int_depth_int_ssflux_sand03_tran3_nowaves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names3, 'sand_03', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Get the depths
# Standard Run (use this for all since it is similar)
z_rho_tran1_nowaves, z_rho_tran2_nowaves, z_rho_tran3_nowaves = get_z_rho_at_transects(file_names3, eta_rho_len, s_rho_len, xi1, xi2, xi3)

# Take the time-averaged z_rho to get average depth for all times
time_avg_z_rho_tran1_nowaves = np.mean(z_rho_tran1_nowaves, axis=0)
time_avg_z_rho_tran2_nowaves = np.mean(z_rho_tran2_nowaves, axis=0)
time_avg_z_rho_tran3_nowaves = np.mean(z_rho_tran3_nowaves, axis=0)


# --- Double Waves Run ---
# Call the functions for the double waves run (copy from above)
# All sed 
time_int_ssflux_tran1_double_waves, time_int_ssflux_tran2_double_waves, time_int_ssflux_tran3_double_waves, time_int_depth_int_ssflux_tran1_double_waves, time_int_depth_int_ssflux_tran2_double_waves, time_int_depth_int_ssflux_tran3_double_waves = get_time_int_depth_int_ss_flux_allsed_3transects_cross_shore(file_names4, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Mud01
time_int_ssflux_mud01_tran1_double_waves, time_int_ssflux_mud01_tran2_double_waves, time_int_ssflux_mud01_tran3_double_waves, time_int_depth_int_ssflux_mud01_tran1_double_waves, time_int_depth_int_ssflux_mud01_tran2_double_waves, time_int_depth_int_ssflux_mud01_tran3_double_waves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names4, 'mud_01', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Mud02
time_int_ssflux_mud02_tran1_double_waves, time_int_ssflux_mud02_tran2_double_waves, time_int_ssflux_mud02_tran3_double_waves, time_int_depth_int_ssflux_mud02_tran1_double_waves, time_int_depth_int_ssflux_mud02_tran2_double_waves, time_int_depth_int_ssflux_mud02_tran3_double_waves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names4, 'mud_02', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand01
time_int_ssflux_sand01_tran1_double_waves, time_int_ssflux_sand01_tran2_double_waves, time_int_ssflux_sand01_tran3_double_waves, time_int_depth_int_ssflux_sand01_tran1_double_waves, time_int_depth_int_ssflux_sand01_tran2_double_waves, time_int_depth_int_ssflux_sand01_tran3_double_waves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names4, 'sand_01', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand02
time_int_ssflux_sand02_tran1_double_waves, time_int_ssflux_sand02_tran2_double_waves, time_int_ssflux_sand02_tran3_double_waves, time_int_depth_int_ssflux_sand02_tran1_double_waves, time_int_depth_int_ssflux_sand02_tran2_double_waves, time_int_depth_int_ssflux_sand02_tran3_double_waves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names4, 'sand_02', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Sand03
time_int_ssflux_sand03_tran1_double_waves, time_int_ssflux_sand03_tran2_double_waves, time_int_ssflux_sand03_tran3_double_waves, time_int_depth_int_ssflux_sand03_tran1_double_waves, time_int_depth_int_ssflux_sand03_tran2_double_waves, time_int_depth_int_ssflux_sand03_tran3_double_waves = get_time_int_depth_int_ss_flux_1sed_3transects_cross_shore(file_names4, 'sand_03', eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# Get the depths
# Standard Run (use this for all since it is similar)
z_rho_tran1_double_waves, z_rho_tran2_double_waves, z_rho_tran3_double_waves = get_z_rho_at_transects(file_names4, eta_rho_len, s_rho_len, xi1, xi2, xi3)

# Take the time-averaged z_rho to get average depth for all times
time_avg_z_rho_tran1_double_waves = np.mean(z_rho_tran1_double_waves, axis=0)
time_avg_z_rho_tran2_double_waves = np.mean(z_rho_tran2_double_waves, axis=0)
time_avg_z_rho_tran3_double_waves = np.mean(z_rho_tran3_double_waves, axis=0)


# Get the time-averaged u currents for each run
# --- Standard Run ---
time_avg_u_tran1_std, time_avg_u_tran2_std, time_avg_u_tran3_std = get_time_avg_u_currents_at_transects(file_names2, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# --- No Waves Run ---
time_avg_u_tran1_nowaves, time_avg_u_tran2_nowaves, time_avg_u_tran3_nowaves = get_time_avg_u_currents_at_transects(file_names3, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)

# --- Double Waves Run ---
time_avg_u_tran1_double_waves, time_avg_u_tran2_double_waves, time_avg_u_tran3_double_waves = get_time_avg_u_currents_at_transects(file_names4, eta_rho_len, s_rho_len, xi1, xi2, xi3, regridder_u2rho)


# More data prep for plotting
# Set the number of cells in the sponge on each open boundary
c_west = 36
c_north = 45
c_east = 36

# Make a fake xy with the right resolution to be able to plot without the angle
x_rho_flat = np.arange(0,750*len(grid.x_rho[0,:]),750)
y_rho_flat = np.arange(0,600*len(grid.y_rho[:,0]),600)

# Prep the data by multiplying by the mask and trimming
# Trim 
x_rho_flat_trimmed = x_rho_flat[c_west:-c_west]

# Multiply bathymetry by mask and trim
h_masked = grid.h.values*grid.mask_rho.values*mask_rho_nan.nudge_mask_rho_nan
h_masked_trimmed = h_masked[:,c_west:-c_west]

# Make it so land will appear
temp_mask = grid.mask_rho.copy()
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)

# Make reverse mask so land can be gray 
reverse_colors = colors.ListedColormap(['darkgray','#FF000000'])
bounds=[0,1,2]
norm = colors.BoundaryNorm(bounds, reverse_colors.N)

# Make new mask so land is gray - NEW
mask_rho2 = grid.mask_rho.values
mask_rho2[mask_rho2==1.0] = 2.0
mask_rho2_trimmed = mask_rho2[:,c_west:-c_west]


# ----------------------------- Plot the Data ---------------------------------

# ---------------------------------------------------------------------------
# ---------------------- Plot 1: Transect Locations ----------------------
# ---------------------------------------------------------------------------
# Make a plot of the transect locations 

# Set the colormap
cmap_bathy = cmocean.cm.deep

# Set the levels
lev_bathy = np.arange(10, 100, 10)

# Set the bathymetry 

# Make the figure
fig1, ax1 = plt.subplots(figsize=(16,8))

# Plot bathymetry contours 
cs1 = ax1.contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       grid.h[:,c_west:-c_west].values, lev_bathy, cmap=cmap_bathy, extend='both')
# Plot land as gray 
cs0a = ax1.contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Make and label the colorbar
cbar1_ax = fig1.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar1 = plt.colorbar(cs1, orientation='vertical', cax=cbar1_ax).set_label(label='Bathymetry (m)', size=fontsize-2, labelpad=15)

# Plot the transects 
ax1.axvline(x=x_rho_flat_trimmed[xi1]/1000, linestyle='-', linewidth=8, color='#D81B60', label='Transect 1')
ax1.axvline(x=x_rho_flat_trimmed[xi2]/1000, linestyle='-', linewidth=8, color='#1E88E5', label='Transect 2')
ax1.axvline(x=x_rho_flat_trimmed[xi3]/1000, linestyle='-', linewidth=8, color='#FFC107', label='Transect 3')
#ax1.plot(x_rho_flat_trimmed[xi1]/1000, y_rho_flat/1000, '-', color='red', label='Transect 1')
#ax1.plot(x_rho_flat_trimmed[xi1]/1000, y_rho_flat/1000, '-', color='green', label='Transect 2')
#ax1.plot(x_rho_flat_trimmed[xi1]/1000, y_rho_flat/1000, '-', color='blue', label='Transect 3')

# Add a legend 
ax1.legend(fontsize=fontsize-2)

# Label the plot
#plt.setp(ax1.get_xticklabels(), visible=False)
ax1.set_ylabel('Y (km)', fontsize=fontsize-2)
ax1.set_xlabel('X (km)', fontsize=fontsize-2)
ax1.set_title('Transect Locations', fontsize=fontsize)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/transect_locations_0002.png', bbox_inches='tight')



# ***** Standard Run *****
# ---------------------------------------------------------------------------
# ---------- Plot 2: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the standard run 

# Set the colormap
cmap_ssflux_std = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_std = np.arange(-100, 100, 1) 
vmin2 = -6000
vmax2 = 6000


# Make the figure
fig2, ax2 = plt.subplots(3, figsize=(22,21))

# Set the title
fig2.suptitle('Time-Integrated Suspended Sediment Flux \n(Standard Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# Transect 1
# Plot the ss flux
cs2 = ax2[0].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran1_std, time_int_ssflux_tran1_std, 
               vmin=vmin2, vmax=vmax2, cmap=cmap_ssflux_std)
# Set the y lim
ax2[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[0].get_xticklabels(), visible=False)
ax2[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax2[0].set_title('Transect 1', fontsize=fontsize-1)
# Make background gray
ax2[0].set_facecolor('lightgray')


# Transect 2
# Plot the ss flux 
cs3 = ax2[1].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran2_std, time_int_ssflux_tran2_std, 
               vmin=vmin2, vmax=vmax2, cmap=cmap_ssflux_std)
# Set the y lim
ax2[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[1].get_xticklabels(), visible=False)
ax2[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax2[1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax2[1].set_facecolor('lightgray')

# Transect 3
# Plot the ss flux 
cs4 = ax2[2].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran3_std, time_int_ssflux_tran3_std, 
               vmin=vmin2, vmax=vmax2, cmap=cmap_ssflux_std)
# Set the y lim
ax2[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax2[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax2[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax2[2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax2[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar2_ax = fig2.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar2 = plt.colorbar(cs4, orientation='vertical', cax=cbar2_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=15)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_standard_run_3transects_0009.png', bbox_inches='tight')

# Print mins and maxes 
print('time_int_ssflux_tran1_std min: ', np.nanmin(time_int_ssflux_tran1_std))
print('time_int_ssflux_tran1_std max: ', np.nanmax(time_int_ssflux_tran1_std))
print('time_int_ssflux_tran1_std min magnitude: ', np.nanmin(abs(time_int_ssflux_tran1_std)))
print('time_int_ssflux_tran1_std max magnitude: ', np.nanmax(abs(time_int_ssflux_tran1_std)))
print('time_int_ssflux_tran2_std min: ', np.nanmin(time_int_ssflux_tran2_std))
print('time_int_ssflux_tran2_std max: ', np.nanmax(time_int_ssflux_tran2_std))
print('time_int_ssflux_tran2_std min magnitude: ', np.nanmin(abs(time_int_ssflux_tran2_std)))
print('time_int_ssflux_tran2_std max magnitude: ', np.nanmax(abs(time_int_ssflux_tran2_std)))
print('time_int_ssflux_tran3_std min: ', np.nanmin(time_int_ssflux_tran3_std))
print('time_int_ssflux_tran3_std max: ', np.nanmax(time_int_ssflux_tran3_std))
print('time_int_ssflux_tran3_std min magnitude: ', np.nanmin(abs(time_int_ssflux_tran3_std)))
print('time_int_ssflux_tran3_std max magnitude: ', np.nanmax(abs(time_int_ssflux_tran3_std)))



# ---------------------------------------------------------------------------
# ---------- Plot 3: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the standard run but with 
# 6 subplots, looking at different depths/regions 

# Set the colormap
cmap_ssflux_std = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_std = np.arange(-100, 100, 1) 
vmin3 = -5000 # -10000
vmax3 = 5000 # 10000


# Make the figure
fig3, ax3 = plt.subplots(3, 2, figsize=(22,16))

# Set the title
fig3.suptitle('Time-Integrated Suspended Sediment Flux \n(Standard Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# Transect 1 - shallow
# Plot the ss flux
cs5 = ax3[0,0].pcolormesh(y_rho_flat[41:110]/1000, time_avg_z_rho_tran1_std[:,41:110], time_int_ssflux_tran1_std[:,41:110], # [41:110]
               vmin=vmin3, vmax=vmax3, cmap=cmap_ssflux_std)
# Set the y lim
#ax3[0,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[0,0].get_xticklabels(), visible=False)
ax3[0,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax3[0,0].set_title('Shallow', fontsize=fontsize-1)
# Make background gray
ax3[0,0].set_facecolor('lightgray')

# Transect 1 - deep
# Plot the ss flux
cs6 = ax3[0,1].pcolormesh(y_rho_flat[110:165]/1000, time_avg_z_rho_tran1_std[:,110:165], time_int_ssflux_tran1_std[:,110:165], 
               vmin=vmin3, vmax=vmax3, cmap=cmap_ssflux_std)
# Set the y lim
#ax3[0,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[0,1].get_xticklabels(), visible=False)
ax3[0,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax3[0,1].set_title('Deep', fontsize=fontsize-1)
# Make background gray
ax3[0,1].set_facecolor('lightgray')


# Transect 2 - shallow
# Plot the ss flux 
cs7 = ax3[1,0].pcolormesh(y_rho_flat[17:95]/1000, time_avg_z_rho_tran2_std[:,17:95], time_int_ssflux_tran2_std[:,17:95], # [17:95]
               vmin=vmin3, vmax=vmax3, cmap=cmap_ssflux_std)
# Set the y lim
#ax3[1,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[1,0].get_xticklabels(), visible=False)
ax3[1,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[1,0].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax3[1,0].set_facecolor('lightgray')

# Transect 2 - deep
# Plot the ss flux 
cs8 = ax3[1,1].pcolormesh(y_rho_flat[95:150]/1000, time_avg_z_rho_tran2_std[:,95:150], time_int_ssflux_tran2_std[:,95:150], 
               vmin=vmin3, vmax=vmax3, cmap=cmap_ssflux_std)
# Set the y lim
#ax3[1,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[1,1].get_xticklabels(), visible=False)
ax3[1,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[1,1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax3[1,1].set_facecolor('lightgray')

# Transect 3 - shallow
# Plot the ss flux 
cs9 = ax3[2,0].pcolormesh(y_rho_flat[12:70]/1000, time_avg_z_rho_tran3_std[:,12:70], time_int_ssflux_tran3_std[:,12:70], # [12:70]
               vmin=vmin3, vmax=vmax3, cmap=cmap_ssflux_std)
# Set the y lim
#ax3[2,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[2,0].get_xticklabels(), visible=False)
ax3[2,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax3[2,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[2,0].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax3[2,0].set_facecolor('lightgray')

# Transect 3 - deep
# Plot the ss flux 
cs10 = ax3[2,1].pcolormesh(y_rho_flat[70:125]/1000, time_avg_z_rho_tran3_std[:,70:125], time_int_ssflux_tran3_std[:,70:125], 
               vmin=vmin3, vmax=vmax3, cmap=cmap_ssflux_std)
# Set the y lim
#ax3[2,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[2,1].get_xticklabels(), visible=False)
ax3[2,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax3[2,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[2,1].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax3[2,1].set_facecolor('lightgray')

# Make and label the colorbar
cbar3_ax = fig3.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar3 = plt.colorbar(cs10, orientation='vertical', cax=cbar3_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.30, wspace=0.22) #0.08

# Add labels for the rows
plt.text(-0.05, 0.769, 'Transect 1', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(-0.05, 0.469, 'Transect 2', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(-0.05, 0.200, 'Transect 3', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_standard_run_3transects_diff_depths_0006.png', bbox_inches='tight')


# -------------------------------------------------------------------------------
# ---- Plot 4: 2D Plot of All Sed Time-Integrated, Depth-Integrated SS Flux ------
# --------------------------------------------------------------------------------
# Make a plot of the time-integrated, depth-integrated suspended sediment flux
# for the standard run with one plot for all sedimnet and the other plots for 
# each sediment class or start with them all on the same y-axis for now 

# Make the figure
fig4, ax4 = plt.subplots(2, figsize=(16,12)) # (22,21)

# Set the title
fig4.suptitle('Time-Integrated, Depth-Integrated Suspended Sediment Flux \n(Standard Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# First subplot = one line for each transect, all sediment 
ax4[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran1_std, color='#D81B60', label='Transect 1')
ax4[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran2_std, color='#1E88E5', label='Transect 2')
ax4[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran3_std, color='#FFC107', label='Transect 3')
# Label the plot
plt.setp(ax4[0].get_xticklabels(), visible=False)
ax4[0].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)
# Add a legend
ax4[0].legend(fontsize=fontsize-2)



# Second subplot = 5 lines for each transect (one for each class)
# Transect 1
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran1_std, linestyle='solid', color='#D81B60', label='Transect 1, Mud01')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran1_std, linestyle='dashed', color='#D81B60', label='Transect 1, Mud02')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran1_std, linestyle='dotted', color='#D81B60', label='Transect 1, Sand01')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran1_std, linestyle='dashdot', color='#D81B60', label='Transect 1, Sand02')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran1_std, linestyle='', marker='+', color='#D81B60', label='Transect 1, Sand03')

# Transect 2
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran2_std, linestyle='solid', color='#1E88E5', label='Transect 2, Mud01')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran2_std, linestyle='dashed', color='#1E88E5', label='Transect 2, Mud02')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran2_std, linestyle='dotted', color='#1E88E5', label='Transect 2, Sand01')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran2_std, linestyle='dashdot', color='#1E88E5', label='Transect 2, Sand02')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran2_std, linestyle='', marker='+', color='#1E88E5', label='Transect 2, Sand03')

# Transect 3
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran3_std, linestyle='solid', color='#FFC107', label='Transect 3, Mud01')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran3_std, linestyle='dashed', color='#FFC107', label='Transect 3, Mud02')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran3_std, linestyle='dotted', color='#FFC107', label='Transect 3, Sand01')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran3_std, linestyle='dashdot', color='#FFC107', label='Transect 3, Sand02')
ax4[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran3_std, linestyle='', marker='+', color='#FFC107', label='Transect 3, Sand03')
# Label the plot
#plt.setp(ax4[1].get_xticklabels(), visible=False)
ax4[1].set_xlabel('Y (km)', fontsize=fontsize-2)
ax4[1].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)

# Add a legend
ax4[1].legend(bbox_to_anchor=(1.0, -0.20), ncol=3, fontsize=fontsize-2) #(0.2, -0.45)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.08) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_standard_run_2d_3transects_byclass_005.png', bbox_inches='tight')



# ---------------------------------------------------------------------------
# ---------- Plot 5: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the standard run but 
# only for top 60 m depth at each transect 

# Set the colormap
cmap_ssflux_std = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_std = np.arange(-100, 100, 1) 
vmin5 = -5000
vmax5 = 5000


# Make the figure
fig5, ax5 = plt.subplots(3, figsize=(22,21))

# Set the title
fig5.suptitle('Time-Integrated Suspended Sediment Flux \n(Standard Run - All Sed)', x=0.5, y=0.93, fontsize=fontsize)

# Transect 1
# Plot the ss flux
cs11 = ax5[0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_std[:,41:165], time_int_ssflux_tran1_std[:,41:165], 
               vmin=vmin5, vmax=vmax5, cmap=cmap_ssflux_std)
# Set the y lim
ax5[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax5[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax5[0].set_title('Transect 1', fontsize=fontsize-1)
# Make background gray
ax5[0].set_facecolor('lightgray')


# Transect 2
# Plot the ss flux 
cs12 = ax5[1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_std[:,17:150], time_int_ssflux_tran2_std[:,17:150], 
               vmin=vmin5, vmax=vmax5, cmap=cmap_ssflux_std)
# Set the y lim
ax5[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
ax5[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax5[1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax5[1].set_facecolor('lightgray')

# Transect 3
# Plot the ss flux 
cs13 = ax5[2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_std[:,12:125], time_int_ssflux_tran3_std[:,12:125], 
               vmin=vmin5, vmax=vmax5, cmap=cmap_ssflux_std)
# Set the y lim
ax5[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax5[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax5[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax5[2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax5[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar5_ax = fig5.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar5 = plt.colorbar(cs13, orientation='vertical', cax=cbar5_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.25) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_standard_run_3transects_60m_0004.png', bbox_inches='tight')


# ***** No Waves Run *****
# ---------------------------------------------------------------------------
# ---------- Plot 6: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the no waves run

# Set the colormap
cmap_ssflux_nowaves = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_std = np.arange(-100, 100, 1) 
vmin6 = -6000
vmax6 = 6000


# Make the figure
fig6, ax6 = plt.subplots(3, figsize=(22,21))

# Set the title
fig6.suptitle('Time-Integrated Suspended Sediment Flux \n(No Waves Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# Transect 1
# Plot the ss flux
cs14 = ax6[0].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran1_nowaves, time_int_ssflux_tran1_nowaves, 
               vmin=vmin6, vmax=vmax6, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax6[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[0].get_xticklabels(), visible=False)
ax6[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax6[0].set_title('Transect 1', fontsize=fontsize-1)
# Make background gray
ax6[0].set_facecolor('lightgray')


# Transect 2
# Plot the ss flux 
cs15 = ax6[1].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran2_nowaves, time_int_ssflux_tran2_nowaves, 
               vmin=vmin6, vmax=vmax6, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax6[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[1].get_xticklabels(), visible=False)
ax6[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax6[1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax6[1].set_facecolor('lightgray')

# Transect 3
# Plot the ss flux 
cs16 = ax6[2].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran3_nowaves, time_int_ssflux_tran3_nowaves, 
               vmin=vmin6, vmax=vmax6, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax6[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax6[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax6[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax6[2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax6[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar6_ax = fig6.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar6 = plt.colorbar(cs16, orientation='vertical', cax=cbar6_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=15)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_nowaves_run_3transects_0002.png', bbox_inches='tight')

# Print mins and maxes 
print('time_int_ssflux_tran1_nowaves min: ', np.nanmin(time_int_ssflux_tran1_nowaves))
print('time_int_ssflux_tran1_nowaves max: ', np.nanmax(time_int_ssflux_tran1_nowaves))
print('time_int_ssflux_tran1_nowaves min magnitude: ', np.nanmin(abs(time_int_ssflux_tran1_nowaves)))
print('time_int_ssflux_tran1_nowaves max magnitude: ', np.nanmax(abs(time_int_ssflux_tran1_nowaves)))
print('time_int_ssflux_tran2_nowaves min: ', np.nanmin(time_int_ssflux_tran2_nowaves))
print('time_int_ssflux_tran2_nowaves max: ', np.nanmax(time_int_ssflux_tran2_nowaves))
print('time_int_ssflux_tran2_nowaves min magnitude: ', np.nanmin(abs(time_int_ssflux_tran2_nowaves)))
print('time_int_ssflux_tran2_nowaves max magnitude: ', np.nanmax(abs(time_int_ssflux_tran2_nowaves)))
print('time_int_ssflux_tran3_nowaves min: ', np.nanmin(time_int_ssflux_tran3_nowaves))
print('time_int_ssflux_tran3_nowaves max: ', np.nanmax(time_int_ssflux_tran3_nowaves))
print('time_int_ssflux_tran3_nowaves min magnitude: ', np.nanmin(abs(time_int_ssflux_tran3_nowaves)))
print('time_int_ssflux_tran3_nowaves max magnitude: ', np.nanmax(abs(time_int_ssflux_tran3_nowaves)))



# ---------------------------------------------------------------------------
# ---------- Plot 7: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the no waves run but with 
# 6 subplots, looking at different depths/regions 

# Set the colormap
cmap_ssflux_nowaves = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_nowaves = np.arange(-100, 100, 1) 
vmin7 = -5000 # -10000
vmax7 = 5000 # 10000


# Make the figure
fig7, ax7 = plt.subplots(3, 2, figsize=(22,16))

# Set the title
fig7.suptitle('Time-Integrated Suspended Sediment Flux \n(No Waves Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# Transect 1 - shallow
# Plot the ss flux
cs17 = ax7[0,0].pcolormesh(y_rho_flat[41:110]/1000, time_avg_z_rho_tran1_nowaves[:,41:110], time_int_ssflux_tran1_nowaves[:,41:110], # [41:110]
               vmin=vmin7, vmax=vmax7, cmap=cmap_ssflux_nowaves)
# Set the y lim
#ax3[0,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[0,0].get_xticklabels(), visible=False)
ax7[0,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax7[0,0].set_title('Shallow', fontsize=fontsize-1)
# Make background gray
ax7[0,0].set_facecolor('lightgray')

# Transect 1 - deep
# Plot the ss flux
cs18 = ax7[0,1].pcolormesh(y_rho_flat[110:165]/1000, time_avg_z_rho_tran1_nowaves[:,110:165], time_int_ssflux_tran1_nowaves[:,110:165], 
               vmin=vmin7, vmax=vmax7, cmap=cmap_ssflux_nowaves)
# Set the y lim
#ax3[0,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[0,1].get_xticklabels(), visible=False)
ax7[0,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax7[0,1].set_title('Deep', fontsize=fontsize-1)
# Make background gray
ax7[0,1].set_facecolor('lightgray')


# Transect 2 - shallow
# Plot the ss flux 
cs19 = ax7[1,0].pcolormesh(y_rho_flat[17:95]/1000, time_avg_z_rho_tran2_nowaves[:,17:95], time_int_ssflux_tran2_nowaves[:,17:95], # [17:95]
               vmin=vmin7, vmax=vmax7, cmap=cmap_ssflux_nowaves)
# Set the y lim
#ax3[1,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[1,0].get_xticklabels(), visible=False)
ax7[1,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[1,0].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax7[1,0].set_facecolor('lightgray')

# Transect 2 - deep
# Plot the ss flux 
cs20 = ax7[1,1].pcolormesh(y_rho_flat[95:150]/1000, time_avg_z_rho_tran2_nowaves[:,95:150], time_int_ssflux_tran2_nowaves[:,95:150], 
               vmin=vmin7, vmax=vmax7, cmap=cmap_ssflux_nowaves)
# Set the y lim
#ax3[1,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[1,1].get_xticklabels(), visible=False)
ax7[1,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[1,1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax7[1,1].set_facecolor('lightgray')

# Transect 3 - shallow
# Plot the ss flux 
cs21 = ax7[2,0].pcolormesh(y_rho_flat[12:70]/1000, time_avg_z_rho_tran3_nowaves[:,12:70], time_int_ssflux_tran3_nowaves[:,12:70], # [12:70]
               vmin=vmin7, vmax=vmax7, cmap=cmap_ssflux_nowaves)
# Set the y lim
#ax3[2,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[2,0].get_xticklabels(), visible=False)
ax7[2,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax7[2,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[2,0].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax7[2,0].set_facecolor('lightgray')

# Transect 3 - deep
# Plot the ss flux 
cs22 = ax7[2,1].pcolormesh(y_rho_flat[70:125]/1000, time_avg_z_rho_tran3_nowaves[:,70:125], time_int_ssflux_tran3_nowaves[:,70:125], 
               vmin=vmin7, vmax=vmax7, cmap=cmap_ssflux_nowaves)
# Set the y lim
#ax3[2,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[2,1].get_xticklabels(), visible=False)
ax7[2,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax7[2,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[2,1].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax7[2,1].set_facecolor('lightgray')

# Make and label the colorbar
cbar7_ax = fig7.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar7 = plt.colorbar(cs22, orientation='vertical', cax=cbar7_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.30, wspace=0.22) #0.08

# Add labels for the rows
plt.text(-0.05, 0.769, 'Transect 1', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(-0.05, 0.469, 'Transect 2', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(-0.05, 0.200, 'Transect 3', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_nowaves_run_3transects_diff_depths_0002.png', bbox_inches='tight')


# -------------------------------------------------------------------------------
# ---- Plot 8: 2D Plot of All Sed Time-Integrated, Depth-Integrated SS Flux ------
# --------------------------------------------------------------------------------
# Make a plot of the time-integrated, depth-integrated suspended sediment flux
# for the no waves run with one plot for all sedimnet and the other plots for 
# each sediment class or start with them all on the same y-axis for now 

# Make the figure
fig8, ax8 = plt.subplots(2, figsize=(16,12)) # (22,21)

# Set the title
fig8.suptitle('Time-Integrated, Depth-Integrated Suspended Sediment Flux \n(No Waves Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# First subplot = one line for each transect, all sediment 
ax8[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran1_nowaves, color='#D81B60', label='Transect 1')
ax8[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran2_nowaves, color='#1E88E5', label='Transect 2')
ax8[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran3_nowaves, color='#FFC107', label='Transect 3')
# Label the plot
plt.setp(ax8[0].get_xticklabels(), visible=False)
ax8[0].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)
# Add a legend
ax8[0].legend(fontsize=fontsize-2)



# Second subplot = 5 lines for each transect (one for each class)
# Transect 1
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran1_nowaves, linestyle='solid', color='#D81B60', label='Transect 1, Mud01')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran1_nowaves, linestyle='dashed', color='#D81B60', label='Transect 1, Mud02')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran1_nowaves, linestyle='dotted', color='#D81B60', label='Transect 1, Sand01')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran1_nowaves, linestyle='dashdot', color='#D81B60', label='Transect 1, Sand02')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran1_nowaves, linestyle='', marker='+', color='#D81B60', label='Transect 1, Sand03')

# Transect 2
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran2_nowaves, linestyle='solid', color='#1E88E5', label='Transect 2, Mud01')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran2_nowaves, linestyle='dashed', color='#1E88E5', label='Transect 2, Mud02')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran2_nowaves, linestyle='dotted', color='#1E88E5', label='Transect 2, Sand01')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran2_nowaves, linestyle='dashdot', color='#1E88E5', label='Transect 2, Sand02')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran2_nowaves, linestyle='', marker='+', color='#1E88E5', label='Transect 2, Sand03')

# Transect 3
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran3_nowaves, linestyle='solid', color='#FFC107', label='Transect 3, Mud01')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran3_nowaves, linestyle='dashed', color='#FFC107', label='Transect 3, Mud02')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran3_nowaves, linestyle='dotted', color='#FFC107', label='Transect 3, Sand01')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran3_nowaves, linestyle='dashdot', color='#FFC107', label='Transect 3, Sand02')
ax8[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran3_nowaves, linestyle='', marker='+', color='#FFC107', label='Transect 3, Sand03')
# Label the plot
#plt.setp(ax4[1].get_xticklabels(), visible=False)
ax8[1].set_xlabel('Y (km)', fontsize=fontsize-2)
ax8[1].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)

# Add a legend
ax8[1].legend(bbox_to_anchor=(1.0, -0.20), ncol=3, fontsize=fontsize-2) #(0.2, -0.45)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.08) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_nowaves_run_2d_3transects_byclass_001.png', bbox_inches='tight')


# ---------------------------------------------------------------------------
# ---------- Plot 9: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the no waves run but 
# only for top 60 m depth at each transect 

# Set the colormap
cmap_ssflux_nowaves = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_nowaves = np.arange(-100, 100, 1) 
vmin9 = -5000
vmax9 = 5000


# Make the figure
fig9, ax9 = plt.subplots(3, figsize=(22,21))

# Set the title
fig9.suptitle('Time-Integrated Suspended Sediment Flux \n(No Waves Run - All Sed)', x=0.5, y=0.93, fontsize=fontsize)

# Transect 1
# Plot the ss flux
cs23 = ax9[0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_nowaves[:,41:165], time_int_ssflux_tran1_nowaves[:,41:165], 
               vmin=vmin9, vmax=vmax9, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax9[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax9[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax9[0].set_title('Transect 1', fontsize=fontsize-1)
# Make background gray
ax9[0].set_facecolor('lightgray')


# Transect 2
# Plot the ss flux 
cs24 = ax9[1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_nowaves[:,17:150], time_int_ssflux_tran2_nowaves[:,17:150], 
               vmin=vmin9, vmax=vmax9, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax9[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
ax9[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax9[1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax9[1].set_facecolor('lightgray')

# Transect 3
# Plot the ss flux 
cs25 = ax9[2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_nowaves[:,12:125], time_int_ssflux_tran3_nowaves[:,12:125], 
               vmin=vmin9, vmax=vmax9, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax9[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax9[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax9[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax9[2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax9[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar9_ax = fig9.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar9 = plt.colorbar(cs25, orientation='vertical', cax=cbar9_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.25) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_nowaves_run_3transects_60m_0002.png', bbox_inches='tight')



# ***** Double Waves Run *****
# ---------------------------------------------------------------------------
# ---------- Plot 10: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the double waves run

# Set the colormap
cmap_ssflux_double_waves = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_double_waves = np.arange(-100, 100, 1) 
vmin10 = -6000
vmax10 = 6000


# Make the figure
fig10, ax10 = plt.subplots(3, figsize=(22,21))

# Set the title
fig10.suptitle('Time-Integrated Suspended Sediment Flux \n(Double Waves Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# Transect 1
# Plot the ss flux
cs26 = ax10[0].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran1_double_waves, time_int_ssflux_tran1_double_waves, 
               vmin=vmin10, vmax=vmax10, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax10[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[0].get_xticklabels(), visible=False)
ax10[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax10[0].set_title('Transect 1', fontsize=fontsize-1)
# Make background gray
ax10[0].set_facecolor('lightgray')


# Transect 2
# Plot the ss flux 
cs27 = ax10[1].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran2_double_waves, time_int_ssflux_tran2_double_waves, 
               vmin=vmin10, vmax=vmax10, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax10[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[1].get_xticklabels(), visible=False)
ax10[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax10[1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax10[1].set_facecolor('lightgray')

# Transect 3
# Plot the ss flux 
cs28 = ax10[2].pcolormesh(y_rho_flat/1000, time_avg_z_rho_tran3_double_waves, time_int_ssflux_tran3_double_waves, 
               vmin=vmin10, vmax=vmax10, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax10[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax10[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax10[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax10[2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax10[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar10_ax = fig10.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar10 = plt.colorbar(cs28, orientation='vertical', cax=cbar10_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=15)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_double_waves_run_3transects_0002.png', bbox_inches='tight')

# Print mins and maxes 
print('time_int_ssflux_tran1_double_waves min: ', np.nanmin(time_int_ssflux_tran1_double_waves))
print('time_int_ssflux_tran1_double_waves max: ', np.nanmax(time_int_ssflux_tran1_double_waves))
print('time_int_ssflux_tran1_double_waves min magnitude: ', np.nanmin(abs(time_int_ssflux_tran1_double_waves)))
print('time_int_ssflux_tran1_double_waves max magnitude: ', np.nanmax(abs(time_int_ssflux_tran1_double_waves)))
print('time_int_ssflux_tran2_double_waves min: ', np.nanmin(time_int_ssflux_tran2_double_waves))
print('time_int_ssflux_tran2_double_waves max: ', np.nanmax(time_int_ssflux_tran2_double_waves))
print('time_int_ssflux_tran2_double_waves min magnitude: ', np.nanmin(abs(time_int_ssflux_tran2_double_waves)))
print('time_int_ssflux_tran2_double_waves max magnitude: ', np.nanmax(abs(time_int_ssflux_tran2_double_waves)))
print('time_int_ssflux_tran3_double_waves min: ', np.nanmin(time_int_ssflux_tran3_double_waves))
print('time_int_ssflux_tran3_double_waves max: ', np.nanmax(time_int_ssflux_tran3_double_waves))
print('time_int_ssflux_tran3_double_waves min magnitude: ', np.nanmin(abs(time_int_ssflux_tran3_double_waves)))
print('time_int_ssflux_tran3_double_waves max magnitude: ', np.nanmax(abs(time_int_ssflux_tran3_double_waves)))



# ---------------------------------------------------------------------------
# ---------- Plot 11: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the double waves run but with 
# 6 subplots, looking at different depths/regions 

# Set the colormap
cmap_ssflux_double_waves = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_double_waves = np.arange(-100, 100, 1) 
vmin11 = -5000 # -10000
vmax11 = 5000 # 10000


# Make the figure
fig11, ax11 = plt.subplots(3, 2, figsize=(22,16))

# Set the title
fig11.suptitle('Time-Integrated Suspended Sediment Flux \n(Double Waves Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# Transect 1 - shallow
# Plot the ss flux
cs29 = ax11[0,0].pcolormesh(y_rho_flat[41:110]/1000, time_avg_z_rho_tran1_double_waves[:,41:110], time_int_ssflux_tran1_double_waves[:,41:110], # [41:110]
               vmin=vmin11, vmax=vmax11, cmap=cmap_ssflux_double_waves)
# Set the y lim
#ax3[0,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[0,0].get_xticklabels(), visible=False)
ax11[0,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax11[0,0].set_title('Shallow', fontsize=fontsize-1)
# Make background gray
ax11[0,0].set_facecolor('lightgray')

# Transect 1 - deep
# Plot the ss flux
cs30 = ax11[0,1].pcolormesh(y_rho_flat[110:165]/1000, time_avg_z_rho_tran1_double_waves[:,110:165], time_int_ssflux_tran1_double_waves[:,110:165], 
               vmin=vmin11, vmax=vmax11, cmap=cmap_ssflux_double_waves)
# Set the y lim
#ax3[0,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[0,1].get_xticklabels(), visible=False)
ax11[0,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax11[0,1].set_title('Deep', fontsize=fontsize-1)
# Make background gray
ax11[0,1].set_facecolor('lightgray')


# Transect 2 - shallow
# Plot the ss flux 
cs31 = ax11[1,0].pcolormesh(y_rho_flat[17:95]/1000, time_avg_z_rho_tran2_double_waves[:,17:95], time_int_ssflux_tran2_double_waves[:,17:95], # [17:95]
               vmin=vmin11, vmax=vmax11, cmap=cmap_ssflux_double_waves)
# Set the y lim
#ax3[1,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[1,0].get_xticklabels(), visible=False)
ax11[1,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[1,0].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax11[1,0].set_facecolor('lightgray')

# Transect 2 - deep
# Plot the ss flux 
cs32 = ax11[1,1].pcolormesh(y_rho_flat[95:150]/1000, time_avg_z_rho_tran2_double_waves[:,95:150], time_int_ssflux_tran2_double_waves[:,95:150], 
               vmin=vmin11, vmax=vmax11, cmap=cmap_ssflux_double_waves)
# Set the y lim
#ax3[1,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[1,1].get_xticklabels(), visible=False)
ax11[1,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[1,1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax11[1,1].set_facecolor('lightgray')

# Transect 3 - shallow
# Plot the ss flux 
cs33 = ax11[2,0].pcolormesh(y_rho_flat[12:70]/1000, time_avg_z_rho_tran3_double_waves[:,12:70], time_int_ssflux_tran3_double_waves[:,12:70], # [12:70]
               vmin=vmin11, vmax=vmax11, cmap=cmap_ssflux_double_waves)
# Set the y lim
#ax3[2,0].set_ylim([-30,0])
# Label the plot
#plt.setp(ax3[2,0].get_xticklabels(), visible=False)
ax11[2,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax11[2,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[2,0].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax11[2,0].set_facecolor('lightgray')

# Transect 3 - deep
# Plot the ss flux 
cs34 = ax11[2,1].pcolormesh(y_rho_flat[70:125]/1000, time_avg_z_rho_tran3_double_waves[:,70:125], time_int_ssflux_tran3_double_waves[:,70:125], 
               vmin=vmin11, vmax=vmax11, cmap=cmap_ssflux_double_waves)
# Set the y lim
#ax3[2,1].set_ylim([-60,-30])
# Label the plot
#plt.setp(ax3[2,1].get_xticklabels(), visible=False)
ax11[2,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax11[2,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax3[2,1].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax11[2,1].set_facecolor('lightgray')

# Make and label the colorbar
cbar11_ax = fig11.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar11 = plt.colorbar(cs34, orientation='vertical', cax=cbar11_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.30, wspace=0.22) #0.08

# Add labels for the rows
plt.text(-0.05, 0.769, 'Transect 1', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(-0.05, 0.469, 'Transect 2', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(-0.05, 0.200, 'Transect 3', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_double_waves_run_3transects_diff_depths_0002.png', bbox_inches='tight')


# -------------------------------------------------------------------------------
# ---- Plot 12: 2D Plot of All Sed Time-Integrated, Depth-Integrated SS Flux ------
# --------------------------------------------------------------------------------
# Make a plot of the time-integrated, depth-integrated suspended sediment flux
# for the double waves run with one plot for all sedimnet and the other plots for 
# each sediment class or start with them all on the same y-axis for now 

# Make the figure
fig12, ax12 = plt.subplots(2, figsize=(16,12)) # (22,21)

# Set the title
fig12.suptitle('Time-Integrated, Depth-Integrated Suspended Sediment Flux \n(Double Waves Run - All Sed)', x=0.5, y=0.95, fontsize=fontsize)

# First subplot = one line for each transect, all sediment 
ax12[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran1_double_waves, color='#D81B60', label='Transect 1')
ax12[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran2_double_waves, color='#1E88E5', label='Transect 2')
ax12[0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran3_double_waves, color='#FFC107', label='Transect 3')
# Label the plot
plt.setp(ax12[0].get_xticklabels(), visible=False)
ax12[0].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)
# Add a legend
ax12[0].legend(fontsize=fontsize-2)



# Second subplot = 5 lines for each transect (one for each class)
# Transect 1
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran1_double_waves, linestyle='solid', color='#D81B60', label='Transect 1, Mud01')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran1_double_waves, linestyle='dashed', color='#D81B60', label='Transect 1, Mud02')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran1_double_waves, linestyle='dotted', color='#D81B60', label='Transect 1, Sand01')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran1_double_waves, linestyle='dashdot', color='#D81B60', label='Transect 1, Sand02')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran1_double_waves, linestyle='', marker='+', color='#D81B60', label='Transect 1, Sand03')

# Transect 2
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran2_double_waves, linestyle='solid', color='#1E88E5', label='Transect 2, Mud01')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran2_double_waves, linestyle='dashed', color='#1E88E5', label='Transect 2, Mud02')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran2_double_waves, linestyle='dotted', color='#1E88E5', label='Transect 2, Sand01')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran2_double_waves, linestyle='dashdot', color='#1E88E5', label='Transect 2, Sand02')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran2_double_waves, linestyle='', marker='+', color='#1E88E5', label='Transect 2, Sand03')

# Transect 3
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran3_double_waves, linestyle='solid', color='#FFC107', label='Transect 3, Mud01')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran3_double_waves, linestyle='dashed', color='#FFC107', label='Transect 3, Mud02')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran3_double_waves, linestyle='dotted', color='#FFC107', label='Transect 3, Sand01')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran3_double_waves, linestyle='dashdot', color='#FFC107', label='Transect 3, Sand02')
ax12[1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran3_double_waves, linestyle='', marker='+', color='#FFC107', label='Transect 3, Sand03')
# Label the plot
#plt.setp(ax4[1].get_xticklabels(), visible=False)
ax12[1].set_xlabel('Y (km)', fontsize=fontsize-2)
ax12[1].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)

# Add a legend
ax12[1].legend(bbox_to_anchor=(1.0, -0.20), ncol=3, fontsize=fontsize-2) #(0.2, -0.45)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.08) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_double_waves_run_2d_3transects_byclass_001.png', bbox_inches='tight')


# ---------------------------------------------------------------------------
# ---------- Plot 13: Transect of All Sed Time-Integrated SS Flux ------------
# ---------------------------------------------------------------------------
# Make suplots of the SS flux for each transect in the double waves run but 
# only for top 60 m depth at each transect 

# Set the colormap
cmap_ssflux_double_waves = cmocean.cm.curl

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_double_waves = np.arange(-100, 100, 1) 
vmin13 = -5000
vmax13 = 5000


# Make the figure
fig13, ax13 = plt.subplots(3, figsize=(22,21))

# Set the title
fig13.suptitle('Time-Integrated Suspended Sediment Flux \n(Double Waves Run - All Sed)', x=0.5, y=0.93, fontsize=fontsize)

# Transect 1
# Plot the ss flux
cs35 = ax13[0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_double_waves[:,41:165], time_int_ssflux_tran1_double_waves[:,41:165], 
               vmin=vmin13, vmax=vmax13, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax13[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax13[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax13[0].set_title('Transect 1', fontsize=fontsize-1)
# Make background gray
ax13[0].set_facecolor('lightgray')


# Transect 2
# Plot the ss flux 
cs36 = ax13[1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_double_waves[:,17:150], time_int_ssflux_tran2_double_waves[:,17:150], 
               vmin=vmin13, vmax=vmax13, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax13[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
ax13[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax13[1].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax13[1].set_facecolor('lightgray')

# Transect 3
# Plot the ss flux 
cs37 = ax13[2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_double_waves[:,12:125], time_int_ssflux_tran3_double_waves[:,12:125], 
               vmin=vmin13, vmax=vmax13, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax13[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax13[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax13[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax13[2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax13[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar13_ax = fig13.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar13 = plt.colorbar(cs37, orientation='vertical', cax=cbar13_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.25) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_double_waves_run_3transects_60m_0002.png', bbox_inches='tight')




# ******** Combo Plots ************
# ---------------------------------------------------------------------------
# ---------- Plot 14: Transect 1 of All Sed Time-Integrated SS Flux -----------
# ------------- Standard and Difference from Senstivity Tests ---------------
# ---------------------------------------------------------------------------
# Make a plots with the time-integrated fluxes for transect 1 with the first 
# subplot being the standard run, the middle being no waves minus standard,
# and the last being double waves minus standard run 

# Set the colormap
cmap_ssflux_diff = cmocean.cm.balance

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_diff = np.arange(-1000, 1000, 10) 
vmin14 = -2000
vmax14 = 2000

# Make the figure
fig14, ax14 = plt.subplots(3, figsize=(22,21))

# Set the title
fig14.suptitle('Time-Integrated Suspended Sediment Flux \nAll Sed - Transect 1', x=0.5, y=0.93, fontsize=fontsize)

# Standard
# Plot the ss flux
cs38 = ax14[0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_std[:,41:165], time_int_ssflux_tran1_std[:,41:165], 
               vmin=vmin5, vmax=vmax5, cmap=cmap_ssflux_std)
# Set the y lim
ax14[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax14[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax14[0].set_title('Standard Run', fontsize=fontsize-1)
# Make background gray
ax14[0].set_facecolor('lightgray')
# Add a colorbar
cbar14_ax = fig14.add_axes([0.92, 0.60, 0.02, 0.23]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar14 = plt.colorbar(cs38, orientation='vertical', cax=cbar14_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# No waves minus standard
# Plot the ss flux 
cs39 = ax14[1].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_nowaves[:,41:165], (time_int_ssflux_tran1_nowaves[:,41:165]-time_int_ssflux_tran1_std[:,41:165]), 
               vmin=vmin14, vmax=vmax14, cmap=cmap_ssflux_diff)
# Set the y lim
ax14[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
ax14[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax14[1].set_title('No Waves Run Minus Standard Run', fontsize=fontsize-1)
# Make background gray
ax14[1].set_facecolor('lightgray')

# Double waves minus standard
# Plot the ss flux 
cs40 = ax14[2].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_double_waves[:,41:165], (time_int_ssflux_tran1_double_waves[:,41:165]-time_int_ssflux_tran1_std[:,41:165]), 
               vmin=vmin14, vmax=vmax14, cmap=cmap_ssflux_diff)
# Set the y lim
ax14[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax14[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax14[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax14[2].set_title('Double Waves Run Minus Standard Run', fontsize=fontsize-1)
# Make background gray
ax14[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar14b_ax = fig14.add_axes([0.92, 0.11, 0.02, 0.45]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar14b = plt.colorbar(cs40, orientation='vertical', cax=cbar14b_ax, extend='both').set_label(label='Difference in Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.25) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_std_diff_tran1_60m_0003.png', bbox_inches='tight')


                          
# ---------------------------------------------------------------------------
# ---------- Plot 15: Transect 2 of All Sed Time-Integrated SS Flux -----------
# ------------- Standard and Difference from Senstivity Tests ---------------
# ---------------------------------------------------------------------------
# Make a plots with the time-integrated fluxes for transect 2 with the first 
# subplot being the standard run, the middle being no waves minus standard,
# and the last being double waves minus standard run 

# Set the colormap
cmap_ssflux_diff = cmocean.cm.balance

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_diff = np.arange(-1000, 1000, 10) 
vmin15 = -2000
vmax15 = 2000

# Make the figure
fig15, ax15 = plt.subplots(3, figsize=(22,21))

# Set the title
fig15.suptitle('Time-Integrated Suspended Sediment Flux \nAll Sed - Transect 2', x=0.5, y=0.93, fontsize=fontsize)

# Standard
# Plot the ss flux
cs41 = ax15[0].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_std[:,17:150], time_int_ssflux_tran2_std[:,17:150], 
               vmin=vmin5, vmax=vmax5, cmap=cmap_ssflux_std)
# Set the y lim
ax15[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax15[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax15[0].set_title('Standard Run', fontsize=fontsize-1)
# Make background gray
ax15[0].set_facecolor('lightgray')
# Add a colorbar
cbar15_ax = fig15.add_axes([0.92, 0.60, 0.02, 0.25]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar15 = plt.colorbar(cs41, orientation='vertical', cax=cbar15_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# No waves minus standard
# Plot the ss flux 
cs42 = ax15[1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_nowaves[:,17:150], (time_int_ssflux_tran2_nowaves[:,17:150]-time_int_ssflux_tran2_std[:,17:150]), 
               vmin=vmin15, vmax=vmax15, cmap=cmap_ssflux_diff)
# Set the y lim
ax15[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
ax15[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax15[1].set_title('No Waves Run Minus Standard Run', fontsize=fontsize-1)
# Make background gray
ax15[1].set_facecolor('lightgray')

# Double waves minus standard
# Plot the ss flux 
cs43 = ax15[2].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_double_waves[:,17:150], (time_int_ssflux_tran2_double_waves[:,17:150]-time_int_ssflux_tran2_std[:,17:150]), 
               vmin=vmin15, vmax=vmax15, cmap=cmap_ssflux_diff)
# Set the y lim
ax15[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax15[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax15[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax15[2].set_title('Double Waves Run Minus Standard Run', fontsize=fontsize-1)
# Make background gray
ax15[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar15b_ax = fig15.add_axes([0.92, 0.11, 0.02, 0.38]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar15b = plt.colorbar(cs43, orientation='vertical', cax=cbar15b_ax, extend='both').set_label(label='Difference in Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.25) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_std_diff_tran2_60m_0003.png', bbox_inches='tight')



# ---------------------------------------------------------------------------
# ---------- Plot 16: Transect 3 of All Sed Time-Integrated SS Flux -----------
# ------------- Standard and Difference from Senstivity Tests ---------------
# ---------------------------------------------------------------------------
# Make a plots with the time-integrated fluxes for transect 3 with the first 
# subplot being the standard run, the middle being no waves minus standard,
# and the last being double waves minus standard run 

# Set the colormap
cmap_ssflux_diff = cmocean.cm.balance

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
lev_ssflux_diff = np.arange(-1000, 1000, 10) 
vmin16 = -2000
vmax16 = 2000

# Make the figure
fig16, ax16 = plt.subplots(3, figsize=(22,21))

# Set the title
fig16.suptitle('Time-Integrated Suspended Sediment Flux \nAll Sed - Transect 3', x=0.5, y=0.93, fontsize=fontsize)

# Standard
# Plot the ss flux
cs44 = ax16[0].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_std[:,12:125], time_int_ssflux_tran3_std[:,12:125], 
               vmin=vmin5, vmax=vmax5, cmap=cmap_ssflux_std)
# Set the y lim
ax16[0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax16[0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax16[0].set_title('Standard Run', fontsize=fontsize-1)
# Make background gray
ax16[0].set_facecolor('lightgray')
# Add a colorbar
cbar16_ax = fig16.add_axes([0.92, 0.62, 0.02, 0.23]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar16 = plt.colorbar(cs44, orientation='vertical', cax=cbar16_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# No waves minus standard
# Plot the ss flux 
cs45 = ax16[1].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_nowaves[:,12:125], (time_int_ssflux_tran3_nowaves[:,12:125]-time_int_ssflux_tran3_std[:,12:125]), 
               vmin=vmin16, vmax=vmax16, cmap=cmap_ssflux_diff)
# Set the y lim
ax16[1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
ax16[1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax16[1].set_title('No Waves Run Minus Standard Run', fontsize=fontsize-1)
# Make background gray
ax16[1].set_facecolor('lightgray')

# Double waves minus standard
# Plot the ss flux 
cs46 = ax16[2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_double_waves[:,12:125], (time_int_ssflux_tran3_double_waves[:,12:125]-time_int_ssflux_tran3_std[:,12:125]), 
               vmin=vmin16, vmax=vmax16, cmap=cmap_ssflux_diff)
# Set the y lim
ax16[2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax16[2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax16[2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax16[2].set_title('Double Waves Run Minus Standard Run', fontsize=fontsize-1)
# Make background gray
ax16[2].set_facecolor('lightgray')

# Make and label the colorbar
cbar16b_ax = fig16.add_axes([0.92, 0.11, 0.02, 0.47]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar16b = plt.colorbar(cs46, orientation='vertical', cax=cbar16b_ax, extend='both').set_label(label='Difference in Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.25) #0.08

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_std_diff_tran3_60m_0003.png', bbox_inches='tight')



# ---------------------------------------------------------------------------
# ---------- Plot 17: All Transects of All Sed Time-Integrated SS Flux -----------
# -------------------------------- All runs -----------------------------------
# ---------------------------------------------------------------------------
# Make a plots with the time-integrated fluxes for all transects (rows) and 
# all runs (columns)

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
vmin17 = -5000
vmax17 = 5000

# Make the figure
fig17, ax17 = plt.subplots(3, 3, figsize=(42,20)) # (35,25)

# Set the title
#fig17.suptitle('Time-Integrated Suspended Sediment Flux', x=0.5, y=0.93, fontsize=fontsize)

# Standard - Transect 1
# Plot the ss flux
cs47 = ax17[0,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_std[:,41:165], time_int_ssflux_tran1_std[:,41:165], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_std)
# Set the y lim
ax17[0,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax17[0,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax17[0,0].set_title('Transect 1', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax17[0,0].set_facecolor('lightgray')

# No waves - Transect 1
# Plot the ss flux
cs50 = ax17[1,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_nowaves[:,41:165], time_int_ssflux_tran1_nowaves[:,41:165], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax17[1,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax17[1,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,0].set_title('No Waves', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax17[1,0].set_facecolor('lightgray')

# Double waves - Transect 1
# Plot the ss flux
cs53 = ax17[2,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_double_waves[:,41:165], time_int_ssflux_tran1_double_waves[:,41:165], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax17[2,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax17[2,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax17[2,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,0].set_title('Double Waves', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax17[2,0].set_facecolor('lightgray')


# Standard - Transect 2
# Plot the ss flux 
cs48 = ax17[0,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_std[:,17:150], time_int_ssflux_tran2_std[:,17:150], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_std)
# Set the y lim
ax17[0,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[1,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax17[0,1].set_title('Transect 2', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax17[0,1].set_facecolor('lightgray')

# No waves - Transect 2
# Plot the ss flux 
cs51 = ax17[1,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_nowaves[:,17:150], time_int_ssflux_tran2_nowaves[:,17:150], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax17[1,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[1,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,2].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax17[1,1].set_facecolor('lightgray')

# Double waves - Transect 2
# Plot the ss flux 
cs54 = ax17[2,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_double_waves[:,17:150], time_int_ssflux_tran2_double_waves[:,17:150], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax17[2,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[2,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax17[2,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,2].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax17[2,1].set_facecolor('lightgray')


# Standard - Transect 3
# Plot the ss flux 
cs49 = ax17[0,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_std[:,12:125], time_int_ssflux_tran3_std[:,12:125], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_std)
# Set the y lim
ax17[0,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax17[0,2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax17[0,2].set_title('Transect 3', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax17[0,2].set_facecolor('lightgray')

# No waves - Transect 3
# Plot the ss flux 
cs52 = ax17[1,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_nowaves[:,12:125], time_int_ssflux_tran3_nowaves[:,12:125], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax17[1,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax17[1,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax17[1,2].set_facecolor('lightgray')

# Double waves - Transect 3
# Plot the ss flux 
cs55 = ax17[2,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_double_waves[:,12:125], time_int_ssflux_tran3_double_waves[:,12:125], 
               vmin=vmin17, vmax=vmax17, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax17[2,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax17[2,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax17[2,2].set_facecolor('lightgray')

# Make and label the colorbar
cbar17_ax = fig17.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar17 = plt.colorbar(cs55, orientation='vertical', cax=cbar17_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Add labels for the rows
plt.text(0.022, 0.769, 'Standard', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.485, 'No Waves', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.218, 'Double Waves', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.18, wspace=0.18) #0.08

# Add subplot labels (bottom left corners)
plt.text(0.131, 0.664, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.664, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.664, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.395, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.395, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.395, 'f)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.127, 'g)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.127, 'h)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.127, 'i)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_allruns_tran3_60m_0009.png', bbox_inches='tight')
                          

    
# -------------------------------------------------------------------------------
# ---- Plot 18: 2D Plot of All Sed Time-Integrated, Depth-Integrated SS Flux ------
# --------------------------------------------------------------------------------
# Make a plot of the time-integrated, depth-integrated suspended sediment flux
# for all runs - with first column being different transects, one color for
# each run and second column being sediment class by run for that transect

# Set the linewidth
lw18 = 4

# Make the figure
fig18, ax18 = plt.subplots(2, 3, figsize=(42,20)) # (22,21)

# Set the title
#fig18.suptitle('Time-Integrated, Depth-Integrated Suspended Sediment Flux', x=0.5, y=0.95, fontsize=fontsize)

# First subplot = one line for each run at Transect 1, all sediment 
ax18[0,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran1_std, linewidth=lw18, color='#17BECF', label='Standard')
ax18[0,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran1_nowaves, linewidth=lw18, color='#8C564B', label='No Waves')
ax18[0,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran1_double_waves, linewidth=lw18, color='#9467BD', label='Double Waves')
ax18[0,0].axhline(y=0, color='black', linewidth=1)
ax18[0,0].legend(fontsize=fontsize-2)
ax18[0,0].set_title('Transect 1', fontsize=fontsize-1, fontweight='bold')
ax18[0,0].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)

# Second subplot = one line for each run at Transect 2, all sediment 
ax18[0,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran2_std, linewidth=lw18, color='#17BECF', label='Standard')
ax18[0,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran2_nowaves, linewidth=lw18, color='#8C564B', label='No Waves')
ax18[0,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran2_double_waves, linewidth=lw18, color='#9467BD', label='Double Waves')
ax18[0,1].axhline(y=0, color='black', linewidth=1)
ax18[0,1].legend(fontsize=fontsize-2)
ax18[0,1].set_title('Transect 2', fontsize=fontsize-1, fontweight='bold')

# Third subplots = one line for each run at Transect 3, all sediment 
ax18[0,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran3_std, linewidth=lw18, color='#17BECF', label='Standard')
ax18[0,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran3_nowaves, linewidth=lw18, color='#8C564B', label='No Waves')
ax18[0,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_tran3_double_waves, linewidth=lw18, color='#9467BD', label='Double Waves')
ax18[0,2].axhline(y=0, color='black', linewidth=1)
ax18[0,2].legend(fontsize=fontsize-2)
ax18[0,2].set_title('Transect 3', fontsize=fontsize-1, fontweight='bold')

# Fourth subplots = one line for each sediment class by run for transect 1
# Standard
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran1_std, linestyle='solid', linewidth=lw18, color='#17BECF', label='Standard, Mud01')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran1_std, linestyle='dashed', linewidth=lw18, color='#17BECF', label='Standard, Mud02')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran1_std, linestyle='dotted', linewidth=lw18, color='#17BECF', label='Standard, Sand01')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran1_std, linestyle='dashdot', linewidth=lw18, color='#17BECF', label='Standard, Sand02')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran1_std, linestyle='', marker='+', linewidth=lw18, color='#17BECF', label='Standard, Sand03')
# No waves
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran1_nowaves, linestyle='solid', linewidth=lw18, color='#8C564B', label='No Waves, Mud01')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran1_nowaves, linestyle='dashed', linewidth=lw18, color='#8C564B', label='No Waves, Mud02')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran1_nowaves, linestyle='dotted', linewidth=lw18, color='#8C564B', label='No Waves, Sand01')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran1_nowaves, linestyle='dashdot', linewidth=lw18, color='#8C564B', label='No Waves, Sand02')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran1_nowaves, linestyle='', marker='+', linewidth=lw18, color='#8C564B', label='No Waves, Sand03')
# Double waves
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran1_double_waves, linestyle='solid', linewidth=lw18, color='#9467BD', label='Double Waves, Mud01')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran1_double_waves, linestyle='dashed', linewidth=lw18, color='#9467BD', label='Double Waves, Mud02')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran1_double_waves, linestyle='dotted', linewidth=lw18, color='#9467BD', label='Double Waves, Sand01')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran1_double_waves, linestyle='dashdot', linewidth=lw18, color='#9467BD', label='Double Waves, Sand02')
ax18[1,0].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran1_double_waves, linestyle='', marker='+', linewidth=lw18, color='#9467BD', label='Double Waves, Sand03')
#ax18[1,0].set_title('Transect 1', fontsize=fontsize-1)
ax18[1,0].axhline(y=0, color='black', linewidth=1)
ax18[1,0].set_xlabel('Y (km)', fontsize=fontsize-2)
ax18[1,0].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)

# Fifth subplots = one line for each sediment class by run for transect 2
# Standard
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran2_std, linestyle='solid', linewidth=lw18, color='#17BECF', label='Standard, Mud01')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran2_std, linestyle='dashed', linewidth=lw18, color='#17BECF', label='Standard, Mud02')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran2_std, linestyle='dotted', linewidth=lw18, color='#17BECF', label='Standard, Sand01')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran2_std, linestyle='dashdot', linewidth=lw18, color='#17BECF', label='Standard, Sand02')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran2_std, linestyle='', marker='+', linewidth=lw18, color='#17BECF', label='Standard, Sand03')
# No waves
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran2_nowaves, linestyle='solid', linewidth=lw18, color='#8C564B', label='No Waves, Mud01')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran2_nowaves, linestyle='dashed', linewidth=lw18, color='#8C564B', label='No Waves, Mud02')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran2_nowaves, linestyle='dotted', linewidth=lw18, color='#8C564B', label='No Waves, Sand01')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran2_nowaves, linestyle='dashdot', linewidth=lw18, color='#8C564B', label='No Waves, Sand02')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran2_nowaves, linestyle='', marker='+', linewidth=lw18, color='#8C564B', label='No Waves, Sand03')
# Double waves
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran2_double_waves, linestyle='solid', linewidth=lw18, color='#9467BD', label='Double Waves, Mud01')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran2_double_waves, linestyle='dashed', linewidth=lw18, color='#9467BD', label='Double Waves, Mud02')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran2_double_waves, linestyle='dotted', linewidth=lw18, color='#9467BD', label='Double Waves, Sand01')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran2_double_waves, linestyle='dashdot', linewidth=lw18, color='#9467BD', label='Double Waves, Sand02')
ax18[1,1].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran2_double_waves, linestyle='', marker='+', linewidth=lw18, color='#9467BD', label='Double Waves, Sand03')
#ax18[1,1].set_title('Transect 2', fontsize=fontsize-1)
ax18[1,1].axhline(y=0, color='black', linewidth=1)
ax18[1,1].set_xlabel('Y (km)', fontsize=fontsize-2)

# Sixth subplots = one line for each sediment class by run for transect 3
# Standard
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran3_std, linestyle='solid', linewidth=lw18, color='#17BECF', label='Standard, Mud01')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran3_std, linestyle='dashed', linewidth=lw18, color='#17BECF', label='Standard, Mud02')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran3_std, linestyle='dotted', linewidth=lw18, color='#17BECF', label='Standard, Sand01')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran3_std, linestyle='dashdot', linewidth=lw18, color='#17BECF', label='Standard, Sand02')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran3_std, linestyle='', marker='+', linewidth=lw18, color='#17BECF', label='Standard, Sand03')
# No waves
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran3_nowaves, linestyle='solid', linewidth=lw18, color='#8C564B', label='No Waves, Mud01')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran3_nowaves, linestyle='dashed', linewidth=lw18, color='#8C564B', label='No Waves, Mud02')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran3_nowaves, linestyle='dotted', linewidth=lw18, color='#8C564B', label='No Waves, Sand01')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran3_nowaves, linestyle='dashdot', linewidth=lw18, color='#8C564B', label='No Waves, Sand02')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran3_nowaves, linestyle='', marker='+', linewidth=lw18, color='#8C564B', label='No Waves, Sand03')
# Double waves
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud01_tran3_double_waves, linestyle='solid', linewidth=lw18, color='#9467BD', label='Double Waves, Mud01')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_mud02_tran3_double_waves, linestyle='dashed', linewidth=lw18, color='#9467BD', label='Double Waves, Mud02')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand01_tran3_double_waves, linestyle='dotted', linewidth=lw18, color='#9467BD', label='Double Waves, Sand01')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand02_tran3_double_waves, linestyle='dashdot', linewidth=lw18, color='#9467BD', label='Double Waves, Sand02')
ax18[1,2].plot(y_rho_flat/1000, time_int_depth_int_ssflux_sand03_tran3_double_waves, linestyle='', marker='+', linewidth=lw18, color='#9467BD', label='Double Waves, Sand03')
#ax18[2,1].set_title('Transect 3', fontsize=fontsize-1)
ax18[1,2].axhline(y=0, color='black', linewidth=1)
ax18[1,2].set_xlabel('Y (km)', fontsize=fontsize-2)

# Label the plot
#plt.setp(ax12[0].get_xticklabels(), visible=False)
#ax12[0].set_ylabel('Time-Integrated, \nDepth-Integrated \nSuspended Sediment \nFlux (kg/m)', fontsize=fontsize-2, rotation=0, va='center', labelpad=120)
# Add a legend
#ax12[0].legend(fontsize=fontsize-2)

# Add a legend
ax18[1,1].legend(bbox_to_anchor=(1.25, -0.15), ncol=3, fontsize=fontsize-2) #(0.2, -0.45)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.22, wspace=0.25) #0.08

# Add subplot labels (bottom right corners)
plt.text(0.331, 0.548, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.610, 0.548, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.887, 0.548, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.331, 0.123, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.610, 0.123, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.887, 0.123, 'f)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_allruns_2d_3transects_byclass_008.png', bbox_inches='tight')


# Print some stats
# 20 km out is 33 or 34 index in x-axis (transect 2)
# 25 km out is 41 or 42 index in x-axis (transect 1)
# Print ratio of unaggregated mud flux divided by aggregated mud flux
# Transect 1
# Find the real index I want 
min_time_int_depth_int_ssflux_mud01_tran1_double_waves = np.nanmin(time_int_depth_int_ssflux_mud01_tran1_double_waves)
print('index of minimum 2d flux transect 1, double waves, mud01: ', np.where(time_int_depth_int_ssflux_mud01_tran1_double_waves == min_time_int_depth_int_ssflux_mud01_tran1_double_waves))
# Double wave run
# 25 km out take A
#print('transect 1 unaggregated time-integrated, depth-integrated flux 25 km out, double waves (index 41): ', time_int_depth_int_ssflux_mud01_tran1_double_waves[41]/time_int_depth_int_ssflux_mud02_tran1_double_waves[41])
# 25 km out take B
print('transect 1 unaggregated time-integrated, depth-integrated flux 25 km out, double waves (index 42): ', time_int_depth_int_ssflux_mud01_tran1_double_waves[42]/time_int_depth_int_ssflux_mud02_tran1_double_waves[42])
# No Waves run 
# 25 km out take A
#print('transect 1 unaggregated time-integrated, depth-integrated flux 25 km out, no waves (index 41): ', time_int_depth_int_ssflux_mud01_tran1_nowaves[41]/time_int_depth_int_ssflux_mud02_tran1_nowaves[41])
# 25 km out take B
print('transect 1 unaggregated time-integrated, depth-integrated flux 25 km out, no waves (index 42): ', time_int_depth_int_ssflux_mud01_tran1_nowaves[42]/time_int_depth_int_ssflux_mud02_tran1_nowaves[42])

# Transect 2
# Find the real index I want 
min_time_int_depth_int_ssflux_mud01_tran2_double_waves = np.nanmin(time_int_depth_int_ssflux_mud01_tran2_double_waves)
print('index of minimum 2d flux transect 2, double waves, mud01: ', np.where(time_int_depth_int_ssflux_mud01_tran2_double_waves == min_time_int_depth_int_ssflux_mud01_tran2_double_waves))
# Double wave run
# 20 km out Correct
print('transect 2 unaggregated time-integrated, depth-integrated flux 20 km out, double waves (index 31): ', time_int_depth_int_ssflux_mud01_tran2_double_waves[31]/time_int_depth_int_ssflux_mud02_tran2_double_waves[31])
# 20 km out take A
#print('transect 2 unaggregated time-integrated, depth-integrated flux 20 km out, double waves (index 33): ', time_int_depth_int_ssflux_mud01_tran2_double_waves[33]/time_int_depth_int_ssflux_mud02_tran2_double_waves[33])
# 20 km out take B
#print('transect 2 unaggregated time-integrated, depth-integrated flux 20 km out, double waves (index 34): ', time_int_depth_int_ssflux_mud01_tran2_double_waves[34]/time_int_depth_int_ssflux_mud02_tran2_double_waves[34])
# No Waves run 
# 20 km out take Correct
print('transect 2 unaggregated time-integrated, depth-integrated flux 20 km out, no waves (index 31): ', time_int_depth_int_ssflux_mud01_tran2_nowaves[31]/time_int_depth_int_ssflux_mud02_tran2_nowaves[31])
# 20 km out take A
#print('transect 2 unaggregated time-integrated, depth-integrated flux 20 km out, no waves (index 33): ', time_int_depth_int_ssflux_mud01_tran2_nowaves[33]/time_int_depth_int_ssflux_mud02_tran2_nowaves[33])
# 20 km out take B
#print('transect 2 unaggregated time-integrated, depth-integrated flux 20 km out, no waves (index 34): ', time_int_depth_int_ssflux_mud01_tran2_nowaves[34]/time_int_depth_int_ssflux_mud02_tran2_nowaves[34])


                          
                          
# ---------------------------------------------------------------------------
# ---------------------------- Plot 19: All Transects of  --------------------
# --------------------- All runs with Time-averaged currents ---------------------
# ---------------------------------------------------------------------------
# Make a plot of the time-averaged u currents at each transect for each type of run
# Then if they all look similar, add the standard run one onto the plot below

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
vmin19 = -0.1
vmax19 = 0.1

# Set the colormap 
#cmap_time_avg_u = cmocean.cm.diff
#cmap_time_avg_u = cmocean.cm.balance
cmap_time_avg_u = plt.get_cmap('PuOr_r')

# Make the figure
fig19, ax19 = plt.subplots(3, 3, figsize=(42,20)) # (35,25)

# Set the title
#fig19.suptitle('Time-Averaged U Currents', x=0.5, y=0.93, fontsize=fontsize)

# Standard - Transect 1
# Plot the time-averaged u
cs56 = ax19[0,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_std[:,41:165], time_avg_u_tran1_std[:,41:165], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[0,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax19[0,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax19[0,0].set_title('Transect 1', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax19[0,0].set_facecolor('lightgray')

# No waves - Transect 1
# Plot the time-averaged u
cs57 = ax19[1,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_nowaves[:,41:165], time_avg_u_tran1_nowaves[:,41:165], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[1,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax19[1,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,0].set_title('No Waves', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax19[1,0].set_facecolor('lightgray')

# Double waves - Transect 1
# Plot the time-averaged u
cs58 = ax19[2,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_double_waves[:,41:165], time_avg_u_tran1_double_waves[:,41:165], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[2,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax19[2,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax19[2,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,0].set_title('Double Waves', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax19[2,0].set_facecolor('lightgray')


# Standard - Transect 2
# Plot the time-averaged u
cs59 = ax19[0,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_std[:,17:150], time_avg_u_tran2_std[:,17:150], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[0,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[1,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
ax19[0,1].set_title('Transect 2', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax19[0,1].set_facecolor('lightgray')

# No waves - Transect 2
# Plot the time-averaged u
cs60 = ax19[1,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_nowaves[:,17:150], time_avg_u_tran2_nowaves[:,17:150], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[1,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[1,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,2].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax19[1,1].set_facecolor('lightgray')

# Double waves - Transect 2
# Plot the time-averaged u
cs61 = ax19[2,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_double_waves[:,17:150], time_avg_u_tran2_double_waves[:,17:150], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[2,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[2,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax19[2,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,2].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax19[2,1].set_facecolor('lightgray')


# Standard - Transect 3
# Plot the time-averaged u
cs62 = ax19[0,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_std[:,12:125], time_avg_u_tran3_std[:,12:125], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[0,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax17[0,2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax19[0,2].set_title('Transect 3', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax19[0,2].set_facecolor('lightgray')

# No waves - Transect 3
# Plot the time-averaged u 
cs63 = ax19[1,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_nowaves[:,12:125], time_avg_u_tran3_nowaves[:,12:125], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[1,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax17[1,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax19[1,2].set_facecolor('lightgray')

# Double waves - Transect 3
# Plot the time-averaged u
cs64 = ax19[2,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_double_waves[:,12:125], time_avg_u_tran3_double_waves[:,12:125], 
               vmin=vmin19, vmax=vmax19, cmap=cmap_time_avg_u)
# Set the y lim
ax19[2,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax19[2,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax19[2,2].set_facecolor('lightgray')

# Make and label the colorbar
cbar19_ax = fig19.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar19 = plt.colorbar(cs56, orientation='vertical', cax=cbar19_ax, extend='both').set_label(label='Current Speed (m/s)', size=fontsize-2, labelpad=10)

# Add labels for the rows
plt.text(0.022, 0.769, 'Standard', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.485, 'No Waves', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.218, 'Double Waves', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.18, wspace=0.18) #0.08

# Add subplot labels (bottom left corners)
plt.text(0.131, 0.664, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.664, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.664, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.395, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.395, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.395, 'f)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.127, 'g)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.127, 'h)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.127, 'i)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/time_avg_u_allruns_tran3_60m_0005.png', bbox_inches='tight')

        

# ---------------------------------------------------------------------------
# ---------- Plot 20: All Transects of All Sed Time-Integrated SS Flux -----------
# --------------------- All runs with Time-averaged currents ---------------------
# ---------------------------------------------------------------------------
# Make a plots with the time-integrated fluxes for all transects (rows) and 
# all runs (columns) and include the time-averaged currents 

# Set the levels 
#lev_ssflux_std = np.arange(-20, 20, 1) 
# Sediment fluxes
vmin20a = -5000
vmax20a = 5000
# Currents
vmin20b = -0.1
vmax20b = 0.1

# Make the figure
fig20, ax20 = plt.subplots(4, 3, figsize=(42,30)) # (35,25)

# Set the title
#fig17.suptitle('Time-Integrated Suspended Sediment Flux', x=0.5, y=0.93, fontsize=fontsize)

# No waves - Transect 1
# Plot the ss flux
cs66 = ax20[0,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_nowaves[:,41:165], time_int_ssflux_tran1_nowaves[:,41:165], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax20[0,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax20[0,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax20[0,0].set_title('Transect 1', fontsize=fontsize+1, fontweight='bold')
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,0].set_title('No Waves', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax20[0,0].set_facecolor('lightgray')

# Standard - Transect 1
# Plot the ss flux
cs65 = ax20[1,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_std[:,41:165], time_int_ssflux_tran1_std[:,41:165], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_std)
# Set the y lim
ax20[1,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[0].get_xticklabels(), visible=False)
ax20[1,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax20[1,0].set_title('Transect 1', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax20[1,0].set_facecolor('lightgray')

# Double waves - Transect 1
# Plot the ss flux
cs67 = ax20[2,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_double_waves[:,41:165], time_int_ssflux_tran1_double_waves[:,41:165], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax20[2,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax20[2,0].get_xticklabels(), visible=False)
ax20[2,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax20[2,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,0].set_title('Double Waves', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax20[2,0].set_facecolor('lightgray')



# No waves - Transect 2
# Plot the ss flux 
cs69 = ax20[0,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_nowaves[:,17:150], time_int_ssflux_tran2_nowaves[:,17:150], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax20[0,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[1,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax20[0,1].set_title('Transect 2', fontsize=fontsize+1, fontweight='bold')
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,2].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax20[0,1].set_facecolor('lightgray')

# Standard - Transect 2
# Plot the ss flux 
cs68 = ax20[1,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_std[:,17:150], time_int_ssflux_tran2_std[:,17:150], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_std)
# Set the y lim
ax20[1,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax17[1,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax2.set_xlabel('Y (km)', fontsize=fontsize-2)
#ax20[1,1].set_title('Transect 2', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax20[1,1].set_facecolor('lightgray')

# Double waves - Transect 2
# Plot the ss flux 
cs70 = ax20[2,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_double_waves[:,17:150], time_int_ssflux_tran2_double_waves[:,17:150], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax20[2,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax20[2,1].get_xticklabels(), visible=False)
#ax17[2,1].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax20[2,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[1,2].set_title('Transect 2', fontsize=fontsize-1)
# Make background gray
ax20[2,1].set_facecolor('lightgray')



# No waves - Transect 3
# Plot the ss flux 
cs72 = ax20[0,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_nowaves[:,12:125], time_int_ssflux_tran3_nowaves[:,12:125], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_nowaves)
# Set the y lim
ax20[0,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax17[1,2].set_xlabel('Y (km)', fontsize=fontsize-2)
ax20[0,2].set_title('Transect 3', fontsize=fontsize+1, fontweight='bold')
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax20[0,2].set_facecolor('lightgray')

# Standard - Transect 3
# Plot the ss flux 
cs71 = ax20[1,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_std[:,12:125], time_int_ssflux_tran3_std[:,12:125], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_std)
# Set the y lim
ax20[1,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax17[0,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax20[1,2].set_title('Transect 3', fontsize=fontsize-1, fontweight='bold')
# Make background gray
ax20[1,2].set_facecolor('lightgray')

# Double waves - Transect 3
# Plot the ss flux 
cs73 = ax20[2,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_double_waves[:,12:125], time_int_ssflux_tran3_double_waves[:,12:125], 
               vmin=vmin20a, vmax=vmax20a, cmap=cmap_ssflux_double_waves)
# Set the y lim
ax20[2,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax20[2,2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
#ax20[2,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax20[2,2].set_facecolor('lightgray')

# Transect 1 Time-Averaged U Currents - std
cs74 = ax20[3,0].pcolormesh(y_rho_flat[41:165]/1000, time_avg_z_rho_tran1_std[:,41:165], time_avg_u_tran1_std[:,41:165], 
               vmin=vmin20b, vmax=vmax20b, cmap=cmap_time_avg_u)
# Set the y lim
ax20[3,0].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
ax20[3,0].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax20[3,0].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax20[3,0].set_facecolor('lightgray')

# Transect 2 Time-Averaged U Currents - std
cs75 = ax20[3,1].pcolormesh(y_rho_flat[17:150]/1000, time_avg_z_rho_tran2_std[:,17:150], time_avg_u_tran2_std[:,17:150], 
               vmin=vmin20b, vmax=vmax20b, cmap=cmap_time_avg_u)
# Set the y lim
ax20[3,1].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax20[3,1].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax20[3,1].set_facecolor('lightgray')

# Transect 3 Time-Averaged U Currents - std
cs76 = ax20[3,2].pcolormesh(y_rho_flat[12:125]/1000, time_avg_z_rho_tran3_std[:,12:125], time_avg_u_tran3_std[:,12:125], 
               vmin=vmin20b, vmax=vmax20b, cmap=cmap_time_avg_u)
# Set the y lim
ax20[3,2].set_ylim([-60,0])
# Label the plot
#plt.setp(ax2[2].get_xticklabels(), visible=False)
#ax17[2,2].set_ylabel('Depth (m)', fontsize=fontsize-2)
ax20[3,2].set_xlabel('Y (km)', fontsize=fontsize-2)
#ax17[2,2].set_title('Transect 3', fontsize=fontsize-1)
# Make background gray
ax20[3,2].set_facecolor('lightgray')


# Make and label the colorbar
cbar20_ax = fig20.add_axes([0.92, 0.32, 0.02, 0.56]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar20 = plt.colorbar(cs65, orientation='vertical', cax=cbar20_ax, extend='both').set_label(label='Suspended Sediment Flux (kg/m\u00b2)', size=fontsize-2, labelpad=10)

# Make and label colorbar for currents 
cbar20b_ax = fig20.add_axes([0.92, 0.11, 0.02, 0.17]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar20b = plt.colorbar(cs74, orientation='vertical', cax=cbar20b_ax, extend='both').set_label(label='Current Speed (m/s)', size=fontsize-2, labelpad=10)

# Add labels for the rows
plt.text(0.022, 0.802, 'No Waves', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.593, 'Standard', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.392, 'Double Waves', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.022, 0.200, 'East-West \nCurrents', fontsize=fontsize-1, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Adjust subplots spacing 
plt.subplots_adjust(hspace=0.18, wspace=0.18) #0.08

# Add subplot labels (bottom left corners)
plt.text(0.131, 0.721, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.721, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.721, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.520, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.520, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.520, 'f)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.322, 'g)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.322, 'h)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.322, 'i)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.131, 0.121, 'j)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.404, 0.121, 'k)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.677, 0.121, 'l)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the figure
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/time_integrated_flow/ssflux_allsed_allruns_tran3_time_avg_u_std_60m_0007.png', bbox_inches='tight')        
        
        

        
                        
# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Make an s_rho 
s_rho = np.arange(0,20,1)

# 3D transect
# Set up the data
roms_time_int_ssflux_transects_allsed = xr.Dataset(
    data_vars=dict(
        time_int_ssflux_tran1_nowaves=(['s_rho','y'], time_int_ssflux_tran1_nowaves),
        time_int_ssflux_tran2_nowaves=(['s_rho','y'], time_int_ssflux_tran2_nowaves),
        time_int_ssflux_tran3_nowaves=(['s_rho','y'], time_int_ssflux_tran3_nowaves),
        time_int_ssflux_tran1_std=(['s_rho','y'], time_int_ssflux_tran1_std),
        time_int_ssflux_tran2_std=(['s_rho','y'], time_int_ssflux_tran2_std),
        time_int_ssflux_tran3_std=(['s_rho','y'], time_int_ssflux_tran3_std),
        time_int_ssflux_tran1_double_waves=(['s_rho','y'], time_int_ssflux_tran1_double_waves),
        time_int_ssflux_tran2_double_waves=(['s_rho','y'], time_int_ssflux_tran2_double_waves),
        time_int_ssflux_tran3_double_waves=(['s_rho','y'], time_int_ssflux_tran3_double_waves),
        time_avg_u_tran1_std=(['s_rho','y'], time_avg_u_tran1_std),
        time_avg_u_tran2_std=(['s_rho','y'], time_avg_u_tran2_std),
        time_avg_u_tran3_std=(['s_rho','y'], time_avg_u_tran3_std),
        z_rho_tran1_nowaves=(['s_rho','y'], time_avg_z_rho_tran1_nowaves),
        z_rho_tran2_nowaves=(['s_rho','y'], time_avg_z_rho_tran2_nowaves),
        z_rho_tran3_nowaves=(['s_rho','y'], time_avg_z_rho_tran3_nowaves),
        z_rho_tran1_std=(['s_rho','y'], time_avg_z_rho_tran1_std),
        z_rho_tran2_std=(['s_rho','y'], time_avg_z_rho_tran2_std),
        z_rho_tran3_std=(['s_rho','y'], time_avg_z_rho_tran3_std),
        z_rho_tran1_double_waves=(['s_rho','y'], time_avg_z_rho_tran1_double_waves),
        z_rho_tran2_double_waves=(['s_rho','y'], time_avg_z_rho_tran2_double_waves),
        z_rho_tran3_double_waves=(['s_rho','y'], time_avg_z_rho_tran3_double_waves)
        ),
    coords=dict(
        y_full=('y', y_rho_flat), 
        s_rho=('s_rho', s_rho)
        ),
    attrs=dict(description='Time-integrated suspended sediment fluxes at three different transects for the no waves, standard, and double waves model runs')) 
# Add more metadata?
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran1_nowaves.name='time-integrated suspended sediment flux at transect 1 for no waves model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran2_nowaves.name='time-integrated suspended sediment flux at transect 2 for no waves model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran3_nowaves.name='time-integrated suspended sediment flux at transect 3 for no waves model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran1_std.name='time-integrated suspended sediment flux at transect 1 for standard model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran2_std.name='time-integrated suspended sediment flux at transect 2 for standard model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran3_std.name='time-integrated suspended sediment flux at transect 3 for standard model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran1_double_waves.name='time-integrated suspended sediment flux at transect 1 for double waves model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran2_double_waves.name='time-integrated suspended sediment flux at transect 2 for double waves model run'
roms_time_int_ssflux_transects_allsed.time_int_ssflux_tran3_double_waves.name='time-integrated suspended sediment flux at transect 3 for double waves model run'
roms_time_int_ssflux_transects_allsed.time_avg_u_tran1_std.name='time-averaged u currents at transect 1 for the standard model run'
roms_time_int_ssflux_transects_allsed.time_avg_u_tran2_std.name='time-averaged u currents at transect 2 for the standard model run'
roms_time_int_ssflux_transects_allsed.time_avg_u_tran3_std.name='time-averaged u currents at transect 3 for the standard model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran1_nowaves.name='z-rho coordinates at transect 1 for no waves model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran2_nowaves.name='z-rho coordinates at transect 2 for no waves model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran3_nowaves.name='z-rho coordinates at transect 3 for no waves model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran1_std.name='z-rho coordinates at transect 1 for standard model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran2_std.name='z-rho coordinates at transect 2 for standard model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran3_std.name='z-rho coordinates at transect 3 for standard model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran1_double_waves.name='z-rho coordinates at transect 1 for double waves model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran2_double_waves.name='z-rho coordinates at transect 2 for double waves model run'
roms_time_int_ssflux_transects_allsed.z_rho_tran3_double_waves.name='z-rho coordinates at transect 3 for double waves model run'
roms_time_int_ssflux_transects_allsed.y_full.name='latitudinal distance (meter)'
roms_time_int_ssflux_transects_allsed.s_rho.name='index for vertical level'

# Save to a netcdf
#roms_time_int_ssflux_transects_allsed.to_netcdf('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/fig11_roms_time_int_ssflux_transects_allsed.nc')


# 2D transect
# Set up the data
roms_time_int_depth_int_ssflux_transects_by_class = xr.Dataset(
    data_vars=dict(
        time_int_depth_int_ssflux_tran1_std=(['y'], time_int_depth_int_ssflux_tran1_std),
        time_int_depth_int_ssflux_tran2_std=(['y'], time_int_depth_int_ssflux_tran2_std),
        time_int_depth_int_ssflux_tran3_std=(['y'], time_int_depth_int_ssflux_tran3_std),
        time_int_depth_int_ssflux_tran1_nowaves=(['y'], time_int_depth_int_ssflux_tran1_nowaves),
        time_int_depth_int_ssflux_tran2_nowaves=(['y'], time_int_depth_int_ssflux_tran2_nowaves),
        time_int_depth_int_ssflux_tran3_nowaves=(['y'], time_int_depth_int_ssflux_tran3_nowaves),
        time_int_depth_int_ssflux_tran1_double_waves=(['y'], time_int_depth_int_ssflux_tran1_double_waves),
        time_int_depth_int_ssflux_tran2_double_waves=(['y'], time_int_depth_int_ssflux_tran2_double_waves),
        time_int_depth_int_ssflux_tran3_double_waves=(['y'], time_int_depth_int_ssflux_tran3_double_waves),
        time_int_depth_int_ssflux_mud01_tran1_std=(['y'], time_int_depth_int_ssflux_mud01_tran1_std),
        time_int_depth_int_ssflux_mud02_tran1_std=(['y'], time_int_depth_int_ssflux_mud02_tran1_std),
        time_int_depth_int_ssflux_sand01_tran1_std=(['y'], time_int_depth_int_ssflux_sand01_tran1_std),
        time_int_depth_int_ssflux_sand02_tran1_std=(['y'], time_int_depth_int_ssflux_sand02_tran1_std),
        time_int_depth_int_ssflux_sand03_tran1_std=(['y'], time_int_depth_int_ssflux_sand03_tran1_std),
        time_int_depth_int_ssflux_mud01_tran1_nowaves=(['y'], time_int_depth_int_ssflux_mud01_tran1_nowaves),
        time_int_depth_int_ssflux_mud02_tran1_nowaves=(['y'], time_int_depth_int_ssflux_mud02_tran1_nowaves),
        time_int_depth_int_ssflux_sand01_tran1_nowaves=(['y'], time_int_depth_int_ssflux_sand01_tran1_nowaves),
        time_int_depth_int_ssflux_sand02_tran1_nowaves=(['y'], time_int_depth_int_ssflux_sand02_tran1_nowaves),
        time_int_depth_int_ssflux_sand03_tran1_nowaves=(['y'], time_int_depth_int_ssflux_sand03_tran1_nowaves),
        time_int_depth_int_ssflux_mud01_tran1_double_waves=(['y'], time_int_depth_int_ssflux_mud01_tran1_double_waves),
        time_int_depth_int_ssflux_mud02_tran1_double_waves=(['y'], time_int_depth_int_ssflux_mud02_tran1_double_waves),
        time_int_depth_int_ssflux_sand01_tran1_double_waves=(['y'], time_int_depth_int_ssflux_sand01_tran1_double_waves),
        time_int_depth_int_ssflux_sand02_tran1_double_waves=(['y'], time_int_depth_int_ssflux_sand02_tran1_double_waves),
        time_int_depth_int_ssflux_sand03_tran1_double_waves=(['y'], time_int_depth_int_ssflux_sand03_tran1_double_waves),
        time_int_depth_int_ssflux_mud01_tran2_std=(['y'], time_int_depth_int_ssflux_mud01_tran2_std),
        time_int_depth_int_ssflux_mud02_tran2_std=(['y'], time_int_depth_int_ssflux_mud02_tran2_std),
        time_int_depth_int_ssflux_sand01_tran2_std=(['y'], time_int_depth_int_ssflux_sand01_tran2_std),
        time_int_depth_int_ssflux_sand02_tran2_std=(['y'], time_int_depth_int_ssflux_sand02_tran2_std),
        time_int_depth_int_ssflux_sand03_tran2_std=(['y'], time_int_depth_int_ssflux_sand03_tran2_std),
        time_int_depth_int_ssflux_mud01_tran2_nowaves=(['y'], time_int_depth_int_ssflux_mud01_tran2_nowaves),
        time_int_depth_int_ssflux_mud02_tran2_nowaves=(['y'], time_int_depth_int_ssflux_mud02_tran2_nowaves),
        time_int_depth_int_ssflux_sand01_tran2_nowaves=(['y'], time_int_depth_int_ssflux_sand01_tran2_nowaves),
        time_int_depth_int_ssflux_sand02_tran2_nowaves=(['y'], time_int_depth_int_ssflux_sand02_tran2_nowaves),
        time_int_depth_int_ssflux_sand03_tran2_nowaves=(['y'], time_int_depth_int_ssflux_sand03_tran2_nowaves),
        time_int_depth_int_ssflux_mud01_tran2_double_waves=(['y'], time_int_depth_int_ssflux_mud01_tran2_double_waves),
        time_int_depth_int_ssflux_mud02_tran2_double_waves=(['y'], time_int_depth_int_ssflux_mud02_tran2_double_waves),
        time_int_depth_int_ssflux_sand01_tran2_double_waves=(['y'], time_int_depth_int_ssflux_sand01_tran2_double_waves),
        time_int_depth_int_ssflux_sand02_tran2_double_waves=(['y'], time_int_depth_int_ssflux_sand02_tran2_double_waves),
        time_int_depth_int_ssflux_sand03_tran2_double_waves=(['y'], time_int_depth_int_ssflux_sand03_tran2_double_waves),
        time_int_depth_int_ssflux_mud01_tran3_std=(['y'], time_int_depth_int_ssflux_mud01_tran3_std),
        time_int_depth_int_ssflux_mud02_tran3_std=(['y'], time_int_depth_int_ssflux_mud02_tran3_std),
        time_int_depth_int_ssflux_sand01_tran3_std=(['y'], time_int_depth_int_ssflux_sand01_tran3_std),
        time_int_depth_int_ssflux_sand02_tran3_std=(['y'], time_int_depth_int_ssflux_sand02_tran3_std),
        time_int_depth_int_ssflux_sand03_tran3_std=(['y'], time_int_depth_int_ssflux_sand03_tran3_std),
        time_int_depth_int_ssflux_mud01_tran3_nowaves=(['y'], time_int_depth_int_ssflux_mud01_tran3_nowaves),
        time_int_depth_int_ssflux_mud02_tran3_nowaves=(['y'], time_int_depth_int_ssflux_mud02_tran3_nowaves),
        time_int_depth_int_ssflux_sand01_tran3_nowaves=(['y'], time_int_depth_int_ssflux_sand01_tran3_nowaves),
        time_int_depth_int_ssflux_sand02_tran3_nowaves=(['y'], time_int_depth_int_ssflux_sand02_tran3_nowaves),
        time_int_depth_int_ssflux_sand03_tran3_nowaves=(['y'], time_int_depth_int_ssflux_sand03_tran3_nowaves),
        time_int_depth_int_ssflux_mud01_tran3_double_waves=(['y'], time_int_depth_int_ssflux_mud01_tran3_double_waves),
        time_int_depth_int_ssflux_mud02_tran3_double_waves=(['y'], time_int_depth_int_ssflux_mud02_tran3_double_waves),
        time_int_depth_int_ssflux_sand01_tran3_double_waves=(['y'], time_int_depth_int_ssflux_sand01_tran3_double_waves),
        time_int_depth_int_ssflux_sand02_tran3_double_waves=(['y'], time_int_depth_int_ssflux_sand02_tran3_double_waves),
        time_int_depth_int_ssflux_sand03_tran3_double_waves=(['y'], time_int_depth_int_ssflux_sand03_tran3_double_waves)
        ),
    coords=dict(
        y_full=('y', y_rho_flat)
        ),
    attrs=dict(description='Time-integrated, depth-integrated suspended sediment fluxes at three different transects for the no waves, standard, and double waves model runs for each sediment class')) 
# Add more metadata?
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran1_std.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the standard run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran2_std.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the standard run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran3_std.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the standard run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran1_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the no waves run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran2_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the no waves run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran3_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the no waves run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran1_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the double waves run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran2_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the double waves run for all sediment classes combined'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_tran3_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the double waves run for all sediment classes combined'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran1_std.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the standard run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran1_std.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the standard run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran1_std.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the standard run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran1_std.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the standard run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran1_std.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the standard run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran2_std.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the standard run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran2_std.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the standard run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran2_std.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the standard run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran2_std.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the standard run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran2_std.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the standard run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran3_std.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the standard run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran3_std.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the standard run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran3_std.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the standard run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran3_std.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the standard run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran3_std.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the standard run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran1_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the no waves run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran1_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the no waves run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran1_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the no waves run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran1_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the no waves run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran1_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the no waves run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran2_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the no waves run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran2_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the no waves run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran2_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the no waves run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran2_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the no waves run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran2_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the no waves run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran3_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the no waves run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran3_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the no waves run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran3_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the no waves run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran3_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the no waves run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran3_nowaves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the no waves run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran1_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the double waves run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran1_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the double waves run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran1_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the double waves run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran1_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the double waves run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran1_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 1 for the double waves run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran2_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the double waves run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran2_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the double waves run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran2_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the double waves run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran2_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the double waves run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran2_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 2 for the double waves run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud01_tran3_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the double waves run for mud01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_mud02_tran3_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the double waves run for mud02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand01_tran3_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the double waves run for sand01'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand02_tran3_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the double waves run for sand02'
roms_time_int_depth_int_ssflux_transects_by_class.time_int_depth_int_ssflux_sand03_tran3_double_waves.name='time-integrated, depth-integrated suspended sediment flux at transect 3 for the double waves run for sand03'

roms_time_int_depth_int_ssflux_transects_by_class.y_full.name='latitudinal distance (meter)'

# Save to a netcdf
#roms_time_int_depth_int_ssflux_transects_by_class.to_netcdf('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/fig12_roms_time_int_depth_int_ssflux_transects_by_class.nc')







