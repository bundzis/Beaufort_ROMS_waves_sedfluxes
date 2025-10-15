#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 14:35:55 2023

@author: brun1463
"""

####################### Violin Plots for ROMS & CODA ###########################
# The purpose of this script is to make violin plots comparing the eastward
# an northward currents from CODA with those from ROMS output. This will also
# be done to compare CODA waves with ERA5 waves used to force ROMS. 
# 
# Notes:
# - This script needs to be run in the xroms env 
################################################################################


# Load in the packages
import xarray as xr
import xroms
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.dates import DateFormatter
import pandas as pd
import numpy as np
import cmocean 
import scipy.io
import seaborn as sns
from datetime import datetime #as dt
from datetime import timedelta
from dateutil.relativedelta import relativedelta
from scipy.stats import pearsonr
from glob import glob

# Set a universal fontsize
fontsize = 20

# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize) 
matplotlib.rc('ytick', labelsize=fontsize)

# Prevent tick labels from overlapping
matplotlib.rcParams['xtick.major.pad'] = 10
matplotlib.rcParams['ytick.major.pad'] = 10


# Load in the CODA data
mooring_data = scipy.io.loadmat('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/CODA_mooring_data/CODAmoorings_masterfile (1).mat')
mooring_data

# ------------------- Prep CODA data -------------------------
# Use a function from the internet to convert matlab time
# (thank you Example #3: https://www.programcreek.com/python/example/100062/datetime.datetime.fromordinal)
def _from_ordinal(x, tz=None):
    ix = int(x)
    dt = datetime.fromordinal(ix)
    remainder = float(x) - ix
    hour, remainder = divmod(24 * remainder, 1)
    minute, remainder = divmod(60 * remainder, 1)
    second, remainder = divmod(60 * remainder, 1)
    microsecond = int(1e6 * remainder)
    if microsecond < 10:
        microsecond = 0  # compensate for rounding errors
    dt = datetime(dt.year, dt.month, dt.day, int(hour), int(minute),
                  int(second), microsecond)
    if tz is not None:
        dt = dt.astimezone(tz)

    if microsecond > 999990:  # compensate for rounding errors
        dt += timedelta(microseconds=1e6 - microsecond)

    return dt


# --- Set up time for both moorings ---
# First convert and slice time for CODA moorings 
# Pull out time
# CODA mooring 2
moor_time1 = mooring_data['S2A1average']['time'][0]
moor_time1_burst = mooring_data['S2A1burst']['time'][0]
# CODA mooring 3
moor_time2 = mooring_data['S3A1average']['time'][0]
moor_time2_burst = mooring_data['S3A1burst']['time'][0]

# Loop through each time step and convert it to datetime
# Make empty lists to hold the dates
moor1_datetime_tmp = []
moor1_datetime_final = []
moor1_datetime_tmp_burst = []
moor1_datetime_final_burst = []
moor2_datetime_tmp = []
moor2_datetime_final = []
moor2_datetime_tmp_burst = []
moor2_datetime_final_burst = []
# Loop through time1 average 
for i in range(len(moor_time1)):
    # Mooring 2
    a = _from_ordinal(moor_time1[i])
    b = a - relativedelta(years=1)
    moor1_datetime_final.append(b)
    
# Loop through time1 burst
for i in range(len(moor_time1_burst)):
    # Mooring 2
    a1 = _from_ordinal(moor_time1_burst[i])
    b1 = a1 - relativedelta(years=1)
    moor1_datetime_final_burst.append(b1)
    
# Loop through time2 average
for i in range(len(moor_time2)):
    # Mooring 3
    a2 = _from_ordinal(moor_time2[i])
    b2 = a2 - relativedelta(years=1)
    moor2_datetime_final.append(b2)
    
# Loop through time2 burst
for i in range(len(moor_time2_burst)):
    # Mooring 3
    a3 = _from_ordinal(moor_time2_burst[i])
    b3 = a3 - relativedelta(years=1)
    moor2_datetime_final_burst.append(b3)
    

# The year 2020 was a leap year and this does not account for this
# so manually subtract a day from all dates with year 2020
# Make a copy to work with 
moor1_datetime_final_leap = moor1_datetime_final.copy()
moor1_datetime_final_leap_burst = moor1_datetime_final_burst.copy()
moor2_datetime_final_leap = moor2_datetime_final.copy()
moor2_datetime_final_leap_burst = moor2_datetime_final_burst.copy()
# Loop through time1 average 
for j in range(len(moor1_datetime_final_leap)):
    # Check the year
    if (moor1_datetime_final_leap[j].year == 2020):
        # Check the month
        if(moor1_datetime_final_leap[j].month > 2):
           # print('pre leap: ', moor_datetime_final_leap[j])
            moor1_datetime_final_leap[j] = moor1_datetime_final_leap[j] - timedelta(days=1)
            #print('post leap: ', moor_datetime_final_leap[j])

# Loop through time1 burst 
for j in range(len(moor1_datetime_final_leap_burst)):
    # Check the year
    if (moor1_datetime_final_leap_burst[j].year == 2020):
        # Check the month
        if(moor1_datetime_final_leap_burst[j].month > 2):
           # print('pre leap: ', moor_datetime_final_leap[j])
            moor1_datetime_final_leap_burst[j] = moor1_datetime_final_leap_burst[j] - timedelta(days=1)
            #print('post leap: ', moor_datetime_final_leap[j])

# Loop through time2 average 
for j in range(len(moor2_datetime_final_leap)):
    # Check the year
    if (moor2_datetime_final_leap[j].year == 2020):
        # Check the month
        if(moor2_datetime_final_leap[j].month > 2):
           # print('pre leap: ', moor_datetime_final_leap[j])
            moor2_datetime_final_leap[j] = moor2_datetime_final_leap[j] - timedelta(days=1)
            #print('post leap: ', moor_datetime_final_leap[j])
            
# Loop through time2 burst
for j in range(len(moor2_datetime_final_leap_burst)):
    # Check the year
    if (moor2_datetime_final_leap_burst[j].year == 2020):
        # Check the month
        if(moor2_datetime_final_leap_burst[j].month > 2):
           # print('pre leap: ', moor_datetime_final_leap[j])
            moor2_datetime_final_leap_burst[j] = moor2_datetime_final_leap_burst[j] - timedelta(days=1)
            #print('post leap: ', moor_datetime_final_leap[j])

# Check to see if this worked (it does!)
#print(moor_datetime_final_leap[0])
#print(moor_datetime_final_leap[-1])
#print(moor_datetime_final_leap[2589])
#print(moor_datetime_final_leap[2588])
#input('press entre to continue...')

# Pull out matching time for moorings
# CODA Mooring 2 (July 1 - September 23)
# Pull out July 1 - September 24, 2020 (until end) in coda data 
moor1_time_julsep_2020 = moor1_datetime_final_leap[33240:]

# CODA Mooring 3 (July 1 - August 1)
moor2_time_julaug_2020 = moor2_datetime_final_leap[32993:]


# --- Set up latitude and longitude for both moorings ---
# CODA mooring 2
# Latitude 
# Check min and max
print('max moor lat: ', np.max(mooring_data['S2A1average']['lat']))
print('min moor lat: ', np.min(mooring_data['S2A1average']['lat']))
# Since these are the same, pull out the first value
moor_lat1 = mooring_data['S2A1average']['lat'][0,0]
print('moor lat1: ', moor_lat1)

# Longitude
print('max moor lon: ', np.max(mooring_data['S2A1average']['lon']))
print('min moor lon: ', np.min(mooring_data['S2A1average']['lon']))
# Since these are the same, pull out the first value
moor_lon1 = mooring_data['S2A1average']['lon'][0,0]
print('moor lon1: ', moor_lon1)

# Pull out the indices of the ROMS grid that match the mooring
# (found by hand)
eta_rho_moor1 = 103
xi_rho_moor1 = 202

# CODA mooring 3
# Check min and max
print('max moor lat: ', np.max(mooring_data['S3A1average']['lat']))
print('min moor lat: ', np.min(mooring_data['S3A1average']['lat']))
# Since these are the same, pull out the first value
moor_lat2 = mooring_data['S3A1average']['lat'][0,0]
print('moor lat2: ', moor_lat2)

# Longitude
print('max moor lon: ', np.max(mooring_data['S3A1average']['lon']))
print('min moor lon: ', np.min(mooring_data['S3A1average']['lon']))
# Since these are the same, pull out the first value
moor_lon2 = mooring_data['S3A1average']['lon'][0,0]
print('moor lon2: ', moor_lon2)

# Pull out the indices of the ROMS grid that match the mooring
# (found by hand)
eta_rho_moor2 = 73
xi_rho_moor2 = 390


# Try making an array and filling it with the values in a shape I want
# Use the full data first
# --- CODA 2 ---
coda2_east = np.empty((len(mooring_data['S2A1average']['time'][0]), len(mooring_data['S2A1average']['z'][0,0][0])))
coda2_east_inv = np.empty((len(mooring_data['S2A1average']['z'][0,0][0]), len(mooring_data['S2A1average']['time'][0])))
coda2_north_inv = np.empty((len(mooring_data['S2A1average']['z'][0,0][0]), len(mooring_data['S2A1average']['time'][0])))

# Loop through time
for i in range(len(mooring_data['S2A1average']['time'][0])):
    # Loop through depth
    for j in range(len(mooring_data['S2A1average']['z'][0,0][0])):
        # Fill the arrays
        coda2_east[i,j] = mooring_data['S2A1average']['east'][0,i][0,j]
        
# Loop through depth
for k in range(len(mooring_data['S2A1average']['z'][0,0][0])):
    # Loop through time
    for l in range(len(mooring_data['S2A1average']['time'][0])):
        # Fill the arrays
        coda2_east_inv[k,l] = mooring_data['S2A1average']['east'][0,l][0,k]
        coda2_north_inv[k,l] = mooring_data['S2A1average']['north'][0,l][0,k]
        
# --- CODA 3 ---
coda3_east = np.empty((len(mooring_data['S3A1average']['time'][0]), len(mooring_data['S3A1average']['z'][0,0][0])))
coda3_east_inv = np.empty((len(mooring_data['S3A1average']['z'][0,0][0]), len(mooring_data['S3A1average']['time'][0])))
coda3_north_inv = np.empty((len(mooring_data['S3A1average']['z'][0,0][0]), len(mooring_data['S3A1average']['time'][0])))

# Loop through time
for i in range(len(mooring_data['S3A1average']['time'][0])):
    # Loop through depth
    for j in range(len(mooring_data['S3A1average']['z'][0,0][0])):
        # Fill the arrays
        coda3_east[i,j] = mooring_data['S3A1average']['east'][0,i][0,j]
        
# Loop through depth
for k in range(len(mooring_data['S3A1average']['z'][0,0][0])):
    # Loop through time
    for l in range(len(mooring_data['S3A1average']['time'][0])):
        # Fill the arrays
        coda3_east_inv[k,l] = mooring_data['S3A1average']['east'][0,l][0,k]
        coda3_north_inv[k,l] = mooring_data['S3A1average']['north'][0,l][0,k]


# ------------------- Prep ROMS data -------------------------
# WORKS FOR 1 FILE 
# =============================================================================
# # Load in the ROMS data
# #ds = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/model_output/Full_run_0003_sponge_swell/ocean_his_biggrid010_gridwindsiniwaves_rivs_si_smooth006_nobulk_chaflaradnudclm_dbsed0007_0001.nc')
# ds = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/model_output/Full_run_0001/ocean_his_biggrid010_gridwindsiniwaves_rivs_si_smooth006_nobulk_chaflaradnudclm_dbsed0004_0001.nc')
# ds['h'] = ds.bath
# ds, xgrid = xroms.roms_dataset(ds)
# 
# # Interpolate the ROMS data to the CODA depths to make things easier to compare
# # Now format these to plot the same things as is in the CODA plots of currents
# # with depth and interpolate onto the CODA depths
# # Interpolate u and v onto the CODA depths at each mooring
# # Set the coda depths 
# coda2_depths = [-3, -5, -7, -9, -11, -13, -15, -17, -19, -21]
# coda2_depths = np.asarray(coda2_depths)
# coda3_depths = [-3.25, -5.25, -7.25, -9.25, -11.25, -13.25, -15.25, -17.25, -19.25, -21.25,
#                 -23.25, -25.25, -27.25, -29.25, -31.25, -33.25]
# coda3_depths = np.asarray(coda3_depths)
# 
# # Set the CODA indices
# # CODA 2
# eta_rho_coda2 = 103
# xi_rho_coda2 = 202
# # CODA 3
# eta_rho_coda3 = 73
# xi_rho_coda3 = 390
# 
# 
# # Now interpolate u and v
# # (easier to do this for whole array, then slice to CODA areas)
# # U 
# u_roms_og = ds.u
# # V
# v_roms_og = ds.v
# 
# # Interpolate all of these
# u_roms_interp_coda2 = xroms.isoslice(u_roms_og, coda2_depths, xgrid)
# v_roms_interp_coda2 = xroms.isoslice(v_roms_og, coda2_depths, xgrid)
# u_roms_interp_coda3 = xroms.isoslice(u_roms_og, coda3_depths, xgrid)
# v_roms_interp_coda3 = xroms.isoslice(v_roms_og, coda3_depths, xgrid)
# =============================================================================

# NEW ATTEMPT WITH LOOP?
# MAke a function to process the output
def load_and_interp_roms_output(filename, coda_depths):
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
    
    # Pull out the ROMS currents  
    # U Currents
    u_roms_og = ds.u
    # V Currents
    v_roms_og = ds.v
    
    # Interopolate onto the given CODA depths
    u_roms_interp_coda = xroms.isoslice(u_roms_og, coda_depths, xgrid)
    v_roms_interp_coda = xroms.isoslice(v_roms_og, coda_depths, xgrid)
    
    # Return these currents 
    return(u_roms_interp_coda, v_roms_interp_coda)
    

# Loop through and pull out all the file names 
# First, get all the output file names
# dbsed0003 (2020)
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
#num_files = len(file_names2)
# We only want to use the first 5 ocean_his files so
# manually set this to 5
num_files = 5

# Set the CODA depths 
coda2_depths = [-3, -5, -7, -9, -11, -13, -15, -17, -19, -21]
coda2_depths = np.asarray(coda2_depths)
coda3_depths = [-3.25, -5.25, -7.25, -9.25, -11.25, -13.25, -15.25, -17.25, -19.25, -21.25,
                -23.25, -25.25, -27.25, -29.25, -31.25, -33.25]
coda3_depths = np.asarray(coda3_depths)

# Set the CODA indices
# CODA 2
eta_rho_coda2 = 103
xi_rho_coda2 = 202
# CODA 3
eta_rho_coda3 = 73
xi_rho_coda3 = 390

# Make an xarray dataset to merge the timeseries into for the different points
# CODA 2
u_roms_interp_coda2_all_output = xr.Dataset()
v_roms_interp_coda2_all_output = xr.Dataset()
# CODA 3
u_roms_interp_coda3_all_output = xr.Dataset()
v_roms_interp_coda3_all_output = xr.Dataset()


# Now loop through each file and call the function to interpolate
# the u and v currents at both moorings
for i in range(num_files):
    # Call the function for CODA 2 
    u_roms_interp_coda2_output1, v_roms_interp_coda2_output1 = load_and_interp_roms_output(file_names2[i], coda2_depths)
    # Call the function for CODA 3
    u_roms_interp_coda3_output1, v_roms_interp_coda3_output1 = load_and_interp_roms_output(file_names2[i], coda3_depths)
    
    
    
    # Add this onto the datasets for each mooring and current 
    # CODA 2 
    u_roms_interp_coda2_all_output = xr.merge([u_roms_interp_coda2_all_output, u_roms_interp_coda2_output1])
    v_roms_interp_coda2_all_output = xr.merge([v_roms_interp_coda2_all_output, v_roms_interp_coda2_output1])
    # CODA 3
    u_roms_interp_coda3_all_output = xr.merge([u_roms_interp_coda3_all_output, u_roms_interp_coda3_output1])
    v_roms_interp_coda3_all_output = xr.merge([v_roms_interp_coda3_all_output, v_roms_interp_coda3_output1])

# Print shapes to see if this worked 
print('u_roms_interp_coda2_all_output: ', u_roms_interp_coda2_all_output)
print('\nv_roms_interp_coda2_all_output: ', v_roms_interp_coda2_all_output)
print('\nu_roms_interp_coda3_all_output: ', u_roms_interp_coda3_all_output)
print('\nv_roms_interp_coda3_all_output: ', v_roms_interp_coda3_all_output)



# Slice to the time period we want 
# (JK use the full file for now since we are just loading in the first ocean_his
# for now)

# ------------------------ Violin Plots for Currents ----------------------------
# Make a dataframe with everything in it 
# CODA2 _east_overlap, ROMS_at_CODA2_east_overlap, CODA2_west_overlap, ROMS_at_CODA2_west_overlap,
# CODA3_east_overlap, ROMS_at_CODA3_east_overlap, CODA3_west_overlap, ROMS_at_CODA3_west_overlap
#df_all =
#


# First start with CODA data
# Combine all of the CODA data into a list
# **Confirm the indexing for the depths of interpolated ROMS data**
vio_data_surf_coda = [coda2_east_inv[0,33240:], coda2_north_inv[0,33240:], coda3_east_inv[0,32993:37442], coda3_north_inv[0,32993:37442]]
# Surface 
vio_data_surf_all = [coda2_east_inv[0,33240:], u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2], coda2_north_inv[0,33240:], 
                     v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2], coda3_east_inv[0,32993:37442], 
                     u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3], coda3_north_inv[0,32993:37442],
                     v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3]]
# Bottom (-1 index gives all nan)
vio_data_bot_all = [coda2_east_inv[-2,33240:], u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2], coda2_north_inv[-2,33240:], 
                     v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2], coda3_east_inv[-3,32993:37442], 
                     u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3], coda3_north_inv[-3,32993:37442],
                     v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3]]

# Surface Currents
# Make the figure 
fig1, ax1 = plt.subplots(figsize=(21,5))
# Set the palette
#palette=['r', 'lightcoral', 'blue', 'dodgerblue', 'orangered', 'lightsalmon', 'royalblue', 'lightskyblue']
palette=['r', 'r', 'dodgerblue', 'dodgerblue', 'lightcoral', 'lightcoral', 'lightskyblue', 'lightskyblue']
s1 = sns.violinplot(data=vio_data_surf_all, ax=ax1, palette=palette)
# CODA
#ax1.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
ax1.set_xticklabels(['CODA2 \nEastward', 'ROMS at CODA2 \nEastward', 'CODA2 \nNorthward', 'ROMS at CODA2 \nNorthward','CODA3 \nEastward', 
                     'ROMS at CODA3 \nEastward', 'CODA3 \nNorthward', 'ROMS at CODA3 \nNorthward'])
ax1.set_title('CODA & ROMS Surface Currents', fontsize=fontsize)
ax1.set_ylabel('Current Speed (m/s)', fontsize=fontsize)


# Bottom Currents
# Make the figure 
fig2, ax2 = plt.subplots(figsize=(21,5))
s2 = sns.violinplot(data=vio_data_bot_all, ax=ax2, palette=palette)
# CODA
#ax1.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
ax2.set_xticklabels(['CODA2 \nEastward', 'ROMS at CODA2 \nEastward', 'CODA2 \nNorthward', 'ROMS at CODA2 \nNorthward','CODA3 \nEastward', 
                     'ROMS at CODA3 \nEastward', 'CODA3 \nNorthward', 'ROMS at CODA3 \nNorthward'])
ax2.set_title('CODA & ROMS Bottom Currents', fontsize=fontsize)
ax2.set_ylabel('Current Speed (m/s)', fontsize=fontsize)



# Depth-Averaged 
# Average everything over depth
coda2_east_inv_depthavg = np.nanmean(coda2_east_inv, axis=0)
coda2_north_inv_depthavg = np.nanmean(coda2_north_inv, axis=0)
coda3_east_inv_depthavg = np.nanmean(coda3_east_inv, axis=0)
coda3_north_inv_depthavg = np.nanmean(coda3_north_inv, axis=0)
u_roms_interp_coda2_all_output_depthavg = u_roms_interp_coda2_all_output.mean(dim='z_rho_u')
v_roms_interp_coda2_all_output_depthavg = v_roms_interp_coda2_all_output.mean(dim='z_rho_v')
u_roms_interp_coda3_all_output_depthavg = u_roms_interp_coda3_all_output.mean(dim='z_rho_u')
v_roms_interp_coda3_all_output_depthavg = v_roms_interp_coda3_all_output.mean(dim='z_rho_v')
vio_data_depthavg_all = [coda2_east_inv_depthavg[33240:], u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2], 
                         coda2_north_inv_depthavg[33240:], v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2], 
                         coda3_east_inv_depthavg[32993:37442], u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3], 
                         coda3_north_inv_depthavg[32993:37442], v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3]]

# Make the figure (m/s)
fig3, ax3 = plt.subplots(figsize=(21,5))
s3 = sns.violinplot(data=vio_data_depthavg_all, ax=ax3, palette=palette)
# CODA
#ax3.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
ax3.set_xticklabels(['CODA2 \nEastward', 'ROMS at CODA2 \nEastward', 'CODA2 \nNorthward', 'ROMS at CODA2 \nNorthward','CODA3 \nEastward', 
                     'ROMS at CODA3 \nEastward', 'CODA3 \nNorthward', 'ROMS at CODA3 \nNorthward'])
ax3.set_title('CODA & ROMS Depth-Averaged Currents', fontsize=fontsize)
ax3.set_ylabel('Current Speed (m/s)', fontsize=fontsize)



# =============================================================================
# # Panel version of all the plots above (m/s)
# fontsize4 = 18
# # Set the tick size for all plots
# matplotlib.rc('xtick', labelsize=fontsize4) 
# matplotlib.rc('ytick', labelsize=fontsize4)
# # Make the figure 
# fig4, ax4 = plt.subplots(3, figsize=(15,9)) # (20,10)
# s4 = sns.violinplot(data=vio_data_surf_all, ax=ax4[0], palette=palette, width=1.15)
# s5 = sns.violinplot(data=vio_data_bot_all, ax=ax4[1], palette=palette, width=1.15)
# s6 = sns.violinplot(data=vio_data_depthavg_all, ax=ax4[2], palette=palette, width=0.99)
# ax4[0].axhline(y=0.0, color='k', linestyle='--', linewidth=1.40)
# ax4[1].axhline(y=0.0, color='k', linestyle='--', linewidth=1.40)
# ax4[2].axhline(y=0.0, color='k', linestyle='--', linewidth=1.40)
# # CODA
# #ax3.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# # All
# #ax4[2].set_xticklabels(['CODA2 \nEastward', 'ROMS at CODA2 \nEastward', 'CODA2 \nNorthward', 'ROMS at CODA2 \nNorthward','CODA3 \nEastward', 
#  #                    'ROMS at CODA3 \nEastward', 'CODA3 \nNorthward', 'ROMS at CODA3 \nNorthward'])
# ax4[2].set_xticklabels(['Observed \nCODA2 \nEastward', 'Modeled \nCODA2 \nEastward', 
#                      'Observed \nCODA2 \nNorthward', 'Modeled \nCODA2 \nNorthward',
#                      'Observed \nCODA3 \nEastward', 'Modeled \nCODA3 \nEastward', 
#                      'Observed \nCODA3 \nNorthward', 'Modeled \nCODA3 \nNorthward'])
# #ax4.set_title('CODA & ROMS Depth-Averaged Currents', fontsize=fontsize)
# plt.setp(ax4[0].get_xticklabels(), visible=False)
# plt.setp(ax4[1].get_xticklabels(), visible=False)
# ax4[0].set_ylabel('Surface \nCurrent \nSpeed \n(m/s)', fontsize=fontsize4, rotation=0,
#                   labelpad=40, va='center')
# ax4[1].set_ylabel('Bottom \nCurrent \nSpeed \n(m/s)', fontsize=fontsize4, rotation=0,
#                   labelpad=40, va='center')
# ax4[2].set_ylabel('Depth-\nAveraged \nCurrent \nSpeed \n(m/s)', fontsize=fontsize4, rotation=0,
#                   labelpad=40, va='center')
# 
# # Edit colors and elements of whisker plot inside violin
# # Top Row
# # First violin plot: whisker, box, dot
# s4.get_children()[1].set_color('yellow')
# s4.get_children()[2].set_color('k')
# s4.get_children()[3].set_color('white')
# # Second violin plot: whisker, box, dot
# s4.get_children()[5].set_color('yellow') 
# s4.get_children()[6].set_color('k')
# s4.get_children()[7].set_color('white')
# # Third violin plot: whisker, box, dot
# s4.get_children()[9].set_color('yellow') 
# s4.get_children()[10].set_color('k')
# s4.get_children()[11].set_color('white')
# # Fourth violin plot: whisker, box, dot
# s4.get_children()[13].set_color('yellow') 
# s4.get_children()[14].set_color('k')
# s4.get_children()[15].set_color('white')
# # Fifth violin plot: whisker, box, dot
# s4.get_children()[17].set_color('yellow') 
# s4.get_children()[18].set_color('k')
# s4.get_children()[19].set_color('white')
# # Sixth violin plot: whisker, box, dot
# s4.get_children()[21].set_color('yellow') 
# s4.get_children()[22].set_color('k')
# s4.get_children()[23].set_color('white')
# # Seventh violin plot: whisker, box, dot
# s4.get_children()[25].set_color('yellow') 
# s4.get_children()[26].set_color('k')
# s4.get_children()[27].set_color('white')
# # Eighth violin plot: whisker, box, dot
# s4.get_children()[29].set_color('yellow') 
# s4.get_children()[30].set_color('k')
# s4.get_children()[31].set_color('white')
# # Middle Row
# # First violin plot: whisker, box, dot
# s5.get_children()[1].set_color('yellow')
# s5.get_children()[2].set_color('k')
# s5.get_children()[3].set_color('white')
# # Second violin plot: whisker, box, dot
# s5.get_children()[5].set_color('yellow') 
# s5.get_children()[6].set_color('k')
# s5.get_children()[7].set_color('white')
# # Third violin plot: whisker, box, dot
# s5.get_children()[9].set_color('yellow') 
# s5.get_children()[10].set_color('k')
# s5.get_children()[11].set_color('white')
# # Fourth violin plot: whisker, box, dot
# s5.get_children()[13].set_color('yellow') 
# s5.get_children()[14].set_color('k')
# s5.get_children()[15].set_color('white')
# # Fifth violin plot: whisker, box, dot
# s5.get_children()[17].set_color('yellow') 
# s5.get_children()[18].set_color('k')
# s5.get_children()[19].set_color('white')
# # Sixth violin plot: whisker, box, dot
# s5.get_children()[21].set_color('yellow') 
# s5.get_children()[22].set_color('k')
# s5.get_children()[23].set_color('white')
# # Seventh violin plot: whisker, box, dot
# s5.get_children()[25].set_color('yellow') 
# s5.get_children()[26].set_color('k')
# s5.get_children()[27].set_color('white')
# # Eighth violin plot: whisker, box, dot
# s5.get_children()[29].set_color('yellow') 
# s5.get_children()[30].set_color('k')
# s5.get_children()[31].set_color('white')
# # Bottom Row
# # First violin plot: whisker, box, dot
# s6.get_children()[1].set_color('yellow')
# s6.get_children()[2].set_color('k')
# s6.get_children()[3].set_color('white')
# # Second violin plot: whisker, box, dot
# s6.get_children()[5].set_color('yellow') 
# s6.get_children()[6].set_color('k')
# s6.get_children()[7].set_color('white')
# # Third violin plot: whisker, box, dot
# s6.get_children()[9].set_color('yellow') 
# s6.get_children()[10].set_color('k')
# s6.get_children()[11].set_color('white')
# # Fourth violin plot: whisker, box, dot
# s6.get_children()[13].set_color('yellow') 
# s6.get_children()[14].set_color('k')
# s6.get_children()[15].set_color('white')
# # Fifth violin plot: whisker, box, dot
# s6.get_children()[17].set_color('yellow') 
# s6.get_children()[18].set_color('k')
# s6.get_children()[19].set_color('white')
# # Sixth violin plot: whisker, box, dot
# s6.get_children()[21].set_color('yellow') 
# s6.get_children()[22].set_color('k')
# s6.get_children()[23].set_color('white')
# # Seventh violin plot: whisker, box, dot
# s6.get_children()[25].set_color('yellow') 
# s6.get_children()[26].set_color('k')
# s6.get_children()[27].set_color('white')
# # Eighth violin plot: whisker, box, dot
# s6.get_children()[29].set_color('yellow') 
# s6.get_children()[30].set_color('k')
# s6.get_children()[31].set_color('white')
# 
# plt.subplots_adjust(hspace=0.1)
# 
# # Label the subplots 
# plt.text(0.879, 0.845, 'a)', fontsize=fontsize4, fontweight='bold', transform=plt.gcf().transFigure)
# plt.text(0.879, 0.579, 'b)', fontsize=fontsize4, fontweight='bold', transform=plt.gcf().transFigure)
# plt.text(0.879, 0.321, 'c)', fontsize=fontsize4, fontweight='bold', transform=plt.gcf().transFigure)
# =============================================================================


# Panel version of all the plots above (cm/s)
# Save data as cm/s
# Surface 
vio_data_surf_all_cms = [coda2_east_inv[0,33240:]*100, u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2]*100, coda2_north_inv[0,33240:]*100, 
                     v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2]*100, coda3_east_inv[0,32993:37442]*100, 
                     u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3]*100, coda3_north_inv[0,32993:37442]*100,
                     v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3]*100]
# Bottom (-1 index gives all nan)
vio_data_bot_all_cms = [coda2_east_inv[-2,33240:]*100, u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2]*100, coda2_north_inv[-2,33240:]*100, 
                     v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2]*100, coda3_east_inv[-3,32993:37442]*100, 
                     u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3]*100, coda3_north_inv[-3,32993:37442]*100,
                     v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3]*100]
# Depth-averaged
vio_data_depthavg_all_cms = [coda2_east_inv_depthavg[33240:]*100, u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2]*100, 
                         coda2_north_inv_depthavg[33240:]*100, v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2]*100, 
                         coda3_east_inv_depthavg[32993:37442]*100, u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3]*100, 
                         coda3_north_inv_depthavg[32993:37442]*100, v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3]*100]


fontsize4 = 18
# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize4) 
matplotlib.rc('ytick', labelsize=fontsize4)
# Make the figure 
fig4, ax4 = plt.subplots(3, figsize=(15,9)) # (20,10)
# Depth-averaged
s6 = sns.violinplot(data=vio_data_depthavg_all_cms, ax=ax4[0], palette=palette, width=0.99)
# Surface 
s4 = sns.violinplot(data=vio_data_surf_all_cms, ax=ax4[1], palette=palette, width=1.15)
# Bottom
s5 = sns.violinplot(data=vio_data_bot_all_cms, ax=ax4[2], palette=palette, width=1.15)
ax4[0].axhline(y=0.0, color='k', linestyle='--', linewidth=1.40)
ax4[1].axhline(y=0.0, color='k', linestyle='--', linewidth=1.40)
ax4[2].axhline(y=0.0, color='k', linestyle='--', linewidth=1.40)
# CODA
#ax3.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
#ax4[2].set_xticklabels(['CODA2 \nEastward', 'ROMS at CODA2 \nEastward', 'CODA2 \nNorthward', 'ROMS at CODA2 \nNorthward','CODA3 \nEastward', 
 #                    'ROMS at CODA3 \nEastward', 'CODA3 \nNorthward', 'ROMS at CODA3 \nNorthward'])
ax4[2].set_xticklabels(['Observed \nCODA2 \nEastward', 'Modeled \nCODA2 \nEastward', 
                     'Observed \nCODA2 \nNorthward', 'Modeled \nCODA2 \nNorthward',
                     'Observed \nCODA3 \nEastward', 'Modeled \nCODA3 \nEastward', 
                     'Observed \nCODA3 \nNorthward', 'Modeled \nCODA3 \nNorthward'])
#ax4.set_title('CODA & ROMS Depth-Averaged Currents', fontsize=fontsize)
plt.setp(ax4[0].get_xticklabels(), visible=False)
plt.setp(ax4[1].get_xticklabels(), visible=False)
ax4[1].set_ylabel('Surface \nCurrent \nSpeed \n(cm/s)', fontsize=fontsize4, rotation=0,
                  labelpad=40, va='center')
ax4[2].set_ylabel('Bottom \nCurrent \nSpeed \n(cm/s)', fontsize=fontsize4, rotation=0,
                  labelpad=40, va='center')
ax4[0].set_ylabel('Depth-\nAveraged \nCurrent \nSpeed \n(cm/s)', fontsize=fontsize4, rotation=0,
                  labelpad=40, va='center')

# Edit colors and elements of whisker plot inside violin
# Top Row
# First violin plot: whisker, box, dot
s4.get_children()[1].set_color('yellow')
s4.get_children()[2].set_color('k')
s4.get_children()[3].set_color('white')
# Second violin plot: whisker, box, dot
s4.get_children()[5].set_color('yellow') 
s4.get_children()[6].set_color('k')
s4.get_children()[7].set_color('white')
# Third violin plot: whisker, box, dot
s4.get_children()[9].set_color('yellow') 
s4.get_children()[10].set_color('k')
s4.get_children()[11].set_color('white')
# Fourth violin plot: whisker, box, dot
s4.get_children()[13].set_color('yellow') 
s4.get_children()[14].set_color('k')
s4.get_children()[15].set_color('white')
# Fifth violin plot: whisker, box, dot
s4.get_children()[17].set_color('yellow') 
s4.get_children()[18].set_color('k')
s4.get_children()[19].set_color('white')
# Sixth violin plot: whisker, box, dot
s4.get_children()[21].set_color('yellow') 
s4.get_children()[22].set_color('k')
s4.get_children()[23].set_color('white')
# Seventh violin plot: whisker, box, dot
s4.get_children()[25].set_color('yellow') 
s4.get_children()[26].set_color('k')
s4.get_children()[27].set_color('white')
# Eighth violin plot: whisker, box, dot
s4.get_children()[29].set_color('yellow') 
s4.get_children()[30].set_color('k')
s4.get_children()[31].set_color('white')
# Middle Row
# First violin plot: whisker, box, dot
s5.get_children()[1].set_color('yellow')
s5.get_children()[2].set_color('k')
s5.get_children()[3].set_color('white')
# Second violin plot: whisker, box, dot
s5.get_children()[5].set_color('yellow') 
s5.get_children()[6].set_color('k')
s5.get_children()[7].set_color('white')
# Third violin plot: whisker, box, dot
s5.get_children()[9].set_color('yellow') 
s5.get_children()[10].set_color('k')
s5.get_children()[11].set_color('white')
# Fourth violin plot: whisker, box, dot
s5.get_children()[13].set_color('yellow') 
s5.get_children()[14].set_color('k')
s5.get_children()[15].set_color('white')
# Fifth violin plot: whisker, box, dot
s5.get_children()[17].set_color('yellow') 
s5.get_children()[18].set_color('k')
s5.get_children()[19].set_color('white')
# Sixth violin plot: whisker, box, dot
s5.get_children()[21].set_color('yellow') 
s5.get_children()[22].set_color('k')
s5.get_children()[23].set_color('white')
# Seventh violin plot: whisker, box, dot
s5.get_children()[25].set_color('yellow') 
s5.get_children()[26].set_color('k')
s5.get_children()[27].set_color('white')
# Eighth violin plot: whisker, box, dot
s5.get_children()[29].set_color('yellow') 
s5.get_children()[30].set_color('k')
s5.get_children()[31].set_color('white')
# Bottom Row
# First violin plot: whisker, box, dot
s6.get_children()[1].set_color('yellow')
s6.get_children()[2].set_color('k')
s6.get_children()[3].set_color('white')
# Second violin plot: whisker, box, dot
s6.get_children()[5].set_color('yellow') 
s6.get_children()[6].set_color('k')
s6.get_children()[7].set_color('white')
# Third violin plot: whisker, box, dot
s6.get_children()[9].set_color('yellow') 
s6.get_children()[10].set_color('k')
s6.get_children()[11].set_color('white')
# Fourth violin plot: whisker, box, dot
s6.get_children()[13].set_color('yellow') 
s6.get_children()[14].set_color('k')
s6.get_children()[15].set_color('white')
# Fifth violin plot: whisker, box, dot
s6.get_children()[17].set_color('yellow') 
s6.get_children()[18].set_color('k')
s6.get_children()[19].set_color('white')
# Sixth violin plot: whisker, box, dot
s6.get_children()[21].set_color('yellow') 
s6.get_children()[22].set_color('k')
s6.get_children()[23].set_color('white')
# Seventh violin plot: whisker, box, dot
s6.get_children()[25].set_color('yellow') 
s6.get_children()[26].set_color('k')
s6.get_children()[27].set_color('white')
# Eighth violin plot: whisker, box, dot
s6.get_children()[29].set_color('yellow') 
s6.get_children()[30].set_color('k')
s6.get_children()[31].set_color('white')

plt.subplots_adjust(hspace=0.1)

# Label the subplots 
plt.text(0.879, 0.845, 'a)', fontsize=fontsize4, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.879, 0.579, 'b)', fontsize=fontsize4, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.879, 0.321, 'c)', fontsize=fontsize4, fontweight='bold', transform=plt.gcf().transFigure)




# Print a bunch of stats from the plot above (m/s)
# =============================================================================
# # Top row (surface)
# # Medians
# # East
# median_surf_cur_coda2_east = np.nanmedian(coda2_east_inv[0,33240:])
# print('Median surf cur coda2 east: ', median_surf_cur_coda2_east)
# median_surf_cur_coda2_east_roms = np.median(u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2])
# print('Median surf cur coda2 east roms: ', median_surf_cur_coda2_east_roms)
# median_surf_cur_coda3_east = np.nanmedian(coda3_east_inv[0,32993:37442])
# print('Median surf cur coda3 east: ', median_surf_cur_coda3_east)
# median_surf_cur_coda3_east_roms = np.median(u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3])
# print('Median ssurf cur coda3 east roms: ', median_surf_cur_coda3_east_roms)
# # North
# median_surf_cur_coda2_north = np.nanmedian(coda2_north_inv[0,33240:])
# print('Median surf cur coda2 north: ', median_surf_cur_coda2_north)
# median_surf_cur_coda2_north_roms = np.median(v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2])
# print('Median surf cur coda2 north roms: ', median_surf_cur_coda2_north_roms)
# median_surf_cur_coda3_north = np.nanmedian(coda3_north_inv[0,32993:37442])
# print('Median surf cur coda3 north: ', median_surf_cur_coda3_north)
# median_surf_cur_coda3_north_roms = np.median(v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3])
# print('Median ssurf cur coda3 north roms: ', median_surf_cur_coda3_north_roms)
# # Means
# # East
# mean_surf_cur_coda2_east = np.nanmean(coda2_east_inv[0,33240:])
# print('mean surf cur coda2 east: ', mean_surf_cur_coda2_east)
# mean_surf_cur_coda2_east_roms = np.nanmean(u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2])
# print('mean surf cur coda2 east roms: ', mean_surf_cur_coda2_east_roms)
# mean_surf_cur_coda3_east = np.nanmean(coda3_east_inv[0,32993:37442])
# print('mean surf cur coda3 east: ', mean_surf_cur_coda3_east)
# mean_surf_cur_coda3_east_roms = np.nanmean(u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3])
# print('mean ssurf cur coda3 east roms: ', mean_surf_cur_coda3_east_roms)
# # North
# mean_surf_cur_coda2_north = np.nanmean(coda2_north_inv[0,33240:])
# print('mean surf cur coda2 north: ', mean_surf_cur_coda2_north)
# mean_surf_cur_coda2_north_roms = np.nanmean(v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2])
# print('mean surf cur coda2 north roms: ', mean_surf_cur_coda2_north_roms)
# mean_surf_cur_coda3_north = np.nanmean(coda3_north_inv[0,32993:37442])
# print('mean surf cur coda3 north: ', mean_surf_cur_coda3_north)
# mean_surf_cur_coda3_north_roms = np.nanmean(v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3])
# print('mean ssurf cur coda3 north roms: ', mean_surf_cur_coda3_north_roms)
# # Standard deviation
# # East
# stddev_surf_cur_coda2_east = np.nanstd(coda2_east_inv[0,33240:])
# print('stddev surf cur coda2 east: ', stddev_surf_cur_coda2_east)
# stddev_surf_cur_coda2_east_roms = np.nanstd(u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2])
# print('stddev surf cur coda2 east roms: ',stddev_surf_cur_coda2_east_roms)
# stddev_surf_cur_coda3_east = np.nanstd(coda3_east_inv[0,32993:37442])
# print('stddev surf cur coda3 east: ', stddev_surf_cur_coda3_east)
# stddev_surf_cur_coda3_east_roms = np.nanstd(u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3])
# print('stddev ssurf cur coda3 east roms: ', stddev_surf_cur_coda3_east_roms)
# # North
# stddev_surf_cur_coda2_north = np.nanstd(coda2_north_inv[0,33240:])
# print('stddev surf cur coda2 north: ', stddev_surf_cur_coda2_north)
# stddev_surf_cur_coda2_north_roms = np.nanstd(v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2])
# print('stddev surf cur coda2 north roms: ', stddev_surf_cur_coda2_north_roms)
# stddev_surf_cur_coda3_north = np.nanstd(coda3_north_inv[0,32993:37442])
# print('stddev surf cur coda3 north: ', stddev_surf_cur_coda3_north)
# stddev_surf_cur_coda3_north_roms = np.nanstd(v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3])
# print('stddev ssurf cur coda3 north roms: ', stddev_surf_cur_coda3_north_roms)
# =============================================================================
# =============================================================================
# # Middle row (bottom)
# # Medians
# # East
# median_bot_cur_coda2_east = np.nanmedian(coda2_east_inv[-2,33240:])
# print('Median bot cur coda2 east: ', median_bot_cur_coda2_east)
# median_bot_cur_coda2_east_roms = np.median(u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2])
# print('Median bot cur coda2 east roms: ', median_bot_cur_coda2_east_roms)
# median_bot_cur_coda3_east = np.nanmedian(coda3_east_inv[-3,32993:37442])
# print('Median bot cur coda3 east: ', median_bot_cur_coda3_east)
# median_bot_cur_coda3_east_roms = np.median(u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3])
# print('Median bot cur coda3 east roms: ', median_bot_cur_coda3_east_roms)
# # North
# median_bot_cur_coda2_north = np.nanmedian(coda2_north_inv[-2,33240:])
# print('Median bot cur coda2 north: ', median_bot_cur_coda2_north)
# median_bot_cur_coda2_north_roms = np.median(v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2])
# print('Median bot cur coda2 north roms: ', median_bot_cur_coda2_north_roms)
# median_bot_cur_coda3_north = np.nanmedian(coda3_north_inv[-3,32993:37442])
# print('Median bot cur coda3 north: ', median_bot_cur_coda3_north)
# median_bot_cur_coda3_north_roms = np.median(v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3])
# print('Median bot cur coda3 north roms: ', median_bot_cur_coda3_north_roms)
# # Means
# # East
# mean_bot_cur_coda2_east = np.nanmean(coda2_east_inv[-2,33240:])
# print('mean bot cur coda2 east: ', mean_bot_cur_coda2_east)
# mean_bot_cur_coda2_east_roms = np.nanmean(u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2])
# print('mean bot cur coda2 east roms: ', mean_bot_cur_coda2_east_roms)
# mean_bot_cur_coda3_east = np.nanmean(coda3_east_inv[-3,32993:37442])
# print('mean bot cur coda3 east: ', mean_bot_cur_coda3_east)
# mean_bot_cur_coda3_east_roms = np.nanmean(u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3])
# print('mean bot cur coda3 east roms: ', mean_bot_cur_coda3_east_roms)
# # North
# mean_bot_cur_coda2_north = np.nanmean(coda2_north_inv[-2,33240:])
# print('mean bot cur coda2 north: ', mean_bot_cur_coda2_north)
# mean_bot_cur_coda2_north_roms = np.nanmean(v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2])
# print('mean bot cur coda2 north roms: ', mean_bot_cur_coda2_north_roms)
# mean_bot_cur_coda3_north = np.nanmean(coda3_north_inv[-3,32993:37442])
# print('mean bot cur coda3 north: ', mean_bot_cur_coda3_north)
# mean_bot_cur_coda3_north_roms = np.nanmean(v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3])
# print('mean bot cur coda3 north roms: ', mean_bot_cur_coda3_north_roms)
# # Standard deviation
# # East
# stddev_bot_cur_coda2_east = np.nanstd(coda2_east_inv[-2,33240:])
# print('stddev bot cur coda2 east: ', stddev_bot_cur_coda2_east)
# stddev_bot_cur_coda2_east_roms = np.nanstd(u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2])
# print('stddev bot cur coda2 east roms: ',stddev_bot_cur_coda2_east_roms)
# stddev_bot_cur_coda3_east = np.nanstd(coda3_east_inv[-3,32993:37442])
# print('stddev bot cur coda3 east: ', stddev_bot_cur_coda3_east)
# stddev_bot_cur_coda3_east_roms = np.nanstd(u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3])
# print('stddev bot cur coda3 east roms: ', stddev_bot_cur_coda3_east_roms)
# # North
# stddev_bot_cur_coda2_north = np.nanstd(coda2_north_inv[-2,33240:])
# print('stddev bot cur coda2 north: ', stddev_bot_cur_coda2_north)
# stddev_bot_cur_coda2_north_roms = np.nanstd(v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2])
# print('stddev bot cur coda2 north roms: ', stddev_bot_cur_coda2_north_roms)
# stddev_bot_cur_coda3_north = np.nanstd(coda3_north_inv[-3,32993:37442])
# print('stddev bot cur coda3 north: ', stddev_bot_cur_coda3_north)
# stddev_bot_cur_coda3_north_roms = np.nanstd(v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3])
# print('stddev bot cur coda3 north roms: ', stddev_bot_cur_coda3_north_roms)
# =============================================================================
# =============================================================================
# # Bottom row (depth-averaged)
# # Medians
# # East
# median_davg_cur_coda2_east = np.nanmedian(coda2_east_inv_depthavg[33240:])
# print('Median davg cur coda2 east: ', median_davg_cur_coda2_east)
# median_davg_cur_coda2_east_roms = np.median(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2])
# print('Median davg cur coda2 east roms: ', median_davg_cur_coda2_east_roms)
# median_davg_cur_coda3_east = np.nanmedian(coda3_east_inv_depthavg[32993:37442])
# print('Median davg cur coda3 east: ', median_davg_cur_coda3_east)
# median_davg_cur_coda3_east_roms = np.median(u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3])
# print('Median davg cur coda3 east roms: ', median_davg_cur_coda3_east_roms)
# # North
# median_davg_cur_coda2_north = np.nanmedian(coda2_north_inv_depthavg[33240:])
# print('Median davg cur coda2 north: ', median_davg_cur_coda2_north)
# median_davg_cur_coda2_north_roms = np.median(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2])
# print('Median davg cur coda2 north roms: ', median_davg_cur_coda2_north_roms)
# median_davg_cur_coda3_north = np.nanmedian(coda3_north_inv_depthavg[32993:37442])
# print('Median davg cur coda3 north: ', median_davg_cur_coda3_north)
# median_davg_cur_coda3_north_roms = np.median(v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3])
# print('Median davg cur coda3 north roms: ', median_davg_cur_coda3_north_roms)
# # Means
# # East velocity
# mean_davg_cur_coda2_east = np.nanmean(coda2_east_inv_depthavg[33240:])
# print('mean davg cur coda2 east: ', mean_davg_cur_coda2_east)
# mean_davg_cur_coda2_east_roms = np.nanmean(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2])
# print('mean davg cur coda2 east roms: ', mean_davg_cur_coda2_east_roms)
# mean_davg_cur_coda3_east = np.nanmean(coda3_east_inv_depthavg[32993:37442])
# print('mean davg cur coda3 east: ', mean_davg_cur_coda3_east)
# mean_davg_cur_coda3_east_roms = np.nanmean(u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3])
# print('mean davg cur coda3 east roms: ', mean_davg_cur_coda3_east_roms)
# # East speed
# mean_davg_cur_coda2_east_spd = np.nanmean(abs(coda2_east_inv_depthavg[33240:]))
# print('mean davg cur coda2 east spd: ', mean_davg_cur_coda2_east_spd)
# mean_davg_cur_coda2_east_roms_spd = np.nanmean(abs(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2]))
# print('mean davg cur coda2 east roms spd: ', mean_davg_cur_coda2_east_roms_spd)
# mean_davg_cur_coda3_east_spd = np.nanmean(abs(coda3_east_inv_depthavg[32993:37442]))
# print('mean davg cur coda3 east spd: ', mean_davg_cur_coda3_east_spd)
# mean_davg_cur_coda3_east_roms_spd = np.nanmean(abs(u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3]))
# print('mean davg cur coda3 east roms spd: ', mean_davg_cur_coda3_east_roms_spd)
# # North velocity
# mean_davg_cur_coda2_north = np.nanmean(coda2_north_inv_depthavg[33240:])
# print('mean davg cur coda2 north: ', mean_davg_cur_coda2_north)
# mean_davg_cur_coda2_north_roms = np.nanmean(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2])
# print('mean davg cur coda2 north roms: ', mean_davg_cur_coda2_north_roms)
# mean_davg_cur_coda3_north = np.nanmean(coda3_north_inv_depthavg[32993:37442])
# print('mean davg cur coda3 north: ', mean_davg_cur_coda3_north)
# mean_davg_cur_coda3_north_roms = np.nanmean(v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3])
# print('mean davg cur coda3 north roms: ', mean_davg_cur_coda3_north_roms)
# # North speed
# mean_davg_cur_coda2_north_spd = np.nanmean(abs(coda2_north_inv_depthavg[33240:]))
# print('mean davg cur coda2 north spd: ', mean_davg_cur_coda2_north_spd)
# mean_davg_cur_coda2_north_roms_spd = np.nanmean(abs(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2]))
# print('mean davg cur coda2 north roms spd: ', mean_davg_cur_coda2_north_roms_spd)
# mean_davg_cur_coda3_north_spd = np.nanmean(abs(coda3_north_inv_depthavg[32993:37442]))
# print('mean davg cur coda3 north spd: ', mean_davg_cur_coda3_north_spd)
# mean_davg_cur_coda3_north_roms_spd = np.nanmean(abs(v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3]))
# print('mean davg cur coda3 north roms spd: ', mean_davg_cur_coda3_north_roms_spd)
# # Standard deviation
# # East
# stddev_davg_cur_coda2_east = np.nanstd(coda2_east_inv_depthavg[33240:])
# print('stddev davg cur coda2 east: ', stddev_davg_cur_coda2_east)
# stddev_davg_cur_coda2_east_roms = np.nanstd(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2])
# print('stddev davg cur coda2 east roms: ',stddev_davg_cur_coda2_east_roms)
# stddev_davg_cur_coda3_east = np.nanstd(coda3_east_inv_depthavg[32993:37442])
# print('stddev davg cur coda3 east: ', stddev_davg_cur_coda3_east)
# stddev_davg_cur_coda3_east_roms = np.nanstd(u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3])
# print('stddev davg cur coda3 east roms: ', stddev_davg_cur_coda3_east_roms)
# # East speed
# stddev_davg_cur_coda2_east_spd = np.nanstd(abs(coda2_east_inv_depthavg[33240:]))
# print('stddev davg cur coda2 east spd: ', stddev_davg_cur_coda2_east_spd)
# stddev_davg_cur_coda2_east_roms_spd = np.nanstd(abs(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2]))
# print('stddev davg cur coda2 east roms spd: ', stddev_davg_cur_coda2_east_roms_spd)
# stddev_davg_cur_coda3_east_spd = np.nanstd(abs(coda3_east_inv_depthavg[32993:37442]))
# print('stddev davg cur coda3 east spd: ', stddev_davg_cur_coda3_east_spd)
# stddev_davg_cur_coda3_east_roms_spd = np.nanstd(abs(u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3]))
# print('stddev davg cur coda3 east roms spd: ', stddev_davg_cur_coda3_east_roms_spd)
# # North
# stddev_davg_cur_coda2_north = np.nanstd(coda2_north_inv_depthavg[33240:])
# print('stddev davg cur coda2 north: ', stddev_davg_cur_coda2_north)
# stddev_davg_cur_coda2_north_roms = np.nanstd(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2])
# print('stddev davg cur coda2 north roms: ', stddev_davg_cur_coda2_north_roms)
# stddev_davg_cur_coda3_north = np.nanstd(coda3_north_inv_depthavg[32993:37442])
# print('stddev davg cur coda3 north: ', stddev_davg_cur_coda3_north)
# stddev_davg_cur_coda3_north_roms = np.nanstd(v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3])
# print('stddev davg cur coda3 north roms: ', stddev_davg_cur_coda3_north_roms)
# # North speed
# stddev_davg_cur_coda2_north_spd = np.nanstd(abs(coda2_north_inv_depthavg[33240:]))
# print('stddev davg cur coda2 north spd: ', stddev_davg_cur_coda2_north_spd)
# stddev_davg_cur_coda2_north_roms_spd = np.nanstd(abs(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2]))
# print('stddev davg cur coda2 north roms spd: ', stddev_davg_cur_coda2_north_roms_spd)
# stddev_davg_cur_coda3_north_spd = np.nanstd(abs(coda3_north_inv_depthavg[32993:37442]))
# print('stddev davg cur coda3 north spd: ', stddev_davg_cur_coda3_north_spd)
# stddev_davg_cur_coda3_north_roms_spd = np.nanstd(abs(v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3]))
# print('stddev davg cur coda3 north roms spd: ', stddev_davg_cur_coda3_north_roms_spd)
# =============================================================================


# ------------------------ Violin Plots for Waves ----------------------------
# Load in the wave input data for ROMS - use the one that was used for dbsed0007 
# which includes swell (ERA5 data)
#roms_wave_2019 = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/wave_forcing_file_kaktovik_shelf_era5_data006.nc')
roms_wave_2020 = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/Forcing_files/wave_forcing_file_kaktovik_shelf_ww3_2020_data002.nc')
#roms_wave_2020 = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Beaufort_Shelf_Rivers_proj_002/Model_Input/wave_forcing_file_kaktovik_shelf_cop_2020_data002.nc') 


# Convert the wave time 
wave_time_2020 = pd.to_datetime(roms_wave_2020.wave_time.values+86400, origin='1999-12-31', unit='s') 

# Significant Wave Height 
# Combine all the data into a list
vio_wave_data_swh_all = [mooring_data['S2A1burst']['sigwaveheight'][0,5543:],
                         roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2], # old index [:685]
                         mooring_data['S3A1burst']['sigwaveheight'][0,5503:], 
                         roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3]] # old index [:260]

# Make the figure 
# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize) 
matplotlib.rc('ytick', labelsize=fontsize)
fig5, ax5 = plt.subplots(figsize=(21,5))
s7 = sns.violinplot(data=vio_wave_data_swh_all, ax=ax5, palette=palette)
# CODA
#ax3.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
ax5.set_xticklabels(['CODA2', 'ERA5 at CODA2', 'CODA3','ERA5 at CODA3'])
ax5.set_title('CODA & ERA5 Significant Wave Height', fontsize=fontsize)
ax5.set_ylabel('Significant Wave Height (m)', fontsize=fontsize)


# Peak Wave Period 
# Combine all the data into a list
vio_wave_data_peak_all = [mooring_data['S2A1burst']['peakwaveperiod'][0,5543:],
                         roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2], # old index [:685]
                         mooring_data['S3A1burst']['peakwaveperiod'][0,5503:], 
                         roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3]] # old index [:260]

# Make the figure 
fig6, ax6 = plt.subplots(figsize=(21,5))
s8 = sns.violinplot(data=vio_wave_data_peak_all, ax=ax6, palette=palette)
# CODA
#ax3.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
ax6.set_xticklabels(['CODA2', 'ERA5 at CODA2', 'CODA3','ERA5 at CODA3'])
ax6.set_title('CODA & ERA5 Peak Wave Period', fontsize=fontsize)
ax6.set_ylabel('Peak Wave Period (s)', fontsize=fontsize)


# Panel version of all the plots above 
# Make the figure 
fontsize7 = 12 #18 (make y-label spacing 60 and 50 if doing 18 fontsize)
# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize7) 
matplotlib.rc('ytick', labelsize=fontsize7)
palette2 = ['orange', 'orange', 'orchid','orchid'] # try better colors
fig7, ax7 = plt.subplots(2, figsize=(10,4)) # 16,6
# Significant wave height 
s9 = sns.violinplot(data=vio_wave_data_swh_all, ax=ax7[0], palette=palette2,
                    linecolor='k')
# Edit colors and elements of whisker plot inside violin
# Whisker on first violin
s9.get_children()[1].set_color('white')
# Box on first violin
s9.get_children()[2].set_color('k')
# Dot on first violin
s9.get_children()[3].set_color('white')
# Whisker on second violin
s9.get_children()[5].set_color('white') 
# Box on second violin
s9.get_children()[6].set_color('k')
# Dot on second violin
s9.get_children()[7].set_color('white')
# Whisker on third violin
s9.get_children()[9].set_color('white') 
# Box on third violin
s9.get_children()[10].set_color('k')
# Dot on third violin
s9.get_children()[11].set_color('white')
# Whisker on fourth violin
s9.get_children()[13].set_color('white') 
# Box on fourth violin
s9.get_children()[14].set_color('k')
# Dot on fourth violin
s9.get_children()[15].set_color('white')

# Peak wave period
s10 = sns.violinplot(data=vio_wave_data_peak_all, ax=ax7[1], palette=palette2)
# Edit colors and elements of whisker plot inside violin
# Whisker on first violin
s10.get_children()[1].set_color('white')
# Box on first violin
s10.get_children()[2].set_color('k')
# Dot on first violin
s10.get_children()[3].set_color('white')
# Whisker on second violin
s10.get_children()[5].set_color('white') 
# Box on second violin
s10.get_children()[6].set_color('k')
# Dot on second violin
s10.get_children()[7].set_color('white')
# Whisker on third violin
s10.get_children()[9].set_color('white') 
# Box on third violin
s10.get_children()[10].set_color('k')
# Dot on third violin
s10.get_children()[11].set_color('white')
# Whisker on fourth violin
s10.get_children()[13].set_color('white') 
# Box on fourth violin
s10.get_children()[14].set_color('k')
# Dot on fourth violin
s10.get_children()[15].set_color('white')
# CODA
#ax3.set_xticklabels(['CODA2 East', 'CODA2 North', 'CODA3 East', 'CODA3 North'])
# All
#ax7[1].set_xticklabels(['CODA2', 'ERA5 at CODA2', 'CODA3','ERA5 at CODA3'])
ax7[1].set_xticklabels(['Observed \nat CODA2', 'Modeled \nat CODA2', 'Observed \nat CODA3','Modeled \nat CODA3'])
#ax4.set_title('CODA & ROMS Depth-Averaged Currents', fontsize=fontsize)
plt.setp(ax7[0].get_xticklabels(), visible=False)
ax7[0].set_ylabel('Significant \nWave \nHeight \n(m)', fontsize=fontsize7, 
                  rotation=0, labelpad=40, va='center')
ax7[1].set_ylabel('Peak \nWave \nPeriod \n(s)', fontsize=fontsize7, rotation=0, 
                  labelpad=30, va='center')

plt.subplots_adjust(hspace=0.1)

# Label the subplots 
plt.text(0.865, 0.815, 'a)', fontsize=fontsize7, fontweight='bold', transform=plt.gcf().transFigure) #(0.875, 0.825)
plt.text(0.865, 0.420, 'b)', fontsize=fontsize7, fontweight='bold', transform=plt.gcf().transFigure) #(0.875, 0.825)


# =============================================================================
# # Print some stats from the above plot
# # Significant wave height
# # Medians
# median_swh_coda2 = np.nanmedian(mooring_data['S2A1burst']['sigwaveheight'][0,5543:]) # 5517 --> 5543; 5477 --> 5503
# print('Median swh coda2: ', median_swh_coda2)
# median_swh_coda2_ww3 = np.median(roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2]) # 685 --> 2030; 260 --> 755
# print('Median swh coda2 ww3: ', median_swh_coda2_ww3)
# median_swh_coda3 = np.nanmedian(mooring_data['S3A1burst']['sigwaveheight'][0,5503:])
# print('Median swh coda3: ', median_swh_coda3)
# median_swh_coda3_ww3 = np.median(roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3])
# print('Median swh coda3 ww3: ', median_swh_coda3_ww3)
# # Means
# mean_swh_coda2 = np.nanmean(mooring_data['S2A1burst']['sigwaveheight'][0,5543:])
# print('Mean swh coda2: ', mean_swh_coda2)
# mean_swh_coda2_ww3 = np.mean(roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2])
# print('Mean swh coda2 ww3: ', mean_swh_coda2_ww3)
# mean_swh_coda3 = np.nanmean(mooring_data['S3A1burst']['sigwaveheight'][0,5503:])
# print('Mean swh coda3: ', mean_swh_coda3)
# mean_swh_coda3_ww3 = np.mean(roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3])
# print('Mean swh coda3 ww3: ', mean_swh_coda3_ww3)
# # Standard deviations
# stddev_swh_coda2 = np.nanstd(mooring_data['S2A1burst']['sigwaveheight'][0,5543:])
# print('stddev swh coda2: ', stddev_swh_coda2)
# stddev_swh_coda2_ww3 = np.nanstd(roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2])
# print('stddev swh coda2 ww3: ', stddev_swh_coda2_ww3)
# stddev_swh_coda3 = np.nanstd(mooring_data['S3A1burst']['sigwaveheight'][0,5503:])
# print('stddev swh coda3: ', stddev_swh_coda3)
# stddev_swh_coda3_ww3 = np.nanstd(roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3])
# print('stddev swh coda3 ww3: ', stddev_swh_coda3_ww3)
# # 25%
# perc_25_swh_coda2 = np.nanpercentile(mooring_data['S2A1burst']['sigwaveheight'][0,5543:], 25)
# print('25th percentile swh coda2: ', perc_25_swh_coda2)
# perc_25_swh_coda2_ww3 = np.nanpercentile(roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2], 25)
# print('25th percentile swh coda2 ww3: ', perc_25_swh_coda2_ww3)
# perc_25_swh_coda3 = np.nanpercentile(mooring_data['S3A1burst']['sigwaveheight'][0,5503:], 25)
# print('25th percentile swh coda3: ', perc_25_swh_coda3)
# perc_25_swh_coda3_ww3 = np.nanpercentile(roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3], 25)
# print('25th percentile swh coda3 ww3: ', perc_25_swh_coda3_ww3)
# # 75%
# perc_75_swh_coda2 = np.nanpercentile(mooring_data['S2A1burst']['sigwaveheight'][0,5543:], 75)
# print('75th percentile swh coda2: ', perc_75_swh_coda2)
# perc_75_swh_coda2_ww3 = np.nanpercentile(roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2], 75)
# print('75th percentile swh coda2 ww3: ', perc_75_swh_coda2_ww3)
# perc_75_swh_coda3 = np.nanpercentile(mooring_data['S3A1burst']['sigwaveheight'][0,5503:], 75)
# print('75th percentile swh coda3: ', perc_75_swh_coda3)
# perc_75_swh_coda3_ww3 = np.nanpercentile(roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3], 75)
# print('75th percentile swh coda3 ww3: ', perc_75_swh_coda3_ww3)
# 
# # Peak Wave Period
# # Medians
# median_pwp_coda2 = np.nanmedian(mooring_data['S2A1burst']['peakwaveperiod'][0,5543:])
# print('Median pwp coda2: ', median_pwp_coda2)
# median_pwp_coda2_ww3 = np.median(roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2])
# print('Median pwp coda2 ww3: ', median_pwp_coda2_ww3)
# median_pwp_coda3 = np.nanmedian(mooring_data['S3A1burst']['peakwaveperiod'][0,5503:])
# print('Median pwp coda3: ', median_pwp_coda3)
# median_pwp_coda3_ww3 = np.median(roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3])
# print('Median pwp coda3 ww3: ', median_pwp_coda3_ww3)
# # Means
# mean_pwp_coda2 = np.nanmean(mooring_data['S2A1burst']['peakwaveperiod'][0,5543:])
# print('Mean pwp coda2: ', mean_pwp_coda2)
# mean_pwp_coda2_ww3 = np.mean(roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2])
# print('Mean pwp coda2 ww3: ', mean_pwp_coda2_ww3)
# mean_pwp_coda3 = np.nanmean(mooring_data['S3A1burst']['peakwaveperiod'][0,5503:])
# print('Mean pwp coda3: ', mean_pwp_coda3)
# mean_pwp_coda3_ww3 = np.mean(roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3])
# print('Mean pwp coda3 ww3: ', mean_pwp_coda3_ww3)
# # Standard deviations
# stddev_pwp_coda2 = np.nanstd(mooring_data['S2A1burst']['peakwaveperiod'][0,5543:])
# print('stddev pwp coda2: ', stddev_pwp_coda2)
# stddev_pwp_coda2_ww3 = np.nanstd(roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2])
# print('stddev pwp coda2 ww3: ', stddev_pwp_coda2_ww3)
# stddev_pwp_coda3 = np.nanstd(mooring_data['S3A1burst']['peakwaveperiod'][0,5503:])
# print('stddev pwp coda3: ', stddev_pwp_coda3)
# stddev_pwp_coda3_ww3 = np.nanstd(roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3])
# print('stddev pwp coda3 ww3: ', stddev_pwp_coda3_ww3)
# # 25%
# perc_25_pwp_coda2 = np.nanpercentile(mooring_data['S2A1burst']['peakwaveperiod'][0,5543:], 25)
# print('25th percentile pwp coda2: ', perc_25_pwp_coda2)
# perc_25_pwp_coda2_ww3 = np.nanpercentile(roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2], 25)
# print('25th percentile pwp coda2 ww3: ', perc_25_pwp_coda2_ww3)
# perc_25_pwp_coda3 = np.nanpercentile(mooring_data['S3A1burst']['peakwaveperiod'][0,5503:], 25)
# print('25th percentile pwp coda3: ', perc_25_pwp_coda3)
# perc_25_pwp_coda3_ww3 = np.nanpercentile(roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3], 25)
# print('25th percentile pwp coda3 ww3: ', perc_25_pwp_coda3_ww3)
# # 75%
# perc_75_pwp_coda2 = np.nanpercentile(mooring_data['S2A1burst']['peakwaveperiod'][0,5543:], 75)
# print('75th percentile pwp coda2: ', perc_75_pwp_coda2)
# perc_75_pwp_coda2_ww3 = np.nanpercentile(roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2], 75)
# print('75th percentile pwp coda2 ww3: ', perc_75_pwp_coda2_ww3)
# perc_75_pwp_coda3 = np.nanpercentile(mooring_data['S3A1burst']['peakwaveperiod'][0,5503:], 75)
# print('75th percentile pwp coda3: ', perc_75_pwp_coda3)
# perc_75_pwp_coda3_ww3 = np.nanpercentile(roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3], 75)
# print('75th percentile pwp coda3 ww3: ', perc_75_pwp_coda3_ww3)
# =============================================================================


# -------------------------------------------------------------------------------
# ----------------------------- More Statistics --------------------------------
# -------------------------------------------------------------------------------
# Calculate and print correlations and RMSE for the different datasets 
from scipy import stats

# Currents 
# Resample the data to match

# Which is coarser in time?

# Pull out matching time for moorings
# CODA Mooring 2 (July 1 - September 23)
# Pull out July 1 - September 24, 2020 (until end) in coda data 
coda2_time_julsep_2020 = moor1_time_julsep_2020.copy()

# CODA Mooring 3 (July 1 - August 1)
coda3_time_julaug_2020 = moor2_time_julaug_2020.copy()

# ROMS current time at 
# CODA 2
roms_time_at_coda2 = u_roms_interp_coda2_all_output.ocean_time[:677]
# CODA 3
roms_time_at_coda3 = u_roms_interp_coda2_all_output.ocean_time[:250]



# ROMS wave time at 
# CODA 2
wave_time_at_coda2 = wave_time_2020[:2030]
# CODA 3
wave_time_at_coda3 = wave_time_2020[:755]


# Coarsest time is model output at four hourly so
# resample current observations to match model output of currents

# Surface 
vio_data_surf_all = [coda2_east_inv[0,33240:], u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2], coda2_north_inv[0,33240:], 
                     v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2], coda3_east_inv[0,32993:37442], 
                     u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3], coda3_north_inv[0,32993:37442],
                     v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3]]
# Bottom (-1 index gives all nan)
vio_data_bot_all = [coda2_east_inv[-2,33240:], u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2], coda2_north_inv[-2,33240:], 
                     v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2], coda3_east_inv[-3,32993:37442], 
                     u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3], coda3_north_inv[-3,32993:37442],
                     v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3]]
# Depth-Averaged 
vio_data_depthavg_all = [coda2_east_inv_depthavg[33240:], u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2], 
                         coda2_north_inv_depthavg[33240:], v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2], 
                         coda3_east_inv_depthavg[32993:37442], u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3], 
                         coda3_north_inv_depthavg[32993:37442], v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3]]


# ------- Currents ------
# Focus on just depth-averaged since we have to pick one...

# Resample current data from moorings 
# CODA 2
coda2_east_davg_pd = pd.Series(coda2_east_inv_depthavg[33240:], index=coda2_time_julsep_2020)
coda2_east_davg_pd_hourly = coda2_east_davg_pd.resample('3H').mean()
coda2_north_davg_pd = pd.Series(coda2_north_inv_depthavg[33240:], index=coda2_time_julsep_2020)
coda2_north_davg_pd_hourly = coda2_north_davg_pd.resample('3H').mean()

# =============================================================================
# coda2_east_davg_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
#                        for elem in coda2_east_davg_pd_hourly],
#                       dtype=float)
# coda2_north_davg_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
#                        for elem in coda2_north_davg_pd_hourly],
#                       dtype=float)
# =============================================================================
# CODA 3
coda3_east_davg_pd = pd.Series(coda3_east_inv_depthavg[32993:37442], index=coda3_time_julaug_2020[25:-50])
coda3_east_davg_pd_hourly = coda3_east_davg_pd.resample('3H').mean()
coda3_north_davg_pd = pd.Series(coda3_north_inv_depthavg[32993:37442], index=coda3_time_julaug_2020[25:-50])
coda3_north_davg_pd_hourly = coda3_north_davg_pd.resample('3H').mean()

# =============================================================================
# coda3_east_davg_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
#                        for elem in coda3_east_davg_pd_hourly],
#                       dtype=float)
# coda3_north_davg_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
#                        for elem in coda3_north_davg_pd_hourly],
#                       dtype=float)
# =============================================================================

# =============================================================================
# # Mask nans in both ROMS and CODA
# # CODA 2
# mask_coda2_east_davg = ~np.isnan(coda2_east_davg_flat[1:]) & ~np.isnan(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2])
# mask_coda2_north_davg = ~np.isnan(coda2_north_davg_flat[1:]) & ~np.isnan(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2])
# # CODA 3
# mask_coda3_east_davg = ~np.isnan(coda3_east_davg_flat) & ~np.isnan(u_roms_interp_coda3_all_output_depthavg.u[:249,eta_rho_coda3,xi_rho_coda3])
# mask_coda3_north_davg = ~np.isnan(coda3_north_davg_flat) & ~np.isnan(v_roms_interp_coda3_all_output_depthavg.v[:249,eta_rho_coda3,xi_rho_coda3])
# 
# # Mask out nans 
# # CODA 2
# coda2_east_davg_flat_clean = coda2_east_davg_flat[1:][mask_coda2_east_davg]
# roms_east_davg_coda2_clean = u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2][mask_coda2_east_davg]
# coda2_north_davg_flat_clean = coda2_north_davg_flat[1:][mask_coda2_north_davg]
# roms_north_davg_coda2_clean = v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2][mask_coda2_north_davg]
# # CODA 3
# coda3_east_davg_flat_clean = coda3_east_davg_flat[mask_coda3_east_davg]
# roms_east_davg_coda3_clean = u_roms_interp_coda3_all_output_depthavg.u[:249,eta_rho_coda3,xi_rho_coda3][mask_coda3_east_davg]
# coda3_north_davg_flat_clean = coda3_north_davg_flat[mask_coda3_north_davg]
# roms_north_davg_coda3_clean = v_roms_interp_coda3_all_output_depthavg.v[:249,eta_rho_coda3,xi_rho_coda3][mask_coda3_north_davg]
# =============================================================================


# Calculate and print correlations
# CODA2
coda2_east_r, coda2_east_p = stats.pearsonr(coda2_east_davg_pd_hourly[1:], u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2])
print('coda2_east_r: ', coda2_east_r)
print('coda2_east_p: ', coda2_east_p)
coda2_north_r, coda2_north_p = stats.pearsonr(coda2_north_davg_pd_hourly[1:], v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2])
print('coda2_north_r: ', coda2_north_r)
print('coda2_north_p: ', coda2_north_p)

# CODA3
coda3_east_r, coda3_east_p = stats.pearsonr(coda3_east_davg_pd_hourly, u_roms_interp_coda3_all_output_depthavg.u[:249,eta_rho_coda3,xi_rho_coda3])
print('coda3_east_r: ', coda3_east_r)
print('coda3_east_p: ', coda3_east_p)
coda3_north_r, coda3_north_p = stats.pearsonr(coda3_north_davg_pd_hourly, v_roms_interp_coda3_all_output_depthavg.v[:249,eta_rho_coda3,xi_rho_coda3])
print('coda3_north_r: ', coda3_north_r)
print('coda3_north_p: ', coda3_north_p)

# Calculate and print RMSE
# CODA 2
coda2_east_rmse = sum((abs(u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2]-coda2_east_davg_pd_hourly[1:])**2)/len(coda2_east_davg_pd_hourly[1:]))**0.5
print('coda2_east_rmse: ', coda2_east_rmse.values)
coda2_north_rmse = sum((abs(v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2]-coda2_north_davg_pd_hourly[1:])**2)/len(coda2_north_davg_pd_hourly[1:]))**0.5
print('coda2_north_rmse: ', coda2_north_rmse.values)
# CODA 3
coda3_east_rmse = sum((abs(u_roms_interp_coda3_all_output_depthavg.u[:249,eta_rho_coda3,xi_rho_coda3]-coda3_east_davg_pd_hourly)**2)/len(coda3_east_davg_pd_hourly))**0.5
print('coda3_east_rmse: ', coda3_east_rmse.values)
coda3_north_rmse = sum((abs(v_roms_interp_coda3_all_output_depthavg.v[:249,eta_rho_coda3,xi_rho_coda3]-coda3_north_davg_pd_hourly)**2)/len(coda3_north_davg_pd_hourly))**0.5
print('coda3_north_rmse: ', coda3_north_rmse.values)

# Calculate and print Briar skill score
# bss=1-np.sum((mod-obs)**2)/np.sum((mod-np.mean(obs))**2)  #Briar skill score
# CODA 2
mod_east_coda2 = u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2]
obs_east_coda2 = coda2_east_davg_pd_hourly[1:]
mod_north_coda2 = v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2]
obs_north_coda2 = coda2_north_davg_pd_hourly[1:]
coda2_east_bss = 1-np.sum((mod_east_coda2-obs_east_coda2)**2)/np.sum((mod_east_coda2-np.mean(obs_east_coda2))**2)
print('coda2_east_bss: ', coda2_east_bss.values)
coda2_north_bss = 1-np.sum((mod_north_coda2-obs_north_coda2)**2)/np.sum((mod_north_coda2-np.mean(obs_north_coda2))**2)
print('coda2_north_bss: ', coda2_north_bss.values)
# CODA 3
mod_east_coda3 = u_roms_interp_coda3_all_output_depthavg.u[:249,eta_rho_coda3,xi_rho_coda3]
obs_east_coda3 = coda3_east_davg_pd_hourly
mod_north_coda3 = v_roms_interp_coda3_all_output_depthavg.v[:249,eta_rho_coda3,xi_rho_coda3]
obs_north_coda3 = coda3_north_davg_pd_hourly
coda3_east_bss = 1-np.sum((mod_east_coda3-obs_east_coda3)**2)/np.sum((mod_east_coda3-np.mean(obs_east_coda3))**2)
print('coda3_east_bss: ', coda3_east_bss.values)
coda3_north_bss = 1-np.sum((mod_north_coda3-obs_north_coda3)**2)/np.sum((mod_north_coda3-np.mean(obs_north_coda3))**2)
print('coda3_north_bss: ', coda3_north_bss.values)



# ------ Waves ------
# Resample wave data from moorings to match wave data of model input
# Jk wave data from CODA is already hourly...sooooo fix thissss
# Significant wave height 

# Flatten matlb data to take correlations 
# CODA 2
coda2_hwave_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
                       for elem in mooring_data['S2A1burst']['sigwaveheight'][0, 5543:]],
                      dtype=float)
coda2_pwp_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
                       for elem in mooring_data['S2A1burst']['peakwaveperiod'][0, 5543:]],
                      dtype=float)
# CODA 3
# resample to fill gaps
coda3_hwave_pd = pd.Series(mooring_data['S3A1burst']['sigwaveheight'][0,5503:], index=moor2_datetime_final_leap_burst[5503:])
coda3_hwave_hourly = coda3_hwave_pd.resample('1H').mean()
coda3_pwp_pd = pd.Series(mooring_data['S3A1burst']['peakwaveperiod'][0,5503:], index=moor2_datetime_final_leap_burst[5503:])
coda3_pwp_hourly = coda3_pwp_pd.resample('1H').mean()

coda3_hwave_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
                       for elem in coda3_hwave_hourly],
                      dtype=float)
coda3_pwp_flat = np.array([elem[0, 0] if isinstance(elem, np.ndarray) else np.nan
                       for elem in coda3_pwp_hourly],
                      dtype=float)

# Mask out nans in both arrays (ROMS and CODA)
# CODA 2
mask_coda2_hwave = ~np.isnan(coda2_hwave_flat) & ~np.isnan(roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2])
mask_coda2_pwp = ~np.isnan(coda2_pwp_flat) & ~np.isnan(roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2])
# CODA 3
mask_coda3_hwave = ~np.isnan(coda3_hwave_flat) & ~np.isnan(roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3])
mask_coda3_pwp = ~np.isnan(coda3_pwp_flat) & ~np.isnan(roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3])

# Mask out the nans 
# CODA 2
coda2_hwave_flat_clean = coda2_hwave_flat[mask_coda2_hwave]
roms_hwave_coda2_clean = roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2][mask_coda2_hwave]
coda2_pwp_flat_clean = coda2_pwp_flat[mask_coda2_pwp]
roms_pwp_coda2_clean = roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2][mask_coda2_pwp]
# CODA 3
coda3_hwave_flat_clean = coda3_hwave_flat[mask_coda3_hwave]
roms_hwave_coda3_clean = roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3][mask_coda3_hwave]
coda3_pwp_flat_clean = coda3_pwp_flat[mask_coda3_pwp]
roms_pwp_coda3_clean = roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3][mask_coda3_pwp]


# Calculate and print correlations 
# CODA2
coda2_hwave_r, coda2_hwave_p = stats.pearsonr(coda2_hwave_flat_clean, roms_hwave_coda2_clean)
print('coda2_hwave_r: ', coda2_hwave_r)
print('coda2_hwave_p: ', coda2_hwave_p)
coda2_pwp_r, coda2_pwp_p = stats.pearsonr(coda2_pwp_flat_clean, roms_pwp_coda2_clean)
print('coda2_pwp_r: ', coda2_pwp_r)
print('coda2_pwp_p: ', coda2_pwp_p)

# CODA3
coda3_hwave_r, coda3_hwave_p = stats.pearsonr(coda3_hwave_flat_clean, roms_hwave_coda3_clean)
print('coda3_hwave_r: ', coda3_hwave_r)
print('coda3_hwave_p: ', coda3_hwave_p)
coda3_pwp_r, coda3_pwp_p = stats.pearsonr(coda3_pwp_flat_clean, roms_pwp_coda3_clean)
print('coda3_pwp_r: ', coda3_pwp_r)
print('coda3_pwp_p: ', coda3_pwp_p)

# Calculate and print RMSE 
#rmse=sum((abs(mod-obs)**2)/len(obs))**0.5 #Root mean square error
# CODA 2
coda2_hwave_rmse = sum((abs(roms_hwave_coda2_clean-coda2_hwave_flat_clean)**2)/len(coda2_hwave_flat_clean))**0.5
print('coda2_hwave_rmse: ', coda2_hwave_rmse.values)
coda2_pwp_rmse = sum((abs(roms_pwp_coda2_clean-coda2_pwp_flat_clean)**2)/len(coda2_pwp_flat_clean))**0.5
print('coda2_pwp_rmse: ', coda2_pwp_rmse.values)
# CODA 3
coda3_hwave_rmse = sum((abs(roms_hwave_coda3_clean-coda3_hwave_flat_clean)**2)/len(coda3_hwave_flat_clean))**0.5
print('coda3_hwave_rmse: ', coda3_hwave_rmse.values)
coda3_pwp_rmse = sum((abs(roms_pwp_coda3_clean-coda3_pwp_flat_clean)**2)/len(coda3_pwp_flat_clean))**0.5
print('coda3_pwp_rmse: ', coda3_pwp_rmse.values)

# Calculate and print Briar skill score
# bss=1-np.sum((mod-obs)**2)/np.sum((mod-np.mean(obs))**2)  #Briar skill score
# CODA 2
mod_hwave_coda2 = roms_hwave_coda2_clean
obs_hwave_coda2 = coda2_hwave_flat_clean
mod_pwp_coda2 = roms_pwp_coda2_clean
obs_pwp_coda2 = coda2_pwp_flat_clean
coda2_hwave_bss = 1-np.sum((mod_hwave_coda2-obs_hwave_coda2)**2)/np.sum((mod_hwave_coda2-np.mean(obs_hwave_coda2))**2)
print('coda2_hwave_bss: ', coda2_hwave_bss.values)
coda2_pwp_bss = 1-np.sum((mod_pwp_coda2-obs_pwp_coda2)**2)/np.sum((mod_pwp_coda2-np.mean(obs_pwp_coda2))**2)
print('coda2_pwp_bss: ', coda2_pwp_bss.values)
# CODA 3
mod_hwave_coda3 = roms_hwave_coda3_clean
obs_hwave_coda3 = coda3_hwave_flat_clean
mod_pwp_coda3 = roms_pwp_coda3_clean
obs_pwp_coda3 = coda3_pwp_flat_clean
coda3_hwave_bss = 1-np.sum((mod_hwave_coda3-obs_hwave_coda3)**2)/np.sum((mod_hwave_coda3-np.mean(obs_hwave_coda3))**2)
print('coda3_hwave_bss: ', coda3_hwave_bss.values)
coda3_pwp_bss = 1-np.sum((mod_pwp_coda3-obs_pwp_coda3)**2)/np.sum((mod_pwp_coda3-np.mean(obs_pwp_coda3))**2)
print('coda3_pwp_bss: ', coda3_pwp_bss.values)



# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Waves
roms_wave_swh_pwp = xr.Dataset(
    data_vars=dict(
        Pwave_top_coda2=(['coda2_time'], roms_wave_2020.Pwave_top[:2030, eta_rho_coda2, xi_rho_coda2].values), # 685 --> 2030
        Pwave_top_coda3=(['coda3_time'], roms_wave_2020.Pwave_top[:755, eta_rho_coda3, xi_rho_coda3].values), # 260 --> 755
        Hwave_top_coda2=(['coda2_time'], roms_wave_2020.Hwave[:2030, eta_rho_coda2, xi_rho_coda2].values),
        Hwave_top_coda3=(['coda3_time'], roms_wave_2020.Hwave[:755, eta_rho_coda3, xi_rho_coda3].values),
        ),
    coords=dict(
        coda2_time=wave_time_2020[:2030],
        coda3_time=wave_time_2020[:755], 
        ),
    attrs=dict(description='ROMS wave properties (1) peak wave period - Pwave (second), (2) significant wave height - Hwave (meter) at CODA moorings 2 and 3'))
# Save this to a netcdf
#roms_wave_swh_pwp.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Paper1_Take2/Data/fig3_pwp_swh_coda_02.nc')

# Currents 
roms_currents_davg_surf_bot = xr.Dataset(
    data_vars=dict(
        surf_u_roms_interp_coda2=(['coda2_time'], u_roms_interp_coda2_all_output.u[:677,0,eta_rho_coda2,xi_rho_coda2].values),
        surf_v_roms_interp_coda2=(['coda2_time'], v_roms_interp_coda2_all_output.v[:677,0,eta_rho_coda2,xi_rho_coda2].values),
        surf_u_roms_interp_coda3=(['coda3_time'], u_roms_interp_coda3_all_output.u[:250,0,eta_rho_coda3,xi_rho_coda3].values),
        surf_v_roms_interp_coda3=(['coda3_time'], v_roms_interp_coda3_all_output.v[:250,0,eta_rho_coda3,xi_rho_coda3].values),
        bot_u_roms_interp_coda2=(['coda2_time'], u_roms_interp_coda2_all_output.u[:677,-2,eta_rho_coda2,xi_rho_coda2].values),
        bot_v_roms_interp_coda2=(['coda2_time'], v_roms_interp_coda2_all_output.v[:677,-2,eta_rho_coda2,xi_rho_coda2].values),
        bot_u_roms_interp_coda3=(['coda3_time'], u_roms_interp_coda3_all_output.u[:250,-3,eta_rho_coda3,xi_rho_coda3].values),
        bot_v_roms_interp_coda3=(['coda3_time'], v_roms_interp_coda3_all_output.v[:250,-3,eta_rho_coda3,xi_rho_coda3].values),
        depth_avg_u_roms_interp_coda2=(['coda2_time'], u_roms_interp_coda2_all_output_depthavg.u[:677,eta_rho_coda2,xi_rho_coda2].values),
        depth_avg_v_roms_interp_coda2=(['coda2_time'], v_roms_interp_coda2_all_output_depthavg.v[:677,eta_rho_coda2,xi_rho_coda2].values),
        depth_avg_u_roms_interp_coda3=(['coda3_time'], u_roms_interp_coda3_all_output_depthavg.u[:250,eta_rho_coda3,xi_rho_coda3].values),
        depth_avg_v_roms_interp_coda3=(['coda3_time'], v_roms_interp_coda3_all_output_depthavg.v[:250,eta_rho_coda3,xi_rho_coda3].values),
        ),
    coords=dict(
        coda2_time=u_roms_interp_coda2_all_output.ocean_time[:677].values,
        coda3_time=v_roms_interp_coda2_all_output.ocean_time[:250].values,
        ),
    attrs=dict(description='ROMS currents interpolated onto CODA moorings 2 and 3; Includes surface, bottom, and depth-averaged u and v currents (meter per second)'))
# Save this to a netcdf
#roms_currents_davg_surf_bot.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Paper1_Take2/Data/fig4_surf_bot_depth_avg_roms_currents_at_coda.nc')





# -------------------------------------------------------------------------------
# --------------- Plot CODA waves and currents at CODA 2 and 3 ------------------ 
# -------------------------------------------------------------------------------
# Plot time series of Significant wave height and east-west currents at 
# CODAs 2 and 3 to see if there is overlap between westward currents and 
# high wave heights 

# =============================================================================
# coda2_hwave_flat_clean
# coda3_hwave_flat_clean
# coda2_east_inv_depthavg[33240:]
# coda3_east_inv_depthavg[32993:37442] 
# 
# # Pull out July 1 - September 24, 2020 (until end) in coda data 
# coda2_time_julsep_2020 = moor1_time_julsep_2020.copy()
# # CODA Mooring 3 (July 1 - August 1)
# coda3_time_julaug_2020 = moor2_time_julaug_2020.copy()
# # ROMS current time at 
# # CODA 2
# roms_time_at_coda2 = u_roms_interp_coda2_all_output.ocean_time[:677]
# # CODA 3
# roms_time_at_coda3 = u_roms_interp_coda2_all_output.ocean_time[:250]
# # ROMS wave time at 
# # CODA 2
# wave_time_at_coda2 = wave_time_2020[:2030]
# # CODA 3
# wave_time_at_coda3 = wave_time_2020[:755]
# =============================================================================

# Make the figure 
fig8, ax8 = plt.subplots(2, 2, figsize=(24,8))

# CODA 2
ax8[0,0].set_title('CODA 2', fontsize=fontsize-5)
ax8[0,0].plot(wave_time_at_coda2, coda2_hwave_flat, color='blue', linewidth=2)
ax8[0,0].set_ylabel('Sig. \nWave. \nHeight (m)', rotation=0, va='center', fontsize=fontsize-8, labelpad=45)
ax8[1,0].plot(coda2_time_julsep_2020, coda2_east_inv_depthavg[33240:], color='orange', linewidth=2)
ax8[1,0].axhline(y=0, color='dimgray', linewidth=0.5)
ax8[1,0].set_ylabel('Eastward \nCurrents \n(m/s)', rotation=0, va='center', fontsize=fontsize-8, labelpad=45)
ax8[1,0].set_xlabel('Time', fontsize=fontsize-6)

# CODA 3
ax8[0,1].set_title('CODA 3', fontsize=fontsize-5)
ax8[0,1].plot(wave_time_at_coda3, coda3_hwave_flat, color='blue', linewidth=2)
ax8[0,1].set_ylabel('Sig. \nWave \nHeight (m)', rotation=0, va='center', fontsize=fontsize-8, labelpad=45)
ax8[1,1].plot(coda3_time_julaug_2020[25:-50], coda3_east_inv_depthavg[32993:37442], color='orange', linewidth=2)
ax8[1,1].axhline(y=0, color='dimgray', linewidth=0.5)
ax8[1,1].set_ylabel('Eastward \nCurrents \n(m/s)', rotation=0, va='center', fontsize=fontsize-8, labelpad=45)
ax8[1,1].set_xlabel('Time', fontsize=fontsize-6)

plt.setp(ax8[0,0].get_xticklabels(), visible=False)
#plt.setp(ax8[1].get_xticklabels(), visible=False)
plt.setp(ax8[0,1].get_xticklabels(), visible=False)

plt.subplots_adjust(hspace=0.1, swpace=0.15)







