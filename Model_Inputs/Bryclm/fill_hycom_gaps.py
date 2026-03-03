#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:38:33 2023

@author: brun1463
"""

########################### Fill HYCOM Gaps ###############################
# The purpose of this script is to fill gaps in the HYCOM data 
# using linear interpolation. This script is the same as 
# fill_hycom_gaps.ipynb on alpine but fils gaps in the surface 
# stress data from HYCOM. The script reads in the original data and 
# writes out the gap-filled data.
#
# Notes:
#
#
###########################################################################


# Load in the packages
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

# Set some plotting things 
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20)

# =============================================================================
# # U Surface Stress
# # Load in the HYCOM data 
# surtx_data = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_surtx/HYCOM_2020_surtx.nc4')
# print(surtx_data)
# 
# # Look at the times in this dataset
# print(surtx_data.time[0].values)
# print(surtx_data.time[1].values)
# print(surtx_data.time[-1].values)
# 
# # Make a copy of the data to work with
# surtx_data_cp1 = surtx_data.copy()
# surtx_data_cp1['time'] = pd.to_datetime(surtx_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 3 hours since this is the frequency 
# # of the HYCOM data
# surtx_data_nogap = surtx_data_cp1.resample(time='3H').interpolate('linear')
# print(surtx_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', surtx_data_nogap.time[0:5].values)
# print(surtx_data_nogap.time[-6:-1].values)
# print(surtx_data_nogap.time[-1].values)
# print(len(surtx_data_nogap.time))
# print('\nOriginal data: ', surtx_data.time[0:5].values)
# print(surtx_data.time[-6:-1].values)
# print(surtx_data.time[-1].values)
# print(len(surtx_data.time))
# 
# # Plot the gap-free data to see if this worked - just look at one location 
# fig1, ax1 = plt.subplots(2,figsize=(22,18))
# #fig1.autofmt_xdate()
# 
# # Plot the original hycom data for this location
# ax1[0].scatter(surtx_data.time[1402:2330].values, surtx_data.surtx[1402:2330,53,133].values, color='red', linewidth=3, label='Original')
# ax1[0].set_title('Original U Surface Stress (Pa)', fontsize=30)
# ax1[0].set_xlabel('Time', fontsize=30)
# ax1[0].set_ylabel('U Surface Stress (Pa)', fontsize=30)
# #ax1[0].legend()
# 
# # Plot the gap free data for this location
# ax1[1].scatter(surtx_data_nogap.time[1444:2444].values, surtx_data_nogap.surtx[1444:2444,53,133].values, color='red', linewidth=3)
# ax1[1].set_title('Gap Free U Surface Stress (Pa)', fontsize=30)
# ax1[1].set_xlabel('Time', fontsize=30)
# ax1[1].set_ylabel('U Surface Stress (Pa)', fontsize=30)
# #fig1.set_facecolor('white')
# #fig1.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# surtx_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_surtx/HYCOM_2020_surtx_nogaps.nc4')
# #surtx_data_nogap.to_netcdf('/Users/brun1463/Desktop/HYCOM_2019_surface_downward_eastward_stress_nogaps.nc4')
# 
# # Remove all of this from memory to save RAM
# del(surtx_data)
# del(surtx_data_cp1)
# del(surtx_data_nogap)
# 
# 
# # V Surface Stress
# # Load in the HYCOM data 
# surty_data = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_surty/HYCOM_2020_surty.nc4')
# print(surty_data)
# 
# # Look at the times in this dataset
# print(surty_data.time[0].values)
# print(surty_data.time[1].values)
# print(surty_data.time[-1].values)
# 
# # Make a copy of the data to work with
# surty_data_cp1 = surty_data.copy()
# surty_data_cp1['time'] = pd.to_datetime(surty_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 3 hours since this is the frequency 
# # of the HYCOM data
# surty_data_nogap = surty_data_cp1.resample(time='3H').interpolate('linear')
# print(surty_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', surty_data_nogap.time[0:5].values)
# print(surty_data_nogap.time[-6:-1].values)
# print(surty_data_nogap.time[-1].values)
# print(len(surty_data_nogap.time))
# print('\nOriginal data: ', surty_data.time[0:5].values)
# print(surty_data.time[-6:-1].values)
# print(surty_data.time[-1].values)
# print(len(surty_data.time))
# 
# # Plot the gap-free data to see if this worked - just look at one location 
# fig2, ax2 = plt.subplots(2,figsize=(22,18))
# #fig1.autofmt_xdate()
# 
# # Plot the original hycom data for this location
# ax2[0].scatter(surty_data.time[1402:2330].values, surty_data.surty[1402:2330,53,133].values, color='green', linewidth=3, label='Original')
# ax2[0].set_title('Original V Surface Stress (Pa)', fontsize=30)
# ax2[0].set_xlabel('Time', fontsize=30)
# ax2[0].set_ylabel('U Surface Stress (Pa)', fontsize=30)
# #ax2[0].legend()
# 
# # Plot the gap free data for this location
# ax2[1].scatter(surty_data_nogap.time[1444:2444].values, surty_data_nogap.surty[1444:2444,53,133].values, color='green', linewidth=3)
# ax2[1].set_title('Gap Free V Surface Stress (Pa)', fontsize=30)
# ax2[1].set_xlabel('Time', fontsize=30)
# ax2[1].set_ylabel('V Surface Stress (Pa)', fontsize=30)
# #fig2.set_facecolor('white')
# #fig2.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# surty_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_surty/HYCOM_2020_surty_nogaps.nc4')
# 
# # Remove all of this from memory to save RAM
# del(surty_data)
# del(surty_data_cp1)
# del(surty_data_nogap)
# =============================================================================


# Qtot = Total Surface Heat Flux
# Load in the HYCOM data 
qtot_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_qtot/HYCOM_2020_qtot.nc4')
print(qtot_data)

# Get rid of duplicate times 
qtot_data = qtot_data.sel(time=~qtot_data.get_index("time").duplicated())

# Look at the times in this dataset
print(qtot_data.time[0].values)
print(qtot_data.time[1].values)
print(qtot_data.time[-1].values)

# Make a copy of the data to work with
qtot_data_cp1 = qtot_data.copy()
qtot_data_cp1['time'] = pd.to_datetime(qtot_data_cp1['time'])

# Resample and interpolate the data over time
# Resampling for every 1 hours since this is the frequency 
# of the HYCOM data
qtot_data_nogap = qtot_data_cp1.resample(time='1H').interpolate('linear') # hourly
#print(qtot_data_nogap)

# Compare the times in this new dataset with the times in the original one
# New data 
print('New gap free data: ', qtot_data_nogap.time[0:5].values)
print(qtot_data_nogap.time[-6:-1].values)
print(qtot_data_nogap.time[-1].values)#
print(len(qtot_data_nogap.time))
print('\nOriginal data: ', qtot_data.time[0:5].values)
print(qtot_data.time[-6:-1].values)
print(qtot_data.time[-1].values)
print(len(qtot_data.time))

# Plot the gap-free data to see if this worked - just look at one location 
fig3, ax3 = plt.subplots(2,figsize=(22,18))
#fig1.autofmt_xdate()

# Plot the original hycom data for this location
ax3[0].scatter(qtot_data.time[1402:2330].values, qtot_data.qtot[1402:2330,53,133].values, color='green', linewidth=3, label='Original')
ax3[0].set_title('Original Qtot (W/m2)', fontsize=30)
ax3[0].set_xlabel('Time', fontsize=30)
ax3[0].set_ylabel('Qtot (W/m2)', fontsize=30)
#ax3[0].legend()

# Plot the gap free data for this location
ax3[1].scatter(qtot_data_nogap.time[1444:2444].values, qtot_data_nogap.qtot[1444:2444,53,133].values, color='green', linewidth=3)
ax3[1].set_title('Gap Free Qtot (W/m2)', fontsize=30)
ax3[1].set_xlabel('Time', fontsize=30)
ax3[1].set_ylabel('Qtot (W/m2)', fontsize=30)
#fig3.set_facecolor('white')
#fig3.savefig('hycom_gapfill_surtx01.png')


# Save this new dataset as the new HYCOM zeta data
qtot_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_qtot/HYCOM_2020_qtot_nogaps.nc4')

# Remove all of this from memory to save RAM
del(qtot_data)
del(qtot_data_cp1)
del(qtot_data_nogap)


# =============================================================================
# # Emp = Total Surface Freshwater Flux
# # Load in the HYCOM data 
# emp_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_emp/HYCOM_2020_emp.nc4')
# print(emp_data)
# 
# # Get rid of duplicate times 
# emp_data = emp_data.sel(time=~emp_data.get_index("time").duplicated())
# 
# # Look at the times in this dataset
# print(emp_data.time[0].values)
# print(emp_data.time[1].values)
# print(emp_data.time[-1].values)
# 
# # Make a copy of the data to work with
# emp_data_cp1 = emp_data.copy()
# emp_data_cp1['time'] = pd.to_datetime(emp_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 1 hours since this is the frequency 
# # of the HYCOM data
# emp_data_nogap = emp_data_cp1.resample(time='1H').interpolate('linear') # hourly
# #print(qtot_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', emp_data_nogap.time[0:5].values)
# print(emp_data_nogap.time[-6:-1].values)
# print(emp_data_nogap.time[-1].values)#
# print(len(emp_data_nogap.time))
# print('\nOriginal data: ', emp_data.time[0:5].values)
# print(emp_data.time[-6:-1].values)
# print(emp_data.time[-1].values)
# print(len(emp_data.time))
# 
# # Plot the gap-free data to see if this worked - just look at one location 
# fig4, ax4 = plt.subplots(2,figsize=(22,18))
# #fig1.autofmt_xdate()
# 
# # Plot the original hycom data for this location
# ax4[0].scatter(emp_data.time[1402:2330].values, emp_data.emp[1402:2330,53,133].values, color='green', linewidth=3, label='Original')
# ax4[0].set_title('Original emp (kg/m2/s)', fontsize=30)
# ax4[0].set_xlabel('Time', fontsize=30)
# ax4[0].set_ylabel('emp (kg/m2/s)', fontsize=30)
# #ax4[0].legend()
# 
# # Plot the gap free data for this location
# ax4[1].scatter(emp_data_nogap.time[1444:2444].values, emp_data_nogap.emp[1444:2444,53,133].values, color='green', linewidth=3)
# ax4[1].set_title('Gap Free emp (kg/m2/s)', fontsize=30)
# ax4[1].set_xlabel('Time', fontsize=30)
# ax4[1].set_ylabel('emp (kg/m2/s)', fontsize=30)
# #fig4.set_facecolor('white')
# #fig4.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# emp_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_emp/HYCOM_2020_emp_nogaps.nc4')
# 
# # Remove all of this from memory to save RAM
# del(emp_data)
# del(emp_data_cp1)
# del(emp_data_nogap)
# =============================================================================


# =============================================================================
# # surf_el = Free surface elevation
# # Load in the HYCOM data 
# surf_el_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_surf_el.nc4')
# print(surf_el_data)
# 
# # Get rid of duplicate times 
# surf_el_data = surf_el_data.sel(time=~surf_el_data.get_index("time").duplicated())
# 
# # Look at the times in this dataset
# print(surf_el_data.time[0].values)
# print(surf_el_data.time[1].values)
# print(surf_el_data.time[-1].values)
# 
# # Make a copy of the data to work with
# surf_el_data_cp1 = surf_el_data.copy()
# surf_el_data_cp1['time'] = pd.to_datetime(surf_el_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 1 hours since this is the frequency 
# # of the HYCOM data
# surf_el_data_nogap = surf_el_data_cp1.resample(time='3H').interpolate('linear') # hourly
# #print(qtot_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', surf_el_data_nogap.time[0:5].values)
# print(surf_el_data_nogap.time[-6:-1].values)
# print(surf_el_data_nogap.time[-1].values)#
# print(len(surf_el_data_nogap.time))
# print('\nOriginal data: ', surf_el_data.time[0:5].values)
# print(surf_el_data.time[-6:-1].values)
# print(surf_el_data.time[-1].values)
# print(len(surf_el_data.time))
# 
# # Plot the gap-free data to see if this worked - just look at one location 
# fig4, ax4 = plt.subplots(2,figsize=(22,18))
# #fig1.autofmt_xdate()
# 
# # Plot the original hycom data for this location
# ax4[0].scatter(surf_el_data.time[1402:2330].values, surf_el_data.surf_el[1402:2330,53,133].values, color='green', linewidth=3, label='Original')
# ax4[0].set_title('Original emp (kg/m2/s)', fontsize=30)
# ax4[0].set_xlabel('Time', fontsize=30)
# ax4[0].set_ylabel('emp (kg/m2/s)', fontsize=30)
# #ax4[0].legend()
# 
# # Plot the gap free data for this location
# ax4[1].scatter(surf_el_data_nogap.time[1444:2444].values, surf_el_data_nogap.surf_el[1444:2444,53,133].values, color='green', linewidth=3)
# ax4[1].set_title('Gap Free emp (kg/m2/s)', fontsize=30)
# ax4[1].set_xlabel('Time', fontsize=30)
# ax4[1].set_ylabel('emp (kg/m2/s)', fontsize=30)
# #fig4.set_facecolor('white')
# #fig4.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# surf_el_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_surf_el_nogaps.nc4')
# 
# # Remove all of this from memory to save RAM
# del(surf_el_data)
# del(surf_el_data_cp1)
# del(surf_el_data_nogap)
# =============================================================================


# =============================================================================
# # water_u = water u currents 
# # Load in the HYCOM data 
# water_u_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data//HYCOM_2020_water_u/HYCOM_2020_water_u.nc4')
# print(water_u_data)
# 
# # Get rid of duplicate times 
# water_u_data = water_u_data.sel(time=~water_u_data.get_index("time").duplicated())
# 
# # Look at the times in this dataset
# print(water_u_data.time[0].values)
# print(water_u_data.time[1].values)
# print(water_u_data.time[-1].values)
# 
# # Make a copy of the data to work with
# water_u_data_cp1 = water_u_data.copy()
# water_u_data_cp1['time'] = pd.to_datetime(water_u_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 1 hours since this is the frequency 
# # of the HYCOM data
# water_u_data_nogap = water_u_data_cp1.resample(time='3H').interpolate('linear') # hourly
# #print(qtot_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', water_u_data_nogap.time[0:5].values)
# print(water_u_data_nogap.time[-6:-1].values)
# print(water_u_data_nogap.time[-1].values)#
# print(len(water_u_data_nogap.time))
# print('\nOriginal data: ', water_u_data.time[0:5].values)
# print(water_u_data.time[-6:-1].values)
# print(water_u_data.time[-1].values)
# print(len(water_u_data.time))
# 
# #input('press enter to continue...')
# 
# # =============================================================================
# # # Plot the gap-free data to see if this worked - just look at one location 
# # fig4, ax4 = plt.subplots(2,figsize=(22,18))
# # #fig1.autofmt_xdate()
# # 
# # # Plot the original hycom data for this location
# # ax4[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax4[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax4[0].set_xlabel('Time', fontsize=30)
# # ax4[0].set_ylabel('water_u (m/s)', fontsize=30)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax4[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax4[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax4[1].set_xlabel('Time', fontsize=30)
# # ax4[1].set_ylabel('water_u (m/s)', fontsize=30)
# # #fig4.set_facecolor('white')
# # #fig4.savefig('hycom_gapfill_surtx01.png')
# # 
# # # Zoom in to the area with weird values
# # fig5, ax5 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax5[0].scatter(water_u_data.time[2150:2230].values, water_u_data.water_u[2150:2230,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax5[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax5[0].set_xlabel('Time', fontsize=30)
# # ax5[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax5[0].set_ylim(-1)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax5[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax5[1].set_title('Gap Free, Large Values Replaced water_u (m/s)', fontsize=30)
# # ax5[1].set_xlabel('Time', fontsize=30)
# # ax5[1].set_ylabel('water_u (m/s)', fontsize=30)
# # =============================================================================
# #fig4.set_facecolor('white')
# #fig4.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Attempt 1
# water_u_data_nogap2 = water_u_data_nogap.copy()
# water_u_data_nogap3 = water_u_data_nogap2.drop([np.datetime64('2020-09-29T00:00:00.000000000')], dim='time') #, '2020-09-29T03:00:00', 
# water_u_data_nogap4 = water_u_data_nogap3.drop([np.datetime64('2020-09-29T03:00:00.000000000')], dim='time')                                        # '2020-09-29T06:00:00', '2020-09-29T09:00:00',
# water_u_data_nogap5 = water_u_data_nogap4.drop([np.datetime64('2020-09-29T06:00:00.000000000')], dim='time')                                          #'2020-09-29T12:00:00', '2020-09-29T15:00:00',
# water_u_data_nogap6 = water_u_data_nogap5.drop([np.datetime64('2020-09-29T09:00:00.000000000')], dim='time')                                        #'2020-09-29T18:00:00', '2020-09-29T21:00:00')])
# water_u_data_nogap7 = water_u_data_nogap6.drop([np.datetime64('2020-09-29T12:00:00.000000000')], dim='time')
# water_u_data_nogap8 = water_u_data_nogap7.drop([np.datetime64('2020-09-29T15:00:00.000000000')], dim='time')
# water_u_data_nogap9 = water_u_data_nogap8.drop([np.datetime64('2020-09-29T18:00:00.000000000')], dim='time')
# water_u_data_nogap10 = water_u_data_nogap9.drop([np.datetime64('2020-09-29T21:00:00.000000000')], dim='time')
# 
# # Resampling for every 3 hours since this is the frequency 
# # of the HYCOM data
# water_u_data_nogap_nolarge = water_u_data_nogap10.resample(time='3H').interpolate('linear') # hourly
# 
# 
# # =============================================================================
# # # Plot this fix to see if it worked 
# # fig6, ax6 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax6[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax6[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax6[0].set_xlabel('Time', fontsize=30)
# # ax6[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[0].set_ylim(-0.4, 0.2)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax6[1].scatter(water_u_data_nogap_nolarge.time[1444:2444].values, water_u_data_nogap_nolarge.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax6[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax6[1].set_xlabel('Time', fontsize=30)
# # ax6[1].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[1].set_ylim(-0.4, 0.2)
# # =============================================================================
# 
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# # No gaps but large values
# #water_u_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_u/HYCOM_2020_water_u_nogaps.nc4')
# # No gaps and large values replaced
# water_u_data_nogap_nolarge.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_u/HYCOM_2020_water_u_nogaps_nolargevalues_02.nc4')
# 
# # Remove all of this from memory to save RAM
# del(water_u_data)
# del(water_u_data_cp1)
# del(water_u_data_nogap)
# =============================================================================



# =============================================================================
# # water_u = water u currents 
# # Load in the HYCOM data 
# water_v_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data//HYCOM_2020_water_v/HYCOM_2020_water_v.nc4')
# print(water_v_data)
# 
# # Get rid of duplicate times 
# water_v_data = water_v_data.sel(time=~water_v_data.get_index("time").duplicated())
# 
# # Look at the times in this dataset
# print(water_v_data.time[0].values)
# print(water_v_data.time[1].values)
# print(water_v_data.time[-1].values)
# 
# # Make a copy of the data to work with
# water_v_data_cp1 = water_v_data.copy()
# water_v_data_cp1['time'] = pd.to_datetime(water_v_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 1 hours since this is the frequency 
# # of the HYCOM data
# water_v_data_nogap = water_v_data_cp1.resample(time='3H').interpolate('linear') # hourly
# #print(qtot_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', water_v_data_nogap.time[0:5].values)
# print(water_v_data_nogap.time[-6:-1].values)
# print(water_v_data_nogap.time[-1].values)#
# print(len(water_v_data_nogap.time))
# print('\nOriginal data: ', water_v_data.time[0:5].values)
# print(water_v_data.time[-6:-1].values)
# print(water_v_data.time[-1].values)
# print(len(water_v_data.time))
# 
# #input('press enter to continue...')
# 
# # =============================================================================
# # # Plot the gap-free data to see if this worked - just look at one location 
# # fig4, ax4 = plt.subplots(2,figsize=(22,18))
# # #fig1.autofmt_xdate()
# # 
# # # Plot the original hycom data for this location
# # ax4[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax4[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax4[0].set_xlabel('Time', fontsize=30)
# # ax4[0].set_ylabel('water_u (m/s)', fontsize=30)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax4[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax4[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax4[1].set_xlabel('Time', fontsize=30)
# # ax4[1].set_ylabel('water_u (m/s)', fontsize=30)
# # #fig4.set_facecolor('white')
# # #fig4.savefig('hycom_gapfill_surtx01.png')
# # 
# # # Zoom in to the area with weird values
# # fig5, ax5 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax5[0].scatter(water_u_data.time[2150:2230].values, water_u_data.water_u[2150:2230,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax5[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax5[0].set_xlabel('Time', fontsize=30)
# # ax5[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax5[0].set_ylim(-1)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax5[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax5[1].set_title('Gap Free, Large Values Replaced water_u (m/s)', fontsize=30)
# # ax5[1].set_xlabel('Time', fontsize=30)
# # ax5[1].set_ylabel('water_u (m/s)', fontsize=30)
# # =============================================================================
# #fig4.set_facecolor('white')
# #fig4.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Attempt 1
# #water_v_data_nogap2 = water_v_data_nogap.copy()
# #water_v_data_nogap3 = water_v_data_nogap2.drop([np.datetime64('2020-09-29T00:00:00.000000000')], dim='time') #, '2020-09-29T03:00:00', 
# #water_v_data_nogap4 = water_v_data_nogap3.drop([np.datetime64('2020-09-29T03:00:00.000000000')], dim='time')                                        # '2020-09-29T06:00:00', '2020-09-29T09:00:00',
# #water_v_data_nogap5 = water_v_data_nogap4.drop([np.datetime64('2020-09-29T06:00:00.000000000')], dim='time')                                          #'2020-09-29T12:00:00', '2020-09-29T15:00:00',
# #water_v_data_nogap6 = water_v_data_nogap5.drop([np.datetime64('2020-09-29T09:00:00.000000000')], dim='time')                                        #'2020-09-29T18:00:00', '2020-09-29T21:00:00')])
# #water_v_data_nogap7 = water_v_data_nogap6.drop([np.datetime64('2020-09-29T12:00:00.000000000')], dim='time')
# #water_v_data_nogap8 = water_v_data_nogap7.drop([np.datetime64('2020-09-29T15:00:00.000000000')], dim='time')
# #water_v_data_nogap9 = water_v_data_nogap8.drop([np.datetime64('2020-09-29T18:00:00.000000000')], dim='time')
# #water_v_data_nogap10 = water_v_data_nogap9.drop([np.datetime64('2020-09-29T21:00:00.000000000')], dim='time')
# 
# # Resampling for every 3 hours since this is the frequency 
# # of the HYCOM data
# #water_v_data_nogap_nolarge = water_v_data_nogap10.resample(time='3H').interpolate('linear') # hourly
# 
# 
# # =============================================================================
# # # Plot this fix to see if it worked 
# # fig6, ax6 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax6[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax6[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax6[0].set_xlabel('Time', fontsize=30)
# # ax6[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[0].set_ylim(-0.4, 0.2)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax6[1].scatter(water_u_data_nogap_nolarge.time[1444:2444].values, water_u_data_nogap_nolarge.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax6[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax6[1].set_xlabel('Time', fontsize=30)
# # ax6[1].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[1].set_ylim(-0.4, 0.2)
# # =============================================================================
# 
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# # No gaps but large values
# water_v_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_v/HYCOM_2020_water_v_nogaps.nc4')
# # No gaps and large values replaced
# #water_v_data_nogap_nolarge.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_u/HYCOM_2020_water_u_nogaps_nolargevalues_02.nc4')
# 
# # Remove all of this from memory to save RAM
# del(water_v_data)
# del(water_v_data_cp1)
# del(water_v_data_nogap)
# =============================================================================



# =============================================================================
# # water_u = water u currents 
# # Load in the HYCOM data 
# water_temp_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data//HYCOM_2020_water_temp/HYCOM_2020_water_temp.nc4')
# print(water_temp_data)
# 
# # Get rid of duplicate times 
# water_temp_data = water_temp_data.sel(time=~water_temp_data.get_index("time").duplicated())
# 
# # Look at the times in this dataset
# print(water_temp_data.time[0].values)
# print(water_temp_data.time[1].values)
# print(water_temp_data.time[-1].values)
# 
# # Make a copy of the data to work with
# water_temp_data_cp1 = water_temp_data.copy()
# water_temp_data_cp1['time'] = pd.to_datetime(water_temp_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 1 hours since this is the frequency 
# # of the HYCOM data
# water_temp_data_nogap = water_temp_data_cp1.resample(time='3H').interpolate('linear') # hourly
# #print(qtot_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', water_temp_data_nogap.time[0:5].values)
# print(water_temp_data_nogap.time[-6:-1].values)
# print(water_temp_data_nogap.time[-1].values)#
# print(len(water_temp_data_nogap.time))
# print('\nOriginal data: ', water_temp_data.time[0:5].values)
# print(water_temp_data.time[-6:-1].values)
# print(water_temp_data.time[-1].values)
# print(len(water_temp_data.time))
# 
# #input('press enter to continue...')
# 
# # =============================================================================
# # # Plot the gap-free data to see if this worked - just look at one location 
# # fig4, ax4 = plt.subplots(2,figsize=(22,18))
# # #fig1.autofmt_xdate()
# # 
# # # Plot the original hycom data for this location
# # ax4[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax4[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax4[0].set_xlabel('Time', fontsize=30)
# # ax4[0].set_ylabel('water_u (m/s)', fontsize=30)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax4[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax4[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax4[1].set_xlabel('Time', fontsize=30)
# # ax4[1].set_ylabel('water_u (m/s)', fontsize=30)
# # #fig4.set_facecolor('white')
# # #fig4.savefig('hycom_gapfill_surtx01.png')
# # 
# # # Zoom in to the area with weird values
# # fig5, ax5 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax5[0].scatter(water_u_data.time[2150:2230].values, water_u_data.water_u[2150:2230,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax5[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax5[0].set_xlabel('Time', fontsize=30)
# # ax5[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax5[0].set_ylim(-1)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax5[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax5[1].set_title('Gap Free, Large Values Replaced water_u (m/s)', fontsize=30)
# # ax5[1].set_xlabel('Time', fontsize=30)
# # ax5[1].set_ylabel('water_u (m/s)', fontsize=30)
# # =============================================================================
# #fig4.set_facecolor('white')
# #fig4.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Attempt 1
# #water_v_data_nogap2 = water_v_data_nogap.copy()
# #water_v_data_nogap3 = water_v_data_nogap2.drop([np.datetime64('2020-09-29T00:00:00.000000000')], dim='time') #, '2020-09-29T03:00:00', 
# #water_v_data_nogap4 = water_v_data_nogap3.drop([np.datetime64('2020-09-29T03:00:00.000000000')], dim='time')                                        # '2020-09-29T06:00:00', '2020-09-29T09:00:00',
# #water_v_data_nogap5 = water_v_data_nogap4.drop([np.datetime64('2020-09-29T06:00:00.000000000')], dim='time')                                          #'2020-09-29T12:00:00', '2020-09-29T15:00:00',
# #water_v_data_nogap6 = water_v_data_nogap5.drop([np.datetime64('2020-09-29T09:00:00.000000000')], dim='time')                                        #'2020-09-29T18:00:00', '2020-09-29T21:00:00')])
# #water_v_data_nogap7 = water_v_data_nogap6.drop([np.datetime64('2020-09-29T12:00:00.000000000')], dim='time')
# #water_v_data_nogap8 = water_v_data_nogap7.drop([np.datetime64('2020-09-29T15:00:00.000000000')], dim='time')
# #water_v_data_nogap9 = water_v_data_nogap8.drop([np.datetime64('2020-09-29T18:00:00.000000000')], dim='time')
# #water_v_data_nogap10 = water_v_data_nogap9.drop([np.datetime64('2020-09-29T21:00:00.000000000')], dim='time')
# 
# # Resampling for every 3 hours since this is the frequency 
# # of the HYCOM data
# #water_v_data_nogap_nolarge = water_v_data_nogap10.resample(time='3H').interpolate('linear') # hourly
# 
# 
# # =============================================================================
# # # Plot this fix to see if it worked 
# # fig6, ax6 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax6[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax6[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax6[0].set_xlabel('Time', fontsize=30)
# # ax6[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[0].set_ylim(-0.4, 0.2)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax6[1].scatter(water_u_data_nogap_nolarge.time[1444:2444].values, water_u_data_nogap_nolarge.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax6[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax6[1].set_xlabel('Time', fontsize=30)
# # ax6[1].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[1].set_ylim(-0.4, 0.2)
# # =============================================================================
# 
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# # No gaps but large values
# water_temp_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_temp/HYCOM_2020_water_temp_nogaps.nc4')
# # No gaps and large values replaced
# #water_v_data_nogap_nolarge.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_u/HYCOM_2020_water_u_nogaps_nolargevalues_02.nc4')
# 
# # Remove all of this from memory to save RAM
# del(water_temp_data)
# del(water_temp_data_cp1)
# del(water_temp_data_nogap)
# =============================================================================



# =============================================================================
# # water_u = water u currents 
# # Load in the HYCOM data 
# salinity_data = xr.open_mfdataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data//HYCOM_2020_salinity/HYCOM_2020_salinity.nc4')
# print(salinity_data)
# 
# # Get rid of duplicate times 
# salinity_data = salinity_data.sel(time=~salinity_data.get_index("time").duplicated())
# 
# # Look at the times in this dataset
# print(salinity_data.time[0].values)
# print(salinity_data.time[1].values)
# print(salinity_data.time[-1].values)
# 
# # Make a copy of the data to work with
# salinity_data_cp1 = salinity_data.copy()
# salinity_data_cp1['time'] = pd.to_datetime(salinity_data_cp1['time'])
# 
# # Resample and interpolate the data over time
# # Resampling for every 1 hours since this is the frequency 
# # of the HYCOM data
# salinity_data_nogap = salinity_data_cp1.resample(time='3H').interpolate('linear') # hourly
# #print(qtot_data_nogap)
# 
# # Compare the times in this new dataset with the times in the original one
# # New data 
# print('New gap free data: ', salinity_data_nogap.time[0:5].values)
# print(salinity_data_nogap.time[-6:-1].values)
# print(salinity_data_nogap.time[-1].values)#
# print(len(salinity_data_nogap.time))
# print('\nOriginal data: ', salinity_data.time[0:5].values)
# print(salinity_data.time[-6:-1].values)
# print(salinity_data.time[-1].values)
# print(len(salinity_data.time))
# 
# #input('press enter to continue...')
# 
# # =============================================================================
# # # Plot the gap-free data to see if this worked - just look at one location 
# # fig4, ax4 = plt.subplots(2,figsize=(22,18))
# # #fig1.autofmt_xdate()
# # 
# # # Plot the original hycom data for this location
# # ax4[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax4[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax4[0].set_xlabel('Time', fontsize=30)
# # ax4[0].set_ylabel('water_u (m/s)', fontsize=30)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax4[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax4[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax4[1].set_xlabel('Time', fontsize=30)
# # ax4[1].set_ylabel('water_u (m/s)', fontsize=30)
# # #fig4.set_facecolor('white')
# # #fig4.savefig('hycom_gapfill_surtx01.png')
# # 
# # # Zoom in to the area with weird values
# # fig5, ax5 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax5[0].scatter(water_u_data.time[2150:2230].values, water_u_data.water_u[2150:2230,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax5[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax5[0].set_xlabel('Time', fontsize=30)
# # ax5[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax5[0].set_ylim(-1)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax5[1].scatter(water_u_data_nogap.time[1444:2444].values, water_u_data_nogap.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax5[1].set_title('Gap Free, Large Values Replaced water_u (m/s)', fontsize=30)
# # ax5[1].set_xlabel('Time', fontsize=30)
# # ax5[1].set_ylabel('water_u (m/s)', fontsize=30)
# # =============================================================================
# #fig4.set_facecolor('white')
# #fig4.savefig('hycom_gapfill_surtx01.png')
# 
# 
# # Attempt 1
# #water_v_data_nogap2 = water_v_data_nogap.copy()
# #water_v_data_nogap3 = water_v_data_nogap2.drop([np.datetime64('2020-09-29T00:00:00.000000000')], dim='time') #, '2020-09-29T03:00:00', 
# #water_v_data_nogap4 = water_v_data_nogap3.drop([np.datetime64('2020-09-29T03:00:00.000000000')], dim='time')                                        # '2020-09-29T06:00:00', '2020-09-29T09:00:00',
# #water_v_data_nogap5 = water_v_data_nogap4.drop([np.datetime64('2020-09-29T06:00:00.000000000')], dim='time')                                          #'2020-09-29T12:00:00', '2020-09-29T15:00:00',
# #water_v_data_nogap6 = water_v_data_nogap5.drop([np.datetime64('2020-09-29T09:00:00.000000000')], dim='time')                                        #'2020-09-29T18:00:00', '2020-09-29T21:00:00')])
# #water_v_data_nogap7 = water_v_data_nogap6.drop([np.datetime64('2020-09-29T12:00:00.000000000')], dim='time')
# #water_v_data_nogap8 = water_v_data_nogap7.drop([np.datetime64('2020-09-29T15:00:00.000000000')], dim='time')
# #water_v_data_nogap9 = water_v_data_nogap8.drop([np.datetime64('2020-09-29T18:00:00.000000000')], dim='time')
# #water_v_data_nogap10 = water_v_data_nogap9.drop([np.datetime64('2020-09-29T21:00:00.000000000')], dim='time')
# 
# # Resampling for every 3 hours since this is the frequency 
# # of the HYCOM data
# #water_v_data_nogap_nolarge = water_v_data_nogap10.resample(time='3H').interpolate('linear') # hourly
# 
# 
# # =============================================================================
# # # Plot this fix to see if it worked 
# # fig6, ax6 = plt.subplots(2,figsize=(22,18))
# # 
# # # Plot the original hycom data for this location
# # ax6[0].scatter(water_u_data.time[1402:2330].values, water_u_data.water_u[1402:2330,0,100,100].values, color='green', linewidth=3, label='Original')
# # ax6[0].set_title('Original water_u (m/s)', fontsize=30)
# # ax6[0].set_xlabel('Time', fontsize=30)
# # ax6[0].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[0].set_ylim(-0.4, 0.2)
# # #ax4[0].legend()
# # 
# # # Plot the gap free data for this location
# # ax6[1].scatter(water_u_data_nogap_nolarge.time[1444:2444].values, water_u_data_nogap_nolarge.water_u[1444:2444,0,100,100].values, color='green', linewidth=3)
# # ax6[1].set_title('Gap Free water_u (m/s)', fontsize=30)
# # ax6[1].set_xlabel('Time', fontsize=30)
# # ax6[1].set_ylabel('water_u (m/s)', fontsize=30)
# # ax6[1].set_ylim(-0.4, 0.2)
# # =============================================================================
# 
# 
# 
# # Save this new dataset as the new HYCOM zeta data
# # No gaps but large values
# salinity_data_nogap.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_salinity/HYCOM_2020_salinity_nogaps.nc4')
# # No gaps and large values replaced
# #water_v_data_nogap_nolarge.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/HYCOM_data/HYCOM_2020_water_u/HYCOM_2020_water_u_nogaps_nolargevalues_02.nc4')
# 
# # Remove all of this from memory to save RAM
# del(salinity_data)
# del(salinity_data_cp1)
# del(salinity_data_nogap)
# =============================================================================


