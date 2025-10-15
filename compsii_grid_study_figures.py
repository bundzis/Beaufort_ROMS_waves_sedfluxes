#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 19:57:50 2023

@author: brun1463
"""

########################### COMPS II Grid Figures ##############################
# The purpose of this script is to make grid figures that also serve as 
# study site figures for my COMPS II paper. There will be two subplots, one
# with shelf batyhemtry from 0 - 200 m with 10 m gray contours, and all
# rivers labeled. The other wll habe the grid cells labeled with some 
# major towns labeled.
#
# Notes:
# - This script needs to be run in the geopand env
#
################################################################################


# Load in the packages 
import numpy as np
import xarray as xr
import pandas as pd
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
import matplotlib
#from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as patches
#from matplotlib import transforms 
import cmocean
#import matplotlib.ticker as tick
#import matplotlib.patches as patches
#from matplotlib.colors import LinearSegmentedColormap
from shapely.geometry import Polygon
#import seaborn as sns
#import missingno as msno
#import openpyxl
#import math
import os
import wget
import geopandas as gpd
import matplotlib.colors as mcolors
from matplotlib import colors

# set a universal fontsize 
fontsize = 20 # 32 #20

# Set the tick size for all plots
matplotlib.rc('xtick', labelsize=fontsize) 
matplotlib.rc('ytick', labelsize=fontsize)

# Prevent tick labels from overlapping
matplotlib.rcParams['xtick.major.pad'] = 12 # 28 # 12
matplotlib.rcParams['ytick.major.pad'] = 12 # 18 # 12

# Load in the model grid
grid = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Grids/KakAKgrd_shelf_big010_smooth006.nc') 

# Load in some output
#output = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/model_output/Full_run_0001/ocean_his_biggrid010_gridwindsiniwaves_rivs_si_smooth006_nobulk_chaflaradnudclm_dbsed0004_0001.nc')
#output = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Output/dbsed0003/ocean_his_beaufort_shelf_2020_dbsed0003_0001.nc')


# Pull out specific transects 
# Find the coordinates in eta and xi for Flaxman Island and Harrison Bay
# Flaxman Island (70.11745, -145.95978)
print('Flaxman: ', np.where(np.round(grid.lon_rho.values, 2) == -145.96))
flax_eta = slice(28, 170)
flax_xi = 388
flax_idx = np.where(np.round(grid.lon_rho.values, 2) == -145.96)

# Harrison Bay (70.41205, -151.50262)
print('Harrison: ', np.where(np.round(grid.lon_rho.values, 2) == -151.50))
har_eta = slice(26, 170) # (26, 170)
har_xi = 105
har_idx = np.where(np.round(grid.lon_rho.values, 2) == -151.50)

# Transect 1
eta1 = slice(21, 160)
xi1 = 320

# Transect 2
eta2 = slice(45, 153)
xi2 = 500

# CODA moorings 
# Mooring 2
eta_rho_moor1 = 103
xi_rho_moor1 = 202
# Mooring 3
eta_rho_moor2 = 73
xi_rho_moor2 = 390

# Point 1 
eta_point1 = 30
xi_point1 = 120

# Point 2
eta_point2 = 70
xi_point2 = 510

# -------------------------------------------------
# ------ Plot 1: Ocean Bathymetry and Grid Cells ------
# -------------------------------------------------
# Make a plot of the ocean bathymetry from 5 - 200 m
# and have every 16th grid line drawn, as well as the good
# land mask

# Make a plot of bathymetry with grid lines
# depths for the colorbar
fig1, ax1 = plt.subplots(figsize=(25,12)) # (15, 12) for horizontal colorbar; (20, 12) for vertical colorbar

# Determine spacing of contours
lev1 = np.arange(grid.h.values.min()-1, 200, 1)
cmap1 = cmocean.cm.deep
cmap1.set_under('darkgray')

lev2 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Plot the bathymetry
cs1 = ax1.contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev1, cmap=cmap1, extend='both')

# Plot a few key towns
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517
sc1 = ax1.scatter(grid.lon_rho[eta_kakt_idx, xi_kakt_idx].values, grid.lat_rho[eta_kakt_idx, xi_kakt_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kaktovik')

# Prudhoe Bay
eta_prud_idx = 37
xi_prud_idx = 279
sc2 = ax1.scatter(grid.lon_rho[eta_prud_idx, xi_prud_idx].values, grid.lat_rho[eta_prud_idx, xi_prud_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Prudhoe Bay')

# Plot grid lines
# This puts the vertical lines through the longitudes over all latitudes 
ax1.plot(grid.lon_rho[:, ::16].values, grid.lat_rho[:, ::16].values, '-', color='k')

# This plots the vertical line on the right boundary of the grid
ax1.plot(grid.lon_rho[:, -1].values, grid.lat_rho[:, -1].values, '-', color='k')

# This plots the upper horizontal boundary of the grid
ax1.plot(grid.lon_rho[-1, :].values, grid.lat_rho[-1, :].values, '-', color='k')

for i in range(12):
    # This plots the horizontal lines (including the bottom boundary)
    ax1.plot(grid.lon_rho[(i*16), :].values, grid.lat_rho[(i*16), :].values, '-', color='k')


# Add title
ax1.set_title('The Grid Cells (every 16th) in the Beaufort Sea', fontsize=fontsize, y=1.08)
 
# Format axes
ax1.set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1.set_ylabel('Latitude (degrees)', fontsize=fontsize)
#plt.grid(True)

ax1.legend(fontsize=fontsize)
 
# Specify colorbar
#cbar6 = fig6.colorbar(cs9, label='Depth (meters)', orientation='horizontal')
cbar1 = fig1.colorbar(cs1, label='Depth (m)', orientation='vertical')
cbar1.set_label(label='Depth (meters)', size=fontsize)



# -------------------------------------------------
# ------ Plot 2: Ocean Bathymetry and Rivers ------
# -------------------------------------------------
# Make ansimilar plot to above but accentuate shelf bathymetry more and
# label all of the rivers in the model

# Make a plot of bathymetry with grid lines
# depths for the colorbar
fig2, ax2 = plt.subplots(figsize=(25,12)) # (15, 12) for horizontal colorbar; (20, 12) for vertical colorbar

# Determine spacing of contours
lev2 = np.arange(grid.h.values.min()-1, 200, 1)
cmap2 = cmocean.cm.deep
cmap2.set_under('darkgray')

lev3 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Plot the bathymetry
cs2 = ax2.contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev2, cmap=cmap1, extend='both')
cs3 = ax2.contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev3, colors='dimgray')

# Plot the 14 rivers in the grid
# Go from West to East
# Kalikpik River
eta_kal_idx = 23 #22
xi_kal_idx = 87
s1 = ax2.scatter(grid.lon_rho[eta_kal_idx, xi_kal_idx].values, grid.lat_rho[eta_kal_idx, xi_kal_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
eta_fis_idx = 20
xi_fis_idx = 117 #116
s2 = ax2.scatter(grid.lon_rho[eta_fis_idx, xi_fis_idx].values, grid.lat_rho[eta_fis_idx, xi_fis_idx].values, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
eta_col_idx = 39
xi_col_idx = 166 #166
s3 = ax2.scatter(grid.lon_rho[eta_col_idx, xi_col_idx].values, grid.lat_rho[eta_col_idx, xi_col_idx].values, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
eta_sak_idx = 46 #45
xi_sak_idx = 234
s4 = ax2.scatter(grid.lon_rho[eta_sak_idx, xi_sak_idx].values, grid.lat_rho[eta_sak_idx, xi_sak_idx].values, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk
eta_kuk_idx = 41 #40
xi_kuk_idx = 239
s5 = ax2.scatter(grid.lon_rho[eta_kuk_idx, xi_kuk_idx].values, grid.lat_rho[eta_kuk_idx, xi_kuk_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Kukpuk')

# Kuparuk
eta_kup_idx = 41 #40
xi_kup_idx = 242
s6 = ax2.scatter(grid.lon_rho[eta_kup_idx, xi_kup_idx].values, grid.lat_rho[eta_kup_idx, xi_kup_idx].values, 
            marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#eta_faw_idx = 44 #43
#xi_faw_idx = 249
#s7 = ax2.scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
eta_put_idx = 28 #27
xi_put_idx = 264
s8 = ax2.scatter(grid.lon_rho[eta_put_idx, xi_put_idx].values, grid.lat_rho[eta_put_idx, xi_put_idx].values, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
eta_sag_idx = 37 #36
xi_sag_idx = 279
s9 = ax2.scatter(grid.lon_rho[eta_sak_idx, xi_sak_idx].values, grid.lat_rho[eta_sak_idx, xi_sak_idx].values, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
eta_sta_idx = 27 #26
xi_sta_idx = 393
s10 = ax2.scatter(grid.lon_rho[eta_sta_idx, xi_sta_idx].values, grid.lat_rho[eta_sta_idx, xi_sta_idx].values, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
eta_can_idx = 20 #19
xi_can_idx = 416
s11 = ax2.scatter(grid.lon_rho[eta_can_idx, xi_can_idx].values, grid.lat_rho[eta_can_idx, xi_can_idx].values, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
eta_kat_idx = 9 #8
xi_kat_idx = 447
s12 = ax2.scatter(grid.lon_rho[eta_kat_idx, xi_kat_idx].values, grid.lat_rho[eta_kat_idx, xi_kat_idx].values, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
eta_hul_idx = 40
xi_hul_idx = 489
s13 = ax2.scatter(grid.lon_rho[eta_hul_idx, xi_hul_idx].values, grid.lat_rho[eta_hul_idx, xi_hul_idx].values, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
eta_jag_idx = 62 #61
xi_jag_idx = 528
s14 = ax2.scatter(grid.lon_rho[eta_jag_idx, xi_jag_idx].values, grid.lat_rho[eta_jag_idx, xi_jag_idx].values, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
eta_sik_idx = 46
xi_sik_idx = 574 #573
s15 = ax2.scatter(grid.lon_rho[eta_sik_idx, xi_sik_idx].values, grid.lat_rho[eta_sik_idx, xi_sik_idx].values, 
            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')

#ax2.legend(fontsize=20, ncol=4, columnspacing=0.3, loc='upper right')

lines=[s1,s2,s3,s4,s5,s6,s8,s9,s10,s11,s12,s13,s14,s15]
include=[0,1,2,3,4,5,6]
legend1 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc='lower left', fontsize=fontsize, ncol=2,
                     columnspacing=0.3)
legend2 = plt.legend([lines[i] for i in [7,8,9,10,11,12,13]],['Sagavanirktok', 'Staines', 'Canning','Katakturuk', 'Hulahula', 'Jago', 
                                                               'Siksik'], loc='upper right', fontsize=fontsize, ncol=2, columnspacing=0.3)
plt.gca().add_artist(legend1)
#l1 = ax2.legend(['Kalikpik', 'Fish Creek', 'Colville', 'Sakonowyak', 'Kukpuk', 'Kuparuk', 'Fawn Creek'], fontsize=20, ncol=2, columnspacing=0.3, loc='lower left')
#l2 = ax2.legend([s8, s9, s10, s11, s12, s13, s14], fontsize=20, ncol=2, columnspacing=0.3, loc='upper right')
#gca().add_artist(l1)
#plt.gca().add_artist(l2)


# add title
ax2.set_title('The Rivers in the Beaufort Sea Model Domain', fontsize=fontsize, y=1.08)
 
# format axes
ax2.set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax2.set_ylabel('Latitude (degrees)', fontsize=fontsize)
#plt.grid(True)
 
# specify colorbar
#cbar6 = fig6.colorbar(cs9, label='Depth (meters)', orientation='horizontal')
cbar2 = fig2.colorbar(cs2, label='Depth (m)', orientation='vertical')
cbar2.set_label(label='Depth (meters)', size=fontsize)


# -----------------------------------------------------------
# ------ Plot 3: Subplot of two figures from above ------
# -------------------------------------------------------------
# Make a subplot with the two figures from above together as 
# one figure (for the paper). Remove the title and have them share
# the colorbar
# (15,30)
# Make the figure
fig3, ax3 = plt.subplots(2, figsize=(17,16)) # (width, height), (15, 12) for horizontal colorbar; (20, 16) for vertical colorbaroutside plots

# Determine spacing of contours
#lev3 = np.arange(grid.h.values.min()-1, 200, 1)
lev3 = np.arange(grid.h.values.min()-1, 100, 1)
cmap3 = cmocean.cm.deep
cmap3.set_under('darkgray')

lev4 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# --- Top subplot - Grid ---
# Plot the bathymetry
cs6 = ax3[0].contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev3, cmap=cmap3, extend='both')

# Plot a few key towns
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517
sc3 = ax3[0].scatter(grid.lon_rho[eta_kakt_idx, xi_kakt_idx].values, grid.lat_rho[eta_kakt_idx, xi_kakt_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kaktovik Weather Station')

# Prudhoe Bay
eta_prud_idx = 37
xi_prud_idx = 279
sc4 = ax3[0].scatter(grid.lon_rho[eta_prud_idx, xi_prud_idx].values, grid.lat_rho[eta_prud_idx, xi_prud_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Prudhoe Bay Weather Station')

# Plot the CODA moorings
# Mooring 2
sc1 = ax3[0].scatter(grid.lon_rho[eta_rho_moor1, xi_rho_moor1].values, grid.lat_rho[eta_rho_moor1, xi_rho_moor1].values, 
            marker='x', s=300, linewidth=4, color='orange', label='CODA Mooring 2')

# Mooring 3
sc2 = ax3[0].scatter(grid.lon_rho[eta_rho_moor2, xi_rho_moor2].values, grid.lat_rho[eta_rho_moor2, xi_rho_moor2].values, 
            marker='x', s=300, linewidth=4, color='purple', label='CODA Mooring 3')

# Point 1 
sc5 = ax3[0].scatter(grid.lon_rho[eta_point1, xi_point1].values, grid.lat_rho[eta_point1, xi_point1].values, 
            marker='+', s=300, linewidth=4, color='darkgoldenrod', label='Harrison Bay')

# Point 2
#sc6 = ax3[0].scatter(grid.lon_rho[eta_point2, xi_point2].values, grid.lat_rho[eta_point2, xi_point2].values, 
 #           marker='+', s=300, linewidth=4, color='lightcoral', label='Point 2')

# Plot grid lines
# This puts the vertical lines through the longitudes over all latitudes 
ax3[0].plot(grid.lon_rho[:, ::16].values, grid.lat_rho[:, ::16].values, '-', color='k')

# This plots the vertical line on the right boundary of the grid
ax3[0].plot(grid.lon_rho[:, -1].values, grid.lat_rho[:, -1].values, '-', color='k')

# This plots the upper horizontal boundary of the grid
ax3[0].plot(grid.lon_rho[-1, :].values, grid.lat_rho[-1, :].values, '-', color='k')

for i in range(12):
    # This plots the horizontal lines (including the bottom boundary)
    ax3[0].plot(grid.lon_rho[(i*16), :].values, grid.lat_rho[(i*16), :].values, '-', color='k')

# Add title
#ax1.set_title('The Grid Cells (every 16th) in the Beaufort Sea', fontsize=30, y=1.08)
 
# Format axes
#ax3[0].set_xlabel('Longitude \n(degrees)', fontsize=30)
ax3[0].set_ylabel('Latitude \n(degrees)', fontsize=fontsize, rotation=0, labelpad=60, va='center') # 118
#plt.grid(True)

ax3[0].legend(fontsize=14, labelspacing=0.2, ncol=2, loc='lower left')

# Add a subset of Alaska
# Old - external vertical colorbar
#left, bottom, width, height = [0.62, 0.75, 0.18, 0.14]
#ax2 = fig3.add_axes([left, bottom, width, height])
# New - internal horizontal colorbar
left, bottom, width, height = [0.77, 0.75, 0.18, 0.14]
ax2 = fig3.add_axes([left, bottom, width, height])


# --- Attempt at getting plot of Alaska ---
# download the data excel file from the USDA website to local folder
#filename = wget.download("https://www.ers.usda.gov/media/rbmpu1zi/mapdata2021.xlsx")
# load the excel file into a pandas dataframe & skip header rows
#df = pd.read_excel(os.getcwd()+'/mapdata2021.xlsx',skiprows=4)
df = pd.read_excel('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/mapdata2021.xlsx',skiprows=4) # TMP
# rename the columns of interest
df = df.rename(columns={'Unnamed: 0':'state',
                       'Percent':'pct_food_insecure'})
# retain only the columns of interest
df = df[['state','pct_food_insecure']]

# take a look at the dataframe for missing values
#msno.matrix(df)

df = df[df.state.str.len()==2]
df.pct_food_insecure.hist()

#wget.download("https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip")
#gdf = gpd.read_file(os.getcwd()+'/cb_2018_us_state_500k')
gdf = gpd.read_file('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/cb_2018_us_state_500k') # TMP
#gdf.head()
gdf = gdf.merge(df,left_on='STUSPS',right_on='state')
#gdf.plot()

world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))

# =============================================================================
# # NOTE: the convention for polygon points is (Long, Lat)....counterintuitive
# polygon = Polygon([(-175,50),(-175,72),(-140, 72),(-140,50)])
# # polygon = Polygon([(-180,0),(-180,90),(-120,90),(-120,0)])
# 
# # polygon=hipolygon
# poly_gdf = gpd.GeoDataFrame( geometry=[polygon], crs=world.crs)
# 
# fig, ax1 = plt.subplots(1, figsize=(8, 18))
# world.plot(ax=ax1)
# poly_gdf.boundary.plot(ax = ax1, color="red")
# ax1.set_title("The red polygon can be used to clip Alaska's western islands", fontsize=20)
# ax1.set_axis_off()
# plt.show()
# =============================================================================

# Add Alaska inset (old - vertical external colorbar)
# Create the polygon outline for Alaska 
#polygon = Polygon([(-170,50),(-170,72),(-140, 72),(-140,50)])
# Plot this Alaska polygon
#world.clip(polygon).plot(color='lightblue', linewidth=0.8, edgecolor='0.8', ax=ax2) 
# Create a Rectangle patch to place around grid site
#rect = patches.Rectangle((-155, 69.8), 14.5, 2.5, linewidth=3, angle=-4, edgecolor='r', facecolor='none')
# Add the patch to the Axes
#ax2.add_patch(rect)
# Add letters to outline
#plt.text(0.70, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # old, vertical external colorbar


# Add Alaska inset (new - horizontal internal colorbar)
# Create the polygon outline for Alaska 
polygon = Polygon([(-170,50),(-170,72),(-140, 72),(-140,50)])
# Plot this Alaska polygon
world.clip(polygon).plot(color='lightblue', linewidth=0.8, edgecolor='0.8', ax=ax2) 
# Create a Rectangle patch to place around grid site
rect = patches.Rectangle((-155, 69.8), 14.5, 2.5, linewidth=3, angle=-4, edgecolor='r', facecolor='none')
# Add the patch to the Axes
ax2.add_patch(rect)
# Add letters to outline
#plt.text(0.70, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # old, vertical external colorbar
plt.text(0.86, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # new, horizontal internal colorbar


# Adjust plot labels and such
ax2.grid(False)
ax2.set_xlim(-170, -135)
ax2.set_ylim(53, 73)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.tick_params(left = False, bottom = False)

# =============================================================================
# # Add in a scale bar 
# import matplotlib_scalebar 
# from matplotlib_scalebar.scalebar import ScaleBar
# from sklearn.metrics.pairwise import haversine_distances
# lat1 = 70.75
# lon1 = -152
# lon2 = -151
# A = [lon1*np.pi/180., lat1*np.pi/180.]
# B = [lon2*np.pi/180., lat1*np.pi/180.]
# scale_dist = (6371000)*haversine_distances([A,B])[0,1]
# ax3[0].add_artist(ScaleBar(dx=scale_dist, units='m', font_properties={'size':'large'}, 
#                            location='lower center'))
# =============================================================================


# --- end of Alaska attempt ---



# --- Bottom subplot - Rivers ---
# Plot the bathymetry
cs4 = ax3[1].contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev3, cmap=cmap3, extend='both')
cs5 = ax3[1].contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev4, colors='dimgray')

# =============================================================================
# # Plot the Flaxman Island and Harrison Bay transects
# # Flaxman Island 
# p1 = ax3[1].plot(grid.lon_rho[flax_eta, flax_xi].values, grid.lat_rho[flax_eta, flax_xi].values, color='red', 
#          linewidth=5, label='Flaxman Island')
# 
# # Harrison Bay
# p2 = ax3[1].plot(grid.lon_rho[har_eta, har_xi].values, grid.lat_rho[har_eta, har_xi].values, color='gold', 
#          linewidth=5, label='Harrison Bay')
# 
# # Plot the other two transects used 
# # Transect 1
# p3 = ax3[1].plot(grid.lon_rho[eta1, xi1].values, grid.lat_rho[eta1, xi1].values, color='purple', 
#          linewidth=5, label='Transect 1')
# 
# # Transect 2
# p4 = ax3[1].plot(grid.lon_rho[eta2, xi2].values, grid.lat_rho[eta2, xi2].values, color='navy', 
#          linewidth=5, label='Transect 2')
# =============================================================================



# Plot the 14 rivers in the grid
# Go from West to East
# Kalikpik River
eta_kal_idx = 23 #22
xi_kal_idx = 87
s1 = ax3[1].scatter(grid.lon_rho[eta_kal_idx, xi_kal_idx].values, grid.lat_rho[eta_kal_idx, xi_kal_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
eta_fis_idx = 20
xi_fis_idx = 117 #116
s2 = ax3[1].scatter(grid.lon_rho[eta_fis_idx, xi_fis_idx].values, grid.lat_rho[eta_fis_idx, xi_fis_idx].values, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
eta_col_idx = 39
xi_col_idx = 166 #166
s3 = ax3[1].scatter(grid.lon_rho[eta_col_idx, xi_col_idx].values, grid.lat_rho[eta_col_idx, xi_col_idx].values, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
eta_sak_idx = 46 #45
xi_sak_idx = 234
s4 = ax3[1].scatter(grid.lon_rho[eta_sak_idx, xi_sak_idx].values, grid.lat_rho[eta_sak_idx, xi_sak_idx].values, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s5 = ax3[1].scatter(grid.lon_rho[eta_kup_idx, xi_kup_idx].values, grid.lat_rho[eta_kup_idx, xi_kup_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax3[1].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#eta_faw_idx = 44 #43
#xi_faw_idx = 249
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
eta_put_idx = 28 #27
xi_put_idx = 264
s8 = ax3[1].scatter(grid.lon_rho[eta_put_idx, xi_put_idx].values, grid.lat_rho[eta_put_idx, xi_put_idx].values, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
eta_sag_idx = 37 #36
xi_sag_idx = 279
s9 = ax3[1].scatter(grid.lon_rho[eta_sag_idx, xi_sag_idx].values, grid.lat_rho[eta_sag_idx, xi_sag_idx].values, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
eta_sta_idx = 27 #26
xi_sta_idx = 393
s10 = ax3[1].scatter(grid.lon_rho[eta_sta_idx, xi_sta_idx].values, grid.lat_rho[eta_sta_idx, xi_sta_idx].values, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
eta_can_idx = 20 #19
xi_can_idx = 416
s11 = ax3[1].scatter(grid.lon_rho[eta_can_idx, xi_can_idx].values, grid.lat_rho[eta_can_idx, xi_can_idx].values, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
eta_kat_idx = 9 #8
xi_kat_idx = 447
s12 = ax3[1].scatter(grid.lon_rho[eta_kat_idx, xi_kat_idx].values, grid.lat_rho[eta_kat_idx, xi_kat_idx].values, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
eta_hul_idx = 40
xi_hul_idx = 489
s13 = ax3[1].scatter(grid.lon_rho[eta_hul_idx, xi_hul_idx].values, grid.lat_rho[eta_hul_idx, xi_hul_idx].values, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
eta_jag_idx = 62 #61
xi_jag_idx = 528
s14 = ax3[1].scatter(grid.lon_rho[eta_jag_idx, xi_jag_idx].values, grid.lat_rho[eta_jag_idx, xi_jag_idx].values, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
eta_sik_idx = 46
xi_sik_idx = 574 #573
s15 = ax3[1].scatter(grid.lon_rho[eta_sik_idx, xi_sik_idx].values, grid.lat_rho[eta_sik_idx, xi_sik_idx].values, 
            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')

#ax2.legend(fontsize=20, ncol=4, columnspacing=0.3, loc='upper right')

# Old but works
# =============================================================================
# lines=[s1,s2,s3,s4,s5,s6,s8,s9,s10,s11,s12,s13,s14,s15]
# include=[0,1,2,3,4,5,6]
# legend3 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc='lower left', fontsize=14, ncol=3,
#                      columnspacing=0.1, labelspacing=0.2)
# legend4 = ax3[1].legend([lines[i] for i in [7,8,9,10,11,12,13]],['Sagavanirktok', 'Staines', 'Canning','Katakturuk', 'Hulahula', 'Jago', 
#                                                                'Siksik'], loc='upper right', fontsize=14, ncol=3, columnspacing=0.1, 
#                         labelspacing=0.2)
# plt.gca().add_artist(legend3)
# =============================================================================
#gca().add_artist(legend3)

# New attempt to include transectsin legend
# =============================================================================
# lines=[s1,s2,s3,s4,s5,s6,s8,s9,s10,s11,s12,s13,s14,s15]
# include=[0,1,2,3,4,5,6]
# #legend3 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc='lower left', fontsize=14, ncol=3,
#  #                    columnspacing=0.1, labelspacing=0.2)
# legend4 = ax3[1].legend([lines[i] for i in [7,8,9,10,11,12,13]],['Sagavanirktok', 'Staines', 'Canning','Katakturuk', 'Hulahula', 'Jago', 
#                                                                'Siksik'], loc='upper right', fontsize=14, ncol=3, columnspacing=0.1, 
#                         labelspacing=0.2)
# curves = [p1,p2,p3,p4]
# include2 = [0,1,2,3]
# #legend5 = ax3[1].legend()
# legend5 = plt.legend([curves[i] for i in [0,1,2,3]], ['Harrison Bay', 'Transect 1', 'Flaxman Island', 'Transect 2'], 
#                         loc='lower right', fontsize=14, ncol=2, columnspacing=0.1, labelspacing=0.2)
# plt.gca().add_artist(legend5)
# #plt.gca().add_artist(legend5)
# =============================================================================

# Temporary not too bad legend
ax3[1].legend(loc='lower left', ncol=4, fontsize=14, columnspacing=0.1, labelspacing=0.1)

# Format axes
ax3[1].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax3[1].set_ylabel('Latitude \n(degrees)', fontsize=fontsize, rotation=0, labelpad=60, va='center') #118

# Specify colorbar
#cbar3 = fig3.colorbar(cs6, label='Depth (m)', orientation='vertical', ax=ax3)
#cbar3.set_label(label='Depth (meters)', size=fontsize)
# If horizontal and in axes
cbar3 = plt.colorbar(cs6, cax=ax3[1].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 20, 40, 60, 80], ax=ax3[1], 
                     orientation='horizontal').set_label(label='Depth (m)', size=fontsize-2)

# Label subplots
# put labels in top right corners
#plt.text(0.715, 0.853, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
#plt.text(0.715, 0.442, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
# put labels in bottom right corners
plt.text(0.870, 0.553, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.870, 0.142, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)

#plt.tight_layout()




# -------------------------------------------------
# ------ Plot 4: Ocean Bathymetry  ------
# -------------------------------------------------
# Make ansimilar plot to above but accentuate shelf bathymetry more and
# label all of the rivers in the model

# Make a plot of bathymetry with grid lines
# depths for the colorbar
fig4, ax4 = plt.subplots(figsize=(25,12)) # (15, 12) for horizontal colorbar; (20, 12) for vertical colorbar

# Determine spacing of contours
lev4 = np.arange(grid.h.values.min()-1, 200, 1)
cmap2 = cmocean.cm.deep
cmap2.set_under('darkgray')

lev5 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Plot the bathymetry
cs5 = ax4.contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev4, cmap=cmap1, extend='both')
cs6 = ax4.contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev5, colors='dimgray')


# add title
#ax2.set_title('The Rivers in the Beaufort Sea Model Domain', fontsize=fontsize, y=1.08)
 
# format axes
ax4.set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax4.set_ylabel('Latitude (degrees)', fontsize=fontsize)
#plt.grid(True)
 
# specify colorbar
#cbar6 = fig6.colorbar(cs9, label='Depth (meters)', orientation='horizontal')
cbar4 = fig4.colorbar(cs4, label='Depth (m)', orientation='vertical')
cbar4.set_label(label='Depth (meters)', size=fontsize)



# -----------------------------------------------------------
# ------ Plot 5: Subplot of two figures from above ------
# -------------------------------------------------------------
# Make a subplot with the two figures from above together as 
# one figure (for the paper). Remove the title and have them share
# the colorbar and make a few other edits from other version
# (15,30)

# Determine spacing of contours
cmap7 = cmocean.cm.deep
lev6 = np.arange(grid.h.values.min()-1, 100, 1)
lev7 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# Make the figure
fig5, ax5 = plt.subplots(2, figsize=(17,16)) # (width, height), (15, 12) for horizontal colorbar; (20, 16) for vertical colorbaroutside plots

# --- Top subplot - Grid ---
# Plot the bathymetry
cs7 = ax5[0].contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev6, cmap=cmap7, extend='max')
# Plot land as gray 
# Make reverse mask so land can be gray 
reverse_colors = colors.ListedColormap(['darkgray','#FF000000'])
bounds=[0,1,2]
norm = colors.BoundaryNorm(bounds, reverse_colors.N)
# Make new mask so land is gray - NEW
mask_rho2 = grid.mask_rho.values
mask_rho2[mask_rho2==1.0] = 2.0
#mask_rho2_trimmed = mask_rho2[:,c_west:-c_west]
# Plot land
cs0f = ax5[0].contourf(grid.lon_rho.values, grid.lat_rho.values, mask_rho2, cmap=reverse_colors, norm=norm)

# Plot a few key towns
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517
sc3 = ax5[0].scatter(grid.lon_rho[eta_kakt_idx, xi_kakt_idx].values, grid.lat_rho[eta_kakt_idx, xi_kakt_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kaktovik Weather Station')

# Prudhoe Bay
eta_prud_idx = 37
xi_prud_idx = 279
sc4 = ax5[0].scatter(grid.lon_rho[eta_prud_idx, xi_prud_idx].values, grid.lat_rho[eta_prud_idx, xi_prud_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Prudhoe Bay Weather Station')

# Plot the CODA moorings
# Mooring 2
sc1 = ax5[0].scatter(grid.lon_rho[eta_rho_moor1, xi_rho_moor1].values, grid.lat_rho[eta_rho_moor1, xi_rho_moor1].values, 
            marker='x', s=300, linewidth=4, color='orange', label='CODA Mooring 2')

# Mooring 3
sc2 = ax5[0].scatter(grid.lon_rho[eta_rho_moor2, xi_rho_moor2].values, grid.lat_rho[eta_rho_moor2, xi_rho_moor2].values, 
            marker='x', s=300, linewidth=4, color='purple', label='CODA Mooring 3')

# Point 1 
sc5 = ax5[0].scatter(grid.lon_rho[eta_point1, xi_point1].values, grid.lat_rho[eta_point1, xi_point1].values, 
            marker='+', s=300, linewidth=4, color='darkgoldenrod', label='Harrison Bay')

# Point 2
#sc6 = ax3[0].scatter(grid.lon_rho[eta_point2, xi_point2].values, grid.lat_rho[eta_point2, xi_point2].values, 
 #           marker='+', s=300, linewidth=4, color='lightcoral', label='Point 2')

# Plot grid lines
# This puts the vertical lines through the longitudes over all latitudes 
ax5[0].plot(grid.lon_rho[:, ::16].values, grid.lat_rho[:, ::16].values, '-', color='k')

# This plots the vertical line on the right boundary of the grid
ax5[0].plot(grid.lon_rho[:, -1].values, grid.lat_rho[:, -1].values, '-', color='k')

# This plots the upper horizontal boundary of the grid
ax5[0].plot(grid.lon_rho[-1, :].values, grid.lat_rho[-1, :].values, '-', color='k')

for i in range(12):
    # This plots the horizontal lines (including the bottom boundary)
    ax5[0].plot(grid.lon_rho[(i*16), :].values, grid.lat_rho[(i*16), :].values, '-', color='k')
    


# Add title
#ax1.set_title('The Grid Cells (every 16th) in the Beaufort Sea', fontsize=30, y=1.08)
 
# Format axes
#ax3[0].set_xlabel('Longitude \n(degrees)', fontsize=30)
ax5[0].set_ylabel('Latitude \n(degrees)', fontsize=fontsize, rotation=0, labelpad=60, va='center') # 118
#plt.grid(True)

ax5[0].legend(fontsize=14, labelspacing=0.2, ncol=2, loc='lower left')

# Add a subset of Alaska
# Old - external vertical colorbar
#left, bottom, width, height = [0.62, 0.75, 0.18, 0.14]
#ax2 = fig3.add_axes([left, bottom, width, height])
# New - internal horizontal colorbar
left, bottom, width, height = [0.77, 0.75, 0.18, 0.14]
ax5b = fig5.add_axes([left, bottom, width, height])
# Add Alaska inset (new - horizontal internal colorbar)
# Create the polygon outline for Alaska 
polygon = Polygon([(-170,50),(-170,72),(-140, 72),(-140,50)])
# Plot this Alaska polygon
world.clip(polygon).plot(color='lightblue', linewidth=0.8, edgecolor='0.8', ax=ax5b) 
# Create a Rectangle patch to place around grid site
rect = patches.Rectangle((-155, 69.8), 14.5, 2.5, linewidth=3, angle=-4, edgecolor='r', facecolor='none')
# Add the patch to the Axes
ax5b.add_patch(rect)
# Add letters to outline
#plt.text(0.70, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # old, vertical external colorbar
plt.text(0.86, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # new, horizontal internal colorbar

# Adjust plot labels and such
ax5b.grid(False)
ax5b.set_xlim(-170, -135)
ax5b.set_ylim(53, 73)
plt.setp(ax5b.get_xticklabels(), visible=False)
plt.setp(ax5b.get_yticklabels(), visible=False)
plt.tick_params(left = False, bottom = False)

# =============================================================================
# # Add in a scale bar 
# import matplotlib_scalebar 
# from matplotlib_scalebar.scalebar import ScaleBar
# from sklearn.metrics.pairwise import haversine_distances
# lat1 = 70.75
# lon1 = -152
# lon2 = -151
# A = [lon1*np.pi/180., lat1*np.pi/180.]
# B = [lon2*np.pi/180., lat1*np.pi/180.]
# scale_dist = (6371000)*haversine_distances([A,B])[0,1]
# ax3[0].add_artist(ScaleBar(dx=scale_dist, units='m', font_properties={'size':'large'}, 
#                            location='lower center'))
# =============================================================================


# --- end of Alaska attempt ---



# --- Bottom subplot - Rivers ---
# Plot the bathymetry
cs8 = ax5[1].contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev6, cmap=cmap7, extend='max')
cs9 = ax5[1].contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev7, colors='dimgray')

# Plot land as gray 
cs0g = ax5[1].contourf(grid.lon_rho.values, grid.lat_rho.values, mask_rho2, cmap=reverse_colors, norm=norm)


# =============================================================================
# # Plot the Flaxman Island and Harrison Bay transects
# # Flaxman Island 
# p1 = ax3[1].plot(grid.lon_rho[flax_eta, flax_xi].values, grid.lat_rho[flax_eta, flax_xi].values, color='red', 
#          linewidth=5, label='Flaxman Island')
# 
# # Harrison Bay
# p2 = ax3[1].plot(grid.lon_rho[har_eta, har_xi].values, grid.lat_rho[har_eta, har_xi].values, color='gold', 
#          linewidth=5, label='Harrison Bay')
# 
# # Plot the other two transects used 
# # Transect 1
# p3 = ax3[1].plot(grid.lon_rho[eta1, xi1].values, grid.lat_rho[eta1, xi1].values, color='purple', 
#          linewidth=5, label='Transect 1')
# 
# # Transect 2
# p4 = ax3[1].plot(grid.lon_rho[eta2, xi2].values, grid.lat_rho[eta2, xi2].values, color='navy', 
#          linewidth=5, label='Transect 2')
# =============================================================================



# Plot the 14 rivers in the grid
# Go from West to East
# Kalikpik River
eta_kal_idx = 23 #22
xi_kal_idx = 87
s1 = ax5[1].scatter(grid.lon_rho[eta_kal_idx, xi_kal_idx].values, grid.lat_rho[eta_kal_idx, xi_kal_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
eta_fis_idx = 20
xi_fis_idx = 117 #116
s2 = ax5[1].scatter(grid.lon_rho[eta_fis_idx, xi_fis_idx].values, grid.lat_rho[eta_fis_idx, xi_fis_idx].values, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
eta_col_idx = 39
xi_col_idx = 166 #166
s3 = ax5[1].scatter(grid.lon_rho[eta_col_idx, xi_col_idx].values, grid.lat_rho[eta_col_idx, xi_col_idx].values, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
eta_sak_idx = 46 #45
xi_sak_idx = 234
s4 = ax5[1].scatter(grid.lon_rho[eta_sak_idx, xi_sak_idx].values, grid.lat_rho[eta_sak_idx, xi_sak_idx].values, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s5 = ax5[1].scatter(grid.lon_rho[eta_kup_idx, xi_kup_idx].values, grid.lat_rho[eta_kup_idx, xi_kup_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax3[1].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#eta_faw_idx = 44 #43
#xi_faw_idx = 249
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
eta_put_idx = 28 #27
xi_put_idx = 264
s8 = ax5[1].scatter(grid.lon_rho[eta_put_idx, xi_put_idx].values, grid.lat_rho[eta_put_idx, xi_put_idx].values, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
eta_sag_idx = 37 #36
xi_sag_idx = 279
s9 = ax5[1].scatter(grid.lon_rho[eta_sag_idx, xi_sag_idx].values, grid.lat_rho[eta_sag_idx, xi_sag_idx].values, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
eta_sta_idx = 27 #26
xi_sta_idx = 393
s10 = ax5[1].scatter(grid.lon_rho[eta_sta_idx, xi_sta_idx].values, grid.lat_rho[eta_sta_idx, xi_sta_idx].values, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
eta_can_idx = 20 #19
xi_can_idx = 416
s11 = ax5[1].scatter(grid.lon_rho[eta_can_idx, xi_can_idx].values, grid.lat_rho[eta_can_idx, xi_can_idx].values, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
eta_kat_idx = 9 #8
xi_kat_idx = 447
s12 = ax5[1].scatter(grid.lon_rho[eta_kat_idx, xi_kat_idx].values, grid.lat_rho[eta_kat_idx, xi_kat_idx].values, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
eta_hul_idx = 40
xi_hul_idx = 489
s13 = ax5[1].scatter(grid.lon_rho[eta_hul_idx, xi_hul_idx].values, grid.lat_rho[eta_hul_idx, xi_hul_idx].values, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
eta_jag_idx = 62 #61
xi_jag_idx = 528
s14 = ax5[1].scatter(grid.lon_rho[eta_jag_idx, xi_jag_idx].values, grid.lat_rho[eta_jag_idx, xi_jag_idx].values, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
eta_sik_idx = 46
xi_sik_idx = 574 #573
s15 = ax5[1].scatter(grid.lon_rho[eta_sik_idx, xi_sik_idx].values, grid.lat_rho[eta_sik_idx, xi_sik_idx].values, 
            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')


# Plot the transects
# Set the indices
xi1 = 152
xi2 = 304
xi3 = 456

# Transect 1
ax5[1].plot(grid.lon_rho[:,xi1].values, grid.lat_rho[:,xi1].values, '-', linewidth=5, color='#D81B60', label='Transect 1')
# Transect 2
ax5[1].plot(grid.lon_rho[:,xi2].values, grid.lat_rho[:,xi2].values, '-', linewidth=5, color='#1E88E5', label='Transect 2')
# Transect 3
ax5[1].plot(grid.lon_rho[:,xi3].values, grid.lat_rho[:,xi3].values, '-', linewidth=5, color='#FFC107', label='Transect 3')


# New attempt to include transects in legend
# =============================================================================
# lines=[s1,s2,s3,s4,s5,s6,s8,s9,s10,s11,s12,s13,s14,s15]
# include=[0,1,2,3,4,5,6]
# #legend3 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc='lower left', fontsize=14, ncol=3,
#  #                    columnspacing=0.1, labelspacing=0.2)
# legend4 = ax3[1].legend([lines[i] for i in [7,8,9,10,11,12,13]],['Sagavanirktok', 'Staines', 'Canning','Katakturuk', 'Hulahula', 'Jago', 
#                                                                'Siksik'], loc='upper right', fontsize=14, ncol=3, columnspacing=0.1, 
#                         labelspacing=0.2)
# curves = [p1,p2,p3,p4]
# include2 = [0,1,2,3]
# #legend5 = ax3[1].legend()
# legend5 = plt.legend([curves[i] for i in [0,1,2,3]], ['Harrison Bay', 'Transect 1', 'Flaxman Island', 'Transect 2'], 
#                         loc='lower right', fontsize=14, ncol=2, columnspacing=0.1, labelspacing=0.2)
# plt.gca().add_artist(legend5)
# #plt.gca().add_artist(legend5)
# =============================================================================

# Temporary not too bad legend
ax5[1].legend(loc='lower left', ncol=4, fontsize=14, columnspacing=0.75, labelspacing=0.1)

# Format axes
ax5[1].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax5[1].set_ylabel('Latitude \n(degrees)', fontsize=fontsize, rotation=0, labelpad=60, va='center') #118

# Specify colorbar
#cbar3 = fig3.colorbar(cs6, label='Depth (m)', orientation='vertical', ax=ax3)
#cbar3.set_label(label='Depth (meters)', size=fontsize)
# If horizontal and in axes
cbar5 = plt.colorbar(cs7, cax=ax5[1].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 20, 40, 60, 80], ax=ax5[1], 
                     orientation='horizontal').set_label(label='Depth (m)', size=fontsize-2)

# Label subplots
# put labels in top right corners
#plt.text(0.715, 0.853, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
#plt.text(0.715, 0.442, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
# put labels in bottom right corners
plt.text(0.870, 0.553, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.870, 0.142, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)

#plt.tight_layout()






# -------------------------------------------------
# ------ Plot 6: Ocean Bathymetry and Grid Cells ------
# -------------------------------------------------
# Make ansimilar plot to above but accentuate shelf bathymetry more and
# put grid cells

# Make a plot of bathymetry with grid lines
# depths for the colorbar
fig6, ax6 = plt.subplots(figsize=(25,12)) # (15, 12) for horizontal colorbar; (20, 12) for vertical colorbar

# Determine spacing of contours
lev6 = np.arange(grid.h.values.min()-1, 200, 1)
cmap6 = cmocean.cm.deep
cmap6.set_under('darkgray')

# Make reverse mask so land can be gray 
reverse_colors = colors.ListedColormap(['darkgray','#FF000000'])
bounds=[0,1,2]
norm = colors.BoundaryNorm(bounds, reverse_colors.N)

# Plot the bathymetry
cs10 = ax6.contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev6, cmap=cmap6, extend='max')
#cs11 = ax6.contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values*grid.mask_rho.values, lev3, colors='dimgray')

# Make new mask so land is gray - NEW
mask_rho2 = grid.mask_rho.values
mask_rho2[mask_rho2==1.0] = 2.0
cs6b = ax6.contourf(grid.lon_rho.values, grid.lat_rho.values, mask_rho2, cmap=reverse_colors, norm=norm)

# Plot grid lines
# This puts the vertical lines through the longitudes over all latitudes 
ax6.plot(grid.lon_rho[:, ::16].values, grid.lat_rho[:, ::16].values, '-', color='k')

# This plots the vertical line on the right boundary of the grid
ax6.plot(grid.lon_rho[:, -1].values, grid.lat_rho[:, -1].values, '-', color='k')

# This plots the upper horizontal boundary of the grid
ax6.plot(grid.lon_rho[-1, :].values, grid.lat_rho[-1, :].values, '-', color='k')

for i in range(12):
    # This plots the horizontal lines (including the bottom boundary)
    ax6.plot(grid.lon_rho[(i*16), :].values, grid.lat_rho[(i*16), :].values, '-', color='k')
 
# format axes
ax6.set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax6.set_ylabel('Latitude (degrees)', fontsize=fontsize)
#plt.grid(True)
 
# specify colorbar
#cbar6 = fig6.colorbar(cs9, label='Depth (meters)', orientation='horizontal')
cbar6 = fig6.colorbar(cs10, label='Depth (m)', orientation='vertical', pad=0.01)
cbar6.set_label(label='Depth (m)', size=fontsize)

# Save the plot
#plt.savefig('/Users/brun1463/Desktop/Year_5/CSDMS/grid_transparent_001.png', transparent=True, bbox_inches='tight', pad_inches=0)





# -----------------------------------------------------------
# ------ Plot 7: Subplot of two figures from above ------
# ------------------ with mud fraction -----------------------
# -------------------------------------------------------------
# Make a subplot with the two figures from above together as 
# one figure (for the paper). Remove the title and made she shading 
# in the second one the percentage of mud 

# Load in the initial conditions file
ini_conds = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/Forcing_files/initial_conds_beaufort_shelf_2020_001.nc') # restrat_mix_ini_kaktovik_biggrid013.nc

# Add together the percentages of the two mud classes
tot_mud_percent = ini_conds.mudfrac_01 + ini_conds.mudfrac_02

# Determine spacing of contours
cmap7 = cmocean.cm.deep
lev6 = np.arange(grid.h.values.min()-1, 100, 1)
lev7 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
lev_mud_frac = np.arange(0,1,0.01)
cmap_sed = cmocean.cm.turbid

# Make the figure
fig7, ax7 = plt.subplots(2, figsize=(17,16)) # (width, height), (15, 12) for horizontal colorbar; (20, 16) for vertical colorbaroutside plots

# --- Top subplot - Grid ---
# Plot the bathymetry
cs11 = ax7[0].contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev6, cmap=cmap7, extend='max')
# Plot land as gray 
# Make reverse mask so land can be gray 
reverse_colors = colors.ListedColormap(['darkgray','#FF000000'])
bounds=[0,1,2]
norm = colors.BoundaryNorm(bounds, reverse_colors.N)
# Make new mask so land is gray - NEW
mask_rho2 = grid.mask_rho.values
mask_rho2[mask_rho2==1.0] = 2.0
#mask_rho2_trimmed = mask_rho2[:,c_west:-c_west]
# Plot land
cs11f = ax7[0].contourf(grid.lon_rho.values, grid.lat_rho.values, mask_rho2, cmap=reverse_colors, norm=norm)

# Plot a few key towns
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517
sc6 = ax7[0].scatter(grid.lon_rho[eta_kakt_idx, xi_kakt_idx].values, grid.lat_rho[eta_kakt_idx, xi_kakt_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kaktovik Weather Station')

# Prudhoe Bay
eta_prud_idx = 37
xi_prud_idx = 279
sc7 = ax7[0].scatter(grid.lon_rho[eta_prud_idx, xi_prud_idx].values, grid.lat_rho[eta_prud_idx, xi_prud_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Prudhoe Bay Weather Station')

# Plot the CODA moorings
# Mooring 2
sc8 = ax7[0].scatter(grid.lon_rho[eta_rho_moor1, xi_rho_moor1].values, grid.lat_rho[eta_rho_moor1, xi_rho_moor1].values, 
            marker='x', s=300, linewidth=4, color='orange', label='CODA Mooring 2')

# Mooring 3
sc9 = ax7[0].scatter(grid.lon_rho[eta_rho_moor2, xi_rho_moor2].values, grid.lat_rho[eta_rho_moor2, xi_rho_moor2].values, 
            marker='x', s=300, linewidth=4, color='purple', label='CODA Mooring 3')

# Point 1 
sc10 = ax7[0].scatter(grid.lon_rho[eta_point1, xi_point1].values, grid.lat_rho[eta_point1, xi_point1].values, 
            marker='+', s=300, linewidth=4, color='darkgoldenrod', label='Harrison Bay')

# Point 2
#sc6 = ax3[0].scatter(grid.lon_rho[eta_point2, xi_point2].values, grid.lat_rho[eta_point2, xi_point2].values, 
 #           marker='+', s=300, linewidth=4, color='lightcoral', label='Point 2')

# Plot grid lines
# This puts the vertical lines through the longitudes over all latitudes 
ax7[0].plot(grid.lon_rho[:, ::16].values, grid.lat_rho[:, ::16].values, '-', color='k')

# This plots the vertical line on the right boundary of the grid
ax7[0].plot(grid.lon_rho[:, -1].values, grid.lat_rho[:, -1].values, '-', color='k')

# This plots the upper horizontal boundary of the grid
ax7[0].plot(grid.lon_rho[-1, :].values, grid.lat_rho[-1, :].values, '-', color='k')

for i in range(12):
    # This plots the horizontal lines (including the bottom boundary)
    ax7[0].plot(grid.lon_rho[(i*16), :].values, grid.lat_rho[(i*16), :].values, '-', color='k')
    


# Add title
#ax1.set_title('The Grid Cells (every 16th) in the Beaufort Sea', fontsize=30, y=1.08)
 
# Format axes
#ax3[0].set_xlabel('Longitude \n(degrees)', fontsize=30)
ax7[0].set_ylabel('Latitude \n(degrees)', fontsize=fontsize, rotation=0, labelpad=60, va='center') # 118
#plt.grid(True)

ax7[0].legend(fontsize=14, labelspacing=0.2, ncol=2, loc='lower left')

# Add a subset of Alaska
# Old - external vertical colorbar
#left, bottom, width, height = [0.62, 0.75, 0.18, 0.14]
#ax2 = fig3.add_axes([left, bottom, width, height])
# New - internal horizontal colorbar
left, bottom, width, height = [0.77, 0.75, 0.18, 0.14]
ax7b = fig7.add_axes([left, bottom, width, height])
# Add Alaska inset (new - horizontal internal colorbar)
# Create the polygon outline for Alaska 
polygon = Polygon([(-170,50),(-170,72),(-140, 72),(-140,50)])
# Plot this Alaska polygon
world.clip(polygon).plot(color='lightblue', linewidth=0.8, edgecolor='0.8', ax=ax7b) 
# Create a Rectangle patch to place around grid site
rect = patches.Rectangle((-155, 69.8), 14.5, 2.5, linewidth=3, angle=-4, edgecolor='r', facecolor='none')
# Add the patch to the Axes
ax7b.add_patch(rect)
# Add letters to outline
#plt.text(0.70, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # old, vertical external colorbar
plt.text(0.86, 0.83, 'AK', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure) # new, horizontal internal colorbar

# Adjust plot labels and such
ax7b.grid(False)
ax7b.set_xlim(-170, -135)
ax7b.set_ylim(53, 73)
plt.setp(ax7b.get_xticklabels(), visible=False)
plt.setp(ax7b.get_yticklabels(), visible=False)
plt.tick_params(left = False, bottom = False)

# =============================================================================
# # Add in a scale bar 
# import matplotlib_scalebar 
# from matplotlib_scalebar.scalebar import ScaleBar
# from sklearn.metrics.pairwise import haversine_distances
# lat1 = 70.75
# lon1 = -152
# lon2 = -151
# A = [lon1*np.pi/180., lat1*np.pi/180.]
# B = [lon2*np.pi/180., lat1*np.pi/180.]
# scale_dist = (6371000)*haversine_distances([A,B])[0,1]
# ax3[0].add_artist(ScaleBar(dx=scale_dist, units='m', font_properties={'size':'large'}, 
#                            location='lower center'))
# =============================================================================


# --- end of Alaska attempt ---



# --- Bottom subplot - Rivers ---
# Plot the bathymetry
#cs12 = ax7[1].contourf(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev6, cmap=cmap7, extend='max')
#cs13 = ax7[1].contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev7, colors='dimgray')

# Plot fraction of mud 
cs12 = ax7[1].contourf(grid.lon_rho.values, grid.lat_rho.values, tot_mud_percent[0,0,:,:], lev_mud_frac, cmap=cmap_sed, extend='max')
cs13 = ax7[1].contour(grid.lon_rho.values, grid.lat_rho.values, grid.h[:,:].values, lev7, colors='white')

# Plot land as gray 
cs12g = ax7[1].contourf(grid.lon_rho.values, grid.lat_rho.values, mask_rho2, cmap=reverse_colors, norm=norm)


# =============================================================================
# # Plot the Flaxman Island and Harrison Bay transects
# # Flaxman Island 
# p1 = ax3[1].plot(grid.lon_rho[flax_eta, flax_xi].values, grid.lat_rho[flax_eta, flax_xi].values, color='red', 
#          linewidth=5, label='Flaxman Island')
# 
# # Harrison Bay
# p2 = ax3[1].plot(grid.lon_rho[har_eta, har_xi].values, grid.lat_rho[har_eta, har_xi].values, color='gold', 
#          linewidth=5, label='Harrison Bay')
# 
# # Plot the other two transects used 
# # Transect 1
# p3 = ax3[1].plot(grid.lon_rho[eta1, xi1].values, grid.lat_rho[eta1, xi1].values, color='purple', 
#          linewidth=5, label='Transect 1')
# 
# # Transect 2
# p4 = ax3[1].plot(grid.lon_rho[eta2, xi2].values, grid.lat_rho[eta2, xi2].values, color='navy', 
#          linewidth=5, label='Transect 2')
# =============================================================================



# Plot the 14 rivers in the grid
# Go from West to East
# Kalikpik River
eta_kal_idx = 23 #22
xi_kal_idx = 87
s11 = ax7[1].scatter(grid.lon_rho[eta_kal_idx, xi_kal_idx].values, grid.lat_rho[eta_kal_idx, xi_kal_idx].values, 
            marker='.', s=300, linewidth=4, color='r', label='Kalikpik')

# Fish Creek
eta_fis_idx = 20
xi_fis_idx = 117 #116
s12 = ax7[1].scatter(grid.lon_rho[eta_fis_idx, xi_fis_idx].values, grid.lat_rho[eta_fis_idx, xi_fis_idx].values, 
            marker='x', s=100, linewidth=4, color='orange', label='Fish Creek')

# Colville River
eta_col_idx = 39
xi_col_idx = 166 #166
s13 = ax7[1].scatter(grid.lon_rho[eta_col_idx, xi_col_idx].values, grid.lat_rho[eta_col_idx, xi_col_idx].values, 
            marker='.', s=300, linewidth=4, color='brown', label='Colville')

# Sakonowyak River
eta_sak_idx = 46 #45
xi_sak_idx = 234
s14 = ax7[1].scatter(grid.lon_rho[eta_sak_idx, xi_sak_idx].values, grid.lat_rho[eta_sak_idx, xi_sak_idx].values, 
            marker='x', s=150, linewidth=6, color='green', label='Sakonowyak')

# Kuparik
# Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
# of the Kuparuk River
s15 = ax7[1].scatter(grid.lon_rho[eta_kup_idx, xi_kup_idx].values, grid.lat_rho[eta_kup_idx, xi_kup_idx].values, 
            marker='.', s=300, linewidth=4, color='b', label='Kuparuk')

# Kuparuk - commented out to move dot onto old Kukpuk since that is the main 
# channel of the Kuparuk 
#s6 = ax3[1].scatter(x_rho_flat_trimmed[xi_kup_idx]/1000, y_rho_flat[eta_kup_idx]/1000, 
 #           marker='x', s=100, linewidth=4, color='pink', label='Kuparuk')

# Fawn Creek
#eta_faw_idx = 44 #43
#xi_faw_idx = 249
#s7 = ax3[1].scatter(grid.lon_rho[eta_faw_idx, xi_faw_idx].values, grid.lat_rho[eta_faw_idx, xi_faw_idx].values, 
 #           marker='.', s=300, linewidth=4, color='darkviolet', label='Fawn Creek')

# Putuligayuk River
eta_put_idx = 28 #27
xi_put_idx = 264
s16 = ax7[1].scatter(grid.lon_rho[eta_put_idx, xi_put_idx].values, grid.lat_rho[eta_put_idx, xi_put_idx].values, 
            marker='x', s=100, linewidth=4, color='dodgerblue', label='Putuligayuk')

# Sagavanirktok River
eta_sag_idx = 37 #36
xi_sag_idx = 279
s17 = ax7[1].scatter(grid.lon_rho[eta_sag_idx, xi_sag_idx].values, grid.lat_rho[eta_sag_idx, xi_sag_idx].values, 
            marker='.', s=300, linewidth=3, color='deepskyblue', label='Sagavanirktok')

# Canning River
# Staines River
eta_sta_idx = 27 #26
xi_sta_idx = 393
s18 = ax7[1].scatter(grid.lon_rho[eta_sta_idx, xi_sta_idx].values, grid.lat_rho[eta_sta_idx, xi_sta_idx].values, 
            marker='x', s=100, linewidth=4, color='gold', label='Staines')

# Canning River
eta_can_idx = 20 #19
xi_can_idx = 416
s19 = ax7[1].scatter(grid.lon_rho[eta_can_idx, xi_can_idx].values, grid.lat_rho[eta_can_idx, xi_can_idx].values, 
            marker='.', s=300, linewidth=4, color='darkorange', label='Canning')

# Katakturuk River
eta_kat_idx = 9 #8
xi_kat_idx = 447
s20 = ax7[1].scatter(grid.lon_rho[eta_kat_idx, xi_kat_idx].values, grid.lat_rho[eta_kat_idx, xi_kat_idx].values, 
            marker='x', s=100, linewidth=4, color='aquamarine', label='Katakturuk')

# Hulahula River
eta_hul_idx = 40
xi_hul_idx = 489
s21 = ax7[1].scatter(grid.lon_rho[eta_hul_idx, xi_hul_idx].values, grid.lat_rho[eta_hul_idx, xi_hul_idx].values, 
            marker='.', s=300, linewidth=4, color='blueviolet', label='Hulahula')

# Jago River
eta_jag_idx = 62 #61
xi_jag_idx = 528
s22 = ax7[1].scatter(grid.lon_rho[eta_jag_idx, xi_jag_idx].values, grid.lat_rho[eta_jag_idx, xi_jag_idx].values, 
            marker='x', s=100, linewidth=4, color='magenta', label='Jago')

# Siksik River
eta_sik_idx = 46
xi_sik_idx = 574 #573
s23 = ax7[1].scatter(grid.lon_rho[eta_sik_idx, xi_sik_idx].values, grid.lat_rho[eta_sik_idx, xi_sik_idx].values, 
            marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')


# Plot the transects
# Set the indices
xi1 = 152
xi2 = 304
xi3 = 456

# Transect 1
ax7[1].plot(grid.lon_rho[41:165,xi1].values, grid.lat_rho[41:165,xi1].values, '-', linewidth=5, color='#D81B60', label='Transect 1')
# Transect 2
ax7[1].plot(grid.lon_rho[17:150,xi2].values, grid.lat_rho[17:150,xi2].values, '-', linewidth=5, color='#1E88E5', label='Transect 2')
# Transect 3
ax7[1].plot(grid.lon_rho[12:125,xi3].values, grid.lat_rho[12:125,xi3].values, '-', linewidth=5, color='#FFC107', label='Transect 3')


# New attempt to include transects in legend
# =============================================================================
# lines=[s1,s2,s3,s4,s5,s6,s8,s9,s10,s11,s12,s13,s14,s15]
# include=[0,1,2,3,4,5,6]
# #legend3 = plt.legend([lines[i] for i in include],[lines[i].get_label() for i in include], loc='lower left', fontsize=14, ncol=3,
#  #                    columnspacing=0.1, labelspacing=0.2)
# legend4 = ax3[1].legend([lines[i] for i in [7,8,9,10,11,12,13]],['Sagavanirktok', 'Staines', 'Canning','Katakturuk', 'Hulahula', 'Jago', 
#                                                                'Siksik'], loc='upper right', fontsize=14, ncol=3, columnspacing=0.1, 
#                         labelspacing=0.2)
# curves = [p1,p2,p3,p4]
# include2 = [0,1,2,3]
# #legend5 = ax3[1].legend()
# legend5 = plt.legend([curves[i] for i in [0,1,2,3]], ['Harrison Bay', 'Transect 1', 'Flaxman Island', 'Transect 2'], 
#                         loc='lower right', fontsize=14, ncol=2, columnspacing=0.1, labelspacing=0.2)
# plt.gca().add_artist(legend5)
# #plt.gca().add_artist(legend5)
# =============================================================================

# Temporary not too bad legend
ax7[1].legend(loc='lower left', ncol=4, fontsize=14, columnspacing=0.75, labelspacing=0.1)

# Format axes
ax7[1].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax7[1].set_ylabel('Latitude \n(degrees)', fontsize=fontsize, rotation=0, labelpad=60, va='center') #118

# Specify colorbar - top plot
#cbar3 = fig3.colorbar(cs6, label='Depth (m)', orientation='vertical', ax=ax3)
#cbar3.set_label(label='Depth (meters)', size=fontsize)
# If horizontal and in axes
cbar7 = plt.colorbar(cs11, cax=ax7[0].inset_axes((0.55, 0.92, 0.32, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 20, 40, 60, 80], ax=ax7[0], 
                     orientation='horizontal').set_label(label='Depth (m)', size=fontsize-2)

# Specify colorbar - bottom plot 
# If horizontal and in axes
cbar8 = plt.colorbar(cs12, cax=ax7[1].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, .2, .4, .6, .8, 1], ax=ax7[1], 
                     orientation='horizontal').set_label(label='Mud Fraction', size=fontsize-2)

# Label subplots
# put labels in top right corners
#plt.text(0.715, 0.853, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
#plt.text(0.715, 0.442, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
# put labels in bottom right corners
plt.text(0.870, 0.553, 'a)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)
plt.text(0.870, 0.142, 'b)', fontsize=fontsize, fontweight='bold', transform=plt.gcf().transFigure)

#plt.tight_layout()

