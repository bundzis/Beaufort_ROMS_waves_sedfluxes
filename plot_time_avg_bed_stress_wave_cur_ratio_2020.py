##################### Time-Averaged Bed Stress #######################
# The purpose of this script is to look at the time-averaged 
# bed shear stress from waves, currents, the two combined, and the 
# ratio of the two.
#
# 
# Notes:
# - The script is being edited to loop through model output
##################################################################################


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
grid = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Grids/KakAKgrd_shelf_big010_smooth006.nc')

# Pull out some dimensions
eta_rho_len = len(grid.eta_rho)
xi_rho_len = len(grid.xi_rho)

# Load in the rho masks 
#mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
#mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')


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



# Make a function to process the model output files
def process_bed_shear_stress_output(filename):
    """
    This function takes a specified model output netcdf
    and does a bunch of processing of the data, then outputs 
    the desired values.
    
    Inputs:
    - filename: string, path to and name of netcdf output
    
    Output:
    - bustrc_mean_tmp: time-averaged current-induced bed stress in u-direction
    - bvstrc_mean_tmp: time-averaged current-induced bed stress in v-direction
    - bustrw_mean_tmp: time-averaged wave-induced bed stress in u-direction
    - bvstrw_mean_tmp: time-averaged wave-induced bed stress in v-direction
    - bstrcwmax_mean_tmp: time-averaged current and wave-induced bed stress
    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Take the temporal average of the bed stress
    # Current-induced = bustrc & bvstrc
    bustrc_mean_tmp = model_output['bustrc'].mean(dim='ocean_time').copy()
    bvstrc_mean_tmp = model_output['bvstrc'].mean(dim='ocean_time').copy()
    # Wave-induced = bustrw & bvstrw
    bustrw_mean_tmp = model_output['bustrw'].mean(dim='ocean_time').copy()
    bvstrw_mean_tmp = model_output['bvstrw'].mean(dim='ocean_time').copy()
    # Combined Wave-Current = bstrcwmax
    bstrcwmax_mean_tmp = model_output['bstrcwmax'].mean(dim='ocean_time').copy()
    
    # Return the means
    return(bustrc_mean_tmp, bvstrc_mean_tmp, bustrw_mean_tmp, bvstrw_mean_tmp, bstrcwmax_mean_tmp)


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
bustrc_means = np.empty((num_files, eta_rho_len, xi_rho_len))
bvstrc_means = np.empty((num_files, eta_rho_len, xi_rho_len))
bustrw_means = np.empty((num_files, eta_rho_len, xi_rho_len))
bvstrw_means = np.empty((num_files, eta_rho_len, xi_rho_len))
bstrcwmax_means = np.empty((num_files, eta_rho_len, xi_rho_len))

# Loop through the model output
for j in range(num_files):
        
    # Call the function to open the output
    bustrc_tmp, bvstrc_tmp, bustrw_tmp, bvstrw_tmp, bstrcwmax_tmp = process_bed_shear_stress_output(file_names2[j])
    
    # Save these to the arrays
    bustrc_means[j,:,:] = bustrc_tmp
    bvstrc_means[j,:,:] = bvstrc_tmp
    bustrw_means[j,:,:] = bustrw_tmp
    bvstrw_means[j,:,:] = bvstrw_tmp
    bstrcwmax_means[j,:,:] = bstrcwmax_tmp
    
# Now that the averages have been found for all outputs,
# combine them into one average
# Be sure to weight them all by the portion of the full time that they represent
# bustrc & bvstrc

# Go through and weight them by the length of that file
# Make some empty arrays to hold final values
bustrc_mean_total = np.zeros((eta_rho_len, xi_rho_len))
bvstrc_mean_total = np.zeros((eta_rho_len, xi_rho_len))
bustrw_mean_total = np.zeros((eta_rho_len, xi_rho_len))
bvstrw_mean_total = np.zeros((eta_rho_len, xi_rho_len))
bstrcwmax_total = np.zeros((eta_rho_len, xi_rho_len))

for g in range(num_files):
    # Multiply each mean by length of time output
    # for that file, divide by total length of
    # time of model run
    # bustrc
    buc_tmp = (bustrc_means[g,:,:]*time_lengths[g])/full_time_len
    bustrc_mean_total = bustrc_mean_total + buc_tmp
    # bvstrc
    bvc_tmp = (bvstrc_means[g,:,:]*time_lengths[g])/full_time_len
    bvstrc_mean_total = bvstrc_mean_total + bvc_tmp
    # bustrw
    buw_tmp = (bustrw_means[g,:,:]*time_lengths[g])/full_time_len
    bustrw_mean_total = bustrw_mean_total + buw_tmp
    # bvstrw
    bvw_tmp = (bvstrw_means[g,:,:]*time_lengths[g])/full_time_len
    bvstrw_mean_total = bvstrw_mean_total + bvw_tmp
    # bstrcwmax
    bstrcwm_tmp = (bstrcwmax_means[g,:,:]*time_lengths[g])/full_time_len
    bstrcwmax_total = bstrcwmax_total + bstrcwm_tmp



# Load in the last model output to use for plotting things
model_output = xr.open_dataset(file_names2[-1])



# -------------------------------------------------------------------
# ---------- Plot 1: Time-Averaged Bed Stress ----------
# -------------------------------------------------------------------
# Make a plot of the time-averaged bed shear stress
# with subplots for current-induced, waves-induced,
# the two combined, and the ratio of the two. Since this 
# has both u and v components, the total 
# magnitude on the rho points must be calculated

# Combine the u and v
# Current-induced 
bstrc_tot_mag = np.sqrt(((bustrc_mean_total)**2) + ((bvstrc_mean_total)**2))
# Wave-induced
bstrw_tot_mag = np.sqrt(((bustrw_mean_total)**2) + ((bvstrw_mean_total)**2))

# Replace the nans with -1 to help with plotting out of bounds data
# Current-induced
temp_mask = model_output.mask_rho.copy()
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)
bstrc_tot_mag_wland = bstrc_tot_mag * temp_mask
bstrc_tot_mag_wland = np.nan_to_num(bstrc_tot_mag_wland, nan=-1)
#bstrc_tot_mag_wland = bstrc_tot_mag_wland.fillna(-1)
# Wave-induced
bstrw_tot_mag_wland = bstrw_tot_mag * temp_mask
bstrw_tot_mag_wland = np.nan_to_num(bstrw_tot_mag_wland, nan=-1)
#bstrw_tot_mag_wland = bstrw_tot_mag_wland.fillna(-1)
# Current and wave induced
bstrcwmax_total_wland = bstrcwmax_total * temp_mask
bstrcwmax_total_wland = np.nan_to_num(bstrcwmax_total_wland, nan=-1)
#bstrcwmax_total_wland = bstrcwmax_total_wland.fillna(-1)


# Make the figure
fig1, ax1 = plt.subplots(4, figsize=(18, 25))  # (18,16) (16,18) (16,21)

# Plot current-induced
# Set the colormap
#cmap = matplotlib.cm.cividis
cmap1=cmocean.cm.amp
cmap1.set_under('darkgray')
#cmap2.set_bad('darkgray')

# Set a colormap and levels to use for all plots
cmap1b = cmocean.cm.tempo
cmap1b.set_under('darkgray')
lev1b = np.arange(0-0.000001, 0.25, 0.01)


# Set the levels 
print(np.max(bstrc_tot_mag), flush=True)
print(np.min(bstrc_tot_mag), flush=True)
lev1  = np.arange(-0.00001, np.nanmax(bstrc_tot_mag), 0.01)

# Plot the bathymetry in 10m contours from 10 - 100 m
lev2 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax1[0].contour(model_output.lon_rho.values, model_output.lat_rho.values,
             model_output.bath[0,:,:].values, lev2, colors='gray')

cs1 = ax1[0].contourf(model_output.lon_rho[:,:].values, model_output.lat_rho[:,:].values, 
                  bstrc_tot_mag_wland[:,:], lev1b, cmap=cmap1b, extend='both')

# Label the plot
#ax1[0].set_title('Time-Averaged Current-Induced Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
plt.setp(ax1[0].get_xticklabels(), visible=False)
#ax1[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[0].set_ylabel('Latitude \n(degrees)', fontsize=fontsize)
#cbar1 = plt.colorbar(cs1, orientation='vertical', ax=ax1[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


# Plot  wave-induced
# Set the colormap
#cmap = matplotlib.cm.cividis
cmap2=cmocean.cm.tempo
cmap2.set_under('darkgray')

# Set the levels 
lev3  = np.arange(-0.00001, np.nanmax(bstrw_tot_mag), 0.01)

# Plot the bathymetry in 10m contours from 10 - 100 m
ax1[1].contour(model_output.lon_rho.values, model_output.lat_rho.values,
             model_output.bath[0,:,:].values, lev2, colors='gray')

cs2 = ax1[1].contourf(model_output.lon_rho[:,:].values, model_output.lat_rho[:,:].values, 
                  bstrw_tot_mag_wland[:,:], lev1b, cmap=cmap1b, extend='both')

# Label the plot
#ax1[1].set_title('Time-Averaged Wave-Induced Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
plt.setp(ax1[1].get_xticklabels(), visible=False)
#ax1[1].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[1].set_ylabel('Latitude \n(degrees)', fontsize=fontsize)
#cbar2 = plt.colorbar(cs2, orientation='vertical', ax=ax1[1]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


# Plot  wave and current induced
# Set the colormap
#cmap = matplotlib.cm.cividis
cmap3=cmocean.cm.dense
cmap3.set_under('darkgray')

# Set the levels 
lev4  = np.arange(-0.00001, np.nanmax(bstrcwmax_total), 0.01)

# Plot the bathymetry in 10m contours from 10 - 100 m
ax1[2].contour(model_output.lon_rho.values, model_output.lat_rho.values,
             model_output.bath[0,:,:].values, lev2, colors='gray')

cs3 = ax1[2].contourf(model_output.lon_rho[:,:].values, model_output.lat_rho[:,:].values, 
                  bstrcwmax_total_wland[:,:], lev1b, cmap=cmap1b, extend='both')

# Label the plot
#ax1[2].set_title('Time-Averaged Total Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
plt.setp(ax1[2].get_xticklabels(), visible=False)
#ax1[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[2].set_ylabel('Latitude \n(degrees)', fontsize=fontsize)

# Adjust colorbar placement
#bottom, top = 0.1, 0.9
#left, right = 0.1, 0.8
#fig1.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
#cbar1_ax = fig1.add_axes([0.85, bottom, 0.05, top-bottom])
#cbar3 = plt.colorbar(cs3, orientation='vertical', cax=cbar1_ax, ax=[ax1[0], ax1[1], ax1[2]]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)



# Plot ratio of current over wave
# Take the ratio
bstrc_over_bstrw = bstrc_tot_mag/bstrw_tot_mag

# Replace the nans with -100 to help with plotting out of bounds data
bstrc_over_bstrw_wland = bstrc_over_bstrw * temp_mask

# TEMP
print('bstrc_tot_mag[0,0]: ', bstrc_tot_mag[0,0])
print('bstrw_tot_mag[0,0]: ', bstrw_tot_mag[0,0])
print('bstrc_tot_mag[-1,-1]: ', bstrc_tot_mag[-1,-1])
print('bstrw_tot_mag[-1,-1]: ', bstrw_tot_mag[-1,-1])
print('bstrc_over_bstrw[0,0]: ', bstrc_over_bstrw[0,0])
print('bstrc_over_bstrw[-1.-1]: ', bstrc_over_bstrw[-1,-1])
# Get rid of infinity 
no_inf1 = np.where(bstrc_over_bstrw_wland == np.inf, 10**6, bstrc_over_bstrw_wland)
# Get rid of nan
no_inf2 = np.nan_to_num(no_inf1, nan=10**-6)
print('no_inf2[0,0]: ', no_inf2[0,0])
print('no_inf2[-1.-1]: ', no_inf2[-1,-1]) 
#bstrc_over_bstrw_wland = np.nan_to_num(bstrc_over_bstrw_wland, nan=10**-6)
#bstrc_over_bstrw_wland = bstrc_over_bstrw_wland.fillna(-1)

# Set the colormap
#cmap = matplotlib.cm.cividis
cmap4=cmocean.cm.matter
cmap4.set_under('darkgray')
cmap4.set_bad('darkgray')

# Set the levels 
print(np.nanmax(bstrc_over_bstrw), flush=True)
#lev5  = np.arange(0, 10, 0.1)
#lev5 = np.arange(10**(-4), 10**4, 0.1)
lev5 = 10.**np.arange(-3,4,1)

# Plot the bathymetry in 10m contours from 10 - 100 m
ax1[3].contour(model_output.lon_rho.values, model_output.lat_rho.values,
             model_output.bath[0,:,:].values, lev2, colors='gray')

#cs4 = ax1[3].contourf(model_output.lon_rho[:,:].values, model_output.lat_rho[:,:].values, 
#                  bstrc_over_bstrw_wland[:,:], locator=LogLocator(), cmap=cmap4, extend='both')

# TEMP
cs4 = ax1[3].contourf(model_output.lon_rho[:,:].values, model_output.lat_rho[:,:].values, 
                  no_inf2[:,:], lev5, locator=LogLocator(), cmap=cmap4, extend='both')

# Label the plot
#ax1[3].set_title('Ratio of Current Over Wave Stress', fontsize=fontsize, y=1.08)
ax1[3].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[3].set_ylabel('Latitude \n(degrees)', fontsize=fontsize)
#cbar4_ax = fig1.add_axes(top=0.25, bottom=0.1, left=0.1, right=0.8)
cbar4_ax = fig1.add_axes([0.82, 0.1, 0.02, 0.19]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar4 = plt.colorbar(cs4, orientation='vertical', format=LogFormatterSciNotation(), cax=cbar4_ax, ax=ax1[3]).set_label(label='Ratio of Current \nStress Over Wave \nStress', size=fontsize)
#cbar4.formatter

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig1.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
cbar1_ax = fig1.add_axes([0.82, 0.32, 0.02, 0.63])
cbar3 = plt.colorbar(cs3, orientation='vertical', cax=cbar1_ax, ax=[ax1[0], ax1[1], ax1[2]]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08)

# Add labels
plt.text(0.770, 0.927, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.770, 0.710, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.770, 0.494, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.770, 0.278, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the plot
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/time_avg_current_wave_combo_ratio_bedstress_nomerge_dbsed0007_0001.png') 

# Check white areas
print('bstrc_over_bstrw_wland[0,0]: ', bstrc_over_bstrw_wland[0,0])
print('bstrc_over_bstrw_wland[-1.-1]: ', bstrc_over_bstrw_wland[-1,-1])


# -------------------------------------------------------------------
# ---------- Plot 2: Time-Averaged Wave-Current Bed Stress ----------
# ------------------------ Masked -----------------------------------
# -------------------------------------------------------------------
# Make a plot of the time-averaged bed shear stress from both waves
# and currents.

# Make it so land will appear
temp_mask = model_output.mask_rho.copy()
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)

# Make it so land appears 
bstrcwmax_total_wland = bstrcwmax_total * temp_mask
bstrcwmax_total_wland = np.nan_to_num(bstrcwmax_total_wland, nan=-5000000)

# Set the number of cells in the sponge on each open boundary
c_west = 36
c_north = 45
c_east = 36

# Prep the data by  ultiplying by the mask and trimming
# Multiply by mask
h_masked = grid.h.values*grid.mask_rho.values*mask_rho_nan.nudge_mask_rho_nan
# Trim 
lon_rho_trimmed = grid.lon_rho[:,c_west:-c_west].values
lat_rho_trimmed = grid.lat_rho[:,c_west:-c_west].values
h_masked_trimmed = h_masked[:,c_west:-c_west]

# Mask, trim, slice
bstrcwmax_total_wland_masked = bstrcwmax_total_wland*mask_rho_nan.nudge_mask_rho_nan
bstrcwmax_total_wland_masked_trimmed = bstrcwmax_total_wland_masked[:,c_west:-c_west]

# Make the figure
fig2, ax2 = plt.subplots(figsize=(26, 12))  # (18,16) (16,18) (16,21)

# Plot current-induced
# Set the colormap
#cmap = matplotlib.cm.cividis
cmap2=cmocean.cm.amp
cmap2.set_under('darkgray')
#cmap2.set_bad('darkgray')

# Set a colormap and levels to use for all plots
cmap2b = cmocean.cm.tempo
cmap2b.set_under('darkgray')
lev2b = np.arange(0-0.000001, 0.25, 0.01)


# Plot the bathymetry in 10m contours from 10 - 100 m
lev3 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax2.contour(lon_rho_trimmed, lat_rho_trimmed,
             h_masked_trimmed, lev3, colors='gray')

cs5 = ax2.contourf(lon_rho_trimmed, lat_rho_trimmed,
                  bstrcwmax_total_wland_masked_trimmed, lev2b, cmap=cmap2b, extend='both')

# Label the plot
#ax1[0].set_title('Time-Averaged Current-Induced Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
#plt.setp(ax1[0].get_xticklabels(), visible=False)
ax2.set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax2.set_ylabel('Latitude \n(degrees)', fontsize=fontsize)
#cbar1 = plt.colorbar(cs1, orientation='vertical', ax=ax1[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


cbar2 = plt.colorbar(cs5, orientation='vertical').set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)
#cbar4.formatter

# Save the plot
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/time_avg_bstrcwmax_masked_dbsed0007_0001.png')



# -------------------------------------------------------------------
# ---------- Plot 3: Time-Averaged Wave-Current Bed Stress ----------
# ---------------------- Masked XY ----------------------------------
# -------------------------------------------------------------------
# Make a plot of the time-averaged bed shear stress from both waves
# and currents.

# Make it so land will appear
temp_mask = model_output.mask_rho.copy()
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)

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

# Prep the data by  ultiplying by the mask and trimming
# Multiply by mask
h_masked = grid.h.values*grid.mask_rho.values*mask_rho_nan.nudge_mask_rho_nan
# Trim 
lon_rho_trimmed = grid.lon_rho[:,c_west:-c_west].values
lat_rho_trimmed = grid.lat_rho[:,c_west:-c_west].values
h_masked_trimmed = h_masked[:,c_west:-c_west]

# Make the figure
fig3, ax3 = plt.subplots(figsize=(26, 8))  # (18,16) (16,18) (16,21)

# Set a colormap and levels to use for all plots
cmap3b = cmocean.cm.tempo
cmap3b.set_under('darkgray')
lev3b = np.arange(0-0.000001, 0.25, 0.01)

# Plot the bathymetry in 10m contours from 10 - 100 m
lev4 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax3.contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev3, colors='gray')

cs6 = ax3.contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  bstrcwmax_total_wland_masked_trimmed, lev3b, cmap=cmap3b, extend='both')

# Label the plot
#ax1[0].set_title('Time-Averaged Current-Induced Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
#plt.setp(ax1[0].get_xticklabels(), visible=False)
ax3.set_xlabel('X (km)', fontsize=fontsize)
ax3.set_ylabel('Y (km)', fontsize=fontsize)
#cbar1 = plt.colorbar(cs1, orientation='vertical', ax=ax1[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


cbar3 = plt.colorbar(cs6, orientation='vertical').set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)
#cbar4.formatter

# Save the plot
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/time_avg_bstrcwmax_masked_xy_dbsed0007_0001.png')



# -------------------------------------------------------------------
# ---------- Plot 4: Time-Averaged Bed Stress ----------
# ------------- Wave, Current, Combo, Ratio Masked XY ---------------
# -------------------------------------------------------------------
# Make a plot of the time-averaged bed shear stress
# with subplots for current-induced, waves-induced,
# the two combined, and the ratio of the two. Since this 
# has both u and v components, the total 
# magnitude on the rho points must be calculated. Do this but masked and
# in xy 


# Mask, trim
bstrc_tot_mag_wland_masked = bstrc_tot_mag_wland*mask_rho_nan.nudge_mask_rho_nan*temp_mask
bstrw_tot_mag_wland_masked = bstrw_tot_mag_wland*mask_rho_nan.nudge_mask_rho_nan*temp_mask
bstrcwmax_total_wland_masked = bstrcwmax_total*mask_rho_nan.nudge_mask_rho_nan*temp_mask
bstrc_tot_mag_wland_masked_trimmed = bstrc_tot_mag_wland_masked[:,c_west:-c_west]
bstrw_tot_mag_wland_masked_trimmed = bstrw_tot_mag_wland_masked[:,c_west:-c_west]
bstrcwmax_total_wland_masked_trimmed = bstrcwmax_total_wland_masked[:,c_west:-c_west]


# Make the figure
fig4, ax4 = plt.subplots(4, figsize=(22, 21))  # (18,16) (16,18) (16,21) (18,25)

# Plot current-induced
# Set the colormap
cmap4=cmocean.cm.amp

# Set a colormap and levels to use for all plots
cmap4b = cmocean.cm.tempo
cmap4b.set_under('darkgray')
lev4b = np.arange(0.0, 0.25, 0.01)


# Set the levels 
lev5  = np.arange(-0.00001, np.nanmax(bstrc_tot_mag), 0.01)
# Shade background 
ax4[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[0].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
lev6 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax4[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')
cs7 = ax4[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                  bstrc_tot_mag_wland_masked_trimmed, lev4b, cmap=cmap4b, extend='max')
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517-36
ax4[0].scatter(x_rho_flat_trimmed[xi_kakt_idx]/1000, y_rho_flat[eta_kakt_idx]/1000, 
            marker='.',  s=900, linewidth=4, color='brown', label='Kaktovik')
ax4[0].legend(fontsize=fontsize)

# Label the plot
#ax1[0].set_title('Time-Averaged Current-Induced Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
plt.setp(ax4[0].get_xticklabels(), visible=False)
#ax1[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax4[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar1 = plt.colorbar(cs1, orientation='vertical', ax=ax1[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


# Plot  wave-induced
# Shade background 
ax4[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[1].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax4[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')

cs8 = ax4[1].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                  bstrw_tot_mag_wland_masked_trimmed, lev4b, cmap=cmap4b, extend='max')

# Label the plot
#ax1[1].set_title('Time-Averaged Wave-Induced Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
plt.setp(ax4[1].get_xticklabels(), visible=False)
#ax1[1].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax4[1].set_ylabel('Y (km)', fontsize=fontsize)
#cbar2 = plt.colorbar(cs2, orientation='vertical', ax=ax1[1]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


# Plot  wave and current induced
# Shade background 
ax4[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[2].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax4[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')

cs9 = ax4[2].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 bstrcwmax_total_wland_masked_trimmed, lev4b, cmap=cmap4b, extend='max')

# Label the plot
#ax1[2].set_title('Time-Averaged Total Bed Stress (N/m\u00b2)', fontsize=fontsize, y=1.08)
plt.setp(ax4[2].get_xticklabels(), visible=False)
#ax1[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax4[2].set_ylabel('Y (km)', fontsize=fontsize)

# Plot ratio of current over wave
# Take the ratio
bstrc_over_bstrw = bstrc_tot_mag/bstrw_tot_mag

# Replace the nans with -100 to help with plotting out of bounds data
bstrc_over_bstrw_wland_masked = bstrc_over_bstrw * temp_mask * mask_rho_nan.nudge_mask_rho_nan
bstrc_over_bstrw_wland_masked_trimmed = bstrc_over_bstrw_wland_masked[:,c_west:-c_west]

# Get rid of infinity 
no_inf1_masked_trimmed = np.where(bstrc_over_bstrw_wland_masked_trimmed == np.inf, 10**6, bstrc_over_bstrw_wland_masked_trimmed)
# Get rid of nan
no_inf2 = np.nan_to_num(no_inf1, nan=10**-6)
no_inf2_masked_trimmed = np.nan_to_num(no_inf1_masked_trimmed, nan=10**-6)
print('no_inf2[0,0]: ', no_inf2[0,0])
print('no_inf2[-1.-1]: ', no_inf2[-1,-1]) 
#bstrc_over_bstrw_wland = np.nan_to_num(bstrc_over_bstrw_wland, nan=10**-6)
#bstrc_over_bstrw_wland = bstrc_over_bstrw_wland.fillna(-1)

# Set the colormap
#cmap = matplotlib.cm.cividis
cmap5=cmocean.cm.balance


# Set the levels 
lev5 = 10.**np.arange(-3,4,1)
# Shade background 
ax4[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[3].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax4[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')

# Plot ratio
cs10 = ax4[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                  no_inf2_masked_trimmed[:,:], lev5, locator=LogLocator(), cmap=cmap5, extend='max')

# Label the plot
#ax1[3].set_title('Ratio of Current Over Wave Stress', fontsize=fontsize, y=1.08)
ax4[3].set_xlabel('X (km)', fontsize=fontsize)
ax4[3].set_ylabel('Y (km)', fontsize=fontsize)
#cbar4_ax = fig1.add_axes(top=0.25, bottom=0.1, left=0.1, right=0.8)
cbar10_ax = fig4.add_axes([0.82, 0.1, 0.02, 0.19]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar10 = plt.colorbar(cs10, orientation='vertical', format=LogFormatterSciNotation(), cax=cbar10_ax, ax=ax4[3]).set_label(label='Ratio of Current Stress \nOver Wave Stress', size=fontsize)
#cbar4.formatter

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig4.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
cbar9_ax = fig4.add_axes([0.82, 0.32, 0.02, 0.63])
cbar9 = plt.colorbar(cs9, orientation='vertical', cax=cbar9_ax, ax=[ax4[0], ax4[1], ax4[2]]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08)

# Add labels
# Top right 
#plt.text(0.770, 0.927, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.770, 0.710, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.770, 0.494, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.770, 0.278, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Bottom right 
plt.text(0.770, 0.757, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.770, 0.540, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.770, 0.324, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.770, 0.108, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the plot
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/time_avg_current_wave_combo_bedstress_masked_xy_dbsed0007_0001.png')




# -------------------------------------------------------------------
# ---------- Plot 5: Time-Averaged Bed Stress Combo ----------
# ------------- & Ratios Masked XY ---------------
# -------------------------------------------------------------------
# Make a plot of the time-averaged bed shear stress
# with subplots for total bed stress, and the ratio of the eaach
# over the total. Since this 
# has both u and v components, the total 
# magnitude on the rho points must be calculated. Do this but masked and
# in xy 

# Make the figure
fig5, ax5 = plt.subplots(3, figsize=(22, 21))  # (18,16) (16,18) (16,21) (18,25)

# Plot total bed stress
# Set the colormap
cmap5=cmocean.cm.amp
# Set a colormap and levels to use for all plots
cmap5b = cmocean.cm.tempo
cmap5b.set_under('darkgray')
lev5b = np.arange(0.0, 0.25, 0.01)
# Set the levels 
lev8  = np.arange(-0.00001, np.nanmax(bstrc_tot_mag), 0.01)
# Shade background 
ax5[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[0].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
lev9 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax5[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev9, colors='darkorchid')
cs11 = ax5[0].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 bstrcwmax_total_wland_masked_trimmed, lev5b, cmap=cmap5b, extend='max')
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517-36
ax5[0].scatter(x_rho_flat_trimmed[xi_kakt_idx]/1000, y_rho_flat[eta_kakt_idx]/1000, 
            marker='.',  s=900, linewidth=4, color='brown', label='Kaktovik')
ax5[0].legend(fontsize=fontsize)
# Label the plot
plt.setp(ax5[0].get_xticklabels(), visible=False)
ax5[0].set_ylabel('Y (km)', fontsize=fontsize)
cbar11_ax = fig5.add_axes([0.82, 0.39, 0.02, 0.55]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar11 = plt.colorbar(cs11, orientation='vertical', cax=cbar11_ax, ax=ax5[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize, labelpad=15)

# Set the colormap
#cmap = matplotlib.cm.cividis
cmap7=cmocean.cm.matter
# Set the levels 
lev10 = np.arange(0,100,0.1)

# Plot  current-induced/total
# Shade background 
ax5[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[1].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax5[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')

cs12 = ax5[1].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                  (bstrc_tot_mag_wland_masked_trimmed/bstrcwmax_total_wland_masked_trimmed)*100, 
                  lev10, cmap=cmap7, extend='max')
# Label the plot
plt.setp(ax5[1].get_xticklabels(), visible=False)
ax5[1].set_ylabel('Y (km)', fontsize=fontsize)


# Plot wave-induced/total
# Shade background 
ax5[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[2].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax5[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')

# Plot ratio
cs13 = ax5[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                  (bstrw_tot_mag_wland_masked_trimmed/bstrcwmax_total_wland_masked_trimmed)*100, 
                  lev10, cmap=cmap7, extend='max') 

# Label the plot
ax5[2].set_xlabel('X (km)', fontsize=fontsize)
ax5[2].set_ylabel('Y (km)', fontsize=fontsize)
# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig5.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
cbar13_ax = fig5.add_axes([0.82, 0.1, 0.02, 0.27]) #[left, bottom, width, height]
cbar13 = plt.colorbar(cs13, orientation='vertical', cax=cbar13_ax, ax=[ax5[1], ax5[2]]).set_label(label='Percent of \nTotal Stress (%)', size=fontsize, labelpad=15)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08)

# Add labels
# Bottom right 
plt.text(0.773, 0.689, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.399, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.113, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Label subplots
plt.text(0.470, 0.620, 'Current Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.470, 0.335, 'Wave Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# -------------------------------------------------------------------
# ---------- Plot 6: Time-Averaged Bed Stress Combo ----------
# ----- & Ratios with Frequency of Resuspension Masked XY -----------
# -------------------------------------------------------------------
# Same as above but with a fourth panel for the frequency of resuspension for
# the critical shear stress for muds and fine sand

# First find the frequency of resuspension
# Load in the rho masks 
#mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
#mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')

# Pull out the critical shear stresses 
# Pull form sediment.in for now since tau_ce in ocean_his has different units 
# than bustr 
tau_crit_mud01 = 0.18 # N/m2 , 0.05
tau_crit_mud02 = 0.18 # N/m2 , 0.05
tau_crit_sand01 = 0.18 # N/m2 , 0.1
tau_crit_sand02 = 0.34 # N/m2 , 0.24
tau_crit_sand03 = 24 # N/m2 # 3783


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


# Make a function to count the time steps above critical
def cnt_time_steps_above_crit(filename, tau_crit_mud01, tau_crit_mud02, tau_crit_sand01, tau_crit_sand02, tau_crit_sand03):
    """
    This function takes a given model output and 
    calculates the number of time steps where the 
    bed shear stress is above the critical shear stress
    for erosion for each sediment class. This function 
    looks at bstrcwmax to get at wave and current 
    induced bed stress.
     
    Inputs:
    - filename: string, path to and name of netcdf output
    - tau_crit_mud01: critical shear stress for mud01 (N/m2)
    - tau_crit_mud02: critical shear stress for mud02 (N/m2)
    - tau_crit_sand01: critical shear stress for sand01 (N/m2)
    - tau_crit_sand02: critical shear stress for sand02 (N/m2)
    - tau_crit_sand03: critical shear stress for sand03 (N/m2)
    
    Outputs:
    - time_above_mud01_cnt_tmp: number of time steps where bed 
       shear stress is above critical value for mud01
    - time_above_mud02_cnt_tmp: number of time steps where bed 
       shear stress is above critical value for mud02
    - time_above_sand01_cnt_tmp: number of time steps where bed 
       shear stress is above critical value for sand01
    - time_above_sand02_cnt_tmp: number of time steps where bed 
       shear stress is above critical value for sand02
    - time_above_sand03_cnt_tmp: number of time steps where bed 
       shear stress is above critical value for sand03
    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Pull out total bottom shear stress
    tot_bottom_stress = model_output.bstrcwmax.copy()
    
    # Save the length of time of this model output
    time_len_tmp = len(model_output.ocean_time)
    
    # Find the number of time steps where the bed stress is above the critical
    # stress for each sediment class
    time_above_mud01_cnt_tmp = (tot_bottom_stress.where(tot_bottom_stress > tau_crit_mud01).groupby('eta_rho',).count(dim='ocean_time'))
    time_above_mud02_cnt_tmp = (tot_bottom_stress.where(tot_bottom_stress > tau_crit_mud02).groupby('eta_rho',).count(dim='ocean_time'))
    time_above_sand01_cnt_tmp = (tot_bottom_stress.where(tot_bottom_stress > tau_crit_sand01).groupby('eta_rho',).count(dim='ocean_time'))
    time_above_sand02_cnt_tmp = (tot_bottom_stress.where(tot_bottom_stress > tau_crit_sand02).groupby('eta_rho',).count(dim='ocean_time'))
    time_above_sand03_cnt_tmp = (tot_bottom_stress.where(tot_bottom_stress > tau_crit_sand03).groupby('eta_rho',).count(dim='ocean_time'))
    
    # Pull out the time for this run
    time_tmp = model_output.ocean_time.copy()
    
    # Return each of the counts
    return(time_len_tmp, time_above_mud01_cnt_tmp, time_above_mud02_cnt_tmp, time_above_sand01_cnt_tmp, time_above_sand02_cnt_tmp, time_above_sand03_cnt_tmp)


# Make some arrays to hold output
time_lengths = np.empty((num_files))
time_above_mud01_cnt = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_mud02_cnt = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand01_cnt = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand02_cnt = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand03_cnt = np.empty((num_files, eta_rho_len, xi_rho_len))

# Loop through the model output
for j in range(num_files):
        
    # Call the function to open the output
    time_len_tmp, time_above_mud01_tmp, time_above_mud02_tmp, time_above_sand01_tmp, time_above_sand02_tmp, time_above_sand03_tmp = cnt_time_steps_above_crit(file_names2[j], tau_crit_mud01, tau_crit_mud02, tau_crit_sand01, tau_crit_sand02, tau_crit_sand03)
    
    # Save these to the arrays
    time_lengths[j] = time_len_tmp
    time_above_mud01_cnt[j,:,:] = time_above_mud01_tmp
    time_above_mud02_cnt[j,:,:] = time_above_mud02_tmp
    time_above_sand01_cnt[j,:,:] = time_above_sand01_tmp
    time_above_sand02_cnt[j,:,:] = time_above_sand02_tmp
    time_above_sand03_cnt[j,:,:] = time_above_sand03_tmp
    
    
# Now that the counts have been found for all outputs,
# combine them into one percent for each location in grid
# Find the total length of time steps
total_time_len = np.sum(time_lengths, axis=0)
time_above_mud01 = (np.sum(time_above_mud01_cnt, axis=0)/total_time_len)
time_above_mud02 = (np.sum(time_above_mud02_cnt, axis=0)/total_time_len)
time_above_sand01 = (np.sum(time_above_sand01_cnt, axis=0)/total_time_len)
time_above_sand02 = (np.sum(time_above_sand02_cnt, axis=0)/total_time_len)
time_above_sand03 = (np.sum(time_above_sand03_cnt, axis=0)/total_time_len)


# Check shapes to see if this worked
print('time_above_mud01_cnt shape: ', np.shape(time_above_mud01_cnt), flush=True)
print('time_above_mud01 shape: ', np.shape(time_above_mud01), flush=True)
print('total_time_len: ', total_time_len, flush=True)

# Load in the last model output to use for plotting things
model_output = xr.open_dataset(file_names2[-1])


# Print these value to the output file
print('Mud01 min percent: ', time_above_mud01.min(), '; Mud01 max percent: ', time_above_mud01.max(), flush=True)
print('Mud02 min percent: ', time_above_mud02.min(), '; Mud02 max percent: ', time_above_mud02.max(), flush=True)
print('Sand01 min percent: ', time_above_sand01.min(), '; Sand01 max percent: ', time_above_sand01.max(), flush=True)
print('Sand02 min percent: ', time_above_sand02.min(), '; Sand02 max percent: ', time_above_sand02.max(), flush=True)
print('Sand03 min percent: ', time_above_sand03.min(), '; Sand03 max percent: ', time_above_sand03.max(), flush=True)





# Make the figure
fig6, ax6 = plt.subplots(4, figsize=(22, 21))  # (18,16) (16,18) (16,21) (18,25)

# Plot total bed stress
# Set the colormap
cmap8=cmocean.cm.amp
# Set a colormap and levels to use for all plots
cmap8b = cmocean.cm.tempo
cmap8b.set_under('darkgray')
lev8b = np.arange(0.0, 0.25, 0.01)
# Set the levels 
lev11  = np.arange(-0.00001, np.nanmax(bstrc_tot_mag), 0.01)
# Shade background 
ax6[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[0].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
lev9 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax6[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev9, colors='darkorchid')
cs14 = ax6[0].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 bstrcwmax_total_wland_masked_trimmed, lev8b, cmap=cmap5b, extend='max')
# Kaktovik
eta_kakt_idx = 58
xi_kakt_idx = 517-36
ax6[0].scatter(x_rho_flat_trimmed[xi_kakt_idx]/1000, y_rho_flat[eta_kakt_idx]/1000, 
            marker='.',  s=900, linewidth=4, color='brown', label='Kaktovik')
ax6[0].legend(fontsize=fontsize)
# Label the plot
plt.setp(ax6[0].get_xticklabels(), visible=False)
ax6[0].set_ylabel('Y (km)', fontsize=fontsize)
cbar14_ax = fig6.add_axes([0.82, 0.75, 0.02, 0.2]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar14 = plt.colorbar(cs14, orientation='vertical', cax=cbar14_ax, ax=ax6[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize, labelpad=15)
# Set the colormap
#cmap = matplotlib.cm.cividis
cmap9=cmocean.cm.dense
# Set the levels 
#lev12 = np.arange(0,100,0.1)
lev12 = np.arange(0,50,0.1)

# Plot  current-induced/total
# Shade background 
ax6[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[1].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax6[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev9, colors='deeppink')

cs15 = ax6[1].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                  (bstrc_tot_mag_wland_masked_trimmed/bstrcwmax_total_wland_masked_trimmed)*100, 
                  lev12, cmap=cmap9, extend='max')
# Label the plot
plt.setp(ax6[1].get_xticklabels(), visible=False)
ax6[1].set_ylabel('Y (km)', fontsize=fontsize)


# Plot wave-induced/total
# Shade background 
ax6[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[2].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax6[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='deeppink')

# Plot ratio
cs16 = ax6[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                  (bstrw_tot_mag_wland_masked_trimmed/bstrcwmax_total_wland_masked_trimmed)*100, 
                  lev12, cmap=cmap9, extend='max') 

# Label the plot
#ax6[2].set_xlabel('X (km)', fontsize=fontsize)
ax6[2].set_ylabel('Y (km)', fontsize=fontsize)
plt.setp(ax6[2].get_xticklabels(), visible=False)
# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig6.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
cbar15_ax = fig6.add_axes([0.82, 0.32, 0.02, 0.42]) #[left, bottom, width, height]
cbar15 = plt.colorbar(cs15, orientation='vertical', cax=cbar15_ax, ax=[ax6[1], ax6[2]]).set_label(label='Percent of \nTotal Stress (%)', size=fontsize, labelpad=15)


# Plot frequency of resuspension
# Make it so land appears 
time_above_mud01_wland = time_above_mud01 * temp_mask
time_above_mud01_wland = np.nan_to_num(time_above_mud01_wland, nan=-5000000)
# Prep the data by  ultiplying by the mask and trimming
# Mask, trim, slicei
time_above_mud01_wland_masked = time_above_mud01*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_wland_masked_trimmed = time_above_mud01_wland_masked[:,c_west:-c_west]
# Set the colormap
cmap10=cmocean.cm.matter
lev13 = np.arange(0, 0.75, 0.01)
# Tau_crit = 0.18 (mud01, mud02, sand01)
# Whole grid
# Make it so land will be gray and ocean white 
ax6[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6[3].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax6[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')
cs16 = ax6[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                     time_above_mud01_wland_masked_trimmed, lev13, cmap=cmap10, extend='max')
# Sagavanirktok River
eta_sag_idx = 37 #36
xi_sag_idx = 279-36
ax6[3].scatter(x_rho_flat_trimmed[xi_sag_idx]/1000, y_rho_flat[eta_sag_idx]/1000, 
            marker='.', s=2100, linewidth=5, color='deepskyblue', label='Sagavanirktok')
#ax6[3].set_title('Tau Critical: ' +str(tau_crit_mud01) + ' (Pa)' , size=30)
ax6[3].set_xlabel('X (km)', size=fontsize)
ax6[3].set_ylabel('Y (km)', size=fontsize)
cbar16_ax = fig6.add_axes([0.82, 0.1, 0.02, 0.2]) #[left, bottom, width, height]
cbar16 = plt.colorbar(cs16, orientation='vertical', cax=cbar16_ax, ax=ax6[3]).set_label(label='Fraction of Time', size=fontsize, labelpad=15)
ax6[3].legend(fontsize=fontsize, loc='upper right')


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08)

# Add labels
# Bottom right 
plt.text(0.773, 0.760, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.542, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.329, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.113, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Label subplots
plt.text(0.470, 0.700, 'Current Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.470, 0.487, 'Wave Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.393, 0.268, 'Tau Critical: ' +str(tau_crit_mud01) + ' (Pa)', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# Make the figure again but slightly different 
fig7, ax7 = plt.subplots(4, figsize=(22, 21))  # (18,16) (16,18) (16,21) (18,25)

# Plot total bed stress
# Set the colormap
cmap8=cmocean.cm.amp
# Set a colormap and levels to use for all plots
cmap8b = cmocean.cm.tempo
cmap8b.set_under('darkgray')
lev8b = np.arange(0.0, 0.25, 0.01)
# Set the levels 
lev11  = np.arange(-0.00001, np.nanmax(bstrc_tot_mag), 0.01)
# Shade background 
ax7[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[0].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
lev9 = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
ax7[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev9, colors='darkorchid')
cs14 = ax7[0].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 bstrcwmax_total_wland_masked_trimmed, lev8b, cmap=cmap5b, extend='max')
# Kaktovik
#eta_kakt_idx = 58
#xi_kakt_idx = 517-36
#ax7[0].scatter(x_rho_flat_trimmed[xi_kakt_idx]/1000, y_rho_flat[eta_kakt_idx]/1000, 
 #           marker='.',  s=900, linewidth=4, color='brown', label='Kaktovik')
#ax7[0].legend(fontsize=fontsize)
# Label the plot
plt.setp(ax7[0].get_xticklabels(), visible=False)
ax7[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar14_ax = fig7.add_axes([0.82, 0.75, 0.02, 0.2]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
#cbar14 = plt.colorbar(cs14, orientation='vertical', cax=cbar14_ax, ax=ax7[0]).set_label(label='Bed Shear \nStress (N/m\u00b2)', size=fontsize, labelpad=15)
# If horizontal and in axes
cbar14 = plt.colorbar(cs14, cax=ax7[0].inset_axes((0.48, 0.92, 0.5, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 0.05, 0.1, 0.15, 0.2], ax=ax7[0], 
                     orientation='horizontal').set_label(label='Bed Shear Stress (N/m\u00b2)', size=fontsize-2)
# Set the colormap
#cmap = matplotlib.cm.cividis
cmap9=cmocean.cm.dense
# Set the levels 
#lev12 = np.arange(0,100,0.1)
lev12 = np.arange(0,50,0.1)

# Plot  current-induced/total
# Shade background 
ax7[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[1].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax7[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev9, colors='deeppink')

cs15 = ax7[1].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                  (bstrc_tot_mag_wland_masked_trimmed/bstrcwmax_total_wland_masked_trimmed)*100, 
                  lev12, cmap=cmap9, extend='max')
# Label the plot
plt.setp(ax7[1].get_xticklabels(), visible=False)
ax7[1].set_ylabel('Y (km)', fontsize=fontsize)


# Plot wave-induced/total
# Shade background 
ax7[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[2].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax7[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='deeppink')

# Plot ratio
cs16 = ax7[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                  (bstrw_tot_mag_wland_masked_trimmed/bstrcwmax_total_wland_masked_trimmed)*100, 
                  lev12, cmap=cmap9, extend='max') 

# Label the plot
#ax6[2].set_xlabel('X (km)', fontsize=fontsize)
ax7[2].set_ylabel('Y (km)', fontsize=fontsize)
plt.setp(ax7[2].get_xticklabels(), visible=False)
# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig7.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
#cbar15_ax = fig7.add_axes([0.82, 0.32, 0.02, 0.42]) #[left, bottom, width, height]
#cbar15 = plt.colorbar(cs15, orientation='vertical', cax=cbar15_ax, ax=[ax7[1], ax7[2]]).set_label(label='Percent of \nTotal Stress (%)', size=fontsize, labelpad=15)
# If horizontal and in axes
cbar15a = plt.colorbar(cs15, cax=ax7[1].inset_axes((0.48, 0.92, 0.5, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 10, 20, 30, 40], ax=ax7[1], 
                     orientation='horizontal').set_label(label='% Current Stress', size=fontsize-2)
# If horizontal and in axes
cbar15b = plt.colorbar(cs15, cax=ax7[2].inset_axes((0.48, 0.92, 0.5, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 10, 20, 30, 40], ax=ax7[2], 
                     orientation='horizontal').set_label(label='% Wave Stress', size=fontsize-2)

# Plot frequency of resuspension
# Make it so land appears 
time_above_mud01_wland = time_above_mud01 * temp_mask
time_above_mud01_wland = np.nan_to_num(time_above_mud01_wland, nan=-5000000)
# Prep the data by  ultiplying by the mask and trimming
# Mask, trim, slicei
time_above_mud01_wland_masked = time_above_mud01*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_wland_masked_trimmed = time_above_mud01_wland_masked[:,c_west:-c_west]
# Set the colormap
cmap10=cmocean.cm.matter
lev13 = np.arange(0, 0.75, 0.01)
# Tau_crit = 0.18 (mud01, mud02, sand01)
# Whole grid
# Make it so land will be gray and ocean white 
ax7[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7[3].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
# Plot the bathymetry in 10m contours from 10 - 100 m
ax7[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000,
             h_masked_trimmed, lev6, colors='darkorchid')
cs16 = ax7[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                     time_above_mud01_wland_masked_trimmed, lev13, cmap=cmap10, extend='max')
# Sagavanirktok River
#eta_sag_idx = 37 #36
#xi_sag_idx = 279-36
#ax7[3].scatter(x_rho_flat_trimmed[xi_sag_idx]/1000, y_rho_flat[eta_sag_idx]/1000, 
 #           marker='.', s=2100, linewidth=5, color='deepskyblue', label='Sagavanirktok')
#ax6[3].set_title('Tau Critical: ' +str(tau_crit_mud01) + ' (Pa)' , size=30)
ax7[3].set_xlabel('X (km)', size=fontsize)
ax7[3].set_ylabel('Y (km)', size=fontsize)
#cbar16_ax = fig7.add_axes([0.82, 0.1, 0.02, 0.2]) #[left, bottom, width, height]
#cbar16 = plt.colorbar(cs16, orientation='vertical', cax=cbar16_ax, ax=ax7[3]).set_label(label='Fraction of Time', size=fontsize, labelpad=15)
#ax7[3].legend(fontsize=fontsize, loc='upper right')
# If horizontal and in axes
cbar16 = plt.colorbar(cs16, cax=ax7[3].inset_axes((0.48, 0.92, 0.50, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, .15, .30, .45, .60], ax=ax7[3], 
                     orientation='horizontal').set_label(label='Fraction of Time', size=fontsize-2)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08)

# Add labels
# Bottom right 
plt.text(0.773, 0.760, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.542, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.329, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.773, 0.113, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)








# Make some calculations and print them
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
# 10 - 20 m depth
h_masked4b = h_masked.copy()
inner_10_20m_mask_rho = masked_array_lowhigh_2dloop(h_masked4b, 10, 20)
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
# 10 - 20 m depth
inner_10_20m_mask_rho_nan_idx = np.where(inner_10_20m_mask_rho == 0.0)
inner_10_20m_mask_rho_nan = inner_10_20m_mask_rho.copy()
inner_10_20m_mask_rho_nan = inner_10_20m_mask_rho_nan.astype('float')
inner_10_20m_mask_rho_nan[inner_10_20m_mask_rho_nan_idx] = np.nan
# 10 - 60 m depth
outer_10_60m_mask_rho_nan_idx = np.where(outer_10_60m_mask_rho == 0.0)
outer_10_60m_mask_rho_nan = outer_10_60m_mask_rho.copy()
outer_10_60m_mask_rho_nan = outer_10_60m_mask_rho_nan.astype('float')
outer_10_60m_mask_rho_nan[outer_10_60m_mask_rho_nan_idx] = np.nan


# Now multiply by the mask to get the different regions 
# For bed shear stress (for now)
# Inner 
bstrcwmax_total_masked_inner = bstrcwmax_total * inner_shelf_mask_rho_nan
# Mid
bstrcwmax_total_masked_mid = bstrcwmax_total * mid_shelf_mask_rho_nan
# Outer
bstrcwmax_total_masked_outer = bstrcwmax_total * outer_shelf_mask_rho_nan
# 0 - 10 m
bstrcwmax_total_masked_10m = bstrcwmax_total * inner_10m_mask_rho_nan
# 10 - 20 m
bstrcwmax_total_masked_10_20m = bstrcwmax_total * inner_10_20m_mask_rho_nan
# 10 - 60 m
bstrcwmax_total_masked_10_60m = bstrcwmax_total * outer_10_60m_mask_rho_nan
# For frequency of resuspension
# Inner 
time_above_mud01_masked_inner = time_above_mud01 * inner_shelf_mask_rho_nan
# Mid
time_above_mud01_masked_mid = time_above_mud01 * mid_shelf_mask_rho_nan
# Outer
time_above_mud01_masked_outer = time_above_mud01 * outer_shelf_mask_rho_nan
# 0 - 10 m
time_above_mud01_masked_10m = time_above_mud01 * inner_10m_mask_rho_nan
# 10 - 20 m
time_above_mud01_masked_10_20m = time_above_mud01 * inner_10_20m_mask_rho_nan
# 10 - 60 m
time_above_mud01_masked_10_60m = time_above_mud01 * outer_10_60m_mask_rho_nan

# Mask and trim these so that we are just looking at the regions in the plot
# Mask, trim, slice
# Mask
# Inner
bstrcwmax_total_masked_inner_masked = bstrcwmax_total_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_masked_inner_masked = time_above_mud01_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Mid
bstrcwmax_total_masked_mid_masked = bstrcwmax_total_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_masked_mid_masked = time_above_mud01_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Outer
bstrcwmax_total_masked_outer_masked = bstrcwmax_total_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_masked_outer_masked = time_above_mud01_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 0 - 10 m
bstrcwmax_total_masked_10m_masked = bstrcwmax_total_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_masked_10m_masked = time_above_mud01_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 20 m
bstrcwmax_total_masked_10_20m_masked = bstrcwmax_total_masked_10_20m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_masked_10_20m_masked = time_above_mud01_masked_10_20m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 60 m
bstrcwmax_total_masked_10_60m_masked = bstrcwmax_total_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
time_above_mud01_masked_10_60m_masked = time_above_mud01_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask

# Trim
# Inner
bstrcwmax_total_masked_inner_masked_trimmed = bstrcwmax_total_masked_inner_masked[:,c_west:-c_west]
time_above_mud01_masked_inner_masked_trimmed = time_above_mud01_masked_inner_masked[:,c_west:-c_west]
# 1 m
bstrcwmax_total_masked_mid_masked_trimmed = bstrcwmax_total_masked_mid_masked[:,c_west:-c_west]
time_above_mud01_masked_mid_masked_trimmed = time_above_mud01_masked_mid_masked[:,c_west:-c_west]
# Depth-averaged 
bstrcwmax_total_masked_outer_masked_trimmed = bstrcwmax_total_masked_outer_masked[:,c_west:-c_west]
time_above_mud01_masked_outer_masked_trimmed = time_above_mud01_masked_outer_masked[:,c_west:-c_west]
# 0 - 10 m
bstrcwmax_total_masked_10m_masked_trimmed = bstrcwmax_total_masked_10m_masked[:,c_west:-c_west]
time_above_mud01_masked_10m_masked_trimmed = time_above_mud01_masked_10m_masked[:,c_west:-c_west]
# 10 - 20 m
bstrcwmax_total_masked_10_20m_masked_trimmed = bstrcwmax_total_masked_10_20m_masked[:,c_west:-c_west]
time_above_mud01_masked_10_20m_masked_trimmed = time_above_mud01_masked_10_20m_masked[:,c_west:-c_west]
# 10 - 60 m
bstrcwmax_total_masked_10_60m_masked_trimmed = bstrcwmax_total_masked_10_60m_masked[:,c_west:-c_west]
time_above_mud01_masked_10_60m_masked_trimmed = time_above_mud01_masked_10_60m_masked[:,c_west:-c_west]

# Print some statistics 
# Mean 
# Bed Stress
bstrcwmax_total_masked_inner_masked_trimmed_avg = np.nanmean(bstrcwmax_total_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged total bed shear stress (N/m2): ', bstrcwmax_total_masked_inner_masked_trimmed_avg)
bstrcwmax_total_masked_mid_masked_trimmed_avg = np.nanmean(bstrcwmax_total_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged total bed shear stress (N/m2): ', bstrcwmax_total_masked_mid_masked_trimmed_avg)
bstrcwmax_total_masked_outer_masked_trimmed_avg = np.nanmean(bstrcwmax_total_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged total bed shear stress (N/m2): ', bstrcwmax_total_masked_outer_masked_trimmed_avg)
bstrcwmax_total_masked_10m_masked_trimmed_avg = np.nanmean(bstrcwmax_total_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean total bed shear stress (N/m2): ', bstrcwmax_total_masked_10m_masked_trimmed_avg)
bstrcwmax_total_masked_10_20m_masked_trimmed_avg = np.nanmean(bstrcwmax_total_masked_10_20m_masked_trimmed, axis=(0,1))
print('10 - 20 m mean total bed shear stress (N/m2): ', bstrcwmax_total_masked_10_20m_masked_trimmed_avg)
bstrcwmax_total_masked_10_60m_masked_trimmed_avg = np.nanmean(bstrcwmax_total_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean total bed shear stress (N/m2): ', bstrcwmax_total_masked_10_60m_masked_trimmed_avg)
# Frequency of resuspension
time_above_mud01_masked_inner_masked_trimmed_avg = np.nanmean(time_above_mud01_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged frequency of resuspension (%): ', time_above_mud01_masked_inner_masked_trimmed_avg)
time_above_mud01_masked_mid_masked_trimmed_avg = np.nanmean(time_above_mud01_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged frequency of resuspension (%): ', time_above_mud01_masked_mid_masked_trimmed_avg)
time_above_mud01_masked_outer_masked_trimmed_avg = np.nanmean(time_above_mud01_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged frequency of resuspension (%): ', time_above_mud01_masked_outer_masked_trimmed_avg)
time_above_mud01_masked_10m_masked_trimmed_avg = np.nanmean(time_above_mud01_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m mean frequency of resuspension (%): ', time_above_mud01_masked_10m_masked_trimmed_avg)
time_above_mud01_masked_10_60m_masked_trimmed_avg = np.nanmean(time_above_mud01_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m mean frequency of resuspension (%): ', time_above_mud01_masked_10_60m_masked_trimmed_avg)

# Standard deviation 
# Bed Stress
bstrcwmax_total_masked_inner_masked_trimmed_std = np.nanstd(bstrcwmax_total_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) std total bed shear stress (N/m2): ', bstrcwmax_total_masked_inner_masked_trimmed_std)
bstrcwmax_total_masked_mid_masked_trimmed_std = np.nanstd(bstrcwmax_total_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) std total bed shear stress (N/m2): ', bstrcwmax_total_masked_mid_masked_trimmed_std)
bstrcwmax_total_masked_outer_masked_trimmed_std = np.nanstd(bstrcwmax_total_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) std total bed shear stress (N/m2): ', bstrcwmax_total_masked_outer_masked_trimmed_std)
bstrcwmax_total_masked_10m_masked_trimmed_std = np.nanstd(bstrcwmax_total_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m std total bed shear stress (N/m2): ', bstrcwmax_total_masked_10m_masked_trimmed_std)
bstrcwmax_total_masked_10_20m_masked_trimmed_std = np.nanstd(bstrcwmax_total_masked_10_20m_masked_trimmed, axis=(0,1))
print('10 - 20 m std total bed shear stress (N/m2): ', bstrcwmax_total_masked_10_20m_masked_trimmed_std)
bstrcwmax_total_masked_10_60m_masked_trimmed_std = np.nanstd(bstrcwmax_total_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m std (m/s) total bed shear stress (N/m2): ', bstrcwmax_total_masked_10_60m_masked_trimmed_std)
# Frequency of resuspension
time_above_mud01_masked_inner_masked_trimmed_std = np.nanstd(time_above_mud01_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) std frequency of resuspension (%): ', time_above_mud01_masked_inner_masked_trimmed_std)
time_above_mud01_masked_mid_masked_trimmed_std = np.nanstd(time_above_mud01_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) std frequency of resuspension (%): ', time_above_mud01_masked_mid_masked_trimmed_std)
time_above_mud01_masked_outer_masked_trimmed_std = np.nanstd(time_above_mud01_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) std frequency of resuspension (%): ', time_above_mud01_masked_outer_masked_trimmed_std)
time_above_mud01_masked_10m_masked_trimmed_std = np.nanstd(time_above_mud01_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m std frequency of resuspension (%): ', time_above_mud01_masked_10m_masked_trimmed_std)
time_above_mud01_masked_10_60m_masked_trimmed_std = np.nanstd(time_above_mud01_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m std (m/s) frequency of resuspension (%): ', time_above_mud01_masked_10_60m_masked_trimmed_std)

# Min
# Bed Stress
bstrcwmax_total_masked_inner_masked_trimmed_min = np.nanmin(bstrcwmax_total_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) min total bed shear stress (N/m2): ', bstrcwmax_total_masked_inner_masked_trimmed_min)
bstrcwmax_total_masked_mid_masked_trimmed_min = np.nanmin(bstrcwmax_total_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) min total bed shear stress (N/m2): ', bstrcwmax_total_masked_mid_masked_trimmed_min)
bstrcwmax_total_masked_outer_masked_trimmed_min = np.nanmin(bstrcwmax_total_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) min total bed shear stress (N/m2): ', bstrcwmax_total_masked_outer_masked_trimmed_min)
bstrcwmax_total_masked_10m_masked_trimmed_min = np.nanmin(bstrcwmax_total_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m min total bed shear stress (N/m2): ', bstrcwmax_total_masked_10m_masked_trimmed_min)
bstrcwmax_total_masked_10_60m_masked_trimmed_min = np.nanmin(bstrcwmax_total_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m min (m/s) total bed shear stress (N/m2): ', bstrcwmax_total_masked_10_60m_masked_trimmed_min)
# Frequency of resuspension
time_above_mud01_masked_inner_masked_trimmed_min = np.nanmin(time_above_mud01_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) min frequency of resuspension (%): ', time_above_mud01_masked_inner_masked_trimmed_min)
time_above_mud01_masked_mid_masked_trimmed_min = np.nanmin(time_above_mud01_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) min frequency of resuspension (%): ', time_above_mud01_masked_mid_masked_trimmed_min)
time_above_mud01_masked_outer_masked_trimmed_min = np.nanmin(time_above_mud01_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) min frequency of resuspension (%): ', time_above_mud01_masked_outer_masked_trimmed_min)
time_above_mud01_masked_10m_masked_trimmed_min = np.nanmin(time_above_mud01_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m min frequency of resuspension (%): ', time_above_mud01_masked_10m_masked_trimmed_min)
time_above_mud01_masked_10_60m_masked_trimmed_min = np.nanmin(time_above_mud01_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m min (m/s) frequency of resuspension (%): ', time_above_mud01_masked_10_60m_masked_trimmed_min)

# Max
# Bed Stress
bstrcwmax_total_masked_inner_masked_trimmed_max = np.nanmax(bstrcwmax_total_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) max total bed shear stress (N/m2): ', bstrcwmax_total_masked_inner_masked_trimmed_max)
bstrcwmax_total_masked_mid_masked_trimmed_max = np.nanmax(bstrcwmax_total_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) max total bed shear stress (N/m2): ', bstrcwmax_total_masked_mid_masked_trimmed_max)
bstrcwmax_total_masked_outer_masked_trimmed_max = np.nanmax(bstrcwmax_total_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) max total bed shear stress (N/m2): ', bstrcwmax_total_masked_outer_masked_trimmed_max)
bstrcwmax_total_masked_10m_masked_trimmed_max = np.nanmax(bstrcwmax_total_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m max total bed shear stress (N/m2): ', bstrcwmax_total_masked_10m_masked_trimmed_max)
bstrcwmax_total_masked_10_60m_masked_trimmed_max = np.nanmax(bstrcwmax_total_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m max (m/s) total bed shear stress (N/m2): ', bstrcwmax_total_masked_10_60m_masked_trimmed_max)
# Frequency of resuspension
time_above_mud01_masked_inner_masked_trimmed_max = np.nanmax(time_above_mud01_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) max frequency of resuspension (%): ', time_above_mud01_masked_inner_masked_trimmed_max)
time_above_mud01_masked_mid_masked_trimmed_max = np.nanmax(time_above_mud01_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) max frequency of resuspension (%): ', time_above_mud01_masked_mid_masked_trimmed_max)
time_above_mud01_masked_outer_masked_trimmed_max = np.nanmax(time_above_mud01_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) max frequency of resuspension (%): ', time_above_mud01_masked_outer_masked_trimmed_max)
time_above_mud01_masked_10m_masked_trimmed_max = np.nanmax(time_above_mud01_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m max frequency of resuspension (%): ', time_above_mud01_masked_10m_masked_trimmed_max)
time_above_mud01_masked_10_60m_masked_trimmed_max = np.nanmax(time_above_mud01_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m max (m/s) frequency of resuspension (%): ', time_above_mud01_masked_10_60m_masked_trimmed_max)



# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Set up the data
roms_bed_stress_cur_wave_freq_resusp = xr.Dataset(
    data_vars=dict(
        roms_bed_shear_stress_time_avg=(['y','x'], bstrcwmax_total_wland_masked_trimmed.values),
        roms_current_bed_stress_over_total_bed_stress_time_avg_percent=(['y','x'], (bstrc_tot_mag_wland_masked_trimmed.values/bstrcwmax_total_wland_masked_trimmed.values)*100),
        roms_wave_bed_stress_over_total_bed_stress_time_avg_percent=(['y','x'], (bstrw_tot_mag_wland_masked_trimmed.values/bstrcwmax_total_wland_masked_trimmed.values)*100),
        roms_freq_resuspension_fraction=(['y','x'], time_above_mud01_wland_masked_trimmed.values),
        ),
    coords=dict(
        x_full=('x', x_rho_flat_trimmed),
        y_full=('y', y_rho_flat), 
        ),
    attrs=dict(description='Time-averaged ROMS output including total bed shear stress (N/m2), percentage of current-induced bed stress over total bed stress, percentage of wave-induced bed stress over total bed stress, and frequency of resuspension for critical stress of 0.18 Pa (fraction of time)'))
# Add more metadata?
roms_bed_stress_cur_wave_freq_resusp.roms_bed_shear_stress_time_avg.name='time-averaged total bed shear stress (N/m2)'
roms_bed_stress_cur_wave_freq_resusp.roms_current_bed_stress_over_total_bed_stress_time_avg_percent.name='time-averaged ratio of current-induced bed stress over total bed stress (percentage)'
roms_bed_stress_cur_wave_freq_resusp.roms_wave_bed_stress_over_total_bed_stress_time_avg_percent.name='time-averaged ratio of wave-induced bed stress over total bed stress (percentage)'
roms_bed_stress_cur_wave_freq_resusp.roms_freq_resuspension_fraction.name='frequency of resuspension for critical stress of 0.18 Pa (fraction of time)'

# Save to a netcdf
#roms_bed_stress_cur_wave_freq_resusp.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Paper1_Take2/Data/fig7_roms_time_avg_bed_stress_wave_cur_freq_resuspension.nc')



