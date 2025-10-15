####################### Plot Frequency of Resuspension and Differences Plots #########################
# The purpose of this script is to plot the frequency of resuspension for the standard run, then 
# calcualte the same values for the no-waves and double waves run, then plot the difference 
# between these sensitivity tests and the standard run. This may change but that is the current
# vision.
#
# Notes:
# - 
######################################################################################################


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

# Load in the rho masks
mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
#mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
#mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')


# ---------------------- Define a bunch of functions -------------------------
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


# ----------------------- End Function Definitions -------------------------

# Loop through model output and call the function
# First, get all the file names 
# ROMS 2020 output
# dbsed0003
file_names = glob('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0003/ocean_his_beaufort_shelf_2020_dbsed0003_*.nc')

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


# Find the frequency of resuspension
# Load in the rho masks 
mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
#mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
#mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')

# Set the critical shear stresses 
tau_crit_mud01 = 0.18 # N/m2 , 0.05
tau_crit_mud02 = 0.18 # N/m2 , 0.05
tau_crit_sand01 = 0.18 # N/m2 , 0.1
tau_crit_sand02 = 0.34 # N/m2 , 0.24
tau_crit_sand03 = 24 # N/m2 # 3783

# Make some arrays to hold output
time_lengths_std = np.empty((num_files))
time_above_mud01_cnt_std = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_mud02_cnt_std = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand01_cnt_std = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand02_cnt_std = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand03_cnt_std = np.empty((num_files, eta_rho_len, xi_rho_len))

# Loop through the model output
for j in range(num_files):
        
    # Call the function to open the output
    time_len_tmp, time_above_mud01_tmp, time_above_mud02_tmp, time_above_sand01_tmp, time_above_sand02_tmp, time_above_sand03_tmp = cnt_time_steps_above_crit(file_names2[j], tau_crit_mud01, tau_crit_mud02, tau_crit_sand01, tau_crit_sand02, tau_crit_sand03)
    
    # Save these to the arrays
    time_lengths_std[j] = time_len_tmp
    time_above_mud01_cnt_std[j,:,:] = time_above_mud01_tmp
    time_above_mud02_cnt_std[j,:,:] = time_above_mud02_tmp
    time_above_sand01_cnt_std[j,:,:] = time_above_sand01_tmp
    time_above_sand02_cnt_std[j,:,:] = time_above_sand02_tmp
    time_above_sand03_cnt_std[j,:,:] = time_above_sand03_tmp
    
    
# Now that the counts have been found for all outputs,
# combine them into one percent for each location in grid
# Find the total length of time steps
total_time_len_std = np.sum(time_lengths_std, axis=0)
time_above_mud01_std = (np.sum(time_above_mud01_cnt_std, axis=0)/total_time_len_std)
time_above_mud02_std = (np.sum(time_above_mud02_cnt_std, axis=0)/total_time_len_std)
time_above_sand01_std = (np.sum(time_above_sand01_cnt_std, axis=0)/total_time_len_std)
time_above_sand02_std = (np.sum(time_above_sand02_cnt_std, axis=0)/total_time_len_std)
time_above_sand03_std = (np.sum(time_above_sand03_cnt_std, axis=0)/total_time_len_std)


# Check shapes to see if this worked
print('time_above_mud01_cnt shape (standard run): ', np.shape(time_above_mud01_cnt_std), flush=True)
print('time_above_mud01 shape (standard run): ', np.shape(time_above_mud01_std), flush=True)
print('total_time_len (standard run): ', total_time_len_std, flush=True)

# Load in the last model output to use for plotting things
#model_output = xr.open_dataset(file_names2[-1])


# Print these value to the output file
print('Mud01 min percent (standard run): ', time_above_mud01_std.min(), '; Mud01 max percent (standard run): ', time_above_mud01_std.max(), flush=True)
print('Mud02 min percent (standard run): ', time_above_mud02_std.min(), '; Mud02 max percent (standard run): ', time_above_mud02_std.max(), flush=True)
print('Sand01 min percent (standard run): ', time_above_sand01_std.min(), '; Sand01 max percent (standard run): ', time_above_sand01_std.max(), flush=True)
print('Sand02 min percent (standard run): ', time_above_sand02_std.min(), '; Sand02 max percent (standard run): ', time_above_sand02_std.max(), flush=True)
print('Sand03 min percent (standard run): ', time_above_sand03_std.min(), '; Sand03 max percent (standard run): ', time_above_sand03_std.max(), flush=True)


# ----------------------- No Waves Run Prep -------------------------------
# Do this all again for the no waves run 

# Loop through model output and call the function
# First, get all the file names 
# ROMS 2020 output
# dbsed0005
file_names_nowaves = glob('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0005_nowaves/ocean_his_beaufort_shelf_2020_dbsed0005_11234017_*.nc')

# Sort them to be in order
file_names3 = sorted(file_names_nowaves)

# Check to see if this worked
print(file_names3[0], flush=True)
print(file_names3[-1], flush=True)

# Pull out the number of files
num_files_nowaves = len(file_names3)

# Pull out the length of time of the full run, the time steps, 
# and the length of time of each output file
full_time_len_nowaves, time_steps_nowaves, time_lengths_nowaves = get_model_time(file_names3, num_files_nowaves)


# Make some arrays to hold output
time_lengths_nowaves = np.empty((num_files_nowaves))
time_above_mud01_cnt_nowaves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_mud02_cnt_nowaves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand01_cnt_nowaves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand02_cnt_nowaves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand03_cnt_nowaves = np.empty((num_files, eta_rho_len, xi_rho_len))

# Loop through the model output
for jj in range(num_files_nowaves):
        
    # Call the function to open the output
    time_len_tmp, time_above_mud01_tmp, time_above_mud02_tmp, time_above_sand01_tmp, time_above_sand02_tmp, time_above_sand03_tmp = cnt_time_steps_above_crit(file_names3[jj], tau_crit_mud01, tau_crit_mud02, tau_crit_sand01, tau_crit_sand02, tau_crit_sand03)
    
    # Save these to the arrays
    time_lengths_nowaves[jj] = time_len_tmp
    time_above_mud01_cnt_nowaves[jj,:,:] = time_above_mud01_tmp
    time_above_mud02_cnt_nowaves[jj,:,:] = time_above_mud02_tmp
    time_above_sand01_cnt_nowaves[jj,:,:] = time_above_sand01_tmp
    time_above_sand02_cnt_nowaves[jj,:,:] = time_above_sand02_tmp
    time_above_sand03_cnt_nowaves[jj,:,:] = time_above_sand03_tmp
    
    
# Now that the counts have been found for all outputs,
# combine them into one percent for each location in grid
# Find the total length of time steps
total_time_len_nowaves = np.sum(time_lengths_nowaves, axis=0)
time_above_mud01_nowaves = (np.sum(time_above_mud01_cnt_nowaves, axis=0)/total_time_len_nowaves)
time_above_mud02_nowaves = (np.sum(time_above_mud02_cnt_nowaves, axis=0)/total_time_len_nowaves)
time_above_sand01_nowaves = (np.sum(time_above_sand01_cnt_nowaves, axis=0)/total_time_len_nowaves)
time_above_sand02_nowaves = (np.sum(time_above_sand02_cnt_nowaves, axis=0)/total_time_len_nowaves)
time_above_sand03_nowaves = (np.sum(time_above_sand03_cnt_nowaves, axis=0)/total_time_len_nowaves)


# Check shapes to see if this worked
print('time_above_mud01_cnt shape (no waves run): ', np.shape(time_above_mud01_cnt_nowaves), flush=True)
print('time_above_mud01 shape (no waves run): ', np.shape(time_above_mud01_nowaves), flush=True)
print('total_time_len (no waves run): ', total_time_len_nowaves, flush=True)

# Load in the last model output to use for plotting things
#model_output = xr.open_dataset(file_names3[-1])


# Print these value to the output file
print('Mud01 min percent (no waves run): ', time_above_mud01_nowaves.min(), '; Mud01 max percent (no waves run): ', time_above_mud01_nowaves.max(), flush=True)
print('Mud02 min percent (no waves run): ', time_above_mud02_nowaves.min(), '; Mud02 max percent (no waves run): ', time_above_mud02_nowaves.max(), flush=True)
print('Sand01 min percent (no waves run): ', time_above_sand01_nowaves.min(), '; Sand01 max percent (no waves run): ', time_above_sand01_nowaves.max(), flush=True)
print('Sand02 min percent (no waves run): ', time_above_sand02_nowaves.min(), '; Sand02 max percent (no waves run): ', time_above_sand02_nowaves.max(), flush=True)
print('Sand03 min percent (no waves run): ', time_above_sand03_nowaves.min(), '; Sand03 max percent (no waves run): ', time_above_sand03_nowaves.max(), flush=True)


# ------------------------- Double Waves Prep ------------------------------
# Do this again for the double waves

# Loop through model output and call the function
# First, get all the file names 
# ROMS 2020 output
# dbsed0006
file_names_double_waves = glob('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0006_double_waves/ocean_his_beaufort_shelf_2020_dbsed0006_*.nc')

# Sort them to be in order
file_names4 = sorted(file_names_double_waves)

# Check to see if this worked
print(file_names4[0], flush=True)
print(file_names4[-1], flush=True)

# Pull out the number of files
num_files_double_waves = len(file_names4)

# Pull out the length of time of the full run, the time steps, 
# and the length of time of each output file
full_time_len_double_waves, time_steps_double_waves, time_lengths_double_waves = get_model_time(file_names4, num_files_double_waves)


# Make some arrays to hold output
time_lengths_double_waves = np.empty((num_files_double_waves))
time_above_mud01_cnt_double_waves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_mud02_cnt_double_waves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand01_cnt_double_waves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand02_cnt_double_waves = np.empty((num_files, eta_rho_len, xi_rho_len))
time_above_sand03_cnt_double_waves = np.empty((num_files, eta_rho_len, xi_rho_len))

# Loop through the model output
for jjj in range(num_files_double_waves):
        
    # Call the function to open the output
    time_len_tmp, time_above_mud01_tmp, time_above_mud02_tmp, time_above_sand01_tmp, time_above_sand02_tmp, time_above_sand03_tmp = cnt_time_steps_above_crit(file_names4[jjj], tau_crit_mud01, tau_crit_mud02, tau_crit_sand01, tau_crit_sand02, tau_crit_sand03)
    
    # Save these to the arrays
    time_lengths_double_waves[jjj] = time_len_tmp
    time_above_mud01_cnt_double_waves[jjj,:,:] = time_above_mud01_tmp
    time_above_mud02_cnt_double_waves[jjj,:,:] = time_above_mud02_tmp
    time_above_sand01_cnt_double_waves[jjj,:,:] = time_above_sand01_tmp
    time_above_sand02_cnt_double_waves[jjj,:,:] = time_above_sand02_tmp
    time_above_sand03_cnt_double_waves[jjj,:,:] = time_above_sand03_tmp
    
    
# Now that the counts have been found for all outputs,
# combine them into one percent for each location in grid
# Find the total length of time steps
total_time_len_double_waves = np.sum(time_lengths_double_waves, axis=0)
time_above_mud01_double_waves = (np.sum(time_above_mud01_cnt_double_waves, axis=0)/total_time_len_double_waves)
time_above_mud02_double_waves = (np.sum(time_above_mud02_cnt_double_waves, axis=0)/total_time_len_double_waves)
time_above_sand01_double_waves = (np.sum(time_above_sand01_cnt_double_waves, axis=0)/total_time_len_double_waves)
time_above_sand02_double_waves = (np.sum(time_above_sand02_cnt_double_waves, axis=0)/total_time_len_double_waves)
time_above_sand03_double_waves = (np.sum(time_above_sand03_cnt_double_waves, axis=0)/total_time_len_double_waves)


# Check shapes to see if this worked
print('time_above_mud01_cnt shape (double waves run): ', np.shape(time_above_mud01_cnt_double_waves), flush=True)
print('time_above_mud01 shape (double waves run): ', np.shape(time_above_mud01_double_waves), flush=True)
print('total_time_len (double waves run): ', total_time_len_double_waves, flush=True)

# Load in the last model output to use for plotting things
model_output = xr.open_dataset(file_names4[-1])

# Print these value to the output file
print('Mud01 min percent (double waves run): ', time_above_mud01_double_waves.min(), '; Mud01 max percent (double waves run): ', time_above_mud01_double_waves.max(), flush=True)
print('Mud02 min percent (double waves run): ', time_above_mud02_double_waves.min(), '; Mud02 max percent (double waves run): ', time_above_mud02_double_waves.max(), flush=True)
print('Sand01 min percent (double waves run): ', time_above_sand01_double_waves.min(), '; Sand01 max percent (double waves run): ', time_above_sand01_double_waves.max(), flush=True)
print('Sand02 min percent (double waves run): ', time_above_sand02_double_waves.min(), '; Sand02 max percent (double waves run): ', time_above_sand02_double_waves.max(), flush=True)
print('Sand03 min percent (double waves run): ', time_above_sand03_double_waves.min(), '; Sand03 max percent (double waves run): ', time_above_sand03_double_waves.max(), flush=True)


# Prep the data for plotting 
# Mask the data 
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

# Make land gray
# Make mask
temp_mask = grid.mask_rho.values
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)

# Mask the data
# Standard run 
time_above_mud01_std_masked = time_above_mud01_std*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_mud02_std_masked = time_above_mud02_std*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand01_std_masked = time_above_sand01_std*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand02_std_masked = time_above_sand02_std*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand03_std_masked = time_above_sand03_std*temp_mask*mask_rho_nan.nudge_mask_rho_nan

# No waves run 
time_above_mud01_nowaves_masked = time_above_mud01_nowaves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_mud02_nowaves_masked = time_above_mud02_nowaves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand01_nowaves_masked = time_above_sand01_nowaves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand02_nowaves_masked = time_above_sand02_nowaves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand03_nowaves_masked = time_above_sand03_nowaves*temp_mask*mask_rho_nan.nudge_mask_rho_nan

# Double waves run 
time_above_mud01_double_waves_masked = time_above_mud01_double_waves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_mud02_double_waves_masked = time_above_mud02_double_waves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand01_double_waves_masked = time_above_sand01_double_waves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand02_double_waves_masked = time_above_sand02_double_waves*temp_mask*mask_rho_nan.nudge_mask_rho_nan
time_above_sand03_double_waves_masked = time_above_sand03_double_waves*temp_mask*mask_rho_nan.nudge_mask_rho_nan

# Trim the data
# Standard run
time_above_mud01_std_masked_trimmed = time_above_mud01_std_masked[:,c_west:-c_west]
time_above_mud02_std_masked_trimmed = time_above_mud02_std_masked[:,c_west:-c_west]
time_above_sand01_std_masked_trimmed = time_above_sand01_std_masked[:,c_west:-c_west]
time_above_sand02_std_masked_trimmed = time_above_sand02_std_masked[:,c_west:-c_west]
time_above_sand03_std_masked_trimmed = time_above_sand03_std_masked[:,c_west:-c_west]

# No waves run 
time_above_mud01_nowaves_masked_trimmed = time_above_mud01_nowaves_masked[:,c_west:-c_west]
time_above_mud02_nowaves_masked_trimmed = time_above_mud02_nowaves_masked[:,c_west:-c_west]
time_above_sand01_nowaves_masked_trimmed = time_above_sand01_nowaves_masked[:,c_west:-c_west]
time_above_sand02_nowaves_masked_trimmed = time_above_sand02_nowaves_masked[:,c_west:-c_west]
time_above_sand03_nowaves_masked_trimmed = time_above_sand03_nowaves_masked[:,c_west:-c_west]

# Double waves run 
time_above_mud01_double_waves_masked_trimmed = time_above_mud01_double_waves_masked[:,c_west:-c_west]
time_above_mud02_double_waves_masked_trimmed = time_above_mud02_double_waves_masked[:,c_west:-c_west]
time_above_sand01_double_waves_masked_trimmed = time_above_sand01_double_waves_masked[:,c_west:-c_west]
time_above_sand02_double_waves_masked_trimmed = time_above_sand02_double_waves_masked[:,c_west:-c_west]
time_above_sand03_double_waves_masked_trimmed = time_above_sand03_double_waves_masked[:,c_west:-c_west]

# Make reverse mask so land can be gray 
reverse_colors = colors.ListedColormap(['darkgray','#FF000000'])
bounds=[0,1,2]
norm = colors.BoundaryNorm(bounds, reverse_colors.N)

# Make new mask so land is gray - NEW
mask_rho2 = grid.mask_rho.values
mask_rho2[mask_rho2==1.0] = 2.0
mask_rho2_trimmed = mask_rho2[:,c_west:-c_west]

# ------------------------------- Plot the Data -----------------------------------

# ---------------------------------------------------------------------------
# ------------- Plot 1: Frequency of Resuspension for Each Run --------------
# ---------------------------------------------------------------------------
# Plot the frequency of resuspension for each run

# Make the figure
fig1, ax1 = plt.subplots(3, figsize=(20, 22))  # (22,21) (18,16) (16,18) (16,21) (18,25)

# Set a colormap and levels to use for all plots
cmap1 = cmocean.cm.matter
#cmap8b.set_under('darkgray')
lev_freq = np.arange(0, 0.75, 0.01)
lev_bathy = np.arange(10,70,10)

# Standard Run
# Plot land as gray 
cs0a = ax1[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs1 = ax1[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs2 = ax1[0].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 time_above_mud01_std_masked_trimmed, lev_freq, cmap=cmap1, extend='max')
# Label the plot
plt.setp(ax1[0].get_xticklabels(), visible=False)
ax1[0].set_ylabel('Y (km)', fontsize=fontsize-2)
ax1[0].set_title('Standard Run', fontsize=fontsize)

# No Waves Run
# Plot land as gray 
cs0b = ax1[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs3 = ax1[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs4 = ax1[1].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 time_above_mud01_nowaves_masked_trimmed, lev_freq, cmap=cmap1, extend='max')
# Label the plot
plt.setp(ax1[1].get_xticklabels(), visible=False)
ax1[1].set_ylabel('Y (km)', fontsize=fontsize-2)
ax1[1].set_title('No Waves Run', fontsize=fontsize)

# Double Waves Run
# Plot land as gray 
cs0c = ax1[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs5 = ax1[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs6 = ax1[2].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 time_above_mud01_double_waves_masked_trimmed, lev_freq, cmap=cmap1, extend='max')
# Label the plot
plt.setp(ax1[2].get_xticklabels(), visible=False)
ax1[2].set_ylabel('Y (km)', fontsize=fontsize-2)
ax1[2].set_title('Double Waves Run', fontsize=fontsize)


# Make and label the colorbar
cbar1_ax = fig1.add_axes([0.92, 0.11, 0.02, 0.77]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
cbar1 = plt.colorbar(cs6, orientation='vertical', cax=cbar1_ax).set_label(label='Fraction of Time', size=fontsize-2, labelpad=15)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.085, left=0.08)

# Add labels
# Bottom right 
plt.text(0.876, 0.648, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.876, 0.393, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.876, 0.122, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Label subplots
#plt.text(0.470, 0.700, 'Current Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.470, 0.487, 'Wave Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.393, 0.268, 'Tau Critical: ' +str(tau_crit_mud01) + ' (Pa)', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the plot
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/freq_of_resusp_standard_no_double_0005.png')



# ---------------------------------------------------------------------------
# ------------- Plot 2: Frequency of Resuspension Difference Plots --------------
# ---------------------------------------------------------------------------
# Plot the frequency of resuspension for each sensitivity test and the difference plot

# Make the figure
fig2, ax2 = plt.subplots(4, figsize=(19, 23))  # (22,21) (18,16) (16,18) (16,21) (18,25) (width, height)

# Set a colormap and levels to use for all plots
cmap1 = cmocean.cm.matter
cmap_diff = cmocean.cm.balance
#cmap8b.set_under('darkgray')
lev_freq = np.arange(0, 0.75, 0.01)
lev_bathy = np.arange(10,70,10)
lev_diff = np.arange(-0.5, 0.5, 0.05)

# --- No waves ---
# - Frequency of resuspension -
# Plot land as gray 
cs0d = ax2[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs7 = ax2[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs8 = ax2[0].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 time_above_mud01_nowaves_masked_trimmed, lev_freq, cmap=cmap1, extend='max')
# Label the plot
plt.setp(ax2[0].get_xticklabels(), visible=False)
ax2[0].set_ylabel('Y (km)', fontsize=fontsize-2)
ax2[0].set_title('No Waves Run', fontsize=fontsize-1)
# Make and label the colorbar
#cbar2_ax = fig2.add_axes([0.82, 0.74, 0.02, 0.2]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
#cbar2 = plt.colorbar(cs8, orientation='vertical', cax=cbar2_ax).set_label(label='Fraction of Time', size=fontsize-2, labelpad=15)
# If horizontal and in axes
cbar2 = plt.colorbar(cs8, cax=ax2[0].inset_axes((0.48, 0.92, 0.5, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 0.15, 0.30, 0.45, 0.60], ax=ax2[0], 
                     orientation='horizontal').set_label(label='Fraction of Time', size=fontsize-2)

# - Difference from standard -
# Plot land as gray 
cs0e = ax2[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs9 = ax2[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs10 = ax2[1].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                (time_above_mud01_nowaves_masked_trimmed-time_above_mud01_std_masked_trimmed), lev_diff, cmap=cmap_diff, extend='both')
# Label the plot
plt.setp(ax2[1].get_xticklabels(), visible=False)
ax2[1].set_ylabel('Y (km)', fontsize=fontsize-2)
ax2[1].set_title('No Waves Run Minus Standard', fontsize=fontsize-1)
# Set the colorbar
# If horizontal and in axes
cbar3 = plt.colorbar(cs10, cax=ax2[1].inset_axes((0.48, 0.92, 0.5, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[-0.50, -0.25, 0.0, 0.25, 0.50], ax=ax2[1], 
                     orientation='horizontal').set_label(label='Difference in Fraction of Time', size=fontsize-2)


# --- Double Waves Run ---
# - Frequency of Resuspension -
# Plot land as gray 
cs0f = ax2[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs11 = ax2[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs12 = ax2[2].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 (time_above_mud01_double_waves_masked_trimmed), lev_freq, cmap=cmap1, extend='max')
# Label the plot
plt.setp(ax2[2].get_xticklabels(), visible=False)
ax2[2].set_ylabel('Y (km)', fontsize=fontsize-2)
ax2[2].set_title('Double Waves Run', fontsize=fontsize-1)
# Make and label the colorbar
#cbar3_ax = fig2.add_axes([0.82, 0.2, 0.02, 0.18]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
#cbar3 = plt.colorbar(cs12, orientation='vertical', cax=cbar3_ax).set_label(label='Difference (Fraction of Time)', size=fontsize-2, labelpad=15)
# If horizontal and in axes
cbar4 = plt.colorbar(cs12, cax=ax2[2].inset_axes((0.48, 0.92, 0.5, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 0.15, 0.30, 0.45, 0.60], ax=ax2[2], 
                     orientation='horizontal').set_label(label='Fraction of Time', size=fontsize-2)

# - Difference -
# Plot land as gray 
cs0f = ax2[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, mask_rho2_trimmed, cmap=reverse_colors, norm=norm)
# Plot bathymetry contours 
cs13 = ax2[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, 
                       h_masked_trimmed, lev_bathy, colors='dimgray')
# Plot the frequency of resuspension
cs14 = ax2[3].contourf(x_rho_flat_trimmed/1000,  y_rho_flat/1000, 
                 (time_above_mud01_double_waves_masked_trimmed-time_above_mud01_std_masked_trimmed), lev_diff, cmap=cmap_diff, extend='both')
# Label the plot
plt.setp(ax2[3].get_xticklabels(), visible=False)
ax2[3].set_ylabel('Y (km)', fontsize=fontsize-2)
ax2[3].set_xlabel('X (km)', size=fontsize-2)
ax2[3].set_title('Double Waves Run Minus Standard', fontsize=fontsize)
# Make and label the colorbar
#cbar3_ax = fig2.add_axes([0.82, 0.2, 0.02, 0.18]) # (0.85, 0.1, 0.02, 0.15) [left, bottom, width, height]
#cbar3 = plt.colorbar(cs12, orientation='vertical', cax=cbar3_ax).set_label(label='Difference (Fraction of Time)', size=fontsize-2, labelpad=15)
# If horizontal and in axes
cbar5 = plt.colorbar(cs14, cax=ax2[3].inset_axes((0.48, 0.92, 0.50, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[-0.50, -0.25, 0.0, 0.25, 0.50], ax=ax2[3], 
                     orientation='horizontal').set_label(label='Difference in Fraction of Time', size=fontsize-2)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.14)

# Add labels
# Add labels
# Bottom right 
plt.text(0.876, 0.718, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.876, 0.518, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.876, 0.317, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.876, 0.115, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Bottom right 
#plt.text(0.843, 0.660, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.853, 0.4402, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.863, 0.129, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Label subplots
#plt.text(0.470, 0.700, 'Current Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.470, 0.487, 'Wave Stress Over Total Stress', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.393, 0.268, 'Tau Critical: ' +str(tau_crit_mud01) + ' (Pa)', fontsize=fontsize+2, fontweight='bold', color='k', transform=plt.gcf().transFigure)

# Save the plot
#plt.savefig('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/freq_of_resusp_diff_0004.png')




# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Set up the data
roms_time_above_tau_crit = xr.Dataset(
    data_vars=dict(
        time_above_tau_crit_std=(['y','x'], time_above_mud01_std_masked_trimmed.values),
        time_above_tau_crit_double_waves=(['y','x'], time_above_mud01_double_waves_masked_trimmed.values),
        time_above_tau_crit_nowaves=(['y','x'], time_above_mud01_nowaves_masked_trimmed.values)
        ),
    coords=dict(
        x_full=('x', x_rho_flat_trimmed),
        y_full=('y', y_rho_flat), 
        ),
    attrs=dict(description='Time above critical shear stress (0.18 Pa) for the standard, no waves, and double waves model runs')) 
# Add more metadata?
roms_time_above_tau_crit.time_above_tau_crit_std.name='fraction of time above critical shear stress (0.18 Pa) in standard model run'
roms_time_above_tau_crit.time_above_tau_crit_double_waves.name='fraction of time above critical shear stress (0.18 Pa) in double waves run'
roms_time_above_tau_crit.time_above_tau_crit_nowaves.name='fraction of time above critical shear stress (0.18 Pa) in no waves run'

# Save to a netcdf
#roms_time_above_tau_crit.to_netcdf('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Scripts/Analysis/fig10_roms_freq_resusp_std_nowaves_double.nc')




