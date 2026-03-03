##################### Time-Averaged Currents #######################
# The purpose of this script is to look at the time-averaged 
# current direction on the shelf in ROMS model output. Tiem-averaged
# currents will be look at for the surface/top sigma layer, the seabed/
# bottom sigma layer, and depth-averaged.
#
# 
# Notes:
# - This script has been updated to use the 2020 ROMS  output
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

# Everything will need to be done once for each output file 
# but setting up the grids for interpolating and finding the weights
# can just be done once so do that first  

# Before we can look at average current directions, we need
# to interpolate ubar and vbar onto rho points

# Set the input and output grids, and sepcify the lat/lon
# Since we are looking at ubar and vbar, we will use lon_u and lat_u as the primary lat/lon for the grid 
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


# Define a function to process the model output
def interpolate_uv_currents(filename, regridder_u2rho, regridder_v2rho):
    """
    This function takes a given model output file and processes
    the current data by first regridding the u and v currents 
    onto rho points.
    
    Inputs:
    - filename: string, path to and name of netcdf output
    - regridder_u2rho: weights from xesmf to regrid from u to rho points
    - regridder_v2rho: weights from xesmf to regrid from v to rho points
    
    Outputs:
    - u_rho_tmp: u current interpolated to rho points
    - v_rho_tmp: v current interpolated to rho points
    - ubar_rho_tmp: depth-averaged u currents interpolated to rho points 
    - vbar_rho_tmp: depth-averaged v currents interpolated to rho points 
    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Use the regridder weights to regrid the currents
    # All depths
    # u
    dr_u2rho_u = model_output['u'].copy()
    u_rho_tmp = regridder_u2rho(dr_u2rho_u)
    
    # v
    dr_v2rho_v = model_output['v'].copy()
    v_rho_tmp = regridder_v2rho(dr_v2rho_v) 
    
    # Depth-averaged 
    # ubar
    dr_u2rho_ubar = model_output['ubar'].copy()
    ubar_rho_tmp = regridder_u2rho(dr_u2rho_ubar) 
    
    # vbar
    dr_v2rho_vbar = model_output['vbar'].copy()
    vbar_rho_tmp = regridder_v2rho(dr_v2rho_vbar) 
    
    # Return the means
    return(u_rho_tmp, v_rho_tmp, ubar_rho_tmp, vbar_rho_tmp)


# Define a function to process the model output
def interpolate_uv_surf_depthavg_currents(filename, regridder_u2rho, regridder_v2rho):
    """
    This function takes a given model output file and processes
    the current data by first regridding the u and v currents 
    onto rho points.
    
    Inputs:
    - filename: string, path to and name of netcdf output
    - regridder_u2rho: weights from xesmf to regrid from u to rho points
    - regridder_v2rho: weights from xesmf to regrid from v to rho points
    
    Outputs:
    - u_rho_tmp: u current interpolated to rho points
    - v_rho_tmp: v current interpolated to rho points
    - ubar_rho_tmp: depth-averaged u currents interpolated to rho points 
    - vbar_rho_tmp: depth-averaged v currents interpolated to rho points 
    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Use the regridder weights to regrid the currents
    # All depths
    # u
    dr_u2rho_u = model_output['u'].copy()
    u_rho_tmp = regridder_u2rho(dr_u2rho_u)
    u_surf_rho_tmp = u_rho_tmp[:,-1,:,:]
    
    # v
    dr_v2rho_v = model_output['v'].copy()
    v_rho_tmp = regridder_v2rho(dr_v2rho_v)
    v_surf_rho_tmp = v_rho_tmp[:,-1,:,:]
    
    # Depth-averaged 
    # ubar
    dr_u2rho_ubar = model_output['ubar'].copy()
    ubar_rho_tmp = regridder_u2rho(dr_u2rho_ubar) 
    
    # vbar
    dr_v2rho_vbar = model_output['vbar'].copy()
    vbar_rho_tmp = regridder_v2rho(dr_v2rho_vbar) 
    
    
    # Return the means
    return(u_surf_rho_tmp, v_surf_rho_tmp, ubar_rho_tmp, vbar_rho_tmp)



# Make a function to pull out the spatial time series of salinity 
def get_salt_time_series(filename):
    """
    This function opens a given model output and pulls out the time series
    for salinity and for depth-averaged salinity.

    Parameters
    ----------
    filename : Path to model output file (string)

    Returns
    -------
    salt_tmp: Salinity time series on rho points
    salt_depthavg_tmp: time series of depth-averaged salinity 

    """
    
    # Load in the model output
    model_output = xr.open_dataset(filename)
    
    # Pull out the bstrcwmax time series
    salt_tmp = model_output.salt[:,:,:,:]
    
    # Take the depth-average and return as time series 
    # Calculate the time-varying thickness of the cells
    cell_thick = abs(model_output.z_w[:,:-1,:,:].values - model_output.z_w[:,1:,:,:].values)
    
    # Multiply each layer by its thickness (PSU*m ) and
    # sum up over all layers to get depth-integrated SSC
    # Multiply by thickness
    salt_times_thick = salt_tmp*cell_thick
    
    # Sum over layers to get depth-integrated SSC (kg/m2)
    salt_depth_int = salt_times_thick.sum(dim='s_rho')
    
    # Divide by bathymetry to get depth-averaged SSC (kg/m3)
    salt_depthavg_tmp = salt_depth_int/model_output.bath[:,:,:].values
    
    # Return these arrays
    return(salt_tmp, salt_depthavg_tmp)




# Make a function to process the output
def interp_roms_output_to1m_rho(filename, regridder_u2rho, regridder_v2rho):
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
    
    # Pull out u, v, and salt
    #u_tmp = ds.u
    #v_tmp = ds.v
    #salt_tmp = ds.salt
    
    # Interpolate these onto the rho points 
    #u_rho_tmp = regridder_u2rho(u_tmp) 
    #v_rho_tmp = regridder_v2rho(v_tmp) 

    # Interopolate onto the given CODA depths
    u_rho_interp_1m = xroms.isoslice(regridder_u2rho(ds.u), depths, xgrid, 
                                         iso_array=height_from_seabed, axis='Z')
    v_rho_interp_1m = xroms.isoslice(regridder_v2rho(ds.v), depths, xgrid, 
                                         iso_array=height_from_seabed, axis='Z')
    #salt_roms_interp_1m = xroms.isoslice(ds.salt, depths, xgrid, 
     #                                    iso_array=height_from_seabed, axis='Z')
     
    
    # Return these currents 
    #return(u_rho_interp_1m, v_rho_interp_1m, salt_roms_interp_1m)
    return(u_rho_interp_1m, v_rho_interp_1m)


# Make a function to process the output 
def interp_roms_currents_to1m_fromsurf_rho(filename, regridder_u2rho, regridder_v2rho):
    """
    This function takes a given ROMS ocean_his file, opens it and 
    interpolates the currents onto the given depths, then returns the interpolated
    data. Right now, it is set up to interpolate u and v currents.
    

    Returns
    -------
    None.

    """
    
    # Load in the ROMS output 
    ds = xr.open_dataset(filename)
    ds['h'] = ds.bath
    ds, xgrid = xroms.roms_dataset(ds)
    
    # Set two depths, 1 and 2 m above seafloor
    #depths = np.asarray([-1.0])
    #depths = np.asarray([-2.0])
    depths = np.asarray([-0.5,-1.5])
    

    # Interopolate onto the given CODA depths
    u_interp_1m_fromsurf = xroms.isoslice(ds.u, depths, xgrid, axis='Z')
    v_interp_1m_fromsurf = xroms.isoslice(ds.v, depths, xgrid, axis='Z')
    #salt_roms_interp_1m = xroms.isoslice(ds.salt, depths, xgrid, 
     #                                    iso_array=height_from_seabed, axis='Z')
    
    # Interpolate onto rho points
    u_rho_interp_1m_fromsurf = regridder_u2rho(u_interp_1m_fromsurf)
    v_rho_interp_1m_fromsurf = regridder_v2rho(v_interp_1m_fromsurf)
    
    print('u_rho_interp_1m_fromsurf shape: ', np.shape(u_rho_interp_1m_fromsurf), flush=True)
    
    # If taking the value over two depths, take the average 
    if i == 6:
        u_rho_interp_1m_fromsurf = np.nanmean(u_rho_interp_1m_fromsurf, axis=0)
        v_rho_interp_1m_fromsurf = np.nanmean(v_rho_interp_1m_fromsurf, axis=0)
    else:
        u_rho_interp_1m_fromsurf = np.nanmean(u_rho_interp_1m_fromsurf, axis=1)
        v_rho_interp_1m_fromsurf = np.nanmean(v_rho_interp_1m_fromsurf, axis=1)
        
    
    # Return these currents 
    #return(u_rho_interp_1m, v_rho_interp_1m, salt_roms_interp_1m)
    return(u_rho_interp_1m_fromsurf, v_rho_interp_1m_fromsurf)


# Make a function to process the output
def interp_salt_to1m(filename):
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
    
    # Set two depths, 1 m from surface and above bottom
    #depths_surf = np.asarray([-1.0])
    #depths_surf = np.asarray([-2.0])
    depths_surf = np.asarray([-0.5,-1.5])
    depths_bot = np.asarray([1.0])

    # Interopolate onto the given CODA depths
    # 1m from surface
    salt_roms_interp_1m_fromsurf = xroms.isoslice(ds.salt, depths_surf, xgrid, axis='Z')
    
    # 1m from bottom
    salt_roms_interp_1m_frombot = xroms.isoslice(ds.salt, depths_bot, xgrid, 
                                         iso_array=height_from_seabed, axis='Z')
    
    # If taking the value over two depths, take the average 
    if i == 6:
        salt_roms_interp_1m_fromsurf = np.nanmean(salt_roms_interp_1m_fromsurf, axis=0)
    else:
        salt_roms_interp_1m_fromsurf = np.nanmean(salt_roms_interp_1m_fromsurf, axis=1)
    
    # Return these currents 
    return(salt_roms_interp_1m_fromsurf, salt_roms_interp_1m_frombot)




# =============================================================================
# 
# # Define a function to process the model output
# def interpolate_and_average_currents(filename, regridder_u2rho, regridder_v2rho):
#     """
#     This function takes a given model output file and processes
#     the current data by first regridding the u and v currents 
#     onto rho points, then pulling out either surface, bottom, or
#     depth-averaged currents, then averages them over time.
#     
#     Inputs:
#     - filename: string, path to and name of netcdf output
#     - regridder_u2rho: weights from xesmf to regrid from u to rho points
#     - regridder_v2rho: weights from xesmf to regrid from v to rho points
#     
#     Outputs:
#     - u_surface_mean_tmp: time-averaged, surface u current interpolated to rho points
#     - u_bottom_mean_tmp: time-averaged, bottom u current interpolated to rho points
#     - v_surface_mean_tmp: time-averaged, surface v current interpolated to rho points
#     - v_bottom_mean_tmp: time-averaged, bottom v current interpolated to rho points
#     - ubar_mean_tmp: time-averaged, depth-averaged u currents interpolated to rho points 
#     - vbar_mean_tmp: time-averaged, depth-averaged v currents interpolated to rho points 
#     """
#     
#     # Load in the model output
#     model_output = xr.open_dataset(filename)
#     
#     # Use the regridder weights to regrid the currents
#     # All depths
#     # u
#     dr_u2rho_u = model_output['u'].copy()
#     dr_out_u2rho_u = regridder_u2rho(dr_u2rho_u)
#     
#     # v
#     dr_v2rho_v = model_output['v'].copy()
#     dr_out_v2rho_v = regridder_v2rho(dr_v2rho_v) 
#     
#     # Depth-averaged 
#     # ubar
#     dr_u2rho_ubar = model_output['ubar'].copy()
#     dr_out_u2rho_ubar = regridder_u2rho(dr_u2rho_ubar) 
#     
#     # vbar
#     dr_v2rho_vbar = model_output['vbar'].copy()
#     dr_out_v2rho_vbar = regridder_v2rho(dr_v2rho_vbar) 
#     
#     # Pull out the surface and bottom currents
#     dr_out_u2rho_u_surface_tmp = dr_out_u2rho_u[:,19,:,:]
#     dr_out_u2rho_u_bottom_tmp = dr_out_u2rho_u[:,0,:,:]
#     dr_out_v2rho_v_surface_tmp = dr_out_v2rho_v[:,19,:,:]
#     dr_out_v2rho_v_bottom_tmp = dr_out_v2rho_v[:,0,:,:]
#     
#     # Take the average over time
#     # u & v
#     u_surface_mean_tmp = dr_out_u2rho_u_surface_tmp.mean(dim='ocean_time')
#     u_bottom_mean_tmp = dr_out_u2rho_u_bottom_tmp.mean(dim='ocean_time')
#     v_surface_mean_tmp = dr_out_v2rho_v_surface_tmp.mean(dim='ocean_time')
#     v_bottom_mean_tmp = dr_out_v2rho_v_bottom_tmp.mean(dim='ocean_time')
#     # ubar & vbar
#     ubar_mean_tmp = dr_out_u2rho_ubar.mean(dim='ocean_time')
#     vbar_mean_tmp = dr_out_v2rho_vbar.mean(dim='ocean_time')
#     
#     # Return the means
#     return(u_surface_mean_tmp, u_bottom_mean_tmp, v_surface_mean_tmp, v_bottom_mean_tmp, ubar_mean_tmp, vbar_mean_tmp)
#     
#     
# # Define a function to process the model output
# def average_east_north_currents(filename):
#     """
#     This function takes a given model output file and processes
#     the current data in the eastward and northward direcitons by 
#     pulling out either surface, bottom, or
#     depth-averaged currents, then averages them over time.
#     
#     Inputs:
#     - filename: string, path to and name of netcdf output
#     
#     Outputs:
#     - u_eastward_surface_mean_tmp: time-averaged, surface u eastward current on rho points
#     - u_eastward_bottom_mean_tmp: time-averaged, bottom u eastward current on rho points
#     - v_northward_surface_mean_tmp: time-averaged, surface v northward current on rho points
#     - v_northward_bottom_mean_tmp: time-averaged, bottom v northward current on rho points
#     - ubar_eastward_mean_tmp: time-averaged, depth-averaged u currents in eastward direction 
#     - vbar_northward_mean_tmp: time-averaged, depth-averaged v currents in northward direction 
#     """
# 
#     # Load in the model output
#     model_output = xr.open_dataset(filename)
# 
#     # Pull out the surface and bottom currents
#     # in the eastward and northward directions 
#     u_eastward_surface_tmp = model_output.u_eastward[:,19,:,:]
#     u_eastward_bottom_tmp = model_output.u_eastward[:,0,:,:]
#     v_northward_surface_tmp = model_output.v_northward[:,19,:,:]
#     v_northward_bottom_tmp = model_output.v_northward[:,0,:,:]
# 
#     # Take the average over time
#     # u_eastward & v_eastward
#     u_eastward_surface_mean_tmp = u_eastward_surface_tmp.mean(dim='ocean_time')
#     u_eastward_bottom_mean_tmp = u_eastward_bottom_tmp.mean(dim='ocean_time')
#     v_northward_surface_mean_tmp = v_northward_surface_tmp.mean(dim='ocean_time')
#     v_northward_bottom_mean_tmp = v_northward_bottom_tmp.mean(dim='ocean_time')
#     # ubar_eastward & vbar_northward
#     ubar_eastward_mean_tmp = model_output.ubar_eastward.mean(dim='ocean_time')
#     vbar_northward_mean_tmp = model_output.vbar_northward.mean(dim='ocean_time')
# 
#     # Return the means
#     return(u_eastward_surface_mean_tmp, u_eastward_bottom_mean_tmp, v_northward_surface_mean_tmp, v_northward_bottom_mean_tmp, ubar_eastward_mean_tmp, vbar_northward_mean_tmp)
# 
# =============================================================================


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
u_surf_rho_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_surf_rho_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
ubar_rho_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
vbar_rho_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))
salt_full = np.empty((full_time_len, s_rho_len, eta_rho_len, xi_rho_len))
salt_depthavg_full = np.empty((full_time_len, eta_rho_len, xi_rho_len))


# Make an xarray dataset to merge dep_net timeseries into 
#dep_nets_xr = xr.Dataset()

# Set a time step to track which time step the loop is on
time_step = 0

# Loop through the model output
for j in range(num_files):
#for j in range(1):

    print('j: ', j)
    
    # Call the function to process the output
    # Get u, v, ubar, vbar on rho points 
    u_surf_rho_tmp, v_surf_rho_tmp, ubar_rho_tmp, vbar_rho_tmp = interpolate_uv_surf_depthavg_currents(file_names2[j], regridder_u2rho, regridder_v2rho)
    
    # Get salinity time series 
    salt_tmp, salt_depthavg_tmp = get_salt_time_series(file_names2[j])
    
    # Save these to the arrays 
    #print('time_step: ', time_step)
    #print('time_step + time_lengths[j]: ', time_step+time_lengths[j])
    start = int(time_step)
    end = int(time_step+time_lengths[j])
    # All sediment classes added together
    u_surf_rho_full[start:end,:,:] = u_surf_rho_tmp
    v_surf_rho_full[start:end,:,:] = v_surf_rho_tmp
    ubar_rho_full[start:end,:,:] = ubar_rho_tmp
    vbar_rho_full[start:end,:,:] = vbar_rho_tmp
    salt_full[start:end,:,:,:] = salt_tmp
    salt_depthavg_full[start:end,:,:] = salt_depthavg_tmp
    
    # Update the base time_step
    time_step = time_step + time_lengths[j]






# Now that w ehave full time series for all variables, split into surface,
# 1 m, and depth-averaged and take time-averages

# OLD using not 1 m from surface but just top grid cell 
# =============================================================================
# # Before averaging over time, find the time series
# # of the current magnitude surface
# cur_mag_surf_rho = np.sqrt((u_surf_rho_full**2)+(v_surf_rho_full**2))
# # Average this over time
# cur_mag_surf_rho_avg = np.mean(cur_mag_surf_rho, axis=0)
# # Surface u, v, salinity
# u_surf_rho_avg = np.mean(u_surf_rho_full[:,:,:], axis=0)
# v_surf_rho_avg = np.mean(v_surf_rho_full[:,:,:], axis=0)
# salt_surf_avg = np.mean(salt_full[:,-1,:,:], axis=0)
# =============================================================================

# Depth-averaged 
# Before averaging over time, find the time series
# of the current magnitude depth-averaged
cur_mag_davg_rho = np.sqrt((ubar_rho_full**2)+(vbar_rho_full**2))
# Average this over time
cur_mag_davg_rho_avg = np.mean(cur_mag_davg_rho, axis=0)
# Get depth-averaged and average over time
ubar_rho_avg = np.mean(ubar_rho_full, axis=0)
vbar_rho_avg = np.mean(vbar_rho_full, axis=0)
salt_depthavg_avg = np.mean(salt_depthavg_full, axis=0)

# 1 meter
# =============================================================================
# # Save full time series to netcdfs (especially since interpoalted to ROMS)
# s_rho_tmp = np.arange(0,20,1)
# # U
# ds_u_rho = xr.DataArray(data=u_rho_full,
#                                       dims=['ocean_time', 's_rho', 'eta_rho', 'xi_rho'],
#                                      coords=dict(ocean_time=(['ocean_time'], time_steps),
#                                                  s_rho=(['s_rho'], s_rho_tmp),
#                                                  eta_rho=(['eta_rho'], grid.eta_rho.values), 
#                                                  xi_rho=(['xi_rho'], grid.xi_rho.values)),
#                                      name='u_rho')
# ds_u_rho.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/Time_Series/u_rho.nc')
# # V
# ds_v_rho = xr.DataArray(data=v_rho_full,
#                                       dims=['ocean_time', 's_rho', 'eta_rho', 'xi_rho'],
#                                      coords=dict(ocean_time=(['ocean_time'], time_steps),
#                                                  s_rho=(['s_rho'], s_rho_tmp),
#                                                  eta_rho=(['eta_rho'], grid.eta_rho.values), 
#                                                  xi_rho=(['xi_rho'], grid.xi_rho.values)),
#                                      name='v_rho')
# ds_v_rho.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/Time_Series/v_rho.nc')
# # Salt
# ds_salt_rho = xr.DataArray(data=salt_full,
#                                       dims=['ocean_time', 's_rho', 'eta_rho', 'xi_rho'],
#                                      coords=dict(ocean_time=(['ocean_time'], time_steps),
#                                                  s_rho=(['s_rho'], s_rho_tmp),
#                                                  eta_rho=(['eta_rho'], grid.eta_rho.values), 
#                                                  xi_rho=(['xi_rho'], grid.xi_rho.values)),
#                                      name='salt_rho')
# ds_salt_rho.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/Time_Series/salt_rho.nc')
# =============================================================================


# Interpolate u, v, and salinity onto 1 m depth
# Load in the netcdf version 
#u_rho_full_nc = xr.open_dataset('/Users//brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/Time_Series/u_rho.nc')
#v_rho_full_nc = xr.open_dataset('/Users//brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/Time_Series/v_rho.nc')
#salt_full_nc = xr.open_dataset('/Users//brun1463/Desktop/Research_Lab/Kaktovik_Alaska/Code/Time_Series/salt_rho.nc')

# Make empty arrays to hold output 
u_1m_fromsurf_rho = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_1m_fromsurf_rho = np.empty((full_time_len, eta_rho_len, xi_rho_len))
u_1m_frombot_rho = np.empty((full_time_len, eta_rho_len, xi_rho_len))
v_1m_frombot_rho = np.empty((full_time_len, eta_rho_len, xi_rho_len))
salt_1m_fromsurf = np.empty((full_time_len, eta_rho_len, xi_rho_len))
salt_1m_frombot = np.empty((full_time_len, eta_rho_len, xi_rho_len))

# Set a time step to track which time step the loop is on
time_step2 = 0

# Loop through the model output
for i in range(num_files):
#for j in range(1):
    print('i: ', i)

    # Call the function to process the output
    # Get u and v 1 m below surface 
    u_roms_1m_fromsurf_tmp, v_roms_1m_fromsurf_tmp = interp_roms_currents_to1m_fromsurf_rho(file_names2[i], regridder_u2rho, regridder_v2rho)
    
    # Get u, v on 1 m above seafloor
    u_roms_1m_frombot_tmp, v_roms_1m_frombot_tmp = interp_roms_output_to1m_rho(file_names2[i], regridder_u2rho, regridder_v2rho)
    
    # Get salt 1 m from surface and from seafloor 
    salt_roms_1m_fromsurf_tmp, salt_roms_1m_frombot_tmp = interp_salt_to1m(file_names2[i])
    
    # Save these to the arrays 
    #print('time_step: ', time_step)
    #print('time_step + time_lengths[i]: ', time_step+time_lengths[i])
    start = int(time_step2)
    end = int(time_step2+time_lengths[i])
    # All sediment classes added together
    u_1m_fromsurf_rho[start:end,:,:] = u_roms_1m_fromsurf_tmp
    v_1m_fromsurf_rho[start:end,:,:] = v_roms_1m_fromsurf_tmp
    u_1m_frombot_rho[start:end,:,:] = u_roms_1m_frombot_tmp
    v_1m_frombot_rho[start:end,:,:] = v_roms_1m_frombot_tmp
    salt_1m_fromsurf[start:end,:,:] = salt_roms_1m_fromsurf_tmp
    salt_1m_frombot[start:end,:,:] = salt_roms_1m_frombot_tmp
    
    # Update the base time_step
    time_step2 = time_step2 + time_lengths[i]

# Load in the inteproalted netcdfs 
#u_rho_1m = 
#v_rho_1m = 
#salt_rho_1m = 


# Before averaging over time, find the time series
# of the current magnitude 1m below surface
cur_mag_surf_rho = np.sqrt((u_1m_fromsurf_rho**2)+(v_1m_fromsurf_rho**2))
# Average this over time
cur_mag_surf_rho_avg = np.mean(cur_mag_surf_rho, axis=0)
# Now take the average of these over time 
u_surf_rho_avg = np.mean(u_1m_fromsurf_rho, axis=0)
v_surf_rho_avg = np.mean(v_1m_fromsurf_rho, axis=0)
salt_1m_fromsurf_avg = np.mean(salt_1m_fromsurf, axis=0)


# Before averaging over time, find the time series
# of the current magnitude 1m above seafloor 
cur_mag_1m_rho = np.sqrt((u_1m_frombot_rho**2)+(v_1m_frombot_rho**2))
# Average this over time
cur_mag_1m_rho_avg = np.mean(cur_mag_1m_rho, axis=0)

# Now take the average of these over time 
u_1m_rho_avg = np.mean(u_1m_frombot_rho, axis=0)
v_1m_rho_avg = np.mean(v_1m_frombot_rho, axis=0)
salt_1m_frombot_avg = np.mean(salt_1m_frombot, axis=0)



#input('press enter to continue...')






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

# Prep the data by multiplying by the mask and trimming
# Multiply by mask
h_masked = grid.h.values*grid.mask_rho.values*mask_rho_nan.nudge_mask_rho_nan
# Trim 
lon_rho_trimmed = grid.lon_rho[:,c_west:-c_west].values
lat_rho_trimmed = grid.lat_rho[:,c_west:-c_west].values
h_masked_trimmed = h_masked[:,c_west:-c_west]

# Mask, trim, slice
# Mask
# Surface 
u_surf_rho_avg_wland_masked = u_surf_rho_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
v_surf_rho_avg_wland_masked = v_surf_rho_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
salt_surf_avg_wland_masked = salt_1m_fromsurf_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 1 m 
u_1m_rho_avg_wland_masked = u_1m_rho_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
v_1m_rho_avg_wland_masked = v_1m_rho_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
salt_1m_avg_wland_masked = salt_1m_frombot_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Depth-averaged 
ubar_rho_avg_wland_masked = ubar_rho_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
vbar_rho_avg_wland_masked = vbar_rho_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask
salt_depthavg_avg_wland_masked = salt_depthavg_avg*mask_rho_nan.nudge_mask_rho_nan*temp_mask

# Trim
# Surface 
u_surf_rho_avg_wland_masked_trimmed = u_surf_rho_avg_wland_masked[:,c_west:-c_west]
v_surf_rho_avg_wland_masked_trimmed = v_surf_rho_avg_wland_masked[:,c_west:-c_west]
salt_surf_avg_wland_masked_trimmed = salt_surf_avg_wland_masked[:,c_west:-c_west]
# 1 m 
u_1m_rho_avg_wland_masked_trimmed = u_1m_rho_avg_wland_masked[:,c_west:-c_west]
v_1m_rho_avg_wland_masked_trimmed = v_1m_rho_avg_wland_masked[:,c_west:-c_west]
salt_1m_avg_wland_masked_trimmed = salt_1m_avg_wland_masked[:,c_west:-c_west]
# Depth-averaged  
ubar_rho_avg_wland_masked_trimmed = ubar_rho_avg_wland_masked[:,c_west:-c_west]
vbar_rho_avg_wland_masked_trimmed = vbar_rho_avg_wland_masked[:,c_west:-c_west]
salt_depthavg_avg_wland_masked_trimmed = salt_depthavg_avg_wland_masked[:,c_west:-c_west]

# Slice
nth_slice = 20 #15
# Surface 
u_surf_rho_avg_wland_masked_trimmed_slice = u_surf_rho_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
v_surf_rho_avg_wland_masked_trimmed_slice = v_surf_rho_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
# 1 m
u_1m_rho_avg_wland_masked_trimmed_slice = u_1m_rho_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
v_1m_rho_avg_wland_masked_trimmed_slice = v_1m_rho_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
# Depth-averaged 
ubar_rho_avg_wland_masked_trimmed_slice = ubar_rho_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
vbar_rho_avg_wland_masked_trimmed_slice = vbar_rho_avg_wland_masked_trimmed[::nth_slice,::nth_slice]
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
cmap1=cmocean.cm.haline
#cmap1.set_under('darkgray')


# Make the figure
fig1, ax1 = plt.subplots(3, figsize=(21, 18)) # (21,23)

# Set colorbar levels for all plots currents
lev1 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots salinity 
lev2 = np.arange(20,32,0.5)

# Surface 
# In the grid's u and v directions
# Plot bathymetry
# Plot surface salinity
ax1[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax1[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs1 = ax1[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_surf_avg_wland_masked_trimmed, lev2, cmap=cmap1, extend='max')
# Plot currents
q1 = ax1[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, u_surf_rho_avg_wland_masked_trimmed_slice, 
                   v_surf_rho_avg_wland_masked_trimmed_slice, color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax1[0].quiverkey(q1, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours
ax1[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev1, colors='orange')

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax1[0].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m currents
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax1[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax1[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs2 = ax1[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_1m_avg_wland_masked_trimmed, lev2, cmap=cmap1, extend='max')
# Plot currents
q2 = ax1[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   u_1m_rho_avg_wland_masked_trimmed_slice, v_1m_rho_avg_wland_masked_trimmed_slice, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax1[1].quiverkey(q2, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax1[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev1, colors='orange')

# Label the plot
plt.setp(ax1[1].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax1[1].set_ylabel('Y (km)', fontsize=fontsize)

# Depth-averaged currents
# In the grid's u and v directions
# Plot depth-avg salinity 
ax1[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax1[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs3 = ax1[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_depthavg_avg_wland_masked_trimmed, lev2, cmap=cmap1, extend='max')
# Plot currents
q3 = ax1[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ubar_rho_avg_wland_masked_trimmed_slice, vbar_rho_avg_wland_masked_trimmed_slice, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax1[2].quiverkey(q3, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax1[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev1, colors='orange')
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
cbar1 = plt.colorbar(cs1, ax=axes, cax=cbar1_ax, orientation='vertical').set_label(label='Salinity (PSU)', size=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# --------------------------------------------------------------------------------------
# ----- Plot 2: Time-Averaged Surface, 1m above seafloor, and depth-averaged ------------
# ------------------- currents and Magnitude, Masked XY ----------------------------------
# --------------------------------------------------------------------------------------
# Plot time-averaged surface, 1m above seafloor, and depth-averaged currents as 
# quivers and current magnitude as contours in background 

# Prep the data by  ultiplying by the mask and trimming
# Mask, trim, slice
# Mask
# Surface 
cur_mag_surf_rho_avg_masked = cur_mag_surf_rho_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 1 m 
cur_mag_1m_rho_avg_masked = cur_mag_1m_rho_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Depth-averaged 
cur_mag_davg_rho_avg_masked = cur_mag_davg_rho_avg*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask

# Trim
# Surface 
cur_mag_surf_rho_avg_masked_trimmed = cur_mag_surf_rho_avg_masked[:,c_west:-c_west]
# 1 m
cur_mag_1m_rho_avg_masked_trimmed = cur_mag_1m_rho_avg_masked[:,c_west:-c_west]
# Depth-averaged 
cur_mag_davg_rho_avg_masked_trimmed = cur_mag_davg_rho_avg_masked[:,c_west:-c_west]


# Set the colormap
cmap2=cmocean.cm.tempo
#cmap1.set_under('darkgray')


# Make the figure
fig2, ax2 = plt.subplots(3, figsize=(21, 18)) # (21,23)

# Set colorbar levels for all plots currents
lev3 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots current magnitude
lev4 = np.arange(0, 0.4, 0.01)

# Surface 
# In the grid's u and v directions
# Plot bathymetry
# Plot surface salinity
ax2[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax2[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs4 = ax2[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_surf_rho_avg_masked_trimmed, lev4, cmap=cmap2, extend='max')
# Plot currents
q4 = ax2[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, u_surf_rho_avg_wland_masked_trimmed_slice, 
                   v_surf_rho_avg_wland_masked_trimmed_slice, color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax2[0].quiverkey(q4, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours
ax2[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev3, colors='orange')

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax2[0].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax2[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m currents
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax2[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax2[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs5 = ax2[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_1m_rho_avg_masked_trimmed, lev4, cmap=cmap2, extend='max')
# Plot currents
q5 = ax2[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   u_1m_rho_avg_wland_masked_trimmed_slice, v_1m_rho_avg_wland_masked_trimmed_slice, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax2[1].quiverkey(q5, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax2[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev3, colors='orange')

# Label the plot
plt.setp(ax2[1].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax2[1].set_ylabel('Y (km)', fontsize=fontsize)

# Depth-averaged currents
# In the grid's u and v directions
# Plot depth-avg salinity 
ax2[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax2[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs6 = ax2[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_davg_rho_avg_masked_trimmed, lev4, cmap=cmap2, extend='max')
# Plot currents
q6 = ax2[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ubar_rho_avg_wland_masked_trimmed_slice, vbar_rho_avg_wland_masked_trimmed_slice, 
                   color='purple', linewidth=2,
                   angles='xy', scale_units='xy', units='xy') #deeppink
ax2[2].quiverkey(q6, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax2[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev3, colors='orange')
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax2[2].get_xticklabels(), visible=True)
ax2[2].set_xlabel('X (km)', fontsize=fontsize)
ax2[2].set_ylabel('Y (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig2.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)

axes = ax2.ravel().tolist()
cbar2_ax = fig2.add_axes([0.85, bottom, 0.05, top-bottom])
cbar2 = plt.colorbar(cs4, ax=axes, cax=cbar2_ax, orientation='vertical').set_label(label='Current Magnitude (m/s)', size=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)



# --------------------------------------------------------------------------------------
# ----- Plot 3: Time-Averaged Surface, 1m above seafloor, and depth-averaged ------------
# ------------------- salinity, Masked XY ----------------------------------
# --------------------------------------------------------------------------------------
# Plot time-averaged surface, 1m above seafloor, and depth-averaged 
# salinity as contours in background 


# Set the colormap
cmap3=cmocean.cm.haline
#cmap1.set_under('darkgray')


# Make the figure
fig3, ax3 = plt.subplots(3, figsize=(21, 18)) # (21,23)

# Set colorbar levels for all plots currents
lev5 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots salinity 
lev6 = np.arange(20,32,0.5)

# Surface 
# In the grid's u and v directions
# Plot bathymetry
# Plot surface salinity
ax3[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax3[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs7 = ax3[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_surf_avg_wland_masked_trimmed, lev6, cmap=cmap3, extend='max')
# Plot bathymetry contours
ax3[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='orange')

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax3[0].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax3[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m currents
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax3[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax3[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs8 = ax3[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_1m_avg_wland_masked_trimmed, lev6, cmap=cmap3, extend='max')
# Plot bathymetry contours 
ax3[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='orange')

# Label the plot
plt.setp(ax3[1].get_xticklabels(), visible=False)
#ax9[2].set_xlabel('Longitude (degrees)', fontsize=fontsize)
ax3[1].set_ylabel('Y (km)', fontsize=fontsize)

# Depth-averaged currents
# In the grid's u and v directions
# Plot depth-avg salinity 
ax3[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax3[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs9 = ax3[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_depthavg_avg_wland_masked_trimmed, lev6, cmap=cmap3, extend='max')
# Plot bathymetry contours 
ax3[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='orange')
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax3[2].get_xticklabels(), visible=True)
ax3[2].set_xlabel('X (km)', fontsize=fontsize)
ax3[2].set_ylabel('Y (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Adjust colorbar placement
bottom, top = 0.1, 0.95
left, right = 0.1, 0.8
fig3.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)

axes = ax3.ravel().tolist()
cbar3_ax = fig3.add_axes([0.85, bottom, 0.05, top-bottom])
cbar3 = plt.colorbar(cs7, ax=axes, cax=cbar3_ax, orientation='vertical').set_label(label='Salinity (PSU)', size=fontsize)

# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
plt.text(0.775, 0.924, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.634, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.775, 0.346, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# --------------------------------------------------------------------------------------
# ----------------------------- Wave Orbital Analysis ----------------------------------
# --------------------------------------------------------------------------------------
# Do some wave orbital analysis to make a plot with that tagged onto the end 
# Read in the wave forcing file
wave_frc = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Model_Inputs/Forcing_files/wave_forcing_file_kaktovik_shelf_ww3_2020_data002.nc')
print(wave_frc)


# Take the average of the significant wave height over time
swh_avg = wave_frc.Hwave.mean(dim=('wave_time'))

# Take the average of the bottom wave orbital velocity over time
uwave_avg = wave_frc.Uwave_rms.mean(dim=('wave_time'))

# Make land show
temp_mask = grid.mask_rho
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)
# Significant weave height 
swh_avg_wland = swh_avg * temp_mask
swh_avg_wland = swh_avg_wland.fillna(-100)
# Bottom wave orbital velocity
uwave_avg_wland = uwave_avg * temp_mask
uwave_avg_wland = uwave_avg_wland.fillna(-100)

# Make a plot with this as a fourth subplot of the current plot 
# --------------------------------------------------------------------------------------
# ----- Plot 4: Time-Averaged Surface, 1m above seafloor, and depth-averaged ------------
# ----------- currents and Magnitude, Masked XY Plue Wave Orbitals --------------------
# --------------------------------------------------------------------------------------
# Plot time-averaged surface, 1m above seafloor, and depth-averaged currents as 
# quivers and current magnitude as contours in background with a fourth subplot
# with the time-averaged wave orbital velocities 

# Make it so land can be masked easily 
# Make it so land will be gray
temp_mask = grid.mask_rho.values
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)
# mud01
uwave_avg_wland2 = uwave_avg * temp_mask

# Prep the data by  ultiplying by the mask and trimming
# Multiply by mask
h_masked = grid.h.values*grid.mask_rho.values*mask_rho_nan.nudge_mask_rho_nan
uwave_avg_wland_masked = uwave_avg_wland*mask_rho_nan.nudge_mask_rho_nan
uwave_avg_wland2_masked = uwave_avg_wland2*mask_rho_nan.nudge_mask_rho_nan
# Trim 
lon_rho_trimmed = grid.lon_rho[:,c_west:-c_west].values
lat_rho_trimmed = grid.lat_rho[:,c_west:-c_west].values
h_masked_trimmed = h_masked[:,c_west:-c_west]
uwave_avg_wland_masked_trimmed = uwave_avg_wland_masked[:,c_west:-c_west]
uwave_avg_wland2_masked_trimmed = uwave_avg_wland2_masked[:,c_west:-c_west]


# Set the colormap
cmap4=cmocean.cm.algae # current magnitude
cmap5 = cmocean.cm.ice_r
#cmap1.set_under('darkgray')


# Make the figure
fig4, ax4 = plt.subplots(5, figsize=(18, 25)) # (21,23) (21,18)

# Set colorbar levels for all plots currents
lev7 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots current magnitude
lev8 = np.arange(0, 0.4, 0.01)
# Set the colorbar for orbital wave velocities 
lev9 = np.arange(0,.16,0.001)

# Surface 
# In the grid's u and v directions
# Plot bathymetry
# Plot surface salinity
ax4[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs10 = ax4[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_surf_rho_avg_masked_trimmed, lev8, cmap=cmap4, extend='max')
# Plot currents
q7 = ax4[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, u_surf_rho_avg_wland_masked_trimmed_slice, 
                   v_surf_rho_avg_wland_masked_trimmed_slice, color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax4[0].quiverkey(q7, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours
ax4[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orange')

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax4[0].get_xticklabels(), visible=False)
#ax4[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m currents
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax4[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs11 = ax4[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_1m_rho_avg_masked_trimmed, lev8, cmap=cmap4, extend='max')
# Plot currents
q8 = ax4[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   u_1m_rho_avg_wland_masked_trimmed_slice, v_1m_rho_avg_wland_masked_trimmed_slice, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax4[1].quiverkey(q8, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax4[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orange')

# Label the plot
plt.setp(ax4[1].get_xticklabels(), visible=False)
#ax4[1].set_ylabel('Y (km)', fontsize=fontsize)

# Depth-averaged currents
# In the grid's u and v directions
# Plot depth-avg salinity 
ax4[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs12 = ax4[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_davg_rho_avg_masked_trimmed, lev8, cmap=cmap4, extend='max')
# Plot currents
q9 = ax4[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ubar_rho_avg_wland_masked_trimmed_slice, vbar_rho_avg_wland_masked_trimmed_slice, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy') #deeppink
ax4[2].quiverkey(q9, 0.65, 0.85, U=0.1, label='0.1 m/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax4[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orange')
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax4[2].get_xticklabels(), visible=False)
ax4[2].set_ylabel('Y (km)', fontsize=fontsize)
#ax4[2].set_xlabel('X (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)
# Adjust colorbar placement
bottom, top = 0.1, 0.9
left, right = 0.1, 0.84
fig4.subplots_adjust(top=top, bottom=bottom, left=left, right=right, wspace=0.2, hspace=0.2)
#axes = ax4.ravel().tolist()
bottom4, top4 = 0.425, 0.9
cbar4_ax = fig4.add_axes([0.85, bottom4, 0.04, top4-bottom4])
cbar4 = fig4.colorbar(cs10, cax=cbar4_ax, ax=[ax4[0],ax4[1],ax4[2]], orientation='vertical', pad=0.015).set_label(label='Current Magnitude (m/s)', size=fontsize, labelpad=15)


# Surface salinity 
# Plot surface salinity
ax4[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[3].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs13 = ax4[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_surf_avg_wland_masked_trimmed, lev6, cmap=cmap3, extend='max')
# Plot bathymetry contours
ax4[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='orange')
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax4[3].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
#ax4[3].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)
bottom5, top5 = 0.265, 0.41
cbar5_ax = fig4.add_axes([0.85, bottom5, 0.04, top5-bottom5])
cbar5 = fig4.colorbar(cs13, cax=cbar5_ax, ax=ax4[3], orientation='vertical', 
                     pad=0.015).set_label(label='Surface \nSalinity (PSU)', size=fontsize, labelpad=15)


# Bottom wave orbital velocities
# Make it so land will be gray and ocean white 
ax4[4].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax4[4].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs14 = ax4[4].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, uwave_avg_wland_masked_trimmed, lev9, cmap=cmap5, extend='max')
cs15 = ax4[4].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orangered')
# add title
#ax1.set_title('Time-Averaged Significant Wave Height', fontsize=fontsize, y=1.08)
# format axes
plt.setp(ax4[4].get_xticklabels(), visible=True)
ax4[4].set_xlabel('X (km)', fontsize=fontsize)
#ax4[4].set_ylabel('Y (km)', fontsize=fontsize)
# specify colorbar
bottom6, top6 = 0.1, 0.25
cbar6_ax = fig4.add_axes([0.85, bottom6, 0.04, top6-bottom6])
cbar6 = fig4.colorbar(cs14, cax=cbar6_ax, ax=ax4[4], label='Bottom Wave Orbital Velocity (m/s)', orientation='vertical', 
                      pad=0.015, ticks=[0,0.032,0.064,0.096,0.128])
cbar6.set_label(label='Bottom Wave \nOrbital Velocity (m/s)', size=fontsize, labelpad=15)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
plt.text(0.810, 0.879, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.810, 0.718, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.810, 0.554, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.810, 0.390, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.810, 0.229, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)




# MRemake plot above but semi differently 
# --------------------------------------------------------------------------------------
# ----- Plot 5: Time-Averaged Surface, 1m above seafloor, and depth-averaged ------------
# ----------- currents and Magnitude, Masked XY Plue Wave Orbitals --------------------
# --------------------------------------------------------------------------------------
# Plot time-averaged surface, 1m above seafloor, and depth-averaged currents as 
# quivers and current magnitude as contours in background with a fourth subplot
# with the time-averaged wave orbital velocities 

# Set the colormap
cmap4=cmocean.cm.algae # current magnitude
cmap5 = cmocean.cm.ice_r
#cmap1.set_under('darkgray')


# Make the figure
fig5, ax5 = plt.subplots(5, figsize=(18, 29)) # (18,25) (21,23) (21,18)

# Set colorbar levels for all plots currents
lev7 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Set colorbar levels for all plots current magnitude
lev8b = np.arange(0, 40, 1)
# Set the colorbar for orbital wave velocities 
lev9b = np.arange(0,16,0.1)

# Depth-averaged currents
# In the grid's u and v directions
# Plot depth-avg salinity 
ax5[0].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[0].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs12 = ax5[0].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_davg_rho_avg_masked_trimmed*100, lev8b, cmap=cmap4, extend='max')
# Plot currents
q9 = ax5[0].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   ubar_rho_avg_wland_masked_trimmed_slice*100, vbar_rho_avg_wland_masked_trimmed_slice*100, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy') #deeppink
ax5[0].quiverkey(q9, 0.93, 0.09, U=10, label='10 cm/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax5[0].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orange')
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax5[0].get_xticklabels(), visible=False)
ax5[0].set_ylabel('Y (km)', fontsize=fontsize)
#ax4[2].set_xlabel('X (km)', fontsize=fontsize)
#cbar10 = plt.colorbar(cs10, orientation='vertical', ax=ax8[2]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)
# If horizontal and in axes
# No depth indication 
#cbar5c = plt.colorbar(cs12, cax=ax5[0].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
 #                    ticks=[0, 10, 20, 30, 40], ax=ax5[0], 
  #                   orientation='horizontal').set_label(label='Current Magnitude (cm/s)', size=fontsize-2)
# With depth indication 
cbar5c = plt.colorbar(cs12, cax=ax5[0].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 10, 20, 30, 40], ax=ax5[0], 
                     orientation='horizontal').set_label(label='Depth-Avg. Cur. Mag. (cm/s)', size=fontsize-2)

# Surface 
# In the grid's u and v directions
# Plot bathymetry
# Plot surface salinity
ax5[1].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[1].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs10 = ax5[1].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_surf_rho_avg_masked_trimmed*100, lev8b, cmap=cmap4, extend='max')
# Plot currents
q7 = ax5[1].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, u_surf_rho_avg_wland_masked_trimmed_slice*100, 
                   v_surf_rho_avg_wland_masked_trimmed_slice*100, color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax5[1].quiverkey(q7, 0.93, 0.09, U=10, label='10 cm/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours
ax5[1].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orange')
ax5[1].set_ylabel('Y (km)', fontsize=fontsize)
# If horizontal and in axes
# No depth indication 
#cbar5a = plt.colorbar(cs10, cax=ax5[1].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
 #                    ticks=[0, 10, 20, 30, 40], ax=ax5[1], 
  #                   orientation='horizontal').set_label(label='Current Magnitude (cm/s)', size=fontsize-2)
# With depth indication 
cbar5a = plt.colorbar(cs10, cax=ax5[1].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 10, 20, 30, 40], ax=ax5[1], 
                     orientation='horizontal').set_label(label='Surface Cur. Mag. (cm/s)', size=fontsize-2)

# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax5[1].get_xticklabels(), visible=False)
#ax4[0].set_ylabel('Y (km)', fontsize=fontsize)
#cbar8 = plt.colorbar(cs8, orientation='vertical', ax=ax8[0]).set_label(label='Current \nMagnitude (m/s)', size=fontsize)

# Bottom  1m currents
# In the grid's u and v directions
# Plot salinity at 1 m above seafloor 
ax5[2].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[2].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs11 = ax5[2].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  cur_mag_1m_rho_avg_masked_trimmed*100, lev8b, cmap=cmap4, extend='max')
# Plot currents
q8 = ax5[2].quiver(x_rho_flat_trimmed_slice/1000, y_rho_flat_slice/1000, 
                   u_1m_rho_avg_wland_masked_trimmed_slice*100, v_1m_rho_avg_wland_masked_trimmed_slice*100, 
                   color='deeppink', linewidth=2,
                   angles='xy', scale_units='xy', units='xy')
ax5[2].quiverkey(q8, 0.93, 0.09, U=10, label='10 cm/s', fontproperties={'size':fontsize-2})
# Plot bathymetry contours 
ax5[2].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orange')
ax5[2].set_ylabel('Y (km)', fontsize=fontsize)
# If horizontal and in axes
# No depth indication 
#cbar5b = plt.colorbar(cs11, cax=ax5[2].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
 #                    ticks=[0, 10, 20, 30, 40], ax=ax5[2], 
  #                   orientation='horizontal').set_label(label='Current Magnitude (cm/s)', size=fontsize-2)
# With depth indication 
cbar5b = plt.colorbar(cs11, cax=ax5[2].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 10, 20, 30, 40], ax=ax5[2], 
                     orientation='horizontal').set_label(label='Bottom Cur. Mag. (cm/s)', size=fontsize-2, loc='right')

# Label the plot
plt.setp(ax5[2].get_xticklabels(), visible=False)
#ax4[1].set_ylabel('Y (km)', fontsize=fontsize)


# Surface salinity 
# Plot surface salinity
ax5[3].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[3].fill_between(x_rho_flat_trimmed/1000, 65 ,120, 
               facecolor ='white', alpha = 0.8)
cs13 = ax5[3].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000,
                  salt_surf_avg_wland_masked_trimmed, lev6, cmap=cmap3, extend='max')
# Plot bathymetry contours
ax5[3].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev5, colors='orange')
# Label the plot
#ax7.set_title('Time-Averaged Surface Currents (m/s)', fontsize=fontsize, y=1.08)
plt.setp(ax5[3].get_xticklabels(), visible=False)
#ax8[0].set_xlabel('Longitude (degrees)', fontsize=fontsize)
#ax4[3].set_ylabel('Y (km)', fontsize=fontsize)
ax5[3].set_ylabel('Y (km)', fontsize=fontsize)
# If horizontal and in axes
cbar8 = plt.colorbar(cs13, cax=ax5[3].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[20, 23, 26, 29, 32], ax=ax5[3], 
                     orientation='horizontal').set_label(label='Surface Salinity (PSU)', size=fontsize-2)


# Bottom wave orbital velocities
# Make it so land will be gray and ocean white 
ax5[4].fill_between(x_rho_flat_trimmed/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax5[4].fill_between(x_rho_flat_trimmed/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs14 = ax5[4].contourf(x_rho_flat_trimmed/1000, y_rho_flat/1000, uwave_avg_wland_masked_trimmed*100, lev9b, cmap=cmap5, extend='max')
cs15 = ax5[4].contour(x_rho_flat_trimmed/1000, y_rho_flat/1000, h_masked_trimmed, lev7, colors='orangered')
# add title
#ax1.set_title('Time-Averaged Significant Wave Height', fontsize=fontsize, y=1.08)
# format axes
plt.setp(ax5[4].get_xticklabels(), visible=True)
ax5[4].set_xlabel('X (km)', fontsize=fontsize)
#ax4[4].set_ylabel('Y (km)', fontsize=fontsize)
ax5[4].set_ylabel('Y (km)', fontsize=fontsize)
# If horizontal and in axes
cbar9 = plt.colorbar(cs14, cax=ax5[4].inset_axes((0.55, 0.92, 0.45, 0.05)),      #(0.38, 0.92, 0.6, 0.05)
                     ticks=[0, 4, 8, 12, 16], ax=ax5[4], 
                     orientation='horizontal').set_label(label='Bottom Wave Orbital Velocity (cm/s)', size=fontsize-2)


# Adjust spacing between subplots
plt.subplots_adjust(hspace=0.08) #0.08

# Add subplot labels
# Top right corner 
#plt.text(0.810, 0.879, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.810, 0.718, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.810, 0.554, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.810, 0.390, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
#plt.text(0.810, 0.229, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
# Bottom left corner 
plt.text(0.135, 0.744, 'a)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.135, 0.589, 'b)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.135, 0.439, 'c)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.135, 0.284, 'd)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)
plt.text(0.135, 0.135, 'e)', fontsize=fontsize, fontweight='bold', color='k', transform=plt.gcf().transFigure)


# --------------------------------------------------------------------------------------
# --------------------------- Plot 6: Lowest salinity ----------------------------------
# --------------------------------------------------------------------------------------
# Find the river plume extent by plotting the lowest salinity over all time
# at each spot the look at where it is below a certain threshold

# Find the lowest salinity

# Get the full surface salinity over all time and space 
surf_salt_all = salt_full[:,-1,:,:]

# Find the minimum value over time 
surf_salt_min_all_time = np.min(surf_salt_all, axis=0)

# Plot this 
fig6, ax6 = plt.subplots(figsize=(22,8))
# Make it so land will be gray and ocean white 
ax6.fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax6.fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs16 = ax6.contourf(x_rho_flat/1000, y_rho_flat/1000, surf_salt_min_all_time, lev6, cmap=cmap3)
cbar10 = plt.colorbar(cs16).set_label('Salinity (PSU)', fontsize=fontsize)
ax6.set_title('Minimum Salinity Over All Time', fontsize=fontsize)
ax6.set_xlabel('X (km)', fontsize=fontsize)
ax6.set_ylabel('Y (km)', fontsize=fontsize)

# Find where the salinity is less than or equal to a threshold
# Set min river salimity
#min_salt_river_plume = 27.0 # PSU

# Find the indices where this is true
#min_salt_river_plume_idx = np.where(surf_salt_min_all_time <= min_salt_river_plume)


# Plot threshold of 27 PSU
# Set colorbar levels for all plots salinity 
lev7 = np.arange(20,27,0.5)

fig7, ax7 = plt.subplots(figsize=(22,8))
# Make it so land will be gray and ocean white 
ax7.fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax7.fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs17 = ax7.contourf(x_rho_flat/1000, y_rho_flat/1000, surf_salt_min_all_time, lev7, cmap=cmap3)
ax7.scatter(x_rho_flat[217]/1000, y_rho_flat[95]/1000, marker='x', color='red', s=250)
ax7.scatter(x_rho_flat[217]/1000, y_rho_flat[64]/1000, marker='x', color='red', s=250)
cbar17 = plt.colorbar(cs17).set_label('Salinity (PSU)', fontsize=fontsize)
ax7.set_title('Minimum Salinity Over All Time', fontsize=fontsize)
ax7.set_xlabel('X (km)', fontsize=fontsize)
ax7.set_ylabel('Y (km)', fontsize=fontsize)


# Plot threshold of 25 PSU
# Set colorbar levels for all plots salinity 
lev8 = np.arange(20,25,0.5)

fig8, ax8 = plt.subplots(figsize=(22,8))
# Make it so land will be gray and ocean white 
ax8.fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax8.fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs18 = ax8.contourf(x_rho_flat/1000, y_rho_flat/1000, surf_salt_min_all_time, lev8, cmap=cmap3)
ax8.scatter(x_rho_flat[217]/1000, y_rho_flat[89]/1000, marker='x', color='red', s=250)
ax8.scatter(x_rho_flat[217]/1000, y_rho_flat[64]/1000, marker='x', color='red', s=250)
cbar18 = plt.colorbar(cs18).set_label('Salinity (PSU)', fontsize=fontsize)
ax8.set_title('Minimum Salinity Over All Time', fontsize=fontsize)
ax8.set_xlabel('X (km)', fontsize=fontsize)
ax8.set_ylabel('Y (km)', fontsize=fontsize)


# Caluclate the distance (dy = 600 m)
dy = 600 # m

# For a threshold of 27 PSU
river_plume_dist = dy*(95-64)
print('Maximum distance of river plume extent (m) (lowest salinity of 27 PSU): ', river_plume_dist)

# For a threshold of 25 PSU
river_plume_dist2 = dy*(89-64)
print('Maximum distance of river plume extent (m) (lowest salinity of 25 PSU): ', river_plume_dist2)


# Make and print some calculations for currents
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
# 10 - 40 m depth
h_masked6 = h_masked.copy()
outer_10_40m_mask_rho = masked_array_lowhigh_2dloop(h_masked6, 10, 40)

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
# 10 - 40 m depth
outer_10_40m_mask_rho_nan_idx = np.where(outer_10_40m_mask_rho == 0.0)
outer_10_40m_mask_rho_nan = outer_10_40m_mask_rho.copy()
outer_10_40m_mask_rho_nan = outer_10_40m_mask_rho_nan.astype('float')
outer_10_40m_mask_rho_nan[outer_10_40m_mask_rho_nan_idx] = np.nan

# Now multiply by the mask to get the different regions 
# We already have average current magnitude so just multiply that by the masks 
# Just do depth-averaged current magnitude for now
# Inner 
cur_mag_davg_rho_masked_inner = cur_mag_davg_rho_avg * inner_shelf_mask_rho_nan
# Mid
cur_mag_davg_rho_masked_mid = cur_mag_davg_rho_avg * mid_shelf_mask_rho_nan
# Outer
cur_mag_davg_rho_masked_outer = cur_mag_davg_rho_avg * outer_shelf_mask_rho_nan
# 0 - 10 m
uwave_rms_avg_masked_10m = uwave_avg * inner_10m_mask_rho_nan
# 10 - 60 m
uwave_rms_avg_masked_10_60m = uwave_avg * outer_10_60m_mask_rho_nan
# 0 - 10 m salinity 
salt_surf_avg_masked_10m = salt_1m_fromsurf_avg * inner_10m_mask_rho_nan
# 10 - 40 m salinity
salt_surf_avg_masked_10_40m = salt_1m_fromsurf_avg * outer_10_40m_mask_rho_nan
# 40 - 60 m salinity
salt_surf_avg_masked_40_60m = salt_1m_fromsurf_avg * outer_shelf_mask_rho_nan


# Mask and trim these so that we are just looking at the regions in the plot
# Mask, trim, slice
# Mask
# Inner
cur_mag_davg_rho_masked_inner_masked = cur_mag_davg_rho_masked_inner*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask 
# Mid
cur_mag_davg_rho_masked_mid_masked = cur_mag_davg_rho_masked_mid*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# Outer
cur_mag_davg_rho_masked_outer_masked = cur_mag_davg_rho_masked_outer*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 0 - 10 m
uwave_rms_avg_masked_10m_masked = uwave_rms_avg_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 60 m
uwave_rms_avg_masked_10_60m_masked = uwave_rms_avg_masked_10_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 0 - 10 m salinity 
salt_surf_avg_masked_10m_masked = salt_surf_avg_masked_10m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 10 - 40 m salinity
salt_surf_avg_masked_10_40m_masked = salt_surf_avg_masked_10_40m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask
# 40 - 60 m salinity
salt_surf_avg_masked_40_60m_masked = salt_surf_avg_masked_40_60m*temp_mask*mask_rho_nan.nudge_mask_rho_nan*temp_mask

# Trim
# Inner
cur_mag_davg_rho_masked_inner_masked_trimmed = cur_mag_davg_rho_masked_inner_masked[:,c_west:-c_west]
# 1 m
cur_mag_davg_rho_masked_mid_masked_trimmed = cur_mag_davg_rho_masked_mid_masked[:,c_west:-c_west]
# Depth-averaged 
cur_mag_davg_rho_masked_outer_masked_trimmed = cur_mag_davg_rho_masked_outer_masked[:,c_west:-c_west]
# 0 - 10 m
uwave_rms_avg_masked_10m_masked_trimmed = uwave_rms_avg_masked_10m_masked[:,c_west:-c_west]
# 10 - 60 m
uwave_rms_avg_masked_10_60m_masked_trimmed = uwave_rms_avg_masked_10_60m_masked[:,c_west:-c_west]
# 0 - 10 m salinity 
salt_surf_avg_masked_10m_masked_trimmed = salt_surf_avg_masked_10m_masked[:,c_west:-c_west]
# 10 - 40 m salinity
salt_surf_avg_masked_10_40m_masked_trimmed = salt_surf_avg_masked_10_40m_masked[:,c_west:-c_west]
# 40 - 60 m salinity
salt_surf_avg_masked_40_60m_masked_trimmed = salt_surf_avg_masked_40_60m_masked[:,c_west:-c_west]


# =============================================================================
# # Replace 0s with nans
# idx_zero = np.where(salt_surf_avg_masked_10m_masked_trimmed == 0)
# salt_surf_avg_masked_10m_masked_trimmed[idx_zero] = np.nan
# =============================================================================


# Print some statistics 
# Mean 
cur_mag_davg_rho_masked_inner_masked_trimmed_avg = np.nanmean(cur_mag_davg_rho_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) averaged current magnitude (m/s): ', cur_mag_davg_rho_masked_inner_masked_trimmed_avg)
cur_mag_davg_rho_masked_mid_masked_trimmed_avg = np.nanmean(cur_mag_davg_rho_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) averaged current magnitude (m/s): ', cur_mag_davg_rho_masked_mid_masked_trimmed_avg)
cur_mag_davg_rho_masked_outer_masked_trimmed_avg = np.nanmean(cur_mag_davg_rho_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) averaged current magnitude (m/s): ', cur_mag_davg_rho_masked_outer_masked_trimmed_avg)
uwave_rms_avg_masked_10m_masked_trimmed_avg = np.nanmean(uwave_rms_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m Uwaverms mean (m/s): ', uwave_rms_avg_masked_10m_masked_trimmed_avg)
uwave_rms_avg_masked_10_60m_masked_trimmed_avg = np.nanmean(uwave_rms_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m Uwaverms mean (m/s): ', uwave_rms_avg_masked_10_60m_masked_trimmed_avg)
salt_surf_avg_masked_10m_masked_trimmed_avg = np.nanmean(salt_surf_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m salinity mean (PSU): ', salt_surf_avg_masked_10m_masked_trimmed_avg)
salt_surf_avg_masked_10_40m_masked_trimmed_avg = np.nanmean(salt_surf_avg_masked_10_40m_masked_trimmed, axis=(0,1))
print('10 - 40 m salinity mean (PSU): ', salt_surf_avg_masked_10_40m_masked_trimmed_avg)
salt_surf_avg_masked_40_60m_masked_trimmed_avg = np.nanmean(salt_surf_avg_masked_40_60m_masked_trimmed, axis=(0,1))
print('40 - 60 m salinity mean (PSU): ', salt_surf_avg_masked_40_60m_masked_trimmed_avg)

# Standard deviation 
cur_mag_davg_rho_masked_inner_masked_trimmed_std = np.nanstd(cur_mag_davg_rho_masked_inner_masked_trimmed, axis=(0,1))
print('Inner shelf (0-20 m) std current magnitude (m/s): ', cur_mag_davg_rho_masked_inner_masked_trimmed_std)
cur_mag_davg_rho_masked_mid_masked_trimmed_std = np.nanstd(cur_mag_davg_rho_masked_mid_masked_trimmed, axis=(0,1))
print('Mid shelf (20-40 m) std current magnitude (m/s): ', cur_mag_davg_rho_masked_mid_masked_trimmed_std)
cur_mag_davg_rho_masked_outer_masked_trimmed_std = np.nanstd(cur_mag_davg_rho_masked_outer_masked_trimmed, axis=(0,1))
print('Outer shelf (40-60 m) stdcurrent magnitude (m/s): ', cur_mag_davg_rho_masked_outer_masked_trimmed_std)
uwave_rms_avg_masked_10m_masked_trimmed_std = np.nanstd(uwave_rms_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m Uwaverms std (m/s): ', uwave_rms_avg_masked_10m_masked_trimmed_std)
uwave_rms_avg_masked_10_60m_masked_trimmed_std = np.nanstd(uwave_rms_avg_masked_10_60m_masked_trimmed, axis=(0,1))
print('10 - 60 m Uwaverms std (m/s): ', uwave_rms_avg_masked_10_60m_masked_trimmed_std)
salt_surf_avg_masked_10m_masked_trimmed_std = np.nanstd(salt_surf_avg_masked_10m_masked_trimmed, axis=(0,1))
print('0 - 10 m salinity std (PSU): ', salt_surf_avg_masked_10m_masked_trimmed_std)
salt_surf_avg_masked_10_40m_masked_trimmed_std = np.nanstd(salt_surf_avg_masked_10_40m_masked_trimmed, axis=(0,1))
print('10 - 40 m salinity std (PSU): ', salt_surf_avg_masked_10_40m_masked_trimmed_std)
salt_surf_avg_masked_40_60m_masked_trimmed_std = np.nanstd(salt_surf_avg_masked_40_60m_masked_trimmed, axis=(0,1))
print('40 - 60 m salinity std (PSU): ', salt_surf_avg_masked_40_60m_masked_trimmed_std)



# -------------------------------------------------------------------------------
# ---- Make a netcdf to hold the output data used for plotting 
# -------------------------------------------------------------------------------
# Set up the data
roms_depth_avg_surf_bot_cur_salt_uwave = xr.Dataset(
    data_vars=dict(
        roms_depth_avg_currents_mag_time_avg=(['y','x'], cur_mag_davg_rho_avg_masked_trimmed.values),
        roms_ubar_time_avg=(['y_slice','x_slice'], ubar_rho_avg_wland_masked_trimmed_slice.values),
        roms_vbar_time_avg=(['y_slice','x_slice'], vbar_rho_avg_wland_masked_trimmed_slice.values),
        roms_surf_currents_mag_time_avg=(['y','x'], cur_mag_surf_rho_avg_masked_trimmed.values),
        roms_surf_u_currents_time_avg=(['y_slice','x_slice'], u_surf_rho_avg_wland_masked_trimmed_slice.values),
        roms_surf_v_currents_time_avg=(['y_slice','x_slice'], v_surf_rho_avg_wland_masked_trimmed_slice.values),
        roms_1m_above_seafloor_currents_mag_time_avg=(['y','x'], cur_mag_1m_rho_avg_masked_trimmed.values),
        roms_1m_above_seafloor_u_currents_time_avg=(['y_slice','x_slice'], u_1m_rho_avg_wland_masked_trimmed_slice.values),
        roms_1m_above_seafloor_v_currents_time_avg=(['y_slice','x_slice'], v_1m_rho_avg_wland_masked_trimmed_slice.values),
        roms_surf_salt_time_avg=(['y','x'], salt_surf_avg_wland_masked_trimmed.values),
        roms_uwave_time_avg=(['y','x'], uwave_avg_wland_masked_trimmed.values)
        ),
    coords=dict(
        x_full=('x', x_rho_flat_trimmed),
        x_slice=('x_slice', x_rho_flat_trimmed_slice),
        y_full=('y', y_rho_flat), 
        y_slice=('y_slice', y_rho_flat_slice)
        ),
    attrs=dict(description='Time-averaged ROMS output including dpeth-averaged, surface, and 1 meter above seafloor current magnitudes and vectors (meter per second), surface salinity (PSU), and bottom wave orbital velocities (meter per second)'))
# Add more metadata?
roms_depth_avg_surf_bot_cur_salt_uwave.roms_depth_avg_currents_mag_time_avg.name='depth-averaged, time-averaged current magnitude'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_ubar_time_avg.name='depth-averaged, time-averaged u water momentum'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_vbar_time_avg.name='depth-averaged, time-averaged v water momentum'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_surf_currents_mag_time_avg.name='time-averaged surface current magnitude'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_surf_u_currents_time_avg.name='time-averaged surface u water momentum'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_surf_v_currents_time_avg.name='time-averaged surface v water momentum'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_1m_above_seafloor_currents_mag_time_avg.name='time-averaged current magnitude 1 m above seafloor'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_1m_above_seafloor_u_currents_time_avg.name='time-averaged u water momentum 1 m above seafloor'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_1m_above_seafloor_v_currents_time_avg.name='time-averaged v water momentum 1 m above seafloor'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_surf_salt_time_avg.name='time-averaged surface salinity'
roms_depth_avg_surf_bot_cur_salt_uwave.roms_uwave_time_avg.name='time-averaged bottom wave orbital velocity'
roms_depth_avg_surf_bot_cur_salt_uwave.x_full.name='longitudinal distance (meter)'
roms_depth_avg_surf_bot_cur_salt_uwave.y_full.name='latitudinal distance (meter)'
roms_depth_avg_surf_bot_cur_salt_uwave.x_slice.name='longitudinal distance (meter), sliced for arrows'
roms_depth_avg_surf_bot_cur_salt_uwave.y_slice.name='latitudinal distance (meter), sliced for arrows'
# Save to a netcdf
#roms_depth_avg_surf_bot_cur_salt_uwave.to_netcdf('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2020/Paper1_Take2/Data/fig6_roms_depth_avg_surf_bot_currents_salt_uwave.nc')



# --------------------------------------------------------------------------------------
# -------------------- Plot 7: Max Current Magnitude  ----------------------------------
# --------------------------------------------------------------------------------------
# Find the maximum surface, depth-averaged, an 1 m above seafloor current
# magnitudes for all time 
# Find the minimum value over time 
surf_curmag_max_all_time = np.nanmax(cur_mag_surf_rho, axis=0)
davg_curmag_max_all_time = np.nanmax(cur_mag_davg_rho, axis=0)
bot_curmag_max_all_time = np.nanmax(cur_mag_1m_rho, axis=0)

# Make a subplot of these with bathymetry contours
# Set levels 
lev9 = np.arange(30,120,1)
cmap9 = cmocean.cm.algae
# Set colorbar levels for all plots currents
lev11 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Plot this 
fig9, ax9 = plt.subplots(3, figsize=(22,16))

# Depth-averaged 
# Make it so land will be gray and ocean white 
ax9[0].fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax9[0].fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs19 = ax9[0].contourf(x_rho_flat/1000, y_rho_flat/1000, davg_curmag_max_all_time*100, lev9, cmap=cmap9,
                       extend='both')
# Plot bathymetry contours
ax9[0].contour(x_rho_flat/1000, y_rho_flat/1000, grid.h.values, lev11, colors='deeppink')
ax9[0].set_title('Depth-Averaged Current Magnitude (cm/s)', fontsize=fontsize)
plt.setp(ax9[0].get_xticklabels(), visible=False)

# Surface 
# Make it so land will be gray and ocean white 
ax9[1].fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax9[1].fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs20 = ax9[1].contourf(x_rho_flat/1000, y_rho_flat/1000, surf_curmag_max_all_time*100, lev9, cmap=cmap9,
                       extend='both')
# Plot bathymetry contours
ax9[1].contour(x_rho_flat/1000, y_rho_flat/1000, grid.h.values, lev11, colors='deeppink')
ax9[1].set_title('Surface Current Magnitude (cm/s)', fontsize=fontsize)
plt.setp(ax9[1].get_xticklabels(), visible=False)

# Bottom
# Make it so land will be gray and ocean white 
ax9[2].fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax9[2].fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs21 = ax9[2].contourf(x_rho_flat/1000, y_rho_flat/1000, bot_curmag_max_all_time*100, lev9, cmap=cmap9,
                       extend='both')
# Plot bathymetry contours
ax9[2].contour(x_rho_flat/1000, y_rho_flat/1000, grid.h.values, lev11, colors='deeppink')
ax9[2].set_title('1 m Above Seafloor Current Magnitude (cm/s)', fontsize=fontsize)
ax9[2].set_xlabel('X (km)', fontsize=fontsize)
ax9[2].set_ylabel('Y (km)', fontsize=fontsize)

cbar13 = plt.colorbar(cs21, ax=[ax9[0], ax9[1], ax9[2]], extend='both').set_label('Current Magnitude (cm/s)', fontsize=fontsize)

#fig9.suptitle('Maximum Current Magnitude Over Time (cm/s)', fontsize=fontsize,)
fig9.text(0.5, 0.92, 'Maximum Current Magnitude Over Time', ha="center", va="center", rotation=0, fontsize=fontsize+5)


# --------------------------------------------------------------------------------------
# ----------------- Plot 8: Standard Deviation Current Magnitude  ----------------------
# --------------------------------------------------------------------------------------
# Do the same as above but the standard deviation of the current magnitude over time
# (hold off unless asked for)

# Find the stadndard deviation value over time 
surf_curmag_std_all_time = np.nanstd(cur_mag_surf_rho, axis=0)
davg_curmag_std_all_time = np.nanstd(cur_mag_davg_rho, axis=0)
bot_curmag_std_all_time = np.nanstd(cur_mag_1m_rho, axis=0)

# Make a subplot of these with bathymetry contours
# Set levels 
lev12= np.arange(0,30,0.1)
cmap9 = cmocean.cm.algae
# Set colorbar levels for all plots currents
lev11 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# Plot this 
fig10, ax10 = plt.subplots(3, figsize=(22,16))

# Depth-averaged 
# Make it so land will be gray and ocean white 
ax10[0].fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax10[0].fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs22 = ax10[0].contourf(x_rho_flat/1000, y_rho_flat/1000, davg_curmag_std_all_time*100, lev12, cmap=cmap9,
                       extend='both')
# Plot bathymetry contours
ax10[0].contour(x_rho_flat/1000, y_rho_flat/1000, grid.h.values, lev11, colors='deeppink')
ax10[0].set_title('Depth-Averaged Current Magnitude (cm/s)', fontsize=fontsize)
plt.setp(ax10[0].get_xticklabels(), visible=False)

# Surface 
# Make it so land will be gray and ocean white 
ax10[1].fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax10[1].fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs23 = ax10[1].contourf(x_rho_flat/1000, y_rho_flat/1000, surf_curmag_std_all_time*100, lev12, cmap=cmap9,
                       extend='both')
# Plot bathymetry contours
ax10[1].contour(x_rho_flat/1000, y_rho_flat/1000, grid.h.values, lev11, colors='deeppink')
ax10[1].set_title('Surface Current Magnitude (cm/s)', fontsize=fontsize)
plt.setp(ax10[1].get_xticklabels(), visible=False)

# Bottom
# Make it so land will be gray and ocean white 
ax10[2].fill_between(x_rho_flat/1000, 0, 65, 
               facecolor ='darkgray', alpha = 0.8)
ax10[2].fill_between(x_rho_flat/1000, 65 , 120, 
               facecolor ='white', alpha = 0.8)
cs24 = ax10[2].contourf(x_rho_flat/1000, y_rho_flat/1000, bot_curmag_std_all_time*100, lev12, cmap=cmap9,
                       extend='both')
# Plot bathymetry contours
ax10[2].contour(x_rho_flat/1000, y_rho_flat/1000, grid.h.values, lev11, colors='deeppink')
ax10[2].set_title('1 m Above Seafloor Current Magnitude (cm/s)', fontsize=fontsize)
ax10[2].set_xlabel('X (km)', fontsize=fontsize)
ax10[2].set_ylabel('Y (km)', fontsize=fontsize)

cbar14 = plt.colorbar(cs24, ax=[ax10[0], ax10[1], ax10[2]], extend='both').set_label('Current Magnitude (cm/s)', fontsize=fontsize)

fig10.text(0.5, 0.92, 'Standard Deviation Current Magnitude Over Time', ha="center", va="center", rotation=0, fontsize=fontsize+5)


