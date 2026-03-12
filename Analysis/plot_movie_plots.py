############################# Plots for Movies! ###################################
# The purpose of this script is to ake various plots for various movies
# and run them in parallel and all that jazz so that we can make ool movies of 
# things.
#
# Notes:
# - This works best in xroms? xesmf_env? conda environment(s)
# 
###################################################################################


# Load in the packages 
import numpy as np
import glob
import gc
import cmocean.cm as cmo
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import LogFormatterExponent
import warnings
import xarray as xr
import xesmf as xe
import os
#import xroms
from joblib import Parallel, delayed
from matplotlib import ticker
import cartopy.mpl.ticker as cticker
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
crs = ccrs.PlateCarree()
#Cartopy
land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                edgecolor='face',
                                facecolor=cfeature.COLORS['land'])


# Load in the model output(s)
model_output_std = xr.open_mfdataset('/pl/active/moriarty_lab/BriannaU/Paper1/2020_version/Output/dbsed0003/ocean_his_beaufort_shelf_2020_dbsed0003_*.nc')

# Load in the model grid
grid = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Include/KakAKgrd_shelf_big010_smooth006.nc')
#grid = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Grids/KakAKgrd_shelf_big010_smooth006.nc') # UPDATE PATH

# Load in the wave forcing file 
wave_frc = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine_2020/Include/wave_forcing_file_kaktovik_shelf_ww3_2020_data002.nc')

# Load in nan masks 
mask_rho_nan = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_ones_nans.nc') # UPDATE PATH
mask_rho_zeros = xr.open_dataset('/projects/brun1463/ROMS/Kakak3_Alpine/Scripts_2/Analysis/Nudge_masks/nudge_mask_rho_zeros_ones.nc')
# mask_rho_nan = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_ones_nans.nc')
# mask_rho_zeros = xr.open_dataset('/Users/brun1463/Desktop/Research_Lab/Kaktovik_Alaska_2019/Code/Nudge_masks/nudge_mask_rho_zeros_ones.nc')

# ------------------------------------------------------------------------------
# ----------------------------- Define Functions -------------------------------
# ------------------------------------------------------------------------------

# Make a function to calculate the sediment flux in the 
# u direction for all sediment classes combined
def calc_u_east_depth_int_ssc_flux_allsed(model_output): #, regridder_u2rho):
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
    
    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
    
    # To collapse to horizontal, multiply each layer by its
    # thickness
    # Calculate the time-varying thickness of the cells
    dz = abs(model_output.z_w[:-1,:,:].values - model_output.z_w[1:,:,:].values)
    
    # Pull out the u velocities at all times, depths, spaces
    u_tmp = model_output.u_eastward
    
    # Interpolate them onto rho points 
    u_tmp_rho = u_tmp
    #u_tmp_rho = regridder_u2rho(u_tmp)
    
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
    depth_avg_ssc_allsed = depth_int_ssc_allsed/model_output.bath[:,:].values
    
    # Return the depth-integrated u flux for all sediment classes
    return(depth_int_ssc_flux_u_allsed, depth_avg_ssc_allsed)


# Make a function to calculate the sediment flux in the 
# v direction for all sediment classes combined
def calc_v_north_depth_int_ssc_flux_allsed(model_output): #, regridder_v2rho):
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
    
    # Add all the sediment classes together
    ssc_allsed_tmp = model_output.mud_01 + model_output.mud_02 + model_output.sand_01 + model_output.sand_02 + model_output.sand_03
    
    # To collapse to horizontal, multiply each layer by its
    # thickness
    # Calculate the time-varying thickness of the cells
    dz = abs(model_output.z_w[:-1,:,:].values - model_output.z_w[1:,:,:].values)
    
    # Pull out the v velocities at all times, depths, spaces
    v_tmp = model_output.v_northward
    
    # Interpolate them onto rho points 
    v_tmp_rho = v_tmp
    #v_tmp_rho = regridder_v2rho(v_tmp)
    
    # Pull out the thickness of the cell in the x direction 
    #dx = 1.0/model_output.pm
    
    # Use all of this to calculate depth-integrated sediment flux
    # First just calculate flux at all times over all space
    ssc_flux_allsed = ((ssc_allsed_tmp*v_tmp_rho)*(dz))
    
    # Then depth-integrated by summing over depth and dividing by dx
    depth_int_ssc_flux_v_allsed = (ssc_flux_allsed.sum(dim='s_rho'))
    
    # Return the depth-integrated v flux for all sediment classes
    return(depth_int_ssc_flux_v_allsed)


# ------------------------------------------------------------------------------
# ------------------------------- Movie Prep -----------------------------------
# ------------------------------------------------------------------------------

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


river_marker_colors = ['#FC440F', '#F5ED00', '#5EF38C', '#26532B', '#F43ECF', '#9C00A8',
                       '#0115F5', '#00A6A6', '#AB64EB', '#D44179', '#08E0E3', '#B27009', '#EA8D40']


# ------------------------------------------------------------------------------
# --------------------------- Make Movie(s)! -----------------------------------
# ------------------------------------------------------------------------------

# Pull out shapes
eta_rho_len = len(grid.eta_rho)
xi_rho_len = len(grid.xi_rho)

# Make land mask
temp_mask = grid.mask_rho.values
temp_mask = np.where(temp_mask==0, np.nan, temp_mask)

# Set the speeds
speed_min = 0 # m/s
speed_max = 0.2 # m/s
cmap_spd = cmo.dense

# Set the SSC stuff
ssc_min = 0 # kg/m3
ssc_max = 0.01 # kg/m3
cmap_ssc = cmo.turbid


# Make the function to plot all the plots
# Now make a movie of this

def plot_bstrcwmax_cur_seds(ds, time_idx):

    #print('time_idx: ', time_idx)
    i = time_idx // 1

    # Define the projection
    crs = ccrs.PlateCarree()

    # Calculate the different things 
    # Current magnitude
    cur_mag = np.sqrt(((ds.ubar_eastward[:,:])**2) + ((ds.vbar_northward[:,:])**2))
    # Total depth-averaged SSC
    # Depth-integrated sediment flux
    u_east_depth_int_ssflux, depth_avg_ssc = calc_u_east_depth_int_ssc_flux_allsed(ds)
    v_north_depth_int_ssflux = calc_v_north_depth_int_ssc_flux_allsed(ds)

    # Set the scales for the currents
    scale_current = 1
    scale_sed_flux = 1

    # Set the number of quivers
    n_quiv = 25

    # Make the figure; add things to make it the right ratios and all that
    nrows = 3; ncols = 1
    fig, ax = plt.subplots(nrows, ncols, figsize=(12, 8), # ((9.75/3)*1.9, 5); (4.5, 5) 
                        #sharex=True, sharey=True, dpi=300,
                        dpi=300, subplot_kw={'projection': crs})

    # Total bed shear stress
    m1 = ax[0].pcolormesh(grid.lon_rho.values, grid.lat_rho.values, ds.bstrcwmax[:,:]*mask_rho_nan.nudge_mask_rho_nan,
                            cmap=cmo.amp, vmin=0, vmax=0.2)
    ax[0].set_title('Bed Shear Stress')
    # Add colorbar
    cbar_ax1 = fig.add_axes([0.30, 0.71, 0.13, 0.015])
    fig.colorbar(m1, orientation='horizontal', label='Bed Shear Stress (N/m\u00b2)', 
                cax=cbar_ax1, extend='max')

    # Depth-averaged current magnitude
    m2 = ax[1].pcolormesh(grid.lon_rho.values, grid.lat_rho.values, cur_mag*mask_rho_nan.nudge_mask_rho_nan, 
                        cmap=cmap_spd, vmin = speed_min, vmax=speed_max)
    ax[1].set_title('Depth-Averaged Current Magnitude')
    # Add colorbar
    cbar_ax2 = fig.add_axes([0.30, 0.445, 0.13, 0.015])
    fig.colorbar(m2, orientation='horizontal', label='Current Mag. (m/s)', 
                cax=cbar_ax2, extend='max')
    #fig.colorbar(m1, ax=ax[0,0], label='[C]', extend='both')
    # Add quivers of the normalized current current 
    q1 = ax[1].quiver(grid.lon_rho[::n_quiv,::n_quiv].values, grid.lat_rho[::n_quiv,::n_quiv].values, 
                   model_output_std.ubar_eastward[time_idx,::n_quiv,::n_quiv]*mask_rho_nan.nudge_mask_rho_nan[::n_quiv,::n_quiv], 
                   model_output_std.vbar_northward[time_idx,::n_quiv,::n_quiv]*mask_rho_nan.nudge_mask_rho_nan[::n_quiv,::n_quiv], 
                   transform=ccrs.PlateCarree(), 
                   color='deeppink', width=0.005, scale=scale_current,
                   angles='xy', scale_units='xy')
    ax[1].quiverkey(q1, 0.85, 0.75, U=0.15, label='0.15 m/s', fontproperties={'size':12})

    # SSC and sediment flux
    m3 = ax[2].pcolormesh(grid.lon_rho.values, grid.lat_rho.values,depth_avg_ssc[:,:]*mask_rho_nan.nudge_mask_rho_nan*1000,
                            cmap=cmap_ssc, vmin=0, vmax=10)
    ax[2].set_title('Depth-Averaged SSC (mg/L)')
    # Add colorbar
    cbar_ax3 = fig.add_axes([0.30, 0.18, 0.13, 0.015])
    fig.colorbar(m3, orientation='horizontal', label='SSC (mg/L)', 
                cax=cbar_ax3, extend='max')
    # Add quivers for sediment flux
    q2 = ax[2].quiver(grid.lon_rho[::n_quiv,::n_quiv].values, grid.lat_rho[::n_quiv,::n_quiv].values, 
                   u_east_depth_int_ssflux[::n_quiv,::n_quiv]*mask_rho_nan.nudge_mask_rho_nan[::n_quiv,::n_quiv], 
                   v_north_depth_int_ssflux[::n_quiv,::n_quiv]*mask_rho_nan.nudge_mask_rho_nan[::n_quiv,::n_quiv], 
                   transform=ccrs.PlateCarree(),
                   color='teal', width=0.005, scale=scale_sed_flux,
                   angles='xy', scale_units='xy')
    ax[2].quiverkey(q2, 0.85, 0.75, U=0.001, label='0.001 kg/m\u00b2s', fontproperties={'size':12})


    # Loop through axes to add features
    for r in range(3):
        # Plot bathymetry 
        #c1 = ax[r,j].contour(grid.lon_rho, grid.lat_rho, grid.h*mask_rho_nan.nudge_mask_rho_nan, lev_bathy, colors='darkslategrey', linewidth=2)
        #ax[r,j].clabel(c1, inline=True, fontsize=10, colors='white')
        
        # Set extent and map features
        #ax[r].set_extent([-153.5,-141.4,69.5, 71.5],ccrs.PlateCarree())
        ax[r].set_extent([-153,-142.2,69.75, 71.5],ccrs.PlateCarree())
        # ax.set_aspect(lat_rad)
        ax[r].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                                    facecolor='darkgray'), linewidth=.7)
        ax[r].coastlines(resolution='10m', linewidth=.7)

        # Set tick locations
        xticks = [-152, -150, -148, -146, -144, -142]
        yticks = [70, 70.5, 71, 71.5]

        ax[r].set_yticks(yticks, crs=ccrs.PlateCarree())
        ax[r].yaxis.set_major_formatter(cticker.LatitudeFormatter())
        ax[r].set_xticks(xticks, crs=ccrs.PlateCarree())
        ax[r].xaxis.set_major_formatter(cticker.LongitudeFormatter())

        # ---- Remove labels for first two panels ----
        if r in [0,1]:
            ax[r].tick_params(labelbottom=False)

        # Fix the aspect ratio
        lon_min, lon_max = -153, -142.2
        lat_min, lat_max = 69.75, 71.5
        mid_lat = (lat_min + lat_max) / 2
        aspect = (lon_max - lon_min) * np.cos(np.deg2rad(mid_lat)) / (lat_max - lat_min)
        ax[r].set_aspect(aspect)

        # Add river markers
        # Go from West to East
        # Kalikpik River
        eta_kal_idx = 23 #22
        xi_kal_idx = 87
        s1 = ax[r].scatter(grid.lon_rho[eta_kal_idx, xi_kal_idx].values, grid.lat_rho[eta_kal_idx, xi_kal_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[0], label='Kalikpik')

        # Fish Creek
        eta_fis_idx = 20
        xi_fis_idx = 117 #116
        s2 = ax[r].scatter(grid.lon_rho[eta_fis_idx, xi_fis_idx].values, grid.lat_rho[eta_fis_idx, xi_fis_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[1], label='Fish Creek')

        # Colville River
        eta_col_idx = 39
        xi_col_idx = 166 #166
        s3 = ax[r].scatter(grid.lon_rho[eta_col_idx, xi_col_idx].values, grid.lat_rho[eta_col_idx, xi_col_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[2], label='Colville')

        # Sakonowyak River
        eta_sak_idx = 46 #45
        xi_sak_idx = 234
        #s4 = ax11.scatter(x_rho_flat[xi_sak_idx]/1000, y_rho_flat[eta_sak_idx]/1000, 
        #           marker='.', s=500, linewidth=1, edgecolors='black', color='green', label='Sakonowyak')

        # Kuparik
        # Kukpuk - Change this to be labeled as the Kuparuk since it  is actually the main channel 
        # of the Kuparuk River
        eta_kuk_idx = 41 #40
        xi_kuk_idx = 239
        eta_kup_idx = 41 #40
        xi_kup_idx = 242
        s5 = ax[r].scatter(grid.lon_rho[eta_kuk_idx, xi_kuk_idx].values, grid.lat_rho[eta_kuk_idx, xi_kuk_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[4], label='Kuparuk')


        # Putuligayuk River
        eta_put_idx = 28 #27
        xi_put_idx = 264
        #s8 = ax11.scatter(x_rho_flat[xi_put_idx]/1000, y_rho_flat[eta_put_idx]/1000, 
        #           marker='1', s=500, linewidth=1, edgecolors='black', color='dodgerblue', label='Putuligayuk')

        # Sagavanirktok River
        eta_sag_idx = 37 #36
        xi_sag_idx = 279
        s9 = ax[r].scatter(grid.lon_rho[eta_sag_idx, xi_sag_idx].values, grid.lat_rho[eta_sag_idx, xi_sag_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[6], label='Sagavanirktok')

        # Canning River
        # Staines River
        eta_sta_idx = 27 #26
        xi_sta_idx = 393
        s10 = ax[r].scatter(grid.lon_rho[eta_sta_idx, xi_sta_idx].values, grid.lat_rho[eta_sta_idx, xi_sta_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[7], label='Staines')

        # Canning River
        eta_can_idx = 20 #19
        xi_can_idx = 416
        s11 = ax[r].scatter(grid.lon_rho[eta_can_idx, xi_can_idx].values, grid.lat_rho[eta_can_idx, xi_can_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[8], label='Canning')

        # Katakturuk River
        eta_kat_idx = 9 #8
        xi_kat_idx = 447
        s12 = ax[r].scatter(grid.lon_rho[eta_kat_idx, xi_kat_idx].values, grid.lat_rho[eta_kat_idx, xi_kat_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[9], label='Katakturuk')

        # Hulahula River
        eta_hul_idx = 40
        xi_hul_idx = 489
        s13 = ax[r].scatter(grid.lon_rho[eta_hul_idx, xi_hul_idx].values, grid.lat_rho[eta_hul_idx, xi_hul_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[10], label='Hulahula')

        # Jago River
        eta_jag_idx = 62 #61
        xi_jag_idx = 528
        s14 = ax[r].scatter(grid.lon_rho[eta_jag_idx, xi_jag_idx].values, grid.lat_rho[eta_jag_idx, xi_jag_idx].values, 
                    marker='.', s=500, linewidth=1, edgecolors='black', color=river_marker_colors[11], label='Jago')

        # Siksik River
        eta_sik_idx = 46
        xi_sik_idx = 574 #573
        #s15 = ax11.scatter(x_rho_flat[xi_sik_idx]/1000, y_rho_flat[eta_sik_idx]/1000, 
        #           marker='.', s=300, linewidth=4, color='deeppink', label='Siksik')

    # Add time stamp as title
    fig.suptitle(str(ds.ocean_time.values)[:10], y=0.94)

    # Add legend for rivers
    ax[2].legend(ncol=3, bbox_to_anchor=(0.90,-0.17))
   
    plt.savefig(
        f'/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Movies/Bed_stress_cur_seds/Plots03/bstrcwmax_depth_avg_cur_depth_avg_ssc_int_flux_std_{i}.png',
        dpi=200, bbox_inches='tight'
    )
    plt.close(fig)  # Close to avoid memory leaks
    gc.collect()


# Call the function above in a loop to make the plots 
# for the movie
print('starting movie', flush=True)
for time_idx in range(0, len(model_output_std.ocean_time)):
    # Load the timestep of data
    ds = model_output_std.isel(ocean_time=time_idx).load()
    #print(f"Processing time index: {time_idx}")
    plot_bstrcwmax_cur_seds(ds, time_idx)

print('dont with movie', flush=True)






