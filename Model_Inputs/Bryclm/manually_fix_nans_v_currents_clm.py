# Manually fix nans in land mask in v_current_clm_002.nc

# Import Functions 
import numpy as np
import xarray as xr

# Load in the data
v_clm = xr.open_dataset('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/Attempt001/v_currents_clm_002.nc')

# Make a copy to edit
v_clm_cp = v_clm.copy()

# Pull out the length of time 
time_len = len(v_clm.v3d_time)

# Use Xarray inteprolate_na one just the level(s) that have nans
# Loop through time
for t in range(time_len):
    # Interpolate over nans 
    # Level 15
    #v_clm_cp.v[t,15,:,:].interpolate_na(dim='xi_v', method='linear', max_gap=None)
    # Level 19
    #v_clm_cp.v[t,19,:,:].interpolate_na(dim='eta_v', method='linear')
    
    # Replace nans with 0
    v_clm_cp.v[t,19,:,:] = v_clm_cp.v[t,19,:,:].fillna(0.0)
    

# Save this to a netcdf
v_clm_cp.to_netcdf('/scratch/alpine/brun1463/ROMS_scratch/Kakak3_Alpine_2020_scratch/Final_bryclm_conds/Attempt001/v_currents_clm_002_fix02.nc')
