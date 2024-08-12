#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:45:37 2024

@author: Meredith

This script compiles the spectra from the headwall imagery to bin into 1 m depths
and plots them to compare.


"""
import rasterio as rio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# imagery file
img_file = 'data/headwall/nano_Rrs_deglint_smooth_mosaic_mean_v2_NS.tif'

# bathymetry file
bathy_file = 'data/lidar/2022_ngs_topobathy_chatham_ma_Job964557_reproj_NS.tif'

# wl file
wl_file = 'data/headwall/wv_v2.csv'
wl = pd.read_csv(wl_file,names=['count','wl'])


'''
Extracting depth specific spectra from imagery:
    
open the tif files 

depth bins in 1 m increments from the bathymetry file 

loop through depth bins from bathymetry file
gather all the spectra in that bin
plot  all the spectra for each depth bin.

take the mean of all the spectra in each depth bin and save to a dataframe

plot a summary plot of all the means from depth bins
'''

# Open the bathymetry file
with rio.open(bathy_file) as src:
    bathymetry = src.read(1, masked=True)  # Read bathymetry data, masking nodata values
    bathymetry_profile = src.profile  # Get the metadata profile

# Calculate depth bins in 1 m increments
min_depth = np.floor(np.nanmin(bathymetry))
max_depth = np.ceil(np.nanmax(bathymetry))
depth_bins = np.arange(min_depth, max_depth + 1, 1)

# Open the imagery file
with rio.open(img_file) as src:
    spectral_data = src.read(masked=True)  # Mask nodata values
    spatial_metadata = src.profile

# Create an empty DataFrame to store the mean spectra for each depth bin
mean_spectra_df = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])

# Loop through depth bins
for i in range(len(depth_bins) - 1):
    # Select pixels within the current depth bin
    pixels_in_depth_bin = spectral_data[:, (bathymetry >= depth_bins[i]) & (bathymetry < depth_bins[i + 1])]
   
    # Calculate the mean spectrum for the current depth bin
    mean_spectrum = np.mean(pixels_in_depth_bin, axis=1)
    
    # Add the mean spectrum to the DataFrame
    mean_spectra_df[f'Depth_{depth_bins[i]}'] = mean_spectrum
    
    # Plot the points on the image
    # Get the XY coordinates corresponding to the selected pixels
    rows, cols = np.where((bathymetry >= depth_bins[i]) & (bathymetry < depth_bins[i + 1]))
    
    # Plot the image
    plt.imshow(spectral_data[50,:,:], cmap='viridis')
    
    # Plot the XY coordinates on the image
    plt.scatter(cols, rows, color='red', s=5)  # Use cols as X and rows as Y
    plot_title = f'Depth Bin = {depth_bins[i]} to {depth_bins[i+1]}m'  
    plt.title(plot_title) 
    plt.show()
    
    # Plot all the spectra separately
    plt.plot(wl['wl'],pixels_in_depth_bin)
    plot_title = f'Depth Bin = {depth_bins[i]} to {depth_bins[i+1]}m'
    plt.title(plot_title)
    plt.show()
    
# Plot the mean spectra for each depth bin
for column in mean_spectra_df.columns:
    plt.plot(wl['wl'], mean_spectra_df[column], label=column)
    plt.vlines(685,0,0.01,color='k')
    plt.vlines(740,0,0.01,color='k')

plt.xlabel('Wavelength')
plt.ylabel('Reflectance')
plt.ylim([0,0.01])
plt.title('Mean Spectra for Each Depth Bin')
plt.legend()
plt.show()



plt.plot(wl['wl'],mean_spectra_df['Depth_-4.0'],label='-4 to -3m')
plt.plot(wl['wl'],mean_spectra_df['Depth_-5.0'],label='-5 to -4m')
plt.plot(wl['wl'],mean_spectra_df['Depth_-6.0'],label='-6 to -5m')
plt.plot(wl['wl'],mean_spectra_df['Depth_-7.0'],label='-7 to -6m')
plt.plot(wl['wl'],mean_spectra_df['Depth_-8.0'],label='-8 to -7m')
plt.ylim([0,0.01])
plt.vlines(680,0,0.01,color='k')
plt.vlines(740,0,0.01,color='k')
plt.legend()
plt.show()
