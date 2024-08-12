#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:38:35 2024

@author: Meredith

This script extracts the spectra for kelp presence/absence, and vegetation
presence/absence via depth bins.
The presence/absence is determined via SACFOR percent cover.
Catagories were divided and a threshold cover was determined for presence/absence
"""
import geopandas as gpd
import pandas as pd
import rasterio as rio
import numpy as np
import matplotlib.pyplot as plt


# %%
# sacfor file with vegetation info
sacfor_file = 'data/dropcam/gooseberries_sacfor_v2_subset_v2.shp'
sacfor = gpd.read_file(sacfor_file)

# add column for non-kelp vegetation
sacfor['veg_corr'] = sacfor['veg #'] - sacfor['kel #']
sacfor['veg_corr_PA'] = sacfor['veg_corr'].isin([0, 1, 2, 3]).astype(int)

# wl file
wl_file = 'data/headwall/wv_v2.csv'
wl = pd.read_csv(wl_file,names=['count','wl'])

# image file
img_file = 'data/headwall/nano_Rrs_deglint_smooth_mosaic_mean_v2_NS.tif'

# %%

# Define a function to extract spectra for each point
def extract_spectra(point, src):
    row, col = src.index(point.x, point.y)
    spectra = src.read(window=((row, row + 1), (col, col + 1))).flatten()
    return spectra

# define depth bins
# Calculate depth bins in 1 m increments
min_depth = np.floor(sacfor['Elev(m)'].min())
max_depth = np.ceil(sacfor['Elev(m)'].max())
depth_bins = np.arange(min_depth, max_depth + 1, 1)

# Create an empty DataFrame to store the mean spectra for each depth bin
mean_kelp_pres_df = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])
mean_kelp_abs_df = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])

# Create an empty DataFrame to store the mean spectra for each depth bin
mean_veg_pres_df = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])
mean_veg_abs_df = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])

mean_kelp_depth = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])
mean_veg_depth = pd.DataFrame(columns=[f'Depth_{depth_bin}' for depth_bin in depth_bins[:-1]])

# Read the hyperspectral image
with rio.open(img_file) as src:
    arr = src.read(masked=True)
    transform = src.transform


    # Iterate over each point in the sacfor GeoDataFrame
    for i in range(len(depth_bins) - 1):
        
        # Select rows within the current depth bin
        sacfor_in_depth_bin = sacfor[(sacfor['Elev(m)'] >= depth_bins[i]) & (sacfor['Elev(m)'] < depth_bins[i + 1])]
        
        # calculate mean depths
        
        
        # create emtpy arrays to gather all the spectra for each class
        kelp_present = np.full((len(wl),len(sacfor_in_depth_bin)),np.nan)
        kelp_absent = np.full((len(wl),len(sacfor_in_depth_bin)),np.nan)
        veg_present = np.full((len(wl),len(sacfor_in_depth_bin)),np.nan)
        veg_absent = np.full((len(wl),len(sacfor_in_depth_bin)),np.nan)
        
        for l, (idx, point) in enumerate(sacfor_in_depth_bin.iterrows()):
            # extract spectra
            spectra = extract_spectra(point.geometry, src)
            
            # check if kelp pres/abs
            if sacfor.iloc[idx]['kelp_pres/'] == 1:
                # add spectra to temporary array before taking mean
                kelp_present[:,l] = spectra
            else:
                kelp_absent[:,l] = spectra
                
            # check if veg pres/abs
            if sacfor.iloc[idx]['veg_corr_PA'] == 1:
                # add spectra to temporary array before taking mean
                veg_present[:,l] = spectra
            else:
                veg_absent[:,l] = spectra
        
        # plt.plot(wl['wl'],kelp_present)
        # plt.title(f'Kelp in Depth = {depth_bins[i]} to {depth_bins[i]+1}m')
        # plt.show()
        
        # Calculate the mean spectra for each depth
        kelp_present_mean = np.nanmean(kelp_present,axis=1)
        kelp_absent_mean = np.nanmean(kelp_absent, axis=1)
        veg_present_mean = np.nanmean(veg_present, axis=1)
        veg_absent_mean = np.nanmean(veg_absent, axis=1)
            
        # add means to appropriate depth bin dataframe
        # Add the mean spectrum to the DataFrame
        mean_kelp_pres_df[f'Depth_{depth_bins[i]}'] = kelp_present_mean
        mean_kelp_abs_df[f'Depth_{depth_bins[i]}'] = kelp_absent_mean
        
        mean_veg_pres_df[f'Depth_{depth_bins[i]}'] = veg_present_mean
        mean_veg_abs_df[f'Depth_{depth_bins[i]}'] = veg_absent_mean
 
# %%
# Plot the mean spectra
plt.figure(figsize=(25, 4))

for i in range(len(depth_bins) - 1): 
    plt.subplot(1,5,i+1)
    plt.plot(wl['wl'], mean_kelp_pres_df[f'Depth_{depth_bins[i]}'],label='Kelp Present')
    plt.plot(wl['wl'], mean_veg_pres_df[f'Depth_{depth_bins[i]}'],label = 'Vegetation Present')
    
    # plt.plot(wl['wl'], mean_kelp_abs_df[f'Depth_{depth_bins[i]}'],label = 'Kelp Absent')
    # plt.plot(wl['wl'], mean_veg_abs_df[f'Depth_{depth_bins[i]}'],label = 'Vegetation Absent')
    
    
        
    plt.title(f'Depth = {depth_bins[i]} to {depth_bins[i]+1}m',fontsize=16)
    plt.ylabel('R$_{rs}$ (sr$^{-1}$)',fontsize=14)
    plt.xlabel('Wavelength',fontsize=14)
    plt.ylim([0,0.008])
    if i == 0:
        plt.legend()
    else:
        continue
    
    
plt.tight_layout()
plt.savefig('output/figs/kelp_veg_present_spectra_with_depth.png')
plt.show()
          
# %% Plot the mean spectra in reverse order

# 
plt.figure(figsize=(20, 4))  # Adjust the figure size for 1 row, 5 column layout

for i in range(len(depth_bins) - 2, -1, -1):  # Reverse the range to go from last to first
    # Adjust the subplot parameters to reflect the reversed order
    plt.subplot(1, 5, len(depth_bins) - 1 - i)  
    plt.plot(wl['wl'], mean_kelp_pres_df[f'Depth_{depth_bins[i]}'], label='Kelp Present')
    plt.plot(wl['wl'], mean_veg_pres_df[f'Depth_{depth_bins[i]}'], label='Vegetation Present')
    
    # Uncomment these lines to include absent spectra
    # plt.plot(wl['wl'], mean_kelp_abs_df[f'Depth_{depth_bins[i]}'], label='Kelp Absent')
    # plt.plot(wl['wl'], mean_veg_abs_df[f'Depth_{depth_bins[i]}'], label='Vegetation Absent')
    
    plt.title(f'Depth = {depth_bins[i]} to {depth_bins[i]+1}m', fontsize=16)
    plt.ylabel('R$_{rs}$ (sr$^{-1}$)', fontsize=14)
    plt.xlabel('Wavelength', fontsize=14)
    plt.ylim([0, 0.008])
    if i == len(depth_bins) - 2:
        plt.legend()
    else:
        continue


plt.tight_layout()
plt.savefig('output/figs/kelp_veg_present_spectra_with_depth.png')
plt.show()
# %%
# Plot the mean spectra
# plt.figure(figsize=(10, 10))

for i in range(len(depth_bins) - 1): 
    plt.subplot(1,1,1)
    plt.plot(wl['wl'], mean_kelp_pres_df[f'Depth_{depth_bins[i]}'],label='Kelp Present')
    plt.plot(wl['wl'], mean_veg_pres_df[f'Depth_{depth_bins[i]}'],label = 'Vegetation Present')
    
    # plt.plot(wl['wl'], mean_kelp_abs_df[f'Depth_{depth_bins[i]}'],label = 'Kelp Absent')
    # plt.plot(wl['wl'], mean_veg_abs_df[f'Depth_{depth_bins[i]}'],label = 'Vegetation Absent')
        
    plt.title(f'Depth = {depth_bins[i]} to {depth_bins[i]+1}m',fontsize=22)
    plt.ylabel('R$_{rs}$ (sr$^{-1}$)',fontsize=18)
    plt.xlabel('Wavelength',fontsize=18)
    plt.ylim([0,0.012])
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(f'output/figs/kelp_veg_present_spectra_{depth_bins[i]}_to_{depth_bins[i]+1}m.png',dpi=300)
    plt.show()

# %%
# plot the depth binned kelp present and veg present 


