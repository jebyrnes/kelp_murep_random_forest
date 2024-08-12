#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 14:11:39 2024

@author: Meredith

This script plots plots the 
"""

import numpy as np
import rasterio as rio
import geopandas as gpd
import pandas as pd
import os
from shapely.geometry import box
import matplotlib.pyplot as plt

os.chdir('/Users/Meredith/Library/CloudStorage/GoogleDrive-meredith.mcpherson@umb.edu/.shortcut-targets-by-id/1PJamAGTknM_kuVq8M53Xo4cAlWzDOk6I/kelp_murep/')

# hyperspectral tif file
img_file = 'rs_imagery_processing/Headwall/output/mmcpherson/mosaic/native_cl_no_1m.tif'

# bathymetry tif file
bathy_reproject_file = 'rs_imagery_processing/classification/random_forest/data/headwall/bathy_reprojected_headwall_norm_1m.tif'

# load the sacfor data shape file 
sacfor = gpd.read_file('sidescan_surveys/gooseberries_sacfor_raster_samples/gooseberries_sacfor_raster_samples_corr.shp')
sacfor['SACFOR_veg'] = sacfor['SACFOR_veg'].replace('s', 'S')

# Define the bounding box to crop to (xmin, ymin, xmax, ymax)
gb_subset_bbox = [258025, 919500, 258150, 919650]

# Set depth threshold
z_thresh = -5

# Crop and mask everything

# Load the hyperspectral data 
with rio.open(img_file) as src:
    arr = src.read()
    
    # HEADWALL masking to nan values
    arr[arr == -9999] = np.nan
    arr[arr < 0] = np.nan
    
    meta = src.meta
    
    # convert the sacfor dataframe to the hyperspectral crs
    sacfor = sacfor.to_crs(src.crs)
    
    # Open the reprojected bathymetry file in read mode and read the data into bathy_arr
    with rio.open(bathy_reproject_file) as reproj_bathy:
        bathy_arr = reproj_bathy.read(1)

        bathy_arr[bathy_arr==-3.4028234663852886e+38] = np.nan

        # ------ HEADWALL BBOX
        # get the bounding corners and create window for the headwall src
        headwall_bbox = src.bounds
        headwall_bbox = [headwall_bbox.left, headwall_bbox.bottom, headwall_bbox.right, headwall_bbox.top]
           
        headwall_window = rio.windows.from_bounds(*headwall_bbox, transform=meta['transform'])

        # crop the bathymetry data to the headwall bbox
        bathy_headwall = reproj_bathy.read(window=headwall_window)[0]
        bathy_headwall[bathy_headwall == -3.4028234663852886e+38] = np.nan

        # crop sacfor to headwall bbox
        sacfor_headwall = sacfor[sacfor.within(box(*headwall_bbox))]

        # Mask the headwall array where the values are greater than -4 or nan
        nanmask = np.isnan(bathy_arr)
        arr_headwall_nanmasked = np.ma.masked_array(arr, mask=np.tile(nanmask, (arr.shape[0], 1, 1)))

        zmask = bathy_arr < z_thresh
        arr_headwall_zmasked = np.ma.masked_array(arr_headwall_nanmasked, mask=np.tile(zmask, (arr.shape[0], 1, 1)))
        
        # ----- GOOSEBERRY SUBSET BBOX
        # create window for subset
        subset_window = rio.windows.from_bounds(*gb_subset_bbox, 
                                                transform=meta['transform'])

        # crop the headwall array to the subset bbox
        arr_subset = src.read(window=subset_window)
        
        # crop the bathymetry data to the subset bbox
        bathy_subset = reproj_bathy.read(window=subset_window)[0]
        bathy_subset[bathy_subset == -3.4028234663852886e+38] = np.nan
        
        
        # Filter sacfor points within the subset bounding box
        sacfor_subset = sacfor[sacfor.within(box(*gb_subset_bbox))]

        # sort out deep elevation in SACFOR
        sacfor_subset_filt = sacfor_subset[sacfor_subset['Elev(m)']>=z_thresh]
        
        

        # Mask the subset array where the values are greater than -4 or nan

        nanmask = np.isnan(bathy_subset)
        arr_subset_nanmasked = np.ma.masked_array(arr_subset, 
                                                  mask=np.tile(nanmask, 
                                                               (arr_subset.shape[0], 1, 1)))
        
        zmask = bathy_subset < z_thresh
        arr_subset_zmasked = np.ma.masked_array(arr_subset_nanmasked, 
                                                mask=np.tile(zmask, 
                                                             (arr_subset.shape[0], 1, 1)))

        #bathy_subset_zmasked = np.ma.masked_array(bathy_subset, mask=np.tile(zmask,(bathy_subset.shape,1)))


# ---- Bin the drop camera validation data into depth bins

bins = [-6,-5,-4,-3,-2,-1,0]
labels = [5,4,3,2,1,0]

# Cut the data into bins and count the occurrences in each bin
# subset area all depths
sacfor_subset['Elev(m)_label'] = pd.cut(sacfor_subset['Elev(m)'], bins=bins, 
                                        labels=labels)
bin_counts = sacfor_subset['Elev(m)_label'].value_counts(sort=False)

# subset area depth filtered
sacfor_subset_filt['Elev(m)_label'] = pd.cut(sacfor_subset_filt['Elev(m)'], bins=bins, 
                                        labels=labels)
bin_counts_filt = sacfor_subset_filt['Elev(m)_label'].value_counts(sort=False)

# entire gooseberry area
sacfor['Elev(m)_label'] = pd.cut(sacfor['Elev(m)'], bins=bins, 
                                        labels=labels)
bin_counts_all = sacfor['Elev(m)_label'].value_counts(sort=False)



# Assuming bin_count_all and bin_count_subset are your Series
# Replace 'labels' with your custom labels

labels = ['5-6','4-5','3-4','2-3','1-2','0-1']
bin_counts.index = labels
bin_counts_all.index = labels
bin_counts_filt.index = labels

# Create a DataFrame with the two bin count variables
df_plot = pd.DataFrame({'bin_count_all': bin_counts_all, 
                        'bin_count_subset': bin_counts}, 
                       index=labels)

# Plotting side by side bar plot
ax = df_plot.plot.bar(rot=0, width=0.8, color=['skyblue', 'orange'])
ax.legend([f'Entire GB (n = {len(sacfor)})',
           f'GB subset (all depths; n = {len(sacfor_subset)})'])
# Customize the plot if needed
plt.title('Bin Counts Comparison')
plt.xlabel('Elevation Bins (m)')
plt.ylabel('Count')

# Show the plot
plt.show()


