#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:37:22 2024

@author: Meredith
"""

import rasterio as rio
import numpy as np
import matplotlib.pyplot as plt

import rioxarray as riox

from pysal.lib import weights
from esda.getisord import G_Local
from esda.getisord import G

# %%


fn = 'data/headwall/nano_Rrs_deglint_smooth_mosaic_mean_v2_NS.tif'

with rio.open(fn) as src:
    arr = src.read()
    mask = arr == -9999.
    
    # Get the affine transformation matrix
    transform = src.transform
    
    # get the mask of arr to mask the coordinates

    # Get the grid of pixel centers in the raster
    ny, nx = arr.shape[1:]
    x, y = np.meshgrid(np.arange(0.5, nx), np.arange(0.5, ny))
    # x = np.ma.masked_array(x, mask=arr.mask[0,:,:])
    # y = np.ma.masked_array(y, mask=arr.mask[0,:,:])

    # Convert pixel coordinates to world coordinates
    xt, yt = transform * (x.flatten(), y.flatten())
    xt = xt.astype('float64')
    yt = yt.astype('float64')

    # Create a flattened array of coordinates
    coordinates = [(xt[i], yt[i]) for i in range(len(xt))]

# using rioxarray to open the data
data_riox = riox.open_rasterio(fn, masked=True)

# %% Run Getis Ord Statistic

# everything must be in float64


# Define a function to calculate Getis Ord statistic
def calculate_getis_ord_statistic(spatial_data, raster_data, threshold):
    
    # create empty arrays to fill 
    lg = np.zeros((raster_data.shape))
    lg_stand = np.zeros((raster_data.shape))
    lg_pvals = np.zeros((raster_data.shape))
    
    
    # # for i in range(len(raster_data)):
    for i in range(1):
    
        # flatten raster data which is already masked when opened
        raster_flat = raster_data[0,:,:].flatten().reshape(len(spatial_data),1)
        raster_flat = raster_flat.astype('float64')
        
        
        # Generate spatial weights using distance threshold
        w = weights.distance.DistanceBand.from_array(np.array(spatial_data), threshold=threshold)
        # w = weights.distance.DistanceBand(spatial_data, threshold=threshold)
    
        # Ensure the weights are symmetric
        # w.transform = 'R'
    
        # Calculate Getis Ord statistic
        g_local = G_Local(raster_flat, w,transform='R')
        
        gband = g_local.Gs.reshape(1,raster_data.shape[1],raster_data.shape[2])
        gband_stand = g_local.Zs.reshape(1,raster_data.shape[1],raster_data.shape[2])
        pval = g_local.p_sim.reshape(1,raster_data.shape[1],raster_data.shape[2])
        
        lg[i,:,:] = gband
        lg_stand[i,:,:] = gband_stand
        lg_pvals[i,:,:] = pval
        
    # Return the results
    return g_local, gband, gband_stand, pval#, lg, lg_stand

# Set the distance threshold for spatial weights calculation
threshold = 11.2  # Adjust as needed

# Call the function to calculate Getis Ord statistic
getis_ord_results = calculate_getis_ord_statistic(coordinates, arr,threshold)

# Print the calculated statistics
# print("Getis Ord statistic values:", getis_ord_results.Gs)
# print("P-values:", getis_ord_results.p_sim)

