#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:38:35 2024

@author: Meredith
"""
import geopandas as gpd
import pandas as pd
import rasterio as rio
import numpy as np
import matplotlib.pyplot as plt

# Function to find the closest pixel to a point and extract spectra
def extract_spectra(point, src, wl_columns):
    row, col = src.index(point.x, point.y)
    spectra = src.read(window=((row, row + 1), (col, col + 1))).flatten().tolist()

    spectra_dict = dict(zip(wl_columns, spectra))
    return spectra_dict

# sacfor file
sacfor_file = 'data/dropcam/gooseberries_sacfor_v2_subset.shp'
sacfor = gpd.read_file(sacfor_file)

# wl file
wl_file = 'data/headwall/wv_v2.csv'
wl = pd.read_csv(wl_file,names=['count','wl'])

# image file
img_file = 'data/headwall/nano_Rrs_deglint_smooth_mosaic_mean_v2_NS.tif'

# get x and y coordiantes from sacfor
sacfor['x'] = sacfor.geometry.x
sacfor['y'] = sacfor.geometry.y

# Create columns for each wavelength
wl_columns = ['wl_' + str(wavelength) for wavelength in wl['wl']]
sacfor[wl_columns] = np.nan

with rio.open(img_file) as src:
    arr = src.read(masked=True)
    plt.imshow(arr[50,:,:])
    
    
    for idx, row in sacfor.iterrows():
        spectra_dict = extract_spectra(row.geometry, src, wl_columns)

        # Update the corresponding columns with spectra values
        for column, value in spectra_dict.items():
            sacfor.at[idx, column] = value
            


