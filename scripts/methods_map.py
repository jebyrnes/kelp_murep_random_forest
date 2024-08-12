#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:07:14 2024

@author: Meredith

This is a map that can be layered with different datasets collected at the GB
islands in 2022. Originally made for OSM 2024 talk
"""

#%% import packages

import geopandas as gpd
import numpy as np
import rasterio as rio
from rasterio.plot import show
import matplotlib.pyplot as plt
import contextily as ctx
from skimage import exposure
# import matplotlib.ticker as ticker

# %% Files

# hyperspectral tif file
img_file = 'data/headwall/nano_Rrs_deglint_smooth_mosaic_mean_v2.tif'

with rio.open(img_file) as src:
    arr = src.read(masked=True)
    
# Bathymetry
bathy_file = 'data/lidar/2022_ngs_topobathy_chatham_ma_Job964557_reproj.tif'

with rio.open(bathy_file) as bathy_src:
    bathy_arr = bathy_src.read(1,masked=True)
    nodata_mask = bathy_arr == bathy_arr.min()
    
    bathy_arr_nanmasked = np.ma.masked_array(bathy_arr,
                                             mask=nodata_mask)
    
    land_mask = bathy_arr > 0
    
    bathy_arr_masked = np.ma.masked_array(bathy_arr_nanmasked,
                                          mask= land_mask)


# sidescan survey area
ss_file = '/Users/Meredith/Library/CloudStorage/GoogleDrive-meredith.mcpherson@umb.edu/.shortcut-targets-by-id/1PJamAGTknM_kuVq8M53Xo4cAlWzDOk6I/kelp_murep/sidescan_surveys/SSS_Coverage/2022_SALEM_Goosberry.shp'
ss_survey = gpd.read_file(ss_file)
ss_survey = ss_survey.to_crs(src.crs)

# drop camera
dropcam_file = 'data/dropcam/gooseberries_sacfor_v2.shp'
dropcam = gpd.read_file(dropcam_file)
dropcam = dropcam.to_crs(src.crs)

# %% Make map

# create a GeoAxes instance with cartopy projection
fig, ax1 = plt.subplots(1, 1, figsize=(20, 20))

# plt the rgb of the 1m resolution imagery
# red (650nm) = band 114
# green (550nm) = band 96
# blue (450nm) = band 24
rgb = arr[[114, 75, 25], :, :]
gamma_corrected = exposure.adjust_gamma(rgb, 
                                        gamma=0.23)  

show(arr[50,:,:], 
      transform=src.transform, 
      zorder=3,
      alpha=1,
      ax=ax1)

# show(bathy_arr_masked,
#       transform=src.transform,
#       zorder = 2,
#       alpha = 0.6,
#       ax=ax1)
# # ax1.set_xlabel('Longitude', 
# #                fontsize = 30)
# # ax1.set_ylabel('Latitude',
# #                fontsize = 30)

# # ax1.tick_params(axis='both',
# #                 labelsize=20)

# # Define a function to format the latitude tick labels
# def lat_formatter(x, pos):
#     return f'{x:.3f}'

# # Set the x-axis tick locator and formatter
# # ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))
# # ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lat_formatter))

# # color the drop cam locations by shallow/deep

# dropcam.plot(ax=ax1, 
#             color = 'w',
#             # edgecolor = 'k', 
#             # label = 'Underwater imagery validation points',
#             #legend=True,
#             zorder=21)

ss_survey.plot(ax=ax1,
               color='grey',
               edgecolor='pink',
               alpha=0.3,
               # label = 'Sidescan Survey Area',
               #legend=True,
               zorder=1)

# plot the basemap of Gooseberry Islands
ctx.add_basemap(ax1, crs=dropcam.crs,source=ctx.providers.Esri.WorldImagery)

plt.legend(loc='lower right',
           fontsize = 20,
           labelcolor='w',
           facecolor = 'black')

ax1.set_xticks([])
ax1.set_yticks([])

fig.tight_layout()

fig.savefig('output/figs/Gooseberry_overview_sidescan_headwall.png',dpi=300)
plt.show()

