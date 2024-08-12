#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 13:49:43 2024

@author: Meredith
"""
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
import rasterio as rio
from rasterio.plot import show
import seaborn as sns

# %% 
# Function to find the closest pixel to a point and extract depth

def extract_depth(point, src):
    row, col = src.index(point.x, point.y)
    depth = src.read(1, window=((row, row + 1), (col, col + 1)))

    # spectra_dict = dict(zip(wl_columns, spectra))
    return depth

# %%
# Files 

# vegetation cover
sacfor_file = 'data/dropcam/gooseberries_sacfor_v2_subset.shp'
sacfor = gpd.read_file(sacfor_file)

# bathymetry file
bathy_file = 'data/lidar/2022_ngs_topobathy_chatham_ma_Job964557_reproj_NS.tif'

# %%
# update geodataframe to get elevation from the bathymetry file

with rio.open(bathy_file) as src:
    if 'Elev(m)' not in sacfor.columns:
        sacfor['Elev(m)'] = np.nan
        
    # Plot the bathymetry data
    fig, ax = plt.subplots(figsize=(10, 8))
    bathy_show = show(src, ax=ax, cmap='viridis', title='Bathymetry Data')
    
    
    # Plot the sacfor points on top of the bathymetry
    sacfor.plot(ax=ax, color='red', markersize=10, label='sacfor points')
    
    
    for idx, row in sacfor.iterrows():
        depth = extract_depth(row.geometry, src)

        # Update the corresponding columns with depth
        sacfor.at[idx, 'Elev(m)'] = depth

print(sacfor.tail())

plt.show()


# %% 
# Add a new column 'bare_pres/abs' based on 'veg #'

sacfor['bare_pres/abs'] = sacfor['veg #'].isin([0, 1, 2, 3]).astype(int)
sacfor['veg_corr'] = sacfor['veg #'] - sacfor['kel #']
sacfor['veg_corr_PA'] = sacfor['veg_corr'].isin([0, 1, 2, 3]).astype(int)

# %%
# plot depth density distribution

# bottom_type = ['veg_pres/a', 'kelp_pres/']#,'bare_pres/abs',]
bottom_type = ['veg_corr_PA', 'kelp_pres/','bare_pres/abs']
colors = ['green','brown','blue']

# Iterate through the five airlines
for b,btype in enumerate(bottom_type):
    # Subset to the airline
    subset = sacfor['Elev(m)'][sacfor[btype]==1]

    # Draw the density plot
    sns.distplot(subset, hist = True, kde = False,
                 bins= int(sacfor['Elev(m)'].max()-sacfor['Elev(m)'].min()/0.3),
                 # hist_kws = {'edgecolor':'white'},
                 # kde_kws = {'linewidth': 1},
                 color=colors[b],
                 label = btype
                 )

# Plot formatting
plt.legend(prop={'size': 10})
plt.title('Bottom Type Distribution')
plt.xlabel('Depth (m)')
plt.ylabel('Density')
plt.savefig('output/substrate_type_hist_NS.png',dpi=300)
plt.show()