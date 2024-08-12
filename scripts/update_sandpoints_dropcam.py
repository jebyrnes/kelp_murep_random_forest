#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 10:06:55 2024

@author: Meredith

Adding manually identified sand points to the sacfor validation shape file


"""

#%% 

# import packages
import os
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
import rasterio as rio

# change working dir
os.chdir('/Users/Meredith/Library/CloudStorage/GoogleDrive-meredith.mcpherson@umb.edu/.shortcut-targets-by-id/1PJamAGTknM_kuVq8M53Xo4cAlWzDOk6I/kelp_murep/rs_imagery_processing/classification/random_forest')

#%% 

# load the shapefiles to edit

# original validation shapefile
val_file = 'data/dropcam/gooseberries_sacfor.shp'

# csv file with new sand points
# there is an X and Y column but no geometry.
# geometry must be created
sand_file = 'data/dropcam/sand_points_manual.csv'

# imagery file
img_file = 'data/headwall/nano_Rrs_deglint_smooth_mosaic_mean_v2.tif'

# open the files as geopandas dataframes
val_gpd = gpd.read_file(val_file)
sand_gpd = gpd.read_file(sand_file)

# convert columns to float (Elev(m), kel #, veg #, and veg #corr)
float_columns = ['Elev(m)', 'kel #', 'veg #', 'veg #corr']
sand_gpd[float_columns] = sand_gpd[float_columns].astype(float)


# sand points were obtained from headwall imagery with certain crs
# get crs of headall image to convert sand points
with rio.open(img_file) as src:
    arr = src.read()
    src_crs = src.crs

# set the sand file crs to src_crs
sand_gpd.crs = src_crs

# get crs of validation 
val_crs = val_gpd.crs

# create point geometry from X and Y for new points
geometry = [Point(xy) for xy in zip(sand_gpd['X'], sand_gpd['Y'])]

# Update geometry column
sand_gpd['geometry'] = geometry

# convert new gpd to validation file crs 
sand_gpd = sand_gpd.to_crs(val_crs)

# add the new points to the original validation gpd as pandas dataframe
val_pd = pd.concat([val_gpd, sand_gpd], ignore_index=True)

# convert back to geodataframe
val_gpd = gpd.GeoDataFrame(val_pd, geometry='geometry', crs=val_crs)

# drop extra columns
val_gpd = val_gpd.drop(['X','Y'],axis=1)

# save the updated gpd
out_file = 'data/dropcam/gooseberries_sacfor_v2.shp'
val_gpd.to_file(out_file)
