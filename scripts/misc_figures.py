#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 10:11:42 2024

@author: Meredith

Testing map making with basempa...
Don't work great for small scale because 
Some good script in here though.
"""


import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

import contextily as ctx
from shapely.geometry import box

# %%
# Create the map of Salem Sound
# Draw a box around study site
# lame kinda

fig,ax = plt.subplots(1,1,figsize=(15,15))
m = Basemap(projection = 'merc',
            resolution = 'h', llcrnrlat =41,
            urcrnrlat = 45, llcrnrlon = -72,
            urcrnrlon = -66, lat_ts = 43,
            ax=ax)
            
m.drawcoastlines()
m.drawmapboundary(fill_color='lightblue')
m.fillcontinents( lake_color='lightblue')
m.drawstates()


parallels = np.arange(41,45+1,2)
m.drawparallels(parallels, labels=[False,True,True,False],size=18)

meridians = np.arange(-72,-66+1,2.)
m.drawmeridians(meridians,labels=[True,False,False,True],size=18)

# x,y = m(-70.8 , 45.5)
# ss = m.plot(x,y,marker='*', color='r',markersize = 25)
plt.savefig('output/figs/GOM.png',dpi=300)
plt.show()

# %%

# Global map with salem sound in center

fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
m = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
# m.shadedrelief()
m.drawcoastlines(linewidth=0.25)
m.drawcountries(linewidth=0.25)
m.fillcontinents(color='lightgrey',lake_color='lightblue')
# draw the edge of the map projection region (the projection limb)
m.drawmapboundary(fill_color='lightblue')
# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0,360,30))
m.drawparallels(np.arange(-90,90,30))

plt.savefig('output/figs/global.png',dpi=300)
plt.show()

# %% 
# Using the Mass shapefile to plot coastline

# mass coastline shapefile
mass = gpd.read_file('/Users/Meredith/Desktop/massgis-coast25k-arc-shapefile/GISDATA_COAST25K_ARC.shp')
# mass.plot()

# get the bbox of salem sound
#llat,llon,ulat,ulon
# bbox = [42.45,-71,42.6,-70.5]
bbox = [-70.9, 42.45, -70.7, 42.6]

mass_subset = mass[mass.within(box(*bbox))]
# mass_subset = mass_subset.to_crs(epsg=3857)
# Load your spatial data (shapefile or GeoDataFrame)

fig,ax = plt.subplots(1,1,figsize=(10,10))
# Plot the data using geopandas

ax = mass_subset.plot(figsize=(10, 10),color='w')

# Add a basemap using contextily
ctx.add_basemap(ax, crs=mass_subset.crs, source=ctx.providers.CartoDB.Positron)
# ctx.add_basemap(ax, crs=mass_subset.crs, source=ctx.providers.Stamen.TerrainBackground,zoom=18)
# ctx.add_basemap(ax,source=ctx.providers.MapTiler.Basic,zoom=10)
# ax.set_xticks([])
# ax.set_yticks([])
plt.savefig('output/figs/salem_sound.png',dpi=300)
plt.show()
























