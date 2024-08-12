#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 11:49:57 2024

@author: Meredith
Figure of mean bottom reflectance data for the different bottom types 
at the GB Islands
Originally used in OSM 2024 talk.

"""
# %%

import pandas as pd
import matplotlib.pyplot as plt

# %%
# Make the plot

kelp = pd.read_csv('/Volumes/Mere/UMB/benthic_reflectance_spectra/refl_means/SL_R.csv',
                   skiprows=1)
browns = pd.read_csv('/Volumes/Mere/UMB/benthic_reflectance_spectra/refl_means/browns_R.csv',
                     skiprows=1)
reds = pd.read_csv('/Volumes/Mere/UMB/benthic_reflectance_spectra/refl_means/reds_R.csv',
                   skiprows=1)
greens = pd.read_csv('/Volumes/Mere/UMB/benthic_reflectance_spectra/refl_means/greens_R.csv',
                     skiprows=1)
sand = pd.read_csv('/Volumes/Mere/UMB/benthic_reflectance_spectra/refl_means/sand_R.csv',
                   skiprows=1)


# wl = kelp['Wavelength'][(kelp['Wavelength'] >= kelp['Wavelength'].min()) & (kelp['Wavelength'] <= 750)]



fig,ax = plt.subplots(1,1,figsize = (10,10))
ax.plot(kelp['Wavelength'],
        kelp['sl_R'],
        label='sugar kelp',
        color='saddlebrown',
        linewidth = 3)

ax.plot(browns['Wavelength'],
        browns['browns_R'],
        label='other brown algae',
        color='saddlebrown',
        linestyle = '--',
        linewidth=3)

ax.plot(reds['Wavelength'],
        reds['reds_R'],
        label='red algae',
        color='black',
        linestyle = ':',
        linewidth=2)

ax.plot(greens['Wavelength'],
        greens['greens_R'],
        label='green algae',
        color='black',
        linestyle = '-.',
        linewidth=2)

ax.plot(sand['Wavelength'],
        sand['WHBsand_R'],
        label='sand',
        color='black',
        linestyle = '-',
        linewidth=2)

ax.set_xlim([400,700])

ax.set_ylim([0.8,1.2])
ax.legend(fontsize=16,loc='upper left')

ax.xaxis.set_tick_params(
                labelsize=16)
ax.yaxis.set_tick_params(
                labelsize=16)
ax.set_xlabel('Wavelength (nm)',fontsize=18)
ax.set_ylabel('Normalized Reflectance',
              fontsize = 18)

fig.savefig('output/figs/norm_mean_reflectance.png',dpi=300)
plt.show()

