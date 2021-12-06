# yzTest.py

#========================================================== 

import numpy as np
#import pandas as pd
#import xarray as xr
import netCDF4 as nc

from grid import Grid

import matplotlib.pyplot as plt

import plotting_tools as ptt

#==========================================================

path = '/Users/mh115/Documents/BAS/data/PISOMIP_001/'
fname = 'stateTheta.nc'

grid = Grid(path)

# THETA(TIME, DEPTH, LATITUDE, LONGITUDE)

ncfile = nc.Dataset(path+fname, 'r')
print(ncfile.variables)
#print(ncfile.variables['DEPTH'][:])

ts = -1; i = 4
data = ncfile.variables['THETA']


data = data[ts, :, :, i]


data = ptt.maskDraftYZ(data, grid, subregion=False)

plt.figure(figsize=(5,2.5), dpi=200)
plt.gca().patch.set_color('.25')
plt.pcolormesh(grid.YC[:,0]/1.e3, grid.RC.squeeze(), data, cmap='jet')
plt.colorbar()
plt.show()

ncfile.close()

