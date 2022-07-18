# slopeCurrent.py

# Look at zonal shelf slope current.

#========================================================== 

import time

import numpy as np
#import pandas as pd
#import xarray as xr
import netCDF4 as nc

from grid import Grid

import matplotlib.pyplot as plt

import plotting as pt

import tools

#import xarray as xr; ds = xr.open_dataset(path+fname)

#========================================================== 

path = '/Users/mh115/Documents/BAS/data/PAS_666/'
#path = '/Users/mh115/Documents/BAS/data/PAS_001/'

grid = Grid(path)

#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5
fname = 'stateUvel.nc'; var = 'UVEL'; vmin = -0.2; vmax = 0.2
#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2


ncfile = nc.Dataset(path+fname, 'r')
data = ncfile.variables[var]

ts = -1; ks = 10



quit()
lon = grid.getIndexFromLongitude(250)
XC = grid.XC
print(XC[0])
print(XC[0, lon])
print(XC[0, lon-1])
print(XC[0, lon+1])

