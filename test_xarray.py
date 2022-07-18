# test.py

# Test that ensures xarray and dependencies are functioning.
# Calls basic xarray functions.
# Testing shows orders of magnitude more time needed with this xarray approach compared to netCDF4 and numpy.

#==========================================================

import time
import sys
import os 

import numpy as np
#import pandas as pd
import xarray as xr
import netCDF4 as nc

import matplotlib.pyplot as plt

#==========================================================

#data_np = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20]})
#data_pd = xr.DataArray(pd.Series(range(3), index=list("abc"), name="foo"))

now = time.time()

path = '/Users/mh115/Documents/BAS/data/PAS_001/'
fname = 'stateTheta.nc'

data = xr.open_dataset(path+fname)
print(data)
print(data.isel(TIME=1))

data = data + 10

data.close()

print(time.time() - now)
