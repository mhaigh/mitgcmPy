import numpy as np

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import matplotlib.pyplot as plt
import matplotlib.colors as cl

import plotting as pt
import plotting_tools as ptt

import tools

from readData import readVariable
from varDict import getPlottingVars, getTrefSref

import time

#==========================================================

g = 9.81
rho0 = 1028.5

path_root = '/data/oceans_output/shelf/pahol/mitgcm/'
run = 'PAS_851/run/'
# With PAS_851, Paul suggests looking at dates 1979-2012.
# Pre 1979 is spin up, post 2012 is cold anomaly.

path = path_root + run

# Load grid
grid = Grid_PAS(path)
depth = - grid.bathy

lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

# Load UVEL.
PHIBOT = readVariable('PHIBOT', path, file_format='nc', meta=True)

ts = 107; te = 502
TIME = u['TIME'][:]
print(time.ctime(TIME[ts]))
print(time.ctime(TIME[te]))

PHIBOT = PHIBOT['PHIBOT']
print(PHIBOT.shape)

Pb = depth * rho0 * g + PHIBOT * rho0

Pb = np.mean(Pb[ts:te+1], axis=0)
Pb = np.ma.filled(Pb, fill_value=0)
np.save('Pbmean_PAS851', Pb)



