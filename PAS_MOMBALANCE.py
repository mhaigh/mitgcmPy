import numpy as np

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import matplotlib.pyplot as plt
import matplotlib.colors as cl

import plotting as pt
import plotting_tools as ptt

import tools

from readData import readVariable
from varDict import getPlottingVars

import time

#==========================================================

# Zonal momentum balance for PAS.
# Output: ice-mediated wind stress, bottom friction, TFS.
# If necessary: full pressure gradient contribution, zonal mom. tendency.

#==========================================================

g = 9.81
rho0 = 1028.5

ts = 107; te = 502

path_root = '/data/oceans_output/shelf/pahol/mitgcm/'
run = 'PAS_851/run/'
# With PAS_851, Paul suggests looking at dates 1979-2012.
# Pre 1979 is spin up, post 2012 is cold anomaly.

path = path_root + run

# Load grid
grid = Grid_PAS(path)

lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

#==========================================================

## Advection ##

u = readVariable('UVEL', path, file_format='nc', meta=False)[ts:te+1]
v = readVariable('VVEL', path, file_format='nc', meta=False)[ts:te+1]

u = tools.interp(u, 'u')
v = tools.interp(v, 'v')

uu = np.mean(u*u, axis=0)
uv = np.mean(u*v, axis=0)

w = readVariable('WVEL', path, file_format='nc', meta=False)[ts:te+1]
uw = np.mean(u*w, axis=0)

conv = - tools.ddx(uu, grid.DXG) - tools.ddy(uv, grid.DYG) - tools.ddz(uw, grid.DRF)

conv = np.sum(conv*grid.hFacC*grid.DRF, axis=0)

conv = np.ma.filled(conv, fill_value=0)
np.save('umom_conv_PAS851', conv)

quit()

#==

## PGF ##

dpdx = readVariable('Um_dPhiX', path, file_format='nc', meta=False)[ts:te+1]
dpdx = np.mean(dpdx, axis=0)

# Take depth integral
dpdx = np.sum(dpdx*grid.hFacW*grid.DRF, axis=0)

dpdx = np.ma.filled(dpdx, fill_value=0)
np.save('umom_dpdx_PAS851', dpdx)

quit()

#==

## Surface stress ##

taux = readVariable('oceTAUX', path, file_format='nc', meta=False)[ts:te+1]
taux = np.mean(taux, axis=0)

taux = np.ma.filled(taux, fill_value=0)
np.save('umom_oceTAUX_PAS851', taux)

#==

## Bottom friction ##

drag = tools.computeBotDragQuadr2(path, grid)[ts:te+1]
drag = np.mean(drag, axis=0)

drag = np.ma.filled(drag, fill_value=0)
np.save('umom_drag_PAS851', drag)


#==

## TFS ##

depth = - grid.bathy
SSH = readVariable('ETAN', path, file_format='nc', meta=False)[ts:te+1]
depth = depth + SSH

# Load PHIBOT.
PHIBOT = readVariable('PHIBOT', path, file_format='nc', meta=False)[ts:te+1]

# Compute full bottom pressure.
Pb = depth * rho0 * g + PHIBOT * rho0
Pb = np.mean(Pb, axis=0)

# Depth gradient
Hx = tools.ddx(depth, grid.DXG)

# TFS
TFS = Pb * Hx

#TFS = np.mean(TFS, axis=0)
TFS = np.ma.filled(TFS, fill_value=0)
np.save('umom_TFS_PAS851', TFS)



