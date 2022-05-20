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


path = '/home/michai/Documents/data/PAS_851/run/'
grid = Grid_PAS(path)


rho0 = 1028.5

drag = np.load(path + 'umom_drag_PAS851.npy')
conv = np.load(path + 'umom_conv_PAS851.npy')*rho0
taux = np.load(path + 'umom_oceTAUX_PAS851.npy')
TFS = np.load(path + 'umom_dpdx_PAS851.npy')
print(TFS.shape)


lats = [-76, -71.5]; lons = [245, 260]#[230, 270]#
latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

drag = ptt.maskBathyXY(drag, grid, 0, timeDep=False)*rho0
conv = ptt.maskBathyXY(conv, grid, 0, timeDep=False)
taux = ptt.maskBathyXY(taux, grid, 0, timeDep=False)
taux = ptt.maskDraftXY(taux, grid, 0, timeDep=False)
TFS = ptt.maskBathyXY(TFS, grid, 0, timeDep=False)

conv = tools.getSubregionXY(conv, latsi, lonsi)
taux = tools.getSubregionXY(taux, latsi, lonsi)
TFS = tools.getSubregionXY(TFS, latsi, lonsi)

print(np.ma.mean(drag))
print(np.ma.mean(conv))
print(np.ma.mean(taux))
print(np.ma.mean(TFS))

#TFS = tools.getSubregionXY(TFS, latsi, lonsi)

vmin = -.5e-1; vmax = -vmin
s=1.

vmin = [vmin, s*vmin, vmin]; vmax = [vmax, s*vmax, vmax]

pt.plot1by3([conv, taux, TFS], vmin=vmin, vmax=vmax)


