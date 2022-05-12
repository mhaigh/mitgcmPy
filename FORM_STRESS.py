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
taux = np.load(path + 'umom_oceTAUX_PAS851.npy')
TFS = np.load(path + 'umom_TFS2_PAS851.npy')


lats = [-76, -71.5]; lons = [245, 260]#[230, 270]#
latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

drag = ptt.maskBathyXY(drag, grid, 0, timeDep=False)*rho0
taux = ptt.maskBathyXY(taux, grid, 0, timeDep=False)
TFS = ptt.maskBathyXY(TFS, grid, 0, timeDep=False)

print(np.mean(drag))
print(np.mean(taux))
print(np.mean(TFS))

TFS = tools.getSubregionXY(TFS, latsi, lonsi)

vmin = -1.e-1; vmax = -vmin
s=1.e7
vmin = [vmin, vmin, s*vmin]; vmax = [vmax, vmax, s*vmax]

pt.plot1by3([drag, taux, TFS], vmin=vmin, vmax=vmax)


