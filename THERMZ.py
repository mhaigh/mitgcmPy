# THERMZ.py

# For a number of MITgcm runs, produce time-mean isotherm height.	

#==========================================================

import numpy as np

from grid import Grid

import tools

from readData import readVariable
from varDict import getPlottingVars, getTrefSref

import time

#==========================================================

rootdir = '/data/oceans_output/shelf/michai/mitgcm/'
#rootdir = '/home/michai/Documents/data/'

outdir = rootdir + 'THERMZnpy/'

runs = ['MCS_'+str(i) for i in [108, 114, 115, 116, 117, 118, 119, 120]]

rho = 1030; Cp = 3974.0
T0 = -2

ts = 108; te = 120
THERM = -0.5; THERMt = 'm05'

for run in runs:
	
	path = rootdir + run + '/run/'

	grid = Grid(path)
	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	bathy = grid.bathy

	T = np.mean(readVariable('THETA', path, meta=False)[ts:te], axis=0)
	HC, HC0 = tools.heatContentShelf(T, grid)
	print(run, HC, HC-HC0)
	
	ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True, timeDep=False)

	np.save(outdir+'ThermZ_'+THERMt+'_'+run, ThermZ)






