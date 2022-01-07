# COMPARE_CS_SIMS.py

# Code for plotting and comparing outputs from contintental shelf (CS) simulations.

#==========================================================

import numpy as np

from grid import Grid

import matplotlib.pyplot as plt

import plotting as pt
import plotting_tools as ptt

import tools

from readData import readVariable
from varDict import getPlottingVars

#==========================================================

# SSH
SSH = True
if SSH:
	
	path1 = '/home/michai/Documents/data/PISOMIP_001/run/' 
	path2 = '/home/michai/Documents/data/PISOMIP_002/run/'

	var = '2D'; VAR = 'ETAN'
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.XC[0,:] / 1000.
	Y = grid1.YC[:,0] / 1000.
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Get time data. This assumes both data have same output freq.
	t = data1.variables['TIME'][:]
	print(t)
	quit()
	
	tscale = 86400. * 30
	tscale_name = 'month'
	t = data.variables['TIME'][:] / tscale; print(t)
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[-10], 'fontdict':{'fontsize':14, 'color':'w'}}

	data = data[VAR][:]

	pt.contourf(data[-1,])
	quit()




