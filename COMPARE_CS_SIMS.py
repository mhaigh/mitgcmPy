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
SSH = 0
if SSH:
	
	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

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

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]

	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'w'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyXY(data1, grid1, zi=0)
	data2 = ptt.maskBathyXY(data2, grid2, zi=0)

	titles = ['SSH (m) periodic', 'SSH (m) shelf wall']
	xlabels = ['x (km)', 'x (km)']
	ylabels = ['y (km)', '']

	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], vmin=[vmin,vmin], vmax=[vmax,vmax], titles=titles, text_data=text_data, xlabels=xlabels, ylabels=ylabels)
	quit()

#==

# Plots of zonal-mean theta in Y,Z-plane.
THETA_YZ = 0
if THETA_YZ:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Theta'; VAR = 'THETA'
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.YC[:,0] / 1000.
	Y = grid1.RC.squeeze()
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]

	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[-2], 'fontdict':{'fontsize':14, 'color':'w'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)

	# Take zonal means after masking.
	data1 = np.mean(data1, axis=2)
	data2 = np.mean(data2, axis=2)

	titles = ['THETA (deg. C) periodic', 'THETA (deg. C) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, vmin=[vmin, vmin], vmax=[vmax, vmax], text_data=text_data, xlabels=xlabels, ylabels=ylabels)
	quit()

#==

# Vertically averaged vertical velocity.
UVEL_zav = 0
if UVEL_zav:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Uvel'; VAR = 'UVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.XC[1,:] / 1000.
	Y = grid1.YC[:,0] / 1000.
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]
	
	
	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)
	
	data1 = tools.depthAverage(data1, grid1.RC.squeeze())
	data2 = tools.depthAverage(data2, grid2.RC.squeeze())
	
	vmin /= 5; vmax /= 5
	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	titles = ['UVEL_zav (m/s) periodic', 'UVEL_zav (m/s) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, vmin=[vmin, vmin], vmax=[vmax, vmax], xlabels=xlabels, ylabels=ylabels)
	quit()

#==

# Vertically averaged vertical velocity.
VVEL_zav = 0
if VVEL_zav:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Vvel'; VAR = 'VVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.XC[1,:] / 1000.
	Y = grid1.YC[:,0] / 1000.
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]
	
	
	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)
	
	data1 = tools.depthAverage(data1, grid1.RC.squeeze())
	data2 = tools.depthAverage(data2, grid2.RC.squeeze())

	vmin *= 5; vmax *= 5
	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	titles = ['VVEL_zav (m/s) periodic', 'VVEL_zav (m/s) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, vmin=[vmin, vmin], vmax=[vmax, vmax], xlabels=xlabels, ylabels=ylabels)
	quit()

#==

# Vertically averaged vertical velocity.
WVEL_zav = 0
if WVEL_zav:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Wvel'; VAR = 'WVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.XC[1,:] / 1000.
	Y = grid1.YC[:,0] / 1000.
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]
	
	
	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)
	
	data1 = tools.depthAverage(data1, grid1.RC.squeeze())
	data2 = tools.depthAverage(data2, grid2.RC.squeeze())

	vmin *= 5; vmax *= 5
	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	titles = ['WVEL_zav (m/s) periodic', 'WVEL_zav (m/s) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, vmin=[vmin, vmin], vmax=[vmax, vmax], xlabels=xlabels, ylabels=ylabels)
	quit()

#==

# Plots of zonal-mean wvel in Y,Z-plane.
UVEL_YZ = 0
if UVEL_YZ:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Uvel'; VAR = 'UVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.YC[:,0] / 1000.
	Y = grid1.RC.squeeze()
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]

	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[-2], 'fontdict':{'fontsize':14, 'color':'w'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)

	# Take zonal means after masking.
	data1 = np.mean(data1, axis=2)
	data2 = np.mean(data2, axis=2)

	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	titles = ['UVEL (m/s) periodic', 'UVEL (m/s) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, vmin=[vmin, vmin], vmax=[vmax, vmax], xlabels=xlabels, ylabels=ylabels)
	quit()

#==
# Plots of zonal-mean wvel in Y,Z-plane.
VVEL_YZ = 1
if VVEL_YZ:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Vvel'; VAR = 'VVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.YC[:,0] / 1000.
	Y = grid1.RC.squeeze()
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]

	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[-2], 'fontdict':{'fontsize':14, 'color':'w'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)

	# Take zonal means after masking.
	data1 = np.mean(data1, axis=2)
	data2 = np.mean(data2, axis=2)

	vmin /= 5; vmax /= 5
	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	titles = ['VVEL (m/s) periodic', 'VVEL (m/s) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, vmin=[vmin, vmin], vmax=[vmax, vmax], xlabels=xlabels, ylabels=ylabels)
	quit()

#==

# Plots of zonal-mean wvel in Y,Z-plane.
WVEL_YZ = 0
if WVEL_YZ:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	var = 'Wvel'; VAR = 'WVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)
	X = grid1.YC[:,0] / 1000.
	Y = grid1.RC.squeeze()
	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]

	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[-2], 'fontdict':{'fontsize':14, 'color':'w'}}, None]

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1)
	data2 = ptt.maskBathyAll(data2, grid2)

	# Take zonal means after masking.
	data1 = np.mean(data1, axis=2)
	data2 = np.mean(data2, axis=2)

	data1 = tools.boundData(data1, vmin, vmax)
	data2 = tools.boundData(data2, vmin, vmax)

	titles = ['WVEL (m/s) periodic', 'WVEL (m/s) shelf wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, vmin=[vmin, vmin], vmax=[vmax, vmax], xlabels=xlabels, ylabels=ylabels)
	quit()

#==

overturningStr = 0
if overturningStr:

	path1 = '/home/michai/Documents/data/MCS_001/run/' 
	path2 = '/home/michai/Documents/data/MCS_002/run/'

	#var = 'Wvel'; VAR = 'WVEL'
	#var = 'Uvel'; VAR = 'UVEL'
	var = 'Vvel'; VAR = 'VVEL'	
	vmin, vmax, cmap, title = getPlottingVars(var)

	# Get grid.
	grid1 = Grid(path1)
	grid2 = Grid(path2)

	
	# Read data
	data1 = readVariable(var, path1, file_format='nc', meta=True)
	data2 = readVariable(var, path2, file_format='nc', meta=True)

	# Extract time and data arrays.
	t = data1.variables['TIME'][:]
	data1 = data1[VAR][:]; data2 = data2[VAR][:]

	# Get time data. This assumes both data have same output freq.
	ts = min(data1.shape[0], data2.shape[0]) - 1
	text = 'month ' + str(ts+1)

	data1 = data1[ts]; data2 = data2[ts]
	data1 = ptt.maskBathyAll(data1, grid1); data2 = ptt.maskBathyAll(data2, grid2)
	
	# Get (Y,Z) slices of meridional streamfunction.
	X1 = grid1.XC[0,:]; X2 = grid2.XC[0,:]
	Z1 = grid1.RC.squeeze(); Z2 = grid2.RC.squeeze()
	data1 = tools.meridStreamfunction(data1, X1, Z1)
	data2 = tools.meridStreamfunction(data2, X2, Z2)	
	
	# Data needs to be masked again.
	data1 = ptt.maskBathyYZ(data1, grid1)
	data2 = ptt.maskBathyYZ(data2, grid2, xi=10)

	X = grid1.YC[:,0] / 1000.
	Y = grid1.RC.squeeze()
	titles = ['Merid. Streamfunction (ms) periodic', 'Merid. Streamfunction (ms) wall']
	xlabels = ['y (km)', 'y (km)']
	ylabels = ['depth (m)', '']
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[-2], 'fontdict':{'fontsize':14, 'color':'w'}}, None]
	
	vmax1 = 8000; vmin1 = -vmax1; vmax2 = 16000; vmin2 = -vmax2
	data1 = tools.boundData(data1, vmin1, vmax1); data2 = tools.boundData(data2, vmin2, vmax2)

	pt.plot1by2([data1, data2], X=[X,X], Y=[Y,Y], titles=titles, text_data=text_data, xlabels=xlabels, ylabels=ylabels, vmin=[vmin1, vmin2], vmax=[vmax1, vmax2])
	quit()



