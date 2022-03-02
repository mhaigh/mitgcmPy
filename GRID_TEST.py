
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

#from io import readData

#==========================================================

# Different variables defined at different parts of the grid/stencil.

# build_land_mask returns 2D boolean array of land.
##	 hfac is fraction of vertical cell occupied by water, so if zero, then its land. Full column=0 -> coastline.
# To get ice shelf, get ocean mask, multiply by (surface) regions where hfac < 1.

# What to do today?
# input code. Can read either nc or mitgcm binaries, function takes filed name as input.
# check they give the same result.
# maybe also optionally return a dictonary with plotting properties such as vmax/vmin and cmap, title...
# animate theta.
# make sure grid looks right from hfacs.


# Load which nc file?
#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
#fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2
# PHIHYD
baroclinicEddies = False
if baroclinicEddies:
	
	path = '/home/michai/Documents/data/MCS_002/run/'
	grid = Grid(path)
	X = grid.XC[1,:]/1000.		
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	ts = 16; level = 24; yi = -5

	VAR = 'THETA'
	theta = readVariable(VAR, path, file_format='nc', meta=False)
	VAR = 'VVEL'
	vvel = readVariable(VAR, path, file_format='nc', meta=False)
	
	# Get time data. This assumes both data have same output freq.
	text = 'month ' + str(ts+1) + '; Z = ' + str(Z[level]) + ' m'
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}, None]

	#vvel_ = vvel[ts, :, -5]
	#
	vvel_ = vvel[ts, level, ]; theta_ = theta[ts, level,]	
	pt.plot1by2([vvel_, theta_], X=[X,X], Y=[Y,Y], xlabels=['X (km)', 'X (km)'], ylabels=['Y (km)', ''], titles=['VVel (m/s)', 'Theta (deg. C)'], text_data=text_data); quit()

	
	vvel_ = vvel[ts, :, yi,]
	theta_ = theta[ts, :, yi] - np.tile(np.mean(theta[ts, :, yi], axis=1), [theta.shape[3],1]).T
	
	vmax = 0.04; vmin = - vmax; vvel_ = tools.boundData(vvel_, vmin, vmax, scale=0.99999)	
	vmax = 0.05; vmin = - vmax; theta_ = tools.boundData(theta_, vmin, vmax, scale=0.99999)

	pt.plot1by2([vvel_, theta_], X=[X,X], Y=[Z,Z], xlabels=['X (km)','X (km)'], ylabels=['Depth (m)', ''], titles=['VVel (m/s)', 'Theta eddy (deg. C)'])

	quit()

#==

TEST_depthAverage = False
if TEST_depthAverage:

	path = '/home/michai/Documents/data/MCS_001/run/'

	var = 'Vvel'; VAR = 'VVEL'
	vmin, vmax, cmap, title = getPlottingourcVars(var)

	# Get grid.
	grid = Grid(path)
	X = grid.XC[0,:] / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()	

	# Read data
	data = readVariable(var, path, file_format='nc', meta=True)
	data = data[VAR][0,]

	#v_zav1 = np.trapz(data, Z, axis=0)
	#data = ptt.maskBathyAll(data, grid)
	#v_zav2 = np.trapz(data,-Z, axis=0)
	#v_zav3 = tools.depthAverage(data, Z, timeDep=False)

	
	vstr1 = tools.meridStreamfunction(data, X, Z)

	data = ptt.maskBathyAll(data, grid)
	vstr2 = tools.meridStreamfunction(data, X, Z)
	vstr1 = ptt.maskBathyYZ(vstr2, grid)
	#vstr2 = np.cumsum(np.sum(data[::-1], axis=2), axis=0) * (Z[0]-Z[1]) * (X[1]-X[0]); vstr2 = vstr2[::-1]

	pt.plot1by2([vstr1, vstr2], X=[Y,Y], Y=[Z,Z])	
	
	quit()

#==

TEST_animate = False
if TEST_animate:

	#path = '/home/michai/Documents/data/MCS_002/run/'

	path = '/home/michai/Documents/data/MCS_036/run/'
	#pathG = '/home/michai/Documents/data/MCS_018/run/'

	grid = Grid(path)
	#grid = Grid_PAS(path)

	X = grid.YC[:,1]/1000.
	Y = grid.RC.squeeze()

	#VAR = 'ETAN'
	#VAR = 'RHOAnoma'
	VAR = 'THETA'
	#VAR = 'PHIHYD'
	#VAR = 'DFrE_TH';
	#VAR = 'WVELTH'#','UVELTH','VVELTH','WVELTH', 'TOTTTEND'
	#VAR = 'SALT'	
	#VAR = 'UVEL'
	#VAR = 'VVEL'
	#VAR = 'WVEL'

	vmin, vmax, cmap, title = getPlottingVars(VAR)

	data = readVariable(VAR, path, file_format='nc', meta=True)
	print(data)

	tscale = 86400.; tscale_name = 'day'
	t = data.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[-10], 'fontdict':{'fontsize':14, 'color':'w'}}
	
	MEAN = False
	if not MEAN:
		data = data[VAR][:]
		for ti in range(data.shape[0]):
			data[ti,] = ptt.maskBathyAll(data[ti,], grid)
		#data = np.ma.mean(data, axis=3)
		data = data[...,1]
	
	else:
		data = np.mean(data[VAR][:], axis=3)
		data = ptt.maskBathyYZ(data, grid, timeDep=True)


	data = tools.boundData(data, vmin, vmax, scale=0.99999)
	#data = tools.removeOscillations(data, 1.e-3)
	
	# PLOT.

	xlabel = 'LATITUDE (km)'; ylabel = 'DEPTH (m)'
	#pt.plot1by1(data[-1], X, Y, xlabel=xlabel, ylabel=ylabel, title=title, cmap=cmap, vmin=vmin, vmax=vmax, mesh=False); quit()

	pt.animate1by1(data, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data)

		
	quit()
	
	#==
	
	# Lat-lon plot.
	pt.plot1by2([var1[-1,0], var2[-1,0]], X=[grid.XC, grid.XC], Y=[grid.YC, grid.YC], titles=['rdmds', 'nc'], xlabels=['lon', 'lon'], ylabels=['lat', 'lat'], show=True)
	
#==

animateSurface = False
if animateSurface:

	#path = '/home/michai/Documents/data/PAS_666/run/'
	path = '/home/michai/Documents/data/MCS_038/run/'
	
	grid = Grid(path)
	#grid = Grid_PAS(path)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.

	#VAR = 'ETAN'
	#VAR = 'RHOAnoma'
	VAR = 'THETA' 
	#VAR = 'PHIHYD'
	#VAR = 'DFrE_TH';
	#VAR = 'WVELTH'#','UVELTH','VVELTH','WVELTH', 'TOTTTEND'
	#VAR = 'SALT'	
	#VAR = 'UVEL'
	#VAR = 'VVEL'
	#VAR = 'VVEL'
	#VAR = 'WVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)

	data = readVariable(VAR, path, file_format='nc', meta=True)
	print(data)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

	data = data[VAR][:]
	
	
	#plt.plot(np.mean(data[-1,], axis=1)); plt.show(); quit()

	# DON'T CHANGE THIS. CHANGE ONE IN IF STATEMENT IF YOU WANT.
	level = 0

	if VAR != 'ETAN':
		level = 26
		data = data[:,level]

	vmax = 1.; vmin = -0.4
	data = tools.boundData(data, vmin, vmax, scale=0.99999)
	data = ptt.maskBathyXY(data, grid, level, timeDep=True)
	

	#data = np.mean(data, axis=3)


	xlabel = 'LON (km)'; ylabel = 'LAT (km)'
	pt.animate1by1(data, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, outname='animate1by1Surf.mp4')

	quit()
	
	#==
	
	# Lat-lon plot.
	pt.plot1by2([var1[-1,0], var2[-1,0]], X=[grid.XC, grid.XC], Y=[grid.YC, grid.YC], titles=['rdmds', 'nc'], xlabels=['lon', 'lon'], ylabels=['lat', 'lat'], show=True)
	
#==

animateQuivers = True
if animateQuivers:

	#path = '/home/michai/Documents/data/PAS_666/run/'
	path = '/home/michai/Documents/data/MCS_036/run/'
	
	level = 24

	grid = Grid(path)
	#grid = Grid_PAS(path)	

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	data1 = readVariable(VAR, path, file_format='nc', meta=True)
	TIME = data1['TIME'][:]
	data1 = data1[VAR]

	VAR = 'VVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	data2 = readVariable(VAR, path, file_format='nc', meta=False)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = TIME / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

	data1 = ptt.maskBathyXY(data1[:, level,], grid, level, timeDep=True)
	data2 = ptt.maskBathyXY(data2[:, level,], grid, level, timeDep=True)

	contour = grid.bathy
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)
	
	# Get subregion.
	#xlims = [100, -100]; ylims = [0, -100]
	#data1 = data1[...,ylims[0]:ylims[1], xlims[0]:xlims[1]]
	#data2 = data2[...,ylims[0]:ylims[1], xlims[0]:xlims[1]]
	#X = X[xlims[0]:xlims[1]]; Y = Y[ylims[0]:ylims[1]]
	#contour = contour[ylims[0]:ylims[1], xlims[0]:xlims[1]]

	# Sample rate & level
	d = 6

	title = '(u, v) & bathymetry; Z = ' + str(Z[level]) + ' m'

	data1 = data1[..., ::d, ::d]; data2 = data2[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]
	
	#pt.plot1by2([data1[0], data2[0]]); quit()

	pt.animate1by1quiver(data1, data2, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, text_data=text_data, title=title)

#==

animateSSHs = False
if animateSSHs:
	# Animate zonal mean SSHs for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michai/Documents/data/'

	# These are cases without CS walls.
	#paths = [path_root+'MCS_013/run/', path_root+'MCS_018/run/', path_root+'MCS_012/run/', path_root+'MCS_016/run/']
	# These are cases with CS walls.	
	paths = [path_root+'MCS_015/run/', path_root+'MCS_019/run/', path_root+'MCS_014/run/', path_root+'MCS_017/run/']
	labels = ['05cos', '025cos', '025sin2', '0125sin2']
	
	grid = Grid(paths[0])
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'zonal mean SSH (m)'

	VAR = 'ETAN'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, paths[0], file_format='nc', meta=True)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}

	# Add first zonal mean SSH to data list.
	data_tmp = ptt.maskBathyXY(data_tmp[VAR][:], grid, 0, timeDep=True)
	data.append(np.mean(data_tmp[:,1:-1], axis=2))
	
	# Now get other SSHs.
	for di in range(len(paths)-1):
		data_tmp = readVariable(VAR, paths[di+1], file_format='nc', meta=False)
		data_tmp = ptt.maskBathyXY(data_tmp, grid, 0, timeDep=True)
		data.append(np.mean(data_tmp[:,1:-1], axis=2))

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data)
	quit()

#==

animateSSHs1sim = False
if animateSSHs1sim:
	# Animate SSH at different longitudes for a single MITgcm run. 

	path_root = '/home/michai/Documents/data/'
	path = path_root + 'MCS_034/run/'
	#xis = [1, 190, 80, 120]
	xis = [1, 5, 230, 239]
	labels = ['x = ' + str(xi*2.5) + ' km' for xi in xis]

	grid = Grid(path)
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'zonal mean SSH (m)'

	VAR = 'ETAN'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, path, file_format='nc', meta=True)


	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}
	
	data_tmp = data_tmp[VAR][:]

	data2 = readVariable('VVEL', path, file_format='nc', meta=True)['VVEL'][:]
	data_tmp = (20. + data_tmp) * data2[:,0]

	# Now get other SSHs.
	for xi in xis:
		data.append(data_tmp[...,1:-1,xi])


	vmin = -2e-1; vmax = -vmin
	pt.animateLine(data, X=Y[1:-1], xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, vmin=vmin, vmax=vmax)
	quit()

#==

animateUvelSurfs = False
if animateUvelSurfs:
	# Animate zonal mean SSHs for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michai/Documents/data/'

	# These are cases without CS walls.
	#paths = [path_root+'MCS_013/run/', path_root+'MCS_018/run/', path_root+'MCS_012/run/', path_root+'MCS_016/run/']
	# These are cases with CS walls.	
	#paths = [path_root+'MCS_015/run/', path_root+'MCS_019/run/', path_root+'MCS_014/run/', path_root+'MCS_017/run/']
	#labels = ['05cos', '025cos', '025sin2', '0125sin2']

	# These are for comparing with and without beta.
	paths = [path_root+'MCS_015/run/', path_root+'MCS_022/run/', path_root+'MCS_014/run/', path_root+'MCS_023/run/']
	labels = ['05cos wall', '05cos wall beta', '025sin2 wall', '025sin2 wall beta']
	
	grid = Grid(paths[0])
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'zonal mean Uvel (m)'

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, paths[0], file_format='nc', meta=True)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}

	# Add first zonal mean SSH to data list.
	data_tmp = ptt.maskBathyXY(data_tmp[VAR][:,0], grid, 0, timeDep=True)
	data.append(np.mean(data_tmp[:,1:-1,-100:], axis=2))
	
	# Now get other SSHs.
	for di in range(len(paths)-1):
		data_tmp = readVariable(VAR, paths[di+1], file_format='nc', meta=False)
		data_tmp = ptt.maskBathyXY(data_tmp[:,0], grid, 0, timeDep=True)
		data.append(np.mean(data_tmp[:,1:-1,-100:], axis=2))

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data)
	quit()


#==

plotBathy = False
if plotBathy:

	#path = '/home/michai/Documents/data/PAS_666/run/'
	path = '/home/michai/Documents/data/MCS_038/run/'
	
	grid = Grid(path)
	bathy = grid.bathy

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	xlabel = 'LON (km)'; ylabel = 'LAT (km)'

	bathy = ptt.maskBathyXY(bathy, grid, 0, timeDep=False)

	pt.plot1by1(bathy, X=X, Y=Y, xlabel=xlabel, ylabel=ylabel, mesh=True); quit()

#==

TEST_troughTransport = False
if TEST_troughTransport:

	# Look at volume transport through window west of trough and possibly east of trough.
	# Compare with volume transport through trough and nearby upwelling.
	# If we want time series of volume transport through these slices, have to take care with grid stencil and
	# variable cell thicnkesses. Take a look at xgcm if needed.
	
	# Remaining to do:
	# w through an x,y slice.
	# u through a second y,z slice.
	# volume transport through these slices. (scale by hfac, dz; use u and v)
	# Time series of these.
	# Onto github, onto HPC.
	
	# Four windows:
	lat1 = [-71.8, -71.2]; lon1 = 360-114.5; depth1 = [-250, -750]
	lat2 = -71.6; lon2 = [360-115, 360-112]; depth2 = [-250, -750]
	lat3 = [-71.8, -71.2]; lon3 = [360-115, 360-112]; depth3 = -500
	lat4 = lat1; lon4 = 360-112; depth4 = depth1

	# Convert physical values to indices.
	lats1 = grid.getIndexFromLat(lat1); depths1 = grid.getIndexFromDepth(depth1)
	xi1= grid.getIndexFromLon(lon1)
	
	lons2 = grid.getIndexFromLon(lon2); depths2 = grid.getIndexFromDepth(depth2)
	yi2 = grid.getIndexFromLat(lat2)
	
	lats3 = grid.getIndexFromLat(lat3); lons3 = grid.getIndexFromLon(lon3)
	zi3 = grid.getIndexFromDepth(depth3)
	
	lats4 = grid.getIndexFromLat(lat4); depths4 = grid.getIndexFromDepth(depth4)
	xi4= grid.getIndexFromLon(lon4)
	
	#==
	
	# Get zonal transport
	fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
	#data = io.readData('Uvel', path,)
	ncfile = nc.Dataset(path+'', 'r'); data = ncfile.variables[var]
	ut = tools.zonalTransport(data, grid)
	
	# Zonal transport through windows 1 and 4.
	ut1 = tools.getSubregionYZ(ut[..., xi1], grid, lats1, depths1)
	ut4 = tools.getSubregionYZ(ut[..., xi4], grid, lats4, depths4)

	# Get meridional transport through window 2.
	fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	vt = tools.meridTransport(data, grid)
	vt2 = tools.getSubregionXZ(vt[:, :, yi2, ], grid, lons2, depths2)
	
	# Get vertical transport through window 3.
	#fname = 'stateWvel.nc'; var = 'WVEL'; cmap = 'coolwarm'; vmax = 0.01; vmin = -vmax
	#ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	wt = tools.meridTransport(data, grid)
	wt3 = tools.getSubregionXY(wt[:, zi3, ], grid, lats3, lons3)
		
	# Transports are already scaled by grid size, so can just sum to get total transport.
	ut1 = np.sum(ut1, axis=(1,2))
	ut4 = np.sum(ut4, axis=(1,2))
	vt2 = np.sum(vt2, axis=(1,2))
	wt3 = np.sum(wt3, axis=(1,2))

	data = [ut1, ut4, vt2, wt3]
	labels = ['ut window 1', 'ut window 4', 'vt window 2', 'wt window 3']
	TIME = ncfile.variables['TIME'][:]
	pt.timeSeries(data, TIME=TIME, labels=labels, ylabel='Volume transport')
	
#==

TEST_manyWindows = False
if TEST_manyWindows:

	ti = -1
	
	# Window 1
	lat1 = [-71.8, -71.2]; lon1 = 360-114.5; depth1 = [-250, -750]
	
	# Window 2
	lat2 = -71.6; lon2 = [360-115, 360-112]; depth2 = [-250, -750]
	
	# Window 3.
	lat3 = [-71.8, -71.2]; lon3 = [360-115, 360-112]; depth3 = -500
	
	# Window 4.
	lat4 = lat1; lon4 = 360-112; depth4 = depth1
	
	#==
	
	# 1. Get zonal flow slice through window 1 and window4
	lats1 = grid.getIndexFromLat(lat1)
	depths1 = grid.getIndexFromDepth(depth1)
	xi1= grid.getIndexFromLon(lon1)
	
	lats4 = grid.getIndexFromLat(lat4)
	depths4 = grid.getIndexFromDepth(depth4)
	xi4= grid.getIndexFromLon(lon4)
	
	fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	u1 = tools.getSubregionYZ(data[ti, ..., xi1], grid, lats1, depths1)
	u4 = tools.getSubregionYZ(data[ti, ..., xi4], grid, lats4, depths4)
	
	# Zonal transport.
	#u_transport = tools.zonalTransport(data, grid)
	
	#==
	
	# 2. Get meridional flow through window 2.
	lons2 = grid.getIndexFromLon(lon2)
	depths2 = grid.getIndexFromDepth(depth2)
	yi = grid.getIndexFromLat(lat2)
	
	fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	v = tools.getSubregionXZ(data[ti, :, yi, ], grid, lons2, depths2, Ydim=False)
	
	#==
	
	# 3. Vertical flow through window 3 (don't have NC file!)
	lats3 = grid.getIndexFromLat(lat3)
	lons3 = grid.getIndexFromLon(lon3)
	zi = grid.getIndexFromDepth(depth3)
	
	#fname = 'stateWvel.nc'; var = 'WVEL'; cmap = 'coolwarm'; vmax = 0.01; vmin = -vmax
	#ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	w = tools.getSubregionXY(data[ti, zi, ], grid, lats3, lons3)

	#==
	
	# Mask bathymetry 
	u1 = ptt.maskBathyYZ(u1, grid, xi=xi1, subregion=True, lats=lats1, depths=depths1)
	u4 = ptt.maskBathyYZ(u4, grid, xi=xi4, subregion=True, lats=lats4, depths=depths4)
	v = ptt.maskBathyXZ(v, grid, yi=yi, subregion=True, lons=lons2, depths=depths2, Ydim=False)
	w = ptt.maskBathyXY(w, grid, zi, subregion=True, lons=lons3, lats=lats3)
		
	# Get 1D subregion arrays. Y,Z for u-plot and X,Z for v-plot. 
	Xsubr = grid.Xsubr1D(lon2)
	Ysubr = grid.Ysubr1D(lat1)
	Zsubr = grid.Zsubr1D(depth1)
	
	pt.plot1by2([u1, v], X=[Ysubr, Xsubr], Y=[Zsubr, Zsubr], titles=['u window1', 'v window2'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lat', 'lon'], ylabels=['depth', ''])

	pt.plot1by2([u1, u4], X=[Ysubr, Ysubr], Y=[Zsubr, Zsubr], titles=['u window1', 'u window4'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lat', 'lat'], ylabels=['depth', ''])
	
	pt.plot1by2([w, v], X=[Xsubr, Xsubr], Y=[Ysubr, Zsubr], titles=['v window3', 'v window2'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lon', 'lon'], ylabels=['lat', 'depth '])
	
	quit()

#==

TEST_getSubregionXYZ = False
if TEST_getSubregionXYZ:
	
	# In this test, produce 3D (allow for additional time dimension) fields in subregion.
	# Apply to zonal flows over PITW shelf break.
	# Volume transport or just zonal velocity?
	# Hoevmoller diagram of U(z), could be averaged over some small y-range.

	None

#==

TEST_getSubregionYZ = False
if TEST_getSubregionYZ:

	# Current mask functions require full 2D slice.
	# Extend this so it can take partial 2D slice (have to get corresponding subregion in grid).
	# And extend so can work with general number of dimensions.

	# Load potential temp. Other options above.
	fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]

	ti = -1
	xi = grid.getIndexFromLon(245.5)
	
	#depth_vals = [0, -2000]; lat_vals = [-72.5, -70.5]
	
	# AssmanEtal2013 vals
	depth_vals = [0, -2000]; lat_vals = [-75, -65]; lon_vals = [360-116, 360-111]
	
	lats = grid.getIndexFromLat(lat_vals, xi=xi)
	depths = grid.getIndexFromDepth(depth_vals)
		
	data1 = tools.getSubregionYZ(data, grid, lats, depths)[ti, ..., xi]
	data2 = tools.getSubregionYZ(data[ti, ..., xi], grid, lats, depths)

	data1 = ptt.maskBathyYZ(data1, grid, xi=xi, subregion=True, lats=lats, depths=depths)
	data2 = ptt.maskBathyYZ(data2, grid, xi=xi, subregion=True, lats=lats, depths=depths)
		
	Ysubr = grid.Ysubr1D(lat_vals)
	Zsubr = grid.Zsubr1D(depth_vals)
	#data1 = ptt.maskBathyYZ(data1, grid, xi=xi)

	pt.plot1by2([data1, data2], X=[Ysubr, Ysubr], Y=[Zsubr, Zsubr], titles=['data1', 'data2'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lat', 'lat'], ylabels=['depth', ''], show=False, save=True)

#==

TEST_getSubregionXY = False
if TEST_getSubregionXY:

	ti = -1
	ki = 10

	lons = (245, 260)
	lats = (-72.5, -70.5)
	data1 = tools.getSubregionXY(data, grid, lons, lats)[ti, ki]
	data2 = tools.getSubregionXY(data[ti, ki, ], grid, lons, lats)

	pt.plot1by2(X, Y, data1, data2, titles=['data1', 'data2'])

#==

TEST_bathy = False
if TEST_bathy:
	bathy = grid.bathy
	draft = grid.draft

	landC = grid.landC
	iceC = grid.iceC
	
	
	pt.plot1by2([draft, grid_KN.draft], titles=['Bathymetry', 'KN Bathymetry'])
	
	pt.plot1by2([bathy, grid_KN.bathy], titles=['Draft', 'KN Draft'])

	quit()
	
	
#==
	
TEST_getIndex = False
if TEST_getIndex:

	depths = [0, -2000]
	lats = [-72.5, -70.5]
	lons = [240, 250]
	
	print(grid.getIndexFromLat(lats, xi=0))
	print(grid.getIndexFromLat(lats[0], xi=0))
	print(grid.getIndexFromLat(lats[1], xi=0))

	print(grid.getIndexFromDepth(depths))
	print(grid.getIndexFromDepth(depths[0]))
	print(grid.getIndexFromDepth(depths[1]))
		
	print(grid.getIndexFromLon(lons, yi=0))
	print(grid.getIndexFromLon(lons[0], yi=0))
	print(grid.getIndexFromLon(lons[1], yi=0))
	
	quit()
	

	
	
#==

TEST_readData = False
if TEST_readData:

	path =  '/Users/mh115/Documents/BAS/data/PISOMIP_001/run/'
	grid = Grid(path)
	
	var1 = readVariable('Theta', path, file_format='rdmds')
	var2 = readVariable('Theta', path, file_format='nc')

	#==
	
	ts = 1
		
	# Lat-depth plot.
	data1 = np.mean(var1[ts], axis=2); data2 = np.mean(var2[ts], axis=2)
	data1 = ptt.maskBathyYZ(data1, grid, xi=10); data2 = ptt.maskBathyYZ(data2, grid, xi=10)
	print(data1.shape)
	data1 = ptt.maskDraftYZ(data1, grid, xi=10); data2 = ptt.maskDraftYZ(data2, grid, xi=10)
	
	pt.plot1by2([data1, data2], X=[grid.YC[:,0], grid.YC[:,0]], Y=[grid.RC.squeeze(), grid.RC.squeeze()], titles=['rdmds', 'nc'], xlabels=['lat', 'lat'], ylabels=['depth', ''], show=True, figsize=(8,3))
		
	quit()
	
	#==
	
	# Lat-lon plot.
	pt.plot1by2([var1[-1,0], var2[-1,0]], X=[grid.XC, grid.XC], Y=[grid.YC, grid.YC], titles=['rdmds', 'nc'], xlabels=['lon', 'lon'], ylabels=['lat', 'lat'], show=True)




