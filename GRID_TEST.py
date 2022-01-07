
import numpy as np
import netCDF4 as nc

from grid_KN import Grid as Grid_KN
from grid import Grid

import matplotlib.pyplot as plt
import matplotlib.colors as cl

import plotting as pt
import plotting_tools as ptt

import tools

from readData import readVariable
from varDict import getPlottingVars

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


TEST_PISO = False
if TEST_PISO:

	print(data.shape)

	data = ptt.maskBathyAll(data[-1], grid)
	u1=np.mean(data,axis=0)
	u2=np.mean(data[:26],axis=0)
	pt.plot1by2([u1,u2], vmin=[vmin,vmin], vmax=[vmax,vmax])


	data = np.mean(data, axis=2)
	#data = ptt.maskDraftYZ(data, grid, xi=10, timeDep=True)
	data = ptt.maskBathyYZ(data, grid, xi=0)
	print(data.shape)

	
	figsize=(5,4)
	dpi=100
	X = grid.YC[:,1]
	Y = grid.RC.squeeze()

		
	Nt = data.shape[0]
	cmap='coolwarm'
	levels = np.linspace(-1,1,9)
	
	#==

	pt.plot1by2([data, data], X=[X,X], Y=[Y,Y], vmin=[vmin, vmin], vmax=[vmax, vmax], mesh=False)
	quit()

#==

TEST_animate = True
if TEST_animate:

	path = '/home/michai/Documents/data/PISOMIP_002/run/'
	#path =  '/Users/mh115/Documents/BAS/data/PISOMIP_001/run/' 
	grid = Grid(path)
	X = grid.YC[:,1]/1000.
	Y = grid.RC.squeeze()

	var = '2D'; VAR = 'ETAN'
	#var = 'Theta'; VAR = 'THETA'
	#var = 'Uvel'; VAR = 'UVEL'
	#var = 'Vvel'; VAR = 'VVEL'
	#var = 'Wvel'; VAR = 'WVEL'
	vmin, vmax, cmap, title = getPlottingVars(var)

	data = readVariable(var, path, file_format='nc', meta=True)
	print(data)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[-10], 'fontdict':{'fontsize':14, 'color':'w'}}

	data = data[VAR][:]

	pt.contourf(data[-1,])
	quit()

	data = tools.boundData(data, vmin, vmax, scale=0.99999)
	
	
	for ti in range(data.shape[0]):
		data[ti,] = ptt.maskBathyAll(data[ti,], grid)
		
	data = np.mean(data, axis=3)

	xlabel = 'LATITUDE (km)'; ylabel = 'DEPTH (m)'
	pt.animate1by1(data, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data)

		
	quit()
	
	#==
	
	# Lat-lon plot.
	pt.plot1by2([var1[-1,0], var2[-1,0]], X=[grid.XC, grid.XC], Y=[grid.YC, grid.YC], titles=['rdmds', 'nc'], xlabels=['lon', 'lon'], ylabels=['lat', 'lat'], show=True)
	

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




