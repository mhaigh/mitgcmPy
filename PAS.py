import numpy as np

from grid_PAS import Grid as Grid_PAS

import matplotlib.pyplot as plt
import matplotlib.colors as cl

import plotting as pt
import plotting_tools as ptt

import PAS_tools
import tools

from readData import * 

#==========================================================

DATADIR = '/home/michael/Documents/data/'
EASlats = [-75.5, -70.5]; EASlons = [245, 262] 

MAIN = True
if MAIN:

	grid = Grid_PAS(DATADIR + 'PAS_851/run/')
	bathy = grid.bathy
	X = grid.XC; Y = grid.YC
	
	pathno = 0
	subr = 1
	if subr:
		# Get subregion of bathymetry.
		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)
		bathy = tools.getSubregionXY(bathy, latsi, lonsi)
		X, Y = grid.XYsubr(EASlons, EASlats)
		pathno = 1
	nY, nX = X.shape

	# Get lat, lon lists of continental slope
	slope_x, slope_y = PAS_tools.getSlopeContour(bathy, X, Y, pathno)
	
	# Get bearing of slope (using a centred differencing, one-sided at ends of slope).
	bearing = np.zeros(len(slope_x))
	bearing[0] = PAS_tools.getBearing(slope_y[0], slope_x[0], slope_y[1], slope_x[1])
	bearing[-1] = PAS_tools.getBearing(slope_y[-2], slope_x[-2], slope_y[-1], slope_x[-1])
	bearing[1:-1] = PAS_tools.getBearing(slope_y[:-2], slope_x[:-2], slope_y[2:], slope_x[2:])
	# bearing[i] from slope[i-1] and slope[i+1] 
	
	slope_x_nX = np.zeros(nX); slope_y_nX = np.zeros(nX)
	slope_xi = np.zeros(nX, dtype=int); slope_yi = np.zeros(nY, dtype=int)
	bearing_nX = np.zeros(nX)
	
	# For each X grid point, get the lat and lon of slope contour's nearest point.
	for i in range(nX):
		slope_xi[i] = int(np.argmin(np.abs(slope_x-X[0,i])))
		slope_x_nX[i] = slope_x[slope_xi[i]]
		slope_y_nX[i] = slope_y[slope_xi[i]]
		bearing_nX[i] = bearing[slope_xi[i]]		
		slope_yi[i] = int(np.argmin(np.argmin(slope_y_nX[i]-Y[:,i])))
	
	# For each timestep, load instances of u, v, rho.
	# Depth average beneath a rho contour. (27.95?)
	# Then take dot prodcut with unit vector in direction of slope's bearing.
	# Surface flow can be just the flow or estimated from SSH for geostrophic part.
	for ti in range(1):
		a = 1	
		
	# Get y-grid point closest to slope_y_nX
	# Read u, v at this point and n grid points north and south.
	# Get component of each in along-slope direction
			
	
		
		# For each y, get t
	plt.subplot(121)
	plt.contourf(X, Y, bathy)
	plt.scatter(slope_x_nX, slope_y_nX)
	plt.subplot(122)
	plt.plot(X[0], bearing_nX)
	plt.show()
	quit()
	
#==

# Define latitude of continental slope for each longitude.
SLOPE_CONTOUR = False
if SLOPE_CONTOUR:

	grid = Grid_PAS(DATADIR + 'PAS_851/run/')
	bathy = grid.bathy
	X = grid.XC; Y = grid.YC
	
	pathno = 0
	subr = 1
	if subr:
		# Get subregion of bathymetry.
		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)
		bathy = tools.getSubregionXY(bathy, latsi, lonsi)
		X, Y = grid.XYsubr(EASlons, EASlats)
		pathno = 1

	plt.subplot(121)
	cs = plt.contour(X, Y, bathy, levels=[-1000])
	x = cs.collections[0].get_paths()[pathno].vertices[:,0][::-1]
	y = cs.collections[0].get_paths()[pathno].vertices[:,1][::-1]
	#plt.scatter(x[:30], y[:30])
	plt.grid()
	
	print(len(x))
	print(X.shape)
	
	bs = []
	for i in range(len(x)-1):
		bs.append(PAS_tools.getBearing(y[i],x[i], y[i+1],x[i+1]))
	plt.subplot(122)
	plt.plot(bs)
	plt.grid()
	plt.show()
	quit() 
	
	
#==

latLonToCartesian1 = True
if latLonToCartesian1:

	grid = Grid_PAS(DATADIR + 'PAS_851/run/')
	bathy = grid.bathy
	X = grid.XC; Y = grid.YC
	
	#x, y, dx, dy = PAS_tools.latlon_to_xy(X,Y)
	b = PAS_tools.getBearing(1,1,1,2)
	print(b); quit()
	
	pt.plot1by2([x,y]);quit()
	
	pt.plot1by1(bathy, X=x, Y=y);quit() 
	cs = plt.contour(x, y, bathy, levels=[-1000])
	x = cs.collections[0].get_paths()[0].vertices[:,0]
	y = cs.collections[0].get_paths()[0].vertices[:,1]
	plt.show()
	plt.clf()
	
#==

latLonToCartesian2 = False
if latLonToCartesian2:

	grid = Grid_PAS(DATADIR + 'PAS_851/run/')
	bathy = grid.bathy
	lon = grid.XC; lat = grid.YC
	
	dA, x, y = PAS_tools.dA_from_latlon (lon, lat, periodic=False, return_edges=True)
	
	
	pt.plot1by2([x, y])
	
	
# OPTIONS:
# 1.
# For each point in scatter plot, get nearest grid point.
# Is there worry about double-counting some u,v with this method?

# 2.
# For each lon, get latitude of shelf slope and nearest grid point for u, v.
# 
	
	
	
	
	
	
	
	
# Original code for getting coordinates of continental slope.
#x = []; y = []
#for item in cs.collections:
#	print(1)
#	for i in item.get_paths():
#		v = i.vertices; x.append(v[:, 0]); y.append(v[:, 1])
