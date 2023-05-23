# PAS_tools.py

# tools.py

import sys

import matplotlib.pyplot as plt

import numpy as np
import math

#====================================================================================

# General helper function to get the area of each cell from latitude and longitude arrays giving the coordinates of the cell centres. Adapted from Kaitlin Naughten's mitgcm_python.
def dA_from_latlon (lon, lat, periodic=False, return_edges=False):

	deg2rad = np.pi / 180.
	rEarth = 6.371e6

	# Make sure they're 2D
	if len(lon.shape) == 1 and len(lat.shape) == 1:
		lon, lat = np.meshgrid(lon, lat)
	# Now make the edges
	# Longitude
	if periodic:
		lon_extend = wrap_periodic(lon, is_lon=True)
		lon_edges = 0.5*(lon_extend[:,:-1] + lon_extend[:,1:])
	else:
		lon_edges_mid = 0.5*(lon[:,:-1] + lon[:,1:])
		# Extrapolate the longitude boundaries
		lon_edges_w = 2*lon_edges_mid[:,0] - lon_edges_mid[:,1]
		lon_edges_e = 2*lon_edges_mid[:,-1] - lon_edges_mid[:,-2]
		lon_edges = np.concatenate((lon_edges_w[:,None], lon_edges_mid, lon_edges_e[:,None]), axis=1)
	dlon = lon_edges[:,1:] - lon_edges[:,:-1] 
	# Latitude
	lat_edges_mid = 0.5*(lat[:-1,:] + lat[1:,:])
	lat_edges_s = 2*lat_edges_mid[0,:] - lat_edges_mid[1,:]
	lat_edges_n = 2*lat_edges_mid[-1,:] - lat_edges_mid[-2,:]
	lat_edges = np.concatenate((np.expand_dims(lat_edges_s,0), lat_edges_mid, np.expand_dims(lat_edges_n,0)), axis=0)
	dlat = lat_edges[1:,:] - lat_edges[:-1,:]
	# Now convert to Cartesian
	dx = rEarth*np.cos(lat*deg2rad)*dlon*deg2rad
	dy = rEarth*dlat*deg2rad
	dA = dx*dy
	if return_edges:
		return dA, lon_edges, lat_edges
	else:
		return dA
		
#==

def latlon_to_xy(X, Y, returnGridSpacing=True, deg2rad=np.pi/180., rEarth=6.371e6):
	'''Convert lat-lon (X,y) grid to Cartesian (x, y) grid using the Haversine formula.'''
	
	X *= deg2rad; Y *= deg2rad
	
	nY, nX = X.shape
	x = np.zeros((nY, nX)); y = np.zeros((nY, nX))
	dx = np.zeros((nY-1, nX-1)); dy = np.zeros((nY-1, nX-1))
	for i in range(1,nX):
		for j in range(1,nY):
			dX = X[j,i] - X[j,i-1]
			dY = Y[j,i] - Y[j-1,i]
			
			# Compute grid for x
			a = np.cos(Y[j,i])**2 * np.sin(dX/2)**2
			dx[j-1,i-1] = 2 * rEarth * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 
			x[j,i] = x[j,i-1] + dx[j-1,i-1]
			
			# Repeat for y
			a = np.sin(dY/2)**2
			dy[j-1,i-1] = 2 * rEarth * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 
			y[j,i] = y[j-1,i] + dy[j-1,i-1] 
	
	if returnGridSpacing:
		return x, y, dx, dy
	else:
		return x, y
		
#==

def getBearing(lat1, lon1, lat2, lon2, deg2rad=np.pi/180.):
	
	lat1 = deg2rad * np.copy(lat1)
	lat2 = deg2rad * np.copy(lat2)
	lon1 = deg2rad * np.copy(lon1)
	lon2 = deg2rad * np.copy(lon2)

	dLon = (lon2 - lon1)

	y = np.sin(dLon) * np.cos(lat2)
	x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dLon)

	brng = np.arctan2(y, x)

	brng = np.rad2deg(brng)

	return brng
	
#==

def getSlopeContour(bathy, X, Y, pathno, level=-1000):

	cs = plt.contour(X, Y, bathy, levels=[level])
	x = cs.collections[0].get_paths()[pathno].vertices[:,0][::-1]
	y = cs.collections[0].get_paths()[pathno].vertices[:,1][::-1]
	plt.clf()
	
	return x, y



			
	
