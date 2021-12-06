
# tools.py

import sys

import numpy as np

#==========================================================
	
def getSubregionXY(data, grid, lons, lats):
	'''Given data with lat-lon dimensions in final two (or only two)
	axes, reduce lateral range to range given by lats and lons.
	Assumes lats and lons are tuples with two indices.'''
		
	return data[..., lons[0]:lons[1]+1, lats[0]:lats[1]+1]
	
#==

def getSubregionXZ(data, grid, lons, depths, Ydim=False):
	'''Similar to above: return subregion of data in (y,z)-plane.
	If data is 3- or 4-dimensional, there has to be y-dim between
	z- and x-dims.
	By default this assumes that data does not have an X-dimension.'''
	
	# Different scenarios based on number of dims in data...
	
	# Assume that only dimensions are z, x.
	if len(data.shape) == 2:
		return data[depths[0]:depths[1]+1, lons[0]:lons[1]+1]
	
	# Else, assume dimensions are (t, z, y, x) or (t, z, x).
	else:
		if Ydim:
			return data[..., depths[0]:depths[1]+1, :, lons[0]:lons[1]+1]
		else:
			return data[..., depths[0]:depths[1]+1, lons[0]:lons[1]+1]
		
#==
	
def getSubregionYZ(data, grid, lats, depths, Xdim=False):
	'''Similar to above: return subregion of data in (y,z)-plane.
	By default this assumes that data does not have a Y-dimension.'''
	
	# Different scenarios based on number of dims in data...
	
	# Assume that only dimensions are z, y.
	if len(data.shape) == 2:
		return data[depths[0]:depths[1]+1, lats[0]:lats[1]+1]
	
	# Else, assume dimensions are (t, z, y, x) or (t, z, y).
	else:
		if Xdim:
			return data[..., depths[0]:depths[1]+1, lats[0]:lats[1]+1, :]
		else:
			return data[..., depths[0]:depths[1]+1, lats[0]:lats[1]+1]
		
#==

def zonalTransport(u, grid):
	'''Return zonal transport in given region, taking care with grid stencil.'''
	
	return u * grid.DYG * grid.DRF * grid.hFacW

#==

def meridTransport(v, grid):
	'''Return meridional transport.'''
	
	return v * grid.DXG * grid.hFacS * grid.DRF
	
#==

def vertTransport(w, grid):
	''''Return vertical transport.'''
	
	return w * grid.RAC
	
#==
	
def zonalTransportTEST(u, grid):
	'''Return zonal transport in given region, taking care with grid stencil.'''
	
	level = 10
	time = 0
	print(u.shape)
	print(grid.DYG.shape)
	print(grid.DRF.shape)
	print(grid.hFacW.shape)

	NT, NZ, NY, NX = u.shape
	ud = np.zeros((NT, NZ, NY, NX))
	#NY = 200; NX = 10
	#for ki in range(NZ):
	#	print(ki)
	#	for yi in range(NY):
	#		for xi in range(NX):
	#			ud[0,ki,yi,xi] = u[0,ki,yi,xi] * grid.DYG[yi,xi]# * grid.DRF[ki,0,0] * grid.hFacW[ki,yi,xi]
	
	ud = np.zeros((NT, NZ, NY, NX))
	for yi in range(NY):
		for xi in range(NX):
			ud[time,level, yi,xi] = u[time,level,yi,xi] * grid.DYG[yi,xi] * grid.DRF[level,0,0] * grid.hFacW[level,yi,xi]			
			

	ub = u * grid.DYG * grid.DRF * grid.hFacW

	
	import plotting as pt
	pt.plot1by2([ud[time,level,]-ub[time,level,], ub[time,level,]], titles=['From loop', 'broadcasting'])
	quit()
	
	return ud, ub
	
	
	
	
	
	
	
	
	
	

