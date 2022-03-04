
# tools.py

import sys

import matplotlib.pyplot as plt

import numpy as np

#==========================================================

def maskEdges(data, size, dims=None):
	'''Given data mask its edges.
	Size of mask in each dimension determined by size and dims.'''
	
	if dims is None and len(data.shape)==2:
		return data[size:-size, size:-size]
	
	else:
		print('Error: tools.maskEdges. Option selected which has not been coded up.')
		sys.exit()
		
#==
	
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

#==

def depthAverage(data, Z, timeDep=True):
	'''Return depth average of data.
	Assumes input is time-dependent (timeDep=True) in which case depth is on second axis.
	If timeDep=False, assumes depth is on first axis.
	Averages over all depths, unless a maxLevel is provided.'''

	if timeDep:
		return np.trapz(data[:,::-1], Z[::-1], axis=1)# / (Z[0] - Z[-1])
	else:
		axis = 0
		return np.trapz(data[::-1], Z[::-1], axis=0)# / (Z[0] - Z[-1])

#==

def depthIntegral(data, grid, timeDep=True, SSH=None, norm=True):
	'''Return depth integral (sum) of data.'''
	
	if timeDep:
		axis = 1
	else:
		axis = 0

	if norm:
		norm = np.sum(grid.DRF * grid.hFacW, axis=0)
	else:
		norm = 1

	if SSH is not None:
		DRF0 = grid.DRF.squeeze()[0] + SSH
		if timeDep:
			T0 = data[:,0] * grid.hFacW[0]  * DRF0
			return	(T0 + np.sum(data[:,1:] * grid.hFacW[1:] * grid.DRF[1:], axis=axis)) / norm
		else: 
			T0 = data[0] * grid.hFacW[0] * DRF0
			return	(T0 + np.sum(data[1:] * grid.hFacW[1:] * grid.DRF[1:], axis=axis)) / norm
			
	else:
		return np.sum(data * grid.hFacW * grid.DRF, axis=axis) / norm

#==

def barotropicStreamfunction(u, grid, timeDep=True, SSH=None):
	'''Given zonal velocity, compute barotropic streamfunction.
	If timeDep, assumes axis ordering is (t, z, y, x).
	If not timeDep, assumes axis ordering is (z, y, x).'''
	
	# First, get depth-integrated mass transport.
	T = depthIntegral(u, grid, timeDep=timeDep, SSH=SSH)

	# Second, integrate from southern boundary to northern boundary.
	if timeDep:
		# Axes are (t, y, x).
		axis = 1
	else:
		# Axes are (y, x).
		axis = 0

	Y = grid.YC
	
	# Option using trapezium rule.
	#strfunc = np.zeros(T.shape)
	#if timeDep:
	#	for yi in range(0,Y.shape[0]):
	#		strfunc[:,yi] = - np.trapz(T[:,:yi+1,], Y[:yi+1,0], axis=axis)

	# Simple option using cumulative sum.
	if timeDep:
		strfunc = - np.cumsum(T * grid.DYG, axis=1)
	
	return strfunc

def meridStreamfunction(v, X, Z, timeDep=True):
	'''Given meridional velocity, return meridional overturning streamfunction.
	If timeDep, assumes axis ordering is (t, z, y, x).
	If not timeDep, assumes axis ordering is (z, y, x).'''

	if timeDep:	
		return

	else:
		vx = np.trapz(v, X, axis=2)
		vxz = np.zeros(vx.shape)
		Nz = vx.shape[0]
		for zi in range(Nz):
			vxz[Nz-zi-1] = np.trapz(vx[::-1][:zi+1], Z[::-1][:zi+1], axis=0) 
		
	return vxz

#==

def boundData(data, vmin, vmax, scale=None):
	'''Limit data array by vmin and vmax using np.where.
	scale used to ensure data less than or greater than (not equal) to bounds.
	This is messy fix for contourf plots.'''

	if scale is not None:
		data = np.ma.where(data>=vmax, scale*vmax, data)
		data = np.ma.where(data<=vmin, scale*vmin, data)
	else:
		data = np.ma.where(data>=vmax, vmax, data)
		data = np.ma.where(data<=vmin, vmin, data)
	return data 

#==

def diff(data, axis=0, dx=1):
	'''Implement central differencing on data along specified axis.'''

	ddata = np.zeros(data.shape)

	if axis == 0:
		ddata[1:-1,] = (data[2:,] - data[0:-2,]) / (2 * dx)
		ddata[0,] = (data[1,] - data[0,]) / dx
		ddata[-1,] = (data[-1,] - data[-2,]) / dx

	elif axis == 1:
		ddata[:, 1:-1,] = (data[:, 2:,] - data[:, 0:-2,]) / (2 * dx)
		ddata[:, 0,] = (data[:, 1,] - data[:, 0,]) / dx
		ddata[:, -1,] = (data[:, -1,] - data[:, -2,]) / dx

	return ddata

#==

def removeOscillations(data, lim):
	
	data = np.where(np.abs(data)<lim, 0, data)
	
	return data

#==

def ddy(f, dy):
	'''A crude meridional derivative function.
	Assumes y-axis is second from last.'''
	
	dfdy = np.zeros(f.shape)
	Ny = f.shape[-2]	

	dfdy[...,1:Ny-1,:] = (f[...,2:Ny,:] - f[...,0:Ny-2,:]) / (2 * dy)
	dfdy[...,0,:] = (f[...,1,:] - f[...,0,:]) / dy
	dfdy[...,Ny-1,:] = (f[...,Ny-1,:] - f[...,Ny-2,:]) / dy

	return dfdy

#==

def bottom(data, grid, cellPos, timeDep=True):
	'''Return value of data at bottom of domain.
	Must specify if at u, v, or h location of grid cells.'''

	if cellPos == 'u':
		hFac = grid.hFacW
	elif cellPos == 'v':
		hFac = grid.hFacS
	elif cellPoss == 'h':
		hFac = grid.hFacC
	else:
		print('Error: tools.bottom. Must provide valid cellPos.')
		sys.exit()

	Nt, Nz, Ny, Nx = data.shape

	# Get depth index of ocean floor.
	bottom = np.zeros((Ny, Nx))
	for k in range(1,Nz):
		bottom += hFac[k] > 0
	bottom = bottom.astype(int)

	# bottom is Ny, Nx array of z grid points.
	Y = [[i]*Nx for i in range(Ny)]
	Y = [i for subX in Y for i in subX]

	X = [[i]*Ny for i in range(Nx)]
	X = [X[j][i] for i in range(Ny) for j in range(Nx)]

	datab = data[:,bottom.flatten(),Y,X].reshape((Nt,Ny,Nx))

	# Commmented out a brute-force approach to above list comprehension solution.
	#datab = np.zeros((Nt,Ny,Nx))
	#for j in range(Ny):
	#	for i in range(Nx):
	#		datab[:,j,i] = data[:,bottom[j,i],j,i]

	return datab

#==


def EOSlinear(rho0, sBeta, tAlpha, S, Sref, T, Tref):

	Tref = np.array(Tref); Sref = np.array(Sref) 
	Tref = np.tile(Tref, (1,1,1)).T
	Sref = np.tile(Sref, (1,1,1)).T

	return rho0 * (sBeta * (S - Sref) - tAlpha * (T - Tref)) 

	
	
	
	
	

