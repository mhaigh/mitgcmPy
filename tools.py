# tools.py

import sys

import matplotlib.pyplot as plt

import numpy as np

#====================================================================================

def maskEdges(data, size, dims=None):
	'''Given data mask its edges.
	Size of mask in each dimension determined by size and dims.'''
	
	if dims is None and len(data.shape)==2:
		return data[size:-size, size:-size]
	
	else:
		print('Error: tools.maskEdges. Option selected which has not been coded up.')
		sys.exit()
		
#==
	
def getSubregionXY(data, lats, lons):
	'''Given data with lat-lon dimensions in final two (or only two)
	axes, reduce lateral range to range given by lats and lons.
	Assumes lats and lons are tuples with two indices.'''
		
	return data[..., lats[0]:lats[1]+1, lons[0]:lons[1]+1]
	
#==

def getSubregionXZ(data, lons, depths, Ydim=False):
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
	
def getSubregionYZ(data, lats, depths, Xdim=False):
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

def heatContentShelf(theta, grid, shelf=500, rho=1030, Cp=3974.0, T0=-2, y0=96, z0=25, troughC=True):
	'''Return scalar heat content on shelf.
	Assumes theta is snapshot or time-mean value.'''
	
	# Make bathy and z 3D fields.
	bathy = grid.bathy
	z = grid.RC
	bathy = np.tile(bathy, (z.shape[0], 1, 1))
	z = np.tile(z, (1, bathy.shape[1], bathy.shape[2]))

	# Heat content times volume elements
	HC = Cp * rho * (theta - T0) * grid.DXG * grid.DYG * grid.hFacC * grid.DRF
	
	# Same for initial heat content
	theta0 = -1.8 * np.ones(HC.shape)
	HC0 = Cp * rho * (theta0 - T0) * grid.DXG * grid.DYG * grid.hFacC * grid.DRF
	
	HC = np.where(z<bathy, 0, HC)
	HC0 = np.where(z<bathy, 0, HC0)
	
	# First get all heat south of certain latitude and above the troughs.
	HCtmp = HC[:z0,:y0,:]
	HCtmp0 = HC0[:z0,:y0,:]
	
	TCH = np.sum(HCtmp)
	TCH0 = np.sum(HCtmp0)
	#print(grid.RC.squeeze()[:z0])
		
	# Add in heat content in central trough.
	if troughC:
		xw = 82; xe = 156
		ys = 1; yn = 87		
		HCtmp = HC[z0:, ys:yn, xw:xe]
		HCtmp0 = HC0[z0:, ys:yn, xw:xe]
		TCH += np.sum(HCtmp)
		TCH0 += np.sum(HCtmp0)	

	return TCH, TCH0
	
#==

def zonalTransport(u, grid):
	'''Return zonal transport in given region, taking care with grid stencil.'''
	
	return u * grid.DYG * grid.DRF * grid.hFacW

#==

def meridTransport(v, grid, yi=None):
	'''Return meridional transport.'''
	
	if yi is None:
		return v * grid.DXG * grid.hFacS * grid.DRF
	
	else:
		return v * grid.DXG[yi] * grid.hFacS[:,yi] * grid.DRF[:,0]
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
	
	hFac = grid.hFacC
	
	if timeDep:
		axis = 1
	else:
		axis = 0

	if norm:
		normval = np.sum(grid.DRF * hFac, axis=0)
	else:
		normval = 1

	if SSH is not None:
		DRF0 = grid.DRF.squeeze()[0] + SSH
		if norm:
			normval = normval + SSH
		if timeDep:
			T0 = data[:,0] * hFac[0] * DRF0
			return	(T0 + np.sum(data[:,1:] * hFac[1:] * grid.DRF[1:], axis=axis)) / normval
		else: 
			T0 = data[0] * hFac[0] * DRF0
			return	(T0 + np.sum(data[1:] * hFac[1:] * grid.DRF[1:], axis=axis)) / normval
			
	else:
		normval = np.where(normval==0.0, 1, normval)
		return np.sum(data * hFac * grid.DRF, axis=axis) / normval

#==

def barotropicStreamfunction(u, grid, timeDep=True, SSH=None, norm=False):
	'''Given zonal velocity, compute barotropic streamfunction.
	If timeDep, assumes axis ordering is (t, z, y, x).
	If not timeDep, assumes axis ordering is (z, y, x).'''
	
	# First, get depth-integrated mass transport.
	T = depthIntegral(u, grid, timeDep=timeDep, SSH=SSH, norm=norm)

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

	# Simple option using cumulative sum. Convert to SV.
	return - 1.e-6 * np.cumsum(T * grid.DYG, axis=axis)

#==

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

def ddx(f, dx=1):
	'''A crude meridional derivative function.
	Assumes x-axis is last.'''
	
	dfdx = np.zeros(f.shape)
	Nx = f.shape[-1]

	if isinstance(dx, (np.ndarray, tuple, list)):
		dx_tmp = dx.copy()
	else:
		# Assume dy is scalar.
		dx_tmp = dx * np.ones(f.shape)

	dfdx[...,1:Nx-1] = (f[...,2:Nx] - f[...,0:Nx-2]) / (2 * dx_tmp[...,1:Nx-1])
	dfdx[...,0] = (f[...,1] - f[...,0]) / dx_tmp[...,0]
	dfdx[...,Nx-1] = (f[...,Nx-1] - f[...,Nx-2]) / dx_tmp[...,Nx-1]

	return dfdx

#==

def ddy(f, dy=1):
	'''A crude meridional derivative function.
	Assumes y-axis is second from last.'''
	
	dfdy = np.zeros(f.shape)
	Ny = f.shape[-2]

	if isinstance(dy, (np.ndarray, tuple, list)):
		dy_tmp = dy.copy()
	else:
		# Assume dy is scalar.
		dy_tmp = dy * np.ones(f.shape)

	dfdy[...,1:Ny-1,:] = (f[...,2:Ny,:] - f[...,0:Ny-2,:]) / (2 * dy_tmp[...,1:Ny-1,:])
	dfdy[...,0,:] = (f[...,1,:] - f[...,0,:]) / dy_tmp[...,0,:]
	dfdy[...,Ny-1,:] = (f[...,Ny-1,:] - f[...,Ny-2,:]) / dy_tmp[...,Ny-1,:]

	return dfdy

#==

def ddz(f, dz=1, axis=-3):
	'''A crude vertical derivative function.
	Assumes z-axis is third from last.
	If axis is elsewhere, override by providing axis.
	This then calls ddx or ddy.'''
	
	# If f doesn't have x or y dependence, axis should be -1 or -2.
	if axis == -1:
		return ddx(f, dz)
	elif axis == -2:
		return ddy(f, dz)

	# f has (x,y,z) dependence.
	elif axis == -3:		
		dfdz = np.zeros(f.shape)
		Nz = f.shape[-3]

		if isinstance(dz, (np.ndarray, tuple, list)):
			dz_tmp = dz.copy()
		else:
			# Assume dy is scalar.
			dz_tmp = dz * np.ones(f.shape)

		dfdz[...,1:Nz-1,:,:] = (f[...,2:Nz,:,:] - f[...,0:Nz-2,:,:]) / (2 * dz_tmp[...,1:Nz-1,:,:])
		dfdz[...,0,:,:] = (f[...,1,:,:] - f[...,0,:,:]) / dz_tmp[...,0,:,:]
		dfdz[...,Nz-1,:,:] = (f[...,Nz-1,:,:] - f[...,Nz-2,:,:]) / dz_tmp[...,Nz-1,:,:]

		return dfdz

	else:
		print('Error: tools.ddz. axis must be -1, -2 or -3.')
		sys.exit()
#==

def bottom(data, grid, cellPos, timeDep=True):
	'''Return value of data at bottom of domain.
	Must specify if at u, v, or h location of grid cells.'''

	if cellPos == 'u':
		hFac = grid.hFacW
	elif cellPos == 'v':
		hFac = grid.hFacS
	elif cellPos == 'h':
		hFac = grid.hFacC
	else:
		print('Error: tools.bottom. Must provide valid cellPos.')
		sys.exit()

	if timeDep:
		Nt, Nz, Ny, Nx = data.shape
	else:
		Nz, Ny, Nx = data.shape

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

	if timeDep:
		datab = data[:,bottom.flatten(),Y,X].reshape((Nt,Ny,Nx))
	else:
		datab = data[bottom.flatten(),Y,X].reshape((Ny,Nx))

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

	
#==

def interp(data, edge):
	'''Interpolate data from cell edges to cell centres.
	Interpolate from either u points or v points.'''

	data_interp = np.zeros(data.shape)

	if edge == 'u':
		print('Interpolating from u points. Averaging between masked/unmasked neighbours OK if no-slip/no-normal flow, provided masked areas set to zero. Assuming periodic in x.' )
		data_interp[...,:-1] = 0.5 * (data[...,:-1] + data[...,1:])
		data_interp[...,-1] = 0.5 * (data[...,-1] + data[...,0])

	elif edge == 'v':
		print('Interpolating from v points. Averaging between masked/unmasked neighbours OK if no-slip, provided masked areas set to zero. Assuming periodic in y.' )
		data_interp[...,:-1,:] = 0.5 * (data[...,:-1,:] + data[...,1:,:])
		data_interp[...,-1,:] = 0.5 * (data[...,-1,:] + data[...,0,:])

	elif edge == 'w':
		print('Interpolating from w points. Bottom value is held fixed.' )
		data_interp[...,:-1,:,:] = 0.5 * (data[...,:-1,:,:] + data[...,1:,:,:]) 
		data_interp[...,-1,:,:] = data[...,-1,:,:]

	else:
		print('Error: tools.interp. Must provide valid edge value.')
		sys.exit()

	return data_interp
	
#==

def get05cosWind(nx, ny):
	'''Return default 05cos wind with unit magnitude.'''

	taux = np.zeros((ny,nx));

	for j in range(1,ny+1):
		taux[j-1,:] = - np.cos(np.pi*(j-1)/(ny-1))

	return taux

#==


def getShiftedWind(nx, ny, shift, d_scale=False):
	'''Return shifted 05cos wind.'''

	taux = np.zeros((ny, nx));

	if d_scale:
		d = (ny - 1 - 2*shift) / (ny - 1 + 2*shift)
	else:
		d = 1

	for j in range(1, int(np.floor(ny/2)+1-shift)):
		taux[j-1,:] = - d * np.cos(np.pi*(j-1)/(ny-1-2*shift))

	for j in range(int(np.floor(ny/2)+1-shift), ny+1):
		taux[j-1,:] = - np.cos(np.pi*(j-1+2*shift)/(ny-1+2*shift))

	return taux

#==

def getTranslatedWind(nx, ny, tr):
	'''Return translated 05cos wind.'''

	taux = np.zeros((ny,nx))
	taux_tr = np.zeros((ny,nx))

	for j in range(1,ny+1):
		taux[j-1,:] = - np.cos(np.pi*(j-1)/(ny-1))

	if tr < 0:
		tr = - tr
		taux_tr[tr:,] = taux[:ny-tr] 
		taux_tr[:tr,] = taux[0]
	else:
		taux_tr[:ny-tr,] = taux[tr:,]
		taux_tr[ny-tr:,] = taux[ny-1]

	return taux_tr

#==

def getIsothermHeight(T, THERM, grid, timeDep=True, interp=True, botVal=np.nan):
	'''Return isotherm height from given temp. data.'''

	if timeDep:
		axis = 1
	else:
		axis = 0
	
	# Get isotherm height
	Tz = np.argmin(np.abs(T-THERM), axis=axis)
	ThermZ = grid.RC.squeeze()[Tz]
	ThermZ_save = ThermZ.copy()

	if interp:

		shape = T.shape
		
		Nz, Ny, Nx = shape[-3], shape[-2], shape[-1]

		# Heights of levels above and below
		ThermZup = grid.RC.squeeze()[Tz-1]
		Tz1 = Tz+1; Tz1 = np.where(Tz1>=Nz, Tz1-1, Tz)
		ThermZdown = grid.RC.squeeze()[Tz1]

		Y = [[i]*Nx for i in range(Ny)]
		Y = [i for subX in Y for i in subX]
		X = [[i]*Ny for i in range(Nx)]
		X = [X[j][i] for i in range(Ny) for j in range(Nx)]

		if timeDep:
			Nt = shape[0]
			for ti in range(T.shape[0]):
				THERMUP = T[ti,Tz[ti].flatten()-1,Y,X].reshape((Ny,Nx))
				THERMDOWN = T[ti,Tz1[ti].flatten(),Y,X].reshape((Ny,Nx))
				THERM_Tz = T[ti,Tz[ti].flatten(),Y,X].reshape((Ny,Nx))
				ThermZ[ti] = ThermZdown[ti] + (ThermZup[ti]-ThermZdown[ti]) * (THERM - THERMDOWN) / (THERMUP-THERMDOWN)

		else:
			THERMUP = T[Tz.flatten()-1,Y,X].reshape((Ny,Nx))
			THERMDOWN = T[Tz1.flatten(),Y,X].reshape((Ny,Nx))
			THERM_Tz = T[Tz.flatten(),Y,X].reshape((Ny,Nx))
			ThermZ = ThermZdown + (ThermZup-ThermZdown) * (THERM - THERMDOWN) / (THERMUP-THERMDOWN)

		#==

	# Mask where isotherm doesn't exist
	ThermZ = np.where(ThermZ_save<grid.bathy, botVal, ThermZ)

	return ThermZ


#==

def computeBotDragQuadr0(path, grid, direction='zonal', Cd=2.5e-3):
	'''Compute quadratic bottom drag in given direction.
	0 refers to value of selectBotDragQuadr (0 in MCS, 2 in PAS).'''

	from readData import readVariable
	
	ub = readVariable('UVEL', path, file_format='nc', meta=False)
	ub = bottom(ub, grid, 'u')

	vb = readVariable('VVEL', path, file_format='nc', meta=False)
	vb = bottom(vb, grid, 'v')

	hFacCb = bottom(grid.hFacC, grid, 'h', timeDep=False)
	hFacSb = bottom(grid.hFacS, grid, 'v', timeDep=False)
	hFacWb = bottom(grid.hFacW, grid, 'u', timeDep=False)

	KE = np.zeros(ub.shape)
	KE[..., 0:-1, 0:-1] = 0.25 * (\
ub[...,0:-1,0:-1]**2 * hFacWb[0:-1,0:-1] + ub[...,0:-1,1:]**2 * hFacWb[0:-1,1:] + \
vb[...,0:-1,0:-1]**2 * hFacSb[0:-1,0:-1] + vb[...,1:,0:-1]**2 * hFacSb[1:,0:-1]) \
/ hFacCb[0:-1,0:-1]

	KE[:,:,1:] = 0.5 * (KE[:,:,1:] + KE[:,:,0:-1])

	# Now do KE at boundaries
	ub = interp(ub, 'u'); vb = interp(vb, 'v')
	KE[:,0,:] = 0.5 * (ub[:,0,:]**2 + vb[:,0,:]**2) * hFacCb[0,:]
	KE[:,-1,:] = 0.5 * (ub[:,-1,:]**2 + vb[:,-1,:]**2) * hFacCb[-1,:]
	KE[:,:,0] = 0.5 * (ub[:,:,0]**2 + vb[:,:,0]**2) * hFacCb[:,0]
	KE[:,:,-1] = 0.5 * (ub[:,:,-1]**2 + vb[:,:,-1]**2) * hFacCb[:,-1]

	drag = - Cd * np.sqrt(2*KE) * ub

	return drag

#==

def computeBotDragQuadr2(path, grid, direction='zonal', Cd=2.5e-3):
	'''Compute quadratic bottom drag in given direction.
	0 refers to value of selectBotDragQuadr (0 in MCS, 2 in PAS).'''

	from readData import readVariable
	
	ub = readVariable('UVEL', path, file_format='nc', meta=False)
	ub = bottom(ub, grid, 'u')

	vb = readVariable('VVEL', path, file_format='nc', meta=False)
	vb = bottom(vb, grid, 'v')

	hFacS = bottom(grid.hFacS, grid, 'v', timeDep=False)
	hFacC = bottom(grid.hFacC, grid, 'h', timeDep=False)

	KE = np.zeros(ub.shape)
	hFac = np.zeros(ub.shape)

	hFac[:, 0:-1, 1:] = hFacS[0:-1, 0:-1] + hFacS[0:-1, 1:] + hFacS[1:, 0:-1] + hFacS[1:, 1:]
		
	usq = ub**2
	KE[:, 0:-1, 1:] = ub[:, 0:-1, 1:]**2 + (vb[:, 0:-1,0:-1]**2 * hFacS[0:-1, 0:-1] \
+ vb[:, 0:-1, 1:]**2 * hFacS[0:-1, 1:] + vb[:, 1:, 0:-1]**2 * hFacS[1:, 0:-1] \
+ vb[:, 1:, 1:]**2 * hFacS[1:, 1:]) / hFac[:,0:-1, 1:]

	KE = np.where(hFac>0, KE, usq)

	# Now do KE at boundaries
	ub = interp(ub, 'u'); vb = interp(vb, 'v')
	KE[:,0,:] = 0.5 * (ub[:,0,:]**2 + vb[:,0,:]**2) * hFacC[0,:]
	KE[:,-1,:] = 0.5 * (ub[:,-1,:]**2 + vb[:,-1,:]**2) * hFacC[-1,:]
	KE[:,:,0] = 0.5 * (ub[:,:,0]**2 + vb[:,:,0]**2) * hFacC[:,0]
	KE[:,:,-1] = 0.5 * (ub[:,:,-1]**2 + vb[:,:,-1]**2) * hFacC[:,-1]

	drag = - Cd * np.sqrt(KE) * ub

	return drag

#==

def computeBotDragQuadr_interp(path, grid, direction='zonal', Cd=2.5e-3):
	'''Compute quadratic bottom drag in given direction.'''

	from readData import readVariable
	
	ub = readVariable('UVEL', path, file_format='nc', meta=False)
	ub = bottom(ub, grid, 'u')

	vb = readVariable('VVEL', path, file_format='nc', meta=False)
	vb = bottom(vb, grid, 'v')

	hFacCb = bottom(grid.hFacC, grid, 'h', timeDep=False)

	ub = interp(ub, 'u')
	vb = interp(vb, 'v'); 

	KE = 0.5 * (ub**2 + vb**2) * hFacCb
	drag = - Cd * np.sqrt(2*KE) * ub

	return drag

#==

def theta_f(grid, S=None, p=None, rho0=1030., g=9.81, nonlinear=True, path=None):
	'''Compute in-situ freezing point of seawater from linear EOS.'''

	if nonlinear:
		a0 = -0.0575
		a1 = 1.710523e-3
		a2 = -2.154996e-4
		b  = -7.53e-4
		c0 = 0.0
	
	else:
		a0 = -0.0575
		a1 = 0.0
		a2 = 0.0
		c0 = 0.0901
		b = -7.61e-4


	# First get pressure

	if not isinstance(p, (np.ndarray)):
		
		# Get basic state hydrostatic pressure in decibars
		if p == None:
			patm = 10.1325
			p = patm - rho0 * g * grid.RC / 10000.
			
		# Compute p from file.
		elif p == 'file':
			print('compute p from file')

	# Else, assume p provided is full p.
		
	#==

	# Next get salinity
	if not isinstance(p, (np.ndarray)):
		if S == None or S == 'file':
			print('Get S from file')
		else:
			print('Error: tools.theta_f. If S array not provided, S must be None or file.')

	# Else, assume S provided.
	
	Tf = S * (a0 + a1 * np.sqrt(S) + a2 * S) + b * p + c0
	Thf = Tf * ((patm)/p)**(2./7.)

	return Tf, Thf

#==

def smooth3(data):
	'''Apply running average over 3 time samples to smooth data.	
	Assumes that time is on zeroth axis.'''

	Nt = data.shape[0]
	data_out = data.copy()

	data_out[1:Nt-1] = (data_out[0:Nt-2] + data_out[1:Nt-1] + data_out[2:]) / 3.

	return data_out

#==

def getTSrel(grid, salttoptop=33.2, salttop=33.8, saltbot=34.7, temptop=-1.8, tempbot=1.0, tclinewid=260., tclinetop=-200, hclinewid=400, hclinetop=-60, shift=0, nztophc=4):
	'''Return T/S relaxation profiles used in MCS simulations.'''

	tclinetop += shift
	tclinebot = tclinetop - tclinewid

	hclinetop += shift
	hclinebot = hclinetop - hclinewid

	zc = grid.RC.squeeze()
	
	botweight = np.minimum(1, np.maximum(0, (zc - tclinetop) / (tclinebot - tclinetop)))
	Trel = (1 - botweight) * temptop + botweight * tempbot

	botweight = np.minimum(1, np.maximum(0, (zc - hclinetop) / (hclinebot - hclinetop)))
	Srel = (1 - botweight) * salttop + botweight * saltbot

	# Do surface values of salinity.
	dstop = (salttop - salttoptop) / nztophc
	for zi in range(1, nztophc+1):
		Srel[zi-1] = salttop - (nztophc - zi+1) * dstop

	return Trel, Srel


