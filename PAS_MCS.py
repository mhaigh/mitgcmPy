# PAS_MCS.py

# Sequence of scripts for comparison on PAS and MCS outputs.

#==========================================================

import numpy as np

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import matplotlib.pyplot as plt
import matplotlib.colors as cl

import plotting as pt
import plotting_tools as ptt

import tools

from readData import readVariable
from varDict import getPlottingVars, getTrefSref

import time

#==========================================================

# (Y,Z)-plots of Theta near continental shelf slope.
ASF = False
if ASF:

	# Window 1
	lat1 = [-75, -70]; lon1 = 248; depth1 = [0, -1000]

	# Window 2
	lat2 = [-75, -70]; lon2 = 249; depth2 = depth1
	
	# Window3
	lat3 = [-72, -64]; lon3 = 252; depth3 = depth1

	# Window4
	lat4 = [-75, -64]; lon4 = 253; depth4 = depth1

	lats = [lat1, lat2, lat3, lat4]
	lons = [lon1, lon2, lon3, lon4]
	depths = [depth1, depth2, depth3, depth4]
	
	#==

	path = '/home/michai/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	T = np.load(path+'Tmean_PAS851.npy')

	X = grid.XC#[1,:]
	Y = grid.YC#[:,1]
	Z = grid.RC.squeeze()

	vmin, vmax, cmap, title = getPlottingVars('THETA')

	#==

	# Plot window for each slice.
	for i in range(len(lats)):

		lats_ = lats[i]
		lons_ = lons[i]
		depths_ = depths[i]

		# Get indices corresponding to given lats, lons, depths.
		latsi = grid.getIndexFromLat(lats_)
		depthsi = grid.getIndexFromDepth(depths_)
		lonsi = grid.getIndexFromLon(lons_)

		# Get data in subregion.
		T1 = tools.getSubregionYZ(T[...,lonsi], latsi, depthsi)
		#T1 = tools.boundData(T1, vmin, vmax, scale=0.99999)
		Ysubr = grid.Ysubr1D(lats_); Zsubr = grid.Zsubr1D(depths_)

		# Mask and plot
		T1 = ptt.maskBathyYZ(T1, grid, xi=lonsi, subregion=True, lats=latsi, depths=depthsi)
		pt.plot1by1(T1, X=Ysubr, Y=Zsubr)#, vmin=vmin, vmax=vmax)
	
	#==

	quit()

#==

thetaHeight = False
if thetaHeight:

	PAS = True
	ANIM = True
	THERM = 0.5

	if PAS:

		if ANIM:
			path = '/data/oceans_output/shelf/pahol/mitgcm/PAS_851/run/'
			grid = Grid_PAS(path)
			X = grid.XC#[1,:]
			Y = grid.YC#[:,1]
			Z = grid.RC.squeeze()
			bathy = grid.bathy

			ts = 0; te = 10
			T = readVariable('THETA', path, file_format='nc', meta=True)[ts:te]
			TIME = T['TIME'][:]
			T = T['THETA'][:]

			zi = grid.getIndexFromDepth(-1000)

			# Get z-indices of level with Theta closest to THERM.
			Tz = np.argmin(np.abs(T[:zi]-THERM),axis=0)
			ThermZ = Z[Tz]
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)
			ThermZ = np.where(ThermZ<bathy, np.nan, ThermZ)

			vmin = -800; vmax = -200
			#vmin = -500; vmax=-200		
			title = 'PAS851 ' + str(THERM) + ' deg. isotherm height'
			text_data = ptt.getTextData(TIME, 'ctime', X[1,1], Y[1,1])
			
		else:
			path = '/home/michai/Documents/data/PAS_851/run/'
			grid = Grid_PAS(path)
			T = np.load(path+'Tmean_PAS851.npy')

			X = grid.XC#[1,:]
			Y = grid.YC#[:,1]
			Z = grid.RC.squeeze()
			bathy = grid.bathy

			zi = grid.getIndexFromDepth(-1000)

			# Get z-indices of level with Theta closest to THERM.
			Tz = np.argmin(np.abs(T[:zi]-THERM),axis=1)
			ThermZ = Z[Tz]
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)
			ThermZ = np.where(ThermZ<bathy, np.nan, ThermZ)

			vmin = -800; vmax = -200
			#vmin = -500; vmax=-200		
			title = 'PAS851 mean ' + str(THERM) + ' deg. isotherm height'

		SUBR = True
		if SUBR:
			#lats = [-75, -71]; lons = [240, 270]
			lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

			ThermZ = tools.getSubregionXY(ThermZ, latsi, lonsi)
			bathy = tools.getSubregionXY(bathy, latsi, lonsi)
			X = grid.Xsubr1D(lons)
			Y = grid.Ysubr1D(lats)

	#==

	else:
		path = '/home/michai/Documents/data/MCS_038/run/'
		grid = Grid(path)
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()
		bathy = grid.bathy

		if ANIM:
			T = readVariable('THETA', path, file_format='nc', meta=True)
			TIME = T['TIME'][:]
			T = T['THETA'][:]

			# Get z-indices of level with Theta closest to THERM.
			Tz = np.argmin(np.abs(T-THERM),axis=1)
			ThermZ = Z[Tz]
			ThermZ = np.where(ThermZ<bathy, np.nan, ThermZ)
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)
			text_data = ptt.getTextData(TIME, 'month', X[1], Y[1])
			title = 'MCS038 ' + str(THERM) + ' deg. isotherm height'

		else:
			T = np.mean(readVariable('THETA', path, file_format='nc', meta=False)[60:], axis=0)
			# Get z-indices of level with Theta closest to THERM.
			Tz = np.argmin(np.abs(T-THERM),axis=0)
			ThermZ = Z[Tz]
			ThermZ = np.where(ThermZ<bathy, np.nan, ThermZ)
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)
			title = 'MCS038 mean ' + str(THERM) + ' deg. isotherm height'

		vmin = -550; vmax = -400
		#vmin = None; vmax = None

		#==

	#==

	xlabel = 'LON (km)'; ylabel = 'LAT (km)'
	cmap = 'jet'
	
	if ANIM:
		pt.animate1by1(ThermZ, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=True, outname='animate1by1Surf.mp4', text_data=text_data)
	
	else:
		pt.plot1by1(ThermZ, X=X, Y=Y, mesh=True, title=title, vmin=vmin, vmax=vmax, cmap='jet'); quit()#

	quit()
#==

quiver = False
if quiver:
	
	path = '/home/michai/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	u = np.load(path+'umean_PAS851.npy')
	v = np.load(path+'vmean_PAS851.npy')

	X = grid.XC[1,:]
	Y = grid.YC[:,1]
	Z = grid.RC.squeeze()

	#level = 0
	depth = -250; level = grid.getIndexFromDepth(depth)
	print(level)
		
	#==

	# Select level and interpolate to cell centres.
	u = u[level]; v = v[level]
	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

	u = ptt.maskBathyXY(u, grid, level, timeDep=False)
	v = ptt.maskBathyXY(v, grid, level, timeDep=False)

	#contour = np.load(path+'Tmean_PAS851.npy')[level; vmin = -1.; vmax = 1 
	contour = grid.bathy; vmin=-1000; vmax=-200
	contour = ptt.maskBathyXY(contour, grid, level, timeDep=False)
	
	title = '(u, v) & bathymetry; Z = ' + str(Z[level]) + ' m'
	xlabel = 'LON'; ylabel = 'LAT'

	SUBR = True
	if SUBR:
		#lats = [-76, -68]; lons = [245, 260]#[230, 270]#
		lats = [-75, -70]; lons = [245, 260]		
		latsi = grid.getIndexFromLat(lats)
		lonsi = grid.getIndexFromLon(lons)

		u = tools.getSubregionXY(u, latsi, lonsi)
		v = tools.getSubregionXY(v, latsi, lonsi)
		contour = tools.getSubregionXY(contour, latsi, lonsi)

		X = grid.Xsubr1D(lons)
		Y = grid.Ysubr1D(lats)


	# Sample rate & level
	d = 2
	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]
	
	#
	
	pt.quiver1by1(u, v, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, title=title, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel)

	quit()

#==
	
btpcStr = False
if btpcStr:

	PAS = True

	if PAS:
		path = '/home/michai/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		u = np.load(path+'umean_PAS851.npy')

		X = grid.XC[1,:]
		Y = grid.YC[:,1]
		Z = grid.RC.squeeze()

		u = ptt.maskBathyAll(u, grid, timeDep=False)
		T = tools.barotropicStreamfunction(u, grid, timeDep=False)
		bathy = grid.bathy
		levels = [-900, -700, -500]

		SUBR = True
		if SUBR:
			lats = [-76, -71]; lons = [245, 260]#[230, 270]#
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			T = tools.getSubregionXY(T, latsi, lonsi); 	bathy = tools.getSubregionXY(bathy, latsi, lonsi);
			X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)


		title = 'PAS851 barotropic streamfunction'

	else:

		path = '/home/michai/Documents/data/MCS_038/run/'
		grid = Grid(path)

		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		u = readVariable('UVEL', path, file_format='nc', meta=False)[-1]
		u = ptt.maskBathyAll(u, grid, timeDep=False)
		T = tools.barotropicStreamfunction(u, grid, timeDep=False)
		bathy = grid.bathy
		levels = [-900, -800, -700, -600, -501]

		
		SUBR = True
		if SUBR:
			lats = [0, 300e3]; lons = [0,600e3]
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			T = tools.getSubregionXY(T, latsi, lonsi); bathy = tools.getSubregionXY(bathy, latsi, lonsi);
			X = grid.Xsubr1D(lons)/1000; Y = grid.Ysubr1D(lats)/1000
		
		title = 'MCS038 barotropic streamfunction'

	vmin = -2; vmax = 2
	cmap = 'turbo'#'gist_rainbow'

	pt.plot1by1(T, X=X, Y=Y, mesh=True, cmap=cmap, contour=bathy, contourlevels=levels, title=title)#, vmin=vmin, vmax=vmax)
	#pt.plot1by2([umean[-1], umean2[-1]], mesh=True, vmin=vmin, vmax=vmax)
	quit()

#==

# Time-mean difference between surface and bottom flow.
brclnc = False
if brclnc:

	PAS = True

	if PAS:
		path = '/home/michai/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		u = np.load(path+'umean_PAS851.npy')
		v = np.load(path+'vmean_PAS851.npy')

		X = grid.XC[1,:]
		Y = grid.YC[:,1]
		Z = grid.RC.squeeze()

		# Transfer to t-points.
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Get bottom val.
		ub = tools.bottom(u, grid, 'u', timeDep=False)
		vb = tools.bottom(v, grid, 'v', timeDep=False)

		# Get 'baroclinicity'
		us = u[0] - ub
		vs = v[0] - vb

		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)
		vmin = -1000; vmax = -500

		title = 'PAS851 (u_t-u_b, v_t-v_b) & bathy'
		xlabel = 'LON (deg)'; ylabel = 'LAT (deg)'

		# Sample rate
		d = 6

		SUBR = True
		if SUBR:
			lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			us = tools.getSubregionXY(us, latsi, lonsi); vs = tools.getSubregionXY(vs, latsi, lonsi)
			contour = tools.getSubregionXY(contour, latsi, lonsi);
			X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)

	else:

		path = '/home/michai/Documents/data/MCS_038/run/'
		grid = Grid(path)
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, file_format='nc', meta=False)[60:]
		v = readVariable('VVEL', path, file_format='nc', meta=False)[60:]
		u = np.mean(u, axis=0); v = np.mean(v, axis=0)
	
		# Transfer to t-points.
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Get bottom val.
		ub = tools.bottom(u, grid, 'u', timeDep=False)
		vb = tools.bottom(v, grid, 'v', timeDep=False)

		# Get 'baroclinicity'
		us = u[0] - ub
		vs = v[0] - vb

		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)
		vmin = -1000; vmax = -500

		title = 'MCS038 (u_t-u_b, v_t-v_b) & bathy'
		xlabel = 'LON (km)'; ylabel = 'LAT (km)'

		# Sample rate
		d = 6

	#==

	# PLOT
 
	us = us[..., ::d, ::d]; vs = vs[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]

	pt.quiver1by1(us, vs, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, title=title, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel)

#==

# Animate velocity vectors and temperature. Each frame is different level.
animateUVTdepth = True
if animateUVTdepth:

	PAS = False
	
	if PAS:
		path = '/home/michai/Documents/data/PAS_8512/run/'
		grid = Grid_PAS(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		#vmin = -800; vmax = -100
		vmin = -600; vmax = -200
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = np.load(path+'umean_PAS851.npy')
		v = np.load(path+'vmean_PAS851.npy')
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Load temperature
		T = np.load(path+'Tmean_PAS851.npy')
		cvmin, cvmax, ccmap, title = getPlottingVars('THETA')
		cmvmin = -2.; cvmax = 2.
		T = tools.boundData(T, cvmin, cvmax, 0.9999)
		# Sample rate
		d = 4

	#==

	else:

		#ts = 60
		ts = 17; te = 18

		path = '/home/michai/Documents/data/MCS_038/run/'
		grid = Grid(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		vmin = -1000; vmax = -500
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, file_format='nc', meta=False)[ts:te]
		v = readVariable('VVEL', path, file_format='nc', meta=False)[ts:te]
		u = np.mean(u, axis=0); v = np.mean(v, axis=0)
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Load temperature
		T = readVariable('THETA', path, file_format='nc', meta=False)[ts:te]
		T = np.mean(T, axis=0)
		cvmin, cvmax, ccmap, title = getPlottingVars('THETA')

		# Sample rate
		d = 6

	for zi in range(len(Z)):
		u[zi] = ptt.maskBathyXY(u[zi], grid, zi, timeDep=False)
		v[zi] = ptt.maskBathyXY(v[zi], grid, zi, timeDep=False)
		T[zi] = ptt.maskBathyXY(T[zi], grid, zi, timeDep=False)

	SUBR = True
	if PAS and SUBR:
		lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
		latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
		u = tools.getSubregionXY(u, latsi, lonsi); v = tools.getSubregionXY(v, latsi, lonsi)
		T = tools.getSubregionXY(T, latsi, lonsi); contour = tools.getSubregionXY(contour, latsi, lonsi); 
		X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)

	title = '(u, v); T; bathy'
	text = ['Z = ' + str(z) + ' m' for z in Z]
	outpath = '/home/michai/Documents/Python/mitgcmPy/images/testingMov/'

	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]; Td = T[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]
	u[:,-1,-1] = 10.; v[:,-1,-1] = 10.;
	Td[:,-1,-1] = cvmin; Td[:,-1,-2] = cvmax
	#u = tools.boundData(u, -0.2, 0.2, 0.999); v = tools.boundData(v, -0.2, 0.2, 0.999)

	#for zi in range(len(Z)):
	#	outname = f'{zi:03}'
	#	text_data = {'text':'Z = ' + str(Z[zi]) + ' m', 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}
		
		#pt.quiver1by1(u[zi], v[zi], Xd, Yd, C=Td[zi], ccmap=ccmap, cvmin=cvmin, cvmax=cvmax, contour=contour, X=X, Y=Y, contourf=False, vmin=vmin, vmax=vmax, text_data=text_data, title=title, show=False, save=True, outname=outname, outpath=outpath, figsize=(7,4))
	
	# Animating in normal way doesn't work.	
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}
	pt.animate1by1quiver(u, v, Xd, Yd, C=Td, contour=contour, X=X, Y=Y, cmap='plasma', contourf=False, vmin=vmin, vmax=vmax, text_data=text_data, title=title, figsize=(7,4), dpi=300)

	quit()

#==

T_transport = False
if T_transport:


	PAS = False
	
	if PAS:
		path = '/home/michai/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		vmin = -0.5; vmax = 0.5
		xlabel = 'LAT (deg)'; ylabel = 'LON (deg)';
		X = grid.XC[1,:]
		Y = grid.YC[:,1]
		Z = grid.RC.squeeze()
		
		dxc = tools.ddx(grid.XC, 1)
		dyc = tools.ddy(grid.YC, 1)
		dzc = tools.ddz(grid.RC, 1)	

		# Load time-mean velocities
		uT = np.load(path+'uT_PAS851.npy')
		vT = np.load(path+'vT_PAS851.npy')
	
		uTx = tools.ddx(uT, dxc)
		vTy = tools.ddy(vT, dyc)
	
		Tconv = - uTx - vTy

		latsi = None; lonsi = None
		SUBR = True
		if PAS and SUBR:
			lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			Tconv = tools.getSubregionXY(Tconv, latsi, lonsi)
			X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)
		
		Tconv_lat_= [0 for i in range(len(Z))]
		for zi in range(len(Z)):
			Tconvzi = tools.boundData(Tconv[zi], vmin, vmax, 0.9999)
			Tconv_lat_[zi] = ptt.maskBathyXY(Tconvzi, grid, zi, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	else:
		
		path = '/home/michai/Documents/data/MCS_038/run/'
		grid = Grid(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		vmin = -5.e-8; vmax = -vmin
		xlabel = 'LAT (km)'; ylabel = 'LON (km)';
		scale = 1.e3
		X = grid.XC[1,:]/scale
		Y = grid.YC[:,1]/scale
		Z = grid.RC.squeeze()
	
		dxc = scale*(X[1]-X[0]); dyc = scale*(Y[1]-Y[0]); dzc = Z[1]-Z[0]

		ts = -5; te = -1
		#ts = 0; te = 10

		# Load variables and define fluxes
		T = readVariable('THETA', path, file_format='nc', meta=False)[ts:te]

		u = readVariable('UVEL', path, file_format='nc', meta=False)[ts:te]
		u = tools.interp(u, 'u'); uT = np.mean(u*T, axis=0);

		u = readVariable('VVEL', path, file_format='nc', meta=False)[ts:te]
		u = tools.interp(u, 'v'); vT = np.mean(u*T, axis=0);

		u = readVariable('WVEL', path, file_format='nc', meta=False)[ts:te]
		u = tools.interp(u, 'w'); wT = np.mean(u*T, axis=0);

		# Compute convergences.		
		Tconv_lat = - tools.ddx(uT, dxc) - tools.ddy(vT, dyc)
		Tconv_vert = - tools.ddy(wT, dzc)

		Tconv_lat_ = [0 for i in range(len(Z))]
		Tconv_vert_ = [0 for i in range(len(Z))]
		for zi in range(len(Z)):
			Tconv_tmp = tools.boundData(Tconv_lat[zi], vmin, vmax, 0.9999)
			Tconv_lat_[zi] = ptt.maskBathyXY(Tconv_tmp, grid, zi, timeDep=False)
			Tconv_tmp = tools.boundData(Tconv_vert[zi], vmin, vmax, 0.9999)
			Tconv_vert_[zi] = ptt.maskBathyXY(Tconv_tmp, grid, zi, timeDep=False)

	#== 

	# PLOT

	cmap = 'coolwarm'
	title = 'Theta flux convergence'
	text = ['Z = ' + str(z) + ' m' for z in Z]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

	pt.animate1by1(Tconv_lat_, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, figsize=(6,4), outname='Tconv_lat.mp4')

	pt.animate1by1(Tconv_vert_, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, figsize=(6,4), outname='Tconv_vert.mp4')

#==

PAS_WIND = False
if PAS_WIND:

	path = '/home/michai/Documents/data/PAS_8512/run/'
	grid = Grid_PAS(path)
	contour = grid.bathy
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

	vmin = -1000; vmax = -100
	xlabel = 'LAT (deg)'; ylabel = 'LON (deg)';
	X = grid.XC[1,:]
	Y = grid.YC[:,1]
	Z = grid.RC.squeeze()
	
	dxc = tools.ddx(grid.XC, 1)
	dyc = tools.ddy(grid.YC, 1)
	dzc = tools.ddz(grid.RC, 1)	

	# Load time-mean velocities
	uw = np.load(path+'uwindmean_PAS851.npy')
	vw = np.load(path+'vwindmean_PAS851.npy')


	SUBR = True
	if SUBR:
		lats = [-76, -68]; lons = [245, 260]#[230, 270]#
		latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
		uw = tools.getSubregionXY(uw, latsi, lonsi); vw = tools.getSubregionXY(vw, latsi, lonsi)
		contour = tools.getSubregionXY(contour, latsi, lonsi); 
		X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)

	title = '(uwind, vwind); bathy'
	text = ['Z = ' + str(z) + ' m' for z in Z]
	outpath = '/home/michai/Documents/Python/mitgcmPy/images/testingMov/'

	d = 6
	uw = uw[..., ::d, ::d]; vw = vw[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]

	pt.quiver1by1(uw, vw, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, vmin=vmin, vmax=vmax, title=title)


