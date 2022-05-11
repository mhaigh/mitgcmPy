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

FORMSTRESS = True
if FORMSTRESS:

	# account for SSH in new Pb calc.

	path = '/home/michai/Documents/data/PAS_8512/run/'
	grid = Grid_PAS(path)
	Pb = np.load(path+'Pb_ETAN_mean_PAS851.npy')

	lats = [-76, -71.5]; lons = [245, 260]#[230, 270]#
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
	X = grid.Xsubr1D(lons)
	Y = grid.Ysubr1D(lats)

	Pb = tools.getSubregionXY(Pb, latsi, lonsi)
	H = - grid.bathy
	Hx = tools.ddx(H, grid.DXG)
	Hy = tools.ddy(H, grid.DYG)
	Hx = tools.getSubregionXY(Hx, latsi, lonsi)

	Pb = ptt.maskBathyXY(Pb, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	Hx = ptt.maskBathyXY(Hx, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	vmin = 0.2e7; vmax = 1.e7
	Pb = tools.boundData(Pb, vmin, vmax, scale=0.9999)

	pt.plot1by1(Pb, X=X, Y=Y, vmin=vmin, vmax=vmax, mesh=True, title='Bottom pressure')
	pt.plot1by1(Hx, X=X, Y=Y, mesh=True, title='Hx', vmin=-0.01, vmax=0.01)
	#pt.plot1by1(H, mesh=True, title='H')
	FS = Hx*Pb
	print(np.sum(FS))
	vmin = -2e5; vmax = -vmin
	pt.plot1by1(FS, X=X, Y=Y, vmin=vmin, vmax=vmax, mesh=True, title='TFS')

	quit()
	
	
#==

Bathy = False
if Bathy:

	path = '/home/michai/Documents/data/PAS_851/run/'
	#path = 	'/home/michai/Documents/data/PISOMIP_001/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy

	X = grid.XC#[1,:]
	Y = grid.YC#[:,1]
	print(Y)

	SUBR = True
	if SUBR:
		lats = [-76, -70.5]; lons = [230, 260]#[230, 270]#
		latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)
		X = grid.Xsubr1D(lons)
		Y = grid.Ysubr1D(lats)
	
		bathy = ptt.maskBathyXY(bathy, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	# PLOT 
	
	#bathy = ptt.maskBathyXY(bathy, grid, 0, timeDep=False, subregion=SUBR, lats=latsi, lons=lonsi)

	pt.plot1by1(1./bathy, X=X, Y=Y, cmap='plasma', mesh=False, vmin=-1./200, vmax=-1./800, contourfNlevels=21, figsize=(5,3), title='1/bathymetry')

	bathy = tools.boundData(bathy, vmin=-899, vmax=-201)

	pt.plot1by1(bathy, X=X, Y=Y, cmap='plasma', mesh=False, vmin=-900, vmax=-200, contourfNlevels=21, xlabel='LON', ylabel='LAT', figsize=(5,3), title='bathymetry')

	quit()

#==

# (Y,Z)-plots of Theta near continental shelf slope.
ASF = False
if ASF:

	# Window 1
	lat1 = [-75, -69]; lon1 = 233; depth1 = [0, -1000]

	# Window 2
	lat2 = [-75, -69]; lon2 = 240; depth2 = depth1
	
	# Window3
	lat3 = [-75, -69]; lon3 = 245; depth3 = depth1

	# Window4
	lat4 = [-75, -69]; lon4 = 250; depth4 = depth1

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

	##plt.plot(T[:20,180,20], Z[:20]); plt.xlim(33.5, 35.); plt.grid(); plt.show();

	#T0 = ptt.maskBathyXY(T[0], grid, 0, timeDep=False)
	#T0 = ptt.maskDraftXY(T0, grid, 0, timeDep=False)
	#fig, ax = plt.subplots()
	#plt.subplot(111); plt.gca().patch.set_color('.25')
	#plt.pcolor(T0, vmin=33.4, vmax=34., cmap='plasma'); plt.colorbar(); plt.show(); quit()

	vmin, vmax, cmap, title = getPlottingVars('THETA')
	xlabel = 'LAT (deg.)'; ylabel = 'DEPTH (m)'
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

		title = 'PAS Theta (lon = ' + str(lons_) + ' deg.)'
		# Mask and plot
		T1 = ptt.maskBathyYZ(T1, grid, xi=lonsi, subregion=True, lats=latsi, depths=depthsi)
		T1 = ptt.maskDraftYZ(T1, grid, xi=lonsi, subregion=True, lats=latsi, depths=depthsi)
		pt.plot1by1(T1, X=Ysubr, Y=Zsubr, title=title, xlabel=xlabel, ylabel=ylabel, mesh=False)#, vmin=vmin, vmax=vmax)
	
	#==

	quit()

#==

thetaHeight = False
if thetaHeight:

	PAS = False
	ANIM = True
	THERM = 0.5

	if PAS:

		if ANIM:
			path = '/home/michai/Documents/data/PAS_8512/run/'
			grid = Grid_PAS(path)
			bathy = grid.bathy

			ThermZ = np.load(path+'Therm05Z_PAS851.npy')
		
			# Subregion used in bsl script.
			lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True, subregion=True, lats=latsi, lons=lonsi)
			ThermZ = np.where(ThermZ<bathy[latsi[0]:latsi[1]+1, lonsi[0]:lonsi[1]+1], np.nan, ThermZ)

			X = grid.Xsubr1D(lons)
			Y = grid.Ysubr1D(lats)

			vmin = -800; vmax = -200
			#vmin = -500; vmax=-200		
			title = 'PAS851 ' + str(THERM) + ' deg. isotherm height'
			
			text_data = None#ptt.getTextData(TIME, 'ctime', X[1,1], Y[1,1])
			
		else:
			path = '/home/michai/Documents/data/PAS_851/run/'
			grid = Grid_PAS(path)
			T = np.load(path+'Tmean_PAS851.npy')

			X = grid.XC#[1,:]
			Y = grid.YC#[:,1]
			Z = grid.RC.squeeze()
			bathy = grid.bathy

			print(bathy.shape)
			print(T.shape)
			zi = grid.getIndexFromDepth(-1000)

			# Get z-indices of level with Theta closest to THERM.
			Tz = np.argmin(np.abs(T[:zi]-THERM),axis=0)
			ThermZ = Z[Tz]
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)
			ThermZ = np.where(ThermZ<bathy, np.nan, ThermZ)

			vmin = -800; vmax = -200
			#vmin = -500; vmax=-200		
			title = 'PAS mean ' + str(THERM) + ' deg. isotherm height'
			xlabel = 'LON'; ylabel = 'LAT'

			SUBR = True
			if SUBR:
				#lats = [-75, -71]; lons = [240, 270]
				lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
				latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

				ThermZ = tools.getSubregionXY(ThermZ, latsi, lonsi)
				bathy = tools.getSubregionXY(bathy, latsi, lonsi)
				X = grid.Xsubr1D(lons)
				Y = grid.Ysubr1D(lats)

	else:


		exp = 'MCS_111'
		path = '/home/michai/Documents/data/'+exp+'/run/'
		grid = Grid(path)
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()
		bathy = grid.bathy

		if ANIM:

			ts = 0
			T = readVariable('THETA', path, file_format='nc', meta=True)
			TIME = T['TIME'][ts:]
			T = T['THETA'][ts:]

			# Get z-indices of level with Theta closest to THERM.
			ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True)
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)

			text_data = ptt.getTextData(TIME, 'month', X[1], Y[1])


		else:
			T = np.mean(readVariable('THETA', path, file_format='nc', meta=False)[58:], axis=0)
			# Get z-indices of level with Theta closest to THERM.
			Tz = np.argmin(np.abs(T-THERM),axis=0)
			ThermZ = Z[Tz]
			ThermZ = np.where(ThermZ<bathy, np.nan, ThermZ)
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)

		title = exp + ' ' + str(THERM) + ' deg. isotherm height'
	
		vmin = -600; vmax = -300	
		#vmin = -1000; vmax = -400
		#vmin = None; vmax = None
		xlabel = 'LON (km)'; ylabel = 'LAT (km)'

		#==

	#==


	cmap = 'jet'
	
	if ANIM:
		pt.animate1by1(ThermZ, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=True, outname='animate1by1Surf.mp4', text_data=text_data)
	
	else:
		pt.plot1by1(ThermZ, X=X, Y=Y, mesh=True, title=title, vmin=vmin, vmax=vmax, cmap='jet', xlabel=xlabel, ylabel=ylabel); quit()#

	quit()
#==

quiver = False
if quiver:
	
	PAS = False

	level = 0
	#depth = -350; level = grid.getIndexFromDepth(depth)
	print(level)


	if PAS:
		path = '/home/michai/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		u = np.load(path+'umean_PAS851.npy')
		v = np.load(path+'vmean_PAS851.npy')

		X = grid.XC[1,:]
		Y = grid.YC[:,1]
		Z = grid.RC.squeeze()

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

	else:

		#path = '/home/michai/Documents/data/PAS_666/run/'
		path = '/home/michai/Documents/data/PISOMIP_004/run/'
		ts = -1
		
		grid = Grid(path)
		bathy = grid.bathy

		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		xlabel = 'LON (km)'; ylabel = 'LAT (km)'

		VAR = 'ETAN'
		vmin, vmax, cmap, title = getPlottingVars(VAR)

		contour = readVariable(VAR, path, file_format='nc', meta=False)[ts]
		u = readVariable('UVEL', path, file_format='nc', meta=False)[ts]
		v = readVariable('VVEL', path, file_format='nc', meta=False)[ts]
		
	#==

	# Select level and interpolate to cell centres.
	u = u[level]; v = v[level]
	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

	u = ptt.maskBathyXY(u, grid, level, timeDep=False)
	v = ptt.maskBathyXY(v, grid, level, timeDep=False)

	#contour = np.load(path+'Tmean_PAS851.npy')[level; vmin = -1.; vmax = 1 
	contour = ptt.maskBathyXY(contour, grid, level, timeDep=False)
	
	title = '(u, v) & SSH'

	# Sample rate & level
	d = 6
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
			lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			T = tools.getSubregionXY(T, latsi, lonsi); 	bathy = tools.getSubregionXY(bathy, latsi, lonsi);
			X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)

			pt.plot1by1(tools.boundData(bathy,-900, -100, 0.9999), X=X, Y=Y, cmap='jet', vmin=-900, vmax=-100, mesh=True, title='PAS bathy'); quit()

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
animateUVTdepth = False
if animateUVTdepth:

	PAS = False
	SUBR = True
	BOT = False

	if PAS:
		path = '/home/michai/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		contour = grid.bathy; vmin = -800; vmax = -100; ctitle = 'bathy'
		#contour = grid.draft; vmin = -600; vmax = -0; ctitle = 'Ice shelf draft'

		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		X = grid.XC[1,:]
		Y = grid.YC[:,1]
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = np.load(path+'umean_PAS851.npy')
		v = np.load(path+'vmean_PAS851.npy')
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		if BOT:
			ubot = tools.bottom(u, grid, 'h', timeDep=False)
			vbot = tools.bottom(v, grid, 'h', timeDep=False)

		COLOUR = 'T'
		if COLOUR == 'T': 
			# Load temperature
			T = np.load(path+'Tmean_PAS851.npy')
			cvmin, cvmax, ccmap, title = getPlottingVars('THETA')
			cmvmin = -2.; cvmax = 2.
		elif COLOUR == 'S':
			# For salt to work, use PAS_8512, SUBR=True, comment out maksing procedure and
			# subregioning of u,v,T.
			T = np.load(path+'Smean_PAS851.npy')
			cvmin, cvmax, ccmap, title = getPlottingVars('SALT')
			#cmvmin = -2.; cvmax = 2.

		T = tools.boundData(T, cvmin, cvmax, 0.9999)
		# Sample rate
		d = 4

	#==

	else:

		#ts = 20
		ts = 60; te = 120

		path = '/home/michai/Documents/data/MCS_110/run/'
		grid = Grid(path)
		contour = grid.bathy; ctitle = 'bathy'
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		vmin = -1000; vmax = -400
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, file_format='nc', meta=False)[ts:]
		v = readVariable('VVEL', path, file_format='nc', meta=False)[ts:]
		u = np.mean(u, axis=0); v = np.mean(v, axis=0)
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Load temperature
		T = readVariable('THETA', path, file_format='nc', meta=False)[ts:]
		T = np.mean(T, axis=0)
		cvmin, cvmax, ccmap, title = getPlottingVars('THETA')

		# Sample rate
		d = 6

	for zi in range(len(Z)):
		u[zi] = ptt.maskBathyXY(u[zi], grid, zi, timeDep=False)
		v[zi] = ptt.maskBathyXY(v[zi], grid, zi, timeDep=False)
		T[zi] = ptt.maskBathyXY(T[zi], grid, zi, timeDep=False)


	if PAS and SUBR:
		lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
		latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
		u = tools.getSubregionXY(u, latsi, lonsi); v = tools.getSubregionXY(v, latsi, lonsi)
		T = tools.getSubregionXY(T, latsi, lonsi)
		contour = tools.getSubregionXY(contour, latsi, lonsi); 
		X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)
		if BOT:
			ubot = tools.getSubregionXY(ubot, latsi, lonsi)
			vbot = tools.getSubregionXY(vbot, latsi, lonsi)

	title = '(u, v); S; ' + ctitle
	#title = '(u, v); T; bathy'
	text = ['Z = ' + str(z) + ' m' for z in Z]
	outpath = '/home/michai/Documents/Python/mitgcmPy/images/testingMov/'

	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]; Td = T[..., ::d, ::d]

	Xd = X[::d]; Yd = Y[::d]
	#u[:,-1,-1] = 10.; v[:,-1,-1] = 10.;
	Td[:,-1,-1] = cvmin; Td[:,-1,-2] = cvmax
	u = tools.boundData(u, -0.1, 0.1, 0.999); v = tools.boundData(v, -0.1, 0.1, 0.999)

	#for zi in range(len(Z)):
	#	outname = f'{zi:03}'
	#	text_data = {'text':'Z = ' + str(Z[zi]) + ' m', 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}
		
	if BOT:
		ubot = ubot[..., ::d, ::d]; vbot = vbot[..., ::d, ::d]
		pt.quiver1by1(ubot, vbot, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, vmin=vmin, vmax=vmax, show=True, figsize=(7,4))
		quit()
	
	cmap = 'YlGn'#'plasma'
	# Animating in normal way doesn't work.	
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}
	pt.animate1by1quiver(u, v, Xd, Yd, C=Td, contour=contour, X=X, Y=Y, cmap=cmap, contourf=False, vmin=vmin, vmax=vmax, text_data=text_data, title=title, figsize=(7,4), dpi=300)

	quit()

#==

# Animate velocity vectors and temperature at fixed level.
animateUVT = False
if animateUVT:

	PAS = False
	
	if PAS:

		path = '/home/michai/Documents/data/PAS_8512/run/'
		grid = Grid_PAS(path)
		contour = grid.bathy; vmin = -800; vmax = -100; ctitle = 'bathy'
		#contour = grid.draft; vmin = -600; vmax = -0; ctitle = 'Ice shelf draft'
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)


		vmin = -800; vmax = -100

		level = 16
		#print(grid.RC.squeeze()[level]); quit()
		z = '455' # '257', '355'
		
		ts = 107; te = 502
	
		# Load time-mean velocities
		u = np.load(path+'u'+z+'_PAS851.npy')[ts:te]
		v = np.load(path+'v'+z+'_PAS851.npy')[ts:te]
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		THETA = False
		if THETA:
			# Load temperature
			T = np.load(path+'T'+z+'_PAS851.npy')[ts:te]
			cvmin, cvmax, ccmap, title = getPlottingVars('THETA')
			cmvmin = -2.; cvmax = 2.
			title = '(u, v); T; ' + ctitle + ' Z = ' + str(grid.RC.squeeze()[level]) + ' m'
		else:
			# Load salinity
			T = np.load(path+'S'+z+'_PAS851.npy')[ts:te]
			cvmin, cvmax, ccmap, title = getPlottingVars('SALT')
			cmvmin = -34.2; cvmax = 34.6
			title = '(u, v); S; ' + ctitle + ' Z = ' + str(grid.RC.squeeze()[level]) + ' m'

		T = tools.boundData(T, cvmin, cvmax, 0.9999)

		# Sample rate
		d = 4

		lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
		latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
		X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)
		contour = tools.getSubregionXY(contour, latsi, lonsi)

		TIME = np.load(path+'TIME_PAS851.npy')[ts:te]
		text_data = ptt.getTextData(TIME, 'ctime', X[1], Y[1])

		for ti in range(T.shape[0]):
			u[ti] = ptt.maskBathyXY(u[ti], grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
			v[ti] = ptt.maskBathyXY(v[ti], grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
			T[ti] = ptt.maskBathyXY(T[ti], grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)


	# IF MCS
	else:

		ts = 0; te = -1

		path = '/home/michai/Documents/data/MCS_113/run/'
		grid = Grid(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		depth = -490; level = grid.getIndexFromDepth(depth)

		vmin = -1000; vmax = -300
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, file_format='nc', meta=False)[ts:, level]
		print(u.shape)
		v = readVariable('VVEL', path, file_format='nc', meta=False)[ts:, level]
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Load temperature
		T = readVariable('THETA', path, file_format='nc', meta=False)[ts:, level]
		cvmin, cvmax, ccmap, title = getPlottingVars('THETA')
		title = '(u, v); T; bathy; Z = ' + str(grid.RC.squeeze()[level]) + ' m'

		text = ['month ' + str(ti) for ti in range(T.shape[0])]
		text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

		# Sample rate
		d = 6

		for ti in range(T.shape[0]):
			u[ti] = ptt.maskBathyXY(u[ti], grid, level, timeDep=False)
			v[ti] = ptt.maskBathyXY(v[ti], grid, level, timeDep=False)
			T[ti] = ptt.maskBathyXY(T[ti], grid, level, timeDep=False)
			
	#==

	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]; Td = T[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]
	#u[:,-1,-1] = 10.; v[:,-1,-1] = 10.;
	Td[:,-1,-1] = cvmin; Td[:,-1,-2] = cvmax
	u = tools.boundData(u, -0.1, 0.1, 0.999); v = tools.boundData(v, -0.1, 0.1, 0.999)
	
	cmap = 'YlGn' #'plasma'


	#plt.pcolor(X, Y, T[0], vmin=33.7, vmax=34.3, cmap='plasma'); plt.colorbar();
	#plt.quiver(Xd, Yd, u[0], v[0]) 
	#plt.show(); quit()

	pt.animate1by1quiver(u, v, Xd, Yd, C=Td, contour=contour, X=X, Y=Y, cmap=cmap, contourf=False, vmin=vmin, vmax=vmax, text_data=text_data, title=title, figsize=(7,4), dpi=300)

	quit()

#==

# Animate velocity shear between two levels
animateUVshear = False
if animateUVshear:

	PAS = False
	
	if PAS:

		path = '/home/michai/Documents/data/PAS_8512/run/'
	
	else:

		ts = 0; te = 112

		path = '/home/michai/Documents/data/MCS_113/run/'
		grid = Grid(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		depth1 = -10; level1 = grid.getIndexFromDepth(depth1)
		depth2 = -490; level2 = grid.getIndexFromDepth(depth2)

		vmin = -1000; vmax = -300
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, file_format='nc', meta=False)[ts:]
		v = readVariable('VVEL', path, file_format='nc', meta=False)[ts:]

		u = u[:,level2] - u[:,level1]
		v = v[:,level2] - v[:,level1]

		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		title = 'u' + str(-depth2) + ' - u' + str(-depth1) + '; bathy'

		text = ['month ' + str(ti) for ti in range(u.shape[0])]
		text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

		# Sample rate
		d = 6

		qlim = 0.2
		u = tools.boundData(u, -qlim, qlim, 0.999); v = tools.boundData(v, -qlim, qlim, 0.999)

		for ti in range(u.shape[0]):
			u[ti] = ptt.maskBathyXY(u[ti], grid, level2, timeDep=False)
			v[ti] = ptt.maskBathyXY(v[ti], grid, level2, timeDep=False)
			
	#==

	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]

	
	cmap = 'YlGn' #'plasma'

	pt.animate1by1quiver(u, v, Xd, Yd, qlim=qlim, contour=contour, X=X, Y=Y, cmap=cmap, contourf=False, vmin=vmin, vmax=vmax, text_data=text_data, title=title, figsize=(7,4), dpi=300)

	quit()

#==

# Animate velocity vectors and temperature at fixed level.
PAS_uvst_season = False
if PAS_uvst_season:


	path = '/home/michai/Documents/data/PAS_8512/run/'
	grid = Grid_PAS(path)
	contour = grid.bathy; vmin = -800; vmax = -100; ctitle = 'bathy'
	#contour = grid.draft; vmin = -600; vmax = -0; ctitle = 'Ice shelf draft'
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

	level = 0
	z = '0' # '257', '355'
	
	ts = 107; te = 502

	# Load time-mean velocities
	u = np.load(path+'u'+z+'_PAS851.npy')[ts:te]
	v = np.load(path+'v'+z+'_PAS851.npy')[ts:te]
	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

	THETA = False
	if THETA:
		# Load temperature
		T = np.load(path+'T'+z+'_PAS851.npy')[ts:te]
		cvmin, cvmax, ccmap, title = getPlottingVars('THETA')
		cmvmin = -2.; cvmax = 2.
		title = '(u, v); T; ' + ctitle + ' Z = ' + str(grid.RC.squeeze()[level]) + ' m'
	else:
		# Load salinity
		T = np.load(path+'S'+z+'_PAS851.npy')[ts:te]
		cvmin, cvmax, ccmap, title = getPlottingVars('SALT')
		cmvmin = -33.4; cvmax = 34.2
		title = '(u, v); S; ' + ctitle + ' Z = ' + str(grid.RC.squeeze()[level]) + ' m'
	T = tools.boundData(T, cvmin, cvmax, 0.9999)

	# Sample rate
	d = 4

	lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
	X = grid.Xsubr1D(lons); Y = grid.Ysubr1D(lats)
	contour = tools.getSubregionXY(contour, latsi, lonsi)

	TIME = np.load(path+'TIME_PAS851.npy')[ts:te]
	text_data = ptt.getTextData(TIME, 'ctime', X[1], Y[1])
	months = [mon[0:3] for mon in text_data['text']] 

	#==

	DJF = ['Dec', 'Jan', 'Feb']; MAM = ['Mar', 'Apr', 'May']
	JJA = ['Jun', 'Jul', 'Aug']; SON = ['Sep', 'Oct', 'Nov']	
	DJFts = []; MAMts = []; JJAts = []; SONts = []

	for ti in range(len(months)):
		if months[ti] in DJF:
			DJFts.append(ti)
		elif months[ti] in MAM:
			MAMts.append(ti)
		elif months[ti] in JJA:
			JJAts.append(ti)
		elif months[ti] in SON:
			SONts.append(ti)
		else:
			print(months[ti])

	# Get seaonal means of u
	uSON = np.mean(u[SONts], axis=0); uDJF = np.mean(u[DJFts], axis=0)
	uMAM = np.mean(u[MAMts], axis=0); uJJA = np.mean(u[JJAts], axis=0)
	uSON = ptt.maskBathyXY(uSON, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uDJF = ptt.maskBathyXY(uDJF, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uMAM = ptt.maskBathyXY(uMAM, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uJJA = ptt.maskBathyXY(uJJA, grid, level, subregion=True, lats=latsi, lons=lonsi)
	udata = [[uSON[::d,::d], uDJF[::d,::d]], [uMAM[::d,::d], uJJA[::d,::d]]]
	
	# Get seaonal means of v
	uSON = np.mean(v[SONts], axis=0); uDJF = np.mean(v[DJFts], axis=0)
	uMAM = np.mean(v[MAMts], axis=0); uJJA = np.mean(v[JJAts], axis=0)
	uSON = ptt.maskBathyXY(uSON, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uDJF = ptt.maskBathyXY(uDJF, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uMAM = ptt.maskBathyXY(uMAM, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uJJA = ptt.maskBathyXY(uJJA, grid, level, subregion=True, lats=latsi, lons=lonsi)
	vdata = [[uSON[::d,::d], uDJF[::d,::d]], [uMAM[::d,::d], uJJA[::d,::d]]]

	# Get seaonal means of T
	uSON = np.mean(T[SONts], axis=0); uDJF = np.mean(T[DJFts], axis=0)
	uMAM = np.mean(T[MAMts], axis=0); uJJA = np.mean(T[JJAts], axis=0)
	uSON = ptt.maskBathyXY(uSON, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uDJF = ptt.maskBathyXY(uDJF, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uMAM = ptt.maskBathyXY(uMAM, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uJJA = ptt.maskBathyXY(uJJA, grid, level, subregion=True, lats=latsi, lons=lonsi)
	uSON = uSON[::d,::d]; uDJF = uDJF[::d,::d]
	uMAM = uMAM[::d,::d]; uJJA = uJJA[::d,::d]
	uSON[-1,-1] = cvmin; uSON[-1,-1] = cvmax; uDJF[-1,-1] = cvmin; uDJF[-1,-1] = cvmax
	uMAM[-1,-1] = cvmin; uMAM[-1,-1] = cvmax; uJJA[-1,-1] = cvmin; uJJA[-1,-1] = cvmax
	Tdata = [[uSON, uDJF], [uMAM, uJJA]]

	#==
	
	Xd = X[::d]; Yd = Y[::d]
	#u = tools.boundData(u, -0.1, 0.1, 0.999); v = tools.boundData(v, -0.1, 0.1, 0.999)
	
	cmap = 'YlGn' #'plasma'

	pt.PAS_2by2_uvs(udata, vdata, Xd, Yd, C=Tdata, contour=contour, X=X, Y=Y, cmap=cmap, vmin=vmin, vmax=vmax, text_data=text_data, title=title, figsize=(7,7), dpi=300)

	quit()

#==

# Animate flux convergences of T at each level.
T_transportDepth = False
if T_transportDepth:


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


# Animate flux convergences of T at fixed level.
T_transport = False
if T_transport:

	PAS = False
	
	if PAS:
		a=1
	else:
		
		path = '/home/michai/Documents/data/MCS_038/run/'
		grid = Grid(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		level = 14

		vmin = -5.e-8; vmax = -vmin
		xlabel = 'LAT (km)'; ylabel = 'LON (km)';
		scale = 1.e3
		X = grid.XC[1,:]/scale
		Y = grid.YC[:,1]/scale
		Z = grid.RC.squeeze()
	
		dxc = scale*(X[1]-X[0]); dyc = scale*(Y[1]-Y[0]); dzc = Z[1]-Z[0]


		# Load variables and define fluxes
		T = readVariable('THETA', path, file_format='nc', meta=True)
		TIME = T['TIME'][:]
		T = T['THETA'][:, level]
		text_data = ptt.getTextData(TIME, 'month', X[1], Y[1])
		print(len(text_data['text']))
		print(T.shape)

		u = readVariable('UVEL', path, file_format='nc', meta=False)
		u = tools.interp(u[:,level], 'u'); uT = u*T

		u = readVariable('VVEL', path, file_format='nc', meta=False)
		u = tools.interp(u[:,level], 'v'); vT = u*T

		u = readVariable('WVEL', path, file_format='nc', meta=False)
		u = tools.interp(u[:,level], 'w'); wT = u*T

		# Compute convergences.		
		Tconv_lat = - tools.ddx(uT, dxc) - tools.ddy(vT, dyc)
		Tconv_vert = - tools.ddy(wT, dzc)

		Tconv_lat_ = [0 for i in range(uT.shape[0])]
		Tconv_vert_ = [0 for i in range(wT.shape[0])]
		for ti in range(uT.shape[0]):
			Tconv_tmp = tools.boundData(Tconv_lat[ti], vmin, vmax, 0.9999)
			Tconv_lat_[ti] = ptt.maskBathyXY(Tconv_tmp, grid, level, timeDep=False)
			Tconv_tmp = tools.boundData(Tconv_vert[ti], vmin, vmax, 0.9999)
			Tconv_vert_[ti] = ptt.maskBathyXY(Tconv_tmp, grid, level, timeDep=False)

	#== 

	# PLOT

	cmap = 'coolwarm'
	title = 'Theta flux convergence'

	pt.animate1by1(Tconv_lat_, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, figsize=(6,4), outname='Tconv_lat.mp4')

	pt.animate1by1(Tconv_vert_, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, figsize=(6,4), outname='Tconv_vert.mp4')

#==


PAS_WIND = True
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

	#uw = np.load(path+'tauxmean_PAS851.npy')
	#vw = np.load(path+'tauymean_PAS851.npy')


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


