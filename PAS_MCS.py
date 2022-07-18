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

HEAT_CONTENT = False
if HEAT_CONTENT:

	EXPS = ['MCS_132', 'MCS_133']
	#path = '/data/oceans_output/shelf/michai/mitgcm/'
	path = '/home/michael/Documents/data/'

	
	#pt.plot1by1(grid.bathy, mesh=True, vmin=-700, vmax=-500); quit()
	
	ts = 0
	
	THC = []
	for exp in EXPS:
	
		path_tmp = path + exp + '/run/'
		print(exp)
		print(path_tmp)
		
		grid = Grid(path_tmp)
		T = readVariable('THETA', path_tmp, meta=False)[ts:]
		
		Nt = T.shape[0]
		THC_tmp = np.zeros(Nt)
		for ti in range(Nt):
			THC_tmp[ti], TCH0 = tools.heatContentShelf(T[ti], grid)

		GJ = 10.e19
		plt.plot((THC_tmp-TCH0)/GJ, label=exp) 
	
	#plt.ylim(0.5e20, 4e20)
	plt.grid()
	plt.legend()
	plt.title(r'On-shelf heat content ($10^{19}$ J)')
	plt.xlabel('Time (months)')
	plt.savefig('heatContent.png')
	plt.show()
	
	
	quit()
			
HEAT_TRANSPORT = False
if HEAT_TRANSPORT:
	
	ts = 0	
	
	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)
	rho0 = 1030. # Units kg / m3
	# --> rho * Cp * v * T --> (kg / m3) * (m2 / s2 C) * (m / s) * C
	# --> kg / S3
	# Area integral --> kg / (m2 s3)

	#path = '/data/oceans_output/shelf/michai/mitgcm/MCS_117/run/'
	path = '/home/michael/Documents/data/MCS_133/run/'
	grid = Grid(path)

	X = grid.XC[1,:] / 1000.
	Y = grid.YC[:,1] / 1000.
	Z = grid.RC.squeeze()
	
	# Subregions for heat transport.

	lat = 95
	print(Y[lat])

	# troughW
	lons1 = [100e3, 200e3]; depth1 = [-10, -600]
	lon1 = grid.getIndexFromLon(lons1)
	depth1 = grid.getIndexFromDepth(depth1)
	label1 = 'trough W'

	# troughE
	lons2 = [415e3, 520e3]; depth2 = [-10, -600]
	lon2 = grid.getIndexFromLon(lons2)
	depth2 = grid.getIndexFromDepth(depth2)
	label2 = 'trough E'

	#lons2 = [200e3, 415e3]; depth2 = [-10, -600]
	#lon2 = grid.getIndexFromLon(lons2)
	#depth2 = grid.getIndexFromDepth(depth2)
	#label2 = 'between W and E'

	# troughE AND sill
	lons3 = [415e3, 590e3]; depth3 = [-10, -600]
	lon3 = grid.getIndexFromLon(lons3)
	depth3 = grid.getIndexFromDepth(depth3)
	label3 = 'trough E + sill'

	# For the rest, just subtract.
	label5 = 'All'
	label4 = 'Uniform lons'
	label6 = 'Uniform lons'
	
	vlines = [lons1[0]/1.e3, lons1[1]/1.e3, lons2[0]/1.e3, lons2[1]/1.e3, lons3[1]/1.e3]
	hlines = [Y[lat]]	
	#pt.plot1by1(grid.bathy, X=X, Y=Y, vmin=-1000, vmax=-300, mesh=True, hlines=hlines, vlines=vlines)
	#quit()
	#==

	T = readVariable('SALT', path, meta=False)[ts,:,:lat+1,:]#[ts:,...,lat:lat+2,:]
	#Tf, Thf = tools.theta_f(grid, S=T)
	Tf = -2. * np.ones(T.shape)
	
	#plt.plot(Tf[:,-2,-2], Z)
	#plt.plot(Thf[:,-2,-2], Z)
	#plt.show()
	#quit()


	T = readVariable('THETA', path, meta=False)[ts:,...,lat,:]

	#HC = rho0 * Cp * (T[:,:,:lat+1,:] - Tf)
	#HC = np.sum((HC * grid.DXG[:lat+1,:] * grid.hFacS[:,:lat+1,:] * grid.DRF), axis=(1,2,3))
	#plt.plot(HC); plt.show(); quit()

	v = tools.interp(readVariable('VVEL', path, meta=False)[ts:,...,lat,:], 'v')
	Tf = Tf[...,lat,:]

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	T = rho0 * Cp * v * (T - Tf)
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)

	T *= area

	# Get various T slices
	T1 = np.ma.sum(T[:, depth1[0]:depth1[1], lon1[0]:lon1[1]], axis=(1,2))
	T2 = np.ma.sum(T[:, depth2[0]:depth2[1], lon2[0]:lon2[1]], axis=(1,2))
	T3 = np.ma.sum(T[:, depth3[0]:depth3[1], lon3[0]:lon3[1]], axis=(1,2))
	T5 = np.ma.sum(T, axis=(1,2))
	T4 = T5 - T1 - T3
	T6 = T5 - T3
	

	# Normalise slices
	norm = True
	if norm:
		title = r'Meridional heat transport per unit area across shelf break (TW / m$^2$)'
		area1 = np.ma.sum(area[depth1[0]:depth1[1], lon1[0]:lon1[1]])
		T1 = T1 / area1
		area2 = np.ma.sum(area[depth2[0]:depth2[1], lon2[0]:lon2[1]])
		T2 = T2 / area2 
		area3 = np.ma.sum(area[depth3[0]:depth3[1], lon3[0]:lon3[1]])
		T3 = T3 / area3
		area5 = np.ma.sum(area)
		T5 = T5 / area5
		area4 = area5 - area1 - area3
		T4 = T4 / area4
		area6 = area5 - area3
		T6 = T6 / area6
	else:
		title = 'Meridional heat transport across shelf break (TW)'

	smooth = True
	if smooth:
		T1 = tools.smooth3(T1)
		T2 = tools.smooth3(T2)
		T3 = tools.smooth3(T3)
		T4 = tools.smooth3(T4)
		T5 = tools.smooth3(T5)
		T6 = tools.smooth3(T6)
		
	# NOW PLOT

	normval = 1.e12

	#Ts = [T1, T2, T3, T4, T5]; labels = [label1, label2, label3, label4, label5]
	Ts = [T2, T3, T5, T6]; labels = [label2, label3, label5, label6]

	vmin = -4.e-8; vmax = 2.e-8

	plt.figure(figsize=(9,4))
	for i in range(len(Ts)):
		plt.plot(Ts[i]/normval, label=labels[i])
		
	plt.ylim([vmin, vmax])
	plt.title(title)
	plt.xlabel('Time (months)')
	plt.grid(); plt.legend(); 
	plt.savefig('heatTransport.png')
	plt.show(); quit()
	
	T = np.mean(T, axis=0)
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=False)

	vmin = -2e5; vmax = -vmin
	T = tools.boundData(T, vmin, vmax, scale=0.9999)

	pt.plot1by2([T, T], X=[X,X], Y=[Z,Z], vmin=[vmin,vmin], vmax=[vmax,vmax])

	quit()	

#==

HEAT_TRANSPORT_PAS = False
if HEAT_TRANSPORT_PAS:
	
	ts = 107; te = 502
	
	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)
	rho0 = 1030. # Units kg / m3
	# --> rho * Cp * v * T --> (kg / m3) * (m2 / s2 C) * (m / s) * C
	# --> kg / S3
	# Area integral --> kg / (m2 s3)

	path = '/data/oceans_output/shelf/pahol/mitgcm/PAS_851/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy
	Y = grid.YC
	X = grid.XC
	
	ylim = [-75.5, -70.5]; xlim = [245, 262]
	X0 = xlim[0]; Lx = xlim[1]-xlim[0]
	
		
	# troughW
	lat1 = -71.5
	lons1 = [245.5, 247.5]; depth1 = [-10, -900]
		
	# troughE
	lat2 = -71.3
	lons2 = [255.5, 258.2]; depth2 = [-10, -900]

	lat1 = grid.getIndexFromLat(lat1)
	lon1 = grid.getIndexFromLon(lons1)
	depth1 = grid.getIndexFromDepth(depth1)
	label1 = 'trough W'
	
	lat2 = grid.getIndexFromLat(lat2)
	lon2 = grid.getIndexFromLon(lons2)
	depth2 = grid.getIndexFromDepth(depth2)
	label2 = 'trough E'
	
	#vlines = [lons1[0]/1.e3, lons1[1]/1.e3, lons2[0]/1.e3, lons2[1]/1.e3, lons3[1]/1.e3]
	vlines = None#[]
	hlines = [Y[lat1,0], Y[lat2,0]]
	
	# Compute xmins and xmaxs
	xmin1 = (X[0,lon1[0]] - X0) / Lx
	xmin2 = (X[0,lon2[0]] - X0) / Lx
	xmax1 = (X[0,lon1[1]] - X0) / Lx
	xmax2 = (X[0,lon2[1]] - X0) / Lx
	
	xmin = [xmin1, xmin2]; xmax = [xmax1, xmax2]
	
	pt.plot1by1(bathy, X=X, Y=Y, vmin=-1000, vmax=-300, mesh=True, hlines=hlines, xmin=xmin, xmax=xmax, vlines=vlines, xlim=xlim, ylim=ylim); quit()
	
	bathy = bathy[lat1, lon1[0]:lon1[1]]
	#plt.plot(bathy); plt.show(); quit()
	
	#==
	
	# Get v * Theta
	Tf = -2
	T = tools.interp(readVariable('VVEL', path, meta=False)[ts:te,], 'v')
	T = T * (readVariable('THETA', path, meta=False)[ts:te,] - Tf)
	T = rho0 * Cp * T
		
	# Compute merid. transport.
	T1 = tools.meridTransport(T[:,:,lat1,:], grid, yi=lat1)
	T1 = ptt.maskBathyXZ(T1, grid, yi=lat1, timeDep=True)

	T2 = tools.meridTransport(T[:,:,lat2,:], grid, yi=lat2)
	T2 = ptt.maskBathyXZ(T2, grid, yi=lat2, timeDep=True)

	#pt.plot1by2([T1[-1, depth1[0]:depth1[1], lon1[0]:lon1[1]], T1[-1]]); quit()
	T1 = np.ma.sum(T1[:, depth1[0]:depth1[1], lon1[0]:lon1[1]], axis=(1,2))
	T2 = np.ma.sum(T2[:, depth2[0]:depth2[1], lon2[0]:lon2[1]], axis=(1,2))
	
	area1 = grid.DXG[lat1] * grid.hFacS[:,lat1] * grid.DRF[:,0]
	area1 = ptt.maskBathyXZ(area1, grid, yi=lat1, timeDep=False)
	
	area2 = grid.DXG[lat2] * grid.hFacS[:,lat2] * grid.DRF[:,0]
	area2 = ptt.maskBathyXZ(area2, grid, yi=lat2, timeDep=False)
	
	# Normalise slices
	norm = True
	if norm:
		title = r'Meridional heat transport per unit area across shelf break (TW / m$^2$)'
		area1 = np.ma.sum(area1[depth1[0]:depth1[1], lon1[0]:lon1[1]])
		T1 = T1 / area1
		area2 = np.ma.sum(area2[depth2[0]:depth2[1], lon2[0]:lon2[1]])
		T2 = T2 / area2
	else:
		title = 'Meridional heat transport across shelf break (TW)'

	smooth = True
	if smooth:
		T1 = tools.smooth3(T1)
		T2 = tools.smooth3(T2)

	# NOW PLOT

	normval = 1.e12
	Ts = [T1, T2]; labels = [label1, label2]

	vmin = -4.e-8; vmax = 2.e-8

	plt.figure(figsize=(9,4))
	for i in range(len(Ts)):
		plt.plot(Ts[i]/normval, label=labels[i])
		
	plt.ylim([vmin, vmax])
	plt.title(title)
	plt.xlabel('Time (months)')
	plt.grid(); plt.legend(); 
	plt.savefig('heatTransport.png')
	
	quit()
#==

FORMSTRESS = False
if FORMSTRESS:

	PAS = 0
	
	if PAS:

		path = '/home/michael/Documents/data/PAS_8512/run/'
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


		Pb = ptt.maskBathyXY(Pb, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

		Hx = tools.getSubregionXY(Hx, latsi, lonsi)
		Hx = ptt.maskBathyXY(Hx, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

		vmin = 0.2e7; vmax = 1.e7
		Pb = tools.boundData(Pb, vmin, vmax, scale=0.9999)

		#pt.plot1by1(Pb, X=X, Y=Y, vmin=vmin, vmax=vmax, mesh=True, title='Bottom pressure')
		#pt.plot1by1(Hx, X=X, Y=Y, mesh=True, title='Hx', vmin=-0.01, vmax=0.01)
		#pt.plot1by1(H, mesh=True, title='H')

		TFS = Hx*Pb
		vmin = -1e6; vmax = -vmin
		pt.plot1by1(TFS, X=X, Y=Y, vmin=vmin, vmax=vmax, mesh=False, title='TFS')


	else:

		rho0 = 1030.
		g = 9.81

		path = '/home/michael/Documents/data/MCS_123/run/'
		grid = Grid(path)

		depth = - grid.bathy
		SSH = readVariable('ETAN', path, meta=False)
		depth = depth + SSH

		botx = readVariable('botTauX', path, meta=False)

		Pb = readVariable('PHIBOT', path, meta=False)
		Pb = depth * rho0 * g + Pb * rho0

		depth = ptt.maskBathyXY(depth, grid, 0, timeDep=True)
		Hx = tools.ddx(depth, grid.DXG)

		Hx = np.where(np.abs(Hx)>0.05, 0, Hx)

		Pb = ptt.maskBathyXY(Pb, grid, 0, timeDep=True)
		botx = ptt.maskBathyXY(botx, grid, 0, timeDep=True)

		TFS = Pb * Hx		

		#vmin = 0.2e7; vmax = 1.e7
		#Pb = tools.boundData(Pb, vmin, vmax, scale=0.9999)

		vmin = -1.e5; vmax = -vmin
		TFS = tools.boundData(TFS, vmin, vmax, scale=0.9999)
		print(np.mean(TFS[-1]))
		print(np.mean(botx[-1]))

		pt.plot1by2([TFS[-1]/rho, botx[-1]], mesh=False, titles=['TFS', 'botTauX'])

	quit()
		
	
#==

Bathy = False
if Bathy:

	path = '/home/michael/Documents/data/PAS_851/run/'
	#path = 	'/home/michael/Documents/data/PISOMIP_001/run/'
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
ASF = 0
if ASF:

	# Window 1
	lat1 = [-75, -69]; lon1 = 250; depth1 = [0, -1000]

	# Window 2
	lat2 = [-75, -69]; lon2 = 260; depth2 = depth1
	
	# Window3
	lat3 = [-75, -69]; lon3 = 258; depth3 = depth1

	# Window4
	lat4 = [-75, -69]; lon4 = 254; depth4 = depth1

	lats = [lat1, lat2, lat3, lat4]
	lons = [lon1, lon2, lon3, lon4]
	depths = [depth1, depth2, depth3, depth4]
	
	#==

	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	T = np.load(path+'Tmean_PAS851.npy')

	X = grid.XC#[1,:]
	Y = grid.YC#[:,1]
	Z = grid.RC.squeeze()

	##plt.plot(T[:20,180,20], Z[:20]); plt.xlim(33.5, 35.); plt.grid(); plt.show();


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

thetaHeight = True
if thetaHeight:

	PAS = False
	ANIM = True
	THERM = -0.5

	if PAS:

		if ANIM:
			path = '/home/michael/Documents/data/PAS_8512/run/'
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
			path = '/home/michael/Documents/data/PAS_851/run/'
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


		exp = 'MCS_125'
		path = '/home/michael/Documents/data/'+exp+'/run/'
		grid = Grid(path)
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()
		bathy = grid.bathy

		if ANIM:

			# Use this if loading from pre-computed npy file.
			if 1:
				d = 2
				ThermZ = np.load(path+'ThermZ_m05_MCS_125.npy')[::d]
				Nt = ThermZ.shape[0] * d
				TIME = np.linspace(1,Nt,Nt)[::d]*86400*30
			else:
				ts = 0
				T = readVariable('THETA', path, meta=True)
				TIME = T['TIME'][ts:]
				T = T['THETA'][ts:]

				# Get z-indices of level with Theta closest to THERM.
				ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True)
				
			#==
			
			#plt.plot(ThermZ[-1,:,120]);plt.show(); quit()
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)

			text_data = ptt.getTextData(TIME, 'month', X[1], Y[1])

		else:
			T = np.mean(readVariable('THETA', path, meta=False)[108:120], axis=0)
			# Get z-indices of level with Theta closest to THERM.
			ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True, timeDep=False)
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)
			ThermZ = np.where(ThermZ<grid.bathy+10, np.nan, ThermZ)

		title = exp + ' ' + str(THERM) + ' deg. isotherm height'
	
		#vmin = -450; vmax = -250	
		vmin = -500; vmax = -100		
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


thetaHeight_timeSeries = False
if thetaHeight_timeSeries:

	THERM = -1.5

	exps = [['MCS_116', (200,100)], ['MCS_117', (200,100)]]
	data = []
	labels = []
	
	for exp in exps:
		path = '/home/michael/Documents/data/'+exp[0]+'/run/'
		grid = Grid(path)
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()
		bathy = grid.bathy

		T = readVariable('THETA', path, meta=False)
		# Get z-indices of level with Theta closest to THERM.
		ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True, timeDep=True)
		ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)
		#ThermZ = np.where(ThermZ<grid.bathy+10, np.nan, ThermZ)
		
		data.append(ThermZ[:, exp[1][1], exp[1][0]])
		labels.append(exp[0])
		
	#==

	for di in range(len(data)):
		plt.plot(data[di], label=labels[di])
	plt.legend()
	plt.show()

	quit()
	
#==


quiver = 0
if quiver:
	
	PAS = False

	level = 0
	#depth = -350; level = grid.getIndexFromDepth(depth)
	print(level)


	if PAS:
		path = '/home/michael/Documents/data/PAS_851/run/'
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

		#path = '/home/michael/Documents/data/PAS_666/run/'
		path = '/home/michael/Documents/data/PISOMIP_004/run/'
		ts = -1
		
		grid = Grid(path)
		bathy = grid.bathy

		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		xlabel = 'LON (km)'; ylabel = 'LAT (km)'

		VAR = 'ETAN'
		vmin, vmax, cmap, title = getPlottingVars(VAR)

		contour = readVariable(VAR, path, meta=False)[ts]
		u = readVariable('UVEL', path, meta=False)[ts]
		v = readVariable('VVEL', path, meta=False)[ts]
		
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
	
btpcStr = True
if btpcStr:

	PAS = 0

	if PAS:
		path = '/home/michael/Documents/data/PAS_851/run/'
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

		exp = 'MCS_132'
		path = '/home/michael/Documents/data/'+exp+'/run/'
		grid = Grid(path)

		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		u = readVariable('UVEL', path, meta=False)
		btpStr = tools.barotropicStreamfunction(u, grid, timeDep=True)
		btpStr = ptt.maskBathyXY(btpStr, grid, zi=0, timeDep=True)
		bathy = grid.bathy
		levels = [-900, -800, -700, -600, -501]

		
		SUBR = False
		if SUBR:
			lats = [0, 300e3]; lons = [0,600e3]
			latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
			btpStr = tools.getSubregionXY(T, latsi, lonsi); bathy = tools.getSubregionXY(bathy, latsi, lonsi);
			X = grid.Xsubr1D(lons)/1000; Y = grid.Ysubr1D(lats)/1000
		
		title = exp + ' barotropic streamfunction'

	#==
	
	vmin = 0; vmax = 4
	#t.plot1by1(btpStr[-1])
	
	btpsStr = tools.boundData(btpStr, vmin, vmax)
	cmap = 'turbo'#'gist_rainbow'
	
	text = ['month ' + str(ti) for ti in range(btpStr.shape[0])]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

	pt.animate1by1(btpStr, X, Y, title=title, cmap=cmap, mesh=True, text_data=text_data, outname='btpStr.mp4', vmin=vmin, vmax=vmax)
	
	#pt.plot1by1(T, X=X, Y=Y, mesh=True, cmap=cmap, contour=bathy, contourlevels=levels, title=title)#, vmin=vmin, vmax=vmax)
	#pt.plot1by2([umean[-1], umean2[-1]], mesh=True, vmin=vmin, vmax=vmax)
	quit()

#==

# Time-mean difference between surface and bottom flow.
brclnc = False
if brclnc:

	PAS = 0

	if PAS:
		path = '/home/michael/Documents/data/PAS_851/run/'
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
		#us = u[0] - ub; vs = v[0] - vb
		# Or leave as bottom vels
		us = ub; vs = vb

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

		path = '/home/michael/Documents/data/MCS_123/run/'
		grid = Grid(path)
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, meta=False)[60:]
		v = readVariable('VVEL', path, meta=False)[60:]
		#u = np.mean(u, axis=0); v = np.mean(v, axis=0)
	
		# Transfer to t-points.
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Get bottom val.
		ub = tools.bottom(u, grid, 'u', timeDep=True)
		vb = tools.bottom(v, grid, 'v', timeDep=True)

		# Get 'baroclinicity'
		#us = u[0] - ub; vs = v[0] - vb
		# Or leave as bottom vels
		us = ub; vs = vb

		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)
		vmin = -1000; vmax = -500

		title = ''#MCS038 (u_t-u_b, v_t-v_b) & bathy'
		xlabel = 'LON (km)'; ylabel = 'LAT (km)'

		# Sample rate
		d = 3

	#==

	# PLOT
 
	text = ['month ' + str(ti) for ti in range(us.shape[0])]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}
		
	us = us[..., ::d, ::d]; vs = vs[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]

	pt.quiver1by1(us[-1], vs[-1], Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, title=title, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, outname='animateBrclnc.mp4')

	pt.animate1by1quiver(us, vs, Xd, Yd, contour=contour, X=X, Y=Y, cmap='YlGn', contourf=False, vmin=vmin, vmax=vmax, text_data=text_data, title=title, figsize=(7,4), dpi=300)
	
#==

# Animate velocity vectors and temperature. Each frame is different level.
animateUVTdepth = False
if animateUVTdepth:

	PAS = 1
	SUBR = 1
	BOT = 0

	if PAS:
		path = '/home/michael/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		contour = grid.bathy; vmin = -600; vmax = -200; ctitle = 'bathy'
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
		d = 8

	#==

	else:

		#ts = 20
		ts = 60; te = 120

		path = '/home/michael/Documents/data/MCS_117/run/'
		grid = Grid(path)
		contour = grid.bathy; ctitle = 'bathy'
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		vmin = -1000; vmax = -400
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, meta=False)[ts:]
		v = readVariable('VVEL', path, meta=False)[ts:]
		u = np.mean(u, axis=0); v = np.mean(v, axis=0)
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Load temperature
		T = readVariable('THETA', path, meta=False)[ts:]
		T = np.mean(T, axis=0)
		cvmin, cvmax, ccmap, title = getPlottingVars('THETA')

		# Sample rate
		d = 6

	for zi in range(len(Z)):
		u[zi] = ptt.maskBathyXY(u[zi], grid, zi, timeDep=False)
		v[zi] = ptt.maskBathyXY(v[zi], grid, zi, timeDep=False)
		T[zi] = ptt.maskBathyXY(T[zi], grid, zi, timeDep=False)


	if PAS and SUBR:
		lats = [-76, -70.5]; lons = [245, 280]#[230, 270]#
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
	outpath = '/home/michael/Documents/Python/mitgcmPy/images/testingMov/'

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
animateUVT = True
if animateUVT:

	PAS = False
	
	if PAS:

		path = '/home/michael/Documents/data/PAS_8512/run/'
		grid = Grid_PAS(path)
		contour = grid.bathy; vmin = -800; vmax = -100; ctitle = 'bathy'
		#contour = grid.draft; vmin = -600; vmax = -0; ctitle = 'Ice shelf draft'
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)


		vmin = -800; vmax = -100

		level = 0
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

		#path = '/data/oceans_output/shelf/michai/mitgcm/MCS_126/run/'
		path = '/home/michael/Documents/data/MCS_133/run/'
		grid = Grid(path)
		contour = grid.bathy
		contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

		#pt.plotMbyN(grid.bathy, mesh=True); quit()

		depth = -490; level = grid.getIndexFromDepth(depth)

		vmin = -600; vmax = -500
		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Load time-mean velocities
		u = readVariable('UVEL', path, meta=False)[ts:, level]
		print(u.shape)
		v = readVariable('VVEL', path, meta=False)[ts:, level]
		u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

		# Load temperature
		T = readVariable('THETA', path, meta=False)[ts:, level]
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

# Rough vorticity budget
vortBudget = True
if vortBudget:

	path = '/home/michael/Documents/data/MCS_117/run/'
	grid = Grid(path)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	dx = (X[1]-X[0])*1000.
	dy = (Y[1]-Y[0])*1000.
	xlabel = 'LON (km)'; ylabel = 'LAT (km)'

	g = 9.81; rho0 = 1028.5
	ts = -24 

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	u = readVariable(VAR, path, file_format='nc', meta=True)
	TIME = u['TIME'][ts:]
	u = u[VAR][ts:,]

	VAR = 'VVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	v = readVariable(VAR, path, file_format='nc', meta=False)[ts:]

	# Advection term
	
	uu = np.mean(u*u, axis=0)
	uv = np.mean(u*v, axis=0)
	vv = np.mean(v*v, axis=0)
	
	uu = tools.depthIntegral(uu, grid, timeDep=False, norm=False)
	uv = tools.depthIntegral(uv, grid, timeDep=False, norm=False)
	vv = tools.depthIntegral(vv, grid, timeDep=False, norm=False)
			
	advu = - tools.ddx(uu, dx) - tools.ddy(uv, dy)
	advv = - tools.ddx(uv, dx) - tools.ddy(vv, dy)

	adv = rho0 * (tools.ddx(advv, dx) - tools.ddy(advu, dy))

	# Beta term
	beta = 6.42e-12
	bV = rho0 * beta * tools.depthIntegral(np.mean(v, axis=(0)), grid, timeDep=False, norm=False)
	
	# Vortex stretching term.
	depth = - grid.bathy
	PHIBOT = readVariable('PHIBOT', path, file_format='nc', meta=False)[ts:]
	Pb = depth * rho0 * g + PHIBOT * rho0
	Pb = np.mean(Pb, axis=0)
	
	J = -tools.ddy(Pb, dy) * tools.ddx(depth, dx) + tools.ddx(Pb, dx) * tools.ddy(depth, dy)
	
	vmin = -1.e-4; vmax = -vmin
	cbar = [False, False, True]
	titles = ['Vort. advection', 'beta term', 'Vort. stretching']

	#==

	adv = ptt.maskBathyXY(adv, grid, 0)
	bV = ptt.maskBathyXY(bV, grid, 0)
	J = ptt.maskBathyXY(J, grid, 0)
		
	pt.plot1by3([adv, bV, J], vmin=vmin, vmax=vmax, mesh=True, cbar=cbar, titles=titles, width_ratios=[1,1,1.1])

	quit()
	
#==
	
# Animate velocity vectors and temperature at fixed level.
animateUVT_npy = False
if animateUVT_npy:


	ts = 0; te = -1

	path = '/home/michael/Documents/data/MCS_127/run/'
	fname = '_MCS_127_lvls_0_24.npy'
	levels = [0, 24]
	leveli = 1
	level = levels[leveli]
	
	grid = Grid(path)
	contour = grid.bathy
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

	vmin = -1000; vmax = -300
	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	# Load time-mean velocities
	u = np.load(path+'UVEL'+fname)[:, leveli]
	v = np.load(path+'VVEL'+fname)[:, leveli]
	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

	# Load temperature
	T = np.load(path+'THETA'+fname)[:, leveli]
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

		path = '/home/michael/Documents/data/PAS_8512/run/'
	
	else:

		ts = 0; te = 112

		path = '/home/michael/Documents/data/MCS_113/run/'
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
		u = readVariable('UVEL', path, meta=False)[ts:]
		v = readVariable('VVEL', path, meta=False)[ts:]

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


	path = '/home/michael/Documents/data/PAS_8512/run/'
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
		path = '/home/michael/Documents/data/PAS_851/run/'
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
		
		path = '/home/michael/Documents/data/MCS_038/run/'
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
		T = readVariable('THETA', path, meta=False)[ts:te]

		u = readVariable('UVEL', path, meta=False)[ts:te]
		u = tools.interp(u, 'u'); uT = np.mean(u*T, axis=0);

		u = readVariable('VVEL', path, meta=False)[ts:te]
		u = tools.interp(u, 'v'); vT = np.mean(u*T, axis=0);

		u = readVariable('WVEL', path, meta=False)[ts:te]
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
		
		path = '/home/michael/Documents/data/MCS_038/run/'
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
		T = readVariable('THETA', path, meta=True)
		TIME = T['TIME'][:]
		T = T['THETA'][:, level]
		text_data = ptt.getTextData(TIME, 'month', X[1], Y[1])
		print(len(text_data['text']))
		print(T.shape)

		u = readVariable('UVEL', path, meta=False)
		u = tools.interp(u[:,level], 'u'); uT = u*T

		u = readVariable('VVEL', path, meta=False)
		u = tools.interp(u[:,level], 'v'); vT = u*T

		u = readVariable('WVEL', path, meta=False)
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


PAS_WIND = False
if PAS_WIND:

	path = '/home/michael/Documents/data/PAS_8512/run/'
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
	outpath = '/home/michael/Documents/Python/mitgcmPy/images/testingMov/'

	d = 6
	uw = uw[..., ::d, ::d]; vw = vw[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]

	pt.quiver1by1(uw, vw, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, vmin=vmin, vmax=vmax, title=title)


