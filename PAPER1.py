# PAPER 1

# Scripts for producing figures for BAS PAPER 1.

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

import iceShelfFront as isf

import time

#==========================================================

print(' I could change all Lats to Lat.s, or to Latitudes')

# Plot ERA5 winds, ocean surface stress and bathymetry in PAS.
FIGURE1a = 0
if FIGURE1a:

	#rho0 = 1028.5; Cd = 0.006
	rho0 = 1.3; Cd = 0.003
	#rho0 = 1.3; Cd = 0.006
	
	# These from MITgcm default params.
	cDrag_1 = 2.70e-3
	cDrag_2 = 0.142e-3
	cDrag_3 = 0.0764e-3
	Cd = cDrag_1 + cDrag_2 + cDrag_3
	
	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	
	bathy = grid.bathy
	bathy = ptt.maskBathyXY(bathy, grid, zi=0)	
	bathy = tools.boundData(bathy, -5000, 0)

	lats = [-75.5, -70.5]; lons = [245, 262] 
	outline = [lons, lats]

	icel, icea = isf.get_ice_shelf_front(grid, grid.iceC, grid.landC)

	X = grid.XC; Y = grid.YC
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]

	taux = np.load(path+'tauxmean_PAS851.npy'); tauy = np.load(path+'tauymean_PAS851.npy')
	windx = np.load(path+'uwindmean_PAS851.npy'); windy = np.load(path+'vwindmean_PAS851.npy')

	speed = (windx**2 + windy**2)**0.5
	windx = rho0 * Cd * speed * windx
	windy = rho0 * Cd * speed * windy

	windx = ptt.maskBathyXY(windx, grid, zi=0); windy = ptt.maskBathyXY(windy, grid, zi=0)
	windx = ptt.maskDraftXY(windx, grid, zi=0); windy = ptt.maskDraftXY(windy, grid, zi=0)	
	taux = ptt.maskBathyXY(taux, grid, zi=0); tauy = ptt.maskBathyXY(tauy, grid, zi=0)	
	taux = ptt.maskDraftXY(taux, grid, zi=0); tauy = ptt.maskDraftXY(tauy, grid, zi=0)	
		
	#pt.plot1by2([taux**2+tauy**2, windx**2+windy**2], vmin=[0,0], vmax=[0.01, 0.01]);

	vecx = [windx, taux]; vecy = [windy, tauy]

	# Sample rate
	d = 20

	paras = [-74, -72, -70, -68, -66, -64]
	merids = [230, 240, 250, 260, 270]

	xlabel = 'LON'; ylabel = 'LAT'
	title = r'Surface stresses (N m$^{-2}$), bathy. (m)'

	pt.quiver1by1Basemap(vecx, vecy, X, Y, d, lat_0, lon_0, contourf=bathy, mesh=True, cmap='plasma', parallels=paras, meridians=merids, isf=icea, outline=outline, figsize=(5,3), title=title)

	quit()

#==

# Zoomed in bathymetry
FIGURE1b = 0
if FIGURE1b:

	path = '/home/michael/Documents/data/PAS_851/run/'
	#path = 	'/home/michael/Documents/data/PISOMIP_001/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy

	lats = [-75.5, -70.5]; lons = [245, 262]
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	X, Y = grid.XYsubr(lons, lats)
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)
	bathy = ptt.maskBathyXY(bathy, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	bathy = tools.boundData(bathy, vmin=-1000, vmax=-200)

	# Extra plot params

	paras = [-75, -74, -73, -72, -71]
	merids = [245, 250, 255, 260]

	# Labels for troughs
	p1 = {'x':(1.91e5,1.536e6), 't':'1'}
	p2 = {'x':(1.025e6, 6.e5), 't':'2'}
	p3 = {'x':(1.386e6, 1.574e6), 't':'3'}
	labelData = [p1, p2, p3]

	# PLOT 

	pt.plot1by1Basemap(bathy, X, Y, lat_0, lon_0, cmap='plasma', mesh=False, vmin=-1000, vmax=-200, contourfNlevels=17, figsize=(4,3), parallels=paras, meridians=merids, title='bathymetry', labelData=labelData)

	quit()

#==

FIGURE2 = 0
if FIGURE2:

	levels = [0,16]	

	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy
	vmin = -800; vmax = -200
	bathy = tools.boundData(bathy, vmin, vmax)

	X = grid.XC; Y = grid.YC
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]

	lats = [-75.5, -70.5]; lons = [245, 262]
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	X, Y = grid.XYsubr(lons, lats)
	Z = grid.RC.squeeze()

	# Load time-mean velocities
	u = np.load(path+'umean_PAS851.npy')
	v = np.load(path+'vmean_PAS851.npy')

	# Load temperature
	T = np.load(path+'Tmean_PAS851.npy')

	#u[levels[1]] = 0.5 * (u[levels[1]]+u[levels[1]+1])
	#v[levels[1]] = 0.5 * (v[levels[1]]+v[levels[1]+1])
	#T[levels[1]] = 0.5 * (T[levels[1]]+T[levels[1]+1])

	# Get two levels
	u = u[levels]
	v = v[levels]
	T = T[levels]

	cvmin = -1.6; cvmax = -cvmin
	#cvmin = -2; cvmax = 2
	T = tools.boundData(T, cvmin, cvmax, 0.9999)

	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')
	u = tools.boundData(u, -0.4, 0.4, 0.9999); v = tools.boundData(v, -0.4, 0.4, 0.9999)

	# Get subregion and mask
	bathy = tools.getSubregionXY(bathy, latsi, lonsi)
	bathy = [bathy, bathy]
	u = tools.getSubregionXY(u, latsi, lonsi)
	v = tools.getSubregionXY(v, latsi, lonsi)
	T = tools.getSubregionXY(T, latsi, lonsi)
	for li in range(len(levels)):
		u[li] = ptt.maskBathyXY(u[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		v[li] = ptt.maskBathyXY(v[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		T[li] = ptt.maskBathyXY(T[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		u[li] = ptt.maskDraftXY(u[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		v[li] = ptt.maskDraftXY(v[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		T[li] = ptt.maskDraftXY(T[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		bathy[li] = ptt.maskBathyXY(bathy[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)
		bathy[li] = ptt.maskDraftXY(bathy[li], grid, levels[li], timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	# Sample rate
	d = 4
	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]; T = T[..., ::d, ::d]
	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd, Xd]; Yd = [Yd, Yd]
	X = [X, X]; Y = [Y, Y]

	T[0,-10,0] = cvmin; T[0,-11,1] = cvmax;
	T[1,-10,0] = cvmin; T[1,-11,1] = cvmax;

	paras = [-75, -74, -73, -72, -71]
	merids = [245, 250, 255, 260]
	
	vmin = [vmin, vmin]; vmax = [vmax, vmax]

	title1 = '(a) Surface flow, Pot. Temp. (deg. C), Bathy. (m)'
	title2 = '(b) Deep flow, Pot. Temp. (deg. C), Bathy. (m)'
	title = [title1, title2]

	fontdict = {'fontsize':12, 'color':'k'}
	text0 = {'text':'z = '+str(Z[levels[0]])+' m', 'xloc':5.e4, 'yloc':1.78e6, 'fontdict':fontdict}
	text1 = {'text':'z = '+str(Z[levels[1]])+' m', 'xloc':5.e4, 'yloc':1.78e6, 'fontdict':fontdict} 
	text_data = [text0, text1]

	# Labels for troughs
	p1 = {'x':(1.91e5,1.536e6), 't':'1'}
	p2 = {'x':(1.025e6, 6.e5), 't':'2'}
	p3 = {'x':(1.386e6, 1.574e6), 't':'3'}
	labelData = [[p1, p2, p3],[]]

	# PLOT

	pt.quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, C=T, ccmap='seismic', contourf=bathy, X=X, Y=Y, mesh=False, contourfNlevels=13, vmin=vmin, vmax=vmax, cmap='YlGn', parallels=paras, meridians=merids, width_ratios=[0.95,1.2], title=title, fontsize=8, figsize=(7.7, 3), text_data=text_data, labelData=labelData)

	quit()
	
#==

# Plot of isotherm height and ASF in PAS.
FIGURE3 = 0
if FIGURE3:

	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy
	contourlevels = [-1000, -600]

	T = np.load(path+'Tmean_PAS851.npy')

	X = grid.XC; Y = grid.YC; Z = grid.RC.squeeze()
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]


	lats = [-75.5, -70.5]; lons = [245, 262]
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	# FIRST PANEL: ISOTHERM HEIGHT

	vmin1 = -500; vmax1 = -200	

	THERM = -0.5
	ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True, timeDep=False)
	ThermZ = tools.getSubregionXY(ThermZ, latsi, lonsi)
	
	#ThermZ = np.where(ThermZ!=ThermZ, np.nan, ThermZ)
	ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	a=np.where(ThermZ!=ThermZ, 1, 0)
	#ThermZ = tools.boundData(ThermZ, vmin1, vmax1)

	# SECOND PANEL: THETA SLICE
	
	vmin2 = -2.0; vmax2 = 2.0

	slice_lon = 252.5
	xi = grid.getIndexFromLon(slice_lon)

	zs = [0, -1000]
	zi = grid.getIndexFromDepth(zs)
	Zsubr = grid.Zsubr1D(zs)

	lats2 = [-78., -69.]
	latsi2 = grid.getIndexFromLat(lats2)
	Ysubr = grid.Ysubr1D(lats2)

	T = tools.getSubregionYZ(T[...,xi], latsi2, zi)
	T = ptt.maskBathyYZ(T, grid, xi=xi, subregion=True, lats=latsi2, depths=zi)
	T = ptt.maskDraftYZ(T, grid, xi=xi, subregion=True, lats=latsi2, depths=zi)
	T = tools.boundData(T, vmin2, vmax2, scale=0.9999)

	#==

	vmin = [vmin1, vmin2]; vmax = [vmax1, vmax2]	

	paras = [-75, -74, -73, -72, -71]
	merids = [245, 250, 255, 260]

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)
	X, Y = grid.XYsubr(lons, lats)
	Xp = [X, Ysubr]; Yp = [Y, Zsubr]

	xlabels = [None, 'Lat (deg.)']; ylabels = [None, 'Depth (m)']
	titles = ['(a) -0.5 deg. isotherm depth (m)', '(b) Pot. Temp. (deg. C)']

	pt.plot1by2Basemap1([ThermZ, T], Xp, Yp, lat_0, lon_0, mesh=False, contour=bathy, contourlevels=contourlevels, contourfNlevels=13, vmin=vmin, vmax=vmax, parallels=paras, meridians=merids, yline=[slice_lon, None], xlabels=xlabels, ylabels=ylabels, fontsize=11, titles=titles, cmaps=['YlOrRd', 'seismic'])

	quit()
	
#==

# First figure of idealised model. Panel 1: bathymetry, wind forcing. Panel 2: T/S relaxation profiles.
FIGURE4 = 0
if FIGURE4:

	path = '/home/michael/Documents/data/MCS_104/run/'
	grid = Grid(path)
	bathy = grid.bathy

	Trel, Srel = tools.getTSrel(grid, salttoptop=33.2, salttop=33.5, saltbot=34.5, temptop=-1.8, tempbot=1.0, tclinewid=160., tclinetop=-200, hclinewid=160, hclinetop=-200, shift=0)
	
	#S0=[3.320e+01,3.328e+01,3.335e+01,3.342e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.356e+01,3.369e+01,3.381e+01,3.394e+01,3.406e+01,3.419e+01,3.431e+01,3.444e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01]
#	plt.plot(Srel-S0); plt.show()
#	quit()

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	X2 = 300

	vmin = -1000; vmax = -500

	# Wind
	taux = 300+200*tools.get05cosWind(grid.Nx,grid.Ny)[:,0]
	
	# Latitude grid points to plot arrows at.
	ys = [20, 40, 60, 80, 120, 140, 160, 180] 
	d = 20
	qx = X2
	qy = Y[ys]
	
	qdx = 0.95*(taux[ys] - 300)
	qdy = 0

	#==

	import matplotlib.gridspec as gridspec
	width_ratios = (2,1.)
	fs = 14

	fig = 	plt.figure(figsize=(10,4), dpi=200)
	gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)

	ax = fig.add_subplot(gs[0,0])

	plt.pcolormesh(X, Y, bathy, vmin=vmin, vmax=vmax, cmap='YlGn')
	plt.colorbar()

	for qi in range(len(qy)):
		plt.arrow(qx, qy[qi], qdx[qi], qdy, length_includes_head=True, head_width=10, color='k')
	plt.axvline(X2, color='k', linestyle='--')

	ax.set_xlabel('Lon (km)', fontsize=fs)
	ax.set_ylabel('Lat (km)', fontsize=fs)

	ax.set_aspect(1)
	plt.plot(taux, Y, color='k')

	plt.title('(a) Bathymetry (m), Wind stress')

	#==

	ax = fig.add_subplot(gs[0,1])
	
	ln1 = ax.plot(Trel, Z, color='r', label='Pot. Temp. (deg. C)')
	ax.tick_params(axis='x', labelcolor='r')
	ax.set_xlim(-2, 1.15)
	ax.set_ylabel('Depth (m)', fontsize=fs)
	ax2 = ax.twiny()
	
	#ax2.plot(Srel, Z, color='r', label='Pot. Temp.')
	ln2 = ax2.plot(Srel, Z, color='k', label='Salinity (g/kg)')
	ax2.tick_params(axis='x', labelcolor='k')

	plt.title('(b) Relaxation profiles')

	lns = ln1+ln2
	labs = [l.get_label() for l in lns]
	ax.legend(lns, labs, loc=3, prop={'size': 9})

	plt.tight_layout()
	#ax2.plot(1, 1, ' ', color='r', label="Extra label on the legend")

	plt.show()
	
	quit()	

#==

# Zonal-mean zonal velocity, temperature and SSH.
FIGURE5 = 0
if FIGURE5:

	path = '/home/michael/Documents/data/MCS_104/run/'
	grid = Grid(path)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	# These for snapshot at end of simulation
	#u = np.mean(readVariable('UVEL', path, meta=False)[-1], axis=-1)
	#T = np.mean(readVariable('THETA', path, meta=False)[-1], axis=-1)
	#h = np.mean(readVariable('ETAN', path, meta=False)[-1], axis=-1)

	ts = -12
	# These for mean over last ts months
	u = np.mean(readVariable('UVEL', path, meta=False)[ts:,], axis=(0,-1))
	T = np.mean(readVariable('THETA', path, meta=False)[ts:], axis=(0,-1))
	h = np.mean(readVariable('ETAN', path, meta=False)[ts:], axis=(0,-1))

	
	u = ptt.maskBathyYZ(u, grid, timeDep=False, xi=120)
	T = ptt.maskBathyYZ(T, grid, timeDep=False, xi=120)

	vmin = -0.2; vmax = 0.2
	u[-1,-1] = vmin
	u = tools.boundData(u, vmin, vmax, scale=0.9999)

	#==

	import matplotlib.gridspec as gridspec
	height_ratios = (1,2.)
	fs = 12
	fs2 = 10
	Tlevels = [-0.5, 0, 0.5]

	fig = plt.figure(figsize=(4,8), dpi=200, constrained_layout=True)
	gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=height_ratios)

	ax = fig.add_subplot(gs[0,0])

	plt.plot(Y[1:-1], h[1:-1], color='k')
	
	ax.axes.xaxis.set_ticklabels([])
	plt.ylabel('SSH (m)', fontsize=fs)

	plt.title('(a) Zonal-mean SSH anomaly', fontsize=fs)
	plt.xlim(Y[0], Y[-1])
	plt.ylim(-0.22, 0.22)
	plt.grid()

	#==

	ax = fig.add_subplot(gs[1,0])
	ax.patch.set_color('.5')

	plt.contourf(Y, Z, u, vmin=vmin, vmax=vmax, cmap='coolwarm', levels=15)
	plt.colorbar()

	plt.contour(Y, Z, T, levels=Tlevels, colors='k', linestyles='solid', linewidths=0.8)

	plt.title('(b) Zonal-mean zonal velocity (m/s) & pot. temp.', fontsize=fs2)
	plt.xlabel('Lat (km)', fontsize=fs)
	plt.ylabel('Depth (m)', fontsize=fs)
	plt.grid()

	plt.show()
	
	quit()

#==

# New walled simulation. Show SSH with bathy, Uvel at two longitudes.
FIGURE6 = 0
if FIGURE6:

	path = '/home/michael/Documents/data/MCS_108/run/'
	grid = Grid(path)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	xis = [0, 120]

	# These for snapshot at end of simulation
	ts = -1#116
	u = readVariable('UVEL', path, meta=False)[ts][...,xis]
	T = readVariable('THETA', path, meta=False)[ts][...,xis]
	h = readVariable('ETAN', path, meta=False)[ts]

	# These for mean over last ts months
	#ts = -12
	#u = np.mean(readVariable('UVEL', path, meta=False)[ts:,...,xis], axis=0)
	#T = np.mean(readVariable('THETA', path, meta=False)[ts:,...,xis], axis=0)
	#h = np.mean(readVariable('ETAN', path, meta=False)[ts:], axis=0)

	h = ptt.maskBathyXY(h, grid, timeDep=False, zi=0)
	u0 = ptt.maskBathyYZ(u[...,0], grid, timeDep=False, xi=xis[0])
	T0 = ptt.maskBathyYZ(T[...,0], grid, timeDep=False, xi=xis[0])
	u1 = ptt.maskBathyYZ(u[...,1], grid, timeDep=False, xi=xis[1])
	T1 = ptt.maskBathyYZ(T[...,1], grid, timeDep=False, xi=xis[1])

	vmin = -0.2; vmax = 0.2
	u0 = tools.boundData(u0, vmin, vmax, scale=0.9999)
	u1 = tools.boundData(u1, vmin, vmax, scale=0.9999)

	vminh = -0.2; vmaxh = 0.2
	h = tools.boundData(h, vminh, vmaxh, scale=0.9999)

	#==
	
	X = [X, Y, Y]
	Y = [Y, Z, Z]

	figsize = (12,3)
	width_ratios = [1,1.02,1]	
	
	data = [h, u0, u1]	
	cmaps = ['jet', 'coolwarm', 'coolwarm']
	contour = [grid.bathy, T0, T1]
	contourlevels = [[-950, -750, -550], [-0.5, 0., 0.5], [-0.5, 0., 0.5]]
	vmin = [vminh, vmin, vmin]; vmax = [vmaxh, vmax, vmax]

	cbar = [True, False, True]
	xlabels = ['Lon (km)', 'Lat (km)', 'Lat (km)']
	ylabels = ['Lat (km)', 'Depth (m)', None]
	titles = ['(a) SSH anomaly (m)', '(b) Zonal velocity (m/s), lon = 0 km', '(c) Zonal velocity (m/s), lon = 300 km']	
	fs = 10

	lonticks = [100,200,300,400,500,500]
	latticks = [0, 100, 200, 300, 400, 500]
	zticks = [-1000, -800, -600, -400, -200]

	xticks = [lonticks, latticks, latticks]
	yticks = [latticks, zticks, zticks]
	yticksvis = [True, True, False]

	pt.plot1by3(data, X=X, Y=Y, cmaps=cmaps, contour=contour, contourlevels=contourlevels, vmin=vmin, vmax=vmax, contourfNlevels=17, cbar=cbar, xlabels=xlabels, ylabels=ylabels, titles=titles, fontsize=fs, figsize=figsize, width_ratios=width_ratios, xticks=xticks, yticks=yticks, yticksvis=yticksvis, save=True, show=False)

	quit()
	
#==

# Show bottom flow for bathy, bathyE, bathyS and bathyES cases.
FIGURE7 = 0
if FIGURE7:

	root_ = '/home/michael/Documents/data/'
	runs = [['MCS_108', 'MCS_116'], ['MCS_117', 'MCS_118']]
	#runs = [runs[0]]

	bathy = [[],[]]	
	uvec = [[],[]]
	vvec = [[],[]]
	
	d = 8
	#level = 25
	level = 22

	M = 2; N = 2
	for col in range(M):
		for row in range(N):
			run = runs[row][col]
			path = root_ + run + '/run/'
			grid = Grid(path)

			ts = 108; te = 120
			u = np.mean(readVariable('UVEL', path, meta=False)[ts:te, level], axis=0)
			v = np.mean(readVariable('VVEL', path, meta=False)[ts:te, level], axis=0)
			
			u = ptt.maskBathyXY(u, grid, zi=level, timeDep=False)
			v = ptt.maskBathyXY(v, grid, zi=level, timeDep=False)

			print(np.max(u**2+v**2)**0.5)

			bathy[row].append(ptt.maskBathyXY(grid.bathy, grid, zi=0, timeDep=False))
			uvec[row].append(u[::d, ::d])
			vvec[row].append(v[::d, ::d])

	#==

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	print(Z[level])
	
	Xd = X[::d]; Yd = Y[::d]
	
	title0 = '(a) Zonally uniform slope'
	titlea = '(b) Eastern trough'
	titleb = '(c) Eastern sill'
	titlec = '(d) Eastern trough and sill'
	title = [[title0, titlea], [titleb, titlec]]
	
	cbar = False
	
	xlabel = 'Lon (km)'
	xlabels = [[None]*2, [xlabel]*2]
	ylabel = 'Lat (km)'
	ylabels = [[ylabel, None]]*2
	
	lonticks = [100,200,300,400,500,500]#[0, 200, 400, 600]
	latticks = [0, 100, 200, 300, 400, 500]
	xticks = [[lonticks, lonticks]]*2
	yticks = [[latticks, latticks]]*2
	xticksvis = [[False]*2, [True]*2]
	yticksvis = [[True, False]]*2
	
	figsize = (7, 6)
	fontsize = 12
	width_ratios = [1,1]

	fontdict = {'fontsize':12, 'color':'k'}
	text = {'text':'z = '+str(Z[level])+' m', 'xloc':X[1], 'yloc':Y[-22], 'fontdict':fontdict}
	text_data = [[text, text], [text, text]]

	cbdata = [[0.825, 0.15, 0.015, 0.7], None]

	pt.quiver2by2(uvec, vvec, Xd, Yd, contourf=bathy, X=X, Y=Y, mesh=True, cmap='YlGn', vmin=-1000, vmax=-400, cbar=cbar, cbarShared=True, cbarSharedData=cbdata, figsize=figsize, xticks=xticks, xticksvis=xticksvis, yticks=yticks, yticksvis=yticksvis, save=True, ylabels=ylabels, xlabels=xlabels, title=title, fontsize=fontsize, text_data=text_data, width_ratios=width_ratios, scale=1)

	#pt.quiver1byN(uvec, vvec, Xd, Yd, contourf=bathy, X=X, Y=Y, mesh=True, cmap='YlGn', vmin=-1000, vmax=-400, cbar=cbar, figsize=figsize, xticks=xticks, xticksvis=xticksvis, yticks=yticks, yticksvis=yticksvis, save=True, ylabel=ylabel, xlabel=xlabel, title=title, fontsize=fontsize, text_data=text_data, width_ratios=width_ratios, scale=1)
	
	quit()

#==

FIGURE8 = 0
if FIGURE8:

	levels = [0,22]	
	ts = 108
	
	path = '/home/michael/Documents/data/MCS_114/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy
	vmin = -1000; vmax = -300
	bathy = tools.boundData(bathy, vmin, vmax)
	bathy = [bathy, bathy]
	
	X = grid.XC/1.e3; Y = grid.YC/1.e3
	Z = grid.RC.squeeze()
	
	# Load data
	u = np.mean(readVariable('UVEL', path, meta=False)[ts:,], axis=0)
	v = np.mean(readVariable('VVEL', path, meta=False)[ts:,], axis=0)
	T = np.mean(readVariable('THETA', path, meta=False)[ts:,], axis=0)
	
	# Get two levels
	u = u[levels]
	v = v[levels]
	T = T[levels]

	cvmin = -1.6; cvmax = 1.6
	T = tools.boundData(T, cvmin, cvmax, 0.9999)

	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')
	u = tools.boundData(u, -0.4, 0.4, 0.9999); v = tools.boundData(v, -0.4, 0.4, 0.9999)


	# Mask data
	for li in range(len(levels)):
		u[li] = ptt.maskBathyXY(u[li], grid, levels[li], timeDep=False)
		v[li] = ptt.maskBathyXY(v[li], grid, levels[li], timeDep=False)
		T[li] = ptt.maskBathyXY(T[li], grid, levels[li], timeDep=False)
		u[li] = ptt.maskDraftXY(u[li], grid, levels[li], timeDep=False)
		v[li] = ptt.maskDraftXY(v[li], grid, levels[li], timeDep=False)
		T[li] = ptt.maskDraftXY(T[li], grid, levels[li], timeDep=False)
		bathy[li] = ptt.maskBathyXY(bathy[li], grid, levels[li], timeDep=False)
		bathy[li] = ptt.maskDraftXY(bathy[0], grid, levels[0], timeDep=False)

	#==
	
	# Text data	
	fontdict = {'fontsize':12, 'color':'k'}
	text0 = {'text':'z = '+str(Z[levels[0]])+' m', 'xloc':X[-17,0], 'yloc':Y[-17,0], 'fontdict':fontdict}
	text1 = {'text':'z = '+str(Z[levels[1]])+' m', 'xloc':X[-17,0], 'yloc':Y[-17,0], 'fontdict':fontdict} 
	text_data = [text0, text1]
	
	# Sample rate
	d = 8
	u = u[..., ::d, ::d]; v = v[..., ::d, ::d]; T = T[..., ::d, ::d]
	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd, Xd]; Yd = [Yd, Yd]
	X = [X, X]; Y = [Y, Y]

	T[0,-1,-1] = cvmin; T[0,-1,-2] = cvmax;
	T[1,-1,-1] = cvmin; T[1,-1,-2] = cvmax;

	vmin = [vmin, vmin]; vmax = [vmax, vmax]

	title1 = '(a) Surface flow, Pot. Temp. (deg. C), Bathy. (m)'
	title2 = '(b) Deep flow, Pot. Temp. (deg. C), Bathy. (m)'
	title = [title1, title2]


	lonticks = [0, 200, 400, 600]#[100,200,300,400,500,500]
	latticks = [0, 100, 200, 300, 400, 500]
	xticks = [lonticks, lonticks]
	yticks = [latticks, latticks]
	xticksvis = [True, True]
	yticksvis = [True, False]
	
	xlabel = 'Lon (km)'
	xlabels = [xlabel, xlabel]
	ylabel = 'Lat (km)'
	ylabels = [ylabel, None]

	scale = [2,1]
	width_ratios = [0.8,1.25]
	
	# PLOT

	pt.quiver1by2(u, v, Xd, Yd, C=T, ccmap='seismic', contourf=bathy, X=X, Y=Y, mesh=True, vmin=vmin, vmax=vmax, cmap='YlGn', width_ratios=width_ratios, title=title, fontsize=12, figsize=(8., 3), xticks=xticks, yticks=yticks, xticksvis=xticksvis, yticksvis=yticksvis, xlabels=xlabels, ylabels=ylabels, scale=scale, save=True, text_data=text_data)

	quit()

#==

# A sequence of isotherm heights from different experiments.
FIGURE9 = False
if FIGURE9:

	path = '/home/michael/Documents/data/'	
	path_THERMZ = path + 'THERMZnpy/'
	#runs = [['MCS_108', 'MCS_120', 'MCS_116'], ['MCS_117', 'MCS_118', 'MCS_114']]
	runs =[['MCS_142', 'MCS_135', 'MCS_132'], ['MCS_133', 'MCS_136', 'MCS_137']]
	titles = [['(a) Uniform shelf', '(b) W', '(c) E'], ['(d) S', '(e) E+S', '(f) W+C+E+S']]
	#HCs = [[1.18, 1.37, 1.41], [2.05, 2.01, 2.23]]
	HCs = [[0.728, 0.920, 0.961], [1.61, 1.56, 1.73]]

	# Heat contents
	# Experiment -- Total heat content -- HC minus initial condition.
	#MCS_108 1.1773118609991408e+20 7.2751814824914076e+19
	#MCS_114 2.2267912740532e+20 1.7331224870024774e+20
	#MCS_115 2.2279558036632755e+20 1.734287016612553e+20
	#MCS_116 1.410387010293125e+20 9.605932975431249e+19
	#MCS_117 2.0540254908284826e+20 1.607943812359065e+20
	#MCS_118 2.0096019693955062e+20 1.563433775096586e+20
	#MCS_119 1.9537381653424665e+20 1.507569971043546e+20
	#MCS_120 1.369755546766797e+20 9.19961834016797e+19

	THERM = -0.5; THERMt = 'm05'

	data = [[], []]
	hc = [[], []]
	text_data = [[], []]
	M = len(runs)
	N = len(runs[0])

	for row in range(M):
		for col in range(N):
		
			run = runs[row][col]
	
			grid = Grid(path+run+'/run/')
			X = grid.XC[1,:] / 1.e3
			Y = grid.YC[:,1] / 1.e3
	
			#ThermZ = np.load(path_THERMZ+'ThermZ_'+THERMt+'_'+run+'.npy')
			ThermZ = np.load(path+run+'/run/ThermZ_'+THERMt+'_'+run+'.npy')
			ThermZ = np.mean(ThermZ[-12:], axis=0)
								
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)
			ThermZ = np.where(ThermZ<grid.bathy, np.nan, ThermZ)

			data[row].append(ThermZ)

			txt = r'$\Delta$(HC) = ' + str(HCs[row][col])
			text = {'text':txt, 'xloc':X[1], 'yloc':Y[-20], 'fontdict':{'fontsize':12, 'color':'k'}}
			text_data[row].append(text)

	#==
	
	X = grid.XC[1,:] / 1.e3
	Y = grid.YC[:,1] / 1.e3
	
	hlines = [[240, None, None], [None]*3]
	
	xlabel = 'Lon (km)'
	xlabels = [[None]*3, [xlabel]*3]
	ylabel = 'Lat (km)'
	ylabels = [[ylabel, None, None]]*2
	
	lonticks = [0, 200, 400, 600]#[100,200,300,400,500,500]
	latticks = [0, 100, 200, 300, 400, 500]
	xticks = [[lonticks, lonticks, lonticks]]*2
	yticks = [[latticks, latticks, latticks]]*2
	xticksvis = [[False]*3, [True]*3]
	yticksvis = [[True, False, False]]*2
	
	cbdata = [[0.825, 0.15, 0.015, 0.7], '-0.5 deg. C isotherm depth (m)']
	
	vmin = -500; vmax=-200
	#vmin = -450; vmax=-250
	cmap = 'YlOrRd'#'jet'
	
	pt.plotMbyN(data, X=X, Y=Y, mesh=False, contourfNlevels=13, vmin=vmin, vmax=vmax, titles=titles, cmap=cmap, cbar=False, cbarShared=True, cbarSharedData=cbdata, xlabels=xlabels, ylabels=ylabels, xticks=xticks, yticks=yticks, xticksvis=xticksvis, yticksvis=yticksvis, width_ratios=[1,1,1], text_data=text_data, grid=True, figsize=(10,4), save=True, hlines=hlines)
	
	quit()
	
#==

# Cross-shelf heat transport plots.
FIGURE10 = True
if FIGURE10:

	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)
	rho0 = 1030. # Units kg / m3
	normval = 1.e12
	
	# First experiment defined here.
	path_root = '/home/michael/Documents/data/'
	path = path_root + 'MCS_132/run/'
	grid = Grid(path)

	# Subregions for heat transport.

	lat = 95

	# troughW
	lonsW = [100e3, 200e3]; depthW = [-10, -800]
	lonW = grid.getIndexFromLon(lonsW)
	depthW = grid.getIndexFromDepth(depthW)
	labelW = 'trough W'

	# troughE
	lonsE = [418e3, 520e3]; depthE = [-10, -800]
	lonE = grid.getIndexFromLon(lonsE)
	depthE = grid.getIndexFromDepth(depthE)
	labelE = 'trough E'

	# troughE AND sill
	lonsES = [415e3, 585e3]; depthES = [-10, -800]
	lonES = grid.getIndexFromLon(lonsES)
	depthES = grid.getIndexFromDepth(depthES)
	labelES = 'trough E + sill'

	labelAll = 'All'
	labelU = 'Uniform lons'

	# Define colours for different lon ranges.
	cAll = 'k'
	cE = 'r'
	cU = 'b'
	cES = 'g'
	cW = 'c' 

	#==
	
	# 1. 
	# Get heat transport for first (bathyE) experiment.

	T = readVariable('THETA', path, meta=False)[...,lat,:]
	Tf = -2. * np.ones(T.shape)
	v = tools.interp(readVariable('VVEL', path, meta=False)[...,lat,:], 'v')

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	T = rho0 * Cp * v * (T - Tf)
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)
	T *= area / normval
	
	# Heat transport for troughE, uniform lons and all lons.
	TE = np.ma.sum(T[:, depthE[0]:depthE[1], lonE[0]:lonE[1]], axis=(1,2))
	TAll = np.ma.sum(T, axis=(1,2))
	TU = TAll - TE
	
	TE = tools.smooth3(TE)
	TAll = tools.smooth3(TAll)
	TU = tools.smooth3(TU)
	
	Ts1 = [TE, TU, TAll]
	labels1 = [labelE, labelU, labelAll]
	colours1 = [cE, cU, cAll]
	title1 = '(a) E'

	#==

	# 2. 
	# Repeat for bathyES experiment.
	exp = 'MCS_141'#'MCS_136'
	path = path_root + exp + '/run/'
	grid = Grid(path)
	
	T = readVariable('THETA', path, meta=False)[...,lat,:]
	v = tools.interp(readVariable('VVEL', path, meta=False)[...,lat,:], 'v')

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	T = rho0 * Cp * v * (T - Tf)
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)
	T *= area / normval
	
	# Heat transport for troughE, uniform lons and all lons.
	TE = np.ma.sum(T[:, depthE[0]:depthE[1], lonE[0]:lonE[1]], axis=(1,2))
	TES = np.ma.sum(T[:, depthES[0]:depthES[1], lonES[0]:lonES[1]], axis=(1,2))
	TAll = np.ma.sum(T, axis=(1,2))
	TU = TAll - TES
	
	TE = tools.smooth3(TE)
	TES = tools.smooth3(TES)
	TAll = tools.smooth3(TAll)
	TU = tools.smooth3(TU)
	
	Ts2 = [TE, TES, TU, TAll]
	labels2 = [labelE, labelES, labelU, labelAll]
	colours2 = [cE, cES, cU, cAll]
	title2 = '(b) E+S'
	
	#==
	
	# 3. 
	# Repeat for bathyES experiment.
	exp = 'MCS_141'
	path = path_root + exp + '/run/'
	grid = Grid(path)
	
	T = readVariable('THETA', path, meta=False)[...,lat,:]
	v = tools.interp(readVariable('VVEL', path, meta=False)[...,lat,:], 'v')

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	T = rho0 * Cp * v * (T - Tf)
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)
	T *= area / normval
	
	# Heat transport for troughE, uniform lons and all lons.
	TW = np.ma.sum(T[:, depthW[0]:depthW[1], lonW[0]:lonW[1]], axis=(1,2))
	TE = np.ma.sum(T[:, depthE[0]:depthE[1], lonE[0]:lonE[1]], axis=(1,2))
	TES = np.ma.sum(T[:, depthES[0]:depthES[1], lonES[0]:lonES[1]], axis=(1,2))
	TAll = np.ma.sum(T, axis=(1,2))
	TU = TAll - TES - TW
	
	TE = tools.smooth3(TE)
	TW = tools.smooth3(TW)
	TES = tools.smooth3(TES)
	TAll = tools.smooth3(TAll)
	TU = tools.smooth3(TU)
	
	Ts3 = [TW, TE, TES, TU, TAll]
	labels3 = [labelW, labelE, labelES, labelU, labelAll]
	colours3 = [cW, cE, cES, cU, cAll]
	title3 = '(c) W+C+E+S'
	
	#==
	
	# Prepare plotting data.
	
	Ts = [Ts1, Ts2, Ts3]
	labels = [labels1, labels2, labels3]
	colours = [colours1, colours2, colours3]
	titles = [title1, title2, title3]
	
	ylabel = 'Heat transport (TW)'
	ylabels = [ylabel, ylabel, ylabel]
	
	xlabel = 'Time (months)'
	xlabels = [None, None, xlabel]
	
	xlims = [0, 240]
	ylims = [-3, 3]
		
	xticks = np.linspace(0, 240, 7)
	xticks = [xticks]*3
	xticksvis = [False, False, True]
	
	yticks = [-3, -2, -1, 0, 1, 2, 3]
	yticks = [yticks]*3
	yticksvis = [True, True, True]
	
	#==

    # NOW PLOT
	pt.line1by3(Ts, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis)

	quit()	




	



