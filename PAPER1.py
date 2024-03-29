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

	#rho0 = 1.3; Cd = 0.003
	rho0 = 1.2
	uMin = 0.5
		
	# These from MITgcm default params.
	cDrag_1 = 2.70e-3
	cDrag_2 = 0.142e-3
	cDrag_3 = 0.0764e-3
	exf_scal_BulkCdn = 1.015
	
	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	
	bathy = grid.bathy
	bathy = ptt.maskBathyXY(bathy, grid, zi=0)	
	bathy = tools.boundData(bathy, -5000, 0)

	lats = [-75.5, -70.5]; lons = [245, 262] 
	outline = [lons, lats]

	icel, icea = isf.get_ice_shelf_front(grid, grid.iceC, grid.landC, neighbourPts=1)
	#gll, gla = isf.get_grounding_line(grid, grid.iceC, grid.landC, neighbourPts=True)
		
	X = grid.XC; Y = grid.YC
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]

	taux = np.load(path+'tauxmean_PAS851.npy'); tauy = np.load(path+'tauymean_PAS851.npy')
	#windx = np.load(path+'uwindmean_PAS851.npy'); windy = np.load(path+'vwindmean_PAS851.npy')
	windx = np.load(path+'uwindstressmean_PAS851.npy')
	windy = np.load(path+'vwindstressmean_PAS851.npy')
	
	taux = tools.interp(taux, 'u'); tauy = tools.interp(tauy, 'v')
	#windx = tools.interp(windx, 'u'); windy = tools.interp(windy, 'v')
	
	#speed = (windx**2 + windy**2)**0.5
	#wsm = np.maximum(speed, uMin)
	#Cd = exf_scal_BulkCdn * (cDrag_1 / wsm + cDrag_2 + cDrag_3 * wsm)
	#windx = rho0 * Cd * windx# * speed
	#windy = rho0 * Cd * windy# * speed

	windx = ptt.maskBathyXY(windx, grid, zi=0); windy = ptt.maskBathyXY(windy, grid, zi=0)
	windx = ptt.maskDraftXY(windx, grid, zi=0); windy = ptt.maskDraftXY(windy, grid, zi=0)	
	taux = ptt.maskBathyXY(taux, grid, zi=0); tauy = ptt.maskBathyXY(tauy, grid, zi=0)	
	taux = ptt.maskDraftXY(taux, grid, zi=0); tauy = ptt.maskDraftXY(tauy, grid, zi=0)	
		
	vecx = [windx, taux]; vecy = [windy, tauy]
	vecx = taux; vecy = tauy

	# Sample rate
	d = 20

	paras = [-74, -72, -70, -68, -66, -64]
	merids = [230, 240, 250, 260, 270]

	xlabel = 'LON'; ylabel = 'LAT'
	title = r'Surface stress (N m$^{-2}$), bathy. (m)'

	pt.quiver1by1Basemap(vecx, vecy, X, Y, d, lat_0, lon_0, contourf=bathy, mesh=True, cmap='plasma', parallels=paras, meridians=merids, isf=icea, outline=outline, figsize=(5,3), title=title, save=True, outname='Figure_1.png')

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
	
	iceC = grid.iceC
	iceC = np.where(iceC==0, np.nan, iceC)
	icel, icea = isf.get_ice_shelf_front(grid, grid.iceC, grid.landC, neighbourPts=False)
		
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
	iceC = tools.getSubregionXY(iceC, latsi, lonsi)
	icea = tools.getSubregionXY(icea, latsi, lonsi)
			
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

	fontdict0 = {'fontsize':12, 'color':'k'}
	fontdictIce = {'fontsize':7, 'color':'w'}
	xloc0 = 5.e4; yloc0 = 1.78e6

	xloc = [xloc0, 1.55e6, 0.9e6, 4.6e5, 1.85e5, 1.52e6, 1.4e6]
	yloc = [yloc0, 1.9e5, 2.e5, 1.95e5, 3.85e5, 8.25e5, 1.2e6]
	text = ['z = '+str(Z[levels[0]])+' m', 'PIG', 'THW', 'CRO', 'DOT', 'COS', 'ABB']
	fontdict = [fontdict0, fontdictIce, fontdictIce, fontdictIce, fontdictIce, fontdictIce, fontdictIce]
	
	text0 = {'text':text, 'xloc':xloc, 'yloc':yloc, 'fontdict':fontdict}
	text1 = {'text':'z = '+str(Z[levels[1]])+' m', 'xloc':xloc0, 'yloc':yloc0, 'fontdict':fontdict0} 
	text_data = [text0, text1]

	# Labels for troughs
	p1 = {'x':(1.91e5,1.536e6), 't':'PITW', 'tx':(1.915e5,1.536e6)}
	#p2 = {'x':(1.025e6, 6.e5), 't':'C'}
	p2 = {'x':(9.e5, 8.e5), 't':'C', 'tx':(9.002e5, 8.e5)}
	#p3 = {'x':(1.37e6, 1.56e6), 't':'E'}
	p3 = {'x':(1.32e6, 1.52e6), 't':'PITE', 'tx':(1.03e6, 1.52e6)}
	p4 = {'x':(1.527e6, 1.612e6), 't':'R', 'tx':(1.53e6, 1.612e6)}
	labelData = [[p1, p2, p3, p4], []]

	# PLOT

	pt.quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, C=T, ccmap='seismic', contourf=bathy, X=X, Y=Y, mesh=False, contourfNlevels=13, vmin=vmin, vmax=vmax, cmap='YlGn', parallels=paras, meridians=merids, isf=[iceC,None], width_ratios=[0.95,1.2], title=title, fontsize=8, figsize=(7.7, 3), text_data=text_data, labelData=labelData, save=True, outname='Figure_2.png')

	quit()
	
#==

# Plot of isotherm height and ASF in PAS.
FIGURE3 = 0
if FIGURE3:

	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy

	iceC = grid.iceC
	iceC = np.where(iceC==0, np.nan, iceC)
	icel, icea = isf.get_ice_shelf_front(grid, grid.iceC, grid.landC, neighbourPts=2)

	T = np.load(path+'Tmean_PAS851.npy')
	S = np.load(path+'Smean_PAS851.npy')

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

	iceC = tools.getSubregionXY(iceC, latsi, lonsi)
	icea = tools.getSubregionXY(icea, latsi, lonsi)
	
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
	S = tools.getSubregionYZ(S[...,xi], latsi2, zi)
	T = ptt.maskBathyYZ(T, grid, xi=xi, subregion=True, lats=latsi2, depths=zi)
	S = ptt.maskBathyYZ(S, grid, xi=xi, subregion=True, lats=latsi2, depths=zi)
	T = ptt.maskDraftYZ(T, grid, xi=xi, subregion=True, lats=latsi2, depths=zi)
	S = ptt.maskDraftYZ(S, grid, xi=xi, subregion=True, lats=latsi2, depths=zi)
	T = tools.boundData(T, vmin2, vmax2, scale=0.9999)

	#==

	vmin = [vmin1, vmin2]; vmax = [vmax1, vmax2]	

	paras = [-75, -74, -73, -72, -71]
	merids = [245, 250, 255, 260]
	
	contourlevels = [[-1000, -600],[34.3, 34.5, 34.7]]
	lws = [0.2, 1.]
	
	bathy = tools.getSubregionXY(bathy, latsi, lonsi)
	X, Y = grid.XYsubr(lons, lats)
	Xp = [X, Ysubr]; Yp = [Y, Zsubr]

	xlabels = [None, 'Lat (deg.)']; ylabels = [None, 'Depth (m)']
	titles = ['(a) -0.5 deg. isotherm depth (m)', '(b) Pot. Temp. (deg. C); Salinity']

	pt.plot1by2Basemap1([ThermZ, T], Xp, Yp, lat_0, lon_0, mesh=False, contour=[bathy,S], contourlevels=contourlevels, lws=lws, contourfNlevels=13, vmin=vmin, vmax=vmax, parallels=paras, meridians=merids, isf=[icea, None], yline=[slice_lon, None], xlabels=xlabels, ylabels=ylabels, fontsize=11, titles=titles, cmaps=['YlOrRd', 'seismic'])

	quit()
	
#==

# First figure of idealised model. Panel 1: bathymetry, wind forcing. Panel 2: T/S relaxation profiles.
FIGURE4 = 0
if FIGURE4:

	# MCS_104 for run with no rel. in south.
	path = '/home/michael/Documents/data/MCS_145/run/'
	grid = Grid(path)
	bathy = grid.bathy

	Trel, Srel = tools.getTSrel(grid, salttoptop=33.5, salttop=33.5, saltbot=34.5, temptop=-1.8, tempbot=1.0, tclinewid=140., tclinetop=-200, hclinewid=140, hclinetop=-200, shift=0)
	
	Nz = Trel.shape[0]
	
	TrelS = -1.8*np.ones(Nz); SrelS = 33.5*np.ones(Nz)
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
	tauxTr = 300+200*tools.getTranslatedWind(grid.Nx,grid.Ny,16)[:,0]
		
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

	fig = 	plt.figure(figsize=(8.8,4), dpi=200)
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
	plt.plot(tauxTr, Y, color='k', linestyle='dotted')
	
	plt.title('(a) Bathymetry (m), Wind stress')

	#==

	ax = fig.add_subplot(gs[0,1])
	
	ln1n = ax.plot(Trel, Z, color='r', label='Pot. Temp. (deg. C), North')
	ln1s = ax.plot(TrelS, Z, color='r', alpha=0.5, linestyle='--', label='Pot. Temp. (deg. C), South')
	ax.tick_params(axis='x', labelcolor='r')
	ax.set_xlim(-2.0, 1.2)
	ax.set_ylabel('Depth (m)', fontsize=fs)
	ax2 = ax.twiny()
	
	#ax2.plot(Srel, Z, color='r', label='Pot. Temp.')
	ln2n = ax2.plot(Srel, Z, color='k', label='Salinity (g/kg), North')
	ln2s = ax2.plot(SrelS, Z, color='k', alpha=0.5, linestyle='--', label='Salinity (g/kg), South')
	ax2.tick_params(axis='x', labelcolor='k')

	plt.title('(b) Relaxation profiles')

	lns = ln1n+ln1s+ln2n+ln2s
	labs = [l.get_label() for l in lns]
	ax.legend(lns, labs, loc=8, prop={'size': 8})

	plt.tight_layout()
	#ax2.plot(1, 1, ' ', color='r', label="Extra label on the legend")

	plt.savefig('Figure_4.png')
	#plt.show()
	
	quit()	

#==

# Zonal-mean zonal velocity, temperature and SSH.
FIGURE5 = 0
if FIGURE5:

	# MCS_104 for run with no rel. in south.
	path = '/home/michael/Documents/data/MCS_145/run/'
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

	fig = plt.figure(figsize=(4,4.5), dpi=200, constrained_layout=True)
	gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=height_ratios)

	ax = fig.add_subplot(gs[0,0])

	plt.plot(Y[1:-1], h[1:-1], color='k')
	
	ax.axes.xaxis.set_ticklabels([])
	plt.ylabel('SSH (m)', fontsize=fs)

	plt.title('(a) Zonal-mean SSH anomaly', fontsize=fs2)
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

	plt.savefig('Figure_5.png')
	plt.show()
	
	quit()

#==

# New walled simulation. Show SSH with bathy, Uvel at two longitudes.
FIGURE6 = 0
if FIGURE6:

	path = '/home/michael/Documents/data/MCS_144/run/'
	grid = Grid(path)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	xis = [0, 120]

	# These for snapshot at end of simulation
	#ts = -12;#116
	#h = readVariable('ETAN', path, meta=False)[ts:]
	#u = readVariable('UVEL', path, meta=False)[ts:][...,xis]
	#T = readVariable('THETA', path, meta=False)[ts:][...,xis]

	# These for mean over last ts months
	ts = -12
	u = np.mean(readVariable('UVEL', path, meta=False)[ts:,...,xis], axis=0)
	T = np.mean(readVariable('THETA', path, meta=False)[ts:,...,xis], axis=0)
	h = np.mean(readVariable('ETAN', path, meta=False)[ts:], axis=0)

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
	cmaps = ['viridis', 'coolwarm', 'coolwarm']
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
	#runs = [['MCS_108', 'MCS_116'], ['MCS_117', 'MCS_118']]
	runs = [['MCS_144', 'MCS_132'], ['MCS_133', 'MCS_136']]
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

			ts = 228
			u = np.mean(readVariable('UVEL', path, meta=False)[ts:, level], axis=0)
			v = np.mean(readVariable('VVEL', path, meta=False)[ts:, level], axis=0)
			
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
	ts = -24; te = -12
	
	path = '/home/michael/Documents/data/MCS_141/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy
	vmin = -1000; vmax = -400
	bathy = tools.boundData(bathy, vmin, vmax)
	bathy = [bathy, bathy]
	
	X = grid.XC/1.e3; Y = grid.YC/1.e3
	Z = grid.RC.squeeze()
	
	# Load data
	u = np.mean(readVariable('UVEL', path, meta=False)[ts:te,], axis=0)
	v = np.mean(readVariable('VVEL', path, meta=False)[ts:te,], axis=0)
	T = np.mean(readVariable('THETA', path, meta=False)[ts:te,], axis=0)
	
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

	pt.quiver1by2(u, v, Xd, Yd, qs=[0.2,0.1], C=T, ccmap='seismic', contourf=bathy, X=X, Y=Y, mesh=True, vmin=vmin, vmax=vmax, cmap='YlGn', width_ratios=width_ratios, title=title, fontsize=12, figsize=(8., 3), xticks=xticks, yticks=yticks, xticksvis=xticksvis, yticksvis=yticksvis, xlabels=xlabels, ylabels=ylabels, scale=scale, save=True, text_data=text_data)

	quit()

#==

# A sequence of isotherm heights from different experiments.
FIGURE9 = 0
if FIGURE9:

	path = '/home/michael/Documents/data/'	
	path_THERMZ = path + 'THERMZnpy/'
		
	runs =[['MCS_159', 'MCS_157', 'MCS_158'], ['MCS_154', 'MCS_155', 'MCS_153']]
	HCs = [[1.021, 1.078, 1.063], [1.507, 1.407, 1.935]]
	
	titles = [['(a) Uniform shelf', '(b) W', '(c) E'], ['(d) S', '(e) E+S', '(f) W+C+E+S']]

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

			txt = r'HC = ' + str(HCs[row][col])
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
FIGURE10 = 0
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
	# Repeat for bathyWCES2 experiment.
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

	# 4. 
	# And for bathyWCES2, 05cos_tr16 experiment.
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
	
	Ts4 = [TW, TE, TES, TU, TAll]
	labels4 = [labelW, labelE, labelES, labelU, labelAll]
	colours4 = [cW, cE, cES, cU, cAll]
	title4 = '(d) W+C+E+S, Shifted Wind'
	
	#==
		
	# Prepare plotting data.
	
	Ts = [Ts1, Ts2, Ts3, Ts4]
	labels = [labels1, labels2, labels3, labels4]
	colours = [colours1, colours2, colours3, colours4]
	titles = [title1, title2, title3, title4]
	
	ylabel = 'Heat transport (TW)'
	ylabels = [ylabel, ylabel, ylabel, ylabel]
	
	xlabel = 'Time (months)'
	xlabels = [None, None, None, xlabel]
	
	xlims = [0, 360]
	ylims = [-3, 3]
		
	xticks = np.linspace(0, 360, 10)
	xticks = [xticks]*4
	xticksvis = [False, False, False, True]
	
	yticks = [-3, -2, -1, 0, 1, 2, 3]
	yticks = [yticks]*4
	yticksvis = [True, True, True, True]
	
	#==

    # NOW PLOT
	pt.line1byN(Ts, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis)

	quit()	

#==

FIGURE11 = 0
if FIGURE11:

	levels = [0,22]	
	ts = -12; te = -1
	
	path_root = '/data/oceans_output/shelf/michai/mitgcm/'
	#path = '/home/michael/Documents/data/'
	
	path = path_root + 'MCS_156/run/'
	
	grid = Grid_PAS(path)
	bathy = grid.bathy
	vmin = -1000; vmax = -400
	bathy = tools.boundData(bathy, vmin, vmax)
	bathy = [bathy, bathy]
	
	X = grid.XC/1.e3; Y = grid.YC/1.e3
	Z = grid.RC.squeeze()
	
	# Load data
	u = np.mean(readVariable('UVEL', path, meta=False)[ts:,], axis=0)
	v = np.mean(readVariable('VVEL', path, meta=False)[ts:,], axis=0)
	T = np.mean(readVariable('THETA', path, meta=False)[ts:,], axis=0)
	
	DIFF = True
	if DIFF:
		path = path_root + 'MCS_153/run/'
		u -= np.mean(readVariable('UVEL', path, meta=False)[ts:,], axis=0)
		v -= np.mean(readVariable('VVEL', path, meta=False)[ts:,], axis=0)
		#T -= np.mean(readVariable('THETA', path, meta=False)[ts:,], axis=0)
		T = 2 * np.ones(u.shape)
	
	# Get two levels
	u = u[levels]
	v = v[levels]
	T = T[levels]

	cvmin = -0.1; cvmax = 0.1
	T = tools.boundData(T, cvmin, cvmax, 0.9999)

	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')
	u = tools.boundData(u, -0.1, 0.1, 0.9999); v = tools.boundData(v, -0.1, 0.1, 0.9999)

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

	scale = [.1,.1]
	width_ratios = [0.8,1.25]
	
	# PLOT

	pt.quiver1by2(u, v, Xd, Yd, C=T, ccmap='seismic', contourf=bathy, X=X, Y=Y, mesh=True, vmin=vmin, vmax=vmax, cmap='YlGn', width_ratios=width_ratios, title=title, fontsize=12, figsize=(8., 3), xticks=xticks, yticks=yticks, xticksvis=xticksvis, yticksvis=yticksvis, xlabels=xlabels, ylabels=ylabels, scale=scale, save=True, outname='Figure_11.png', text_data=text_data)

	quit()	
	
#==

# A sequence of isotherm heights from different experiments.
FIGURE12 = 0
if FIGURE12:

	path = '/home/michael/Documents/data/'	
	path_THERMZ = path + 'THERMZnpy/'
	runs = [['MCS_108', 'MCS_120', 'MCS_116'], ['MCS_117', 'MCS_118', 'MCS_160']]
	HCs = [[0.728, 0.920, 0.961], [1.607, 1.561, 1.646]]
	
	titles = [['(a) Uniform shelf', '(b) W', '(c) E'], ['(d) S', '(e) E+S', '(f) W+C+E+S']]

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
			ThermZ = np.mean(ThermZ, axis=0)
								
			ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)
			ThermZ = np.where(ThermZ<grid.bathy, np.nan, ThermZ)

			data[row].append(ThermZ)

			txt = r'HC = ' + str(HCs[row][col])
			text = {'text':txt, 'xloc':X[1], 'yloc':Y[-20], 'fontdict':{'fontsize':12, 'color':'k'}}
			text_data[row].append(text)

	#==
	
	X = grid.XC[1,:] / 1.e3
	Y = grid.YC[:,1] / 1.e3
	
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
	
	pt.plotMbyN(data, X=X, Y=Y, mesh=False, contourfNlevels=13, vmin=vmin, vmax=vmax, titles=titles, cmap=cmap, cbar=False, cbarShared=True, cbarSharedData=cbdata, xlabels=xlabels, ylabels=ylabels, xticks=xticks, yticks=yticks, xticksvis=xticksvis, yticksvis=yticksvis, width_ratios=[1,1,1], text_data=text_data, grid=True, figsize=(10,4), save=True, outname='Figure_12.png')
	
	quit()

#==

# Across-shelf heat flux for some simulations.
FIGURE13 = 0
if FIGURE13:

	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)
	rho0 = 1030. # Units kg / m3
	normval = 1.e12
	T0 = -1.8

	ts = -2
	norm = 1.e5
	maxDepth = -700
	
	# First experiment defined here.

	paths = ['MCS_132', 'MCS_132']
	path_root = '/home/michael/Documents/data/'
	vT = []

	lat = 95
	
	#==
	
	# 1. 
	# Get heat transport for first (bathyE) experiment.

	path = path_root + paths[0] + '/run/'
	grid = Grid(path)
	zb = grid.getIndexFromDepth(maxDepth)

	T = readVariable('THETA', path, meta=False)[ts:,:,lat,:]
	Tf = T0 * np.ones(T.shape)
	#v = tools.interp(readVariable('VVEL', path, meta=False)[ts:,:,lat,:], 'v')

	T = np.mean(rho0 * Cp * 1 * (T - Tf), axis=0) / norm
	vT.append(ptt.maskBathyXZ(T, grid, yi=lat, timeDep=False)[:zb])


	#==

	vT.append(vT[0])	
	#path = path_root + paths[1] + '/run/'
	#grid = Grid(path)

	#T = readVariable('THETA', path, meta=False)[ts:,:,lat,:]
	#Tf = T0 * np.ones(T.shape)
	#v = tools.interp(readVariable('VVEL', path, meta=False)[ts:,:,lat,:], 'v')

	#T = np.mean(rho0 * Cp * v * (T - Tf), axis=0) / norm
	#T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=False)
	
	#==
	
	figsize = (7,3)
	outname = 'Figure_11.png'
	fontsize = 12
	
	cmap = 'coolwarm'
	vmin = -4; vmax = -vmin
	
	X = grid.XC[0,:] / 1.e3
	Y = grid.RC.squeeze()[:zb]
	X = [X, X]; Y = [Y, Y]
	
	cbar = [False, True]
	cbarLabel = [None, 'Cross shelf heat flux']
	width_ratios = [1, 1.2]
	
	xticks = [0, 200, 400, 600]
	xticks = [xticks, xticks]
	yticks = [-600, -500, -400, -300, -200, -100]
	yticks = [yticks, yticks]
	yticksvis = [True, False]
	
	xlabels = ['Lon (km)', 'Lon (km)']
	ylabels = ['Depth (m)', None]
	titles = ['(a) Original wind', '(b) Shifted wind']
	
	
	pt.plot1by2(vT, X=X, Y=Y, figsize=figsize, fontsize=fontsize, width_ratios=width_ratios, mesh=False, titles=titles, cmaps=cmap, vmin=vmin, vmax=vmax, xlabels=xlabels, ylabels=ylabels, xticks=xticks, yticks=yticks, yticksvis=yticksvis, cbar=cbar, cbarLabel=cbarLabel, show=True, save=True, outname=outname)
	
	#==

	#quit()



