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

# Plot ERA5 winds, ocean surface stress and bathymetry in PAS.
FIGURE1a = 0
if FIGURE1a:

	rho0 = 1028.5; Cd = 0.006
	rho0 = 1.3; Cd= 0.006
	
	# These from MITgcm default params.
	cDrag_1 = 2.70e-3
	cDrag_2 = 0.142e-3
	cDrag_3 = 0.0764e-3
	#Cd = cDrag_1 + cDrag_2 + cDrag_3
	
	path = '/home/michai/Documents/data/PAS_8512/run/'
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

	pt.quiver1by1Basemap(vecx, vecy, X, Y, d, lat_0, lon_0, contourf=bathy, mesh=True, cmap='plasma', parallels=paras, meridians=merids, isf=icea, outline=outline, figsize=(5,3))

	quit()

#==

# Zoomed in bathymetry
FIGURE1b = 0
if FIGURE1b:

	path = '/home/michai/Documents/data/PAS_851/run/'
	#path = 	'/home/michai/Documents/data/PISOMIP_001/run/'
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

#==

#
FIGURE2 = 0
if FIGURE2:

	levels = [0,15]	

	path = '/home/michai/Documents/data/PAS_851/run/'
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

	# Get two levels
	u = u[levels]
	v = v[levels]
	T = T[levels]

	cvmin = -2; cvmax = 2
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

	title = ['(a) Surface flow, Pot. Temp., Bathymetry', '(b) Deep flow, Pot. Temp., Bathymetry']

	fontdict = {'fontsize':12, 'color':'k'}
	text0 = {'text':'Z = '+str(Z[levels[0]])+' m', 'xloc':5.e4, 'yloc':1.78e6, 'fontdict':fontdict}
	text1 = {'text':'Z = '+str(Z[levels[1]])+' m', 'xloc':5.e4, 'yloc':1.78e6, 'fontdict':fontdict} 
	text_data = [text0, text1]

	# Labels for troughs
	p1 = {'x':(1.91e5,1.536e6), 't':'1'}
	p2 = {'x':(1.025e6, 6.e5), 't':'2'}
	p3 = {'x':(1.386e6, 1.574e6), 't':'3'}
	labelData = [[p1, p2, p3],[]]

	# PLOT

	pt.quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, C=T, ccmap='seismic', contourf=bathy, X=X, Y=Y, mesh=False, contourfNlevels=13, vmin=vmin, vmax=vmax, cmap='YlGn', parallels=paras, meridians=merids, width_ratios=[0.95,1.2], title=title, fontsize=8, figsize=(7.7, 3), text_data=text_data, labelData=labelData)

#==

# Plot of isotherm height and ASF in PAS.
FIGURE3 = 1
if FIGURE3:

	path = '/home/michai/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	bathy = grid.bathy

	T = np.load(path+'Tmean_PAS851.npy')

	X = grid.XC; Y = grid.YC; Z = grid.RC.squeeze()
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]


	lats = [-75.5, -70.5]; lons = [245, 262]
	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	# FIRST PANEL: ISOTHERM HEIGHT

	vmin1 = -600; vmax1 = -200	

	THERM = -0.5
	ThermZ = tools.getIsothermHeight(T, THERM, grid, interp=True, timeDep=False)
	ThermZ = tools.getSubregionXY(ThermZ, latsi, lonsi)

	ThermZ = np.where(ThermZ!=ThermZ, 0, ThermZ)
	#zi = grid.getIndexFromDepth(-600)
	ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)
	#ThermZ = ptt.maskBathyXY(ThermZ, grid, zi, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	
	ThermZ = tools.boundData(ThermZ, vmin1, vmax1)

	# SECOND PANEL: THETA SLICE
	
	vmin2 = -2.; vmax2 = 2.

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

	xlabels = [None, 'LAT (deg.)']; ylabels = [None, 'DEPTH (m)']
	titles = ['(a) -0.5 deg. isotherm height', '(b) Pot. Temp.']

	pt.plot1by2Basemap1([ThermZ, T], Xp, Yp, lat_0, lon_0, mesh=False, contourfNlevels=17, vmin=vmin, vmax=vmax, parallels=paras, meridians=merids, yline=[slice_lon, None], xlabels=xlabels, ylabels=ylabels, fontsize=11, titles=titles, cmaps=['jet', 'seismic'])

#==

todo = 'in fig3 mask only where land is. Weird things happening with mask. use npwhere that ive added, recreate mask. make sure that where mask is (but where it shouldnt exist) the isotherm doesnt exist. check this using ASF in PAS_MCS.')
print('todo')





