# IPO_FIGS.py



# Scripts for producing figures for BAS IPO paper.



#==========================================================



import numpy as np

from scipy.stats import pearsonr



from grid import Grid

from grid_PAS import Grid as Grid_PAS



import matplotlib.pyplot as plt

import matplotlib.colors as cl



import plotting as pt

import plotting_tools as ptt



import tools

import PAS_tools



from readData import readVariable

from varDict import getPlottingVars, getTrefSref



import iceShelfFront as isf



import time



#==========================================================



DATADIR = '/home/michael/Documents/data/'

PASDIR = DATADIR + 'PAS_851/run/'

	

# Lats for Amundsen Sea subregion.

EASlons = [235, 262] 

EASlats = [-75.5, -70]



# Start/end month in PAS

ts_shift = -2; te_shift = -3

t_start = 24*12 + 5*12 + ts_shift; t_end = -11 + te_shift

IPO_t_start = 14*12+1 + ts_shift; IPO_t_end = te_shift

	

# These are indices fcorresponding to longitude of bathy features, in coords of undercurrent.

westGetz = [1, 60]

Getz = [60, 65]

westPITW = [65, 105]

PITW = [105, 120]

westPITE = [120, 180]

PITE = [180, 230]

eastPITE = [230, 250]

sections = {'westGetz':westGetz, 'Getz':Getz, 'westPITW':westPITW, 'PITW':PITW, 'westPITE':westPITE, 'PITE':PITE, 'eastPITE':eastPITE}



# Size of running-mean window length in months

nn = 60



# Plot which figures?

FIGURE1 = 0

FIGURE2 = 0

FIGURE3 = 1



#==



if FIGURE1:



	level = 16



	grid = Grid_PAS(PASDIR)

	

	bathy = grid.bathy



	iceC = grid.iceC

	iceC = np.where(iceC==0, np.nan, iceC)

	icel, icea = isf.get_ice_shelf_front(grid, grid.iceC, grid.landC, neighbourPts=False)



	X = grid.XC; Y = grid.YC

	Xl = X[0,:]; Yl = Y[:,0]

	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2

	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]



	lats = EASlats; lons = EASlons

	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)



	X, Y = grid.XYsubr(lons, lats)

	Z = grid.RC.squeeze()



	# Load time-mean velocities & temp.

	u = np.load(path+'umean_PAS851_z'+str(level)+'.npy')

	v = np.load(path+'vmean_PAS851_z'+str(level)+'.npy')

	T = np.load(path+'Tmean_PAS851_z'+str(level)+'.npy')



	vmin = -2.; vmax = -vmin

	T = tools.boundData(T, vmin, vmax, 0.9999)



	u = tools.interp(u, 'u'); v = tools.interp(v, 'v')

	lim = 0.1

	u = tools.boundData(u, -lim, lim, 0.9999); v = tools.boundData(v, -lim, lim, 0.9999)



	# Get subregion and mask

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	u = tools.getSubregionXY(u, latsi, lonsi)

	v = tools.getSubregionXY(v, latsi, lonsi)

	T = tools.getSubregionXY(T, latsi, lonsi)

	iceC = tools.getSubregionXY(iceC, latsi, lonsi)

	icea = tools.getSubregionXY(icea, latsi, lonsi)

	

	#iceC = np.where(bathy<Z[level],iceC,0)

	#iceC = np.ma.array(iceC, mask=bathy==0)

	

	u = ptt.maskDraftXY(u, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	v = ptt.maskDraftXY(v, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskDraftXY(T, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	u = ptt.maskBathyXY(u, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	v = ptt.maskBathyXY(v, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskBathyXY(T, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	bathy = ptt.maskBathyXY(bathy, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	bathy = ptt.maskDraftXY(bathy, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	bathyC = ptt.makeBathyContour(bathy, X, Y)

	

	# Sample rate

	d = 6



	T[-10,0] = vmin; T[-11,1] = vmax;



	paras = [-75,  -73,  -71]

	merids = [240, 250, 260]

	

	title = 'Flow & Pot. Temp. (deg. C) at Z = -455 m'



	# List for ice shelves

	r = 1.1e6

	xloc1 = [2.66e5, 1.55e6+r, 0.9e6+r, 4.6e5+r, 1.85e5+r, 1.52e6+r, 1.4e6+r]

	yloc1 = [4.e5, 1.9e5, 2.e5, 1.95e5, 3.85e5, 8.25e5, 1.2e6]

	text1 = ['GETZ', 'PIG', 'THW', 'CRO', 'DOT', 'COS', 'ABB']

	fontdictIce = {'fontsize':7, 'color':'w'}

	fontdict1 = [fontdictIce]*7

	

	# List for contour labels

	xloc2 = []

	yloc2 = []

	text2 = []

	fontdict2 = []

	

	text0 = {'text':text1+text2, 'xloc':xloc1+xloc2, 'yloc':yloc1+yloc2, 'fontdict':fontdict1+fontdict2}

	text_data = text0



	# PLOT



	pt.quiver1by1Basemap(u, v, X, Y, d, lat_0, lon_0, contourf=T, mesh=True, contourfNlevels=11, vmin=vmin, vmax=vmax, cmap='coolwarm', contour=[-bathyC,-bathy], contourlevels=[[1000], [500]], contourColours=['silver','white'], parallels=paras, meridians=merids, isf=iceC, title=title, text_data=text_data, fontsize=7, figsize=(5, 4), show=True, save=False, outname='Figure_1.png')



#==



if FIGURE2:



	from TPI import TPI_unfiltered_PSL as IPO

	IPO = np.array(IPO[IPO_t_start:IPO_t_end])

	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

	IPO /= 100*np.max(np.abs(IPO))

	

	grid = Grid_PAS(PASDIR)

	

	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'



	uvFile = 'slope_uv_max.npy'

	surfuvFile = 'surf_uv_av.npy'

	SHIfwFlxFile = 'SHIfwFlx.npy'

	timeFile = 'PAS_time.npy'

	

	slope_uv = np.load(path+uvFile)

	surf_uv = np.load(path+surfuvFile)

	melt = -np.load(path+SHIfwFlxFile) # With minus sign, so +ve into ocean.

	time = np.load(path+timeFile)

	

	PIbayLats = [-75.5, -74]; PIbayLons = [251, 262] # PI and Thwaites glaciers

	latsi = grid.getIndexFromLat(PIbayLats); lonsi = grid.getIndexFromLon(PIbayLons)

	melt = tools.getSubregionXY(melt[t_start:t_end], latsi, lonsi)

		

	rhoFresh = 1.e3; yearSeconds = 86400. * 365.25

	melt = ptt.maskDraftXY(melt, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi, inverse=True)

	melt = yearSeconds * np.ma.mean(melt, axis=(1,2)) / rhoFresh

	

	#===



	# Get along-slope currents.



	slope_uv = slope_uv[t_start:t_end]

	surf_uv = surf_uv[t_start:t_end]

	year = PAS_tools.getDecimalTime(time, t_start=t_start, t_end=t_end)

	

	#plotSections = ['westGetz', 'westPITW', 'westPITE']

	plotSections = ['westPITW', 'westPITE']



	uv_sum = np.zeros(len(year))

	surf_uv_sum = np.zeros(len(year))

	nx = 0

	for section in plotSections:

		iw = sections[section][0]; ie = sections[section][1]

		uv_sum += np.sum(slope_uv[:,iw:ie], axis=1)

		surf_uv_sum += np.sum(surf_uv[:,iw:ie], axis=1)

		nx += ie-iw

		

	uv_mean = uv_sum / nx

	surf_uv_mean = surf_uv_sum / nx



	#==



	# Detrend

	uv_mean = PAS_tools.detrend(uv_mean, year)

	surf_uv_mean = PAS_tools.detrend(surf_uv_mean, year)

	

	melt = PAS_tools.detrend(melt, year)

	#melt = PAS_tools.deseason(melt)

		

	uv_mean = PAS_tools.windowAv(uv_mean, n=nn)[nn//2:-nn//2+1]

	surf_uv_mean = PAS_tools.windowAv(surf_uv_mean, n=nn)[nn//2:-nn//2+1]

	melt = PAS_tools.windowAv(melt, n=nn)[nn//2:-nn//2+1]

	year = year[nn//2:-nn//2+1]

	

	plt.plot(year, uv_mean, color='r', label='Undercurrent (m/s)')

	plt.plot(year, surf_uv_mean, color='b', label='Surface current (m/s)')

	plt.plot(year, melt/100., color='darkgreen', label='PIG + THW basal melt (100 m/y)')

	plt.plot(year, IPO, color='k', linestyle='--', label='IPO')

	

	plt.xlim(min(year), max(year))

	plt.grid()

	plt.legend()

	plt.show()

	

#==



if FIGURE3:



	fnames = ['FWflx.npy', 'wk.npy']

	plotSections = ['westPITW', 'westPITE']

	

	#==

	

	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

		

	LATS = EASlats; LONS = EASlons

	LATS[1] = -68.; LONS[0] = 230.

	latsi = grid.getIndexFromLat(LATS); lonsi = grid.getIndexFromLon(LONS)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	X, Y = grid.XYsubr(LONS, LATS)

			

	Xl = X[0,:]; Yl = Y[:,0]

	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2

	lat_0 = Yl[Ny2]; lon_0 = Xl[Nx2]

	

	#==

	

	# Get along-slope velocity timeseries

	uvfile = 'slope_uv_max.npy'

	uv = np.load(path+uvfile)[t_start:t_end]

	

	uv_sum = np.zeros(uv.shape[0])	

	nx = 0

	for section in plotSections:

		iw = sections[section][0]; ie = sections[section][1]

		uv_sum += np.sum(uv[:,iw:ie], axis=1)

		nx += ie-iw		

	uv_mean = uv_sum / nx



	uv_mean = PAS_tools.detrend(uv_mean, None)

	uv_mean = PAS_tools.windowAv(uv_mean, n=nn)[nn//2:-nn//2+1]

	

	#==

	

	# Initialise correlations and pvals.

	nF = len(fnames); nY, nX = Y.shape

	corr = []

	pval = []

	

	for fi in range(nF):

	

		fname = fnames[fi]

		print(fname)

		

		corr_tmp = np.zeros((nY, nX))

		pval_tmp = np.zeros((nY, nX))

		

		data = np.load(path+fname)[t_start:t_end]

		data = tools.getSubregionXY(data, latsi, lonsi)

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

		#data = PAS_tools.detrendXY(data, None)

		data = PAS_tools.windowAv(data, n=nn)[nn//2:-nn//2+1]

		

		for j in range(nY):

			for i in range(nX):

				if data[0,j,i] != 0:

					corr_tmp[j,i], pval_tmp[j,i] = pearsonr(data[:,j,i], uv_mean)



		corr_tmp = ptt.maskBathyXY(corr_tmp, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		corr.append(ptt.maskDraftXY(corr_tmp, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi))

		

		pval_tmp = np.where(bathy>=0, 1, pval_tmp)

		pval.append(np.where(draft<0, 1, pval_tmp))

		

	#==

	

	bathyC = ptt.makeBathyContour(bathy, X, Y)

	

	cbarTicks = [[-1., -.5, 0., .5, 1.]]*2

	paras = [[-75,  -73,  -71]]*2

	merids = [[240, 250, 260]]*2

	

	pt.plot1by2Basemap(corr, X, Y, lat_0, lon_0, contourfNlevels=21, stippling=pval, stipData=[0.05, 4,4, .2], cmap='coolwarm', vmin=-1., vmax=1., contour=[bathyC,bathyC], contourlevels=[[-1000],[-1000]], cbarTicks=cbarTicks, landscape=False)



	

#==

	
