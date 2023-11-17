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

import PAS_tools as ptools



from readData import readVariable, readnp

from varDict import getPlottingVars, getTrefSref



from TPI import getIPO, getSAM, getASL



import iceShelfFront as isf



import time



#==========================================================



DATADIR = '/home/michael/Documents/data/'

PASDIR = DATADIR + 'PAS_851/run/'

NPYDIR = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'	

	

black = 'k'; blue = 'cornflowerblue'; red = 'indianred'; green = 'seagreen'; grey = 'gray'

cmapMean = 'plasma'; cmapAnom = 'coolwarm'



# Density/mass constants

rho0 = 1028.5

rhoFresh = 1000.

gt = 1.e12



# Time constants

secs_per_year = 86400. * 365.25

secs_per_month = secs_per_year / 12.



# Lats for Amundsen Sea subregion.

EASlons = [235, 262] 

EASlats = [-75.5, -70]



LATS = [-75.5, -68]

LONS = [225, 262]



# Start/end month in PAS

ts_shift = -2#+24

te_shift = -3#-48

t_start = 24*12 + 5*12 + ts_shift; t_end = -11 + te_shift

IPO_t_start = 14*12+1 + ts_shift; IPO_t_end = te_shift



# Size of running-mean window length in months

nn = 60



# These are indices fcorresponding to longitude of bathy features, in coords of undercurrent.

westGetz = [1, 60]

Getz = [60, 65]

westPITW = [65, 105]

PITW = [105, 120]

westPITE = [120, 180]

PITE = [180, 230]

eastPITE = [230, 250]

sections = {'westGetz':westGetz, 'Getz':Getz, 'westPITW':westPITW, 'PITW':PITW, 'westPITE':westPITE, 'PITE':PITE, 'eastPITE':eastPITE}



timefile = 'PAS_time.npy'

t = np.load(NPYDIR+timefile)

year = ptools.getDecimalTime(t, t_start=t_start, t_end=t_end)

year = ptools.windowAv(year, n=nn, av=False)

	

# IPO

IPO = getIPO(IPO_t_start=IPO_t_start, IPO_t_end=IPO_t_end)

SAM = getSAM(SAM_t_start=IPO_t_start, SAM_t_end=IPO_t_end)

# Options for ASL are COLS = {'lon':1, 'lat':2, 'actP':3, 'P':4}

ASLlon = getASL(ASL_t_start=IPO_t_start, ASL_t_end=IPO_t_end, COL='lon',)

ASLlat = getASL(ASL_t_start=IPO_t_start, ASL_t_end=IPO_t_end, COL='lat')

ASLactP = getASL(ASL_t_start=IPO_t_start, ASL_t_end=IPO_t_end, COL='actP')

ASLP = getASL(ASL_t_start=IPO_t_start, ASL_t_end=IPO_t_end, COL='P')

ASLrelP = getASL(ASL_t_start=IPO_t_start, ASL_t_end=IPO_t_end, COL='relP')

#plt.plot(IPO,label='1'); plt.plot(IPO2,label='2'); plt.plot(IPO3,label='3'); plt.legend(); plt.show()





#plt.subplot(311); plt.plot(year, ASLlon, label='ASL lon'); plt.ylabel('ASL lon')

#plt.subplot(312); plt.plot(year, ASLlat, label='ASL lat'); plt.ylabel('ASL lat')

#plt.subplot(313); plt.plot(year, ASLactP, label='ASL rel. P'); plt.ylabel('ASl actual P')

#plt.show(); quit()



# Plot which figures?

FIGURE1 = 0

FIGURE2 = 0

FIGURE3 = 0

FIGURE4 = 0

FIGURE5 = 0

FIGURE6 = 1

FIGURE7 = 0

FIGURE8 = 0

FIGURE9 = 0



#==



print('make fig 7 again')



if FIGURE1:



	level = 16



	grid = Grid_PAS(PASDIR)	

	bathy = grid.bathy



	iceC = grid.iceC

	iceC = np.where(iceC==0, np.nan, iceC)



	X = grid.XC; Y = grid.YC

	Xl = X[0,:]; Yl = Y[:,0]

	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2

	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]



	lats = EASlats; lons = EASlons

	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)



	X, Y = grid.XYsubr(lons, lats)

	Z = grid.RC.squeeze()



	# Load time-mean velocities & temp.

	u = np.load(PASDIR+'umean_PAS851_z'+str(level)+'.npy')

	v = np.load(PASDIR+'vmean_PAS851_z'+str(level)+'.npy')

	T = np.load(PASDIR+'Tmean_PAS851_z'+str(level)+'.npy')



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

	land = np.where(bathy<0,0,1)

	

	# This creates land mask with appropriate shade of grey masking.

	nY, nX = bathy.shape

	land = np.zeros((nY, nX, 4))

	for j in range(nY):	

		for i in range(nX):

			if bathy[j,i] == 0:

				land[j,i,:3] = 0.4

			else:

				land[j,i,3] = 1.

				

	#iceC = np.where(bathy<Z[level],iceC,0)

	#iceC = np.ma.array(iceC, mask=bathy==0)

	

	u = ptt.maskDraftXY(u, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	v = ptt.maskDraftXY(v, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskDraftXY(T, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	u = ptt.maskBathyXY(u, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	v = ptt.maskBathyXY(v, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskBathyXY(T, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskBathyXY(T, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	iceC = ptt.maskBathyXY(iceC, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	bathy = ptt.maskBathyXY(bathy, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	bathy = ptt.maskDraftXY(bathy, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	bathyC = ptt.makeBathyContour(bathy, X, Y)



	iceMask = np.ma.getmask(iceC)

	Tmask = np.ma.getmask(bathy)

	#for j in range(Y.shape[0]):

	#	for i in range(Y.shape[1]):

			#if not iceMask[j,i] and iceC[j,i] != 1 and grid.bathy[j,i]>=-455:

	#		if not iceMask[j,i] and Tmask[j,i] and iceC[j,i]!=1:

	#			iceC[j,i]=.25

	#pt.plot1by1(iceC, vmin=0, vmax=1); quit()

	

	# Sample rate

	d = 6



	T[-10,0] = vmin; T[-11,1] = vmax;



	paras = [-75,  -73,  -71]

	merids = [240, 250, 260]

	

	title = r'Time-mean flow & pot. temp. ($^{\circ}$C) at Z = -455 m'



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



	pt.quiver1by1Basemap(u, v, X, Y, d, lat_0, lon_0, contourf=T, mesh=False, contourfNlevels=21, vmin=vmin, vmax=vmax, cmap=cmapMean, contour=[-bathyC,-bathy], contourLevels=[[1000], [500]], contourColours=['silver','white'], parallels=paras, meridians=merids, isf=iceC, land=land, title=title, text_data=text_data, fontsize=10, figsize=(5, 4), cbarTicks=[-2,-1,0,1,2], maskColour='.8', AntarcticInsetData=True, show=True, save=True, outname='Figure_1.png')



#==



if FIGURE2:



	plotSections = ['westGetz', 'westPITW', 'westPITE']

	#plotSections = ['westPITW', 'westPITE']

	

	grid = Grid_PAS(PASDIR)

	areaIS = grid.RAC * (1-grid.hFacW[0]) * (1-grid.hFacS[0])

	

	SHIfwFlxFile = 'SHIfwFlx.npy'

	melt = -np.load(NPYDIR+SHIfwFlxFile)[t_start:t_end] # With minus sign, so +ve into ocean.

	

	SHIfw = -np.load(NPYDIR+'SHIfwFlx.npy')[t_start:t_end]

	shelves = ['PIG', 'THW', 'CRO', 'DOT', 'GETZe', 'COS', 'ABB']#, 'GETZw']

	melts = {}

	IStot = np.zeros(len(year))

	for shelf in shelves:

		melt_ = grid.getIceShelf(SHIfw.copy(), shelf)

		# This one for m/yr

		melt_ = secs_per_year * np.ma.mean(melt_, axis=(1,2)) / rhoFresh

		# This one for Gt/yr

		#melt_ = secs_per_year*np.ma.sum(melt_*areaIS, axis=(1,2))

		melts[shelf] = ptools.windowAv(ptools.detrend(melt_, None), n=nn, nBuffer=nn)

		IStot += melts[shelf]

		

	uv_mean = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections)

	surf_uv_mean = ptools.getSurfCurrent(t_start=t_start, t_end=t_end, sections=plotSections)

	

	wind = np.load(NPYDIR+'slope_uw_av.npy')[t_start:t_end]

	wind = ptools.avSlopeSections(wind, sections=plotSections)

	wind = ptools.detrend(wind, None)

	wind = ptools.windowAv(wind, n=nn)

	

	PIbayLats = [-75.5, -74]; PIbayLons = [251, 262] # PI and Thwaites glaciers

	latsi = grid.getIndexFromLat(PIbayLats); lonsi = grid.getIndexFromLon(PIbayLons)

	melt = tools.getSubregionXY(melt, latsi, lonsi)

		

	melt = ptt.maskDraftXY(melt, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi, inverse=True)

	melt = secs_per_year * np.ma.mean(melt, axis=(1,2)) / rhoFresh

	melt = ptools.detrend(melt, None)

	melt = ptools.windowAv(melt, n=nn)

	

	#===



	print(pearsonr(np.linspace(1,100,100), np.linspace(1,100,100)))

	print('IPO')

	print('Surf current: ' + str(pearsonr(IPO, surf_uv_mean)))

	print('Undercurrent: ' + str(pearsonr(IPO, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(IPO, uv_mean-surf_uv_mean)))

	print('SAM')

	print('Surf current: ' + str(pearsonr(SAM, surf_uv_mean)))

	print('Undercurrent: ' + str(pearsonr(SAM, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(SAM, uv_mean-surf_uv_mean)))

	print('ASLlon')

	print('Surf current: ' + str(pearsonr(ASLlon, surf_uv_mean)))

	print('Undercurrent: ' + str(pearsonr(ASLlon, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(ASLlon, uv_mean-surf_uv_mean)))

	print('ASLlat')

	print('Surf current: ' + str(pearsonr(ASLlat, surf_uv_mean)))	

	print('Undercurrent: ' + str(pearsonr(ASLlat, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(ASLlat, uv_mean-surf_uv_mean)))

	print('ASLactP')

	print('Surf current: ' + str(pearsonr(ASLactP, surf_uv_mean)))	

	print('Undercurrent: ' + str(pearsonr(ASLactP, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(ASLactP, uv_mean-surf_uv_mean)))

	print('ASLP')

	print('Surf current: ' + str(pearsonr(ASLP, surf_uv_mean)))

	print('Undercurrent: ' + str(pearsonr(ASLP, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(ASLP, uv_mean-surf_uv_mean)))

	print('ASLrelP')

	print('Surf current: ' + str(pearsonr(ASLrelP, surf_uv_mean)))

	print('Undercurrent: ' + str(pearsonr(ASLrelP, uv_mean)))

	print('Baroclinicity: ' + str(pearsonr(ASLrelP, uv_mean-surf_uv_mean)))

	

	IPO /= np.max(np.abs(IPO))

	SAM = SAM / np.max(np.abs(SAM))

	ASLlon = ASLlon / np.max(np.abs(ASLlon))

	std = np.mean(IPO**2)**0.5

	

	plt.figure(figsize=(9,4), dpi=200)

	plt.plot(year, 1.e2*uv_mean, color=red, label=r'Undercurrent ($10^{-2}$ m s$^{-1}$)')

	plt.plot(year, 1.e2*surf_uv_mean, color=blue, label=r'Surface current ($10^{-2}$ m s$^{-1}$)')

	#plt.plot(year, 1.e2*(uv_mean-surf_uv_mean), color=grey, label=r'Baroclinicty ($10^{-2}$ m s$^{-1}$)')

	#plt.plot(year, melt, color=green, label='PIG + THW basal melt (m y$^{-1}$)')

	plt.plot(year, IStot/2, color=green, label='Ice-shelf basal melt (2 m yr$^{-1}$)')

	plt.plot(year, wind, color='grey', label='Along-slope wind (m s$^{-1}$)')

	plt.plot(year, IPO, color=black, linestyle='--', label='IPO')

	#plt.plot(year, SAM, color='lightgrey', linestyle='--', label='SAM')

	#plt.plot(year, ASLlon, color='darkgrey', linestyle='--', label='ASL')

	#plt.axhline(std, color='k', linestyle='--'); plt.axhline(-std, color='k', linestyle='--')

	

	plt.ylim(-1.35,1.35)

	plt.title('Along-slope flow, along-slope wind and ice-shelf melt anomalies')

	plt.xlabel('Year')

	plt.ylabel('Current speed / wind speed / basal melt')

	plt.xlim(min(year), max(year))

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)

	plt.legend(prop={'size': 8}, ncol=3)

	plt.savefig('Figure_2.png')

	plt.show()

	

#==



if FIGURE3:



	slopeSections = ['westGetz','westPITW','westPITE']

	slopeSection_indices = ptools.getSectionIndices(slopeSections)



	allSections = ['westGetz','Getz','westPITW','PITW','westPITE']

	#allSections = slopeSections

	allSection_indices = ptools.getSectionIndices(allSections)

	

	thresholdScale = .5

	nn = 60



	ts_shift = -2 + 24

	te_shift = -3 -50

	t_start = 24*12 + 5*12 + ts_shift; t_end = -11 + te_shift



	IPO = ptools.getUndercurrent(t_start=t_start, t_end=t_end); compstr = 'undercurrent'



	# Get grid data

	grid = Grid_PAS(PASDIR)

	yslim = 50; ynlim = 50; nY = ynlim + yslim

	dY = 3.6e0;	Y = np.linspace(-yslim*dY, ynlim*dY, nY)

	zs = [0, -600]; zz = grid.getIndexFromDepth(zs)

	Z = grid.RC.squeeze()[zz[0]:zz[1]]



	# Time.

	time = np.load(NPYDIR+'PAS_time.npy')

	year = ptools.getDecimalTime(time, t_start=t_start, t_end=t_end)

	year = ptools.windowAv(year, n=nn, av=False)



	#==



	# Load and treat density

	S = np.load(NPYDIR+'slopeRho.npy')[t_start:t_end,zz[0]:zz[1],:,allSection_indices].mean(axis=-1)

	varName = 'Rho (kg/m^3)'

	vmin = rho0 - 1.3; vmax = rho0 + 1.5

	vminAnom = -.03; vmaxAnom = .03

	

	Smean = rho0 + np.ma.mean(S, axis=0)

	S = ptools.detrendXY(S, None)#, interceptFlag=0)

	S = ptools.windowAv(S, n=nn)

	compPos, compNeg = ptools.computeComposites(S, IPO, thresholdScale=thresholdScale)



	#==

			

	# Load along-slope velocity.

	uv = np.load(NPYDIR+'slope_uv_xyz.npy')[t_start:t_end,zz[0]:zz[1],:,slopeSection_indices].mean(axis=-1)

	uvmin = -0.04; uvmax = 0.04

	uvminAnom = -0.01; uvmaxAnom = 0.01



	uvmean = np.ma.mean(uv, axis=0)

	uv = tools.subtractSurfFlow(uv)

	uv = ptools.demean(uv)

	#uv = ptools.detrendXY(uv, None)

	uv = ptools.windowAv(uv, n=nn)

	uvCompPos, uvCompNeg = ptools.computeComposites(uv, IPO, thresholdScale=thresholdScale)

		

	#==

	

	# Create bathymetry mask.

	nT, nZ, nY = S.shape

	bathy = np.load(NPYDIR+'slopeBathy.npy')

	bathyMask = np.zeros((S.shape))

	for yi in range(nY):

		bathyMask[:, :, yi] = np.max(bathy[yi, allSection_indices], axis=-1)

	for zi in range(len(Z)):

		bathyMask[:,zi,:] -= Z[zi]



	#==

	

	# Plot

		

	cmap = 'jet'; d = 1.e-4

	xlabel = 'Dist. from slope (km)'; ylabel = 'Depth (m)'

	title = varName + ' ('+str(nn)+'-month running mean)'



	Smean = np.ma.masked_where(bathyMask[0]>0, Smean)

	compPos = np.ma.masked_where(bathyMask[0]>0, compPos)

	compNeg = np.ma.masked_where(bathyMask[0]>0, compNeg)

	Smean = tools.boundData(Smean, vmin, vmax, d=d)

	compPos = tools.boundData(compPos, vminAnom, vmaxAnom, d=d)

	compNeg = tools.boundData(compNeg, vminAnom, vmaxAnom, d=d)

		

	uvmean = np.ma.masked_where(bathyMask[0]>0, uvmean)

	uvCompPos = np.ma.masked_where(bathyMask[0]>0, uvCompPos)

	uvCompNeg = np.ma.masked_where(bathyMask[0]>0, uvCompNeg)	



	#==

	

	# Plots

	titles = [r'(a) Time-mean density (kg m$^{-3}$) & along-slope vel.', r'(b) Pos. ' + compstr + ' density & slope vel. composite', r'(c) Neg. ' + compstr + ' density \& slope vel. composite']

	titles = [r'(a) Time-mean density & along-slope vel.', r'(b) Pos. ' + compstr + ' composite', r'(c) Neg. ' + compstr + ' composite']

	

	#contourLevels = [[-.02,0.,.02]]*2

	cl1 = [-0.04,-.03,-0.02,-0.01,0.,.01,0.02,.03,0.04]

	cl2 = [-.002,-.001,0.,.001,.002]

	contourLevels = [cl1, cl2, cl2]

	

	vmin = [vmin, vminAnom, vminAnom]

	vmax = [vmax, vmaxAnom, vmaxAnom]

	

	yticks = [[-500,-400,-300,-200,-100]]*3

	yticksvis = [True, False, False]

	ylabels = ['Depth (m)', '', '']

	xlabels = ['$y_{slope}$ (km)']*3

	cmaps = [cmapMean, cmapAnom, cmapAnom]

	

	pt.plot1by3([Smean, compPos, compNeg], X=Y, Y=Z, vmin=vmin, vmax=vmax, contourfNlevels=17, titles=titles, contour=[uvmean, uvCompPos, uvCompNeg], contourLevels=contourLevels, DASH_NEG_CONTOURS=True, cmaps=cmaps, xlabels=xlabels, ylabels=ylabels, yticks=yticks, yticksvis=yticksvis, fstitle=7.5, fontsize=10, extend=['both']*3, gridls='dashed', gridlw=0.3, gridlc='k', save=True, outname='Figure_3.png')

	

	quit()

	

#==



if FIGURE4:



	fnames = ['FWflx.npy', 'wk.npy']

	#fnames = ['FWflx.npy', 'SHIfwFlx.npy']

	plotSections = ['westGetz', 'westPITW', 'westPITE']#, 'westPITE']



	lag = None # For testing lags between variables. Default is None.

	

	#==

	

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

	X = grid.XC

	

	# Coastal and area masks

	lonsN=[225, 250]; lonsS=[233.5,258.5]; lonsN=lonsS; Nrange=63 # These should match values for computations elsewhere.

	areaS, areaN = ptools.maskSouthNorthSlope(np.ones(bathy.shape), grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange, FOR_PLOTTING=True)

	coast = ptools.coastalMask(grid, maskRad=15, omitBurkeIsland=True)

	coast = np.ma.masked_where(X<lonsS[0], coast)

	coast = np.ma.masked_where(X>lonsS[1], coast)

	

	iceC = grid.iceC

	iceC = np.where(iceC==0, np.nan, iceC)

		

	latsi = grid.getIndexFromLat(LATS); lonsi = grid.getIndexFromLon(LONS)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	iceC = tools.getSubregionXY(iceC, latsi, lonsi)

	X, Y = grid.XYsubr(LONS, LATS)

	areaS = tools.getSubregionXY(areaS, latsi, lonsi)

	areaN = tools.getSubregionXY(areaN, latsi, lonsi)

	coast = tools.getSubregionXY(coast, latsi, lonsi)



	areaS = ptt.maskBathyXY(areaS, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	areaS = ptt.maskDraftXY(areaS, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	Xl = X[0,:]; Yl = Y[:,0]

	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2

	lat_0 = Yl[Ny2]; lon_0 = Xl[Nx2]

	

	#==

	

	# Get along-slope velocity timeseries

	uv_mean = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections)[lag:]; corrstr = 'undercurrent'

	#uv_mean = ptools.getBarocl(t_start=t_start, t_end=t_end, sections=plotSections)[lag:]; corrstr = 'baroclinicty'

	#uv_mean = -IPO; print('Over-riding regressor to IPO.'); corrstr = '-IPO'

	

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

		

		data = np.load(NPYDIR+fname)[t_start:t_end]#[:-lag]

		data = tools.getSubregionXY(data, latsi, lonsi)

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

		data = ptools.detrendXY(data, None)

		data = ptools.windowAv(data, n=nn)

		

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

	paras = [[-75,  -73,  -71, -69]]*2

	merids = [[230, 240, 250, 260]]*2

	titles = ['(a) Corr(sea-ice FW flux, '+corrstr+')', '(b) Corr(surface stress curl, '+corrstr+')']



	areaNplot = np.where(np.ma.getmask(areaN), 0, areaN); contour2 = [areaNplot, None, None]

	areaSplot = np.where(np.ma.getmask(areaS), 0, areaS); contour3 = [areaSplot, None, None]

	coastPlot = np.where(np.ma.getmask(coast), 0, coast); contour4 = [coastPlot, None, None]

	lw = 1.; cc2 = 'w'; cc3 = 'w'; cc4 = cc3; ls4 = 'dashed'

	

	pt.plot1byNbasemap(corr, X, Y, lat_0, lon_0, contourfNlevels=21, stippling=pval, stipData=[0.01, 6, 6, .2], cmap='coolwarm', vmin=-1., vmax=1., contour=[bathyC,bathyC], contourLevels=[[-1000],[-1000]], contour2=contour2, contourLevels2=[[0]], contour3=contour3, contourLevels3=[[0]], contour4=contour4, contourLevels4=[[0]], lw=0.8, lw2=lw, lw3=lw, lw4=lw, cc2=cc2, cc3=cc3, cc4=cc4, ls4=ls4, cbarTicks=cbarTicks, cbarShrink=.9, landscape=False, isf=[iceC]*2, parallels=paras, meridians=merids, fontsize=9, fstitle=[9]*2, titles=titles, figsize=(5,6), save=True, outname='Figure_4.png')



#==



if FIGURE5:



	si = 4 # For default time params, seasons are [winter, spring, summer, autumn, all]	

	AV = True

	thresholdScale = 0.5

	

	SEASONAL = True; ns = 4;

	SUBREGION = True; lonsi = None; latsi = None



	fname1 = 'FWflx_noPr.npy'; fname2 = 'SHIfwFlx.npy'

	SHIscale = 4

	

	plotSections = ['westGetz', 'westPITW', 'westPITE']

	uv = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections)

	

	predictor = uv; compstr = 'undercurrent '

	#predictor = IPO; compstr = 'IPO '

	

	#==



	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

	X = grid.XC

	

	# Coastal and area masks

	lonsN=[225, 250]; lonsS=[233.5,258.5]; lonsN=lonsS; Nrange=63 # These should match values for computations elsewhere.

	areaS, areaN = ptools.maskSouthNorthSlope(np.ones(bathy.shape), grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange, FOR_PLOTTING=True)

	coast = ptools.coastalMask(grid, maskRad=15, omitBurkeIsland=True)

	coast = np.ma.masked_where(X<lonsS[0], coast)

	coast = np.ma.masked_where(X>lonsS[1], coast)

	

	iceC = grid.iceC

	iceC = np.where(iceC==0, np.nan, iceC)

			

	timefile = 'PAS_time.npy'

	t = np.load(NPYDIR+timefile)

	year = ptools.getDecimalTime(t)

	year = ptools.windowAv(year[t_start:t_end], n=nn, av=False)



	# FW flux data

	data = np.load(NPYDIR+fname1)[t_start:t_end] - np.load(NPYDIR+fname2)[t_start:t_end]

	dataAv = np.mean(data, axis=0)

	data = ptools.detrendXY(data, None)

	data = ptools.windowAv(data, n=nn, nBuffer=nn, av=AV)

	#data = ptools.demean(data)



	# Wind data

	u = np.load(NPYDIR+'EXFuwind.npy')[t_start:t_end]

	v = np.load(NPYDIR+'EXFvwind.npy')[t_start:t_end]

	uAv = np.mean(u, axis=0); vAv = np.mean(v, axis=0)

	u = ptools.detrendXY(u, None); v = ptools.detrendXY(v, None)

	u = ptools.windowAv(u, n=nn, nBuffer=nn, av=AV); v = ptools.windowAv(v, n=nn, nBuffer=nn, av=AV)

	qs = [4.,0.2,0.2]; scale = [100.,4.,4.]; qunits = 'm/s'

	

	# Subregions.

	latsi = grid.getIndexFromLat(LATS); lonsi = grid.getIndexFromLon(LONS)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	iceC = tools.getSubregionXY(iceC, latsi, lonsi)

	data = tools.getSubregionXY(data, latsi, lonsi)

	dataAv = tools.getSubregionXY(dataAv, latsi, lonsi)

	u = tools.getSubregionXY(u, latsi, lonsi)

	uAv = tools.getSubregionXY(uAv, latsi, lonsi)

	v = tools.getSubregionXY(v, latsi, lonsi)

	vAv = tools.getSubregionXY(vAv, latsi, lonsi)

	areaS = tools.getSubregionXY(areaS, latsi, lonsi)

	areaN = tools.getSubregionXY(areaN, latsi, lonsi)

	coast = tools.getSubregionXY(coast, latsi, lonsi)

	X, Y = grid.XYsubr(LONS, LATS)

	

	Xl = X[0,:]; Yl = Y[:,0]

	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2

	lat_0 = Yl[Ny2]; lon_0 = Xl[Nx2]

		

	#dataSeas, seasons = ptools.seasonalDataIPO(data, year, predictor, thresholdScale=thresholdScale)

	dataPosIPO, dataNegIPO = ptools.computeComposites(data, predictor, thresholdScale=thresholdScale)

	uPosIPO, uNegIPO = ptools.computeComposites(u, predictor, thresholdScale=thresholdScale)

	vPosIPO, vNegIPO = ptools.computeComposites(v, predictor, thresholdScale=thresholdScale)

	

	dataPosIPO = np.where(draft<0, dataPosIPO/SHIscale, dataPosIPO)

	dataNegIPO = np.where(draft<0, dataNegIPO/SHIscale, dataNegIPO)

	

	try:

		season = seasons[si].lower() + ' '

	except:

		season = ''

	#print(seasons[si])

	

	# Apply masks

	areaS = ptt.maskBathyXY(areaS, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	areaS = ptt.maskDraftXY(areaS, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	dataAv = ptt.maskBathyXY(dataAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	dataAv = ptt.maskBathyXY(dataAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	#dataAv = ptt.maskDraftXY(dataAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uAv = ptt.maskBathyXY(uAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uAv = ptt.maskDraftXY(uAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vAv = ptt.maskBathyXY(vAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vAv = ptt.maskDraftXY(vAv, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	#dataPosIPO = ptt.maskBathyXY(dataSeas[si,1], grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	dataPosIPO = ptt.maskBathyXY(dataPosIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	#dataPosIPO = ptt.maskDraftXY(dataPosIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	#dataNegIPO = ptt.maskBathyXY(dataSeas[si,2], grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	dataNegIPO = ptt.maskBathyXY(dataNegIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	#dataNegIPO = ptt.maskDraftXY(dataNegIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)



	uPosIPO = ptt.maskBathyXY(uPosIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uPosIPO = ptt.maskDraftXY(uPosIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uNegIPO = ptt.maskBathyXY(uNegIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uNegIPO = ptt.maskDraftXY(uNegIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vPosIPO = ptt.maskBathyXY(vPosIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vPosIPO = ptt.maskDraftXY(vPosIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vNegIPO = ptt.maskBathyXY(vNegIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vNegIPO = ptt.maskDraftXY(vNegIPO, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	#==

	

	# Plotting

	

	bathyC = ptt.makeBathyContour(bathy, X, Y)

	

	XX = [X, X, X]; YY = [Y, Y, Y]

	

	#vmax = 1.e-5; vmin = - vmax

	#vmaxAv = 1.e-4

	#cbarTicks = [[-1.e-5, -.5e-5, 0., .5e-5, 1.e-5]]*2



	vmax = 1; vmin = - vmax; vmaxAv = 2

	cbarTicks = [[vmin, vmin/2, 0., vmax/2, vmax]]*2

	dataAv *= 1.e4; dataPosIPO *= 1.e5; dataNegIPO *= 1.e5



	paras = [[-75,  -73,  -71, -69]]*3

	merids = [[230, 240, 250, 260]]*3

	titles = [r'(a) Pos. ' + compstr + season + 'FW flux composite ($10^{-5}$ kg m$^{-2}$ s$^{-1}$)', r'(b) Neg. ' + compstr + season + 'FW flux composite ($10^{-5}$ kg m$^{-2}$ s$^{-1}$)']

	

	#pt.plot1by2Basemap([dataPosIPO, dataNegIPO], X, Y, lat_0, lon_0, mesh=False, contourfNlevels=21, cbarTicks=cbarTicks, cmap='coolwarm', vmin=vmin, vmax=vmax, contour=[bathyC]*2, contourLevels=[[-1000],[-1000]], cbarShrink=.95, extend=['both']*2, landscape=False, isf=[iceC]*2, parallels=paras, meridians=merids, fontsize=9, fstitle=[8]*2, titles=titles, figsize=(5,6), save=True, outname='Figure_4.png')



	#==

	

	# Figure with extra panel for mean FW flux 

	

	d = 14

	u = [uAv[::d,::d], uPosIPO[::d,::d], uNegIPO[::d,::d]]; v = [vAv[::d,::d], vPosIPO[::d,::d], vNegIPO[::d,::d]]

	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd]*3; Yd = [Yd]*3

	

	isf = [None]*3

	

	contour2 = [None]*3; contour3 = [None]*3; contour4 = [None]*3

	lw = 1.3; cc2 = 'k'; cc3 = None; cc4 = None; ls4 = None

	AREAS = False

	if AREAS:

		areaNplot = np.where(np.ma.getmask(areaN), 0, areaN); contour2 = [areaNplot, None, None]

		areaSplot = np.where(np.ma.getmask(areaS), 0, areaS); contour3 = [areaSplot, None, None]

		coastPlot = np.where(np.ma.getmask(coast), 0, coast); contour4 = [coastPlot, None, None]

		lw = 0.8; cc2 = 'darkgrey'; cc3 = 'w'; cc4 = cc3; ls4 = 'dashed'; contourLevels2=[[0],[0],[0]]

		isf = [iceC]*3

	

	contour2 = [draft]*3; contourLevels2=[[-1],[-1],[-1]]



	fstitle = 7

	titles = ['(a) Time-mean FW flux ($10^{-4}$ kg m$^{-2}$ s$^{-1}$)', r'(b) Pos. ' + compstr + season + 'FW flux composite ($10^{-5}$ kg m$^{-2}$ s$^{-1}$)', r'(c) Neg. ' + compstr + season + 'FW flux composite ($10^{-5}$ kg m$^{-2}$ s$^{-1}$)']

	vmin = [-vmaxAv, vmin, vmin]; vmax = [vmaxAv, vmax, vmax]

	cbarTicks = [[-vmaxAv, -vmaxAv/2, 0., vmaxAv/2., vmaxAv], cbarTicks[0], cbarTicks[1]]

	

	title1 = r'(a) Time-mean winds & FW fluxes ($10^{-4}$ kg m$^{-2}$ s$^{-1}$)'

	title2 = r'(b) Pos. undercurrent composite'+'\n'+r'SI FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1}$) & ice shelf FW flux ($4\times 10^{-5}$ kg m$^{-2}$ s$^{-1}$)'

	title3 = r'(c) Neg. undercurrent composite'+'\n'+r'SI FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1}$) & ice shelf FW flux ($4\times 10^{-5}$ kg m$^{-2}$ s$^{-1}$)'

	titles = [title1, title2, title3]

	

	#pt.plot1byNbasemap([dataAv, dataPosIPO, dataNegIPO], X, Y, lat_0, lon_0, mesh=False, contourfNlevels=21, cbarTicks=cbarTicks, cmap='coolwarm', vmin=vmin, vmax=vmax, contour=[bathyC]*3, contourLevels=[[-1000], [-1000],[-1000]], contour2=contour2, contourLevels2=contourLevels2, contour3=contour3, contourLevels3=[[0]], contour4=contour4, contourLevels4=[[0]], lw=0.8, lw2=lw, lw3=lw, lw4=lw, cc2=cc2, cc3=cc3, cc4=cc4, ls4=ls4, cbarShrink=.95, extend=['both']*3, landscape=False, isf=isf, parallels=paras, meridians=merids, fontsize=9, fstitle=[8]*3, titles=titles, width_ratios=[1,1,1], figsize=(5,8), save=True, outname='Figure_4.png')

	cmap = [cmapMean, cmapAnom, cmapAnom]

	qcolor = ['w', 'k', 'k']

	

	pt.quiver1byN_Basemap(u, v, Xd, Yd, lat_0, lon_0, contourf=[dataAv, dataPosIPO, dataNegIPO], X=XX, Y=YY, mesh=False, contourfNlevels=11, cbarTicks=cbarTicks, cmap=cmap, vmin=vmin, vmax=vmax, contour=[bathyC]*3, contourLevels=[[-1000], [-1000],[-1000]], contour2=contour2, contourLevels2=contourLevels2, contour3=contour3, contourLevels3=[[0]], contour4=contour4, contourLevels4=[[0]], lw=0.8, lw2=lw, lw3=lw, lw4=lw, cc2=cc2, cc3=cc3, cc4=cc4, ls4=ls4, cbarShrink=.95, extend=['both']*3, landscape=False, isf=isf, parallels=paras, meridians=merids, fontsize=9, fstitle=[fstitle]*3, titles=titles, width_ratios=[1,1,1], figsize=(4.5,8), qcolor=qcolor, qs=qs, scale=scale, qunits=qunits, labelpos='W', qlabelx=0.35, save=True, outname='Figure_5.png')

	

#==



if FIGURE6:



	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

	areaIS = grid.RAC * (1-grid.hFacW[0]) * (1-grid.hFacS[0])

	area = grid.RAC * grid.hFacW[0] * grid.hFacS[0]

	coast = ptools.coastalMask(grid, maskRad=15, omitBurkeIsland=True)

	X = grid.XC; Y = grid.YC



	Zbot = 500

	zi = grid.getIndexFromDepth(-Zbot)

	Zbot = int(-grid.RC.squeeze()[zi-1])



	# Replace volume of mean-density water with M gt of freshwater (for pos FW flux anom.)

	# M gt of FW (rhoFresh (kg/m3)) takes up a certain volume V: V = M / rhoFresh

	# Mass of mean density (which is denser) water in same volume is rhoAv * V.

	# Resulting mass difference is rhoAv * V - M = rhoAv * M / rhoFresh - M = M (rhoAv/rhoFresh - 1).



	plotSections = ['westGetz', 'westPITW', 'westPITE']

	uv_mean = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections)

	

	time = np.load(NPYDIR+'PAS_time.npy')

	year = ptools.getDecimalTime(time, t_start=t_start, t_end=t_end)

	year = ptools.windowAv(year, n=nn, av=False)



	fname = 'RHO_zIntAbs_505.npy'

	rho = np.load(NPYDIR+fname)[t_start:t_end]

	

	lonsN=[225, 250]; lonsS=[233.5,258.5]; lonsN=lonsS; Nrange=63

	lonsN=[227.5, 252.5]; lonsS=[233.5,258.5];

	rhoS, rhoN = ptools.maskSouthNorthSlope(rho, grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange)

	areaS, areaN = ptools.maskSouthNorthSlope(area, grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange)

	coast = np.ma.masked_where(X<lonsS[0], coast)

	coast = np.ma.masked_where(X>lonsS[1], coast)

	

	LATS = [-75.5, -68]; LONS = [225, 262]

	#LONS = [235, 259]

	latsi = grid.getIndexFromLat(LATS); lonsi = grid.getIndexFromLon(LONS)

	bathy = ptt.maskBathyXY(tools.getSubregionXY(bathy, latsi, lonsi), grid, zi=0, subregion=True, lats=latsi, lons=lonsi)

	rhoS = ptt.maskBathyXY(tools.getSubregionXY(rhoS, latsi, lonsi), grid, zi=0, timeDep=True, subregion=True, lats=latsi, lons=lonsi)

	areaS = tools.getSubregionXY(areaS, latsi, lonsi)

	areaN = tools.getSubregionXY(areaN, latsi, lonsi)

	coast = tools.getSubregionXY(coast, latsi, lonsi)

	hFacC = tools.getSubregionXY(grid.hFacC, latsi, lonsi)

	X, Y = grid.XYsubr(LONS, LATS)



	#==

	

	# Compute volume, mass and density on shelf.

	volS = np.ma.sum(areaS * grid.DRF[:zi] * hFacC[:zi])

	mass = np.ma.sum(rhoS*areaS, axis=(1,2))

	rhoS = mass / volS



	scaleFacDens = (rhoFresh - rhoS.mean()) / (volS * rhoFresh)

	mass = ptools.movingAv(ptools.detrend(mass, None), n=nn)

	rhoS = ptools.movingAv(ptools.detrend(rhoS, None), n=nn)

	dt_rhoS = 12 * tools.ddt(rhoS, dt=1) 

	

	del rhoN



	# Surface FW fluxes

	FW = np.load(NPYDIR+'FWflx.npy')[t_start:t_end]

	FW = ptools.detrendXY(FW, None)

	FW = ptools.movingAv(FW, n=nn)

	

	# Get south/north FW fluxes in subregion.

	FWS, FWN = ptools.maskSouthNorthSlope(FW, grid, lonsN=lonsN, lonsS=lonsS, Nrange=Nrange)

	del FW

	

	FWS = ptt.maskBathyXY(tools.getSubregionXY(FWS, latsi, lonsi), grid, zi=0, timeDep=True, subregion=True, lats=latsi, lons=lonsi)

	FWcoast = FWS.copy()

	FWS = secs_per_year * np.ma.sum(FWS*areaS, axis=(1,2))

	FWN = ptt.maskBathyXY(tools.getSubregionXY(FWN, latsi, lonsi), grid, zi=0, timeDep=True, subregion=True, lats=latsi, lons=lonsi)

	FWN = secs_per_year * np.ma.sum(FWN*areaN, axis=(1,2))

	

	# Coastal FW flux

	coastT = tools.addTimeDim(coast, FWcoast.shape[0])	

	FWcoast = secs_per_year * np.ma.sum(np.ma.masked_where(coastT<1, FWcoast*areaS), axis=(1,2))

		

	# ICE SHELF MELT

	SHIfw = -np.load(NPYDIR+'SHIfwFlx.npy')[t_start:t_end]

	shelves = ['PIG', 'THW', 'CRO', 'DOT', 'GETZe', 'COS', 'ABB']#, 'GETZw']

	#shelves = ['PIG', 'THW', 'GETZe', 'ABB']

	melts = {}

	IStot = np.zeros(len(year))

	for shelf in shelves:

		melt = grid.getIceShelf(SHIfw.copy(), shelf)

		melt = secs_per_year*np.ma.sum(melt*areaIS, axis=(1,2))

		melts[shelf] = ptools.windowAv(ptools.detrend(melt, None), n=nn, nBuffer=nn)

		IStot += melts[shelf]

	

	#==

	

	# Print some correlations before plotting.

	print('Correlations with undercurrent:')

	print('Ice-shelf melt: ' + str(pearsonr(IStot, uv_mean)))

	print('Off-shelf SI melt: ' + str(pearsonr(FWN, uv_mean)))

	print('On-shelf SI melt: ' + str(pearsonr(FWS, uv_mean)))

	print('On-shelf minus off-shelf: ' + str(pearsonr(FWS-FWN, uv_mean)))

	

	print('Correlations with density & density tendency:')

	print('Density & ice-shelf melt: ' + str(pearsonr(IStot, rhoS)))

	print('Density tend. & ice-shelf melt: ' + str(pearsonr(IStot, dt_rhoS)))

	print('Density & SI melt: ' + str(pearsonr(FWS, rhoS)))

	print('Density tend. & SI melt: ' + str(pearsonr(FWS, dt_rhoS)))

	print('Density & SI coast melt: ' + str(pearsonr(FWcoast, rhoS)))

	print('Density tend. & SI coast melt: ' + str(pearsonr(FWcoast, dt_rhoS)))

	

	#==

	

	# PLOT

		

	fs = 9; lw=1.

	

	plt.figure(figsize=(7,7), dpi=200)

	

	plt.subplot(211)

	#plt.plot(year, 80*IPO/np.max(np.abs(IPO)), label='IPO', color=black, linestyle='dashed')

	plt.plot(year, -FWN/gt, label='-1 x off-shelf sea-ice FW flux', color=blue, linewidth=lw)

	plt.plot(year, FWS/gt, label='On-shelf sea-ice FW flux', color=red, linewidth=lw)

	plt.plot(year, (FWS-FWN)/gt, label='On-shelf minus off-shelf sea-ice FW flux', color='grey', linewidth=lw)

	plt.plot(year, IStot/gt, label='Ice-shelf FW flux', color=green, linewidth=lw)

	plt.title(r'(a) Area integrated sea-ice & ice-shelf FW flux anomalies (Gt yr$^{-1}$)', fontsize=fs)

	plt.ylabel('FW flux anomaly (Gt yr$^{-1}$)', fontsize=fs)

	plt.xlim(year[0], year[-1])

	plt.xticks(ticks=[1990,1995,2000,2005,2010,2015], labels=[], fontsize=fs)

	plt.yticks(fontsize=fs)

	plt.legend(prop={'size': 6})

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)



	plt.subplot(212)

	plt.plot(year, scaleFacDens*FWS, label=r'On-shelf sea-ice FW flux (kg m$^{-3}$ yr$^{-1}$)', color=red, linewidth=lw)

	plt.plot(year, scaleFacDens*FWcoast, label='Coastal sea-ice FW flux (kg m$^{-3}$ yr$^{-1}$)', color=red, linestyle='dashed', linewidth=lw)

	plt.plot(year, scaleFacDens*IStot, label='Ice-shelf FW flux (kg m$^{-3}$ yr$^{-1}$)', color=green, linewidth=lw)

	plt.plot(year, rhoS, label=r'On-shelf density (kg m$^{-3}$)', color=black, linewidth=lw)

	plt.plot(year, dt_rhoS, label=r'On-shelf density tendency (kg m$^{-3}$ yr$^{-1}$)', color=black, linestyle='dashed', linewidth=lw)

	plt.legend(prop={'size': 5.}, ncol=2)

	plt.title('(b) On-shelf density (kg m$^{-3}$) & FW flux density tendencies (kg m$^{-3}$ yr$^{-1})$', fontsize=fs)

	plt.xlabel('Year', fontsize=fs)

	plt.xticks(ticks=[1990,1995,2000,2005,2010,2015], fontsize=fs)

	plt.yticks(fontsize=fs)

	plt.ylabel(r'Density / density tendency', fontsize=fs)

	#plt.ylim(-.035,.035)

	plt.xlim(year[0], year[-1])

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)

	

	plt.savefig('Figure_6.png')

	plt.show()

	

#==



if FIGURE7:



	SUBREGION = False; latsi = None; lonsi = None # So that we may consider subregion if desired in future.

	

	si = 0

	

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	

	iceC = grid.iceC

	iceC = np.where(iceC==0, np.nan, iceC)

	

	X = grid.XC; Y = grid.YC

	

	timefile = 'PAS_time.npy'

	t = np.load(NPYDIR+timefile)

	year = ptools.getDecimalTime(t)

	year = year[t_start:t_end][nn//2:-nn//2+1]



	uwind = np.load(NPYDIR+'EXFuwind.npy')[t_start:t_end][nn//2:-nn//2+1]

	vwind = np.load(NPYDIR+'EXFvwind.npy')[t_start:t_end][nn//2:-nn//2+1]

	atemp = np.load(NPYDIR+'EXFatemp.npy')[t_start:t_end][nn//2:-nn//2+1]

	aqh = np.load(NPYDIR+'EXFaqh.npy')[t_start:t_end][nn//2:-nn//2+1]

	lw = np.load(NPYDIR+'EXFlwdn.npy')[t_start:t_end][nn//2:-nn//2+1]

	

	#atemp = ptools.demean(atemp)

	uwind = ptools.detrendXY(uwind, None)

	vwind = ptools.detrendXY(vwind, None)

	atemp = ptools.detrendXY(atemp, None)

	aqh = ptools.detrendXY(aqh, None)

			

	uwindSeas, seasons = ptools.seasonalDataIPO(uwind, year, IPO)

	vwindSeas, seasons = ptools.seasonalDataIPO(vwind, year, IPO)

	atempSeas, seasons = ptools.seasonalDataIPO(atemp, year, IPO)

	aqhSeas, seasons = ptools.seasonalDataIPO(aqh, year, IPO)

	lwSeas, seasons = ptools.seasonalDataIPO(lw, year, IPO)

			

	print(seasons[si])

	

	d = 28

	

	uPos = ptt.maskBathyXY(uwindSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	uNeg = ptt.maskBathyXY(uwindSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	uPos = ptt.maskDraftXY(uPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

	uNeg = ptt.maskDraftXY(uNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

	

	vPos = ptt.maskBathyXY(vwindSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	vNeg = ptt.maskBathyXY(vwindSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	vPos = ptt.maskDraftXY(vPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

	vNeg = ptt.maskDraftXY(vNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		

	atempPos = ptt.maskBathyXY(atempSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	atempNeg = ptt.maskBathyXY(atempSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	atempPos = ptt.maskDraftXY(atempPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	atempNeg = ptt.maskDraftXY(atempNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	

	aqhPos = ptt.maskBathyXY(aqhSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	aqhNeg = ptt.maskBathyXY(aqhSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	aqhPos = ptt.maskDraftXY(aqhPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	aqhNeg = ptt.maskDraftXY(aqhNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)



	lwPos = ptt.maskBathyXY(lwSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	lwNeg = ptt.maskBathyXY(lwSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	lwPos = ptt.maskDraftXY(lwPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	lwNeg = ptt.maskDraftXY(lwNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	

	#==

	

	Ny2 = len(Y[:,0])//2; Nx2 = len(X[0,:])//2

	lat_0=Y[Ny2,0]; lon_0=X[0,Nx2]



	bathyC = ptt.makeBathyContour(bathy, X, Y)

	contour = [[bathyC, bathyC], [bathyC, bathyC]]

	cl = [-1000]

	contourLevels=[[cl, cl], [cl, cl]]

	

	contourf = [[atempPos, atempNeg], [lwPos, lwNeg]]

	extend = [['both', 'both']*2]*2



	contour2 = [[None, None], [aqhPos, aqhNeg]]

	contourLevels2 = [[], [np.linspace(-2.e-4, 2.e-4, 9)]*2]

	DASH_NEG_CONTOURS = [False, True]

	

	isf = [[iceC]*2]*2

	

	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [[Xd, Xd],[X, X]]; Yd = [[Yd, Yd],[Y, Y]]

	X = [[X, X], [X, X]]; Y = [[Y, Y], [Y, Y]]

	

	vmaxT = 1.5; vminT = -vmaxT

	vmaxLW = 6.0; vminLW = -vmaxLW

	vmax = [[vmaxT, vmaxT], [vmaxLW, vmaxLW]]

	vmin = [[vminT, vminT], [vminLW, vminLW]]

	

	titles = [[r'(a) Pos. IPO winter composite: wind & T ($^{\circ}$C) anoms.', r'(b) Neg. IPO winter composite: wind & T ($^{\circ}$C) anoms.'], [r'(c) Pos. IPO winter composite: LW (W m$^{-2}$) & q anoms.', r'(d) Neg. IPO winter composite: LW (W m$^{-2}$) & q anoms.']]

	

	scale = [[10]*2,[]]

	qs = [[1]*2,[]]

	

	paras = [[[-74, -72, -70, -68, -66, -64]]*2]*2

	merids = [[[230, 240, 250, 260, 270]]*2]*2

	

	u = [[uPos, uNeg], [None, None]]

	v = [[vPos, vNeg], [None, None]]

	

	pt.quiverMbyN_Basemap(u, v, Xd, Yd, lat_0, lon_0, contourf=contourf, X=X, Y=Y, vmin=vmin, vmax=vmax, contourfNlevels=13, extend=extend, isf=isf, contour=contour, contourLevels=contourLevels, lw=0.8, contour2=contour2, contourLevels2=contourLevels2, cc2='w', lw2=1, DASH_NEG_CONTOURS=DASH_NEG_CONTOURS, parallels=paras, meridians=merids, scale=scale, qs=qs, cmap='coolwarm', cbarShrink=0.85, titles=titles, fontsize=6.5, qlabelx=.79, qlabely=0.04, save=True, outname='Figure_7.png', landscape=False, figsize=(7,4.3), w_pad=-0.1, h_pad=-0.1)

		



#==



# Idealised model. Need to show:

# Bathy, wind, rel at north & surface. surf FW flux. 

# Time series of surf flow, deep flow, FW flux.

# Time-mean zonal flow and density, pos/neg composites.

if FIGURE8:



	path = DATADIR + '/MCS_417/run/'

	grid = Grid(path)

	bathy = ptt.maskBathyXY(grid.bathy, grid, 0)



	# Get FW flux from bin data. Take mean over one half period.

	nz = 120; ny = 200; nx = 240; dims = (nz, ny, nx)

	FW = 1.e8 * readnp(path+'precipGaussNsinS05_t.bin', dims, rec=0, dtype='>f8')[0:60].mean(axis=0)

	FW = ptt.maskBathyXY(FW, grid, 0)

	

	# Get T/S relaxation profiles from bin data.

	# Dims of relaxation bin files.

	nz = 50; ny = 200; nx = 240; dims = (nz, ny, nx)

	Trel = readnp(path+'RBCt_100_coldS.bin', dims, rec=0, dtype='>f8')[:,-3,-3]

	Srel = readnp(path+'RBCs_allSurf.bin', dims, rec=0, dtype='>f8')[:,-3,-3]

	

	#==

	

	# Grid and plotting data

	

	Nz = Trel.shape[0]

	X = grid.XC[1,:]/1000.

	Y = grid.YC[:,1]/1000.

	Z = grid.RC.squeeze()

	X2 = 300

	

	xticks = [0, 200, 400, 600]

	yticks = [0, 100, 200, 300, 400, 500]



	vmin = -1000; vmax = -400



	# Wind

	taux = 300+200*tools.get05cosWind(grid.Nx,grid.Ny)[:,0]

	# Latitude grid points to plot arrows at.

	ys = [20, 40, 60, 80, 120, 140, 160, 180] 

	d = 20

	qx = X2

	qy = Y[ys]

	

	qdx = 0.92*(taux[ys] - 300)

	qdy = 0



	#==

	

	# PLOT

	

	import matplotlib.gridspec as gridspec

	width_ratios = (2.4,2.4,1.5)

	fs = 8.5

	fstitle = 8.5



	fig = plt.figure(figsize=(10,3), dpi=200)

	gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig, width_ratios=width_ratios)



	#==

	

	# Subplot 1

	

	ax = fig.add_subplot(gs[0,0])

	ax.patch.set_color('.6')



	levels = np.linspace(vmin, vmax, 13)

	plt.contourf(X, Y, bathy, cmap=cmapMean, levels=levels)#vmin=vmin, vmax=vmax,

	plt.colorbar()

	plt.contour(X, Y, bathy, levels=[-600,-500], colors='k', linestyles='solid', linewidths=0.8)

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)



	for qi in range(len(qy)):

		plt.arrow(qx, qy[qi], qdx[qi], qdy, length_includes_head=True, head_width=10, color='k', zorder=4)

	plt.axvline(X2, color='k', linestyle='--', zorder=4)



	ax.set_xlabel('Lon (km)', fontsize=fs)

	ax.set_ylabel('Lat (km)', fontsize=fs)



	#ax.set_aspect(1)

	plt.plot(taux, Y, color='k')

	ptt.doTicks(xticks, True, yticks, True)

	ax.tick_params(axis='both', which='minor', labelsize=1)

	plt.xticks(fontsize=fs)

	plt.yticks(fontsize=fs)

	

	plt.title('(a) Bathymetry (m) & wind', fontsize=fstitle)



	#==

	

	# Subplot 2

	

	ax = fig.add_subplot(gs[0,1])

	ax.patch.set_color('.6')

	

	levels = np.linspace(-2,2,21)

	plt.contourf(X, Y, FW, cmap=cmapAnom, levels=levels, extend='min')

	plt.colorbar(ticks=[-2,-1,0,1,2])

	plt.contour(X, Y, bathy, levels=[-600,-500], colors='k', linestyles='solid', linewidths=0.8)

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)

	

	#ax.set_aspect(1)

	ax.set_xlabel('Lon (km)', fontsize=fs)

	ax.set_ylabel('Lat (km)', fontsize=fs)

	ptt.doTicks(xticks, True, yticks, True)

	plt.xticks(fontsize=fs)

	plt.yticks(fontsize=fs)

	

	plt.title(r'(b) Imposed FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1})$', fontsize=fstitle)

	plt.tight_layout(w_pad=-3)

	

	#==



	# Subplot 3

	

	ax = fig.add_subplot(gs[0,2])

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)

	ln1n = ax.plot(Trel, Z, color='r', label='Pot. Temp. (deg. C)')

	#ln1s = ax.plot(TrelS, Z, color='r', alpha=0.5, linestyle='--', label='Pot. Temp. (deg. C), South')

	ax.tick_params(axis='x', labelcolor='r', color='r')

	plt.xticks([-2,-1,0,1])

	

	plt.xticks(fontsize=fs)

	plt.yticks(fontsize=fs)

	

	ax.set_xlim(-2.0, 1.2)

	ax.set_ylabel('Depth (m)', fontsize=fs)

	ax2 = ax.twiny()

	

	#ax2.plot(Srel, Z, color='r', label='Pot. Temp.')

	ln2n = ax2.plot(Srel+0.015, Z, color='k', label='Salinity (g/kg)') # Add tiny value to S so both T/S lines visible.

	#ln2s = ax2.plot(SrelS, Z, color='k', alpha=0.5, linestyle='--', label='Salinity (g/kg), South')

	plt.xticks([33., 33.5, 34., 34.5])

	ax2.tick_params(axis='x', labelcolor='k')

	ax2.set_xlim(33, 34.6)

	ax.spines['bottom'].set_color('r')

	ax2.spines['bottom'].set_color('r')

	

	plt.ylim(-1000, 0)

	plt.title(r'(b) $\theta$/$S$ relaxation profiles', fontsize=fstitle)



	lns = ln1n+ln2n

	labs = [l.get_label() for l in lns]

	ax.legend(lns, labs, loc=8, prop={'size': 7})



	plt.xticks(fontsize=fs)

	plt.yticks(fontsize=fs)



	#==

	

	plt.tight_layout()

	plt.savefig('Figure_8.png')

	

	#==

	

#==



if FIGURE9:



	path = DATADIR + 'MCS_417/run/'

	grid = Grid(path)



	X = grid.XC[1,:]/1000.

	Y = grid.YC[:,1]/1000.

	Z = grid.RC.squeeze()



	ts = -12

	# These for mean over last ts months

	u = np.mean(readVariable('UVEL', path, meta=False)[ts:,], axis=(0,-1))

	#T = np.mean(readVariable('THETA', path, meta=False)[ts:], axis=(0,-1))



	

	u = ptt.maskBathyYZ(u, grid, timeDep=False, xi=120)

	T = ptt.maskBathyYZ(T, grid, timeDep=False, xi=120)



	vmin = -0.2; vmax = 0.2

	u[-1,-1] = vmin

	u = tools.boundData(u, vmin, vmax, scale=0.9999)



	#==



	import matplotlib.gridspec as gridspec

	height_ratios = (1,2.)

	fs = 1

	fs2 = 8

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



	plt.title(r'(b) Zonal-mean zonal velocity (m s$^{-1}$) & pot. temp. (deg. C)', fontsize=fs2)

	plt.xlabel('Lat (km)', fontsize=fs)

	plt.ylabel('Depth (m)', fontsize=fs)

	plt.grid()



	plt.savefig('Figure_5.png')

	
