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



from readData import readVariable, readnp, readVar

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

FIGURE1 = 1

FIGURE2 = 0

FIGURE3 = 0

FIGURE4 = 0



#==

#==



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



	EASlons[0] = 234

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

	

	title = r'(a) Time-mean flow & pot. temp. ($^{\circ}$C) at Z = -455 m'



	# List for ice shelves

	#r = 1.1e6

	r = 1.2e6

	#xloc1 = [2.66e5, 1.55e6+r, 0.9e6+r, 4.6e5+r, 1.85e5+r, 1.52e6+r, 1.4e6+r]

	xloc1 = [2.66e5, 1.55e6+r, 0.93e6+r, 4.6e5+r, 1.9e5+r, 1.52e6+r, 1.43e6+r]

	yloc1 = [4.e5, 1.9e5, 2.e5, 1.95e5, 3.85e5, 8.25e5, 1.2e6]

	text1 = ['GETZ', 'PIG', 'THW', 'CRO', 'DOT', 'COS', 'ABB']

	fontdictIce = {'fontsize':12, 'color':'w'}

	fontdict1 = [fontdictIce]*7

	

	# List for contour labels

	xloc2 = []

	yloc2 = []

	text2 = []

	fontdict2 = []

	

	text0 = {'text':text1+text2, 'xloc':xloc1+xloc2, 'yloc':yloc1+yloc2, 'fontdict':fontdict1+fontdict2}

	text_data = text0

	

	lw = 2.5

	xw = -125; ysw = -72.85; ynw = -71.85

	xe = -108; yse = -71.6; yne = -70.6

	plotvlines = [(xw, ysw, ynw, 'blue', lw), (xe, yse, yne, 'blue', lw)]



	#==

	

	# Now do computations for timeseries



	interceptFlag = 0 # Set this to zero to not remove time-mean when detrending.

	plotSections = ['westGetz', 'westPITW', 'westPITE']

	areaIS = grid.RAC * (1-grid.hFacW[0]) * (1-grid.hFacS[0])

	

	SHIfwFlxFile = 'SHIfwFlx.npy'

	##melt = -np.load(NPYDIR+SHIfwFlxFile)[t_start:t_end] # With minus sign, so +ve into ocean.

	

	SHIfw = -np.load(NPYDIR+'SHIfwFlx.npy')[t_start:t_end]

	shelves = ['PIG', 'THW', 'CRO', 'DOT', 'GETZe', 'COS', 'ABB']#, 'GETZw']

	melts = {}

	IStot = np.zeros(len(year))

	#IStot_ = np.zeros(len(year))

	for shelf in shelves:

		melt_ = grid.getIceShelf(SHIfw.copy()*areaIS, shelf)

		# This one for m/yr

		#melt_ = secs_per_year * np.ma.mean(melt_, axis=(1,2)) / rhoFresh

		# This one for Gt/yr

		melt_ = secs_per_year*np.ma.sum(melt_, axis=(1,2)) / gt

		#melts[shelf] = ptools.windowAv(melt_, n=nn, nBuffer=nn)

		melts[shelf] = ptools.windowAv(ptools.detrend(melt_, None, interceptFlag=interceptFlag), n=nn, nBuffer=nn)

		IStot += melts[shelf]

		#melts[shelf] = ptools.windowAv(ptools.detrend(melt_, None), n=nn, nBuffer=nn)

		#IStot_ += melts[shelf]	

		

	#IStotAv = np.mean(IStot)

	#std = np.mean((IStot - IStotAv)**2)**0.5

	#print(IStotAv, std)

	

	uv = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections, interceptFlag=interceptFlag, fname='slope_uv_av.npy')

	#uv = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections, interceptFlag=interceptFlag)

	surf_uv = ptools.getSurfCurrent(t_start=t_start, t_end=t_end, sections=plotSections, interceptFlag=interceptFlag)

	

	wind = np.load(NPYDIR+'slope_uw_av.npy')[t_start:t_end]

	wind = ptools.avSlopeSections(wind, sections=plotSections)

	wind = ptools.detrend(wind, None, interceptFlag=interceptFlag)

	wind = ptools.windowAv(wind, n=nn)

	

	PIbayLats = [-75.5, -74]; PIbayLons = [251, 262] # PI and Thwaites glaciers

	latsi = grid.getIndexFromLat(PIbayLats); lonsi = grid.getIndexFromLon(PIbayLons)

	##melt = tools.getSubregionXY(melt, latsi, lonsi)

		

	##melt = ptt.maskDraftXY(melt, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi, inverse=True)

	##melt = secs_per_year * np.ma.mean(melt, axis=(1,2)) / rhoFresh

	##melt = ptools.detrend(melt, None)

	##melt = ptools.windowAv(melt, n=nn)

	

	#===

	binwidth = 40

	print('IPO')

	print('Surf current: ' + str(pearsonr(IPO, surf_uv)))

	#print(ptools.corress(IPO, surf_uv, binwidth, 0))

	print('Undercurrent: ' + str(pearsonr(IPO, uv)))

	print('Baroclinicity: ' + str(pearsonr(IPO, uv-surf_uv)))

	print('SAM')

	print('Surf current: ' + str(pearsonr(SAM, surf_uv)))

	print('Undercurrent: ' + str(pearsonr(SAM, uv)))

	print('Baroclinicity: ' + str(pearsonr(SAM, uv-surf_uv)))

	print('ASLlon')

	print('Surf current: ' + str(pearsonr(ASLlon, surf_uv)))

	print('Undercurrent: ' + str(pearsonr(ASLlon, uv)))

	print('Baroclinicity: ' + str(pearsonr(ASLlon, uv-surf_uv)))

	print('ASLlat')

	print('Surf current: ' + str(pearsonr(ASLlat, surf_uv)))	

	print('Undercurrent: ' + str(pearsonr(ASLlat, uv)))

	print('Baroclinicity: ' + str(pearsonr(ASLlat, uv-surf_uv)))

	print('ASLactP')

	print('Surf current: ' + str(pearsonr(ASLactP, surf_uv)))	

	print('Undercurrent: ' + str(pearsonr(ASLactP, uv)))

	print('Baroclinicity: ' + str(pearsonr(ASLactP, uv-surf_uv)))

	print('ASLP')

	print('Surf current: ' + str(pearsonr(ASLP, surf_uv)))

	print('Undercurrent: ' + str(pearsonr(ASLP, uv)))

	print('Baroclinicity: ' + str(pearsonr(ASLP, uv-surf_uv)))

	print('ASLrelP')

	print('Surf current: ' + str(pearsonr(ASLrelP, surf_uv)))

	print('Undercurrent: ' + str(pearsonr(ASLrelP, uv)))

	print('Baroclinicity: ' + str(pearsonr(ASLrelP, uv-surf_uv)))

	

	IPO = getIPO(IPO_t_start=IPO_t_start, IPO_t_end=IPO_t_end, interceptFlag=interceptFlag)

	IPO /= np.max(np.abs(IPO))

	SAM = SAM / np.max(np.abs(SAM))

	ASLlon = ASLlon / np.max(np.abs(ASLlon))

	std = np.mean(IPO**2)**0.5

	

	#==

	

	# This option for non-normalised data.

	#ts_data = [1.e2*uv, 1.e2*surf_uv, IStot/50, wind, IPO]

	#ts_ylim = (-1.35, 1.35)

	#ts_labels = [r'Undercurrent ($10^{-2}$ m s$^{-1}$)', r'Surface current ($10^{-2}$ m s$^{-1}$)', 'Ice-shelf basal melt (50 Gt yr$^{-1}$)', 'Along-slope wind (m s$^{-1}$)', 'IPO']

	#ts_legendLoc = 'lower right'

	#ts_legendFontsize = 8. 

	#ts_legendNcol = 3

	

	#==

	

	# This option for normalised data.

	# Standardise, print standard deviations and means.

	uv, AV, STD = ptools.normaliseSTD(1.e2*uv, printVarName='Undercurrent', outputStats=True)

	uvt = r'Undercurrent, mean=%.2f cm s$^{-1}$, std. dev.=%.2f cm s$^{-1}$' % (AV, STD)

	uvt = r'Undercurrent, (mean,s.d.)=(%.2f,%.2f) cm s$^{-1}$' % (AV, STD)

	

	surf_uv, AV, STD = ptools.normaliseSTD(1.e2*surf_uv, printVarName='Surface current', outputStats=True)

	surf_uvt = r'Surface current mean=%.2f cm s$^{-1}$, std. dev.=%.2f cm s$^{-1}$' % (AV, STD)

	surf_uvt = r'Surf. current (mean,s.d.)=(%.2f,%.2f) cm s$^{-1}$' % (AV, STD)

	

	IStot, AV, STD = ptools.normaliseSTD(IStot, printVarName='Ice-shelf melt', outputStats=True)

	IStott = r'Ice-shelf basal melt, mean=%.2f Gt yr$^{-1}$, std. dev.=%.2f Gt yr$^{-1}$' % (AV, STD)

	IStott = r'Ice-shelf melt (mean,s.d.)=(%.2f,%.2f) Gt yr$^{-1}$' % (AV, STD)

	

	wind, AV, STD = ptools.normaliseSTD(wind, printVarName='Winds', outputStats=True)

	windt = r'Along-slope wind, mean=%.2f m s$^{-1}$, std. dev.=%.2f m s$^{-1}$' % (AV, STD)

	windt = r'Wind (mean,s.d.)=(%.2f,%.2f) m s$^{-1}$' % (AV, STD)

	

	IPO, AV, STD = ptools.normaliseSTD(IPO, printVarName='IPO', outputStats=True)

	IPOt = r'IPO, mean=%.2f$^{\circ}$C, std. dev.=%.2f$^{\circ}$C' % (AV, STD)

	IPOt = r'IPO (mean,s.d.)=(%.2f,%.2f)$^{\circ}$C' % (AV, STD)

	

	ts_data = [IPO, wind, IStot, uv, surf_uv]

	ts_labels = [IPOt, windt, IStott, uvt, surf_uvt]

	ts_ylim = (-2.,3.)

	ts_yticks = np.linspace(-2., 3., 11)

	ts_legendLoc = 'upper right'

	ts_legendFontsize = 7.4

	ts_legendNcol = 2

	

	#==

	

	ts_time = [year]*len(ts_data)

	

	solid = 'solid'; dashed = 'dashed'

	ts_lines = [(black, solid), (grey, solid), (green, solid), (red, solid), (blue, solid)]

	ts_title = '(b) Eastward along-slope flow, wind and ice-shelf melt anomalies (units std. dev.)'

	ts_xlabel = 'Year'

	ts_ylabel = 'Current/wind speed & basal melt'

	ts_xlim = (min(year), max(year))

	

	#==

	

	# PLOT

	

	contourColours = ['silver', 'white']

	contourColours = ['cyan', 'white']

	

	dxArrow = 0.8e5; dxArrowDiag = 5.67e4

	GetzArrow = [.65e6, 1.50e6, dxArrowDiag, -dxArrowDiag, 'r']

	arrows = [[1.58e6, 8.2e5, dxArrow, 0, 'r'], [.93e6, 1.28e6, -dxArrowDiag, dxArrowDiag, 'b'], [1.18e6, 1.55e6, dxArrow, 0, 'b']]

	#arrows = [[-112., -74., 1, 0, 'r']



	pt.quiver1by1Basemap_timeseries(u, v, X, Y, d, lat_0, lon_0, contourf=T, mesh=False, contourfNlevels=21, vmin=vmin, vmax=vmax, cmap=cmapMean, contour=[-bathyC,-bathy], lw=lw, plotvlines=plotvlines, contourLevels=[[1000], [500]], contourColours=contourColours, parallels=paras, meridians=merids, isf=iceC, land=land, title=title, text_data=text_data, fontsize=12, fstitle=14, cbarTicks=[-2,-1,0,1,2], maskColour='.8', AntarcticInsetData=True, arrows=arrows, \

	ts_time=ts_time, ts_data=ts_data, ts_labels=ts_labels, ts_lines=ts_lines, ts_title=ts_title, ts_xlabel=ts_xlabel, ts_ylabel=ts_ylabel, ts_xlim=ts_xlim, ts_ylim=ts_ylim, ts_yticks=ts_yticks, ts_legendLoc=ts_legendLoc, ts_legendFontsize=ts_legendFontsize, ts_legendNcol=ts_legendNcol, \

	show=False, save=True, outname='Figure_1.png', height_ratios=[1.8,1], figsize=(9,9.2))

	

#==



if FIGURE2:



	slopeSections = ['westGetz','westPITW','westPITE']

	slopeSection_indices = ptools.getSectionIndices(slopeSections)



	allSections = ['westGetz','Getz','westPITW','PITW','westPITE']

	allSection_indices = ptools.getSectionIndices(allSections)

	

	# Get grid data

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

	

	latsi = grid.getIndexFromLat(LATS); lonsi = grid.getIndexFromLon(LONS)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	X, Y = grid.XYsubr(LONS, LATS)

	

	Xl = X[0,:]; Yl = Y[:,0]

	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2

	lat_0 = Yl[Ny2]; lon_0 = Xl[Nx2]

	

	yslim = 50; ynlim = 50; nY = ynlim + yslim

	dY = 3.6e0;	Yslope = np.linspace(-yslim*dY, ynlim*dY, nY)

	zs = [0, -600]; zz = grid.getIndexFromDepth(zs)

	Z = grid.RC.squeeze()[zz[0]:zz[1]]

	

	#==



	# Load along-slong vel. and density data & mask.

	uvmean, uvCompPos, uvCompNeg = np.load(PASDIR+'uvComp.npy')

	uvmin = -0.04; uvmax = 0.04

	uvminAnom = -0.01; uvmaxAnom = 0.01

	

	RhoMean, RhoCompPos, RhoCompNeg = np.load(PASDIR+'rhoComp.npy')

	vmin = rho0 - 1.3; vmax = rho0 + 1.5

	vminAnom = -.03; vmaxAnom = .03

		

	# Create bathymetry mask.

	nZ, nY = RhoMean.shape

	slopeBathy = np.load(NPYDIR+'slopeBathy.npy')

	bathyMask = np.zeros(RhoMean.shape)

	for yi in range(nY):

		bathyMask[:, yi] = np.max(slopeBathy[yi, allSection_indices], axis=-1)

	for zi in range(len(Z)):

		bathyMask[zi,:] -= Z[zi]

		

	RhoMean = np.ma.masked_where(bathyMask>0, RhoMean)

	RhoCompPos = np.ma.masked_where(bathyMask>0, RhoCompPos)

	RhoCompNeg = np.ma.masked_where(bathyMask>0, RhoCompNeg)

	

	d = 1.e-4

	RhoMean = tools.boundData(RhoMean, vmin, vmax, d=d)

	RhoCompPos = tools.boundData(RhoCompPos, vminAnom, vmaxAnom, d=d)

	RhoCompNeg = tools.boundData(RhoCompNeg, vminAnom, vmaxAnom, d=d)

		

	uvmean = np.ma.masked_where(bathyMask>0, uvmean)

	uvCompPos = np.ma.masked_where(bathyMask>0, uvCompPos)

	uvCompNeg = np.ma.masked_where(bathyMask>0, uvCompNeg)

		

	#==

	

	# Load FW flux and wind data

	uwindMean, uwindPos, uwindNeg = np.load(PASDIR+'uwindComp.npy')

	vwindMean, vwindPos, vwindNeg = np.load(PASDIR+'vwindComp.npy')

	fwMean, fwPos, fwNeg = np.load(PASDIR+'fwComp.npy')



	uwindMean = ptt.maskBathyXY(uwindMean, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uwindMean = ptt.maskDraftXY(uwindMean, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uwindPos = ptt.maskBathyXY(uwindPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uwindPos = ptt.maskDraftXY(uwindPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uwindNeg = ptt.maskBathyXY(uwindNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	uwindNeg = ptt.maskDraftXY(uwindNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	vwindMean = ptt.maskBathyXY(vwindMean, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vwindMean = ptt.maskDraftXY(vwindMean, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vwindPos = ptt.maskBathyXY(vwindPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vwindPos = ptt.maskDraftXY(vwindPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vwindNeg = ptt.maskBathyXY(vwindNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vwindNeg = ptt.maskDraftXY(vwindNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	fwMean = ptt.maskBathyXY(fwMean, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	fwPos = ptt.maskBathyXY(fwPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	fwNeg = ptt.maskBathyXY(fwNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	#==

	

	# Plotting details for sections

	sec_titles = [r'(a) Time-mean density (kg m$^{-3}$)'+'\n'+'& eastward along-slope vel.', r'(c) Density (kg m$^{-3}$) &'+'\n'+'eastward along-slope vel. (+ comp)', r'(e) Density (kg m$^{-3}$) &'+'\n'+'eastward along-slope vel. (- comp)']

	

	#contourLevels = [[-.02,0.,.02]]*2

	cl1 = [-0.04,-.03,-0.02,-0.01,0.,.01,0.02,.03,0.04]

	cl2 = [-.002,-.001,0.,.001,.002]

	sec_contourLevels = [cl1, cl2, cl2]

	

	sec_vmin = [vmin, vminAnom, vminAnom]

	sec_vmax = [vmax, vmaxAnom, vmaxAnom]

	

	sec_yticks = [[-500,-400,-300,-200,-100]]*3

	sec_yticksvis = [True]*3

	sec_ylabels = ['Depth (m)']*3

	sec_xlabels = ['','','$y_{slope}$ (km)']

	sec_cmaps = [cmapMean, cmapAnom, cmapAnom]

	

	#pt.plot1by3([RhoMean, RhoCompPos, RhoCompNeg], X=Yslope, Y=Z, vmin=sec_vmin, vmax=sec_vmax, contourfNlevels=17, titles=sec_titles, contour=[uvmean, uvCompPos, uvCompNeg], contourLevels=sec_contourLevels, DASH_NEG_CONTOURS=True, cmaps=sec_cmaps, xlabels=sec_xlabels, ylabels=sec_ylabels, yticks=sec_yticks, yticksvis=sec_yticksvis, fstitle=7.5, fontsize=10, extend=['both']*3, gridls='dashed', gridlw=0.3, gridlc='k', save=True, outname='Figure_3.png')

	

	#==

	

	# Plotting details for Basemap plot

	

	XX = [X]*3; YY = [Y]*3



	bathyC = ptt.makeBathyContour(bathy, X, Y)

	paras = [[-75,  -73,  -71, -69]]*3

	merids = [[230, 240, 250, 260]]*3

	

	d = 14

	u = [uwindMean[::d,::d], uwindPos[::d,::d], uwindNeg[::d,::d]]; v = [vwindMean[::d,::d], vwindPos[::d,::d], vwindNeg[::d,::d]]

	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd]*3; Yd = [Yd]*3

	

	isf = [None]*3

	

	contour2 = [None]*3; contour3 = [None]*3; contour4 = [None]*3

	lw = 1.3; cc2 = 'k'; cc3 = None; cc4 = None; ls4 = None

	contour2 = [draft]*3; contourLevels2=[[-1],[-1],[-1]]

	

	vmax = 1; vmin = - vmax; vmaxAv = 2

	cbarTicks = [[vmin, vmin/2, 0., vmax/2, vmax]]*2

	fwMean *= 1.e4; fwPos *= 1.e5; fwNeg *= 1.e5

	

	fstitle = 8

	vmin = [-vmaxAv, vmin, vmin]; vmax = [vmaxAv, vmax, vmax]

	cbarTicks = [[-vmaxAv, -vmaxAv/2, 0., vmaxAv/2., vmaxAv], cbarTicks[0], cbarTicks[1]]

	cbarLabels = ['TIME MEANS', 'POS. UNDERCURRENT COMPS.', 'NEG. UNDERCURRENT COMPS.']

	

	title1 = r'(b) Time-mean winds & FW fluxes ($10^{-4}$ kg m$^{-2}$ s$^{-1}$)'

	title2 = r'(d) Pos. undercurrent composite'+'\n'+r'SI FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1}$) & ice-shelf FW flux ($4\times 10^{-5}$ kg m$^{-2}$ s$^{-1}$)'

	title3 = r'(f) Neg. undercurrent composite'+'\n'+r'SI FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1}$) & ice-shelf FW flux ($4\times 10^{-5}$ kg m$^{-2}$ s$^{-1}$)'

	

	title1 = r'(b) Time-mean winds & FW fluxes ($10^{-4}$ kg m$^{-2}$ s$^{-1}$)'

	title2 = r'(d) Winds, sea-ice FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1}$) &'+'\n'+r'ice-shelf FW flux ($4\times 10^{-5}$ kg m$^{-2}$ s$^{-1}$) (+ comp)'

	title3 = r'(f) Winds, sea-ice FW flux ($10^{-5}$ kg m$^{-2}$ s$^{-1}$) &'+'\n'+r'ice-shelf FW flux ($4\times 10^{-5}$ kg m$^{-2}$ s$^{-1}$) (- comp)'

	

	titles = [title1, title2, title3]

	

	cmap = [cmapMean, cmapAnom, cmapAnom]

	qcolor = ['w', 'k', 'k']

	qs = [4.,0.2,0.2]; scale = [100.,4.,4.]; qunits = 'm/s'

	

	#width_ratios=[1.1,2], save=True, figsize=(7.5,8)

	

	# PLOT

	pt.quiver1byNBasemapAndSections(u, v, Xd, Yd, lat_0, lon_0, contourf=[fwMean, fwPos, fwNeg], X=XX, Y=YY, mesh=False, contourfNlevels=11, cbarTicks=cbarTicks, cbarLabels=cbarLabels, cmap=cmap, vmin=vmin, vmax=vmax, contour=[bathyC]*3, contourLevels=[[-1000], [-1000],[-1000]], contour2=contour2, contourLevels2=contourLevels2, contour3=contour3, contourLevels3=[[0]], contour4=contour4, contourLevels4=[[0]], lw=0.8, lw2=lw, lw3=lw, lw4=lw, cc2=cc2, cc3=cc3, cc4=cc4, ls4=ls4, cbarShrink=.95, extend=['both']*3, isf=isf, parallels=paras, meridians=merids, fontsize=9, fstitle=[fstitle]*3, titles=titles, qcolor=qcolor, qs=qs, scale=scale, qunits=qunits, labelpos='W', qlabelx=0.36, \

	sec_data=[RhoMean, RhoCompPos, RhoCompNeg], sec_X=Yslope, sec_Y=Z, sec_vmin=sec_vmin, sec_vmax=sec_vmax, sec_contourfNlevels=[15,13,13], sec_titles=sec_titles, sec_contour=[uvmean, uvCompPos, uvCompNeg], sec_contourLevels=sec_contourLevels, sec_DASH_NEG_CONTOURS=True, sec_cmaps=sec_cmaps, sec_xlabels=sec_xlabels, sec_ylabels=sec_ylabels, sec_yticks=sec_yticks, sec_yticksvis=sec_yticksvis, sec_fstitle=fstitle, sec_fontsize=9, sec_extend=['both']*3, sec_gridls='dashed', sec_gridlw=0.3, sec_gridlc='k', \

show=False, width_ratios=[1.2,2], save=True, figsize=(7.7,7.9), outname='Figure_2.png')



#==



if FIGURE3:

	

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

	X = grid.XC

	

	# Coastal and area masks

	#lonsN=[225, 250]; lonsS=[233.5,258.5]; lonsN=lonsS; Nrange=63 # These should match values for computations elsewhere.

	lonsN=[225.1, 258.5]; lonsS=[225.1,258.5]; Nrange=58; # lonsN=lonsS;

	

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

	

	corr = np.load(PASDIR+'corr_FWflx.npy')

	pval = np.load(PASDIR+'pval_FWflx.npy')

	

	corr = ptt.maskBathyXY(corr, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	pval = np.where(bathy>=0, 1, pval)

	#pval = np.where(draft<0, 1, pval)

	

	#=

	#=

	

	# Plotting data for correlation panel

	

	title = '(a) Correlation: freshwater fluxes & undercurrent'



	cbarTicks = [-1., -.5, 0., .5, 1.]

	paras = [-75,  -73,  -71, -69]

	merids = [230, 240, 250, 260]



	areaNplot = np.where(np.ma.getmask(areaN), 0, areaN)

	areaSplot = np.where(np.ma.getmask(areaS), 0, 1)

	nY,nX = areaSplot.shape

	for j in range(nY):	

		for i in range(nX):

			if areaSplot[j,i] == 1:

				#if areaSplot[j,i+1] == 0 or areaSplot[j,i-1] == 0 or areaSplot[j+1,i] == 0 or areaSplot[j-1,i] == 0:

				if areaSplot[j,i+1] != 0 and areaSplot[j,i-1] != 0 and areaSplot[j+1,i] != 0 and areaSplot[j-1,i] != 0:

					areaSplot[j,i] = 2

	#for j in range(nY):	

	#	for i in range(nX):

	#		if areaSplot[j,i] == 2:

	#			#if areaSplot[j,i+1] == 0 or areaSplot[j,i-1] == 0 or areaSplot[j+1,i] == 0 or areaSplot[j-1,i] == 0:

	#			if areaSplot[j,i+1] != 1 and areaSplot[j,i-1] != 1 and areaSplot[j+1,i] != 1 and areaSplot[j-1,i] != 1:

	#				print('here')

	#				areaSplot[j,i] = 3

	coastPlot = np.where(np.ma.getmask(coast), 0, coast); contour4 = [coastPlot, None, None]

	

	bathyC = ptt.makeBathyContour(bathy, X, Y)

	contourLW = [2.,2.8,2.8,2.8]

	solid = 'solid'; dashed = 'dashed'



	# For on-shelf area outline over ice-shelf draft outline

	contour = [bathyC, draft, areaNplot, areaSplot]

	contourLS = [solid, solid, dashed, dashed]

	contourLevels = [[-1000], [-1], [0], [1]]

	contourColours = ['k', 'k', 'w', 'w']

	# For ice-shelf draft outline over on-shelf area outline.

	contour = [bathyC, areaNplot, areaSplot, draft]

	contourLS = [solid, dashed, dashed, solid]

	contourLevels = [[-1000], [0], [1], [-1]]

	contourColours = ['k', 'w', 'w', 'k']

	

	vlines = [(-125., -74.7, -74.0, 'k', 2.)]

	

	#=

	#=

	

	# Plotting data for timeseries panel

	

	FWN = np.load(PASDIR+'FWN.npy')

	FWS = np.load(PASDIR+'FWS.npy')

	IStot = np.load(PASDIR+'IStot.npy')

	

	gt = 1.e12

	ts_data = [-FWN/gt, FWS/gt, (FWS-FWN)/gt, IStot/gt]

	ts_time = [year]*len(ts_data)

	ts_labels = ['-1 x off-shelf sea-ice FW flux', 'On-shelf sea-ice FW flux', 'On-shelf minus off-shelf sea-ice FW flux', 'Ice-shelf FW flux']	

	

	ts_lines = [(blue, solid), (red, solid), ('grey', solid), (green, solid)]

	ts_title = '(b) Area integrated sea-ice & ice-shelf FW flux anomalies (Gt yr$^{-1}$)'

	ts_xlabel = 'Year'

	ts_ylabel = 'FW flux anomaly (Gt yr$^{-1}$)'

	ts_xlim = (min(year), max(year))

	ts_ylim = (-75, 125)

	ts_fslegend = 9

	#ts_xticks = [1990,1995,2000,2005,2010,2015]

	

	#=

	#=

	

	# Plot

	

	pt.plot1by1Basemap_timeseries(corr, X, Y, lat_0, lon_0, contourfNlevels=21, stippling=pval, stipData=[0.01, 6, 6, 2.], cmap='coolwarm', vmin=-1., vmax=1., cbarTicks=cbarTicks, parallels=paras, meridians=merids, fontsize=12, fstitle=14, title=title, contour=contour, contourLevels=contourLevels, contourColours=contourColours, contourLW=contourLW, contourLS=contourLS, mesh=False, plotvlines=vlines, \

	ts_time=ts_time, ts_data=ts_data, ts_labels=ts_labels, ts_lines=ts_lines, ts_title=ts_title, ts_xlabel=ts_xlabel, ts_ylabel=ts_ylabel, ts_xlim=ts_xlim, ts_ylim=ts_ylim, ts_fslegend=ts_fslegend, \

	show=False, save=True, outname='Figure_3.png', figsize=(9,9.2), height_ratios=[1.8,1])

	

	

	

	

	

	

	

	

	
