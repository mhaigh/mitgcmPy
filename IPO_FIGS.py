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



from readData import readVariable

from varDict import getPlottingVars, getTrefSref



from TPI import getIPO, getSAM, getASL



import iceShelfFront as isf



import time



#==========================================================



DATADIR = '/home/michael/Documents/data/'

PASDIR = DATADIR + 'PAS_851/run/'

NPYDIR = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'	

	

black = 'k'; blue = 'cornflowerblue'; red = 'indianred'; green = 'seagreen'; grey = 'gray'



rho0 = 1028.5



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

FIGURE3 = 1

FIGURE4 = 0

FIGURE5 = 0

FIGURE6 = 0



#==



print('Add along slope winds to fig1?')

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

	

	#iceC = np.where(bathy<Z[level],iceC,0)

	#iceC = np.ma.array(iceC, mask=bathy==0)

	

	u = ptt.maskDraftXY(u, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	v = ptt.maskDraftXY(v, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskDraftXY(T, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	u = ptt.maskBathyXY(u, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	v = ptt.maskBathyXY(v, grid, level, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	T = ptt.maskBathyXY(T, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

	iceC = ptt.maskBathyXY(iceC, grid, 0, timeDep=False, subregion=True, lats=latsi, lons=lonsi)

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



	pt.quiver1by1Basemap(u, v, X, Y, d, lat_0, lon_0, contourf=T, mesh=True, contourfNlevels=11, vmin=vmin, vmax=vmax, cmap='coolwarm', contour=[-bathyC,-bathy], contourLevels=[[1000], [500]], contourColours=['silver','white'], parallels=paras, meridians=merids, isf=iceC, title=title, text_data=text_data, fontsize=7, figsize=(5, 4), show=True, save=True, outname='Figure_1.png')



#==



if FIGURE2:



	plotSections = ['westGetz', 'westPITW', 'westPITE']

	#plotSections = ['westPITW', 'westPITE']

	

	grid = Grid_PAS(PASDIR)

	

	SHIfwFlxFile = 'SHIfwFlx.npy'

	melt = -np.load(NPYDIR+SHIfwFlxFile)[t_start:t_end] # With minus sign, so +ve into ocean.

	

	uv_mean = ptools.getUndercurrent(t_start=t_start, t_end=t_end, sections=plotSections)

	surf_uv_mean = ptools.getSurfCurrent(t_start=t_start, t_end=t_end, sections=plotSections)

	

	wind = np.load(NPYDIR+'slope_uw_av.npy')[t_start:t_end]

	wind = ptools.avSlopeSections(wind, sections=plotSections)

	wind = ptools.detrend(wind, None)

	wind = ptools.windowAv(wind, n=nn)

	

	PIbayLats = [-75.5, -74]; PIbayLons = [251, 262] # PI and Thwaites glaciers

	latsi = grid.getIndexFromLat(PIbayLats); lonsi = grid.getIndexFromLon(PIbayLons)

	melt = tools.getSubregionXY(melt, latsi, lonsi)

		

	rhoFresh = 1.e3; yearSeconds = 86400. * 365.25

	melt = ptt.maskDraftXY(melt, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi, inverse=True)

	melt = yearSeconds * np.ma.mean(melt, axis=(1,2)) / rhoFresh

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

	

	IPO /= 100*np.max(np.abs(IPO))

	SAM = 0.01 * SAM / np.max(np.abs(SAM))

	ASLlon = 0.01 * ASLlon / np.max(np.abs(ASLlon))

	std = np.mean(IPO**2)**0.5

	

	plt.figure(figsize=(9,4), dpi=200)

	plt.plot(year, uv_mean, color=red, label='Undercurrent (m/s)')

	plt.plot(year, surf_uv_mean, color=blue, label='Surface current (m/s)')

	plt.plot(year, uv_mean-surf_uv_mean, color=grey, label='Baroclinicty (m/s)')

	plt.plot(year, melt/100., color=green, label='PIG + THW basal melt (100 m/y)')

	plt.plot(year, IPO, color=black, linestyle='--', label='IPO')

	plt.plot(year, wind/100, color='darkgrey', linestyle='--', label='Along-slope wind (100 m/s)')

	#plt.plot(year, SAM, color='lightgrey', linestyle='--', label='SAM')

	#plt.plot(year, ASLlon, color='darkgrey', linestyle='--', label='ASL')

	#plt.axhline(std, color='k', linestyle='--'); plt.axhline(-std, color='k', linestyle='--')

	

	plt.title('Along-slope flow and ice shelf melt anomalies')

	plt.xlim(min(year), max(year))

	plt.grid()

	plt.legend(prop={'size': 8})

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

	#uv = tools.subtractSurfFlow(uv)

	uv = ptools.demean(uv)

	uv = ptools.detrendXY(uv, None)

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

	

	pt.plot1by3([Smean, compPos, compNeg], X=Y, Y=Z, vmin=vmin, vmax=vmax, contourfNlevels=17, titles=titles, contour=[uvmean, uvCompPos, uvCompNeg], contourLevels=contourLevels, DASH_NEG_CONTOURS=True, cmaps='coolwarm', xlabels=xlabels, ylabels=ylabels, yticks=yticks, yticksvis=yticksvis, fstitle=7.5, fontsize=10, extend=['both']*3, gridls='dashed', gridlw=0.3, gridlc='k', save=True, outname='Figure_3.png')

	

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

	titles = ['(a) Corr(surface FW flux, '+corrstr+')', '(b) Corr(surface stress curl, '+corrstr+')']



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

	

	d = 12

	u = [uAv[::d,::d], uPosIPO[::d,::d], uNegIPO[::d,::d]]; v = [vAv[::d,::d], vPosIPO[::d,::d], vNegIPO[::d,::d]]

	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd]*3; Yd = [Yd]*3

	

	isf = [None]*3

	

	contour2 = [None]*3; contour3 = [None]*3; contour4 = [None]*3

	lw = 1.; cc2 = 'k'; cc3 = None; cc4 = None; ls4 = None

	AREAS = False

	if AREAS:

		areaNplot = np.where(np.ma.getmask(areaN), 0, areaN); contour2 = [areaNplot, None, None]

		areaSplot = np.where(np.ma.getmask(areaS), 0, areaS); contour3 = [areaSplot, None, None]

		coastPlot = np.where(np.ma.getmask(coast), 0, coast); contour4 = [coastPlot, None, None]

		lw = 0.8; cc2 = 'darkgrey'; cc3 = 'w'; cc4 = cc3; ls4 = 'dashed'; contourLevels2=[[0],[0],[0]]

		isf = [iceC]*3

	

	contour2 = [draft]*3; contourLevels2=[[-1],[-1],[-1]]



	titles = ['(a) Time-mean FW flux ($10^{-4}$ kg m$^{-2}$ s$^{-1}$)', r'(b) Pos. ' + compstr + season + 'FW flux composite ($10^{-5}$ kg m$^{-2}$ s$^{-1}$)', r'(c) Neg. ' + compstr + season + 'FW flux composite ($10^{-5}$ kg m$^{-2}$ s$^{-1}$)']

	vmin = [-vmaxAv, vmin, vmin]; vmax = [vmaxAv, vmax, vmax]

	cbarTicks = [[-vmaxAv, -vmaxAv/2, 0., vmaxAv/2., vmaxAv], cbarTicks[0], cbarTicks[1]]

	

	#pt.plot1byNbasemap([dataAv, dataPosIPO, dataNegIPO], X, Y, lat_0, lon_0, mesh=False, contourfNlevels=21, cbarTicks=cbarTicks, cmap='coolwarm', vmin=vmin, vmax=vmax, contour=[bathyC]*3, contourLevels=[[-1000], [-1000],[-1000]], contour2=contour2, contourLevels2=contourLevels2, contour3=contour3, contourLevels3=[[0]], contour4=contour4, contourLevels4=[[0]], lw=0.8, lw2=lw, lw3=lw, lw4=lw, cc2=cc2, cc3=cc3, cc4=cc4, ls4=ls4, cbarShrink=.95, extend=['both']*3, landscape=False, isf=isf, parallels=paras, meridians=merids, fontsize=9, fstitle=[8]*3, titles=titles, width_ratios=[1,1,1], figsize=(5,8), save=True, outname='Figure_4.png')

	

	pt.quiver1byN_Basemap(u, v, Xd, Yd, lat_0, lon_0, contourf=[dataAv, dataPosIPO, dataNegIPO], X=XX, Y=YY, mesh=False, contourfNlevels=11, cbarTicks=cbarTicks, cmap='coolwarm', vmin=vmin, vmax=vmax, contour=[bathyC]*3, contourLevels=[[-1000], [-1000],[-1000]], contour2=contour2, contourLevels2=contourLevels2, contour3=contour3, contourLevels3=[[0]], contour4=contour4, contourLevels4=[[0]], lw=0.8, lw2=lw, lw3=lw, lw4=lw, cc2=cc2, cc3=cc3, cc4=cc4, ls4=ls4, cbarShrink=.95, extend=['both']*3, landscape=False, isf=isf, parallels=paras, meridians=merids, fontsize=9, fstitle=[8]*3, titles=titles, width_ratios=[1,1,1], figsize=(4.7,8), qs=qs, scale=scale, qunits=qunits, labelpos='W', save=True, outname='Figure_5.png')

	

#==



if FIGURE6:



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



	#atemp = ptools.demean(atemp)

	#uwind = ptools.detrendXY(uwind, None)

	#vwind = ptools.detrendXY(vwind, None)

	#atemp = ptools.detrendXY(atemp, None)

	#aqh = ptools.detrendXY(aqh, None)

			

	uwindSeas, seasons = ptools.seasonalDataIPO(uwind, year, IPO)

	vwindSeas, seasons = ptools.seasonalDataIPO(vwind, year, IPO)

	atempSeas, seasons = ptools.seasonalDataIPO(atemp, year, IPO)

	aqhSeas, seasons = ptools.seasonalDataIPO(aqh, year, IPO)

		

	print(seasons[si])

	del uwind, vwind, atemp, aqh

	

	d = 32

	

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

	

	#==

	

	Ny2 = len(Y[:,0])//2; Nx2 = len(X[0,:])//2

	lat_0=Y[Ny2,0]; lon_0=X[0,Nx2]



	bathyC = ptt.makeBathyContour(bathy, X, Y)

	contour = [bathyC, bathyC]

	contourLevels=[[-1000]]*2



	contour2 = [aqhPos, aqhNeg]

	contourLevels2 = [np.linspace(-0.0002, 0.0002, 5)]*2

		

	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd, Xd]; Yd = [Yd, Yd]

	X = [X, X]; Y = [Y, Y]

	

	vmax = 2.0; vmin = -vmax

	titles = ['(a) Pos IPO wind/T/q composite', '(b) Neg IPO wind/T/q composite']

	

	scale = [10]*2

	qs = [1]*2

	

	paras = [-74, -72, -70, -68, -66, -64]*2

	merids = [230, 240, 250, 260, 270]*2

	

	print('change arrow style so like in figure 1.')

	

	pt.quiver1byN_Basemap([uPos, uNeg], [vPos, vNeg], Xd, Yd, lat_0, lon_0, contourf=[atempPos, atempNeg], X=X, Y=Y, vmin=vmin, vmax=vmax, contourfNlevels=9, extend=['max', 'neither'], isf=[iceC]*2, contour=contour, contourLevels=contourLevels, contour2=[aqhPos, aqhNeg], contourLevels2=contourLevels2, contourLinewidths=0.8, parallels=paras, meridians=merids, scale=scale, qs=qs, cmap='coolwarm', cbarShrink=0.85, title=titles, fontsize=9, save=True, outname='Figure_6.png', landscape=False, figsize=(5,6))

		

	

#==



	
