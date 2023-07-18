import numpy as np



from time import ctime



from grid_PAS import Grid as Grid_PAS



import matplotlib.pyplot as plt

import matplotlib.colors as cl

import matplotlib.animation as animation

import matplotlib.ticker as ticker



import plotting as pt

import plotting_tools as ptt



import PAS_tools

import tools



from readData import * 



#==========================================================



DATADIR = '/home/michael/Documents/data/'

PASDIR = DATADIR + 'PAS_851/run/'

EASlats = [-75.5, -70.5]; EASlons = [235, 260] 

EASlats = [-75.5, -68]



westGetz = [1, 60]

Getz = [60, 65]

westPITW = [65, 105]

PITW = [105, 120]

westPITE = [120, 180]

PITE = [180, 230]

eastPITE = [230, 250]

	

sections = {'westGetz':westGetz, 'Getz':Getz, 'westPITW':westPITW, 'PITW':PITW, 'westPITE':westPITE, 'PITE':PITE, 'eastPITE':eastPITE}

	



MAIN = False

if MAIN:



	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy



	X = grid.XC; Y = grid.YC

	Z = grid.RC.squeeze()



	pathno = 0; xshift = 0; yshift = 0; latsi = []; lonsi = []

	subr = 1

	if subr:

		# Get subregion of bathymetry.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		pathno = 0

		xshift = lonsi[0]

		yshift = latsi[0]

	nY, nX = X.shape



	# Get lat, lon lists of continental slope

	slope_x, slope_y = PAS_tools.getSlopeContour(bathy, X, Y, pathno)

	np.save('slope_x', slope_x); np.save('slope_y', slope_y)



	bathy = ptt.maskBathyXY(bathy, grid, zi=0, subregion=True, lats=latsi, lons=lonsi)



	bearing = PAS_tools.getBearing(slope_x, slope_y)



	slope_x_nX = np.zeros(nX); slope_y_nX = np.zeros(nX)

	slope_xi = np.zeros(nX, dtype=int); slope_yi = np.zeros(nX, dtype=int)

	bearing_nX = np.zeros(nX)



    # For each X grid point, get the lat and lon of slope contour's nearest point.

	for xi in range(nX):

		slope_xi[xi] = int(np.argmin(np.abs(slope_x-X[0,xi])))

		slope_x_nX[xi] = slope_x[slope_xi[xi]]

		slope_y_nX[xi] = slope_y[slope_xi[xi]]

		bearing_nX[xi] = bearing[slope_xi[xi]]

		slope_yi[xi] = int(np.argmin(np.abs(slope_y_nX[xi]-Y[:,xi])))



    # Plot bathy, check contour points and bearing.

    #plt.subplot(121); plt.contourf(X, Y, bathy); plt.scatter(slope_x_nX, slope_y_nX)

    #plt.subplot(122); plt.plot(X[0], bearing_nX); plt.show(); quit()



    # For each timestep, load instances of u, v, rho.

    # Depth average beneath a rho contour. (27.95?)

    # Then take dot prodcut with unit vector in direction of slope's bearing.

    # Surface flow can be just the flow or estimated from SSH for geostrophic part.



	zs = [0, -1000]

	zz = grid.getIndexFromDepth(zs)



	ufile = 'UVEL'; vfile = 'VVEL'; rhofile = 'RHOAnoma'

	uwindfile = 'EXFuwind'; vwindfile = 'EXFvwind'

	ustressfile = 'oceTAUX'; vstressfile = 'oceTAUY'



	#t_start = 0; t_end = 10

	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779

	surf_uv_max = np.zeros((t_end-t_start, nX))

	surf_uv_av = np.zeros((t_end-t_start, nX))

	slope_uv_av = np.zeros((t_end-t_start, nX))

	slope_uv_max = np.zeros((t_end-t_start, nX))

	slope_uw_av = np.zeros((t_end-t_start, nX))

	slope_us_av = np.zeros((t_end-t_start, nX))

	for ti in range(t_start, t_end):

		print(ti)

		print(ctime(readVariable(rhofile, PASDIR, file_format='nc', meta=True)['TIME'][ti]))

		for xi in range(nX):



			yy = [yshift+slope_yi[xi]-3, yshift+slope_yi[xi]+1]

			rho = 1028.5 + readVariable(rhofile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			rho = ptt.maskBathyYZ(rho, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			rho = tools.boundData(rho, 1026, 1040)



			# Load velocities

			u = readVariable(ufile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			u = ptt.maskBathyYZ(u, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			v = readVariable(vfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			v = ptt.maskBathyYZ(v, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)



			# Load wind and stress

			uw = readVariable(uwindfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, var2D=True)

			vw = readVariable(vwindfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, var2D=True)

			us = readVariable(ustressfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, var2D=True)

			vs = readVariable(vstressfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, var2D=True)



			#==



			# Get along-slope surface and undercurrent flow speeds. 



			slope_uv = u * np.sin(bearing_nX[xi]*np.pi/180.) + v * np.cos(bearing_nX[xi]*np.pi/180.)

			# Mask along-slope velocity before computing average.

			surf_uv_av[ti-t_start,xi] = np.ma.mean(slope_uv[0])

			surf_uv_max[ti-t_start,xi] = np.ma.max(slope_uv[0])

			slope_uv = np.ma.masked_where(rho<1028., slope_uv)

			ZZZ = np.tile(grid.RC.squeeze()[zz[0]:zz[1]], (slope_uv.shape[1],1)).T

			slope_uv = np.ma.masked_where(ZZZ<-800, slope_uv)



			slope_uv_av[ti-t_start,xi] = np.ma.mean(slope_uv)

			slope_uv_max[ti-t_start,xi] = np.ma.max(slope_uv)

			#lim = 0.09; slope_uv = tools.boundData(slope_uv, -lim, lim)

			#slope_uv[0,0] = -lim; slope_uv[0,1] = lim

			#pt.plot1by2([bathy, slope_uv], contour=[None, rho], X=[X,Y[yy[0]:yy[1],0]], Y=[Y,Z[zz[0]:zz[1]]], vlines=[[X[0,xi]],None])

			# Get along-slope wind and wind stress.

			slope_uw = uw * np.sin(bearing_nX[xi]*np.pi/180.) + vw * np.cos(bearing_nX[xi]*np.pi/180.)

			slope_us = us * np.sin(bearing_nX[xi]*np.pi/180.) + vs * np.cos(bearing_nX[xi]*np.pi/180.)

			slope_uw_av[ti-t_start,xi] = np.ma.mean(slope_uw)

			slope_us_av[ti-t_start,xi] = np.ma.mean(slope_us)



		#plt.contourf(bathy)

		#plt.plot(3.e3*slope_uv_av[ti]); plt.show()



	#==



	np.save('surf_uv_max', surf_uv_max)

	np.save('surf_uv_av', surf_uv_av)

	np.save('slope_uv_av', slope_uv_av)

	np.save('slope_uv_max', slope_uv_max)

	np.save('slope_uw_av', slope_uw_av)

	np.save('slope_us_av', slope_us_av)



	#for ti in [1,3,5,7,9,11]:#range(nT):

	#       plt.plot(slope_uv_av[ti], label=ti)

	#plt.legend()

	#plt.show()



	quit()

	# Get y-grid point closest to slope_y_nX

	# Read u, v at this point and n grid points north and south.

	# Get component of each in along-slope direction





#==



# This block computes and saves the spatial-mean wind stress curl evaluated over the continental shelf.

onShelfWindCurl = True

if onShelfWindCurl:



	ustressfile = 'oceTAUX'; vstressfile = 'oceTAUY'



	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy



	X = grid.XC; Y = grid.YC

	dx = grid.DXG; dy = grid.DYG



	pathno = 0; xshift = 0; yshift = 0; latsi = []; lonsi = []

	subr = 1

	if subr:

		# Get subregion of bathymetry.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		dx = tools.getSubregionXY(dx, latsi, lonsi)

		dy = tools.getSubregionXY(dy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		pathno = 0

		xshift = lonsi[0]

		yshift = latsi[0]

	nY, nX = X.shape



	# Get lat, lon lists of continental slope

	slope_x, slope_y = PAS_tools.getSlopeContour(bathy, X, Y, pathno)

	np.save('slope_x', slope_x); np.save('slope_y', slope_y)



	bathy = ptt.maskBathyXY(bathy, grid, zi=0, subregion=True, lats=latsi, lons=lonsi)



	bearing = PAS_tools.getBearing(slope_x, slope_y)



	slope_x_nX = np.zeros(nX); slope_y_nX = np.zeros(nX)

	slope_xi = np.zeros(nX, dtype=int); slope_yi = np.zeros(nX, dtype=int)

	bearing_nX = np.zeros(nX)



    # For each X grid point, get the lat and lon of slope contour's nearest point.

	for xi in range(nX):

		slope_xi[xi] = int(np.argmin(np.abs(slope_x-X[0,xi])))

		slope_x_nX[xi] = slope_x[slope_xi[xi]]

		slope_y_nX[xi] = slope_y[slope_xi[xi]]

		bearing_nX[xi] = bearing[slope_xi[xi]]

		slope_yi[xi] = int(np.argmin(np.abs(slope_y_nX[xi]-Y[:,xi])))

	

	# Load wind and stress

	#us = readVariable(ustressfile, PASDIR, file_format='nc', meta=False, var2D=True)

	#vs = readVariable(vstressfile, PASDIR, file_format='nc', meta=False, var2D=True)

	

	us = bathy.copy()

	vs = bathy.copy()

	

	# Get wind stress curl

	curl = tools.ddx(vs, dx) - tools.ddy(us, dy)



	# Array demarking shelf region. 1s for shelf, 0s for deep ocean or land.

	shelf = np.where(bathy>=0, 0, 1)

	for xi in range(nX):

		shelf[slope_yi[xi]:, xi] = 0

		

	pt.plot1by1(shelf)

	

	quit()

	

		



#==



readSlopeUVfile = False

if readSlopeUVfile:



	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779

	

	nl2 = int(15*12); nl = 2*nl2 + 1

	lags = np.linspace(-nl2,nl2,nl)

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	uvfile = 'slope_uv_av.npy'

	surfuvfile = 'surf_uv_max.npy'

	timefile = 'PAS_time.npy'

	slopexfile = 'slope_x.npy'

	slopeyfile = 'slope_y.npy'

	

	slope_uv = np.load(path+uvfile)

	surf_uv = np.load(path+surfuvfile)

	time = np.load(path+timefile)

	slope_x = np.load(path+slopexfile)

	slope_y = np.load(path+slopeyfile)

	print(slope_uv.shape)

	#plt.scatter(slope_x, slope_y); plt.show(); quit()



	

	#===

	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

	monthDec = {}

	for mi in range(len(months)):

		monthDec[months[mi]] = mi/12.

	year = [float(ctime(int(time_))[-4:])-15+monthDec[ctime(int(time_))[4:7]] for time_ in time[t_start:t_end]]#

	

	start = 24*12 + 5*12; end=-11

	year = year[start:end] 



	slope_uv = slope_uv[start:end]

	surf_uv = surf_uv[start:end]

		

	#year_ = year[nn//2:-nn//2+1]

	# Silvano years: 1987 - 2017 (+-2.5)

	#ysSil = np.argmin(np.abs(np.array(year)-1987.+2.5)); yeSil = np.argmin(np.abs(np.array(year)-2017.-2.5))

	#ys_Sil = np.argmin(np.abs(np.array(year_)-1987.+2.5)); ye_Sil = np.argmin(np.abs(np.array(year_)-2017.-2.5))

	#ys = ysSil; ye = yeSil

	#ys_ = ys_Sil; ye_ = ye_Sil

	

	plotSections = ['westGetz', 'westPITW', 'westPITE']

	#plotSections = ['westPITW']



	# Plot average evaluated over all sections.

	if 1:

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

		

		#uv_mean = uv_mean[ys:ye]; surf_uv_mean = surf_uv_mean[ys:ye]; year = year[ys:ye]	



		#==

		

		# Demean

		uv_mean = PAS_tools.demean(uv_mean)

		surf_uv_mean = PAS_tools.demean(surf_uv_mean)

		

		# Deseason

		uv_mean = PAS_tools.deseason(uv_mean)

		surf_uv_mean = PAS_tools.deseason(surf_uv_mean)

		

		# Detrend

		uv_mean = PAS_tools.detrend(uv_mean, year)

		surf_uv_mean = PAS_tools.detrend(surf_uv_mean, year)

		

		#==

		

		corruv = PAS_tools.crossCorr(uv_mean, surf_uv_mean, nl, nl2, P_VALUE=False)

		corruu = PAS_tools.crossCorr(uv_mean, uv_mean, nl, nl2, P_VALUE=False)

		corrvv = PAS_tools.crossCorr(surf_uv_mean, surf_uv_mean, nl, nl2, P_VALUE=False)

		Fcuv, freq, period = PAS_tools.fft(corruv)

		Fcuu, freq, period = PAS_tools.fft(corruu)

		Fcvv, freq, period = PAS_tools.fft(corrvv)

		

		Puv = PAS_tools.fftPower(Fcuv, nl2)

		Puu = PAS_tools.fftPower(Fcuu, nl2)

		Pvv = PAS_tools.fftPower(Fcvv, nl2)

		coh = Puv**2 / (Puu * Pvv)

		

		# Get correlations after a range of running means applied to data.

		nw = 61

		window_lengths = np.linspace(1,nw,nw)

		corrw, tw, uw, vw = PAS_tools.windowCorr(uv_mean, surf_uv_mean, window_lengths, year, return_uv=True)

		corrwl, p_value = PAS_tools.crossCorrWindowAv(uv_mean, surf_uv_mean, window_lengths, nl, nl2)



		lim = 0.8; corrwlp = tools.boundData(corrwl, -lim, lim); corrwlp[-1,-1] = -lim; corrwlp[-1,-2] = lim

		plt.contourf(lags/12, window_lengths/12, corrwlp, cmap='bwr', vmin=-lim, vmax=lim); plt.colorbar()

		plt.xlabel('Lag (years)'); plt.ylabel('Running-mean window length (years)')

		plt.title('corr(surface current, undercurrent)')

		plt.grid();	plt.show()

		

		plt.pcolor(lags/12, window_lengths/12, p_value, cmap='bwr', vmin=0.0, vmax=1); plt.colorbar()

		plt.xlabel('Lag (years)'); plt.ylabel('Running-mean window length (years)')

		plt.title('corr(surface current, undercurrent)')

		plt.grid();	plt.show()

		

		quit()

				

		# Plot correlations for velocities averaged with range of running mean windows.

		wis = [0,12,18,60]

		plt.figure(figsize=(16,17))

		for wi in range(len(wis)):

			plt.subplot(2,2,wi+1)

			plt.plot(tw[wis[wi]], uw[wis[wi]], color='r', label='deep vel')

			plt.plot(tw[wis[wi]], vw[wis[wi]], color='k', label='surf vel')

			plt.grid()

			title = str(wis[wi]/12) + '-year running mean, corr = {:.2f}'.format(corrw[wis[wi]])

			plt.title('u, v time series, ' + title)

			plt.legend()

		plt.show()

		

		#==

		

		nn = 60

		uv_mean = PAS_tools.windowAv(uv_mean, n=nn)[nn//2:-nn//2+1]

		surf_uv_mean = PAS_tools.windowAv(surf_uv_mean, n=nn)[nn//2:-nn//2+1]

		year = PAS_tools.windowAv(year, n=nn)[nn//2:-nn//2+1]

		

		plt.figure(figsize=(16,17))

		plt.subplot(221)

		plt.plot(year, uv_mean, label='deep vel', color='r')

		plt.plot(year, surf_uv_mean, label='surf vel', color='k')

		plt.title('Slope surf/deep velocity time series'); plt.legend(); plt.grid()

		plt.subplot(222)

		plt.plot((window_lengths-1)/12, corrw);  plt.xlabel('running mean window size (years)')

		plt.title('corr(u,v) vs running mean window size'); plt.grid()

		

		plt.subplot(223)

		for wi in range(len(wis)):

			plt.plot(lags/12., corrwl[wis[wi],:], label=str(wis[wi]) + ' month window')

		plt.grid(); plt.xlabel('Lag (years)'); plt.title('corr(u,v) vs lag'); plt.legend()

		plt.subplot(224)

		#plt.plot(period, np.real(Fcuv), label='Real power'); plt.gca().set_xscale('symlog', base=2)

		#plt.plot(period, np.imag(Fcuv), label='Imag power'); plt.legend()

		plt.plot(period[nl2:],Puv); plt.gca().set_xscale('log', base=2)

		plt.xticks([2,4,12,24,36,48,60,120])

		plt.gca().get_xaxis().set_major_formatter(ticker.ScalarFormatter())

		plt.xlabel('Period (months)');

		plt.title('corr(u,v) Fourier power spectrum'); plt.grid();

		plt.show()

	

		#==

		

		plt.plot(year, uv_mean, label='Deep along-slope flow', color='r')

		plt.plot(year, surf_uv_mean, label='Surface along-slope flow', color='k')		

		#plt.plot(year_[ys_:ye_], uv_mean[ys_:ye_], label='Deep along-slope flow')

		#plt.plot(year_[ys_:ye_], surf_uv_mean[ys_:ye_], label='Surface along-slope flow')

		plt.title('Along-slope average current speed (m/s)')

		#plt.ylim(0, 0.15)

		plt.grid()

		plt.legend()

		plt.show()

		

	# Plot average for each section.

	if 0:

		for section in plotSections:

			iw = sections[section][0]; ie = sections[section][1]

			uv = np.mean(slope_uv[:,iw:ie], axis=1)

			uv = window_average(uv, n=nn)[nn//2:-nn//2+1]

			plt.plot(year_, uv, label=section)

		plt.title('Along-slope average current speed (m/s)')

		#plt.ylim(0, 0.15)

		plt.grid()

		plt.legend()

		plt.show()

	

	quit()

	

#==



readAllSlopeFiles = True

if readAllSlopeFiles:

	

	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779

	

	nl2 = int(15*12); nl = 2*nl2 + 1

	lags = np.linspace(-nl2,nl2,nl)

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	uvfile = 'slope_uv_max.npy'

	surfuvfile = 'surf_uv_max.npy'

	uwfile = 'slope_uw_av.npy'

	usfile = 'slope_us_av.npy'

	timefile = 'PAS_time.npy'

	slopexfile = 'slope_x.npy'

	slopeyfile = 'slope_y.npy'

	

	slope_uv = np.load(path+uvfile)

	surf_uv = np.load(path+surfuvfile)

	uw = np.load(path+uwfile)

	us = np.load(path+usfile)

	time = np.load(path+timefile)

	slope_x = np.load(path+slopexfile)

	slope_y = np.load(path+slopeyfile)

	print(slope_uv.shape)

	#plt.scatter(slope_x, slope_y); plt.show(); quit()



	## This code for test of correlation functions.

	#for ti in range(slope_uv.shape[0]):

	#	slope_uv[ti,:] = np.sin(2*np.pi*ti/120.)

	#	surf_uv[ti,:] = np.sin(2*np.pi*(ti+15)/120.) 

	

	#===

	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

	monthDec = {}

	for mi in range(len(months)):

		monthDec[months[mi]] = mi/12.

	year = [float(ctime(int(time_))[-4:])-15+monthDec[ctime(int(time_))[4:7]] for time_ in time[t_start:t_end]]#

	

	start = 24*12 + 5*12; end=-11

	year = year[start:end] 



	slope_uv = slope_uv[start:end]

	surf_uv = surf_uv[start:end]

	uw = uw[start:end]

	us = us[start:end]

	

	#===============

	#===============



	# Get zonal sums and averages evaluated over sections of interest.

	plotSections = ['westGetz', 'westPITW', 'westPITE']

	#plotSections = ['westGetz']

		

	uv_sum = np.zeros(len(year))

	surf_uv_sum = np.zeros(len(year))

	uw_sum = np.zeros(len(year))

	us_sum = np.zeros(len(year))

	nx = 0

	for section in plotSections:

		iw = sections[section][0]; ie = sections[section][1]

		uv_sum += np.sum(slope_uv[:,iw:ie], axis=1)

		surf_uv_sum += np.sum(surf_uv[:,iw:ie], axis=1)

		uw_sum += np.sum(uw[:,iw:ie], axis=1)

		us_sum += np.sum(us[:,iw:ie], axis=1)

		nx += ie-iw

	uv_mean = uv_sum / nx

	surf_uv_mean = surf_uv_sum / nx

	uw_mean = uw_sum / nx

	us_mean = us_sum / nx

	

	#==

	

	# Demean, deseason, detrend.

	

	# Demean

	uv_mean = PAS_tools.demean(uv_mean)

	surf_uv_mean = PAS_tools.demean(surf_uv_mean)

	uw_mean = PAS_tools.demean(uw_mean)

	us_mean = PAS_tools.demean(us_mean)

	

	# Deseason

	uv_mean = PAS_tools.deseason(uv_mean)

	surf_uv_mean = PAS_tools.deseason(surf_uv_mean)

	uw_mean = PAS_tools.deseason(uw_mean)

	us_mean = PAS_tools.deseason(us_mean)

	

	# Detrend

	uv_mean = PAS_tools.detrend(uv_mean, year)

	surf_uv_mean = PAS_tools.detrend(surf_uv_mean, year)

	uw_mean = PAS_tools.detrend(uw_mean, year)

	us_mean = PAS_tools.detrend(us_mean, year)

	

	#==

	

	# Get correlations after a range of running means applied to data.

	nw = 61

	window_lengths = np.linspace(1,nw,nw)

	corr_surfDeep, p_value_surfDeep = PAS_tools.crossCorrWindowAv(surf_uv_mean, uv_mean, window_lengths, nl, nl2)

	corr_deepWind, p_value_deepWind = PAS_tools.crossCorrWindowAv(uv_mean, uw_mean, window_lengths, nl, nl2)

	corr_surfWind, p_value_surfWind = PAS_tools.crossCorrWindowAv(surf_uv_mean, uw_mean, window_lengths, nl, nl2)

	corr_deepStress, p_value_deepStress = PAS_tools.crossCorrWindowAv(uv_mean, us_mean, window_lengths, nl, nl2)

	corr_surfStress, p_value_surfStress = PAS_tools.crossCorrWindowAv(surf_uv_mean, us_mean, window_lengths, nl, nl2)

	

	corrs = [corr_surfDeep, corr_surfWind, corr_deepWind, corr_deepStress]#, corr_surfStress]

	pvals = [p_value_surfDeep, p_value_surfWind, p_value_deepWind, p_value_deepStress]#, p_value_surfStress]

	titles = ['corr(surface current, undercurrent)', 'corr(surface current, surface wind speed)', 'corr(undercurrent, surface wind speed)', 'corr(undercurrent, surface stress)']#, 'corr(surface current, surface stress)']

	xlabels = ['', '', 'Lag (years)', 'Lag (years)']

	ylabels = ['Running-mean window length (years)', '', 'Running-mean window length (years)', '']



	plt.figure(figsize=(10,9), dpi=200)

	if len(plotSections)>1:

		plt.gcf().suptitle('All sections', fontsize=9)

	else: 

		plt.gcf().suptitle(section, fontsize=9)

	lim = .8;

	for ci, corr in enumerate(corrs):

		corr_ = tools.boundData(corr, -lim, lim); corr_[-1,-1] = -lim; corr_[-1,-2] = lim

		plt.subplot(2,2,ci+1)

		plt.contourf(lags/12, window_lengths/12, corr_, cmap='bwr', vmin=-lim, vmax=lim); plt.colorbar()

		#CS = plt.contour(lags/12, window_lengths/12, pvals[ci], levels=[0.01,0.05], colors='k')

		#plt.gca().clabel(CS, inline=True, fontsize=10)

		ptt.doStippling(pvals[ci], dataLim=0.1, X=lags/12, Y=window_lengths/12, dx=5, dy=2, s=0.5)		

		plt.xlabel(xlabels[ci], fontsize=8); plt.ylabel(ylabels[ci], fontsize=6)

		plt.xticks(fontsize=5); plt.yticks(fontsize=5)

		plt.title(titles[ci], fontsize=8)

		plt.grid(linewidth=0.5)

	#plt.tight_layout()

	plt.show()

	

	#==

	

	nn = 60

	uv_mean = PAS_tools.windowAv(uv_mean, n=nn)[nn//2:-nn//2+1]

	surf_uv_mean = PAS_tools.windowAv(surf_uv_mean, n=nn)[nn//2:-nn//2+1]

	uw_mean = PAS_tools.windowAv(uw_mean, n=nn)[nn//2:-nn//2+1]

	us_mean = PAS_tools.windowAv(us_mean, n=nn)[nn//2:-nn//2+1]

	year = PAS_tools.windowAv(year, n=nn)[nn//2:-nn//2+1]

		

	uv_mean /= np.max(np.abs(uv_mean))

	surf_uv_mean /= np.max(np.abs(surf_uv_mean))

	uw_mean /= np.max(np.abs(uw_mean))

	us_mean /= np.max(np.abs(us_mean))

	

	plt.plot(year, uv_mean, label='Deep along-slope flow', color='r')

	plt.plot(year, surf_uv_mean, label='Surface along-slope flow', color='k')		

	plt.plot(year, uw_mean, label='Along-slope wind')

	plt.plot(year, us_mean, label='Along-slope stress')

	

	plt.title('Along-slope average current speed (m/s)')

	#plt.ylim(0, 0.15)

	plt.grid()

	plt.legend()

	plt.show()

		

	quit()

	

#==



saveTime = False

if saveTime:



	rhofile = 'RHOAnoma'

	rho = readVariable(rhofile, PASDIR, file_format='nc', meta=True)

	time = np.array(rho['TIME'][:])

	np.save('PAS_time', time)

	print(time.shape)

	print(time)

	quit()



#==



animateX = False

if animateX:



	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	X = grid.XC; Y = grid.YC

	Z = grid.RC.squeeze()

	

	pathno = 0;	xshift = 0; yshift = 0; latsi = []; lonsi = []

	subr = 1

	if subr:

		# Get subregion of bathymetry.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		pathno = 0

		xshift = lonsi[0]

		yshift = latsi[0]

	nY, nX = X.shape



	# Get lat, lon lists of continental slope

	slope_x, slope_y = PAS_tools.getSlopeContour(bathy, X, Y, pathno)

	bathy = ptt.maskBathyXY(bathy, grid, zi=0, subregion=True, lats=latsi, lons=lonsi)

	

	# Get bearing of slope (using a centred differencing, one-sided at ends of slope).

	bearing = np.zeros(len(slope_x))

	bearing[0] = PAS_tools.getBearing(slope_y[0], slope_x[0], slope_y[1], slope_x[1])

	bearing[-1] = PAS_tools.getBearing(slope_y[-2], slope_x[-2], slope_y[-1], slope_x[-1])

	bearing[1:-1] = PAS_tools.getBearing(slope_y[:-2], slope_x[:-2], slope_y[2:], slope_x[2:])

	# bearing[i] from slope[i-1] and slope[i+1] 

	

	slope_x_nX = np.zeros(nX); slope_y_nX = np.zeros(nX)

	slope_xi = np.zeros(nX, dtype=int); slope_yi = np.zeros(nX, dtype=int)

	bearing_nX = np.zeros(nX)

	

	# For each X grid point, get the lat and lon of slope contour's nearest point.

	for xi in range(nX):

		slope_xi[xi] = int(np.argmin(np.abs(slope_x-X[0,xi])))

		slope_x_nX[xi] = slope_x[slope_xi[xi]]

		slope_y_nX[xi] = slope_y[slope_xi[xi]]

		bearing_nX[xi] = bearing[slope_xi[xi]]		

		slope_yi[xi] = int(np.argmin(np.abs(slope_y_nX[xi]-Y[:,xi])))



	# Plot bathy, check contour points and bearing.

	#plt.subplot(121); plt.contourf(X, Y, bathy); plt.scatter(slope_x_nX, slope_y_nX)

	#plt.subplot(122); plt.plot(X[0], bearing_nX); plt.show(); quit()

	

	# For each timestep, load instances of u, v, rho.

	# Depth average beneath a rho contour. (27.95?)

	# Then take dot prodcut with unit vector in direction of slope's bearing.

	# Surface flow can be just the flow or estimated from SSH for geostrophic part.

	

	zs = [0, -1000]

	zz = grid.getIndexFromDepth(zs)



	imgpath = '/home/michael/Documents/Python/mitgcmPy/tmpimg/'

	ufile = 'UVEL'; vfile = 'VVEL';	rhofile = 'RHOAnoma'

	for ti in range(1):

		

		for xi in range(nX):

			print(xi)

			

			yy = [yshift+slope_yi[xi]-3, yshift+slope_yi[xi]+3]

			print(yy)

			rho = 1028.5 + readVariable(rhofile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			rho = ptt.maskBathyYZ(rho, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			rho = tools.boundData(rho, 1026, 1040)

			

			u = readVariable(ufile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			u = ptt.maskBathyYZ(u, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			v = readVariable(vfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			v = ptt.maskBathyYZ(v, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			print(yy)

			slope_uv = u * np.sin(bearing_nX[xi]*np.pi/180.) + v * np.cos(bearing_nX[xi]*np.pi/180.)

			

			bathy = tools.boundData(bathy, -4.e3, -4.e2)

			lim = 0.06;	slope_uv = tools.boundData(slope_uv, -lim, lim)

			slope_uv[0,0] = -lim; slope_uv[0,1] = lim

			

			title0 = 'Bathymetry (m), xi='+str(xi) + ', ' + f'bearing = {bearing_nX[xi]:.1f}'

			pt.plot1by2([bathy, slope_uv], contour=[None, rho], X=[X,Y[yy[0]-yshift:yy[1]-yshift,0]], Y=[Y,Z[zz[0]:zz[1]]], vlines=[[X[0,xi]],None], titles=[title0, 'Along-slope flow speed (m/s)'], show=False, save=True, outname=imgpath+f"{xi:03}", cmaps=['viridis', 'coolwarm'])

			

		#anim = animation.FuncAnimation(fig, animate, interval=50, frames=nX)

		#anim.save('PAS_animateX',metadata={'artist':'Guido'},writer='ffmpeg',fps=8,bitrate=-1)

		#plt.draw()

		#quit()

	# Get y-grid point closest to slope_y_nX

	# Read u, v at this point and n grid points north and south.

	# Get component of each in along-slope direction

	

	import cv2

	import os



	imgpath = '/home/michael/Documents/Python/mitgcmPy/tmpimg/'

	video_name = 'video.avi'



	images = [img for img in os.listdir(imgpath) if img.endswith(".png")]

	images.sort()

	frame = cv2.imread(os.path.join(imgpath, images[0]))

	height, width, layers = frame.shape



	video = cv2.VideoWriter(video_name, 0, 8, (width,height))



	for image in images:

		video.write(cv2.imread(os.path.join(imgpath, image)))



	cv2.destroyAllWindows()

	video.release()



	quit()

	

	

#==



# Define latitude of continental slope for each longitude.

SLOPE_CONTOUR = False

if SLOPE_CONTOUR:



	grid = Grid_PAS(DATADIR + 'PAS_851/run/')

	bathy = grid.bathy

	X = grid.XC; Y = grid.YC

	

	pathno = 0

	subr = 1

	if subr:

		# Get subregion of bathymetry.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		pathno = 1



	plt.subplot(121)

	cs = plt.contour(X, Y, bathy, levels=[-1000])

	x = cs.collections[0].get_paths()[pathno].vertices[:,0][::-1]

	y = cs.collections[0].get_paths()[pathno].vertices[:,1][::-1]

	#plt.scatter(x[:30], y[:30])

	plt.grid()

	

	print(len(x))

	print(X.shape)

	

	bs = []

	for i in range(len(x)-1):

		bs.append(PAS_tools.getBearing(y[i],x[i], y[i+1],x[i+1]))

	plt.subplot(122)

	plt.plot(bs)

	plt.grid()

	plt.show()

	quit() 

	

	

#==



latLonToCartesian1 = False

if latLonToCartesian1:



	grid = Grid_PAS(DATADIR + 'PAS_851/run/')

	bathy = grid.bathy

	X = grid.XC; Y = grid.YC

	

	#x, y, dx, dy = PAS_tools.latlon_to_xy(X,Y)

	b = PAS_tools.getBearing(1,1,1,2)

	print(b); quit()

	

	pt.plot1by2([x,y]);quit()

	

	pt.plot1by1(bathy, X=x, Y=y);quit() 

	cs = plt.contour(x, y, bathy, levels=[-1000])

	x = cs.collections[0].get_paths()[0].vertices[:,0]

	y = cs.collections[0].get_paths()[0].vertices[:,1]

	plt.show()

	plt.clf()

	

#==



plotSections = False

if plotSections:



	t_start = 107; t_end = 502	

	path = '/home/michael/Documents/Python/mitgcmPy/movies/31_05_23/'

	path = path + str(t_start) + '_' + str(t_end) + '/'

	

	slopexfile = 'slope_x.npy';	slopeyfile = 'slope_y.npy'

	slopexifile = 'slope_xi.npy';	slopeyifile = 'slope_yi.npy'

	slope_x = np.load(path+slopexfile);	slope_y = np.load(path+slopeyfile)

	slope_xi = np.load(path+slopexifile); slope_yi = np.load(path+slopeyifile)

	

	grid = Grid_PAS(PASDIR)

	X = grid.XC; Y = grid.YC

	bathy = grid.bathy

	bathy = ptt.maskBathyXY(bathy, grid, zi=0)

	bathy = tools.boundData(bathy, -4000, -400)

	

	#plt.scatter(slope_xi, slope_yi); plt.show(); 

	

	subr = 1

	if subr:

		# Get subregion of bathymetry.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		bathy = ptt.maskBathyXY(bathy, grid, zi=0, subregion=True, lats=latsi, lons=lonsi)

		xshift = lonsi[0]; yshift = latsi[0]



	plotSections = ['westGetz', 'westPITW', 'westPITE']

	plt.contourf(X, Y, bathy)

	for sec in plotSections:

		iw = sections[sec][0]; ie = sections[sec][1]

		plt.plot([X[0,iw], X[0,iw]], [Y[slope_yi[iw]-10,0], Y[slope_yi[iw]+10,0]], linestyle='dashed', color='k')

		plt.plot([X[0,ie], X[0,ie]], [Y[slope_yi[ie]-10,0], Y[slope_yi[ie]+10,0]], linestyle='dashed', color='k')

	plt.colorbar()

	plt.tight_layout()

	plt.show()

	

	

#==



latLonToCartesian2 = False

if latLonToCartesian2:



	grid = Grid_PAS(DATADIR + 'PAS_851/run/')

	bathy = grid.bathy

	lon = grid.XC; lat = grid.YC

	

	dA, x, y = PAS_tools.dA_from_latlon (lon, lat, periodic=False, return_edges=True)

	

	

	pt.plot1by2([x, y])

	

#==



coherence = False

if coherence:



	nt = 20*12

	nl2 = 10*12

	nl = 2*nl2 + 1

	

	t = np.linspace(1,nt,nt)

	lags = np.linspace(-nl2,nl2,nl)

	

	rand1 = np.random.random(nt); rand2 = np.random.random(nt); rand3 = np.random.random(nt); rand4 = np.random.random(nt)

	rand1 = 1; rand2 = 1; rand3 = 1; rand4 = 1

	p1 = 4*12; p2 = 4*12; p3 = 2.6

	u = .5*np.sin(2*np.pi*t/p1) + 0*rand1*np.cos(2*np.pi*t/(rand2*p3))

	v = -.5*np.cos(2*np.pi*t/p2) + 0*rand3*np.cos(2*np.pi*t/(rand4*p3))

	

	#plt.plot(u); plt.plot(v); plt.show()



	#==

	

	# Get correlations

	

	# Compute coherence.

	# 1. Get auto and cross correlations as functions of lag.

	# 2. Fourier transform correlations.



	corruv = PAS_tools.crossCorr(u, v, nl, nl2)

	corruu = PAS_tools.crossCorr(u, u, nl, nl2)

	corrvv = PAS_tools.crossCorr(v, v, nl, nl2)

	corrvu = PAS_tools.crossCorr(v, u, nl, nl2)

	#coh = PAS_tools.coherence(corruv, corruu, corrvv)

	

	#plt.plot(lags, corruv); plt.plot(lags, corrvu); plt.grid(); plt.show(); quit()

	

	#==

	

	# Get correlations after a range of running means applied to data.

	nw = 60

	window_lengths = np.linspace(1,nw,nw)

	corrw, tw, uw, vw = PAS_tools.windowCorr(u, v, window_lengths, t, return_uv=True)

	#==

	

	wis = [0,12,18,48]

	plt.figure(figsize=(16,17))

	for wi in range(len(wis)):

		plt.subplot(2,2,wi+1)

		plt.plot(tw[wis[wi]], uw[wis[wi]])

		plt.plot(tw[wis[wi]], vw[wis[wi]])

		plt.xlabel('time (years)'); plt.grid()

		title = str(wis[wi]/12) + '-year running mean, corr = {:.2f}'.format(corrw[wis[wi]])

		plt.title('u, v time series, ' + title)

	plt.show()



	#Fv = np.fft.fftshift(np.fft.fft(u)); freqt = 1./np.fft.fftshift(np.fft.fftfreq(nt))

	#plt.plot(freqt, np.real(Fv), label='Real'); plt.plot(freqt, np.imag(Fv), label='Imag'); plt.legend(); plt.show(); quit()

	

	Fcorruv = np.fft.fftshift(np.fft.fft(corruv))

	freq = np.fft.fftshift(np.fft.fftfreq(nl))

	

	plt.figure(figsize=(16,17))

	plt.subplot(221)

	plt.plot(t[:]/12, u[:], label='u = sin'+str(p1/12) + '+ cos'+str(p3/12))

	plt.plot(t[:]/12, v[:], label='v = -sin'+str(p1/12) + '+ cos'+str(p3/12))

	plt.xlabel('time (years)'); plt.title('u, v time series'); plt.legend(); plt.grid()

	plt.subplot(222)

	plt.plot((window_lengths-1)/12, corrw);  plt.xlabel('running mean window size (years)')

	plt.title('corr(u,v) vs running mean window size'); plt.grid()

	

	plt.subplot(223)

	plt.plot(lags/12., corruv); plt.grid(); plt.xlabel('Lag (years)'); plt.title('corr(u,v) vs lag');

	plt.subplot(224)

	plt.plot(1./12./freq, np.real(Fcorruv), label='Real power')

	plt.plot(1./12./freq, np.imag(Fcorruv), label='Imag power'); plt.legend()

	plt.xlabel('Period (years)'); plt.title('corr(u,v) Fourier power spectrum'); plt.grid();



	plt.show()

	

	

	quit()

	

# OPTIONS:

# 1.

# For each point in scatter plot, get nearest grid point.

# Is there worry about double-counting some u,v with this method?



# 2.

# For each lon, get latitude of shelf slope and nearest grid point for u, v.

# 

	

	

	

	

	

	

	

	

# Original code for getting coordinates of continental slope.

#x = []; y = []

#for item in cs.collections:

#	print(1)

#	for i in item.get_paths():

#		v = i.vertices; x.append(v[:, 0]); y.append(v[:, 1])
