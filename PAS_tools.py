# PAS_tools.py



# tools.py



import sys



import matplotlib.pyplot as plt



from time import ctime



import numpy as np

from sklearn.linear_model import LinearRegression

from scipy.stats import pearsonr

	

#====================================================================================



# General helper function to get the area of each cell from latitude and longitude arrays giving the coordinates of the cell centres. Adapted from Kaitlin Naughten's mitgcm_python.

def dA_from_latlon (lon, lat, periodic=False, return_edges=False):



	deg2rad = np.pi / 180.

	rEarth = 6.371e6



	# Make sure they're 2D

	if len(lon.shape) == 1 and len(lat.shape) == 1:

		lon, lat = np.meshgrid(lon, lat)

	# Now make the edges

	# Longitude

	if periodic:

		lon_extend = wrap_periodic(lon, is_lon=True)

		lon_edges = 0.5*(lon_extend[:,:-1] + lon_extend[:,1:])

	else:

		lon_edges_mid = 0.5*(lon[:,:-1] + lon[:,1:])

		# Extrapolate the longitude boundaries

		lon_edges_w = 2*lon_edges_mid[:,0] - lon_edges_mid[:,1]

		lon_edges_e = 2*lon_edges_mid[:,-1] - lon_edges_mid[:,-2]

		lon_edges = np.concatenate((lon_edges_w[:,None], lon_edges_mid, lon_edges_e[:,None]), axis=1)

	dlon = lon_edges[:,1:] - lon_edges[:,:-1] 

	# Latitude

	lat_edges_mid = 0.5*(lat[:-1,:] + lat[1:,:])

	lat_edges_s = 2*lat_edges_mid[0,:] - lat_edges_mid[1,:]

	lat_edges_n = 2*lat_edges_mid[-1,:] - lat_edges_mid[-2,:]

	lat_edges = np.concatenate((np.expand_dims(lat_edges_s,0), lat_edges_mid, np.expand_dims(lat_edges_n,0)), axis=0)

	dlat = lat_edges[1:,:] - lat_edges[:-1,:]

	# Now convert to Cartesian

	dx = rEarth*np.cos(lat*deg2rad)*dlon*deg2rad

	dy = rEarth*dlat*deg2rad

	dA = dx*dy

	if return_edges:

		return dA, lon_edges, lat_edges

	else:

		return dA

		

#==



def latlon_to_xy(X, Y, returnGridSpacing=True, deg2rad=np.pi/180., rEarth=6.371e6):

	'''Convert lat-lon (X,y) grid to Cartesian (x, y) grid using the Haversine formula.'''

	

	X *= deg2rad; Y *= deg2rad

	

	nY, nX = X.shape

	x = np.zeros((nY, nX)); y = np.zeros((nY, nX))

	dx = np.zeros((nY-1, nX-1)); dy = np.zeros((nY-1, nX-1))

	for i in range(1,nX):

		for j in range(1,nY):

			dX = X[j,i] - X[j,i-1]

			dY = Y[j,i] - Y[j-1,i]

			

			# Compute grid for x

			a = np.cos(Y[j,i])**2 * np.sin(dX/2)**2

			dx[j-1,i-1] = 2 * rEarth * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 

			x[j,i] = x[j,i-1] + dx[j-1,i-1]

			

			# Repeat for y

			a = np.sin(dY/2)**2

			dy[j-1,i-1] = 2 * rEarth * np.arctan2(np.sqrt(a), np.sqrt(1-a)) 

			y[j,i] = y[j-1,i] + dy[j-1,i-1] 

	

	if returnGridSpacing:

		return x, y, dx, dy

	else:

		return x, y

		

#==



def getBearing_(lat1, lon1, lat2, lon2, deg2rad=np.pi/180.):

	

	lat1 = deg2rad * np.copy(lat1)

	lat2 = deg2rad * np.copy(lat2)

	lon1 = deg2rad * np.copy(lon1)

	lon2 = deg2rad * np.copy(lon2)



	dLon = (lon2 - lon1)



	y = np.sin(dLon) * np.cos(lat2)

	x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dLon)



	brng = np.arctan2(y, x)



	brng = np.rad2deg(brng)



	return brng

	

#==



def getBearing(slope_x, slope_y):



	bearing = np.zeros(len(slope_x))

	bearing[0] = getBearing_(slope_y[0], slope_x[0], slope_y[1], slope_x[1])

	bearing[-1] = getBearing_(slope_y[-2], slope_x[-2], slope_y[-1], slope_x[-1])

	bearing[1:-1] = getBearing_(slope_y[:-2], slope_x[:-2], slope_y[2:], slope_x[2:])

	

	return bearing

	

#==



def getSlopeContour(bathy, X, Y, pathno, level=-1000):



	cs = plt.contour(X, Y, bathy, levels=[level])

	x = cs.collections[0].get_paths()[pathno].vertices[:,0][::-1]

	y = cs.collections[0].get_paths()[pathno].vertices[:,1][::-1]

	plt.clf()

	

	return x, y



#==



def getDecimalTime(t, t_start=None, t_end=None, PAS_correction=15):

	'''Given time array t of seconds since epoch, return decimalised time in units years.'''

	

	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

	monthDec = {}

	

	for mi in range(len(months)):

		monthDec[months[mi]] = mi/12.

	

	year = [float(ctime(int(t_))[-4:])-PAS_correction + monthDec[ctime(int(t_))[4:7]] for t_ in t[t_start:t_end]]



	return year



#==



def demean(f, axis=0):

	'''Remove the mean from a field along the given axis.'''



	u = f.copy()

	

	up = np.zeros((u.shape))

	dim = up.shape[axis]

	

	mean = np.mean(u, axis=axis)



	if axis == 0:

		for i in range(dim):

			up[i, ] = u[i, ] - mean



	elif axis == 1:

		for i in range(dim):

			up[:, i, ] = u[:, i, ] - mean

	

	elif axis == 2:

		for i in range(dim):

			up[:, :, i] = u[:, :, i] - mean



	return up

	

#==



def desurf(f, axis=1):

	'''Subtract surface flow from all levels of field f.'''

	

	u = f.copy()

	

	up = np.zeros((u.shape))

	dim = up.shape[axis]



	if axis == 0:

		for i in range(dim):

			up[i, ] = u[i, ] - u[0, ]



	elif axis == 1:

		for i in range(dim):

			up[:, i, ] = u[:, i, ] - u[:, 0, ]

	

	elif axis == 2:

		for i in range(dim):

			up[:, :, i] = u[:, :, i] - u[:, :, 0]

	

	return up

	

#==

	

def corr(u, v, axis=0):

	'''Get the correlation between two fields along axis.'''



	# corr = cov(u) * cov(v) / (sd(u) * sd(v))



	ud = demean(u, axis)

	vd = demean(v, axis)



	return np.mean(ud * vd, axis=axis) / (np.mean(ud**2, axis=axis) * np.mean(vd**2, axis=axis))**0.5



#==



def pearson(u, v):

	'''Use scipy.pearson to compute the correlation.'''

	

	corr, pval = pearsonr(u, v)

	

	return corr



#==



def corr2(u, v, axis=0):



	uav = np.mean(u, axis=axis)

	vav = np.mean(v, axis=axis)

	

	return np.mean((u-uav) * (v-vav), axis=axis) / (np.mean((u-uav)**2, axis=axis) *  np.mean((v-vav)**2, axis=axis))**0.5

	

#==



def cov(u, v, axis=0):

	'''Get the covariance between two fields along axis.'''



	ud = demean(u, axis)

	vd = demean(v, axis)



	return np.mean(ud * vd, axis=axis)

	

#==



def crossCorr(u, v, nl, nl2, P_VALUE=True):

	'''Get cross correlations for number of lags nl.'''

	

	if nl != 2 * nl2 + 1:

		print('PAS_tools.crossCorr error: nl and nl2 incompatible.')

		quit()

		

	corruv = np.zeros(nl)

	p_value = np.zeros(nl)

	

	#corruv[nl2] = corr(u,v)

	corruv[nl2], p_value[nl2] = pearsonr(u,v)

	

	for li in range(1,nl2+1):

		# Shift u ahead of v

		#corruv[nl2-li] = corr(u[li:], v[:-li])

		corruv[nl2-li], p_value[nl2-li] = pearsonr(u[li:], v[:-li])

		# Shift v ahead of u (u, v[li:])

		#corruv[nl2+li] = corr(u[:-li], v[li:])	

		corruv[nl2+li], p_value[nl2+li] = pearsonr(u[:-li], v[li:])	

	

	if P_VALUE:

		return corruv, p_value

		

	else:

		return corruv

	

#==



def movingAv(a, n=3):

	ret = np.ma.cumsum(a.copy(), dtype=float, axis=0)

	ret[n:] = ret[n:] - ret[:-n]

	return ret[n - 1:] / n

	

#==



def windowAv(a, n):

		b = np.copy(a)

		b = np.ma.masked_where(b!=b, b)

		if n % 2 == 0:

			n += 1

		for i in range((n-1)//2, len(b)-(n-1)//2+1):

			b[i] = np.ma.mean(b[i-(n-1)//2:i+(n-1)//2], axis=0)

		return b 

		

#==



def windowAvLoop(a, n):



	b = np.zeros(a.shape)

	nt, ny, nx = b.shape

	for j in range(ny):

		print(j)

		for i in range(nx):

			b[:,j,i] = windowAv(a[:,j,i], n)

	

	return b



#==



def windowCorr(u, v, window_lengths, t, return_uv=False):

	'''Compute correlations between time series u and v average applying 

	running means for a range of running-mean window lengths.'''



	uw = []; vw = []; tw = []

	corrw = np.zeros(len(window_lengths))

	

	tw.append(t.copy()); uw.append(u); vw.append(v)

	corrw[0] = corr(u, v)

	

	for wi in range(1, len(window_lengths)):

		utmp = movingAv(u, n=wi)

		vtmp = movingAv(v, n=wi)

		

		corrw[wi] = corr(utmp, vtmp)

		

		uw.append(utmp); vw.append(vtmp)

		tw.append(movingAv(t, n=wi))

		

	if return_uv:

		return corrw, tw, uw, vw		

	else:

		return corrw, tw

		

#==



def fft(u):

	'''Return fast Fourier transform of u and corresponding frequencies.'''

	

	Fu = np.fft.fftshift(np.fft.fft(u))

	freq = np.fft.fftshift(np.fft.fftfreq(len(u)))

	period = 1. / freq

	

	return Fu, freq, period



#==



def fftPower(Fu, nl2):

	'''From fft of u, get power for each absolute frequency.'''

	

	F = np.abs(np.copy(Fu))



	power = np.real(F[nl2:])

	power[1:] += np.real(F[:nl2][::-1])

	

	return power

	

#==



def crossCorrWindowAv(uin, vin, window_lengths, nl, nl2):

	'''Plot correlation between (u,v) for a range of lags AND running-mean window lengths.'''

	

	# For each running-mean window length, compute cross correlation for all lags.

	

	u = np.copy(uin); v = np.copy(vin)

	

	corrw = np.zeros((len(window_lengths), nl))

	p_value = np.zeros((len(window_lengths), nl))

	

	corrw[0], p_value[0] = crossCorr(u, v, nl, nl2)

	

	for wi in range(1, len(window_lengths)):

		utmp = movingAv(u, n=wi)

		vtmp = movingAv(v, n=wi)

		

		corrw[wi], p_value[wi] = crossCorr(utmp, vtmp, nl, nl2)

		

	return corrw, p_value



#==



def seasonalData(data, year):

	'''Return dictionary of timeseries data averaged over four seasons.'''

	

	jani = np.argmin(year[:12]-np.floor(year[:12]))

	if (jani - 1) % 3 != 0:

		print('Adjust start so that data starts with 3-months of a season')

		quit()

	if len(year) % 12 != 0:

		print('Adjust end so that data spans integer number of years.')

		quit()



	nyears = int(len(year)/12)

	seasons = ['Summer', 'Autumn', 'Winter', 'Spring']

	shift = int((jani-1)/3)

	seasons = seasons[-shift:] + seasons[:-shift]

	

	if len(data.shape) > 1:

		# Assume data has x,y dims

		nt, ny, nx = data.shape

		season0 = np.zeros((nyears, ny, nx)); season1 = np.zeros((nyears, ny, nx))

		season2 = np.zeros((nyears, ny, nx)); season3 = np.zeros((nyears, ny, nx))

	else:

		season0 = np.zeros(nyears); season1 = np.zeros(nyears)

		season2 = np.zeros(nyears); season3 = np.zeros(nyears)

		

	for yi in range(nyears):

		season0[yi] = np.mean(data[12*yi+0:12*yi+3], axis=0)

		season1[yi] = np.mean(data[12*yi+3:12*yi+6], axis=0)

		season2[yi] = np.mean(data[12*yi+6:12*yi+9], axis=0)

		season3[yi] = np.mean(data[12*yi+9:12*yi+12], axis=0)

	

	data_seasonal = {seasons[0]:season0, seasons[1]:season1, seasons[2]:season2, seasons[3]:season3}

	

	return data_seasonal



#==



def seasonalDataIPO(data, year, IPO):

	'''Return dictionary of timeseries data averaged over four seasons.

	Also return seasonal timeseries averaged during positive and negative IPO.'''



	std = np.mean(IPO**2)**0.5

			

	jani = np.argmin(year[:12]-np.floor(year[:12]))

	if (jani - 1) % 3 != 0:

		print('Adjust start so that data starts with 3-months of a season')

		quit()

	if len(year) % 12 != 0:

		print('Adjust end so that data spans integer number of years.')

		quit()



	nyears = int(len(year)/12)

	seasons = ['Summer', 'Autumn', 'Winter', 'Spring']

	shift = int((jani-1)/3)

	seasons = seasons[-shift:] + seasons[:-shift]

	

	# Assume data has x,y dims

	nt, ny, nx = data.shape

	seasonData = np.zeros((5, 3, ny, nx))

	nPos = np.zeros(4);	nNeg = np.zeros(4)

		

	for yi in range(nyears):

		for si in range(4):

			for mi in range(3):

				seasonData[si, 0, ] += data[12*yi + 3*si + mi, ]

				if IPO[12*yi + 3*si + mi] > std:

					seasonData[si, 1, ] += data[12*yi + 3*si + mi, ]

					nPos[si] += 1

				elif IPO[12*yi + 3*si + mi] < - std:

					seasonData[si, 2, ] += data[12*yi + 3*si + mi, ]

					nNeg[si] += 1

					

	seasonData[:,0] /= (0.25 * nt)



	seasonData[4,0] = np.mean(seasonData[:,0], axis=0)

	seasonData[4,1] = np.sum(seasonData[:,1], axis=0) / np.sum(nPos) - seasonData[4,0]

	seasonData[4,2] = np.sum(seasonData[:,2], axis=0) / np.sum(nNeg) - seasonData[4,0]

	

	for si in range(4):

		seasonData[si,1] /= nPos[si]

		seasonData[si,1] -= seasonData[si,0]

		#seasonData[si,1] -= np.mean(seasonData[:,0], axis=0)

		

		seasonData[si,2] /= nNeg[si]

		seasonData[si,2] -= seasonData[si,0]

		#seasonData[si,2] -= np.mean(seasonData[:,0], axis=0)

	

	return seasonData, seasons

	

#==



def monthlyDataIPO(data, year, IPO):

	'''Return dictionary of timeseries data averaged over each month.

	Also return monthly timeseries averaged during positive and negative IPO.'''



	std = np.mean(IPO**2)**0.5

			

	jani = np.argmin(year[:12]-np.floor(year[:12]))

	if len(year) % 12 != 0:

		print('Adjust end so that data spans integer number of years.')

		quit()



	nyears = int(len(year)/12)

	months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'July', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

	months = months[-int(jani):] + months[:-int(jani)]

	

	# Assume data has x,y dims

	nt, ny, nx = data.shape

	monthData = np.zeros((13, 3, ny, nx))

	nPos = np.zeros(12); nNeg = np.zeros(12)

		

	# For each year, for each month, add data to all-month sum.

	# Check if IPO pos/neg, and data to corresponding sum.

	for yi in range(nyears):

		for mi in range(12):

			monthData[mi, 0, ] += data[12*yi + mi, ]

			if IPO[12*yi + mi] > std:

				monthData[mi, 1, ] += data[12*yi + mi, ]

				nPos[mi] += 1

			elif IPO[12*yi + mi] < - std:

				monthData[mi, 2, ] += data[12*yi + mi, ]

				nNeg[mi] += 1

					

	# Get mean by dividing by number of years.

	monthData[:,0] /= (nt / 12.)



	# Compute data averaged over all moonths.

	monthData[12,0] = np.mean(monthData[:,0], axis=0)

	monthData[12,1] = np.sum(monthData[:,1], axis=0) / np.sum(nPos) - monthData[12,0]

	monthData[12,2] = np.sum(monthData[:,2], axis=0) / np.sum(nNeg) - monthData[12,0]

	

	# Normalise monthly data, and remove monthly mean.

	for mi in range(12):

		monthData[mi,1] /= nPos[mi]

		monthData[mi,1] -= monthData[mi,0]

		

		monthData[mi,2] /= nNeg[mi]

		monthData[mi,2] -= monthData[mi,0]

	

	return monthData, months

	

#==		



def surfaceRadiation(qa, Ua, Ta, LW, SW, SSS, hi, hs, bathy, draft, nTiS=1000):

	'''Return all terms in surface radiation balance.

	Requires computation of surface temperature, also returned.'''

	

	qa = qa.copy()

	Ua = Ua.copy()

	Ta = Ta.copy()

	LW = LW.copy()

	SW = SW.copy()

	SSS = SSS.copy()

	hi = hi.copy()

	hs = hs.copy()

	

	aa1 = 2663.5

	aa2 = 12.537

	bb1 = 0.622

	bb2 = 1.0 - bb1

	Ppascals = 100000.0

	lnTEN = np.log(10)

	

	SEAICE_lhEvap = 2.5e6

	SEAICE_lhFusion = 3.34e5

	Ls = SEAICE_lhEvap + SEAICE_lhFusion # lhSublim

	

	rhoa = 1.2 # SEAICE_rhoAir

	rhos = 3.3e2 # SEAICE_rhoSnow

	ca = 1.005e3 # SEAICE_cpAir

	CDi = 1.75e-3 # SEAICE_dalton

	

	eps_s = 9.5e-1 # Sea-ice/snow emissivity. (Same val. for snow/ice.)

	sigma = 5.67e-8 # SEAICE_boltzmann

	

	SEAICE_dTempFrz_dS = -5.75e-2

	SEAICE_tempFrz0 = 9.01e-2

	celsius2K = 273.15



	Tf = SEAICE_dTempFrz_dS * SSS + SEAICE_tempFrz0 + celsius2K

	

	nY, nX = qa.shape

	TiS = np.zeros((nTiS, nY, nX))

	for yi in range(nY):

		for xi in range(nX):

			TiS[:,yi,xi] = np.linspace(250, 270, nTiS)

	

	#==

	

	# LATENT

	

	mm_log10pi = -aa1 / TiS + aa2

	mm_pi = np.exp(mm_log10pi * lnTEN)

	qsat = bb1 * mm_pi / ( Ppascals - (1.0 - bb1) * mm_pi ) # qhice

	# qhice=qsat: saturation vapor pressure of snow/ice surface

	

	LATENT = rhoa * Ls * CDi * Ua * (qsat - qa)



	#==

	

	# SENS



	SENS = rhoa * ca * CDi * Ua * (TiS - Ta)



	#==

	

	# BLACKBODY

	

	BB = eps_s * sigma * TiS**4

	

	#==

	

	# LONGWAVE

	LWtmp = - eps_s * LW

	LW = np.zeros((nTiS, nY, nX))

	for ti in range(0,nTiS):

		LW[ti] = LWtmp

	

	#==

	

	# SHORTWAVE

	

	# DESIGNED TO WORK FOR WINTER ONLY, WHEN SEA ICE IS SNOW-COVERED.

	penetSWFrac = 0.0

	

	# These are values for DRY snow/ice, suitable for winter when Ta < surface melt temp.

	ALB_ICE = 7.5e-1

	ALB_SNOW = 8.4e-1

	

	# For albedo should compute a linear transition between ALB_SNOW and ALB_ICE where hsnow>HCUT=0.15.

	# But just using ALB_SNOW for now.

	ALB = ALB_SNOW

     

	SWtmp = - (1.0 - ALB) * (1.0 - penetSWFrac) * SW    

	SW = np.zeros((nTiS, nY, nX))

	for ti in range(0,nTiS):

		SW[ti] = SWtmp

		

	#==

	

	# CONDUCTIVE HEAT FLUX

	

	ki = 2.16560

	ks = 3.10000e-1

	

	FC = ki * ks * (Tf - TiS) / (ki * hs + ks * hi)

	

	#==

	

	# Find TiS that balances surface radiation:

	# Take argmin of poly and linearly interpolate.

	

	poly = LATENT + SENS + LW + BB + SW - FC

	

	SW_out = SW[0]; LW_out = LW[0]

	

	LATENT_out = np.zeros((nY, nX))

	SENS_out = np.zeros((nY, nX))

	BB_out = np.zeros((nY, nX))

	FC_out = np.zeros((nY, nX))

	

	tsurf = np.zeros((nY, nX))

	for yi in range(nY):	

		for xi in range(nX):

			p = poly[:, yi, xi]

			T = TiS[:, yi, xi]

			t_arg = np.argmin(np.abs(p))

			

			if bathy[yi,xi] < 0 and draft[yi,xi] == 0:

			

				if np.min(p)*np.max(p)>0:

					print('Polynomial has no solutions at (xi, yi) = ' +str(xi,yi))

									

				else:

					tsurf[yi,xi] = T[t_arg]

					LATENT_out[yi, xi] = LATENT[t_arg, yi, xi]

					SENS_out[yi, xi] = SENS[t_arg, yi, xi]

					BB_out[yi, xi] = BB[t_arg, yi, xi]

					FC_out[yi, xi] = FC[t_arg, yi, xi]



	# End loops

	#==

			

	return {'LATENT':LATENT_out, 'SENS':SENS_out, 'BB':BB_out, 'LW':LW_out, 'SW':SW_out, 'FC':FC_out, 'TiS':tsurf}



#==



def surfaceRadiation_testCode():



	SW_out = SW[0]; LW_out = LW[0]

	

	LATENT_out = np.zeros((nY, nX))

	SENS_out = np.zeros((nY, nX))

	BB_out = np.zeros((nY, nX))

	FC_out = np.zeros((nY, nX))

	

	tsurf = np.zeros((nY, nX))

	tsurf2 = np.zeros((nY, nX))

	for yi in range(nY):	

		for xi in range(nX):

			p = poly[:, yi, xi]

			T = TiS[:, yi, xi]

			

			if bathy[yi,xi] < 0 and draft[yi,xi] == 0:

			

				if np.min(p)*np.max(p)>0:

					print(hs[yi,xi])

					print(hi[yi,xi])

					

					plt.subplot(121)

					plt.plot(T, p, label='poly')

					plt.plot(T, FC[:,yi,xi], label='FC')

					#plt.plot(T, SENS[:,yi,xi], label='SENS')

					#plt.plot(T, LW[:,yi,xi], label='LW')

					#plt.plot(T, SW[:,yi,xi], label='SW')

					#plt.plot(T, BB[:,yi,xi], label='BB')

					#plt.plot(T, LATENT[:,yi,xi], label='LATENT')

					plt.grid();	plt.legend()

					plt.subplot(122)

					plt.contourf(bathy, vmin=-1000, vmax=0, mesh=True); plt.plot(xi,yi, 'ro')

					plt.show()

									

				else:



					t_arg = np.argmin(np.abs(p))

					if p[t_arg] > 0:

						tsurf[yi,xi] = T[t_arg] + (T[t_arg] - T[t_arg-1]) / (p[t_arg] - p[t_arg-1])

					else:

						tsurf[yi,xi] = T[t_arg] - p[t_arg] * (T[t_arg] - T[t_arg+1]) / (p[t_arg] - p[t_arg+1])

					tsurf2[yi,xi] = T[t_arg]

					

					LATENT_out[yi, xi] = LATENT[t_arg, yi, xi]

					SENS_out[yi, xi] = SENS[t_arg, yi, xi]

					BB_out[yi, xi] = BB[t_arg, yi, xi]

					FC_out[yi, xi] = FC[t_arg, yi, xi]

			

	#pt.plot1by2([tsurf, tsurf2], vmin=[250]*2, vmax=[300]*2, mesh=True)

	

	RADS2 = {'LATENT':LATENT_out, 'SENS':SENS_out, 'BB':BB_out, 'LW':LW_out, 'SW':SW_out, 'FC':FC_out, 'TiS':tsurf2}



	return

	

#==



def deseason(uin):

	'''Deseasons monthly data. Checks if data has 12*n data points.'''

	

	u = np.copy(uin)

	

	nm = len(u)

	if nm % 12 != 0:

		print('Warning: PAS_tools.deseason. Length of uin should be multiple of 12.')

	

	cycle = np.zeros(12)

	for mi in range(12):

		cycle[mi] = np.mean(u[mi:None:12])

	

	cycle = np.tile(cycle, nm//12)

	

	if len(u.shape) > 1:

		return np.transpose(np.transpose(u, (2,1,0)) - cycle, (2,1,0))

	else:

		return u - cycle

		

	

#==



def detrend(uin, tin, interceptFlag=1):

	'''Use linear regression to detrend uin timeseries.'''

	

	if tin is None:

		tin = np.linspace(1,len(uin),len(uin))

	

	u = np.copy(uin)

	t = np.copy(tin).reshape(-1,1)

	

	reg = LinearRegression().fit(t, u)

	

	return uin - (reg.coef_*tin + interceptFlag * reg.intercept_)

	

#==



def detrendXY(uin, tin, interceptFlag=1):

	'''Use linear regression to detrend uin timeseries,

	where uin has dimensions of (time, lat, lon).'''	

	

	nt, ny, nx = uin.shape

	

	uout = uin.copy()

	

	if tin is None:

		tin = np.linspace(1,nt,nt)

	t = np.copy(tin).reshape(-1,1)

	

	for j in range(ny):

		for i in range(nx):

			reg = LinearRegression().fit(t, uin[:,j,i])

			uout[:,j,i] = uin[:,j,i] - (reg.coef_*tin + interceptFlag * reg.intercept_)

			

	return uout

	

#==



def getCoriolis(lats, Omega=2*np.pi/86400., nx=None):

	'''Get Coriolis coefficient for spherical polar grid.

	Input lats to be in degrees.'''

	

	if len(lats.shape) == 2 or nx is None:

		f = 2. * Omega * np.sin(lats * np.pi / 180.)	

			

	elif len(lats.shape) == 1:

		f = np.zeros((lats.shape[0],nx))

		for i in range(nx):

			f[:,i] = 2. * Omega * np.sin(lats * np.pi / 180.)

	

	else:

		print('PAS_tools.getCoriolis error: lats must have 1 or 2 dimensions.')

		quit()

			

	return f



	

	

	

	

	



