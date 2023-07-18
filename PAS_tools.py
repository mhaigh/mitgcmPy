# PAS_tools.py



# tools.py



import sys



import matplotlib.pyplot as plt



import numpy as np

from sklearn.linear_model import LinearRegression



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

	

def corr(u, v, axis=0):

	'''Get the correlation between two fields along axis.'''



	# corr = cov(u) * cov(v) / (sd(u) * sd(v))



	ud = demean(u, axis)

	vd = demean(v, axis)



	return np.mean(ud * vd, axis=axis) / (np.mean(ud**2, axis=axis) * np.mean(vd**2, axis=axis))**0.5



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

	

	from scipy.stats import pearsonr

	

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

	

def cross_correlation(x, y, lag_range):

    """

    Computes the cross-correlation between two time series for a range of lags.



    Args:

        x (ndarray): The first time series.

        y (ndarray): The second time series.

        lag_range (range or list): The range of lags to compute the cross-correlation.



    Returns:

        ndarray: An array of cross-correlation values for each lag.

    """

    cross_corr = np.zeros(len(lag_range))



    for i, lag in enumerate(lag_range):

        if lag >= 0:

            cross_corr[i] = np.correlate(x[lag:], y[:-lag])

        else:

            cross_corr[i] = np.correlate(x[:lag], y[-lag:])



    return cross_corr

	

#==



def movingAv(a, n=3):

	ret = np.cumsum(a.copy(), dtype=float)

	ret[n:] = ret[n:] - ret[:-n]

	return ret[n - 1:] / n

	

#==



def windowAv(a, n):

		b = np.copy(a)

		if n % 2 == 0:

			n += 1

		for i in range((n-1)//2, len(a)-(n-1)//2+1):

			b[i] = np.mean(a[i-(n-1)//2:i+(n-1)//2])

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

	

	return u - cycle

	

#==



def detrend(uin, tin):

	'''Use linear regression to detrend uin timeseries.'''

	

	u = np.copy(uin)

	t = np.copy(tin).reshape(-1,1)

	

	reg = LinearRegression().fit(t, u)

	

	return uin - (reg.coef_*tin + reg.intercept_)

	

	

	

	

	

	



