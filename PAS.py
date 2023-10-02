import numpy as np



from time import ctime



from scipy.stats import pearsonr



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



#==





MAIN = False

if MAIN:



	OUTPATH = 'slopeCurrent/'



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



	bearing = PAS_tools.getBearing(slope_x, slope_y)

	plt.plot(bearing); plt.show(); quit()

	

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

		

	np.save(OUTPATH+'slope_x', slope_x); np.save(OUTPATH+'slope_y', slope_y)

	np.save(OUTPATH+'slope_xi', slope_xi); np.save(OUTPATH+'slope_yi', slope_yi)



	# Plot bathy, check contour points and bearing.

	#plt.subplot(121); plt.contourf(X, Y, bathy); plt.scatter(slope_x_nX, slope_y_nX)

	#plt.subplot(122); plt.plot(X[0], bearing_nX); plt.show(); quit()

	

	# For each timestep, load instances of u, v, rho.

	# Depth average beneath a rho contour. (27.95?)

	# Then take dot prodcut with unit vector in direction of slope's bearing.

	# Surface flow can be just the flow or estimated from SSH for geostrophic part.

	

	zs = [0, -1000]

	zz = grid.getIndexFromDepth(zs)



	ufile = 'UVEL'; vfile = 'VVEL';	rhofile = 'RHOAnoma'

	uwindfile = 'EXFuwind'; vwindfile = 'EXFvwind'

	ustressfile = 'oceTAUX'; vstressfile = 'oceTAUY'	



	# How many grid points either side of slope do we search for max in along-slope velocity?

	yslim = 3; ynlim = 1



	#t_start = 0; t_end = 10

	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779



	surf_uv_max = np.zeros((t_end-t_start, nX))

	surf_uv_av = np.zeros((t_end-t_start, nX))

	slope_uv_av = np.zeros((t_end-t_start, nX))

	slope_uv_max = np.zeros((t_end-t_start, nX))	

	slope_uw_av = np.zeros((t_end-t_start, nX))

	slope_us_av = np.zeros((t_end-t_start, nX))

	surf_uv_maxj = np.zeros((t_end-t_start, nX))

	slope_uv_maxj = np.zeros((t_end-t_start, nX))

	slope_uw_maxj = np.zeros((t_end-t_start, nX))

	slope_us_maxj = np.zeros((t_end-t_start, nX))

	maxj_out = np.zeros((t_end-t_start, nX), dtype=int)

	for ti in range(t_start, t_end):

		print(ti)

		print(ctime(readVariable(rhofile, PASDIR, file_format='nc', meta=True)['TIME'][ti]))

		for xi in range(nX):

			

			yy = [yshift+slope_yi[xi]-yslim, yshift+slope_yi[xi]+ynlim]

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

			surf_uv = slope_uv[0]

			surf_uv_av[ti-t_start,xi] = np.ma.mean(surf_uv)

			surf_uv_max[ti-t_start,xi] = np.ma.max(surf_uv)

			slope_uv = np.ma.masked_where(rho<1028., slope_uv)

			ZZZ = np.tile(grid.RC.squeeze()[zz[0]:zz[1]], (slope_uv.shape[1],1)).T

			slope_uv = np.ma.masked_where(ZZZ<-800, slope_uv)

	

			slope_uv_av[ti-t_start,xi] = np.ma.mean(slope_uv)

			slope_uv_max[ti-t_start,xi] = np.ma.max(slope_uv)		

					

			maxz, maxj = np.unravel_index(slope_uv.argmax(), slope_uv.shape)

			surf_uv_maxj[ti-t_start,xi] = surf_uv[maxj]

			slope_uv_maxj[ti-t_start,xi] = slope_uv[maxz,maxj]

			maxj_out[ti-t_start,xi] = maxj + slope_yi[xi] - yslim



			#lim = 0.09; slope_uv = tools.boundData(slope_uv, -lim, lim)

			#slope_uv[0,0] = -lim; slope_uv[0,1] = lim

			#pt.plot1by2([bathy, slope_uv], contour=[None, rho], X=[X,Y[yy[0]:yy[1],0]], Y=[Y,Z[zz[0]:zz[1]]], vlines=[[X[0,xi]],None])

			# Get along-slope wind and wind stress.

			slope_uw = uw * np.sin(bearing_nX[xi]*np.pi/180.) + vw * np.cos(bearing_nX[xi]*np.pi/180.)

			slope_us = us * np.sin(bearing_nX[xi]*np.pi/180.) + vs * np.cos(bearing_nX[xi]*np.pi/180.)

			slope_uw_av[ti-t_start,xi] = np.ma.mean(slope_uw)

			slope_us_av[ti-t_start,xi] = np.ma.mean(slope_us)



			slope_uw_maxj[ti-t_start,xi] = slope_uw[maxj]

			slope_us_maxj[ti-t_start,xi] = slope_us[maxj]



			#==

		

		#plt.contourf(bathy)

		#plt.plot(3.e3*slope_uv_av[ti]); plt.show()

			

	#==

	

	OUTPATH = 'slopeCurrent/'

	np.save(OUTPATH+'surf_uv_max', surf_uv_max)

	np.save(OUTPATH+'surf_uv_av', surf_uv_av)

	np.save(OUTPATH+'surf_uv_maxj', surf_uv_maxj)

	np.save(OUTPATH+'slope_uv_av', slope_uv_av)

	np.save(OUTPATH+'slope_uv_max', slope_uv_max)

	np.save(OUTPATH+'slope_uv_maxj', slope_uv_maxj)

	np.save(OUTPATH+'slope_uw_av', slope_uw_av)

	np.save(OUTPATH+'slope_us_av', slope_us_av)

	np.save(OUTPATH+'slope_uw_maxj', slope_uw_maxj)

	np.save(OUTPATH+'slope_us_maxj', slope_us_maxj)

	np.save(OUTPATH+'maxj', maxj_out)





	#for ti in [1,3,5,7,9,11]:#range(nT):

	#	plt.plot(slope_uv_av[ti], label=ti)

	#plt.legend()

	#plt.show()

			

	quit()



#==



# Computes vertical profiles of fields of interest.

# E.g., cross-slope salinity gradient, temp, velocity.

MAINZ = False

if MAINZ:



	OUTPATH = 'slopeCurrent/'



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

		X, Y, dX, dY = grid.XYsubr(EASlons, EASlats, dxdy=True)

		pathno = 0

		xshift = lonsi[0]

		yshift = latsi[0]

	nY, nX = X.shape



	# Get lat, lon lists of continental slope

	slope_x, slope_y = PAS_tools.getSlopeContour(bathy, X, Y, pathno)



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

		

	zs = [0, -1000]

	zz = grid.getIndexFromDepth(zs)

	nZ = zz[1] - zz[0]



	ufile = 'UVEL'; vfile = 'VVEL';	rhofile = 'RHOAnoma'

	Sfile = 'SALT'; Tfile = 'THETA'

	Pxfile = 'Um_dPhiX'; Pyfile = 'Vm_dPhiY'

	

	# How many grid points either side of slope do we search for max in along-slope velocity?

	yslim = 3; ynlim = 1



	#t_start = 0; t_end = 10

	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779



	uv_z = np.zeros((t_end-t_start, nZ, nX))

	Sgrad = np.zeros((t_end-t_start, nZ, nX))

	Tgrad = np.zeros((t_end-t_start, nZ, nX))

	Pgrad = np.zeros((t_end-t_start, nZ, nX))

	maxj_out = np.zeros((t_end-t_start, nX))

	maxz_out = np.zeros((t_end-t_start, nX))

	for ti in range(t_start, t_end):

		print(ti)

		print(ctime(readVariable(rhofile, PASDIR, file_format='nc', meta=True)['TIME'][ti]))

		for xi in range(nX):

			

			ys = yshift+slope_yi[xi]-yslim

			yn = yshift+slope_yi[xi]+ynlim

			yy = [ys, yn]

			

			rho = 1028.5 + readVariable(rhofile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			rho = ptt.maskBathyYZ(rho, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			rho = tools.boundData(rho, 1026, 1040)

			

			# Load velocities

			u = readVariable(ufile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			u = ptt.maskBathyYZ(u, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			v = readVariable(vfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			v = ptt.maskBathyYZ(v, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)



			#==

		

			# Get along-slope surface and undercurrent flow speeds. 

			slope_uv = u * np.sin(bearing_nX[xi]*np.pi/180.) + v * np.cos(bearing_nX[xi]*np.pi/180.)

			

			# Mask along-slope velocity before computing average.

			slope_uv = np.ma.masked_where(rho<1028., slope_uv)

			ZZZ = np.tile(grid.RC.squeeze()[zz[0]:zz[1]], (slope_uv.shape[1],1)).T

			slope_uv = np.ma.masked_where(ZZZ<-800, slope_uv)



			# Find latitude (maxj) where slope_uv is maximised.	

			maxz, maxj = np.unravel_index(slope_uv.argmax(), slope_uv.shape)

			uv_z[ti-t_start,:,xi] = slope_uv[:,maxj]

			

			maxj_out[ti-t_start,xi] = maxj

			maxz_out[ti-t_start,xi] = maxz

						

			# Get centred difference for S, T and P gradients. 

			#Sn = readVariable(Sfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=ys+maxj+1, zz=zz)

			#Ss = readVariable(Sfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=ys+maxj-1, zz=zz)

			#Se = readVariable(Sfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift+1, yy=ys+maxj, zz=zz)

			#Sw = readVariable(Sfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift-1, yy=ys+maxj, zz=zz)

			#Sy = 0.5 * (Sn - Ss) / dY[ys-yshift+maxj, xi]

			#Sx = 0.5 * (Se - Sw) / dX[ys-yshift+maxj, xi]

			

			#Tn = readVariable(Tfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=ys+maxj+1, zz=zz)

			#Ts = readVariable(Tfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=ys+maxj-1, zz=zz)

			#Te = readVariable(Tfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift+1, yy=ys+maxj, zz=zz)

			#Tw = readVariable(Tfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift-1, yy=ys+maxj, zz=zz)

			#Ty = 0.5 * (Tn - Ts) / dY[ys-yshift+maxj, xi]

			#Tx = 0.5 * (Te - Tw) / dX[ys-yshift+maxj, xi]

			

			Px = readVariable(Pxfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=ys+maxj, zz=zz)

			Py = readVariable(Pyfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=ys+maxj, zz=zz)

			

			#lim = 0.09; slope_uv = tools.boundData(slope_uv, -lim, lim)

			#slope_uv[0,0] = -lim; slope_uv[0,1] = lim

			#pt.plot1by2([bathy, slope_uv], contour=[None, rho], X=[X,Y[yy[0]:yy[1],0]], Y=[Y,Z[zz[0]:zz[1]]], vlines=[[X[0,xi]],None])

			

			# Get along-slope wind and wind stress.

			#Sgrad[ti,:,xi] = Sy * np.sin(bearing_nX[xi]*np.pi/180.) - Sx * np.cos(bearing_nX[xi]*np.pi/180.)

			#Tgrad[ti,:,xi] = Ty * np.sin(bearing_nX[xi]*np.pi/180.) - Tx * np.cos(bearing_nX[xi]*np.pi/180.)

			Pgrad[ti,:,xi] = Py * np.sin(bearing_nX[xi]*np.pi/180.) - Px * np.cos(bearing_nX[xi]*np.pi/180.)



			#==

		

	#==

	

	OUTPATH = 'slopeCurrent/'

	np.save('uv_z', np.ma.filled(uv_z, fill_value=0))

	#np.save('Tgrad', np.ma.filled(Tgrad, fill_value=0))

	#np.save('Sgrad', np.ma.filled(Sgrad, fill_value=0))

	np.save('Pgrad', np.ma.filled(Pgrad, fill_value=0))



	np.save('maxj', maxj)

	np.save('maxz', maxz)



	quit()



#==



slopeSalinity = False

if slopeSalinity:



	OUTPATH = 'slopeCurrent/'



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



	slope_x_nX = np.zeros(nX); slope_y_nX = np.zeros(nX)

	slope_xi = np.zeros(nX, dtype=int); slope_yi = np.zeros(nX, dtype=int)

	

	# For each X grid point, get the lat and lon of slope contour's nearest point.

	for xi in range(nX):

		slope_xi[xi] = int(np.argmin(np.abs(slope_x-X[0,xi])))

		slope_x_nX[xi] = slope_x[slope_xi[xi]]

		slope_y_nX[xi] = slope_y[slope_xi[xi]]

		slope_yi[xi] = int(np.argmin(np.abs(slope_y_nX[xi]-Y[:,xi])))

		

	# For each timestep, load instances salinity, S.

	# For each x, load S over given y-range centred on continental slope.

	# Put loaded S into 2-array whose new y axis is centred on slope.

	# Repeat for bathy.

		

	Sfile = 'SALT'	

	Tfile = 'THETA'	

	

	# How many grid points either side of slope?

	yslim = 50; ynlim = 20

	nY = ynlim + yslim



	zs = [0, -1000]

	zz = grid.getIndexFromDepth(zs)

	nZ = zz[1] - zz[0]

	

	#t_start = 0; t_end = 10

	t_start = 0; t_end = 779



	slopeSalt = np.zeros((t_end-t_start, nZ, nY, nX))

	slopeTheta = np.zeros((t_end-t_start, nZ, nY, nX))

	for ti in range(t_start, t_end):

		print(ti)

		for xi in range(nX):

			

			yy = [yshift+slope_yi[xi]-yslim, yshift+slope_yi[xi]+ynlim]

			

			S = readVariable(Sfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			S = ptt.maskBathyYZ(S, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			slopeSalt[ti, :, :, xi] = S[:, :]

			

			T = readVariable(Tfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			T = ptt.maskBathyYZ(T, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			slopeTheta[ti, :, :, xi] = T[:, :]

	

	#==

	

	# Now get bathy centred at slope. 

	slopeBathy = np.zeros((nY, nX))

	for xi in range(nX):

		slopeBathy[:, xi] = bathy[slope_yi[xi]-yslim:slope_yi[xi]+ynlim, xi]

		

	# Save

	np.save(OUTPATH+'slopeSalt.npy', slopeSalt)

	np.save(OUTPATH+'slopeTheta.npy', slopeTheta)

	np.save(OUTPATH+'slopeBathy.npy', slopeBathy)



	quit()



#==



# This block computes and saves the spatial-mean wind stress curl evaluated over the continental shelf.

onShelfWindCurl = False

if onShelfWindCurl:



	ustressfile = 'oceTAUX'; vstressfile = 'oceTAUY'

	SIfile = 'SIarea'; FWflxfile = 'oceFWflx'



	iw = 0; ie = westPITE[1]



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

		X, Y = grid.XYsubr(EASlons, EASlats)

		pathno = 0

		xshift = lonsi[0]

		yshift = latsi[0]

	nY, nX = X.shape



	# Get lat, lon lists of continental slope

	slope_x, slope_y = PAS_tools.getSlopeContour(bathy, X, Y, pathno)



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

	

	#SI = readVariable(SIfile, PASDIR, file_format='nc', meta=False, var2D=True)

	#FWflx = readVariable(FWflxfile, PASDIR, file_format='nc', meta=False, var2D=True)



	us = np.zeros((2, grid.bathy.shape[0], grid.bathy.shape[1]))

	us[0] = grid.bathy; us[1] = grid.bathy

	vs = us.copy(); SI = us.copy(); FWflx = us.copy()

	

	# Get wind stress curl

	curl = tools.ddx(vs, dx) - tools.ddy(us, dy)

	

	# Array demarking shelf region. 0s land or far deep ocean, 1s for shelf, 2s for near deep ocean.

	shelf = np.zeros(bathy.shape)

	shelfSI = np.zeros(bathy.shape)

	for xi in range(nX):

		shelf[slope_yi[xi]-30:slope_yi[xi], xi] = 1

		shelf[slope_yi[xi]:slope_yi[xi]+80, xi] = 2

		shelfSI[:slope_yi[xi], xi] = 1

		shelfSI[slope_yi[xi]:, xi] = 2

	shelf = np.where(bathy>=0, 0, shelf)

	shelfSI = np.where(bathy>=0, 0, shelfSI)

	coast = tools.coastalMask(grid)

	

	if subr:

		curl = tools.getSubregionXY(curl, latsi, lonsi)

		SI = tools.getSubregionXY(SI, latsi, lonsi)

		FWflx = tools.getSubregionXY(FWflx, latsi, lonsi)

		coast = tools.getSubregionXY(coast, latsi, lonsi)

				

	# Convert into Ekman

	rho0 = 1028.5

	f = PAS_tools.getCoriolis(Y)

	wk = curl / (rho0 * f)



	shelf = shelf[...,iw:ie]

	shelfSI = shelfSI[...,iw:ie]

	coast = coast[...,iw:ie]

	curl = curl[...,iw:ie]

	wk = wk[...,iw:ie]

	SI = SI[...,iw:ie]

	FWflx = FWflx[...,iw:ie]

	

	pt.plot1by2([shelfSI, coast])

	np.save('shelf.npy', shelf)

	quit()

	

	# Make mask

	mask = np.zeros((curl.shape))

	maskSI = np.zeros((SI.shape))

	maskCoast = np.zeros((wk.shape))

	for ti in range(curl.shape[0]):

		mask[ti] = shelf

		maskSI[ti] = shelfSI

		maskCoast[ti] = coast

	

	#==

	

	curl_shelf = np.ma.array(curl.copy(), mask=mask!=1)

	curl_deep = np.ma.array(curl.copy(), mask=mask!=2)

	wk_shelf = np.ma.array(wk.copy(), mask=mask!=1)

	wk_deep = np.ma.array(wk.copy(), mask=mask!=2) 

	wk_coast = np.ma.array(wk.copy(), mask=maskCoast!=1)

	SI_shelf = np.ma.array(SI.copy(), mask=maskSI!=1)

	SI_deep = np.ma.array(SI.copy(), mask=maskSI!=2)

	SI_coast = np.ma.array(SI.copy(), mask=maskCoast!=1)

	FWflx_shelf = np.ma.array(FWflx.copy(), mask=maskSI!=1)

	FWflx_deep = np.ma.array(FWflx.copy(), mask=maskSI!=2)

	FWflx_coast = np.ma.array(FWflx.copy(), mask=maskCoast!=1)

	   	

	curl_shelf = np.ma.mean(curl_shelf, axis=(1,2))

	curl_deep = np.ma.mean(curl_deep, axis=(1,2))	

	wk_shelf = np.ma.mean(wk_shelf, axis=(1,2))

	wk_deep = np.ma.mean(wk_deep, axis=(1,2))

	wk_coast = np.ma.mean(wk_coast, axis=(1,2))

	SI_shelf = np.ma.mean(SI_shelf, axis=(1,2))

	SI_deep = np.ma.mean(SI_deep, axis=(1,2))

	SI_coast = np.ma.mean(SI_coast, axis=(1,2))

	FWflx_shelf = np.ma.mean(FWflx_shelf, axis=(1,2))

	FWflx_deep = np.ma.mean(FWflx_deep, axis=(1,2))

	FWflx_coast = np.ma.mean(FWflx_coast, axis=(1,2))

	

	np.save('windStressCurl', np.ma.filled(curl, fill_value=0))

	np.save('curl_shelf', np.ma.filled(curl_shelf, fill_value=0))

	np.save('curl_deep', np.ma.filled(curl_deep, fill_value=0))

	np.save('wk', np.ma.filled(wk, fill_value=0))

	np.save('wk_shelf', np.ma.filled(wk_shelf, fill_value=0))

	np.save('wk_deep', np.ma.filled(wk_deep, fill_value=0))

	np.save('wk_coast', np.ma.filled(wk_coast, fill_value=0))

	np.save('SI_shelf', np.ma.filled(SI_shelf, fill_value=0))

	np.save('SI_deep', np.ma.filled(SI_deep, fill_value=0))

	np.save('SI_coast', np.ma.filled(SI_coast, fill_value=0))

	np.save('FWflx_shelf', np.ma.filled(FWflx_shelf, fill_value=0))

	np.save('FWflx_deep', np.ma.filled(FWflx_deep, fill_value=0))

	np.save('FWflx_coast', np.ma.filled(FWflx_coast, fill_value=0))

	

	quit()

	

#==



# Save 2D files such as vertical Ekman, sea ice area, and surface freshwater fluxes. 

save2Dfiles = False

if save2Dfiles:

	

	SIfile = 'SIarea'; SIufile = 'SIuice'; SIvfile = 'SIvice'

	SIhefffile = 'SIheff'; SIhsnowfile = 'SIhsnow'

	FWflxfile = 'oceFWflx'; SHIfwFlxfile = 'SHIfwFlx'

	oceQnetfile = 'oceQnet'

	

	EXFuwindile = 'EXFuwind'; EXFvwindfile = 'EXFvwind'

	EXFpressfile = 'EXFpress'; EXFatempfile = 'EXFatemp'

	EXFprecifile = 'EXFpreci'; EXFrofffile = 'EXFroff'

	EXFswdnfile = 'EXFswdn '; EXFlwdnfile = 'EXFlwdn'; EXFaqhfile = 'EXFaqh'

	

	fnames = [EXFuwindile, EXFvwindfile, EXFpressfile, EXFatempfile, EXFprecifile, EXFrofffile, EXFswdnfile, EXFlwdnfile, EXFaqhfile]



	for fname in fnames:

		data = readVariable(fname, PASDIR, file_format='nc', meta=False, var2D=True)

		np.save(fname+'.npy', np.ma.filled(SHIfwFlx, fill_value=np.nan))

		

	quit()

	

#==



saveEkman = False

if saveEkman:



	ustressfile = 'oceTAUX'; vstressfile = 'oceTAUY'

	grid = Grid_PAS(PASDIR)

	Y = grid.YC

	dx = grid.DXG; dy = grid.DYG

	

	# Get Ekman

	us = readVariable(ustressfile, PASDIR, file_format='nc', meta=False, var2D=True)

	vs = readVariable(vstressfile, PASDIR, file_format='nc', meta=False, var2D=True)

	curl = tools.ddx(vs, dx) - tools.ddy(us, dy)

	rho0 = 1028.5

	f = PAS_tools.getCoriolis(Y)

	wk = curl / (rho0 * f)

	

	np.save('wk', np.ma.filled(wk, fill_value=np.nan))

	

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

	

	uvfile = 'slope_uv_max.npy'

	surfuvfile = 'surf_uv_av.npy'

	timefile = 'PAS_time.npy'

	slopexfile = 'slope_x.npy'

	slopeyfile = 'slope_y.npy'

	curldeepfile = 'curl_deep.npy'

	curlshelffile = 'curl_shelf.npy'

	

	slope_uv = np.load(path+uvfile)

	surf_uv = np.load(path+surfuvfile)

	time = np.load(path+timefile)

	slope_x = np.load(path+slopexfile)

	slope_y = np.load(path+slopeyfile)

	curl_deep = np.load(path+curldeepfile)

	curl_shelf = np.load(path+curlshelffile)

	

	#===

	

	year = PAS_tools.getDecimalTime(time, t_start=t_start, t_end=t_end)

	

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

	

	#plotSections = ['westGetz', 'westPITW', 'westPITE']

	plotSections = ['westPITW', 'westPITE']



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

		

		from TPI import TPI_unfiltered_PSL as IPO

		IPO_start = 14*12+1; IPO_end = None

		IPO = np.array(IPO[IPO_start:IPO_end])

		IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

		IPO /= 100*np.max(np.abs(IPO))

		

		plt.figure(figsize=(8,3))

		plt.subplot(121)

		plt.plot(year, uv_mean, label='Deep along-slope flow', color='r')

		plt.plot(year, surf_uv_mean, label='Surface along-slope flow', color='k')

		plt.plot(year, IPO, '--k', label='IPO')				

		#plt.plot(year_[ys_:ye_], uv_mean[ys_:ye_], label='Deep along-slope flow')

		#plt.plot(year_[ys_:ye_], surf_uv_mean[ys_:ye_], label='Surface along-slope flow')

		plt.title('Along-slope average current speed (m/s)')

		#plt.ylim(0, 0.15)

		plt.grid();	plt.legend()

	

		plt.subplot(122)

		plt.plot((window_lengths-1)/12, corrw);  plt.xlabel('running mean window size (years)')

		plt.title('corr(u,v) vs running mean window size'); plt.grid()

		

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



readAllSlopeFiles = False

if readAllSlopeFiles:

	

	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779

	#t_start = 0; t_end = 10

	

	#DEMEAN = False;	DESEASON = False; DETREND = False; NORM = False

	DEMEAN = False; DESEASON = True; DETREND = True; NORM = True

	

	nl2 = int(15*12); nl = 2*nl2 + 1

	lags = np.linspace(-nl2,nl2,nl)

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	uvfile = 'slope_uv_max.npy'; uvmaxjfile = 'slope_uv_maxj.npy'

	surfuvfile = 'surf_uv_av.npy'; surfuvmaxjfile = 'surf_uv_maxj.npy'

	uwfile = 'slope_uw_av.npy'; uwmaxjfile = 'slope_uw_maxj.npy'

	usfile = 'slope_us_av.npy'; usmaxjfile = 'slope_us_maxj.npy'

	timefile = 'PAS_time.npy'

	maxjfile = 'maxj.npy'

	slopexfile = 'slope_x.npy';	slopeyfile = 'slope_y.npy'	

	curldeepfile = 'curl_deep.npy';	curlshelffile = 'curl_shelf.npy'

	wkdeepfile = 'wk_deep.npy';	wkshelffile = 'wk_shelf.npy'; wkcoastfile = 'wk_coast.npy'

	SIdeepfile = 'SI_deep.npy';	SIshelffile = 'SI_shelf.npy'; SIcoastfile = 'SI_coast.npy'

	FWflxdeepfile = 'FWflx_deep.npy'; FWflxshelffile = 'FWflx_shelf.npy'; FWflxcoastfile = 'FWflx_coast.npy'

	

	slope_uv = np.load(path+uvfile); slope_uv_maxj = np.load(path+uvmaxjfile)

	surf_uv = np.load(path+surfuvfile); surf_uv_maxj = np.load(path+surfuvmaxjfile)

	uw = np.load(path+uwfile); uw_maxj = np.load(path+uwmaxjfile)

	us = np.load(path+usfile); us_maxj = np.load(path+usmaxjfile)

	time = np.load(path+timefile)

	maxj = np.load(path+maxjfile)

	slope_x = np.load(path+slopexfile)

	slope_y = np.load(path+slopeyfile)

	curl_deep = np.load(path+curldeepfile)

	curl_shelf = np.load(path+curlshelffile)

	wk_deep = np.load(path+wkdeepfile);	wk_shelf = np.load(path+wkshelffile); wk_coast = np.load(path+wkcoastfile)

	SI_deep = np.load(path+SIdeepfile);	SI_shelf = np.load(path+SIshelffile); SI_coast = np.load(path+SIcoastfile)

	FWflx_deep = np.load(path+FWflxdeepfile); FWflx_shelf = np.load(path+FWflxshelffile)

	FWflx_coast = np.load(path+FWflxcoastfile)



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

	curl_deep = curl_deep[start:end]

	curl_shelf = curl_shelf[start:end]

	wk_deep = wk_deep[start:end]

	wk_shelf = wk_shelf[start:end]

	wk_coast = wk_coast[start:end]

	SI_deep = SI_deep[start:end]

	SI_shelf = SI_shelf[start:end]

	SI_coast = SI_coast[start:end]

	FWflx_deep = FWflx_deep[start:end]

	FWflx_shelf = FWflx_shelf[start:end]

	FWflx_coast = FWflx_coast[start:end]

	

	print(slope_uv.shape)

	barocl = slope_uv - surf_uv	

	

	#===============

	#===============



	# Get zonal sums and averages evaluated over sections of interest.

	plotSections = ['westPITW', 'westPITE']

	#plotSections = ['westGetz']

		

	uv_sum = np.zeros(len(year))

	surf_uv_sum = np.zeros(len(year))

	uw_sum = np.zeros(len(year))

	us_sum = np.zeros(len(year))

	barocl_sum = np.zeros(len(year))

	nx = 0

	for section in plotSections:

		iw = sections[section][0]; ie = sections[section][1]

		uv_sum += np.sum(slope_uv[:,iw:ie], axis=1)

		surf_uv_sum += np.sum(surf_uv[:,iw:ie], axis=1)

		uw_sum += np.sum(uw[:,iw:ie], axis=1)

		us_sum += np.sum(us[:,iw:ie], axis=1)

		barocl_sum += np.sum(barocl[:,iw:ie], axis=1)

		nx += ie-iw

		

	uv_mean = uv_sum / nx

	surf_uv_mean = surf_uv_sum / nx

	uw_mean = uw_sum / nx

	us_mean = us_sum / nx

	barocl_mean = barocl_sum / nx

	

	#==

	

	# 135 

	

	# Demean, deseason, detrend.

	

	# Demean

	if DEMEAN:

		uv_mean = PAS_tools.demean(uv_mean)

		surf_uv_mean = PAS_tools.demean(surf_uv_mean)

		barocl_mean = PAS_tools.demean(barocl_mean)

		uw_mean = PAS_tools.demean(uw_mean)

		us_mean = PAS_tools.demean(us_mean)

		curl_deep = PAS_tools.demean(curl_deep)

		curl_shelf = PAS_tools.demean(curl_shelf)

		wk_deep = PAS_tools.demean(wk_deep)

		wk_shelf = PAS_tools.demean(wk_shelf)

		wk_coast = PAS_tools.demean(wk_coast)

		SI_deep = PAS_tools.demean(SI_deep)

		SI_shelf = PAS_tools.demean(SI_shelf)

		SI_coast = PAS_tools.demean(SI_coast)

		FWflx_deep = PAS_tools.demean(FWflx_deep)

		FWflx_shelf = PAS_tools.demean(FWflx_shelf)

		FWflx_coast = PAS_tools.demean(FWflx_coast)

			

	# Deseason

	if DESEASON:

		uv_mean = PAS_tools.deseason(uv_mean)

		surf_uv_mean = PAS_tools.deseason(surf_uv_mean)

		barocl_mean = PAS_tools.deseason(barocl_mean)

		uw_mean = PAS_tools.deseason(uw_mean)

		us_mean = PAS_tools.deseason(us_mean)

		curl_deep = PAS_tools.deseason(curl_deep)

		curl_shelf = PAS_tools.deseason(curl_shelf)

		wk_deep = PAS_tools.deseason(wk_deep)

		wk_shelf = PAS_tools.deseason(wk_shelf)

		wk_coast = PAS_tools.deseason(wk_coast)

		SI_deep = PAS_tools.deseason(SI_deep)

		SI_shelf = PAS_tools.deseason(SI_shelf)

		SI_coast = PAS_tools.deseason(SI_coast)

		FWflx_deep = PAS_tools.deseason(FWflx_deep)

		FWflx_shelf = PAS_tools.deseason(FWflx_shelf)

		FWflx_coast = PAS_tools.deseason(FWflx_coast)

		

	# Detrend

	if DETREND:

		uv_mean = PAS_tools.detrend(uv_mean, year)

		surf_uv_mean = PAS_tools.detrend(surf_uv_mean, year)

		barocl_mean = PAS_tools.detrend(barocl_mean, year)

		uw_mean = PAS_tools.detrend(uw_mean, year)

		us_mean = PAS_tools.detrend(us_mean, year)

		curl_deep = PAS_tools.detrend(curl_deep, year)

		curl_shelf = PAS_tools.detrend(curl_shelf, year)

		wk_deep = PAS_tools.detrend(wk_deep, year)

		wk_shelf = PAS_tools.detrend(wk_shelf, year)

		wk_coast = PAS_tools.detrend(wk_coast, year)

		SI_deep = PAS_tools.detrend(SI_deep, year)

		SI_shelf = PAS_tools.detrend(SI_shelf, year)

		SI_coast = PAS_tools.detrend(SI_coast, year)

		FWflx_deep = PAS_tools.detrend(FWflx_deep, year)

		FWflx_shelf = PAS_tools.detrend(FWflx_shelf, year)

		FWflx_coast = PAS_tools.detrend(FWflx_coast, year)

		

	#==



	# Get correlations after a range of running means applied to data.

	nw = 61

	window_lengths = np.linspace(1,nw,nw)

	corr_surfDeep, p_value_surfDeep = PAS_tools.crossCorrWindowAv(surf_uv_mean, uv_mean, window_lengths, nl, nl2)

	corr_deepWind, p_value_deepWind = PAS_tools.crossCorrWindowAv(uv_mean, uw_mean, window_lengths, nl, nl2)

	corr_surfWind, p_value_surfWind = PAS_tools.crossCorrWindowAv(surf_uv_mean, uw_mean, window_lengths, nl, nl2)

	corr_deepStress, p_value_deepStress = PAS_tools.crossCorrWindowAv(uv_mean, us_mean, window_lengths, nl, nl2)

	corr_surfStress, p_value_surfStress = PAS_tools.crossCorrWindowAv(surf_uv_mean, us_mean, window_lengths, nl, nl2)

	corr_deepDeepCurl, p_value_deepDeepCurl = PAS_tools.crossCorrWindowAv(uv_mean, curl_deep, window_lengths, nl, nl2)

	corr_deepShCurl, p_value_deepShCurl = PAS_tools.crossCorrWindowAv(uv_mean, curl_shelf, window_lengths, nl, nl2)

	corr_deepWkDeep, p_value_deepWkDeep = PAS_tools.crossCorrWindowAv(uv_mean, wk_deep, window_lengths, nl, nl2)

	corr_deepWkShelf, p_value_deepWkShelf = PAS_tools.crossCorrWindowAv(uv_mean, wk_shelf, window_lengths, nl, nl2)

	corr_stressWkShelf, p_value_stressWkShelf = PAS_tools.crossCorrWindowAv(us_mean, wk_shelf, window_lengths, nl, nl2)

	corr_stressWkDeep, p_value_stressWkDeep = PAS_tools.crossCorrWindowAv(us_mean, wk_deep, window_lengths, nl, nl2)

	corr_deepSIdeep, p_value_deepSIdeep = PAS_tools.crossCorrWindowAv(uv_mean, SI_deep, window_lengths, nl, nl2)

	corr_deepSIshelf, p_value_deepSIshelf = PAS_tools.crossCorrWindowAv(uv_mean, SI_shelf, window_lengths, nl, nl2)

	corr_surfSIdeep, p_value_surfSIdeep = PAS_tools.crossCorrWindowAv(surf_uv_mean, SI_deep, window_lengths, nl, nl2)

	corr_surfSIshelf, p_value_surfSIshelf = PAS_tools.crossCorrWindowAv(surf_uv_mean, SI_shelf, window_lengths, nl, nl2)

	

	#corrs = [corr_surfDeep, corr_surfWind, corr_deepWind, corr_deepStress]#, corr_surfStress]

	#pvals = [p_value_surfDeep, p_value_surfWind, p_value_deepWind, p_value_deepStress]#, p_value_surfStress]

	#titles = ['corr(surface current, undercurrent)', 'corr(surface current, surface wind speed)', 'corr(undercurrent, surface wind speed)', 'corr(undercurrent, surface stress)']#, 'corr(surface current, surface stress)']

	

	#corrs = [corr_deepDeepCurl, corr_deepShCurl, corr_deepWind, corr_deepStress]#, corr_surfStress]

	#pvals = [p_value_deepDeepCurl, p_value_deepShCurl, p_value_deepWind, p_value_deepStress]#, p_value_surfStress]

	#titles = ['corr(undercurrent, deep ocean wind stress curl)', 'corr(undercurrent, shelf wind stress curl)', 'corr(undercurrent, surface wind speed)', 'corr(undercurrent, surface stress)']#, 'corr(surface current, surface stress)']

	

	#corrs = [corr_deepDeepCurl, corr_deepShCurl, corr_deepWkDeep, corr_deepWkShelf]

	#pvals = [p_value_deepDeepCurl, p_value_deepShCurl, p_value_deepWkDeep, p_value_deepWkShelf]

	#titles = ['corr(undercurrent, deep ocean wind stress curl)', 'corr(undercurrent, shelf wind stress curl)', 'corr(undercurrent, wk deep ocean)', 'corr(undercurrent, wk shelf)']



	corrs = [corr_stressWkDeep, corr_stressWkShelf, corr_deepWkDeep, corr_deepWkShelf]

	pvals = [p_value_stressWkDeep, p_value_stressWkShelf, p_value_deepWkDeep, p_value_deepWkShelf]

	titles = ['corr(surface stress, wk deep ocean)', 'corr(surface stress, wk shelf)', 'corr(undercurrent, wk deep ocean)', 'corr(undercurrent, wk shelf)']



	#corrs = [corr_surfSIdeep, corr_deepSIdeep, corr_surfSIshelf, corr_deepSIshelf]

	#pvals = [p_value_surfSIdeep, p_value_deepSIdeep, p_value_surfSIshelf, p_value_deepSIshelf]

	#titles = ['corr(surface current, deep ocean SI area)', 'corr(undercurrent, deep ocean SI area)', 'corr(surface current, shelf SI area)', 'corr(undercurrent, shelf SI area)']

	

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

	#year = PAS_tools.windowAv(year, n=nn)[nn//2:-nn//2+1]

	year = year[nn//2:-nn//2+1]

	uv_mean = PAS_tools.windowAv(uv_mean, n=nn)[nn//2:-nn//2+1]

	surf_uv_mean = PAS_tools.windowAv(surf_uv_mean, n=nn)[nn//2:-nn//2+1]

	barocl_mean = PAS_tools.windowAv(barocl_mean, n=nn)[nn//2:-nn//2+1]

	uw_mean = PAS_tools.windowAv(uw_mean, n=nn)[nn//2:-nn//2+1]

	us_mean = PAS_tools.windowAv(us_mean, n=nn)[nn//2:-nn//2+1]

	curl_deep = PAS_tools.windowAv(curl_deep, n=nn)[nn//2:-nn//2+1]	

	curl_shelf = PAS_tools.windowAv(curl_shelf, n=nn)[nn//2:-nn//2+1]

	wk_deep = PAS_tools.windowAv(wk_deep, n=nn)[nn//2:-nn//2+1]	

	wk_shelf = PAS_tools.windowAv(wk_shelf, n=nn)[nn//2:-nn//2+1]

	wk_coast = PAS_tools.windowAv(wk_coast, n=nn)[nn//2:-nn//2+1]

	SI_deep = PAS_tools.windowAv(SI_deep, n=nn)[nn//2:-nn//2+1]	

	SI_shelf = PAS_tools.windowAv(SI_shelf, n=nn)[nn//2:-nn//2+1]

	SI_coast = PAS_tools.windowAv(SI_coast, n=nn)[nn//2:-nn//2+1]

	FWflx_deep = PAS_tools.windowAv(FWflx_deep, n=nn)[nn//2:-nn//2+1]	

	FWflx_shelf = PAS_tools.windowAv(FWflx_shelf, n=nn)[nn//2:-nn//2+1]

	FWflx_coast = PAS_tools.windowAv(FWflx_coast, n=nn)[nn//2:-nn//2+1]

	

	if NORM:

		uv_mean /= np.max(np.abs(uv_mean))

		surf_uv_mean /= np.max(np.abs(surf_uv_mean))

		barocl_mean /= np.max(np.abs(barocl_mean))

		uw_mean /= np.max(np.abs(uw_mean))

		us_mean /= np.max(np.abs(us_mean))

		curl_deep /= np.max(np.abs(curl_deep))

		curl_shelf /= np.max(np.abs(curl_shelf))

		wk_deep /= np.max(np.abs(wk_deep))

		wk_shelf /= np.max(np.abs(wk_shelf))

		wk_coast /= np.max(np.abs(wk_coast))

		SI_deep /= np.max(np.abs(SI_deep))

		SI_shelf /= np.max(np.abs(SI_shelf))

		SI_coast /= np.max(np.abs(SI_coast))

		FWflx_deep /= np.max(np.abs(FWflx_deep))

		FWflx_shelf /= np.max(np.abs(FWflx_shelf))

		FWflx_coast /= np.max(np.abs(FWflx_coast))

	

	plt.plot(year, uv_mean, label='Deep along-slope flow', color='r')

	plt.plot(year, surf_uv_mean, label='Surface along-slope flow', color='k')	

	#plt.plot(year, barocl_mean, label='Baroclinicity')			

	plt.plot(year, uw_mean, label='Along-slope wind')

	plt.plot(year, us_mean, label='Along-slope stress')

	#plt.plot(year, curl_deep, label='Deep ocean wind stress curl')

	#plt.plot(year, curl_shelf, label='Shelf wind stress curl')

	#plt.plot(year, wk_deep, label='Deep ocean Ekman')

	#plt.plot(year, wk_shelf, label='Shelf Ekman')

	#plt.plot(year, wk_coast, label='Coastal Ekman')		

	#plt.plot(year, SI_deep, label='Deep ocean SI area')

	#plt.plot(year, SI_shelf, label='Shelf SI area')

	#plt.plot(year, SI_coast, label='Coastal SI area')

	#plt.plot(year, FWflx_deep, label='Deep ocean FW flux')

	#plt.plot(year, FWflx_shelf, label='Shelf FW flux')

	#plt.plot(year, FWflx_coast, label='Coastal FW flux')	



	plt.title('Along-slope average current speed (m/s)')

	#plt.ylim(0, 0.15)

	plt.grid()

	plt.legend()

	plt.show()

		

	quit()

	

#==



seasonalBarocl = False

if seasonalBarocl:



	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779

	#t_start = 0; t_end = 10

	

	DEMEAN = False; DESEASON = False; DETREND = True; NORM = False

	#DEMEAN = True; DESEASON = True; DETREND = True; NORM = True

	start = 24*12 + 5*12 - 2; end=-1

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	uvfile = 'slope_uv_max.npy'

	surfuvfile = 'surf_uv_av.npy'

	timefile = 'PAS_time.npy'

	

	slope_uv = np.load(path+uvfile)

	surf_uv = np.load(path+surfuvfile)

	

	t = np.load(path+timefile)

	year = PAS_tools.getDecimalTime(t)

	year = year[start:end]



	slope_uv = slope_uv[start:end]

	surf_uv = surf_uv[start:end]	

	barocl = slope_uv - surf_uv

	

	#==

	

	# Get zonal sums and averages evaluated over sections of interest.

	#plotSections = ['westPITW', 'westPITE']

	plotSections = ['westPITE']

		

	uv_sum = np.zeros(len(year))

	surf_uv_sum = np.zeros(len(year))

	uw_sum = np.zeros(len(year))

	us_sum = np.zeros(len(year))

	barocl_sum = np.zeros(len(year))

	nx = 0

	for section in plotSections:

		iw = sections[section][0]; ie = sections[section][1]

		uv_sum += np.sum(slope_uv[:,iw:ie], axis=1)

		surf_uv_sum += np.sum(surf_uv[:,iw:ie], axis=1)

		barocl_sum += np.sum(barocl[:,iw:ie], axis=1)

		nx += ie-iw

	uv_mean = uv_sum / nx

	surf_uv_mean = surf_uv_sum / nx

	barocl_mean = barocl_sum / nx

	

	# Demean

	if DEMEAN:

		uv_mean = PAS_tools.demean(uv_mean)

		surf_uv_mean = PAS_tools.demean(surf_uv_mean)

		barocl_mean = PAS_tools.demean(barocl_mean)

			

	# Deseason

	if DESEASON:

		uv_mean = PAS_tools.deseason(uv_mean)

		surf_uv_mean = PAS_tools.deseason(surf_uv_mean)

		barocl_mean = PAS_tools.deseason(barocl_mean)

		

	# Detrend

	if DETREND:

		uv_mean = PAS_tools.detrend(uv_mean, year)

		surf_uv_mean = PAS_tools.detrend(surf_uv_mean, year)

		barocl_mean = PAS_tools.detrend(barocl_mean, year)

	

	baroclSeasonal = PAS_tools.seasonalData(barocl_mean, year)

		

	nn = 60

	barocl_mean = PAS_tools.windowAv(barocl_mean, n=nn)[nn//2:-nn//2+1]

	year_ = year[nn//2:-nn//2+1]

	

	for season in baroclSeasonal:

		tmp = tools.smooth3(baroclSeasonal[season])

		yr = tools.smooth3(np.array(year[::12]))

		plt.plot(yr, tmp, label=season)

	plt.plot(year_, barocl_mean)	



	plt.title('Along-slope baroclinicty (m/s)')

	plt.grid()

	plt.legend()

	plt.show()

		

	quit()



#==



seasonalData = False

if seasonalData:



	SEASONAL = True; ns = 4;

	SUBREGION = True; lonsi = None; latsi = None



	#SEASONAL = False; ns = 12 # If False, computes monthly instead.

	#SUBREGION = False

	

	

	grid = Grid_PAS(PASDIR)

	

	nn = 60	

	ts = -2; te = - 3

		

	from TPI import TPI_unfiltered_PSL as IPO

	IPO_start = 14*12+1 + ts; IPO_end = te

	IPO = np.array(IPO[IPO_start:IPO_end])

	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]



	#t_start = 107; t_end = 622

	t_start = 0; t_end = 779

	#t_start = 0; t_end = 10

	

	start = 24*12 + 5*12 + ts; end= -11 + te

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

		

	timefile = 'PAS_time.npy'

	t = np.load(path+timefile)

	year = PAS_tools.getDecimalTime(t)

	year = year[start:end]



	#fname = 'wk.npy'; title = 'Ekman (m/s)'; vmax1 = 2.e-6; vmax2 = vmax1

	#fname = 'FWflx.npy'; title = 'FWflx (kg/m^2/s)'; vmax1 = 2.e-4; vmax2 = 5.e-5; c0 = 0#4.e-6

	#fname = 'SST.npy'; title = 'SST (deg. C)'; vmax1 = .1; vmax2 = .1 

	#fname = 'SSS.npy'; title = 'SSS (g/kg)'; vmax1 = 1; vmax2 = .1; c0 = 33.8; vmax3 = 0.5

	#fname = 'SSH.npy'; title = 'SSH (m)'; vmax1 = .1; vmax2 = .1 

	#fname = 'SIv.npy'; title = 'SIv (m)'; vmax1 = .01; vmax2 = .01

	#fname = 'ua.npy'; title = 'Wind speed (m/s)'; vmax1 = 4.; vmax2 = 1.; c0 = 4.; vmax3 = 2. 

	#fname = 'SIheff.npy'; title = 'heff'; vmax1 = .4; vmax2 = .2; c0 = 1.5; vmax3 = 1.3 # vmax1=.4;vmac2=.2 for monthly

	#fname = 'SIhsnow.npy'; title = 'hsnow'; vmax1 = .1; vmax2 = .1; c0 = .25; vmax3 = .25

	#fname = 'EXFuwind.npy'; title = 'EXFuwind (m/s)'; vmax1 = 0.6; vmax2 = 0.6

	#fname = 'EXFvwind.npy'; title = 'EXFvwind (m/s)'; vmax1 = 0.2; vmax2 = 0.2

	#fname = 'SHIfwFlx.npy'; title = 'SHIfwFlx'

	#fname = 'SIconv'; vmax1 = 2.e-1; vmax2 = vmax1

	#fname = 'EXFpreci.npy'; title = 'EXFpreci (m/s)'; vmax1 = 1.e-8; vmax2 = vmax1

	#fname = 'EXFatemp.npy'; title = 'EXFatemp (degK)'; vmax1 = 8.; vmax2 = 2.5; c0 = 260. 

	#fname = 'EXFaqh.npy'; title = 'EXFaqh (kg/kg)'; vmax1 = 1.e-3; vmax2 = 4.e-4; c0 = 1.5e-3

	#fname = 'EXFpress.npy'; title = 'EXFpress (N/m^2)'; vmax1 = 4.e1; vmax2 = 1.5e2    

	#fname = 'EXFroff.npy'; title = 'EXFroff (m/s)'; vmax1 = 1.e-56; vmax2 = vmax1

	#fname = 'EXFlwdn.npy'; title = 'EXFlwdn (W/m^2)'; vmax1 = 4.e1; vmax2 = 1.e1; c0 = 230.  

	fname = 'EXFswdn.npy'; title = 'EXFswdn (W/m^2)'; vmax1 = 2.e2; vmax2 = 1.e1; c0 = .9; vmax3 = .9

	#fname = 'oceQnet.npy'; title = 'net surf. heat flux (W/m^2)'; vmax1 = 8.e1; vmax2 = 3.e1

	

	SAVE = False

	

	try:

		vmax3

	except:

		vmax3 = vmax1

	

	vmin1 = -vmax1; vmin2 = - vmax2; vmin3 = - vmax3



	data = np.load(path+fname)

	data = data[start:end]

	data = data[nn//2:-nn//2+1]

	year = year[nn//2:-nn//2+1]

	

	if SUBREGION:

		EASlons[0] = 225.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		try:

			data = tools.getSubregionXY(data, latsi, lonsi)

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

		except:

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	else:

		X = grid.XC; Y = grid.YC

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True)

	

	#data = PAS_tools.demean(data)

	#data = PAS_tools.detrendXY(data, None); c0 = 0

	# Use one below when plotting absolute values.

	#data = PAS_tools.detrendXY(data, None, interceptFlag=0)

		

	if SEASONAL:

		dataSeas, seasons = PAS_tools.seasonalDataIPO(data, year, IPO)

	else:

		dataSeas, seasons = PAS_tools.monthlyDataIPO(data, year, IPO)



	#s = 1.e5	

	#vmax1 = 1.e-4; vmin1 = -vmax1; vmax2 = 4.e-5; vmin2 = -vmax2

	#vmax1 = 2.e-4; vmin1 = -vmax1; vmax2 = 2.e-4; vmin2 = -vmax2

	

	vmin = [vmin1, vmin2, vmin2]

	vmax = [vmax1, vmax2, vmax2]

	

	vminc = [c0+vmin3, c0+vmin3, c0+vmin3]

	vmaxc = [c0+vmax3, c0+vmax3, c0+vmax3]

	

	X = [X, X, X]; Y = [Y, Y, Y]

	contour = [bathy, bathy, bathy]

	contourlevels=[[-1000]]*3

	savenum = ''

	

	ny = data.shape[-2]; nx = data.shape[-1]

	ipopos = np.zeros((ny,nx))

	

	for si in range(1):

		titles = [seasons[si] + ' ' + title, 'IPO pos comp', 'IPO neg comp']

		print(seasons[si])

		data1 = ptt.maskBathyXY(dataSeas[si,0], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data2 = ptt.maskBathyXY(dataSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data3 = ptt.maskBathyXY(dataSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data1 = ptt.maskDraftXY(data1, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data2 = ptt.maskDraftXY(data2, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data3 = ptt.maskDraftXY(data3, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		

		#data2 += data1; data3 += data1

		if not SEASONAL:

			savenum = str(si+1) 

		outname = 'tmpimg/' + fname[:-4] + savenum + seasons[si]

		

		# For plotting anomalies.

		#pt.plot1by3([data1, data2, data3], X=X, Y=Y, mesh=True, vmin=vmin, vmax=vmax, figsize=(13,3), titles=titles, fontsize=9, contour=contour, contourlevels=contourlevels, save=SAVE, outname=outname, show=SEASONAL)

	

		# For plotting absolute values.

		pt.plot1by3([data1, data2+data1, data3+data1], X=X, Y=Y, mesh=False, contourfNlevels=17, vmin=vminc, vmax=vmaxc, figsize=(13,3), titles=titles, fontsize=9, contour=contour, contourlevels=contourlevels, save=SAVE, outname=outname)

			

		# For checking colorbar limits.		

		#pt.plot1by3([data1, data2+data1, data3+data1], X=X, Y=Y, mesh=False, contourfNlevels=17, figsize=(13,3), titles=titles, fontsize=9, contour=contour, contourlevels=contourlevels, save=False, outname=outname)

		

	#pt.plotMbyN([[dataSeas[0][0], summer], [autumn, winter]], titles=titles, vmin=vmin, vmax=vmax)



	# Seasonailty & IPO

	# For each instance in time, check season

	# If IPO



	quit()



#==



seasonalQuiver = False

if seasonalQuiver:



	SUBREGION = False; lonsi = None; latsi = None

	#SUBREGION = True

	SAVE = True

		

	grid = Grid_PAS(PASDIR)

	

	nn = 60	

	ts = -2; te = - 3

		

	from TPI import TPI_unfiltered_PSL as IPO

	IPO_start = 14*12+1 + ts; IPO_end = te

	IPO = np.array(IPO[IPO_start:IPO_end])

	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]



	t_start = 0; t_end = 779

	

	start = 24*12 + 5*12 + ts; end= -11 + te

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

		

	#start += 99; end -= 201

	

	timefile = 'PAS_time.npy'

	t = np.load(path+timefile)

	year = PAS_tools.getDecimalTime(t)

	year = year[start:end]



	if False: # SI vel. and heff

		u = np.load(path+'SIu.npy'); v = np.load(path+'SIv.npy')

		data = np.load(path+'SIheff.npy'); vmax1 = 1.; vmax2 = 0.4

		title = 'SI vel. & heff'; qs = 0.05; scale = .3; qunits = 'm/s'

		savename = 'SIvel_heff'

	if True:

		u = np.load(path+'EXFuwind.npy'); v = np.load(path+'EXFvwind.npy')

		data = np.load(path+'wk.npy'); vmax1 = 1.e-5; vmax2 = 1.e-5

		title = 'wind & wk'; qs = .5; scale = 10; qunits = 'm/s'

		savename = 'wind_wk'	

	

	data = data[start:end];	data = data[nn//2:-nn//2+1]

	u = u[start:end]; u = u[nn//2:-nn//2+1]

	v = v[start:end]; v = v[nn//2:-nn//2+1]

	year = year[nn//2:-nn//2+1]

	print(u.shape)

	

	#EASlons[0] = 225.

	#EASlons[0] = 240.

	#latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

	#X, Y = grid.XYsubr(EASlons, EASlats)

	#data = tools.getSubregionXY(data, latsi, lonsi)

	#u = tools.getSubregionXY(u, latsi, lonsi)

	#v = tools.getSubregionXY(v, latsi, lonsi)

	#bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	#bathy = ptt.maskBathyXY(bathy, grid, 0, subregion=True, lons=lonsi, lats=latsi)

	

	if SUBREGION:

		EASlons[0] = 225.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

		try:

			data = tools.getSubregionXY(data, latsi, lonsi)

			u = tools.getSubregionXY(u, latsi, lonsi)

			v = tools.getSubregionXY(v, latsi, lonsi)

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

		except:

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	else:

		X = grid.XC; Y = grid.YC

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True)

		

	#data = PAS_tools.demean(data)

	data = PAS_tools.detrendXY(data, None)

	u = PAS_tools.detrendXY(u, None)

	v = PAS_tools.detrendXY(v, None)

			

	dataSeas, seasons = PAS_tools.seasonalDataIPO(data, year, IPO)

	uSeas, seasons = PAS_tools.seasonalDataIPO(u, year, IPO)

	vSeas, seasons = PAS_tools.seasonalDataIPO(v, year, IPO)

		

	del data, u, v

	

	vmin = [-vmax1, -vmax2, -vmax2]

	vmax = [vmax1, vmax2, vmax2]



	scale = [scale]*3

	qs = [qs]*3



	Ny2 = len(Y[:,0])//2; Nx2 = len(X[0,:])//2

	lat_0=Y[Ny2,0]; lon_0=X[0,Nx2]

	paras = [-74, -72, -70, -68, -66, -64]*3

	merids = [230, 240, 250, 260, 270]*3

	

	d = 16; Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	Xd = [Xd, Xd, Xd]; Yd = [Yd, Yd, Yd]

	X = [X, X, X]; Y = [Y, Y, Y]

	

	contour = [bathy, bathy, bathy]

	contourlevels=[[-1000]]*3

	

#	ny = data.shape[-2]; nx = data.shape[-1]

#	ipopos = np.zeros((ny,nx))

	for si in range(4):

		titles = [seasons[si] + ' ' + title, 'IPO pos comp', 'IPO neg comp']

		print(seasons[si])

		data1 = ptt.maskBathyXY(dataSeas[si,0], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data2 = ptt.maskBathyXY(dataSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data3 = ptt.maskBathyXY(dataSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data1 = ptt.maskDraftXY(data1, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data2 = ptt.maskDraftXY(data2, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		data3 = ptt.maskDraftXY(data3, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		

		u1 = ptt.maskBathyXY(uSeas[si,0], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		u2 = ptt.maskBathyXY(uSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		u3 = ptt.maskBathyXY(uSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		u1 = ptt.maskDraftXY(u1, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		u2 = ptt.maskDraftXY(u2, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		u3 = ptt.maskDraftXY(u3, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		

		v1 = ptt.maskBathyXY(vSeas[si,0], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		v2 = ptt.maskBathyXY(vSeas[si,1], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		v3 = ptt.maskBathyXY(vSeas[si,2], grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

		v1 = ptt.maskDraftXY(v1, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		v2 = ptt.maskDraftXY(v2, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		v3 = ptt.maskDraftXY(v3, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)[::d, ::d]

		

		#data2 += data1; data3 += data1

		outname = 'tmpimg/a' + savename + seasons[si]

		#pt.quiver1byN([u1,u2,u3], [v1,v2,v3], Xd, Yd, contourf=[data1,data2,data3], X=X, Y=Y, vmin=vmin, vmax=vmax, contour=contour, contourLevels=contourlevels, scale=scale, qs=qs, mesh=True, title=titles, fontsize=10, save=SAVE, outname=outname, figsize=(12,3))

		

		pt.quiver1byN_Basemap([u1,u2,u3], [v1,v2,v3], Xd, Yd, lat_0, lon_0, contourf=[data1,data2,data3], X=X, Y=Y, vmin=vmin, vmax=vmax, contour=contour, contourLevels=contourlevels, parallels=paras, meridians=merids, scale=scale, qs=qs, mesh=True, cbarShrink=0.5, title=titles, fontsize=7, save=SAVE, outname=outname, figsize=(12,3))

		

		#pt.plot1by3([data1, data2, data3], X=X, Y=Y, mesh=True, vmin=vmin, vmax=vmax, figsize=(13,3), titles=titles, fontsize=9, contour=contour, contourlevels=contourlevels, save=SAVE, outname=outname)

	

	#pt.plotMbyN([[dataSeas[0][0], summer], [autumn, winter]], titles=titles, vmin=vmin, vmax=vmax)



	# Seasonailty & IPO

	# For each instance in time, check season

	# If IPO



	quit()

	

#==



seaIce = False

if seaIce:



	# Load sea ice data. 

	# Figure out what correct metric for volume is.

	# Take running mean and plot.

	ts = -2; te = - 2

		

	nn = 60

	from TPI import TPI_unfiltered_PSL as IPO

	IPO_start = 14*12+1 + ts; IPO_end = te

	IPO = np.array(IPO[IPO_start:IPO_end])

	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

	

	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	grid = Grid_PAS(PASDIR)

	

	start = 24*12 + 5*12 + ts; end=-11 + te

	t = np.load(path+'PAS_time.npy')

	t = t[start:end]

	year = PAS_tools.getDecimalTime(t)[nn//2:-nn//2+1]

	print(len(year))

	

	# https://mitgcm.readthedocs.io/en/latest/phys_pkgs/seaice.html

	heff = np.load(path+'SIheff.npy') * grid.RAC; title = 'SI heff conv' # * grid.RAC

	print(heff.shape)

	heff = heff[start:end]

	

	heff = PAS_tools.deseason(heff)

	#heff = PAS_tools.detrendXY(heff, None)

	

	heff = PAS_tools.windowAv(heff, n=nn)[nn//2:-nn//2+1]

	heff = np.sum(heff, axis=(1,2))

	

	heff /= np.max(np.abs(heff))

	IPO /= np.max(np.abs(IPO))

	

	plt.plot(year, IPO, color='k', linestyle='--', label='IPO')

	plt.plot(year, heff, color='b', label='Sea-ice volume anomaly')

	#plt.plot(year, -heff, color='b', label='-Sea-ice volume')

	plt.grid()

	plt.legend()

	plt.show()

	quit()

	

#==



plotZ_ST = False

if plotZ_ST:



	from TPI import TPI_unfiltered_PSL as IPO

	IPO_start = 14*12+1; IPO_end = None

	IPO = np.array(IPO[IPO_start:IPO_end])

	

	grid = Grid_PAS(PASDIR)



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

		

	start = 24*12 + 5*12; end=-11

	

	fname = 'Pgrad.npy'

	

	t = np.load(path+'PAS_time.npy')

	Sgrad = np.load(path + fname)

	bathy = np.load(path + 'slopeBathy.npy')

	

	section = westPITE

		

	Sgrad = Sgrad[start:end]

	t = t[start:end]

		

	X = grid.Xsubr1D(grid.getIndexFromLon(EASlons))

	zs = [0, -600]; zz = grid.getIndexFromDepth(zs)

	Z = grid.RC.squeeze()[zz[0]:zz[1]]

	print(Z)

	dZ = grid.DRF.squeeze()[zz[0]:zz[1]]

	yslim = 50; ynlim = 20; nY = ynlim + yslim

	dy = 3.6e0

	Y = np.linspace(-yslim*dy, ynlim*dy, nY)

	

	#==



	Sgrad = np.ma.mean(Sgrad[...,section[0]:section[1]], axis=-1)

	

	nn = 60

	Sgrad = PAS_tools.windowAv(Sgrad[:,zz[0]:zz[1]], n=nn)[nn//2:-nn//2+1]

	

	Sgrad = PAS_tools.demean(Sgrad)

	#Sgrad = PAS_tools.desurf(Sgrad, axis=1)



	t = t[nn//2:-nn//2+1]

	year = PAS_tools.getDecimalTime(t)

	nT = len(t)

	text_data = ptt.getTextData(t, 'ctime', Y[1], Z[-1], PAS=True, color='w')



	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

	IPO = PAS_tools.detrend(IPO, tin=t)

	

	#==



	pt.plot1by1(Sgrad.T, X=year, Y=Z, title='Cross-slope pressure gradient', xlabel='Time', ylabel='depth')

	

	if fname == 'Sgrad.npy':

		for zi in range(1, len(dZ)):

			plt.plot(year, np.sum(Sgrad[:,:zi]*dZ[:zi], axis=1), label='depth=' + str(Z[zi]) + 'm')

			plt.legend()

			plt.show()

	else:

		Sgrad = -Sgrad

		Sgrad = PAS_tools.desurf(Sgrad, axis=1)

		#for i in range(12, len(dZ)):

		#	plt.plot(year, PAS_tools.detrend(Sgrad[:,i],None), label='level=' + str(i) + ', depth=' + str(Z[i]) + 'm')

		#	plt.legend(); plt.show()

		

		zs = [0,2,5,9,11,12,13,14,15,16]

		#zs = [9,11,12,13,14,15]

		tmp_max = 0

		for zi in zs:

			tmp = PAS_tools.detrend(Sgrad[:,zi],None); ylabel = 'PGF'

			#tmp = Z[zi]+1.e7*tmp; ylabel = '1.e7 * PGF + depth'

			corr = np.round(PAS_tools.pearson(tmp, IPO), 3)

			label = 'depth=' + str(Z[zi]) + 'm, corr='+str(corr)

			plt.plot(year, tmp, label=label)

			tmp_max = max(tmp_max, np.max(np.abs(tmp)))

		plt.title('Cross-slope pressure gradient')

		plt.ylabel(ylabel)

		plt.plot(year, tmp_max*IPO/np.max(np.abs(IPO)), 'k--', label='IPO')

		plt.grid(); plt.legend()

		plt.show()

	#plt.plot(year, Sgrad[:,9]); plt.plot(year, Sgrad[:,10]); plt.show()

	

	quit()

	

#==



plotZ_uv = False

if plotZ_uv:



	grid = Grid_PAS(PASDIR)



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

		

	#section = westPITW; sec = 'westPITW'

	#section = westPITE; sec = 'westPITE'

	section = westGetz; sec = 'westGetz'

	

	start = 24*12 + 5*12; end=-11

	

	hals_labels = ['34.0', '34.1', '34.2', '34.3', '34.4', '34.5', '34.6']

	hals = np.load(path+'depthS_'+sec+'.npy')

	

	t = np.load(path+'PAS_time.npy')

	uv = np.load(path + 'uv_z.npy')

	bathy = np.load(path + 'slopeBathy.npy')

	

	uv = uv[start:end]

	t = t[start:end]

		

	X = grid.Xsubr1D(grid.getIndexFromLon(EASlons))

	zs = [0, -500]; zz = grid.getIndexFromDepth(zs)

	Z = grid.RC.squeeze()[zz[0]:zz[1]]

	dZ = grid.DRF.squeeze()[zz[0]:zz[1]]

	yslim = 50; ynlim = 20; nY = ynlim + yslim

	dy = 3.6e0

	Y = np.linspace(-yslim*dy, ynlim*dy, nY)

	

	#==



	# Get Y,Z slice, mask, take running mean.	

	#xi = 10; S = S[...,xi]



	uv = np.ma.mean(uv[...,section[0]:section[1]], axis=-1)

	

	nn = 60

	uv = PAS_tools.windowAv(uv[:,zz[0]:zz[1]], n=nn)[nn//2:-nn//2+1]

	

	t = t[nn//2:-nn//2+1]

	year = PAS_tools.getDecimalTime(t)

	nT = len(t)

	text_data = ptt.getTextData(t, 'ctime', Y[1], Z[-1], PAS=True, color='w')



	nT, nZ = uv.shape

	uvb = np.zeros((nT, nZ))

	corr = np.zeros(nZ)

	for zi in range(nZ):

		uvb[:,zi] = uv[:,zi] - uv[:,0]

		corr[zi] = PAS_tools.pearson(uv[:,0], uv[:,zi])

	

	#==

	

	from matplotlib.cm import get_cmap

	cmap = get_cmap('jet')

	

	plt.plot(corr, Z, color='k'); plt.grid()

	plt.title('corr[u(Z), u(Z=surface)]; ' + sec)

	nh = len(hals)

	for hi, hal in enumerate(hals):

		plt.axhline(hal, label=hals_labels[hi], color=cmap(hi/nh))

	plt.ylabel('Z (m)')

	plt.legend()

	plt.show()

	quit()

		

	pt.plot1by1(uvb.T, X=year, Y=Z, contourfNlevels=15)



	pt.plot1by1(np.cumsum(uvb, axis=1).T, X=year, Y=Z, contourfNlevels=15)



	plt.plot(year, np.max(uvb, axis=1))

	plt.show()	



	quit()

	

#==



plotSlopeSalinity = False

if plotSlopeSalinity:

	

	grid = Grid_PAS(PASDIR)



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

		

	start = 24*12 + 5*12; end=-11

	

	t = np.load(path+'PAS_time.npy')

	S = np.load(path + 'slopeSalt.npy'); outname = 'slopeSalinity.mp4'; vmin = 33.8; vmax = 34.8;

	#S = np.load(path + 'slopeTheta.npy'); outname = 'slopeTheta.mp4'; vmin = -2.; vmax = 2

	bathy = np.load(path + 'slopeBathy.npy')

	

	section = westGetz

		

	S = S[start:end]

	t = t[start:end]

		

	X = grid.Xsubr1D(grid.getIndexFromLon(EASlons))

	zs = [0, -1000]; zz = grid.getIndexFromDepth(zs)

	Z = grid.RC.squeeze()[zz[0]:zz[1]]

	yslim = 50; ynlim = 20; nY = ynlim + yslim

	dy = 3.6e0

	Y = np.linspace(-yslim*dy, ynlim*dy, nY)

	

	#==



	# Get Y,Z slice, mask, take running mean.	

	#xi = 10; S = S[...,xi]

	print(S.shape)

	

	S = np.ma.mean(S[...,section[0]:section[1]], axis=-1)

	

	nn = 60

	S = PAS_tools.windowAv(S, n=nn)[nn//2:-nn//2+1]

	

	t = t[nn//2:-nn//2+1]

	nT = len(t)

	text_data = ptt.getTextData(t, 'ctime', Y[1], Z[-1], PAS=True, color='w')



	#==



	bathyMask = np.zeros((S.shape))

	for yi in range(nY):

		bathyMask[:, :, yi] = np.max(bathy[yi, section[0]:section[1]], axis=-1)

	for zi in range(len(Z)):

		bathyMask[:,zi,:] -= Z[zi]

	

	cmap = 'jet'; d = 1.e-4

	xlabel = 'Dist. from slope (km)'; ylabel = 'Depth (m)'

	title = 'S (60-month running mean)'

	#S = np.ma.masked_where(S<20., S)

	S = np.ma.masked_where(bathyMask>0, S)

	S = tools.boundData(S, vmin, vmax, d=d)

		

	pt.animate1by1(S, Y, Z, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, contourfNlevels=17, figsize=(4,4), outname=outname)

	

	quit()



#==



regressions = False

if regressions:



	# Which surface field?

	fname = 'wk.npy'; title = 'Ekman'

	#fname = 'FWflx.npy'; title = 'FWflx'

	#fname = 'SHIfwFlx.npy'; title = 'SHIfwFlx'

		

	# Which undercurrent section?

	#plotSections = ['westGetz', 'westPITW', 'westPITE']; title = title + '& full undercurrent', 

	#plotSections = ['westPITE']; title = title + ' & westPITE undercurrent'

	#plotSections = ['westPITW']; title = title + ' & westPITW undercurrent'

	#plotSections = ['westGetz']; title = title + ' & westGetz undercurrent'

	plotSections = ['westPITW', 'westPITE']; title = title + ' & undercurrent'

	

	#==

	

	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	draft = grid.draft

		

	data = np.load(path+fname)

	start = 24*12 + 5*12; end=-11

	data = data[start:end]	

	

	EASlons[0] = 230.

	latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	X = grid.Xsubr1D(EASlons)

	Y = grid.Ysubr1D(EASlats)

	print(bathy.shape)

	print(data.shape)

			

	try:

		data = tools.getSubregionXY(data, latsi, lonsi)

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	except:

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)



	#==

	

	# Get along-slope velocity timeseries

	uvfile = 'slope_uv_max.npy'

	surfuvfile = 'surf_uv_av.npy'

	

	slope_uv = np.load(path+uvfile)

	surf_uv = np.load(path+surfuvfile)

	

	uv = slope_uv - surf_uv	

	uv = slope_uv[start:end]



	uv_sum = np.zeros(uv.shape[0])	

	nx = 0

	for section in plotSections:

		iw = sections[section][0]; ie = sections[section][1]

		uv_sum += np.sum(uv[:,iw:ie], axis=1)

		nx += ie-iw		

	uv_mean = uv_sum / nx



	nt = len(uv_mean); t = np.linspace(0,nt,nt)

	uv_mean = PAS_tools.detrend(uv_mean, tin=t)

	

	nn = 60

	data = PAS_tools.windowAv(data, n=nn)[nn//2:-nn//2+1]

	uv_mean = PAS_tools.windowAv(uv_mean, n=nn)[nn//2:-nn//2+1]

	

	slope_y = np.load(path+'slope_y.npy')

	

	#==

	

	nt, ny, nx = data.shape	



	li = 6

	# Instantaneous correlation

	

	LAG = 0

	if LAG:

	

		nl2 = int(6); nl = 2*nl2 + 1

		lags = np.linspace(-nl2,nl2,nl)

		vmin2 = -nl2; vmax2 = nl2

	

		# Correlation for range of lags

		corr = np.zeros((nl, ny, nx))

		for j in range(0,ny//4):

			print(j)

			for i in range(nx):

				tmp = PAS_tools.detrend(data[:,j,i], tin=t)

				corr[:,j,i] = PAS_tools.crossCorr(data[:,j,i], uv_mean, nl, nl2, P_VALUE=False)

				 

		#for li in range(nl):

		#	print(li); print(lags[li])

		#	pt.plot1by1(corr[li])

		

		#corr[:nl2] = 0; vmin2 = 0; vmax2 = nl2

		

		corrmax = np.max(corr, axis=0)

		corrmin = np.min(corr, axis=0)

		data = np.where(corrmax>np.abs(corrmin), corrmax, corrmin)

				

		argmax = np.argmax(corr, axis=0)

		argmin = np.argmin(corr, axis=0)

		arg = lags[np.where(corrmax>np.abs(corrmin), argmax, argmin)]

		

		data = ptt.maskBathyXY(data, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		#data = ptt.maskDraftXY(data, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		arg = ptt.maskBathyXY(arg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		#arg = ptt.maskDraftXY(arg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		

		#vmax = np.max(np.abs(data)); vmin = -vmax

		vmax = 0.5

		vmin = [-vmax, vmin2-.5]; vmax = [vmax, vmax2+.5]

		

		titles = ['Largest corr: ' + title, 'Lag (months)']

		cmaps = ['bwr', 'jet']

		

		ymin = min(Y); ymax = max(Y)

		yscale = ymax - ymin

		vline1 = [X[iw], (slope_y[iw]-0-ymin)/yscale, (slope_y[iw]+1-ymin)/yscale]

		vline2 = [X[ie], (slope_y[ie]-0-ymin)/yscale, (slope_y[ie]+1-ymin)/yscale]

		vlines = [[vline1, vline2], [vline1, vline2]]

		

		pt.plot1by2([data, arg], X=[X,X], Y=[Y,Y], contour=[bathy,bathy], contourlevels=[[-1000],[-1000]], contourfNlevels=[11,14], titles=titles, cmaps=cmaps, vmin=vmin, vmax=vmax, vlines=vlines)

		

	#==

	

	else:

		corr = np.zeros((ny, nx))

		pval = np.zeros((ny, nx))

		for j in range(ny):

			for i in range(nx):

				#corr[j,i] = np.correlate(data[:,j,i], uv_mean)

				#corr[j,i], tmp = pearsonr(data[:,j,i], uv_mean)

				if data[0,j,i] != 0:

					#corr[j,i] = np.corrcoef(data[:,j,i], uv_mean)[0,1]

					corr[j,i], pval[j,i] = pearsonr(data[:,j,i], uv_mean)

				#corr[j,i] = np.correlate(data[li:,j,i], uv_mean[:-li])

		#		corr[j,i] = np.correlate(data[:-li,j,i], uv_mean[li:])

		

	

		corr = ptt.maskBathyXY(corr, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		corr = ptt.maskDraftXY(corr, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		pval = np.where(bathy>=0, 1, pval);	pval = np.where(draft<0, 1, pval)

				

		pt.plot1by1(corr, X=X, Y=Y, stippling=pval, stipData=[0.05, 4,4, .2], cmap='bwr', vmin=-1., vmax=1., contour=bathy, contourlevels=[-1000])



	quit()



#==



regressionsIPO = False

if regressionsIPO:



	# TPI_unfiltered_PSL, TPI_filtered_PSL

	from TPI import TPI_unfiltered_PSL as IPO



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	

	# Which surface field?

	#fname = 'wk.npy'; title = 'Ekman'

	#fname = 'FWflx.npy'; title = 'FWflx'

	#fname = 'SHIfwFlx.npy'; title = 'SHIfwFlx'

	fname = 'SIconv';

	

	#==

	

	if fname == 'SIconv':

		# https://mitgcm.readthedocs.io/en/latest/phys_pkgs/seaice.html

		heff = np.load(path+'SIheff.npy'); title = 'SI heff conv' # * grid.RAC

		#hsnow = np.load(path+'SIhsnow.npy'); title = 'SI hsnow conv'

		u = np.load(path+'SIu.npy')

		v = np.load(path+'SIv.npy')

		data = - tools.ddx(u*heff, grid.DXG) - tools.ddy(v*heff, grid.DYG)

		# POS IPO correlates with sea-ice convergence in region of interest.

		# Sea-ice conv. implies larger FW flux.

		# I've linked larger FW fluxes with low baroclinicity.

	else:

		data = np.load(path+fname)

	

	start = 24*12 + 5*12; end=-11

	data = data[start:end]	

	

	# IPO data is 1970-2017. 

	# Default start/end is Feb 1984 - Jan 2019

	IPO_start = 14*12+1; IPO_end = None

	IPO = np.array(IPO[IPO_start:IPO_end])

	

	t = np.load(path+'PAS_time.npy')

	t = t[start:end]

	year = PAS_tools.getDecimalTime(t)

	

	latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	X = grid.Xsubr1D(EASlons)

	Y = grid.Ysubr1D(EASlats)

			

	try:

		data = tools.getSubregionXY(data, latsi, lonsi)

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	except:

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)



	nn = 60

	data = PAS_tools.windowAv(data, n=nn)[nn//2:-nn//2+1]

	year = year[nn//2:-nn//2+1]

	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

	#IPO = IPO[nn//2:-nn//2+1]

	nt, ny, nx = data.shape	

	

	nl2 = int(0); nl = 2*nl2 + 1

	lags = np.linspace(-nl2,nl2,nl)

	vmin2 = -nl2; vmax2 = nl2

	

	nt = len(IPO); t = np.linspace(0,nt,nt)

	IPO = PAS_tools.detrend(IPO, tin=t)

	#plt.plot(year, IPO); plt.grid(); plt.show(); quit()

	

	if nl2 > 0:

		# Correlation for range of lags

		corr = np.zeros((nl, ny, nx))

		for j in range(0,ny):

			print(j)

			for i in range(nx):

				tmp = PAS_tools.detrend(data[:,j,i], tin=t)

				corr[:,j,i] = PAS_tools.crossCorr(tmp, IPO, nl, nl2, P_VALUE=False)

				 

		corrmax = np.max(corr, axis=0)

		corrmin = np.min(corr, axis=0)

		data = np.where(corrmax>np.abs(corrmin), corrmax, corrmin)

				

		argmax = np.argmax(corr, axis=0)

		argmin = np.argmin(corr, axis=0)

		arg = lags[np.where(corrmax>np.abs(corrmin), argmax, argmin)]

		

		data = ptt.maskBathyXY(data, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		#data = ptt.maskDraftXY(data, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		arg = ptt.maskBathyXY(arg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		#arg = ptt.maskDraftXY(arg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

		

		#vmax = np.max(np.abs(data)); vmin = -vmax

		vmax = 1.

		vmin = [-vmax, vmin2-.5]; vmax = [vmax, vmax2+.5]

		

		titles = ['Largest corr: ' + title, 'Lag (months)']

		cmaps = ['bwr', 'jet']

		

		#ymin = min(Y); ymax = max(Y)

		#yscale = ymax - ymin

		#vline1 = [X[iw], (slope_y[iw]-0-ymin)/yscale, (slope_y[iw]+1-ymin)/yscale]

		#vline2 = [X[ie], (slope_y[ie]-0-ymin)/yscale, (slope_y[ie]+1-ymin)/yscale]

		#vlines = [[vline1, vline2], [vline1, vline2]]

		

		pt.plot1by2([data, arg], X=[X,X], Y=[Y,Y], contour=[bathy,bathy], contourlevels=[[-1000],[-1000]], contourfNlevels=[11,14], titles=titles, cmaps=cmaps, vmin=vmin, vmax=vmax)

	

	else:

	

		corr = np.zeros((ny, nx))

		for j in range(0,ny):

			for i in range(nx):

				tmp = PAS_tools.detrend(data[:,j,i], tin=t)

				corr[j,i] = PAS_tools.pearson(tmp, IPO)

			

		title = 'corr[IPO, ' + title + ']'

		pt.plot1by1(corr, X=X, Y=Y, contour=bathy, contourlevels=[-1000], contourfNlevels=11, title=title)

	

	#==



#==



compositeIPO = False

if compositeIPO:



	# TPI_unfiltered_PSL, TPI_filtered_PSL

	from TPI import TPI_unfiltered_PSL as IPO



	#SUBREGION = False; lonsi = None; latsi = None

	SUBREGION = True

	

	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	

	# Which surface field?

	#fname = 'wk.npy'; title = 'Ekman (m/s)'; vmax1 = 2.e-6; vmax2 = vmax1

	#fname = 'FWflx.npy'; title = 'FWflx (kg/m^2/s)'; vmax1 = 1.e-5; vmax2 = 1.e-5#4.e-6

	#fname = 'SST.npy'; title = 'SST (deg. C)'; vmax1 = .1; vmax2 = .1 

	#fname = 'SSH.npy'; title = 'SSH (m)'; vmax1 = .1; vmax2 = .1 

	#fname = 'SIv.npy'; title = 'SIv (m)'; vmax1 = .01; vmax2 = .01 

	#fname = 'SIheff.npy'; title = 'heff'; vmax1 = .1; vmax2 = .1 

	#fname = 'EXFuwind.npy'; title = 'EXFuwind (m/s)'; vmax1 = 0.6; vmax2 = 0.6

	#fname = 'EXFvwind.npy'; title = 'EXFvwind (m/s)'; vmax1 = 0.2; vmax2 = 0.2

	#fname = 'SHIfwFlx.npy'; title = 'SHIfwFlx'

	#fname = 'SIconv'; vmax1 = 2.e-1; vmax2 = vmax1

	#fname = 'EXFpreci.npy'; title = 'EXFpreci (m/s)'; vmax1 = 8.e-10; vmax2 = vmax1

	#fname = 'EXFatemp.npy'; title = 'EXFatemp (degK)'; vmax1 = 3.e-1; vmax2 = vmax1  

	#fname = 'EXFaqh.npy'; title = 'EXFaqh (kg/kg)'; vmax1 = 1.e-4; vmax2 = vmax1  

	#fname = 'EXFlwdn.npy'; title = 'EXFlwdn (W/m^2)'; vmax1 = 4.e2; vmax2 = vmax1    

	#fname = 'EXFswdn.npy'; title = 'EXFswdn (W/m^2)'; vmax1 = 1.e3; vmax2 = vmax1

	#fname = 'oceQnet.npy'; title = 'oceQnet (W/m^2)'; vmax1 = 1.e2; vmax2 = vmax1

	#fname = 'EXFpress.npy'; title = 'EXFpress (N/m^2)'; vmax1 = 4.e1; vmax2 = 1.5e2    

	#fname = 'EXFroff.npy'; title = 'EXFroff (m/s)'; vmax1 = 1.e-56; vmax2 = vmax1



	#==

	

	if fname == 'SIconv':

		# https://mitgcm.readthedocs.io/en/latest/phys_pkgs/seaice.html

		SI = np.load(path+'SIheff.npy') * grid.RAC * np.load(path+'SI.npy'); title = 'SI heff conv' # 

		#SI = np.load(path+'SIhsnow.npy') * grid.RAC; title = 'SI hsnow conv'

		u = np.load(path+'SIu.npy')*SI

		v = np.load(path+'SIv.npy')*SI

		del SI

		data = - tools.ddx(u, grid.DXG) - tools.ddy(v, grid.DYG)

		# POS IPO correlates with sea-ice convergence in region of interest.

		# Sea-ice conv. implies larger FW flux.

		# I've linked larger FW fluxes with low baroclinicity.

	else:

		data = np.load(path+fname)

		

		

	#data = PAS_tools.demean(data)

	#pt.plot1by1(data[400], vmin=-vmax1, vmax=vmax1); quit()

	

	start = 24*12 + 5*12 - 2; end = -11 - 3

	data = data[start:end]	



	pt.plot1by1(np.mean(data,axis=0), vmin=-vmax1, vmax=vmax1, mesh=True); quit()



	# IPO data is 1970-2017. 

	# Default start/end is Feb 1984 - Jan 2019

	IPO_start = 14*12+1 - 2; IPO_end = -3#None

	IPO = np.array(IPO[IPO_start:IPO_end])

	

	t = np.load(path+'PAS_time.npy')

	t = t[start:end]

	year = PAS_tools.getDecimalTime(t)

	

	if SUBREGION:

		EASlons[0] = 230.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X = grid.Xsubr1D(EASlons)

		Y = grid.Ysubr1D(EASlats)		

		try:

			data = tools.getSubregionXY(data, latsi, lonsi)

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

		except:

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	else:

		X = grid.XC[0,:]; Y = grid.YC[:,0]

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True)



	nn = 60

	WINDOWAV = 1

	if WINDOWAV:

		data = PAS_tools.windowAv(data, n=nn)[nn//2:-nn//2+1]

		year = year[nn//2:-nn//2+1]

		IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

	else:

		data = data[nn//2:-nn//2+1]

		year = year[nn//2:-nn//2+1]

		IPO = IPO[nn//2:-nn//2+1]

		vmax1 *= .5e1; vmax2 *= .2e1	

	

	nt, ny, nx = data.shape

	data = PAS_tools.deseason(data)

	data = PAS_tools.demean(data)

	data = PAS_tools.detrendXY(data, None)

	

	IPO = PAS_tools.detrend(IPO, None)

	IPO = PAS_tools.demean(IPO)

	

	std = np.mean(IPO**2)**0.5



	plt.plot(year, IPO); plt.grid(); plt.show()



	compPos = np.zeros((ny, nx)); nPos = 0

	compNeg = np.zeros((ny, nx)); nNeg = 0

	for ti in range(nt):

		if IPO[ti] > std:

			compPos	+= data[ti]

			nPos += 1

		elif IPO[ti] < -std:

			compNeg += data[ti]

			nNeg +=1

			

	compPos /= nPos

	compNeg /= nNeg

	

	#==

	

	d = 1.e-9

	compPos = tools.boundData(compPos, -vmax1, vmax1, d=d)

	compNeg = tools.boundData(compNeg, -vmax2, vmax2, d=d)

	compPos = ptt.maskBathyXY(compPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	compNeg = ptt.maskBathyXY(compNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	compPos = ptt.maskDraftXY(compPos, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	compNeg = ptt.maskDraftXY(compNeg, grid, 0, timeDep=False, subregion=SUBREGION, lons=lonsi, lats=latsi)

	

	#from scipy.io import savemat

	#tmp = np.array(tools.smooth3(np.mean(compPos,axis=1))[23:-21])

	#savemat('EXFuwind_compPos.mat', {'data':tmp})

	#plt.contourf(X,Y, bathy); plt.plot(246+10*tmp, Y[23:-21]); plt.show()

	pt.plot1by1(compPos, vmin=-vmax1, vmax=vmax1)

	contour = [bathy, bathy]

	contourlevels = [[-1000], [-1000]]

	

	vmin = [-vmax1, -vmax2]; vmax = [vmax1, vmax2]

	titles = ['composite(IPO pos, '+title+')', 'composite(IPO neg, '+title+')']

	

	X = [X, X]; Y = [Y, Y]

	

	pt.plot1by2([compPos, compNeg], X=X, Y=Y, vmin=vmin, vmax=vmax, titles=titles, cmaps='bwr',  contourlevels=contourlevels, figsize=(9,3), contour=contour)

	

	quit()

	

	#==



#==



compositeIPOuv = False

if compositeIPOuv:



	# TPI_unfiltered_PSL, TPI_filtered_PSL

	from TPI import TPI_unfiltered_PSL as IPO



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	

	# Which surface field?

	#fname = 'wk.npy'; title1 = 'Ekman (m/s)'; vmax1 = 2.e-6; vmax2 = vmax1

	#fname = 'FWflx.npy'; title1 = 'FWflx (kg/m^2/s)'; vmax1 = 1.e-5; vmax2 = 1.e-5#4.e-6

	#fname = 'SHIfwFlx.npy'; title1 = 'SHIfwFlx'

	fname = 'SIconv'; vmax1 = .1; vmax2 = .4

	#fname = 'EXFpreci.npy'; title1 = 'EXFpreci (m/s)'; vmax1 = 8.e-10; vmax2 = vmax1

	#fname = 'EXFatemp.npy'; title1 = 'EXFatemp (degK)'; vmax1 = 3.e-1; vmax2 = vmax1  

	#fname = 'EXFlwdn.npy'; title1 = 'EXFlwdn (W/m^2)'; vmax1 = 2.e0; vmax2 = vmax1    

	#fname = 'EXFswdn.npy'; title1 = 'EXFswdn (W/m^2)'; vmax1 = 2.e0; vmax2 = vmax1

	#fname = 'EXFpress.npy'; title1 = 'EXFpress (N/m^2)'; vmax1 = 4.e1; vmax2 = 1.5e2    

	#fname = 'EXFroff.npy'; title1 = 'EXFroff (m/s)'; vmax1 = 1.e-56; vmax2 = vmax1

	

	#==

	

	if fname == 'SIconv':

		# https://mitgcm.readthedocs.io/en/latest/phys_pkgs/seaice.html

		#data = np.load(path+'SIheff.npy'); title1 = 'SI heff' # 

		#SI = np.load(path+'SIhsnow.npy') * grid.RAC; title1 = 'SI hsnow conv'

		u = np.load(path+'SIu.npy')

		v = np.load(path+'SIv.npy')

		umax = 6.; title2 = 'SI vel.'; qs = 0.05; scale = .5; qunits = 'm/s'

		data = - tools.ddx(u, grid.DXG) - tools.ddy(v, grid.DYG); title1 = 'SI vel conv'; vmax1=1.e-7;vmax2=vmax1

		#del SI

		# POS IPO correlates with sea-ice convergence in region of interest.

		# Sea-ice conv. implies larger FW flux.

		# I've linked larger FW fluxes with low baroclinicity.

	else:

		#u = np.load(path+'EXFuwind.npy'); v = np.load(path+'EXFvwind.npy'); umax = 6.; title2 = 'EXFwind'; qs = 0.2; scale = 4.; qunits = 'm/s'

		u = np.load(path+'SIu.npy'); v = np.load(path+'SIv.npy'); umax = 6.; title2 = 'SI vel.'; qs = 0.05; scale = .5; qunits = 'm/s'

		#u = np.load(path+'oceTAUX.npy'); v = np.load(path+'oceTAUY.npy'); umax = 1.e0; title2 = 'surf. stress'; qs = 1.e-2; scale = .1; qunits = 'N/m^2'

		data = np.load(path+fname)

		

	#data = PAS_tools.demean(data)

	#pt.plot1by1(data[400]); quit()

	

	start = 24*12 + 5*12; end = -11

	data = data[start:end]	

	u = u[start:end]

	v = v[start:end]

	

	# IPO data is 1970-2017. 

	# Default start/end is Feb 1984 - Jan 2019

	IPO_start = 14*12+1; IPO_end = None

	IPO = np.array(IPO[IPO_start:IPO_end])

	

	t = np.load(path+'PAS_time.npy')

	t = t[start:end]

	year = PAS_tools.getDecimalTime(t)

	

	EASlons[0] = 230.

	

	latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	#X1D = grid.Xsubr1D(EASlons)

	#Y1D = grid.Ysubr1D(EASlats)

	X, Y = grid.XYsubr(EASlons, EASlats)

			

	try:

		data = tools.getSubregionXY(data, latsi, lonsi)

		u = tools.getSubregionXY(u, latsi, lonsi)

		v = tools.getSubregionXY(v, latsi, lonsi)

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	except:

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	

	nn = 60

	WINDOWAV = 0

	if WINDOWAV:

		data = PAS_tools.windowAv(data, n=nn)[nn//2:-nn//2+1]

		u = PAS_tools.windowAv(u, n=nn)[nn//2:-nn//2+1]

		v = PAS_tools.windowAv(v, n=nn)[nn//2:-nn//2+1]

		year = year[nn//2:-nn//2+1]

		IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]

	else:

		data = data[nn//2:-nn//2+1]

		u = u[nn//2:-nn//2+1]

		v = v[nn//2:-nn//2+1]

		year = year[nn//2:-nn//2+1]

		IPO = IPO[nn//2:-nn//2+1]

		vmax1 *= 1.e0; vmax2 *= 1.e0	

		

	nt, ny, nx = data.shape

	

	data = PAS_tools.demean(data)

	data = PAS_tools.detrendXY(data, None)

	u = PAS_tools.demean(u)

	u = PAS_tools.detrendXY(u, None)

	v = PAS_tools.demean(v)

	v = PAS_tools.detrendXY(v, None)

	

	IPO = PAS_tools.detrend(IPO, None)

	IPO = PAS_tools.demean(IPO)

	

	std = np.mean(IPO**2)**0.5

	

	plt.plot(year, IPO); plt.grid(); plt.show()



	compPos = np.zeros((ny, nx)); nPos = 0

	compNeg = np.zeros((ny, nx)); nNeg = 0

	ucompPos = np.zeros((ny, nx)); ucompNeg = np.zeros((ny, nx))

	vcompPos = np.zeros((ny, nx)); vcompNeg = np.zeros((ny, nx))

	

	print(u.shape)

	print(v.shape)

	print(ucompPos.shape)

	for ti in range(nt):

		if IPO[ti] > std:

			compPos	+= data[ti]

			ucompPos += u[ti]

			vcompPos += v[ti]

			nPos += 1

		elif IPO[ti] < -std:

			compNeg += data[ti]

			ucompNeg += u[ti]

			vcompNeg += v[ti]

			nNeg +=1

			

	compPos /= nPos; compNeg /= nNeg

	ucompPos /= nPos; ucompNeg /= nNeg

	vcompPos /= nPos; vcompNeg /= nNeg

	

	pt.plot1by2([compPos,compNeg])

	

	#==

	

	d = 1.e-9

	compPos = tools.boundData(compPos, -vmax1, vmax1, d=d)

	compNeg = tools.boundData(compNeg, -vmax2, vmax2, d=d)

	compPos = ptt.maskBathyXY(compPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	compNeg = ptt.maskBathyXY(compNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	compPos = ptt.maskDraftXY(compPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	compNeg = ptt.maskDraftXY(compNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	ucompPos = tools.boundData(ucompPos, -umax, umax, d=d)

	ucompNeg = tools.boundData(ucompNeg, -umax, umax, d=d)

	ucompPos = ptt.maskBathyXY(ucompPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	ucompNeg = ptt.maskBathyXY(ucompNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	ucompPos = ptt.maskDraftXY(ucompPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	ucompNeg = ptt.maskDraftXY(ucompNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	vcompPos = tools.boundData(vcompPos, -umax, umax, d=d)

	vcompNeg = tools.boundData(vcompNeg, -umax, umax, d=d)

	vcompPos = ptt.maskBathyXY(vcompPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vcompNeg = ptt.maskBathyXY(vcompNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vcompPos = ptt.maskDraftXY(vcompPos, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	vcompNeg = ptt.maskDraftXY(vcompNeg, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	#==

	

	contour = [bathy, bathy]

	contourLevels = [[-1000], [-1000]]

	contourf = [compPos, compNeg]

	XX = [X, X]; YY = [Y, Y]

	vmin = [-vmax1, -vmax2]; vmax = [vmax1, vmax2]

	

	vmin = [-vmax1, -vmax2]; vmax = [vmax1, vmax2]

	titles = ['composite(IPO pos, '+title1+'); ' + title2, 'composite(IPO neg, '+title1+'); ' + title2]

	

	d = 12

	ucompPos = ucompPos[::d, ::d]; ucompNeg = ucompNeg[::d, ::d]

	vcompPos = vcompPos[::d, ::d]; vcompNeg = vcompNeg[::d, ::d]

	Xd = X[::d, ::d]; Yd = Y[::d, ::d]

	qs = [qs, qs]

	

	u = [ucompPos, ucompNeg]; v = [vcompPos, vcompNeg]

	Xd = [Xd, Xd]; Yd = [Yd, Yd]

	

	pt.quiver1byN(u, v, Xd, Yd, contourf=contourf, X=XX, Y=YY, vmin=vmin, vmax=vmax, contourfNlevels=13, cmap='bwr', contour=contour, contourLevels=contourLevels, scale=scale, title=titles, fontsize=7, qs=qs, qunits=qunits)

	

	Ny2 = len(Y[:,0])//2; Nx2 = len(X[0,:])//2

	lat_0=Y[Ny2,0]; lon_0=X[0,Nx2]

	

	pt.quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, X=XX, Y=YY, contourf=contourf, vmin=vmin, vmax=vmax, contourfNlevels=13, scale=scale, contour=contour, contourLevels=contourLevels, title=titles, fontsize=7, qs=qs, qunits=qunits)

	

	quit()



#==



isotherm = False 

if isotherm:

	

	grid = Grid_PAS(PASDIR)



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	#fname = 'ThermZ_0.npy'; vmin = -1000; vmax=0; cmap = 'YlOrRd'; d = 0.1; ylabel = 'Isotherm depth'

	fname = 'HalZ_342_S.npy'; vmin = -500; vmax=0; cmap = 'jet'; d = 0.01; ylabel = 'Isohaline depth'

	

	start = 24*12 + 5*12; end=-11



	X = grid.XC[1,:]; Y = grid.YC[:,1]

	bathy = grid.bathy; draft = grid.draft

	

	data = np.load(path+fname)

	t = np.load(path+'PAS_time.npy')

	slope_x = np.load(path+'slope_x.npy'); slope_y = np.load(path+'slope_y.npy')

	slope_xi = np.load(path+'slope_xi.npy'); slope_yi = np.load(path+'slope_yi.npy')

	maxj = np.load(path+'maxj.npy')



	data = data[start:end]

	t = t[start:end]

	year = PAS_tools.getDecimalTime(t)



	#pt.plot1by1(data[-1], vmin=vmin, vmax=vmax, mesh=True); quit()



	lats = EASlats; lons = EASlons

	lats[1] = - 69.

	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	lonsi[1] = westPITE[1] + lonsi[0]

	

	data = tools.getSubregionXY(data, latsi, lonsi)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	X = X[lonsi[0]:lonsi[1]+1]; Y = Y[latsi[0]:latsi[1]+1]

	

	nX = len(X); nY = len(Y)

	slope_yi = slope_yi[:nX]

	

	data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	

	#data = np.where(draft<0, np.nan, data)

	data = np.where(data!=data, vmax-d, data) 

	data = tools.boundData(data, vmin+d, vmax-d)

	

	nn = 60

	year = year[nn//2:-nn//2+1]

	data1 = PAS_tools.movingAv(data.copy(), n=nn)#[nn//2:-nn//2+1]

	

	maskYrange = 20

	buffer_ = 2

	# yref is maxj or slope_yi

	yref = slope_yi

	maskN = np.zeros(data1.shape)

	maskS = np.zeros(data1.shape)

	for i in range(nX):

		maskN[:,1+yref[i]+buffer_:1+yref[i]+buffer_+maskYrange,i] = 1

		maskS[:,yref[i]-buffer_-maskYrange:yref[i]-buffer_,i] = 1	



	#==

	

	data1N = np.ma.mean(np.ma.masked_where(maskN!=1, data1.copy()), axis=(1,2))

	data1S = np.ma.mean(np.ma.masked_where(maskS!=1, data1.copy()), axis=(1,2))

	

	plt.plot(year, data1N, label='north')

	plt.plot(year, data1S, label='south')	

	plt.plot(year, data1N-data1S, label='north-south')	

	plt.title(fname)

	plt.ylabel(ylabel); plt.xlabel('Time')

	plt.legend(); plt.grid()

	plt.show()

	

	#==

	

	tmp = np.ma.masked_where(maskS[-1]==1, data1[-1].copy())

	tmp = np.ma.masked_where(maskN[-1]==1, tmp)

	tmp = ptt.maskBathyXY(tmp, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	plt.contourf(X, Y, tmp)

	#plt.colorbar()

	plt.contour(X, Y, bathy, levels=[-1000])

	plt.gca().patch.set_color('.25')

	#plt.plot(maxj[0])

	plt.scatter(X, Y[slope_yi], s=1)

	#plt.scatter(slope_xi, slope_yi, s=0.5)

	plt.title('Continental slope north/south masks')

	

	plt.show()



	quit()

	

#==





allIsotherms = False

if allIsotherms:

	

	grid = Grid_PAS(PASDIR)



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'



	#files = ['ThermZ_p05_S.npy', 'ThermZ_0.npy', 'ThermZ_m05.npy', 'HalZ_343_S.npy', 'HalZ_344_S.npy', 'HalZ_345_S.npy', 'HalZ_346_S.npy']

	#files = ['HalZ_340_S.npy', 'HalZ_341_S.npy', 'HalZ_342_S.npy', 'HalZ_343_S.npy', 'HalZ_344_S.npy', 'HalZ_345_S.npy', 'HalZ_346_S.npy']

	files = ['HalZ_341_S.npy', 'HalZ_343_S.npy', 'HalZ_346_S.npy']

	#files = ['HalZ_340_S.npy']

	

	#section = westPITW; sec = 'westPITW'

	section = westPITE; sec = 'westPITE'

	#section = westGetz; sec = 'westGetz'

	

	maskYrange = 10; buffer_ = 2 # Default vals

	#maskYrange = 2; buffer_ = 1

		

	start = 24*12 + 5*12; end=-11

	

	vmin = -1000; vmax=0; d = 1.e-4

	

	from matplotlib.cm import get_cmap

	cmap = get_cmap('jet')



	X = grid.XC[1,:]; Y = grid.YC[:,1]

	bathy = grid.bathy; draft = grid.draft

	

	t = np.load(path+'PAS_time.npy')

	slope_x = np.load(path+'slope_x.npy'); slope_y = np.load(path+'slope_y.npy')

	slope_xi = np.load(path+'slope_xi.npy'); slope_yi = np.load(path+'slope_yi.npy')

	maxj = np.load(path+'maxj.npy')

	

	t = t[start:end]

	year = PAS_tools.getDecimalTime(t)



	lats = EASlats; lons = EASlons

	lats[1] = - 69.

	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	lonsi[1] = westPITE[1] + lonsi[0]

	X = X[lonsi[0]:lonsi[1]+1]; Y = Y[latsi[0]:latsi[1]+1]

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	draft = tools.getSubregionXY(draft, latsi, lonsi)

	nX = len(X); nY = len(Y)

	slope_yi = slope_yi[:nX]

		

	# yref is maxj or slope_yi

	yref = slope_yi

	

	nn = 60

	year = year[nn//2:-nn//2+1]

	nT = len(year)

	

	#==

	

	maskN = None

	barocls = {}

	

	nf = len(files)

	depthS = np.zeros(nf)

	depthN = np.zeros(nf)

	mean_depth = np.zeros(nf)

	

	for fi, fname in enumerate(files):

	

		print(fname)

	

		data = np.load(path+fname)

		data = data[start:end]



		data = tools.getSubregionXY(data, latsi, lonsi)

		data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	

		#data = np.where(draft<0, np.nan, data)

		data = np.where(data!=data, vmax-d, data) 

		data = tools.boundData(data, vmin+d, vmax-d)

	

		#data1 = PAS_tools.movingAv(data.copy(), n=nn)#[nn//2:-nn//2+1]

		data1 = PAS_tools.windowAv(data, n=nn)[nn//2:-nn//2+1]

		

		xw = 0; xe = nX

		xw = section[0]; xe = section[1]

		if maskN is None:

			maskN = np.zeros(data1.shape)

			maskS = np.zeros(data1.shape)

			for i in range(xw, xe):

				maskN[:,1+yref[i]+buffer_:1+yref[i]+buffer_+maskYrange,i] = 1

				maskS[:,yref[i]-buffer_-maskYrange:yref[i]-buffer_,i] = 1

		

		data1N = np.ma.mean(np.ma.masked_where(maskN!=1, data1.copy()), axis=(1,2))

		data1S = np.ma.mean(np.ma.masked_where(maskS!=1, data1.copy()), axis=(1,2))

	

		barocls[fname] = data1N - data1S

		

		depthS[fi] = np.mean(data1S)

		depthN[fi] = np.mean(data1N)

		mean_depth[fi] = 0.5 * np.mean(data1N + data1S)



	#==



	np.save(path+'depthS_'+sec, depthS)

	np.save(path+'depthN_'+sec, depthN)

	np.save(path+'meanDepth_'+sec, mean_depth)

	

	tot = np.zeros(nT)

	nb = len(barocls.keys())

	for bi, b  in enumerate(barocls.keys()):

		tmp = barocls[b]-np.mean(barocls[b])

		tot += tmp

		#tmp += mean_depth[bi]

		plt.plot(year, tmp, label=b, color=cmap(bi/nb))

			

	#plt.plot(year, tot/len(barocls), label='mean', color='k')

	plt.xlabel('Time')

	plt.legend(); plt.grid()

	plt.title('Cross-slope isotherm drop - ' + sec)

	plt.show()



	quit()



#==



computeIsotherm = False

if computeIsotherm:



	VAR = 'THETA'

	THERM = -0.5

	

	grid = Grid_PAS(PASDIR)

	

	T = readVariable(VAR, PASDIR, file_format='nc', meta=False)

	

	ThermZ = tools.getIsothermHeight(T, THERM, grid)

	

	print(ThermZ.shape)

	np.save('ThermZ', ThermZ)

	

	quit()



#==



computeIsotherm2 = False

if computeIsotherm2:



        VAR = 'THETA'

        THERM = -0.5



        grid = Grid_PAS(PASDIR)

        Z = grid.RC.squeeze()

        ny = grid.Ny; nx = grid.Nx



        nt = 779



        ThermZ = np.zeros((nt, ny, nx))



        for ti in range(nt):

                print(ti)

                T = readVariable(VAR, PASDIR, file_format='nc', meta=False, tt=ti)

                ThermZ[ti] = tools.getIsothermDepth3(T, Z, THERM)



        np.save('ThermZ', ThermZ)

        quit()



#==



computeIsohaline = False

if computeIsohaline:



        VAR = 'SALT'

        SAL = 34.3



        grid = Grid_PAS(PASDIR)

        Z = grid.RC.squeeze()

        ny = grid.Ny; nx = grid.Nx



        nt = 779



        HalZ = np.zeros((nt, ny, nx))



        for ti in range(nt):

                print(ti)

                S = readVariable(VAR, PASDIR, file_format='nc', meta=False, tt=ti)

                HalZ[ti] = tools.getIsohalineDepth3(S, Z, SAL)



        np.save('HalZ', HalZ)

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



animateSI = False

if animateSI:



	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	X = grid.XC; Y = grid.YC



	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	

	SIu = np.load(path+'SIu.npy')

	SIv = np.load(path+'SIv.npy')

	t = np.load(path+'PAS_time.npy')



	start = 24*12 + 5*12; end=-11

	SIu = SIu[start:end]

	SIv = SIv[start:end]

	t = t[start:end]



	# Subregion and mask		

	lats = [-76, -69]; lons = [225, 260]#[230, 270]#

	latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)

	SIu = tools.getSubregionXY(SIu, latsi, lonsi)

	SIv = tools.getSubregionXY(SIv, latsi, lonsi)

	bathy = tools.getSubregionXY(bathy, latsi, lonsi)

	X = grid.Xsubr1D(lons);	Y = grid.Ysubr1D(lats)

	SIu = ptt.maskBathyXY(SIu, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	SIv = ptt.maskBathyXY(SIv, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	SIu = ptt.maskDraftXY(SIu, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	SIv = ptt.maskDraftXY(SIv, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

	#

		

	nn = 60

	SIu = PAS_tools.windowAv(SIu, n=nn)[nn//2:-nn//2+1]

	SIv = PAS_tools.windowAv(SIv, n=nn)[nn//2:-nn//2+1]

	

	SIu = PAS_tools.demean(SIu)

	SIv = PAS_tools.demean(SIv)

	

	t = t[nn//2:-nn//2+1]

	text_data = ptt.getTextData(t, 'ctime', X[1], Y[1], PAS=True, color='k')

	

	#==

	

	d = 8

	SIu = SIu[:,::d,::d]; SIv = SIv[:,::d,::d]

	Xd = X[::d]; Yd = Y[::d]

	

	bathy = tools.boundData(bathy, -1000, 0)

	bathy = ptt.maskBathyXY(bathy, grid, 0, timeDep=False, subregion=True, lons=lonsi, lats=latsi)

	

	pt.animate1by1quiver(SIu, SIv, Xd, Yd, contour=bathy, X=X, Y=Y, contourf=True, text_data=text_data, title='Sea ice vel.', figsize=(7,6), cmap='jet')

		

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

	

#==



# Radiative snow-covered sea-ice budget for full PAS domain.

SIrad_all = True

if SIrad_all:



	SEASONAL = True; ns = 4;

	SUBREGION = True; lonsi = None; latsi = None



	#SEASONAL = False; ns = 12 # If False, computes monthly instead.

	#SUBREGION = False

	

	si = 0 # Season of choice. 0 for winter with default params.

	

	nn = 60	

	ts = -2; te = - 3

		

	from TPI import TPI_unfiltered_PSL as IPO

	IPO_start = 14*12+1 + ts; IPO_end = te

	IPO = np.array(IPO[IPO_start:IPO_end])

	IPO = PAS_tools.windowAv(IPO, n=nn)[nn//2:-nn//2+1]



	t_start = 0; t_end = 779

	start = 24*12 + 5*12 + ts; end= -11 + te

	

	path = '/home/michael/Documents/data/slopeCurrent/'

	path = path + str(t_start) + '_' + str(t_end) + '_y1/'

	

	grid = Grid_PAS(PASDIR)

	bathy = grid.bathy

	

	if SUBREGION:

		EASlons[0] = 225.

		latsi = grid.getIndexFromLat(EASlats); lonsi = grid.getIndexFromLon(EASlons)

		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X, Y = grid.XYsubr(EASlons, EASlats)

	else:

		X = grid.XC; Y = grid.YC

		

	timefile = 'PAS_time.npy'

	t = np.load(path+timefile)

	year = PAS_tools.getDecimalTime(t)[nn//2+start:end-nn//2+1]



	#==

	# For data in fnames, load data, remove start/end years, get subregion, detrend, get winter mean/IPO variability.

	

	fnames = ['ua', 'SSS', 'SIheff', 'SIhsnow', 'EXFatemp', 'EXFaqh', 'EXFlwdn', 'EXFswdn', 'FWflx']	

	dataDict = {}



	for fname in fnames:

	

		print('Processing ' + fname)

		

		data = np.load(path+fname+'.npy')[nn//2+start:end-nn//2+1]

		

		if SUBREGION:

			try:

				data = tools.getSubregionXY(data, latsi, lonsi)

				data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

			except:

				data = ptt.maskBathyXY(data, grid, 0, timeDep=True, subregion=True, lons=lonsi, lats=latsi)

		else:

			data = ptt.maskBathyXY(data, grid, 0, timeDep=True)

	

		#data = PAS_tools.detrendXY(data, None); c0 = 0

		# Use one below when plotting absolute values.

		#data = PAS_tools.detrendXY(data, None, interceptFlag=0)

		

		# Get 

		dataSeas, seasons = PAS_tools.seasonalDataIPO(data, year, IPO, )

		data = dataSeas[si]; #print(seasons[si])

		

		# Add season mean and absolute values of data during pos/neg IPO.

		dataDict[fname] = (data[0], data[1]+data[0], data[2]+data[0])

		

		#pt.plot1by3([dataDict[fname][0], dataDict[fname][1], dataDict[fname][2]], titles=[fname]*3)

		

	# End loop

	#==

	

	# Need to compute

	# Sens, blackbody radiation, latent, LW, SW, vertical flux

	

	# Compute surface radation balance for average winter. Involves solving for TiS.

	# Recompute (solve for TiS) given new LS, SW, qsat...

	# What info can we tell from plots of just radiative fluxes.

	

	IPO_index = 0 # 0 for season mean, 1 for positive IPO mean, 2 for negative IPO mean.

	

	# Reload all data from dataDict for season mean in given IPO phase.

	qa = dataDict['EXFaqh'][IPO_index]; Ua = dataDict['ua'][IPO_index]; Ta = dataDict['EXFatemp'][IPO_index]

	LW = dataDict['EXFlwdn'][IPO_index]; SW = dataDict['EXFswdn'][IPO_index]

	hi = dataDict['SIheff'][IPO_index]; hs = dataDict['SIhsnow'][IPO_index]; SSS = dataDict['SSS'][IPO_index]

	

	#==

	 

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



	Tf = SEAICE_dTempFrz_dS* SSS + SEAICE_tempFrz0 + celsius2K

	TiS = np.linspace(250, 280, 100)

	

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

	

	LW = - eps_s * LW

	

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

     

	SW =  (1.0 - ALB) * (1.0 - penetSWFrac) * SW    

	

	#==

	

	# CONDUCTIVE HEAT FLUX

	

	ki = 2.16560

	ks = 3.10000e-1

	

	FC = ki * ks * (Tf - TiS) / (ki * hs + ks * hi)

	

	#==

	

	poly = LATENT + SENS + LW + BB + SW - FC

	

	plt.plot(TiS, poly, label='poly')

	plt.plot(TiS, FC, label='FC')

	plt.plot(TiS, SENS, label='SENS')

	plt.plot(TiS, LW*np.ones(len(TiS)), label='LW')

	plt.plot(TiS, SW*np.ones(len(TiS)), label='SW')

	plt.plot(TiS, BB, label='BB')

	plt.plot(TiS, LATENT, label='LATENT')

	plt.grid()

	plt.legend()

	plt.show()

	

	quit()





#==



# Radiative snow-covered sea-ice budget.

SIrad = False

if SIrad:



	# Need to compute

	# Sens, blackbody radiation, latent, LW, SW, vertical flux

	

	# Compute surface radation balance for average winter. Involves solving for TiS.

	# Recompute (solve for TiS) given new LS, SW, qsat...

	# What info can we tell from plots of just radiative fluxes.

	

	path = '/home/michael/Documents/data/slopeCurrent/0_779_y1/'

	

	#==

	 

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



	

	#qa = np.load(path+'EXFaqh.npy')

	#Ua = (np.load(path+'EXFuwind.npy')**2 + np.load(path+'EXFvwind.npy')**2)**0.5

	#LW = np.load(path+'EXFlwdn.npy')

	#SW = np.load(path+'EXFswdn.npy')

	#Ta = np.load(path+'EXFatemp.npy')

	

	# Pos IPO vals # Ts = 262.8

	posIPOvals = {'qa':0.00182, 'Ua':3.2, 'Ta':263.5, 'LW':247.5, 'SW':0, 'hi':0.525, 'hs':0.15, 'SSS':33.95}

	# Neg IPO vals # Ts = 259.1

	negIPOvals = {'qa':0.0015, 'Ua':2.5, 'Ta':260., 'LW':232., 'SW':0, 'hi':0.7, 'hs':0.21, 'SSS':33.925}

	

	vals = negIPOvals

	qa = vals['qa']; Ua = vals['Ua']; Ta = vals['Ta']; LW = vals['LW']

	SW = vals['SW']; hi = vals['hi']; hs = vals['hs']; SSS = vals['SSS']



	Tf = SEAICE_dTempFrz_dS* SSS + SEAICE_tempFrz0 + celsius2K

	TiS = np.linspace(250, 280, 100)

	

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

	

	LW = - eps_s * LW

	

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

     

	SW =  (1.0 - ALB) * (1.0 - penetSWFrac) * SW    

	

	#==

	

	# CONDUCTIVE HEAT FLUX

	

	ki = 2.16560

	ks = 3.10000e-1

	

	FC = ki * ks * (Tf - TiS) / (ki * hs + ks * hi)

	

	#==

	

	poly = LATENT + SENS + LW + BB + SW - FC

	

	plt.plot(TiS, poly, label='poly')

	plt.plot(TiS, FC, label='FC')

	plt.plot(TiS, SENS, label='SENS')

	plt.plot(TiS, LW*np.ones(len(TiS)), label='LW')

	plt.plot(TiS, SW*np.ones(len(TiS)), label='SW')

	plt.plot(TiS, BB, label='BB')

	plt.plot(TiS, LATENT, label='LATENT')

	plt.grid()

	plt.legend()

	plt.show()

	

	quit()



#==



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
