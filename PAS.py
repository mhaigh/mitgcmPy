import numpy as np



from grid_PAS import Grid as Grid_PAS



import matplotlib.pyplot as plt

import matplotlib.colors as cl

import matplotlib.animation as animation



from time import ctime



import plotting as pt

import plotting_tools as ptt



import PAS_tools

import tools



from readData import * 



#==========================================================



DATADIR = '/data/oceans_output/shelf/pahol/mitgcm/'

#DATADIR = '/home/michael/Documents/data/'

PASDIR = DATADIR + 'PAS_851/run/'

EASlats = [-75.5, -70.5]; EASlons = [235, 260] 

EASlats = [-75.5, -68]



MAIN = True

if MAIN:



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

	np.save('slope_x', slope_x); np.save('slope_y', slope_y)



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



	ufile = 'UVEL'; vfile = 'VVEL';	rhofile = 'RHOAnoma'

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



			#==

		

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

	#	plt.plot(slope_uv_av[ti], label=ti)

	#plt.legend()

	#plt.show()

			

	quit()

	# Get y-grid point closest to slope_y_nX

	# Read u, v at this point and n grid points north and south.

	# Get component of each in along-slope direction



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

		pathno = 1

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

	

	slice_lon = 252.5

	xi = grid.getIndexFromLon(slice_lon)

	zs = [0, -1000]

	zz = grid.getIndexFromDepth(zs)



	imgpath = '/home/michael/Documents/Python/mitgcmPy/tmpimg/'

	ufile = 'UVEL'; vfile = 'VVEL';	rhofile = 'RHOAnoma'

	for ti in range(1):

		

		for xi in range(nX):

			

			yy = [yshift+slope_yi[xi]-3, yshift+slope_yi[xi]+3]

			print(yy)

			rho = 1028.5 + readVariable(rhofile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			rho = ptt.maskBathyYZ(rho, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			rho = tools.boundData(rho, 1026, 1040)

			

			u = readVariable(ufile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			u = ptt.maskBathyYZ(u, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			v = readVariable(vfile, PASDIR, file_format='nc', meta=False, tt=ti, xx=xi+xshift, yy=yy, zz=zz)

			v = ptt.maskBathyYZ(v, grid, xi=xi+xshift, subregion=True, lats=yy, depths=zz, sameIasNC=True)

			

			slope_uv = u * np.sin(bearing_nX[xi]*np.pi/180.) + v * np.cos(bearing_nX[xi]*np.pi/180.)

			lim = 0.09;	slope_uv = tools.boundData(slope_uv, -lim, lim)

			slope_uv[0,0] = -lim; slope_uv[0,1] = lim

			

			pt.plot1by2([bathy, slope_uv], contour=[None, rho], X=[X,Y[yy[0]:yy[1],0]], Y=[Y,Z[zz[0]:zz[1]]], vlines=[[X[0,xi]],None], titles=['Bathymetry (m)', 'Along-slope flow speed (m/s)'], show=False, save=True, outname=imgpath+f"{xi:03}")

			

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



saveTime = True

if saveTime:



	rhofile = 'RHOAnoma'

	rho = readVariable(rhofile, PASDIR, file_format='nc', meta=True)

	time = np.array(rho['TIME'][:])

	np.save('PAS_time', time)

	print(time.shape)

	print(time)

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



latLonToCartesian2 = False

if latLonToCartesian2:



	grid = Grid_PAS(DATADIR + 'PAS_851/run/')

	bathy = grid.bathy

	lon = grid.XC; lat = grid.YC

	

	dA, x, y = PAS_tools.dA_from_latlon (lon, lat, periodic=False, return_edges=True)

	

	

	pt.plot1by2([x, y])

	

	

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
