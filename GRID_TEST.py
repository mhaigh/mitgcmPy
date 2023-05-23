import numpy as np

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import matplotlib.pyplot as plt
import matplotlib.colors as cl

import plotting as pt
import plotting_tools as ptt

import tools

from readData import * 
from varDict import getPlottingVars, getTrefSref

import time

#from io import readData

#==========================================================

# Different variables defined at different parts of the grid/stencil.

# build_land_mask returns 2D boolean array of land.
##	 hfac is fraction of vertical cell occupied by water, so if zero, then its land. Full column=0 -> coastline.
# To get ice shelf, get ocean mask, multiply by (surface) regions where hfac < 1.

# What to do today?
# input code. Can read either nc or mitgcm binaries, function takes filed name as input.
# check they give the same result.
# maybe also optionally return a dictonary with plotting properties such as vmax/vmin and cmap, title...
# animate theta.
# make sure grid looks right from hfacs.


# Load which nc file?
#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
#fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2

#==

# SIZE 12x10 tiled grid 480318175/
# SIZE 24x10 tiled grid 490234079

#==

surfaceArea = 0
if surfaceArea:

	path_root = '/home/michael/Documents/data/'
	run = 'PISOMIP_002'	
	path = path_root + run + '/run/'
	
	grid = Grid(path)
	bathy = grid.bathy

	dx = grid.DXG
	dy = grid.DYG
	hfac = grid.hFacC[0] * grid.hFacS[0] * grid.hFacW[0]
	
	pt.plot1by1(grid.hFacS[0])
	quit()
	
	y2 = 100
	print('Southern half surface area:')
	print(np.sum(dx[:y2,:]*dy[:y2,:]*hfac[:y2,:]))
	print('Northern half surface area:')
	print(np.sum(dx[y2:,:]*dy[y2:,:]*hfac[y2:,:]))
	
	quit()

#==

IBCSOnc = 0
if IBCSOnc:

	import matplotlib.cm as cmaps
	from mpl_toolkits.basemap import Basemap
	from netCDF4 import Dataset
	import copy
	
	my_cmap = copy.copy(plt.cm.get_cmap('Greys')) # get a copy of the gray color map
	my_cmap.set_bad(alpha=0) # set how the colormap handles 'bad' values

	
	#==
	
	# First get PAS bathymetry.
	
	PASpath = '/home/michael/Documents/data/IdealisedAmundsenSea/data/PAS_851/run/'
	#W = -130; E = -90; S = -76; N = -64;
	
			
	grid = Grid_PAS(PASpath)
	ice = grid.iceC
	land_ice = np.where(grid.bathy<0, 1, 0)
	#land_ice = np.ma.masked_where(land_ice==0, land_ice)
	land_ice = np.where(grid.iceC>0, 0.5, land_ice)	
	land_ice = np.where(land_ice==1, np.nan, land_ice)
	land_ice += 0.5
	
	bathy = grid.bathy
	
	alpha = np.where(grid.bathy>0, 1, 0)

	X = grid.XC - 360
	Y = grid.YC

	W = np.min(X[0,:]); E = np.max(X[0,:])
	S = np.min(Y[:,0]); N = np.max(Y[:,0])
	
	
	labelData = []
	figsize = (5,4)
	paras = [-74, -72, -70, -68, -66, -64]
	merids = [230, 240, 250, 260, 270] 
	tid_types = [11, 17, 41, 45, 70, 40]
	cmap = 'gist_rainbow'

	SUBREGION = True
	if SUBREGION:
		figsize = (3,2.5)
		tid_types = [11, 17, 41, 45, 70]
		W = -115; E = -94.9; S = -75.5; N = -70.5
		xw = np.argmin(np.abs(X[0,:]-W)); xe = np.argmin(np.abs(X[0,:]-E))
		ys = np.argmin(np.abs(Y[:,0]-S)); yn = np.argmin(np.abs(Y[:,0]-N))
		X = X[ys:yn, xw:xe]; Y = Y[ys:yn, xw:xe]; land_ice = land_ice[ys:yn, xw:xe]
		bathy = bathy[ys:yn, xw:xe]; alpha = alpha[ys:yn, xw:xe]
		paras = [-74, -72, -70]
		merids = [245, 250, 255, 260, 265]
		# Labels for trough
		ts = 12
		p1 = {'x':(1.91e5,1.536e6), 't':'PITW', 'tx':(1.915e5,1.536e6), 'tc':'k', 'ts':ts}
		p2 = {'x':(9.e5, 8.e5), 't':'C', 'tx':(9.002e5, 8.e5), 'tc':'k', 'ts':ts}
		#p3 = {'x':(1.28e6, 1.52e6), 't':'PITE', 'tx':(.99e6, 1.52e6), 'tc':'k', 'ts':ts}
		#p4 = {'x':(1.507e6, 1.612e6), 't':'R', 'tx':(1.40e6, 1.620e6), 'tc':'k', 'ts':ts}
		p3 = {'x':(1.28e6, 1.52e6), 't':'PITE', 'tx':(.91e6, 1.52e6), 'tc':'k', 'ts':ts}
		p4 = {'x':(1.507e6, 1.612e6), 't':'R', 'tx':(1.375e6, 1.620e6), 'tc':'k', 'ts':ts}
		labelData = [p1, p2, p3, p4]
		cmap = 'hsv'
		#cmap = 'gist_rainbow'
		
	Xl = X[0,:]; Yl = Y[:,0]
	Ny2 = len(Yl)//2; Nx2 = len(Xl)//2
	lat_0=Yl[Ny2]; lon_0=Xl[Nx2]
	
	#==
	
	RID = False
	if RID:
		fname = '/home/michael/Documents/data/IBCSO_v2_RID_WGS84.nc'
		data =  Dataset(fname, 'r')	
		rid = data['rid'][:]
	else:
		fname = '/home/michael/Documents/data/IBCSO_v2_TID_WGS84.nc'
		data =  Dataset(fname, 'r')	
		rid = data['tid'][:]
	
	print(data.variables)

	#rid = data['rid']

	lon = data['lon'][:]
	lat = data['lat'][:]	
	xw = np.argmin(np.abs(lon-W)); xe = np.argmin(np.abs(lon-E))
	ys = np.argmin(np.abs(lat-S)); yn = np.argmin(np.abs(lat-N))
	
	lon = lon[xw:xe]
	lat = lat[ys:yn]
	rid = rid[ys:yn, xw:xe]
	
	#rid = np.ma.masked_where(rid!=70, rid)
	#plt.contourf(lon, lat, rid)#, cmap=cmap)#, levels=levels)
	#plt.colorbar()
	#plt.show(); quit()
	
	col = (np.random.random(), np.random.random(), np.random.random())

	if RID:
		#plt.plot(rid[100,:]);plt.show();quit()
		rid_min = 10800.; rid_max = np.max(rid)
		rid = np.ma.masked_where(rid<rid_min, rid)
		rid = np.ma.masked_where(rid!=11066, rid)
	
		print(rid_max)
		print(np.min(np.where(rid<10000,1.e10,rid)))
		levels = np.linspace(rid_min, rid_max, 50)
	else:
		#rid = np.ma.masked_where(rid!=71,rid)
		# 11 - multibeam
		# 17 - combination of direct measurement methods
		# 41 - Interpolated based on a computer algorithm
		# 45 - Predicted bathymetry based on flight-derived gravity data
		# 70 - Depth value is taken from a pre-generated grid that is based on mixed source data types, e.g. single beam, multibeam, interpolation etc.
		# 40 - Predicted bathymetry based on satellite-derived gravity data
		# tid_types = [11, 17, 41, 45, 70]
		
		# Make all values of interest negative
		for tid_type in tid_types:
			rid = np.where(rid==tid_type, -rid, rid)
		# Then mask all positive values.
		rid = np.ma.masked_where(rid>0, rid)
		# Then map each tid_type to 1,2,3,4...
		for ti, tid_type in enumerate(tid_types):
			rid = np.where(rid==-tid_type, ti+1, rid)
				
		levels = [0.5+ti for ti in range(len(tid_types)+1)]
		cmap = cmaps.get_cmap(cmap)
		colours = [cmap((levels[ti]+levels[ti+1]-1)/(2*len(tid_types))) for ti in range(len(tid_types))]

	#==
	

	
	#==
	
	# PLOT
	
	dpi = 300
	
	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)
	
	fig = plt.figure(figsize=figsize, dpi=dpi)
	plt.subplot(111)
	ax = plt.gca()
	
	X0, Y0 = m(X,Y)
	lon, lat = np.meshgrid(lon, lat)
	lon0, lat0 = m(lon,lat)
	
	d = 1

	for ci, colour in enumerate(colours):
		m.scatter(lon0[0,0],lat0[0,0],color=colour, s=10, label='Type ' + str(ci+1), marker=',')
	ax = plt.gca()
	#plt.legend();plt.show();quit()

	m.contourf(lon0[::d], lat0[::d], rid[::d, ::d], levels=levels, cmap=cmap)
	#plt.contourf(X, Y, land_ice, cmap=my_cmap, zorder=10)
	m.contourf(X0, Y0, land_ice, cmap=my_cmap, vmin=0, vmax=1)
	#plt.imshow(land_ice[::-1], extent=[W, E, S, N], zorder=10, cmap=my_cmap)
	#plt.gca().set_aspect(3)
	
	#plt.scatter(-102, -71.5, s=4, color='r')
	m.contour(X0, Y0, bathy, levels=[-1000, -600], colors='k', linestyles='solid', linewidths=0.3)
	m.drawparallels(paras,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	m.drawmeridians(merids,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	ax.legend(loc=3, fontsize='x-small')
	
	for li in labelData:
		plt.scatter(li['x'][0], li['x'][1], s=1, color='w')
		plt.annotate(li['t'], li['tx'], color=li['tc'], fontsize=li['ts'])
		
	plt.title('Bedmachine/IBCSO bathymetry data source', fontsize=8)
	plt.grid()
	plt.tight_layout()
	plt.savefig('Figure_2.jpg')
	plt.show()
	quit()
	
#==

# ['_MutableMapping__marker', '__abstractmethods__', '__class__', '__class_getitem__', '__contains__', '__delattr__', '__delitem__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__geo_interface__', '__getattribute__', '__getitem__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__reversed__', '__setattr__', '__setitem__', '__sizeof__', '__slots__', '__str__', '__subclasshook__', '__weakref__', '_abc_impl', '_data', '_delegate', '_delegated_properties', '_props', 'clear', 'coordinates', 'from_dict', 'geometries', 'get', 'items', 'keys', 'pop', 'popitem', 'setdefault', 'type', 'update', 'values']
gkpg = 0
if gkpg:

	import fiona		
	
	fname = '/home/michael/Documents/data/IBCSO_v2_coverage.gpkg'
	
	with fiona.open(fname, layer='IBCSO_coverage') as layer:
		print(layer.schema)
		for feature in layer:
			print(feature['geometry'].values())
			
	quit()
			
#==

animBin = 0
if animBin:
	
	#==
	
	path = '/home/michael/Documents/data/PISOMIP_001/run/'
	
	#VAR = 'ADJtaux'; mult_gencost = 1.e9; reverse = True; vmax = 1.e-3; vmin = -vmax
	#title=VAR + ' (deg. C/(N m$^{-2}$)), OBJ=on-shelf heat, days 90-120'
		
	VAR = 'stateUvel'; mult_gencost = 1.; reverse = False; vmax = .2; vmin = -vmax
	title = VAR
	
	#VAR = 'm_boxmean_theta.0000000000.data'; mult_gencost = 1; reverse=False
	#vmax = 0; vmin = -1.e-4; title = VAR
		
	nz = 50; ny = 200; nx = 240; dims = (nz, ny, nx)
	#nx = 240; ny = 200;	dims = (ny, nx)
	
	data = readAllnp(VAR, path, dims, reverse=reverse) / mult_gencost
	#data = readnp(path+VAR, dims, rec=-1)
	print(data.shape)
	nt = data.shape[0]
	
	
	#plt.plot(np.mean(np.abs(data),axis=(1,2))); plt.show(); quit()

	#==
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	data = ptt.maskBathyXY(data, grid, 0, timeDep=True)
			
	#pt.plot1by1(data[0], X=X, Y=Y)
	if reverse:
		time = 2*86400*np.linspace(nt,1,nt)
	else:
		time = 2*86400*np.linspace(0,nt,nt+1)
	text_data = ptt.getTextData(time, 'day', X[1,1], Y[1,1], color='k')
		
	contour = grid.bathy
	levels = [-600, -501, -401]
	
	#for i in range(data.shape[0]):
	#	pt.plot1by1(data[i])
			
	#plt.plot(np.sum(data, axis=(1,2))); plt.show(); quit()	
	
	pt.plot1by1(data[0,], X, Y, mesh=True, cmap='bwr', vmin=0.4*vmin, vmax=0.4*vmax, title='HC sens. to zonal wind, 60 day lag (deg. C/(N m$^{-2}$))', xlabel='Lon (km)', ylabel='Lat (km)', fontsize=10, contour=contour); quit()
	
	pt.animate1by1(data, X, Y, vmin=vmin, vmax=vmax, cmap='bwr', title=title, text_data=text_data, fontsize=9, contour=contour, contourLevels=levels)
	
	quit()

#==

binReadTest = 0
if binReadTest:

	path = '/home/michael/Documents/data/MCS/run/'

	# = 'adxx_bottomdrag.0000000000.data'
	VAR = 'stateTheta.0000000001.data'
	#VAR = 'm_boxmean_theta.0000000000.data'
	#VAR = 'adxx_theta.0000000000.data'
	#VAR = 'ADJtaux.0000002304.data'
	#VAR = 'state2D.0000008640.data'

	nx = 240; ny = 200
	dims = (nx, ny)
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.

	#plt.plot(Y[:,100], grid.bathy[:,100]); plt.grid(); plt.show(); quit()
	#pt.plot1by1(grid.bathy, mesh=True); quit()
	
	#nx = 600; ny = 384
	data = readnp(path+VAR, dims, rec=-1)# / 1.e9
	print(data.shape)

	#data = data[0]
	#data = np.mean(data, axis=0)
	
	#vmax = 1.e6; vmin = -vmax
	#vmax = 1; vmin = -vmax
	vmax = None; vmin = None
	 		
	#data[:, 1:][:, ::2] = data[:,::2]
	
	#data = np.mean(data,axis=0); data = np.where(data<0, data, 0)
	#pt.plot1by1(data[24], mesh=True)#, vmax=1, vmin=-1)
	#plt.plot(data[-1, :]); plt.show(); plt.title(VAR)
 	
	#for i in range(data.shape[0]):
	#	pt.plot1by1(data[i])#, X=X, Y=Y, vmin=vmin, vmax=vmax, mesh=True, title=VAR)

	quit()
	
#==

readWindBin = 0
if readWindBin:

	path = '/home/michael/Documents/data/MCSwinds/'
	
	fnames = ['taux_05cos.bin', 'taux_05cos_tr16.bin', 'taux_05cos_tr-16.bin']
	labels = ['Default wind', 'Southward shift', 'Northward shift']
	nx = 240; ny = 200
	dims = (ny, nx)
	
	pathBathy =  '/home/michael/Documents/data/MCS_308/run/'
	grid = Grid(pathBathy)
	contour = grid.bathy
	
	CURL = True
	
	y = np.linspace(0,ny*2.5,ny)
	x = np.linspace(0,nx*2.5,nx)
	
	for fi, fname in enumerate(fnames):
		data = readnp(path+fnames[fi], dims, rec=0, dtype='>f8')[...,0]

		plt.contourf(x, y, contour)

		if CURL:		
			plt.plot(300+1.6e9*tools.ddx(data, 2.5e3), y, label=labels[fi])
			plt.xlabel('Wind stress curl')
		else:
			plt.plot(1.e4*data+300, y, label=labels[fi])
			plt.xlabel('Wind stress')
			
	plt.colorbar()	
	plt.ylabel('Lat')
	plt.grid()
	plt.legend()
	plt.show()
	quit()
	
#==
	
# Read MITgcm binary output using rdmds.
mitgcmutilsTest = 0
if mitgcmutilsTest:

	from MITgcmutils import rdmds
	from MITgcmutils.mds import readmeta
	
	path = '/home/michael/Documents/data/MCStheta/run/'
	#VAR = 'adxx_tauu.effective.0000000000'
	#VAR = 'Eta.0000000576'
	#VAR = 'stateExf.0000000001'
	VAR = 'm_boxmean_theta.0000000000'
	#VAR = 'ADJtaux.0000008928'
	
	data = rdmds(path+VAR)
	
	for i in range(1):
		data = rdmds(path+VAR, returnmeta=False, rec=i)
		print(data.shape)
				
		pt.plot1by1(data)
	quit()
	
	
# A selection of binary-reading functions for testing.		
binReadTest = 0
if binReadTest:

	from scipy.io import FortranFile

	path = '/home/michael/Documents/data/MCS/run/'
	VAR = 'adxx_tauu.effective.0000000000.data'
	#VAR = 'm_boxmean_eta.0000000000.data'
	#VAR = 'state2D.0000000576.data'
	
	def readnp(inputfilename, nx, ny):

		dtype = '>f'
		count = 1*nx*ny
		size = np.dtype(dtype).itemsize
		f = open(inputfilename, 'rb')
		#f.seek(count*size)
		data = np.fromfile(f, dtype=dtype, count=-1)
		data = data.reshape((int(data.size/(ny*nx)),ny,nx))

		return data

	def readF(inputfilename):
			f = FortranFile(inputfilename, 'r')
			f = f.read_reals(dtype='float64')
			return f
	
	def readslice(inputfilename,nx,ny,timeslice):
		f = open(inputfilename,'rb')
		f.seek(8*timeslice*nx*ny)
		field = np.fromfile(f,dtype='float64',count=nx*ny)
		field = np.reshape(field,(ny,nx))
		f.close()
		return field

	#for i in range(3):
		#field = readslice(path+VAR, 240,200,i)
	#	field = readnp(path+VAR, 240, 200, i)
		#print(np.mean(field))
		#print(field)	
	#	pt.plot1by1(field, vmin=-1.e-7, vmax=1.e-7)

	data = readnp(path+VAR, 240, 200)
	
	print(data.shape)
	pt.plot1by1(data[0]); quit()
	
	for i in range(data.shape[0]):
		pt.plot1by1(data[i])
    	
	quit()


xmitgcmTest = 0
if xmitgcmTest:

	from xmitgcm import open_mdsdataset
	path = '/home/michael/Documents/data/MCS/run/'
	
	VAR = 'adxx_tauu.0000000000.data'
	VAR = 'm_boxmean_eta.0000000000'
	
	data = open_mdsdataset(path, prefix=['adxx_tauu'], read_grid=False, ignore_unknown_vars=True)
	
	print(data)
	data = np.array(data.to_array(dtype=float))
	print(data.shape)
	
	quit()

volTest = False
if volTest:

	path = '/home/michael/Documents/data/MCS_141/run/'
	grid = Grid(path)
	 
	xlims = None#[98750,198750]
	ylims = [0, 238750]
	zlims = [0, -500]
	
	vol = grid.volume(xlims=xlims, ylims=ylims, zlims=zlims)
		
	print(vol)
	
	quit()
		
#==

fbeta = False
if fbeta:
	Omega = 2 * np.pi / 86400.
	a = 6371.e3
	
	# Get f0, beta from lat
	if 0:
		lat = -74
		f0 = 2 * Omega * np.sin(lat * np.pi / 180)
		beta = 2 * Omega * np.cos(lat * np.pi / 180) / a
		print(f0)
		print(beta)
	
	# Get lat, beta from f0
	else:
		f0 = -1.4e-4
		lat = 180 * np.arcsin(f0 / (2*Omega)) / np.pi
		beta = 2 * Omega * np.cos(lat * np.pi / 180) / a
		print(f0)
		print(beta)
		print(lat)
	
	
	quit()

#==

# PHIHYD
baroclinicEddies = False
if baroclinicEddies:
	
	path = '/home/michael/Documents/data/MCS_002/run/'
	grid = Grid(path)
	X = grid.XC[1,:]/1000.		
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	ts = 16; level = 24; yi = -5

	VAR = 'THETA'
	theta = readVariable(VAR, path, file_format='nc', meta=False)
	VAR = 'VVEL'
	vvel = readVariable(VAR, path, file_format='nc', meta=False)
	
	# Get time data. This assumes both data have same output freq.
	text = 'month ' + str(ts+1) + '; Z = ' + str(Z[level]) + ' m'
	text_data = [{'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}, None]

	#vvel_ = vvel[ts, :, -5]
	#
	vvel_ = vvel[ts, level, ]; theta_ = theta[ts, level,]	
	pt.plot1by2([vvel_, theta_], X=[X,X], Y=[Y,Y], xlabels=['X (km)', 'X (km)'], ylabels=['Y (km)', ''], titles=['VVel (m/s)', 'Theta (deg. C)'], text_data=text_data); quit()

	
	vvel_ = vvel[ts, :, yi,]
	theta_ = theta[ts, :, yi] - np.tile(np.mean(theta[ts, :, yi], axis=1), [theta.shape[3],1]).T
	
	vmax = 0.04; vmin = - vmax; vvel_ = tools.boundData(vvel_, vmin, vmax, scale=0.99999)	
	vmax = 0.05; vmin = - vmax; theta_ = tools.boundData(theta_, vmin, vmax, scale=0.99999)

	pt.plot1by2([vvel_, theta_], X=[X,X], Y=[Z,Z], xlabels=['X (km)','X (km)'], ylabels=['Depth (m)', ''], titles=['VVel (m/s)', 'Theta eddy (deg. C)'])

	quit()

#==

TEST_thetaHeight = False
if TEST_thetaHeight:

	PAS = True
	THERM = -0.4

	if PAS:
		path = '/home/michael/Documents/data/PAS_851/run/'
		grid = Grid_PAS(path)
		T = np.load(path+'Tmean_PAS851.npy')

		X = grid.XC#[1,:]
		Y = grid.YC#[:,1]
		Z = grid.RC.squeeze()

		# Get z-indices of level with Theta closest to THERM.
		Tz = np.argmin(np.abs(T-THERM),axis=0)
		ThermZ = Z[Tz]
		ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=False)

		bathy = ptt.maskBathyXY(grid.bathy, grid, 0, timeDep=False)

		title = 'PAS851 mean - 0.4 deg. isotherm height'
 
	else:
		path = '/home/michael/Documents/data/MCS_038/run/'
		grid = Grid(path)
		T = readVariable('THETA', path, file_format='nc', meta=False)[-2:]

		X = grid.XC[1,:]/1000.
		Y = grid.YC[:,1]/1000.
		Z = grid.RC.squeeze()

		# Get z-indices of level with Theta closest to THERM.
		Tz = np.argmin(np.abs(T-THERM),axis=1)
		ThermZ = Z[Tz]
		ThermZ = ptt.maskBathyXY(ThermZ, grid, 0, timeDep=True)

	SUBR = True
	if SUBR:
		lats = [-75, -68]; lons = [230, 270]
		latsi = grid.getIndexFromLat(lats)
		lonsi = grid.getIndexFromLon(lons)

		ThermZ = tools.getSubregionXY(ThermZ, latsi, lonsi)
		bathy = tools.getSubregionXY(bathy, latsi, lonsi)

		X = grid.Xsubr1D(lons)
		Y = grid.Ysubr1D(lats)

	pt.plot1by1(bathy, X=X, Y=Y, mesh=True, vmin=-1000, vmax=0, title='PAS851 bathymetry')
	pt.plot1by1(ThermZ, X=X, Y=Y, mesh=True, vmin=-800, vmax=-100, title=title); quit()#

	#Tt = -1.8; Tb = 1.; print(0.5*(Tt+Tb))

#==

TEST_rho = False
if TEST_rho:

	#path = '/home/michael/Documents/data/MCS_038/run/'
	path = '/home/michael/Documents/data/IdealisedAmundsenSea_master/IdealisedAmundsenSea1/data/MCS_161/run/'
	grid = Grid(path)
	
	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	Rho1 = readVariable('RHOAnoma', path, file_format='nc', meta=False, tt=-1)
	T = readVariable('THETA', path, file_format='nc', meta=False, tt=-1)
	S = readVariable('SALT', path, file_format='nc', meta=False, tt=-1)

	rho0 = 1030.
	tAlpha=3.90e-5
	sBeta =7.41e-4

	Tref, Sref = getTrefSref()
	Rho2 = tools.EOSlinear(rho0, sBeta, tAlpha, S, Sref, T, Tref)

	vmin = -0.2; vmax = -vmin
	pt.plot1by2([Rho1[...,-1], Rho2[...,-1]], X=[Y,Y], Y=[Z,Z], vmin=[vmin,vmin], vmax=[vmax,vmax])
	
	quit()

#==

# Plot T/S at northern boundary to ensure rel. profile correct.
TEST_rel = False
if TEST_rel:

	path = '/home/michael/Documents/data/MCS_116/run/'
	
	grid = Grid(path)
	
	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	T = readVariable('THETA', path, file_format='nc', meta=False)[-1]
	S = readVariable('SALT', path, file_format='nc', meta=False)[-1]

	Tref, Sref = getTrefSref()
	
	T = np.mean(T[...,-2,:], axis=-1)
	S = np.mean(S[...,-2,:], axis=-1)

	plt.subplot(121)
	plt.plot(T, Z); plt.grid()
	plt.subplot(122)
	plt.plot(S, Z); plt.grid()
	plt.show()

#==

TEST_brclnc = False
if TEST_brclnc:
	
	#path = '/home/michael/Documents/data/PAS_666/run/'
	path = '/home/michael/Documents/data/MCS_038/run/'
	
	level = 24

	grid = Grid(path)
	#grid = Grid_PAS(path)	

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	dY = 1.e3 * (Y[1] - Y[0])

	u = readVariable('UVEL', path, file_format='nc', meta=False)[-10:]
	v = readVariable('VVEL', path, file_format='nc', meta=False)[-10:]

	ub = tools.bottom(u, grid, 'u')
	vb = tools.bottom(v, grid, 'v')

	us = u[:,0] - ub
	vs = v[:,0] - vb

	vmin = -0.1; vmax = 0.1; vmin = [vmin, vmin]; vmax = [vmax, vmax]
	pt.plot1by2([us[-1],vs[-1]], vmin=vmin, vmax=vmax)

	contour = grid.bathy
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

	# Sample rate & level
	d = 6

	title = '(u_top-u_bot, v_top-v_bot) & bathymetry'

	us = us[..., ::d, ::d]; vs = vs[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]
	
	#pt.plot1by2([data1[0], data2[0]]); quit()

	pt.animate1by1quiver(us, vs, Xd, Yd, contour=contour, X=X, Y=Y, contourf=False, title=title)
	quit()
	

#==

TEST_baroStr = False
if TEST_baroStr:

	#path = '/home/michael/Documents/data/PAS_666/run/'
	path = '/home/michael/Documents/data/MCS_038/run/'
	
	level = 24

	grid = Grid(path)
	#grid = Grid_PAS(path)	

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	dY = 1.e3 * (Y[1] - Y[0])

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	data = readVariable(VAR, path, file_format='nc', meta=True)
	u = data[VAR][-10:]

	SSH = readVariable('ETAN', path, file_format='nc', meta=False)[-10:]

	for ti in range(u.shape[0]):
		u[ti,] = ptt.maskBathyAll(u[ti,], grid)
		SSH[ti,] = ptt.maskBathyXY(SSH[ti,], grid, 0)

	#T = tools.barotropicStreamfunction(u, grid)
	T = tools.barotropicStreamfunction(u, grid)

	T2 = tools.barotropicStreamfunction(u, grid)

	umean = - tools.ddy(T, dY)
	umean2 = tools.depthIntegral(u, grid)
	umean = ptt.maskBathyXY(umean, grid, 0, timeDep=True)
	
	vmin = -0.1; vmax = 0.1
	vmin = [vmin, vmin]; vmax = [vmax, vmax]
	g = 9.81; f0 = -1.4e-4	
	SSHstr = g*SSH/f0

	Tmean = np.mean(T, axis=(1,2)); T -= np.tile(Tmean,(1,1,1)).T
	
	#print(np.mean(SSHstr,axis=(1,2)))
	#print(np.mean(SSHstr**2,axis=(1,2))**0.5)
	#print(np.mean(T, axis=(1,2)))
	#print(np.mean(T**2,axis=(1,2))**0.5)

	cmap = 'gist_rainbow'
	pt.plot1by2([T[-1], T[-1]-T2[-1]], mesh=True, cmap=cmap)
	#pt.plot1by2([umean[-1], umean2[-1]], mesh=True, vmin=vmin, vmax=vmax)
	quit()

#==

TEST_depthAverage = False
if TEST_depthAverage:

	path = '/home/michael/Documents/data/MCS_001/run/'

	var = 'Vvel'; VAR = 'VVEL'
	vmin, vmax, cmap, title = getPlottingourcVars(var)

	# Get grid.
	grid = Grid(path)
	X = grid.XC[0,:] / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()	

	# Read data
	data = readVariable(var, path, file_format='nc', meta=True)
	data = data[VAR][0,]

	#v_zav1 = np.trapz(data, Z, axis=0)
	#data = ptt.maskBathyAll(data, grid)
	#v_zav2 = np.trapz(data,-Z, axis=0)
	#v_zav3 = tools.depthAverage(data, Z, timeDep=False)

	
	vstr1 = tools.meridStreamfunction(data, X, Z)

	data = ptt.maskBathyAll(data, grid)
	vstr2 = tools.meridStreamfunction(data, X, Z)
	vstr1 = ptt.maskBathyYZ(vstr2, grid)
	#vstr2 = np.cumsum(np.sum(data[::-1], axis=2), axis=0) * (Z[0]-Z[1]) * (X[1]-X[0]); vstr2 = vstr2[::-1]

	pt.plot1by2([vstr1, vstr2], X=[Y,Y], Y=[Z,Z])	
	
	quit()

#==

TEST_animate = False
if TEST_animate:

	#path = '/home/michael/Documents/data/MCS_002/run/'

	path = '/home/michael/Documents/data/PISOMIP_002/run/'

	grid = Grid(path)

	X = grid.YC[:,1]/1000.
	Y = grid.RC.squeeze()

	#VAR = 'ETAN'
	#VAR = 'THETA'
	VAR = 'SALT'	
	#VAR = 'UVEL'	
	
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	#vmin = 33.32; vmax = 34.5
	data = readVariable(VAR, path, file_format='nc', meta=True, tt=[0,None])
	text_data = ptt.getTextData(data['TIME'][:], 'month', X[1], Y[-2], color='k')
	data = data[VAR][:]	
		
	MEAN = False
	ASYM = False
	if not MEAN:
		for ti in range(data.shape[0]):
			data[ti,] = ptt.maskBathyAll(data[ti,], grid)
		#data = np.ma.mean(data, axis=3)
		data = data[...,120]
		#data = np.mean(data[...,1:40], axis=-1)
	
	else:
		if ASYM:
			for ti in range(data.shape[0]):
				data[ti,] = ptt.maskBathyAll(data[ti,], grid)

		data = np.ma.mean(data, axis=3)

		if not ASYM:
			data = ptt.maskBathyYZ(data, grid, timeDep=True, xi=120)

	#plt.pcolor(np.mean(data[48:,], axis=0), vmin=vmin, vmax=vmax, cmap=cmap); plt.colorbar(); plt.show(); quit()

	data = tools.boundData(data, vmin, vmax, scale=0.99999)

	# PLOT.

	#vmin = None; vmax = None
	xlabel = 'LATITUDE (km)'; ylabel = 'DEPTH (m)'
	#pt.plot1by1(data[-1], X, Y, xlabel=xlabel, ylabel=ylabel, title=title, cmap=cmap, vmin=vmin, vmax=vmax, mesh=False); quit()

	pt.animate1by1(data, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data)
		
	quit()
	
	#==

#==

TEST_animateX = False
if TEST_animateX:

	path = '/home/michael/Documents/data/MCS_103/run/'

	grid = Grid(path)
	X = grid.YC[:,1]/1000.
	Y = grid.RC.squeeze()

	#VAR = 'ETAN'
	#VAR = 'THETA'
	VAR = 'UVEL'	
	#VAR = 'VVEL'
	#VAR = 'WVEL'

	vmin, vmax, cmap, title = getPlottingVars(VAR)
	data = readVariable(VAR, path, file_format='nc', meta=True)

	text_data = ptt.getTextData(grid.XC[1,:], 'X', X[1], Y[-2], color='k')
	data = np.mean(data[VAR][-12:], axis=0)

	data = ptt.maskBathyAll(data, grid)
	data = np.transpose(data, (2,0,1))	

	data = tools.boundData(data, vmin, vmax, scale=0.99999)
	
	# PLOT.
	xlabel = 'LATITUDE (km)'; ylabel = 'DEPTH (m)'
	pt.animate1by1(data, X, Y, vmin=vmin, vmax=vmax, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data)

	quit()

	
#==

animateSurface = True	
if animateSurface:

	#path = '/home/michael/Documents/data/PAS_666/run/'
	#path = '/home/michael/Documents/data/PISOMIP_003/run/'
	path_root = '/home/michael/Documents/data/'
	run = 'MCS_320'
	
	path = path_root + run + '/run/'
	grid = Grid(path)
	#grid = Grid_PAS(path)
	bathy = grid.bathy
	#pt.plot1by1(bathy, vmin=-1000, vmax=-300, mesh=True); quit()

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	xlabel = 'LON (km)'; ylabel = 'LAT (km)'
	ny, nx = bathy.shape
	
	#VAR = 'ETAN'
	#VAR = 'PHIBOT'
	#VAR = 'RHOAnoma'
	#VAR = 'THETA' 
	#VAR = 'PHIHYD'
	#VAR = 'WVELTH';
	#VAR = 'WVELTH'#','UVELTH','VVELTH','WVELTH', 'TOTTTEND'
	#VAR = 'SALT'	
	#VAR = 'UVEL'
	#VAR = 'VVEL'
	#VAR = 'VVEL'
	#VAR = 'WVEL'
	#VAR = 'botTauX'
	#VAR = 'ADTAUU'
	#VAR = 'ISOTHERM'
	#VAR = 'oceSflux'
	VAR = 'oceTAUX'

	flatVars = ['ETAN', 'botTauX', 'PHIBOT', 'ISOTHERM', 'oceTAUX', 'oceSflux']

	if VAR == 'ISOTHERM':
	
		data = readVariable('ETAN', path, file_format='nc', meta=True)
		text_data = ptt.getTextData(data.variables['TIME'][:], 'month', X[1], Y[1], color='k')
		
		data = np.load(path+'ThermZ_m05_'+run+'.npy')
		vmin = -600; vmax=-150;	cmap = 'YlOrRd'
		#vmin = -500; vmax = -350; cmap = 'jet'
		
		title = run + ' -0.5 deg. C isotherm depth'
		
	else:
			
		vmin, vmax, cmap, title = getPlottingVars(VAR)

		data = readVariable(VAR, path, file_format='nc', meta=True)
		print(data.variables)

		text_data = ptt.getTextData(data.variables['TIME'][:], 'month', X[1], Y[1], color='k')
		data = data[VAR][:]
		#plt.plot(data[-1,:,1]);plt.grid(); plt.show(); quit()

	#==

	# DON'T CHANGE THIS. CHANGE ONE IN IF STATEMENT IF YOU WANT.
	level = 0

	if VAR not in flatVars:
		level = 0
		data = data[:,level]
		print('Z = ' + str(grid.RC.squeeze()[level]))

	#data = tools.boundData(data, vmin, vmax, scale=0.99999)
	data = ptt.maskBathyXY(data, grid, level, timeDep=True)
	hfac=grid.hFacC[0]*grid.hFacW[0]*grid.hFacS[0]
	print(data.shape)
	print(np.ma.sum(data[0]))
	print(np.ma.sum(data[0,:100]*hfac[:100]))
	print(np.ma.sum(data[0,100:]*hfac[100:]))
	
	#vmin = None; vmax = None
	#vmin = -0.025; vmax = -vmin # For SSS
	
	pt.animate1by1(data, X, Y, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, outname='animate1by1Surf.mp4', vmin=vmin, vmax=vmax, contour=grid.bathy)

	quit()
	
	
#==

animateCONV = False
if animateCONV:

	path = '/home/michael/Documents/data/MCS_123/run/'

	grid = Grid(path)
	#grid = Grid_PAS(path)
	bathy = grid.bathy
	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	
	u = readVariable('UVEL', path, file_format='nc', meta=False)
	v = readVariable('VVEL', path, file_format='nc', meta=False)

	u = tools.interp(u, 'u')
	v = tools.interp(v, 'v')

	level = -1
	if level == -1:
		#u = np.mean(u, axis=1)
		#v = np.mean(v, axis=1)
		
		u = np.sum(u, axis=1)
		v = np.sum(v, axis=1)
		#u = tools.depthIntegral(u, grid, norm=False)
		#v = tools.depthIntegral(v, grid, norm=False)
	else:
		u = u[:,level]
		v = v[:,level]

	uu = u*u
	uv = u*v

	conv = - tools.ddx(uu, grid.DXG) - tools.ddy(uv, grid.DYG)
	conv = ptt.maskBathy(conv, grid, zi=0, timeDep=True)

	#==
	
	text = ['month ' + str(ti) for ti in range(conv.shape[0])]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}
	
	#pt.plot1by1(conv[-1])
	
	vmin = -1.e-6; vmax = -vmin
	
	pt.plot1by1(20*conv[-1]); quit()
	
	pt.animate1by1(conv, X, Y, cmap='bwr', mesh=True, text_data=text_data, outname='animateConv.mp4', vmin=vmin, vmax=vmax)
	
	quit()
	
#==
	
animateRV = False
if animateRV:

	#path = '/home/michael/Documents/data/PAS_666/run/'
	path = '/home/michael/Documents/data/MCS_123/run/'
	#path = '/data/oceans_output/shelf/pahol/mitgcm/PAS_851/run/'

	grid = Grid(path)
	#grid = Grid_PAS(path)	

	#level = 17 # 17 -> 350m depth in MCS.
	level = grid.getIndexFromDepth(-490)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	dx = (X[1]-X[0])*1000.
	dy = (Y[1]-Y[0])*1000.
	xlabel = 'LON (km)'; ylabel = 'LAT (km)'

	ts = 0

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	u = readVariable(VAR, path, file_format='nc', meta=True)
	TIME = u['TIME'][ts:]
	u = u[VAR][ts:, level]

	VAR = 'VVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	v = readVariable(VAR, path, file_format='nc', meta=False)[ts:, level]

	RV = tools.ddy(u, dy) - tools.ddx(v, dx)
	
	tscale = 86400. * 30
	tscale_name = 'month'
	t = TIME / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

	contour = grid.bathy
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)
	RV = ptt.maskBathyXY(RV, grid, 0, timeDep=True)
	
	vmin = -1.e-5; vmax = -vmin
	RV = tools.boundData(RV, vmin, vmax)
	
	pt.animate1by1(RV, X, Y, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, outname='animate1by1Surf.mp4', vmin=vmin, vmax=vmax)

	quit()

#==

animateQuivers = False
if animateQuivers:

	#path = '/home/michael/Documents/data/PAS_666/run/'
	path = '/home/michael/Documents/data/MCS_315/run/'
	#path = '/data/oceans_output/shelf/pahol/mitgcm/PAS_851/run/'

	# Sample rate
	d = 8

	grid = Grid(path)
	#grid = Grid_PAS(path)	

	#level = 17 # 17 -> 350m depth in MCS.
	level = grid.getIndexFromDepth(-450)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	data1 = readVariable(VAR, path, file_format='nc', meta=True)
	TIME = data1['TIME'][:]
	data1 = data1[VAR]

	VAR = 'VVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	data2 = readVariable(VAR, path, file_format='nc', meta=False)

	T = readVariable('THETA', path, file_format='nc', meta=False)
	cvmin, cvmax, cmap, title = getPlottingVars(VAR)
	
	tscale = 86400. * 30
	tscale_name = 'month'
	t = TIME / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}

	contour = grid.bathy
	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)
	
	INTERVAL = False
	if INTERVAL:
		ts = 0; te = 10		
		data1 = data1[ts:te, level]; data2 = data2[ts:te, level]; T = T[ts:te, level]
	else:
		data1 = data1[:, level]; data2 = data2[:, level]; T = T[:, level]


	SUBR = False
	if SUBR:
		#lats = [-75, -71]; lons = [240, 270]
		lats = [-76, -70.5]; lons = [245, 260]#[230, 270]#
		latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)
		ThermZ = tools.getSubregionXY(ThermZ, latsi, lonsi)
		bathy = tools.getSubregionXY(bathy, latsi, lonsi)
		X = grid.Xsubr1D(lons)
		Y = grid.Ysubr1D(lats)

	data1 = tools.interp(data1, 'u'); data2 = tools.interp(data2, 'v')

	data1 = ptt.maskBathyXY(data1, grid, level, timeDep=True)
	data2 = ptt.maskBathyXY(data2, grid, level, timeDep=True)

	title = '(u, v) & bathymetry; Z = ' + str(Z[level]) + ' m'

	data1 = data1[..., ::d, ::d]; data2 = data2[..., ::d, ::d]
	Xd = X[::d]; Yd = Y[::d]
	Td = T[..., ::d, ::d]; Td[:,-1,-1] = cvmin; Td[:,-1,-2] = cvmax
	
	#pt.plot1by2([data1[0], data2[0]]); quit()

	pt.animate1by1quiver(data1, data2, Xd, Yd, C=Td, contour=contour, X=X, Y=Y, contourf=False, text_data=text_data, title=title, figsize=(7,6), cmap='YlGn')
	
	quit()

#==

animateSSHs = False
if animateSSHs:
	# Animate zonal mean SSHs for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michael/Documents/data/'

	# These are cases without CS walls.
	#paths = [path_root+'MCS_013/run/', path_root+'MCS_018/run/', path_root+'MCS_012/run/', path_root+'MCS_016/run/']
	# These are cases with CS walls.	
	paths = [path_root+'MCS_101/run/', path_root+'MCS_104/run/', path_root+'MCS_102/run/', path_root+'MCS_103/run/']
	labels = ['wind0 rel300', 'wind0 rel200', 'wind16 rel300', 'wind16 rel200']
	
	grid = Grid(paths[0])
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'zonal mean SSH (m)'

	VAR = 'ETAN'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, paths[0], file_format='nc', meta=True)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}

	# Add first zonal mean SSH to data list.
	data_tmp = ptt.maskBathyXY(data_tmp[VAR][:], grid, 0, timeDep=True)
	data.append(np.mean(data_tmp[:,1:-1], axis=2))
	#data.append(data_tmp[:,1:-1,1])

	# Now get other SSHs.
	for di in range(len(paths)-1):
		data_tmp = readVariable(VAR, paths[di+1], file_format='nc', meta=False)
		data_tmp = ptt.maskBathyXY(data_tmp, grid, 0, timeDep=True)
		data.append(np.mean(data_tmp[:,1:-1], axis=2))
		#data.append(data_tmp[:,1:-1,1])

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data)
	quit()

#==

animateSSHs1sim = False
if animateSSHs1sim:
	# Animate SSH at different longitudes for a single MITgcm run. 

	path_root = '/home/michael/Documents/data/'
	path = path_root + 'MCS_034/run/'
	#xis = [1, 190, 80, 120]
	xis = [1, 5, 230, 239]
	labels = ['x = ' + str(xi*2.5) + ' km' for xi in xis]

	grid = Grid(path)
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'zonal mean SSH (m)'

	VAR = 'ETAN'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, path, file_format='nc', meta=True)


	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}
	
	data_tmp = data_tmp[VAR][:]

	data2 = readVariable('VVEL', path, file_format='nc', meta=True)['VVEL'][:]
	data_tmp = (20. + data_tmp) * data2[:,0]

	# Now get other SSHs.
	for xi in xis:
		data.append(data_tmp[...,1:-1,xi])


	vmin = -2e-1; vmax = -vmin
	pt.animateLine(data, X=Y[1:-1], xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, vmin=vmin, vmax=vmax)
	quit()

#==

animateIsotherms = False
if animateIsotherms:
	# Animate zonal mean Isotherm Height for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michael/Documents/data/'

	# These are cases with CS walls.	
	#paths = [path_root+'MCS_105/run/', path_root+'MCS_106/run/']; minus_dirs = 1
	#labels = ['wind16 rel300 x = 0', 'wind16 rel300 x = 120', 'wind16 rel200 x = 0', 'wind16 rel200 x = 120']

	# These are fully periodic sims.	
	paths = [path_root+'MCS_101/run/', path_root+'MCS_104/run/', path_root+'MCS_102/run/', path_root+'MCS_103/run/']; minus_dirs = 1
	labels = ['wind0 rel300', 'wind0 rel200', 'wind16 rel300', 'wind16 rel200']

	grid = Grid(paths[0])
	bathy = grid.bathy
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()
	xlabel = 'Y (km)'
	ylabel = '-0.5 deg. isotherm depth (m)'

	THERM = -0.5
	VAR = 'THETA'
	vmin = -500; vmax = -200
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, paths[0], file_format='nc', meta=True)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}
	
	data_tmp = data_tmp[VAR][:]

	# Get z-indices of level with Theta closest to THERM.

	# Two zonal slices.
	if minus_dirs == 3:	
		ThermZ = tools.getIsothermHeight(data_tmp, THERM, grid, interp=True)
		data.append(ThermZ[...,0]); data.append(ThermZ[...,120]);
	elif minus_dirs == 1:
		Nx = data_tmp.shape[-1]
		data_tmp = np.mean(data_tmp, axis=-1)
		data_tmp = np.tile(data_tmp, (Nx,1,1,1)).transpose([1,2,3,0])
		ThermZ = tools.getIsothermHeight(data_tmp, THERM, grid, interp=True)
		data.append(ThermZ[...,0])

	# Now get other SSHs.
	for di in range(len(paths)-minus_dirs):
		data_tmp = readVariable(VAR, paths[di+1], file_format='nc', meta=False)

		# Two zonal slices.
		if minus_dirs == 3:	
			ThermZ = tools.getIsothermHeight(data_tmp, THERM, grid, interp=True)
			data.append(ThermZ[...,0]); data.append(ThermZ[...,120]);
		elif minus_dirs == 1:
			Nx = data_tmp.shape[-1]
			data_tmp = np.mean(data_tmp, axis=-1)
			data_tmp = np.tile(data_tmp, (Nx,1,1,1)).transpose([1,2,3,0])
			ThermZ = tools.getIsothermHeight(data_tmp, THERM, grid, interp=True)
			data.append(ThermZ[...,0])

	pt.animateLine(data, X=Y, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data)
	quit()

#==

animateUbots = False
if animateUbots:
	# Animate zonal mean Ubots for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michael/Documents/data/'

	# These are cases with CS walls.	
	paths = [path_root+'MCS_103/run/', path_root+'MCS_104/run/', path_root+'MCS_105/run/', path_root+'MCS_106/run/']
	labels = ['wind16 rel300 x = 0', 'wind16 rel300 x = 120', 'wind16 rel200 x = 0', 'wind16 rel200 x = 120']
	
	grid = Grid(paths[0])
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'bottom zonal vel. (m/s)'

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	u = readVariable(VAR, paths[0], file_format='nc', meta=True)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = u.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}

	# Get bottom val
	ub = tools.bottom(u[VAR][:], grid, 'u')
	print(ub.shape)

	# COMMENT OUT AS APPROPRIATE

	# Two bottom velocities
	data.append(ub[:,1:-1,0]); data.append(ub[:,1:-1,120])

	# Zonal-mean bottom vel.
	#data.append(np.mean(ub[:,1:-1,:], axis=-1))

	# Now get other SSHs.
	for di in range(len(paths)-3):
		ub = readVariable(VAR, paths[di+1], file_format='nc', meta=False)
		ub = tools.bottom(ub, grid, 'u')
		print(data_tmp.shape)
		#data.append(np.mean(ub[:,1:-1,:], axis=-1))
		data.append(ub[:,1:-1,0]); data.append(ub[:,1:-1,120])

	print(data[-1].shape)

	vmin *= vmax; vmax *= vmax

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data)
	quit()


#==


animateDragBot = True
if animateDragBot:
	# Animate zonal mean Ubots for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michael/Documents/data/'

	# These are cases with CS walls.	
	#paths = [path_root+'MCS_113/run/', path_root+'MCS_104/run/']
	#labels = ['wind16 botStress', 'wind0 botStress', 'wind0 x = 0', 'wind0 x = 120']
	
	paths = [path_root+'MCS_104/run/', path_root+'MCS_104/run/']
	labels = ['wind0 botStress', 'wind0 botStress', 'wind0 x = 0', 'wind0 x = 120']

	grid = Grid(paths[0])
	Y = grid.YC
	xlabel = 'Y (km)'
	ylabel = 'bottom stress'
	Ny, Nx = Y.shape

	Y = Y[:,1]/1000.

	Cd = 2.5e-3
	rho0 = 1030.0
	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	u = readVariable(VAR, paths[0], file_format='nc', meta=True)

	u16 = 0.025 * tools.getTranslatedWind(Nx, Ny, 16) / rho0
	u0 = 0.025 * tools.getTranslatedWind(Nx, Ny, 0) / rho0
	#u_16 = 0.025 * tools.getTranslatedWind(Nx, Ny, -16) / rho0
	#plt.plot(u16[:,0], Y, label='16');plt.plot(u_16[:,0], Y, label='-16'); plt.legend(); plt.show(); quit()

	tscale = 86400. * 30
	tscale_name = 'month'
	t = u.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}

	drag = -tools.computeBotDragQuadr0(paths[0], grid)
	#drag = -readVariable('botTauX', paths[0], file_format='nc', meta=False)/rho0

	u16 = ptt.maskBathyXY(u16, grid, 0, timeDep=False)
	drag = ptt.maskBathyXY(drag, grid, 0, timeDep=True)

	# COMMENT OUT AS APPROPRIATE

	# Two quadr drags
	#data.append(drag[:,1:-1,0]); data.append(drag[:,1:-1,120])

	# Zonal mean
	data.append(np.ma.mean(drag[:,1:-1,:], axis=-1))

	# Now get other SSHs.
	for di in range(len(paths)-1):

		print(paths[di+1])
		drag = -tools.computeBotDragQuadr0(paths[di+1], grid)

		# COMMENT OUT AS APPROPRIATE
		# Two quadr drags
		#data.append(drag[:,1:-1,0]); data.append(drag[:,1:-1,120])
		# Zonal mean
		data.append(np.ma.mean(drag[:,1:-1,:], axis=-1))

	#==

	print(len(data))
	vmin=-4e-5; vmax=-vmin
	constLineLabel = ['wind16', 'wind0']

	u0 = np.ma.sum(u0, axis=-1) / Nx
	u16 = np.ma.sum(u16, axis=-1) / Nx

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[u16[1:-1], u0[1:-1]], constLineLabel=constLineLabel)
	quit()


#==

animateUvelSurfs = False
if animateUvelSurfs:
	# Animate zonal mean SSHs for multiple MITgcm runs.
	# Assumes that all use the same grid and have the same time dependencies.

	path_root = '/home/michael/Documents/data/'

	# These are cases without CS walls.
	#paths = [path_root+'MCS_013/run/', path_root+'MCS_018/run/', path_root+'MCS_012/run/', path_root+'MCS_016/run/']
	# These are cases with CS walls.	
	#paths = [path_root+'MCS_015/run/', path_root+'MCS_019/run/', path_root+'MCS_014/run/', path_root+'MCS_017/run/']
	#labels = ['05cos', '025cos', '025sin2', '0125sin2']

	# These are for comparing with and without beta.
	paths = [path_root+'MCS_015/run/', path_root+'MCS_022/run/', path_root+'MCS_014/run/', path_root+'MCS_023/run/']
	labels = ['05cos wall', '05cos wall beta', '025sin2 wall', '025sin2 wall beta']
	
	grid = Grid(paths[0])
	Y = grid.YC[:,1]/1000.
	xlabel = 'Y (km)'
	ylabel = 'zonal mean Uvel (m)'

	VAR = 'UVEL'
	vmin, vmax, cmap, title = getPlottingVars(VAR)
	
	data = []

	# Get metadata from first run.
	data_tmp = readVariable(VAR, paths[0], file_format='nc', meta=True)

	tscale = 86400. * 30
	tscale_name = 'month'
	t = data_tmp.variables['TIME'][:] / tscale
	text = [tscale_name + ' ' + str(int(tt)) for tt in t]
	text_data = {'text':text, 'xloc':Y[1], 'yloc':vmin, 'fontdict':{'fontsize':14, 'color':'k'}}

	# Add first zonal mean SSH to data list.
	data_tmp = ptt.maskBathyXY(data_tmp[VAR][:,0], grid, 0, timeDep=True)
	data.append(np.mean(data_tmp[:,1:-1,-100:], axis=2))
	
	# Now get other SSHs.
	for di in range(len(paths)-1):
		data_tmp = readVariable(VAR, paths[di+1], file_format='nc', meta=False)
		data_tmp = ptt.maskBathyXY(data_tmp[:,0], grid, 0, timeDep=True)
		data.append(np.mean(data_tmp[:,1:-1,-100:], axis=2))

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data)
	quit()


#==

#==

# yzSlice. Contour and contourf of two fields, e.g., Uvel and Theta.
yzSlice = False
if yzSlice:

	#path = '/home/michael/Documents/data/MCS_002/run/'

	path = '/home/michael/Documents/data/MCS_036/run/'
	#path = '/home/michael/Documents/data/PISOMIP_003/run/'
	#pathG = '/home/michael/Documents/data/MCS_018/run/'

	grid = Grid(path)
	#grid = Grid_PAS(path)
	X = grid.YC[:,1]/1000.
	Y = grid.RC.squeeze()
	
	ts = 30
	VAR1 = 'UVEL' # Contourf
	VAR2 = 'THETA' # Contour

	vmin, vmax, cmap, title = getPlottingVars(VAR1)

	data1 = readVariable(VAR1, path, file_format='nc', meta=True)
	data2 = readVariable(VAR2, path, file_format='nc', meta=False)
	print(data2.shape)

	text_data = ptt.getTextData(data1.variables['TIME'][:], 'month', X[1], Y[-2], color='w')

	MEAN = False
	if not MEAN:
		data1 = data1[VAR1][:]
		#data1 = ptt.maskBathyAll(data1[ts,], grid)
		#data2 = ptt.maskBathyAll(data2[ts,], grid)
		data1 = ptt.maskBathyAll(np.mean(data1[ts:,], axis=0), grid)
		data2 = ptt.maskBathyAll(np.mean(data2[ts:,], axis=0), grid)
		#data = np.ma.mean(data, axis=3)
		data1 = data1[...,0]
		data2 = data2[...,0]
		#data = np.mean(data[...,1:40], axis=-1)
	
	else:
		data1 = np.mean(data1[VAR1][:], axis=3)
		data1 = ptt.maskBathyYZ(data1, grid, timeDep=True, xi=120)
		data2 = np.mean(data2[VAR2][:], axis=3)
		data2 = ptt.maskBathyYZ(data2, grid, timeDep=True, xi=120)

	data1 = tools.boundData(data1, vmin, vmax, scale=0.99999)
	data2 = tools.boundData(data2, vmin, vmax, scale=0.99999)
	
	# PLOT.

	#vmin = None; vmax = None
	xlabel = 'LATITUDE (km)'; ylabel = 'DEPTH (m)'
	title = VAR1 + ' & ' + VAR2
	#pt.plot1by1(data[-1], X, Y, xlabel=xlabel, ylabel=ylabel, title=title, cmap=cmap, vmin=vmin, vmax=vmax, mesh=False); quit()
	
	# YZ plot.
	pt.plot1by1(data1, X=X, Y=Y, title=title, xlabel=xlabel, ylabel=ylabel, cmap=cmap, contour=data2, contourlevels=3, show=True)
	
	quit()

#==

plotBathy = False
if plotBathy:

	#path = '/home/michael/Documents/data/PAS_666/run/'
	path = '/home/michael/Documents/data/MCS_038/run/'
	
	grid = Grid(path)
	bathy = grid.bathy

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	xlabel = 'LON (km)'; ylabel = 'LAT (km)'

	bathy = ptt.maskBathyXY(bathy, grid, 0, timeDep=False)

	pt.plot1by1(bathy, X=X, Y=Y, xlabel=xlabel, ylabel=ylabel, mesh=True); quit()

#==

TEST_troughTransport = True
if TEST_troughTransport:

	# Look at volume transport through window west of trough and possibly east of trough.
	# Compare with volume transport through trough and nearby upwelling.
	# If we want time series of volume transport through these slices, have to take care with grid stencil and
	# variable cell thicnkesses. Take a look at xgcm if needed.
	
	# Remaining to do:
	# w through an x,y slice.
	# u through a second y,z slice.
	# volume transport through these slices. (scale by hfac, dz; use u and v)
	# Time series of these.
	# Onto github, onto HPC.
	
	path = '/home/michael/Documents/data/MCS_123/run/'
	grid = Grid(path)
	
	#path = '/home/michael/Documents/data/PAS_851/run/'
	#grid = Grid_PAS(path)
	
	# Four windows:
	lat1 = [-71.8, -71.2]; lon1 = 360-114.5; depth1 = [-250, -750]
	lat2 = -71.6; lon2 = [360-115, 360-112]; depth2 = [-250, -750]
	lat3 = [-71.8, -71.2]; lon3 = [360-115, 360-112]; depth3 = -500
	lat4 = lat1; lon4 = 360-112; depth4 = depth1

	# Convert physical values to indices.
	lats1 = grid.getIndexFromLat(lat1); depths1 = grid.getIndexFromDepth(depth1)
	xi1= grid.getIndexFromLon(lon1)
	
	lons2 = grid.getIndexFromLon(lon2); depths2 = grid.getIndexFromDepth(depth2)
	yi2 = grid.getIndexFromLat(lat2)
	
	lats3 = grid.getIndexFromLat(lat3); lons3 = grid.getIndexFromLon(lon3)
	zi3 = grid.getIndexFromDepth(depth3)
	
	lats4 = grid.getIndexFromLat(lat4); depths4 = grid.getIndexFromDepth(depth4)
	xi4= grid.getIndexFromLon(lon4)
	
	#==
	
	# Get zonal transport
	u = readVariable('UVEL', path, meta=False)[ts:]
	ut = tools.zonalTransport(data, grid)
	
	# Zonal transport through windows 1 and 4.
	ut1 = tools.getSubregionYZ(ut[..., xi1], grid, lats1, depths1)
	ut4 = tools.getSubregionYZ(ut[..., xi4], grid, lats4, depths4)

	# Get meridional transport through window 2.
	fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	vt = tools.meridTransport(data, grid)
	vt2 = tools.getSubregionXZ(vt[:, :, yi2, ], grid, lons2, depths2)
	
	# Get vertical transport through window 3.
	#fname = 'stateWvel.nc'; var = 'WVEL'; cmap = 'coolwarm'; vmax = 0.01; vmin = -vmax
	#ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	wt = tools.meridTransport(data, grid)
	wt3 = tools.getSubregionXY(wt[:, zi3, ], grid, lats3, lons3)
		
	# Transports are already scaled by grid size, so can just sum to get total transport.
	ut1 = np.sum(ut1, axis=(1,2))
	ut4 = np.sum(ut4, axis=(1,2))
	vt2 = np.sum(vt2, axis=(1,2))
	wt3 = np.sum(wt3, axis=(1,2))

	data = [ut1, ut4, vt2, wt3]
	labels = ['ut window 1', 'ut window 4', 'vt window 2', 'wt window 3']
	TIME = ncfile.variables['TIME'][:]
	pt.timeSeries(data, TIME=TIME, labels=labels, ylabel='Volume transport')
	
#==

TEST_manyWindows = False
if TEST_manyWindows:

	ti = -1
	
	# Window 1
	lat1 = [-71.8, -71.2]; lon1 = 360-114.5; depth1 = [-250, -750]
	
	# Window 2
	lat2 = -71.6; lon2 = [360-115, 360-112]; depth2 = [-250, -750]
	
	# Window 3.
	lat3 = [-71.8, -71.2]; lon3 = [360-115, 360-112]; depth3 = -500
	
	# Window 4.
	lat4 = lat1; lon4 = 360-112; depth4 = depth1
	
	#==
	
	# 1. Get zonal flow slice through window 1 and window4
	lats1 = grid.getIndexFromLat(lat1)
	depths1 = grid.getIndexFromDepth(depth1)
	xi1= grid.getIndexFromLon(lon1)
	
	lats4 = grid.getIndexFromLat(lat4)
	depths4 = grid.getIndexFromDepth(depth4)
	xi4= grid.getIndexFromLon(lon4)
	
	fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	u1 = tools.getSubregionYZ(data[ti, ..., xi1], grid, lats1, depths1)
	u4 = tools.getSubregionYZ(data[ti, ..., xi4], grid, lats4, depths4)
	
	# Zonal transport.
	#u_transport = tools.zonalTransport(data, grid)
	
	#==
	
	# 2. Get meridional flow through window 2.
	lons2 = grid.getIndexFromLon(lon2)
	depths2 = grid.getIndexFromDepth(depth2)
	yi = grid.getIndexFromLat(lat2)
	
	fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	v = tools.getSubregionXZ(data[ti, :, yi, ], grid, lons2, depths2, Ydim=False)
	
	#==
	
	# 3. Vertical flow through window 3 (don't have NC file!)
	lats3 = grid.getIndexFromLat(lat3)
	lons3 = grid.getIndexFromLon(lon3)
	zi = grid.getIndexFromDepth(depth3)
	
	#fname = 'stateWvel.nc'; var = 'WVEL'; cmap = 'coolwarm'; vmax = 0.01; vmin = -vmax
	#ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]
	w = tools.getSubregionXY(data[ti, zi, ], grid, lats3, lons3)

	#==
	
	# Mask bathymetry 
	u1 = ptt.maskBathyYZ(u1, grid, xi=xi1, subregion=True, lats=lats1, depths=depths1)
	u4 = ptt.maskBathyYZ(u4, grid, xi=xi4, subregion=True, lats=lats4, depths=depths4)
	v = ptt.maskBathyXZ(v, grid, yi=yi, subregion=True, lons=lons2, depths=depths2, Ydim=False)
	w = ptt.maskBathyXY(w, grid, zi, subregion=True, lons=lons3, lats=lats3)
		
	# Get 1D subregion arrays. Y,Z for u-plot and X,Z for v-plot. 
	Xsubr = grid.Xsubr1D(lon2)
	Ysubr = grid.Ysubr1D(lat1)
	Zsubr = grid.Zsubr1D(depth1)
	
	pt.plot1by2([u1, v], X=[Ysubr, Xsubr], Y=[Zsubr, Zsubr], titles=['u window1', 'v window2'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lat', 'lon'], ylabels=['depth', ''])

	pt.plot1by2([u1, u4], X=[Ysubr, Ysubr], Y=[Zsubr, Zsubr], titles=['u window1', 'u window4'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lat', 'lat'], ylabels=['depth', ''])
	
	pt.plot1by2([w, v], X=[Xsubr, Xsubr], Y=[Ysubr, Zsubr], titles=['v window3', 'v window2'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lon', 'lon'], ylabels=['lat', 'depth '])
	
	quit()

#==

TEST_getSubregionXYZ = False
if TEST_getSubregionXYZ:
	
	# In this test, produce 3D (allow for additional time dimension) fields in subregion.
	# Apply to zonal flows over PITW shelf break.
	# Volume transport or just zonal velocity?
	# Hoevmoller diagram of U(z), could be averaged over some small y-range.

	None

#==

TEST_getSubregionYZ = False
if TEST_getSubregionYZ:

	# Current mask functions require full 2D slice.
	# Extend this so it can take partial 2D slice (have to get corresponding subregion in grid).
	# And extend so can work with general number of dimensions.

	# Load potential temp. Other options above.
	fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
	ncfile = nc.Dataset(path+fname, 'r'); data = ncfile.variables[var]

	ti = -1
	xi = grid.getIndexFromLon(245.5)
	
	#depth_vals = [0, -2000]; lat_vals = [-72.5, -70.5]
	
	# AssmanEtal2013 vals
	depth_vals = [0, -2000]; lat_vals = [-75, -65]; lon_vals = [360-116, 360-111]
	
	lats = grid.getIndexFromLat(lat_vals, xi=xi)
	depths = grid.getIndexFromDepth(depth_vals)
		
	data1 = tools.getSubregionYZ(data, grid, lats, depths)[ti, ..., xi]
	data2 = tools.getSubregionYZ(data[ti, ..., xi], grid, lats, depths)

	data1 = ptt.maskBathyYZ(data1, grid, xi=xi, subregion=True, lats=lats, depths=depths)
	data2 = ptt.maskBathyYZ(data2, grid, xi=xi, subregion=True, lats=lats, depths=depths)
		
	Ysubr = grid.Ysubr1D(lat_vals)
	Zsubr = grid.Zsubr1D(depth_vals)
	#data1 = ptt.maskBathyYZ(data1, grid, xi=xi)

	pt.plot1by2([data1, data2], X=[Ysubr, Ysubr], Y=[Zsubr, Zsubr], titles=['data1', 'data2'], cmap=cmap, vmin=vmin, vmax=vmax, xlabels=['lat', 'lat'], ylabels=['depth', ''], show=False, save=True)

#==

TEST_getSubregiondiXY = False
if TEST_getSubregionXY:

	ti = -1
	ki = 10

	lons = (245, 260)
	lats = (-72.5, -70.5)
	data1 = tools.getSubregionXY(data, grid, lons, lats)[ti, ki]
	data2 = tools.getSubregionXY(data[ti, ki, ], grid, lons, lats)

	pt.plot1by2(X, Y, data1, data2, titles=['data1', 'data2'])

#==

TEST_btpcStrKN = False
if TEST_btpcStrKN:

	import utils 
	from grid_KN import Grid as Grid_KN
	
	path = '/home/michael/Documents/data/PAS_851/run/'
	grid = Grid_PAS(path)
	u = np.load(path+'umean_PAS851.npy')
	
	#path = '/home/michael/Documents/data/MCS_038/run/'
	#grid = Grid(path)

	gridK = Grid_KN(path)

	#u = readVariable('UVEL', path, file_format='nc', meta=False)[-10:]
	#u = np.mean(u, axis=0)

	X = grid.XC[1,:]
	Y = grid.YC[:,1]
	Z = grid.RC.squeeze()

	u = ptt.maskBathyAll(u, grid, timeDep=False)

	#T = tools.barotropicStreamfunction(u, grid)
	T = - 1.e-6 * tools.barotropicStreamfunction(u, grid, timeDep=False, norm=False)
	TK = utils.barotropic_streamfunction(u, gridK)

	T = ptt.maskBathyXY(T, grid, 0, timeDep=False)

	pt.plot1by2([T, TK], mesh=True)
	quit()

#==

TEST_bathy = False
if TEST_bathy:

	from grid_KN import Grid as Grid_KN
 
	path = '/home/michael/Documents/data/MCS_038/run/'
	gridP = Grid_PAS(path)
	grid = Grid(path)
	gridK = Grid_KN(path)

	#plt.plot(gridP.RF.squeeze()); plt.plot(grid.RF.squeeze()); plt.show(); quit()

	bathy = grid.bathy; draft = grid.draft
	bathyP = gridP.bathy; draftP = gridP.draft
	bathyK = gridK.bathy; draftK = gridK.draft

	landC = grid.landC
	iceC = grid.iceC

	pt.plot1by2([bathyP, bathyK-bathyP], titles=['Bathymetry', 'KN Bathymetry'], mesh=True)
	pt.plot1by2([draftP, draftK-draftP] , titles=['Draft', 'KN Draft'], mesh=True)

	quit()
	
	
#==
	
TEST_getIndex = False
if TEST_getIndex:

	depths = [0, -2000]
	lats = [-72.5, -70.5]
	lons = [240, 250]
	
	print(grid.getIndexFromLat(lats, xi=0))
	print(grid.getIndexFromLat(lats[0], xi=0))
	print(grid.getIndexFromLat(lats[1], xi=0))

	print(grid.getIndexFromDepth(depths))
	print(grid.getIndexFromDepth(depths[0]))
	print(grid.getIndexFromDepth(depths[1]))
		
	print(grid.getIndexFromLon(lons, yi=0))
	print(grid.getIndexFromLon(lons[0], yi=0))
	print(grid.getIndexFromLon(lons[1], yi=0))
	
	quit()
	

	
	
#==

TEST_readData = False
if TEST_readData:

	path =  '/Users/mh115/Documents/BAS/data/PISOMIP_001/run/'
	grid = Grid(path)
	
	var1 = readVariable('Theta', path, file_format='rdmds')
	var2 = readVariable('Theta', path, file_format='nc')

	#==
	
	ts = 1
		
	# Lat-depth plot.
	data1 = np.mean(var1[ts], axis=2); data2 = np.mean(var2[ts], axis=2)
	data1 = ptt.maskBathyYZ(data1, grid, xi=10); data2 = ptt.maskBathyYZ(data2, grid, xi=10)
	print(data1.shape)
	data1 = ptt.maskDraftYZ(data1, grid, xi=10); data2 = ptt.maskDraftYZ(data2, grid, xi=10)
	
	pt.plot1by2([data1, data2], X=[grid.YC[:,0], grid.YC[:,0]], Y=[grid.RC.squeeze(), grid.RC.squeeze()], titles=['rdmds', 'nc'], xlabels=['lat', 'lat'], ylabels=['depth', ''], show=True, figsize=(8,3))
		
	quit()
	
	#==
	
	# Lat-lon plot.
	pt.plot1by2([var1[-1,0], var2[-1,0]], X=[grid.XC, grid.XC], Y=[grid.YC, grid.YC], titles=['rdmds', 'nc'], xlabels=['lon', 'lon'], ylabels=['lat', 'lat'], show=True)




