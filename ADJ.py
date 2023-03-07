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

#==========================================================


objTest = 0
if objTest:

	path = '/home/michael/Documents/data/MCS/run/'

	VAR1 = 'stateTheta.0000000000.data'; norm1 = 1.e9
	VAR2 = 'm_boxmean_theta.0000000000.data'; norm2 = 1.e0

	nx = 240; ny = 200
	dims = (nx, ny)
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	
	vol = grid.hFacC * grid.DRF * grid.RAC
	
	data1 = readnp(path+VAR1, dims, rec=-1)# / norm1
	
	# Integrate over volume and normalise.
	data1 = np.sum(vol * data1, axis=0)	
	totvol = np.sum(vol)
	data1 = data1 / totvol	
	
	sum1 = np.sum(data1)#*norm1

	data2 = readnp(path+VAR2, dims, rec=0)# / norm2
	sum2 = np.sum(data2)#*norm2
	
	print(sum1)
	print(sum2)
	print(sum1/sum2)

	pt.plot1by2([data1, data2-data1])
	
	quit()
	
#==

OBJtimeseries = 0
if OBJtimeseries:

	path_root = '/home/michael/Documents/data/'
	paths = [path_root+'MCS_202/run/', path_root+'MCS_201/run/']
	
	d = 1; ylim = 92	
	LOAD = True

	grid = Grid(paths[0])		
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	
	vol = grid.hFacC * grid.DRF * grid.RAC; vol = vol[:,:ylim,:]
	totvol = np.sum(vol)
	
	VAR = 'THETA'
	
	#== 
	
	data = None
	for pi in range(len(paths)):
		
		# Load or compute timeseries?
		

		if LOAD:
			if data is None:
				data = np.load(paths[pi]+'OBJtimeseries.npy')
			else:
				data = np.append(data, np.load(paths[pi]+'OBJtimeseries.npy'))

		else:
			tmp = readVariable(VAR, paths[pi], file_format='nc', meta=False)[::d,:,:ylim,:]
			if data is None:
				data = np.sum(vol * tmp, axis=(1,2,3))
				#np.save(paths[pi]+'OBJtimeseries', np.ma.filled(data, fill_value=0)); quit()
			else:
				data = np.append(data, np.sum(vol * tmp, axis=(1,2,3)))

	#==
		
	data /= totvol

	time = np.linspace(0, d*len(data), len(data)) / 12.
	
	
	data = data[30*12:]; time = time[30*12:]
	
	
	title = 'OBJ: on-shelf depth-averaged temperature (deg. C)'
	plt.figure(figsize=(8,5))
	plt.plot(time, data)
	plt.title(title)
	plt.xlabel('Time (years)')
	plt.grid()
	plt.tight_layout()
	plt.show()
	quit()	
	
#==

kinDynSens = 0
if kinDynSens:

	# Load sensitivities to T/S for computing kin. and dyn. sensitivities.

	path = '/home/michael/Documents/data/MCS/run/'
	ndays = 5
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.

	VAR1 = 'ADJtheta'
	VAR2 = 'ADJsalt'
	
	tAlpha = 3.90e-5; sBeta = 7.41e-4
	mult_gencost = 1.e9
 
	nx = 240; ny = 200; nz = 20
	dims = (nz, ny, nx)
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	
	vol = grid.hFacC * grid.DRF * grid.RAC
	
	# Compute kin. and dyn. sensitivities.
	Fdyn = - readAllnp(VAR2, path, dims, reverse=False) * tAlpha / sBeta
	Fkin = readAllnp(VAR1, path, dims, reverse=False) - Fdyn
	
	Fdyn /= mult_gencost; Fkin /= mult_gencost

	# Plot variance timeseries of sens. at each level.
	depthVar = 1
	if depthVar:
		d = 1
		kbot = 20
		Fkinmean = np.mean(Fkin, axis=(0,2,3))[:kbot][::d]
		Fdynmean = np.mean(Fdyn, axis=(0,2,3))[:kbot][::d]
		t = np.linspace(ndays*Fdyn.shape[0],1,Fdyn.shape[0])
		tlag = - t + 90
		
		for k in range(kbot//d):	
			plt.plot(tlag, np.sum((Fkin[:,::d][:,k]-Fkinmean[k])**2,axis=(1,2)), label=str(d*k))
		plt.title('Level-wise kin. sens. variance')
		plt.xlabel('Lag (days)')
		plt.grid()
		plt.legend()
		plt.savefig('LevelKinSensVar')
		plt.show()	

		for k in range(kbot//d):	
			plt.plot(tlag, np.sum((Fdyn[:,::d][:,k]-Fdynmean[k])**2,axis=(1,2)), label=str(d*k))
		plt.title('Level-wise dyn. sens. variance')
		plt.xlabel('Lag (days)')
		plt.grid()
		plt.legend()
		plt.savefig('LevelDynSensVar')
		plt.show()	
		
		quit()
		
	Fdyn_zint = tools.depthIntegral(Fdyn, grid)
	Fkin_zint = tools.depthIntegral(Fkin, grid)
		
	#==
	# PLOT
	#==

	Fkin_zint = ptt.maskBathyXY(Fkin_zint, grid, 0, timeDep=True)
	Fdyn_zint = ptt.maskBathyXY(Fdyn_zint, grid, 0, timeDep=True)
			
	time = ndays*86400*np.linspace(0,Fdyn.shape[0]-1,Fdyn.shape[0])
	text_data = ptt.getTextData(time, 'day', X[-12,1], Y[-12,1], color='k', lag=True, lagTime=90)
	
	#==

	# Depth int. dyn. sens. anim.
	if 1:
		vmin = -1.e-3; vmax = -vmin
		levels = [-580, -501, -420]
		title = 'Depth-integrated dyn. sens. (deg C. / deg C.)'
		outname = 'Fdyn_zint.mp4'
		contourLevels = 9
		
		pt.animate1by1(Fdyn_zint[::-1], X, Y, vmin=vmin, vmax=vmax, cmap='bwr', title=title, text_data=text_data, fontsize=9, contour=grid.bathy, contourLevels=levels, outname=outname)
		
	#==

	# Depth int. kin. sens. anim.
	if 1:
		vmin = 0.4e-5; vmax = 2.8e-5
		levels = [-580, -501, -420]
		title = 'Depth-integrated kin. sens. (deg C. / deg C.)'
		outname = 'Fkin_zint.mp4'
		contourLevels = 9
		
		pt.animate1by1(Fkin_zint[::-1], X, Y, vmin=vmin, vmax=vmax, cmap='viridis', title=title, text_data=text_data, fontsize=9, contour=grid.bathy, contourLevels=levels, outname=outname)
	
	# Fixed-depth dyn. sens. anim.
	if 0:
		k = 12
		Fdyn = ptt.maskBathyXY(Fdyn[:,k], grid, k, timeDep=True)
		vmin = -1.e-5; vmax = -vmin
		levels = [-580, -501, -420]
		title = 'Z = ' +str(grid.RC.squeeze()[k])+ ' m, dyn. sens. (deg C. / deg C.)'
		outname = 'Fdyn_z'+str(grid.RC.squeeze()[k])+'.mp4'
		contourLevels = 9
		
		pt.animate1by1(Fdyn[::-1], X, Y, vmin=vmin, vmax=vmax, cmap='bwr', title=title, text_data=text_data, fontsize=9, contour=grid.bathy, contourLevels=levels, outname=outname)
		
	#==

	# Fixed-depth kin. sens. anim.
	if 0:
		k = 12
		Fkin = ptt.maskBathyXY(Fkin[:,k], grid, k, timeDep=True)
		vmin = 0.4e-6; vmax = 2.8e-6
		levels = [-580, -501, -420]
		title = 'Depth-integrated kin. sens. (deg C. / deg C.)'
		outname = 'Fkin_z'+str(grid.RC.squeeze()[k])+'.mp4'
		contourLevels = 9
		
		pt.animate1by1(Fkin[::-1], X, Y, vmin=vmin, vmax=vmax, cmap='viridis', title=title, text_data=text_data, fontsize=9, contour=grid.bathy, contourLevels=levels, outname=outname)
		
	quit()
	
	
#==


animBin = 1
if animBin:
	
	#==
	
	path = '/home/michael/Documents/data/MCS/run/'
	
	VAR = 'ADJtaux'; mult_gencost = 1.e9; reverse = True; vmax = 1.e0; vmin = -vmax
	title=VAR + ' (deg. C/(N m$^{-2}$)), OBJ=on-shelf heat, days 1710-1800'
	ndays = 5
		
	nx = 240; ny = 200
	dims = (ny, nx)
	
	data = readAllnp(VAR, path, dims, reverse=reverse) / mult_gencost
	print(data.shape)
	nt = data.shape[0]

	#==
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	data = ptt.maskBathyXY(data, grid, 0, timeDep=True)
			
	time = ndays*86400*np.linspace(0, nt-1, nt)
	text_data = ptt.getTextData(time, 'day', X[-12,1], Y[-12,1], color='k', lag=True, lagTime=90)
		
	contour = grid.bathy
	levels = [-600, -501, -401]
	
	var = np.var(data, axis=(1,2))
	dtvar = np.diff(var)
	dt2var = np.diff(dtvar)
	#plt.plot(var); plt.show(); quit()
		
	#pt.plot1by1(data[-30,], X, Y, mesh=True, cmap='bwr', vmin=0.4*vmin, vmax=0.4*vmax, title='HC sens. to zonal wind, 60 day lag (deg. C/(N m$^{-2}$))', xlabel='Lon (km)', ylabel='Lat (km)', fontsize=10, contour=contour); quit()
	
	pt.animate1by1varCbar(data, X, Y, vmin=vmin, vmax=vmax, cmap='bwr', title=title, text_data=text_data, fontsize=9, contour=contour, contourLevels=levels)
	
	quit()

