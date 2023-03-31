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

# Meridional limits of mask:
# YC>=2.21250e5

objTest = 1
if objTest:

	path1 = '/home/michael/Documents/data/MCS/run/'
	path2 = '/home/michael/Documents/data/MCS_ad2/run/'
	
	VAR1 = 'stateTheta.0000008640.data'; norm1 = 1.e9
	VAR2 = 'm_boxmean_theta.0000000000.data'; norm2 = 1.e0

	ny = 200; nx = 240
	dims = (ny, nx)
	LAT = 92
	
	grid = Grid(path1)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	
	vol = grid.hFacC * grid.DRF * grid.RAC
	
	data1 = readnp(path1+VAR1, dims, rec=-1)# / norm1
	print(data1.shape)
	
	data1[:,LAT:] = 0
	vol[:,LAT:] = 0
	
	# Integrate over volume and normalise.
	data1 = np.sum(vol * data1, axis=0)	
	
	totvol = np.sum(vol)
	data1 = data1 / totvol	
	
	sum1 = np.sum(data1)#*norm1

	data2 = readnp(path2+VAR2, dims, rec=-1)
	print(data2.shape)
	
	# This gives the output in STDOUT (boxmean#horflux) before scaling by time mask.
	print(np.sum(data2,axis=(1,2)))
	
	data2 = np.mean(data2, axis=0)
	sum2 = np.sum(data2)#*norm2
	
	print(sum1)
	print(sum2)
	print(sum1/sum2)
	
	pt.plot1by2([data1, data2])#, vmin=-2.e-5, vmax=0.)
	
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

animBin = 0
if animBin:
	
	#==
	
	path = '/home/michael/Documents/data/MCS/run/'
	
	VAR = 'ADJtaux'; mult_gencost = 1.e9; reverse = True; vmax = 0.0002; vmin = -vmax
	title = VAR + ' (deg. C/(N m$^{-2}$)), OBJ=on-shelf heat, days 290-360'
	ndays = 5
	
	INTEGRATE = True
	dt = ndays * 86400.
	
	NORM = True
	norm = 1./0.025
	
	nx = 240; ny = 200
	dims = (ny, nx)
	
	data = readAllnp(VAR, path, dims, reverse=reverse) / mult_gencost
	print(data.shape)
	nt = data.shape[0]
	
	if NORM:
		data /= norm
		title = VAR + ' (nondim), OBJ=on-shelf heat, days 290-360'
	#==
	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC / 1000.
	data = ptt.maskBathyXY(data, grid, 0, timeDep=True)
			
	time = ndays*86400*np.linspace(0, nt-1, nt)
	text_data = ptt.getTextData(time, 'day', X[-12,1], Y[-12,1], color='k', lag=True, lagTime=90)
		
	contour = grid.bathy
	levels = [-600, -501, -401]
	
	# PLOT VARIANCE
	#var = np.var(data, axis=(1,2))
	#dtvar = np.diff(var)
	#dt2var = np.diff(dtvar)
	#plt.plot(var); plt.show(); quit()
		
	if INTEGRATE:
		data = dt*np.cumsum(data,axis=0)
		title = VAR + ' int. (nondim), OBJ=on-shelf heat, days 290-360'
				
	pt.animate1by1varCbar(data, X, Y, cmap='bwr', title=title, text_data=text_data, fontsize=9, contour=contour, vmin=vmin, vmax=vmax)
	
	#pt.animate1by1(data, X, Y, cmap='bwr', title=title, text_data=text_data, fontsize=9, contour=contour, vmin=vmin, vmax=vmax, mesh=True)
	
	quit()
	
#==

adjTimings = 0
if adjTimings:

	times = [
	[(48,1,1), 71.218007087707520],\
	[(48,2,1), 129.52633404731750],\
	[(48,2,2), 272.00419402122498],\
	[(48,8,8), 4790.6619710922241],\
	[(48,16,8), 10857.532526016235],\
	[(48,16,16), 23602.149661064148],\
	[(48,32,16), 54387.048799991608],\
	[(48,32,32), 154051.76470588235]\
	]
	
	scale = 1
	Nts = []
	runtimes = []
	for time in times:
		
		chklevs = time[0]
		Nts.append(chklevs[0] * chklevs[1] * chklevs[2])
		runtimes.append(time[1]/scale)
		
	refx = [2**n for n in range(6,18)]
	refy = [2**n/scale for n in range(6,18)]	
	
	fig, ax = plt.subplots()
	plt.plot(Nts, runtimes)
	plt.plot(refx,refy,linestyle='dashed', color='k')
	plt.grid()
	ax.set_xscale('symlog', base=2)
	ax.set_yscale('symlog', base=2)
	plt.xlabel('n timesteps')
	plt.ylabel('Adjoint runtime (s)')

	plt.show()







