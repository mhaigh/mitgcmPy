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

LAT = 88

animYZ = 0
if animYZ:
	
	#==
	
	run = 'MCS_308'; xx = [50,120]; xi = 120#301 and 302
	#run = 'MCS_303'; xx = [0,None]; xi=120#304 and 305
	
	path = '/home/michael/Documents/data/'+run+'/run/'
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()
	
	VAR = 'THETA'	
	data = readVariable(VAR, path, file_format='nc', meta=False, xx=xx, tt=[0,None])
	vmin, vmax, cmap, title = getPlottingVars(VAR)	

	#==
	
	data = np.mean(data, axis=-1)
	data = ptt.maskBathyYZ(data, grid, 0, timeDep=True, xi=120)
			
	nt = data.shape[0]; ndays = 30;
	time = ndays*86400*np.linspace(0, nt-1, nt)
	text_data = ptt.getTextData(time, 'month', Y[1], Z[-1], color='k')
		
	hline = grid.YC[92,0]/1000.
	
	data = tools.boundData(data, vmin, vmax)
	
	pt.animate1by1(data, Y, Z, vmin=vmin, vmax=vmax, cmap='coolwarm', text_data=text_data, mesh=False, hline=hline, outname=VAR+'_'+run+'.mp4')
	
	quit()
	
#==

surfBotFlow = 0
if surfBotFlow:

	path_root = '/home/michael/Documents/data/'
	paths = ['MCS_154', 'MCS_301', 'MCS_302']
	
	path = path_root+ paths[0] + '/run/'	
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()

	xx = [50,120]

	VAR = 'UVEL'	
	vmin, vmax, cmap, title = getPlottingVars(VAR)	

	data = [[],[]]
	for path_tmp in paths:
	
		tmp = readVariable(VAR, path_root+path_tmp+'/run/', file_format='nc', meta=False,  tt=[-12,None])
		
		# Time mean surface flow
		data[0].append(np.mean(tmp[:,0,:,xx[0]:xx[1]+1], axis=(0,-1)))
		
		# Time mean bottom flow
		ub = tools.bottom(tmp, grid, 'u', shiftUp=1)
		data[1].append(np.mean(ub[...,xx[0]:xx[1]+1], axis=(0,-1)))
	
	plt.figure()
	plt.subplot(121)
	for di, d in enumerate(data[0]):
		plt.plot(d[1:-1], Y[1:-1], label=paths[di])
	plt.grid()
	plt.subplot(122)
	for di, d in enumerate(data[1]):
		plt.plot(d[1:-1], Y[1:-1], label=paths[di])
	plt.grid()
	plt.legend()
	plt.show()
	
	quit()
	
#==

surfBotFlowAnim = 0
if surfBotFlowAnim:

	path_root = '/home/michael/Documents/data/'

	#paths = ['MCS_154', ['MCS_306', 'MCS_302']]; bathyName = 'bathyS'; 	xx = [50,120]
	paths = ['MCS_303', ['MCS_304', 'MCS_305']]; bathyName = 'bathyUniform'; xx = [0,None]
	labels = ['Southward shift', 'Northward shift']
		
	ttref = [-120, None]
	tt = [0, 200]
		
	# First get reference flow from IC run.
	path = path_root + paths[0] + '/run/'
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()
	
	VAR = 'UVEL'	
	vmin, vmax, cmap, title = getPlottingVars(VAR)	
	tmp = readVariable(VAR, path, file_format='nc', meta=False,  tt=ttref)
	
	# Time mean surface flow
	surfRef = np.mean(tmp[:,0,:,xx[0]:xx[1]], axis=(0,-1))
	
	# Time mean bottom flow
	ub = tools.bottom(tmp, grid, 'u', shiftUp=1)
	botRef = np.mean(ub[...,xx[0]:xx[1]], axis=(0,-1))

	#==
	
	# Now get time-dep. surface and bottom flows from perturbation experiments.
	
	data = [[],[]]
	for path_tmp in paths[1]:
	
		print(path_tmp)
	
		tmp = readVariable(VAR, path_root+path_tmp+'/run/', file_format='nc', meta=False, tt=tt)
		
		# Surface flow
		surftmp = np.mean(tmp[:,0,1:-1,xx[0]:xx[1]], axis=-1)
		surftmp = tools.smooth3(surftmp)
		data[0].append(surftmp)

		# Bottom flow
		bottmp = tools.bottom(tmp, grid, 'u', shiftUp=2)
		bottmp = np.mean(bottmp[...,1:-1,xx[0]:xx[1]], axis=-1)
		bottmp = tools.smooth3(bottmp)
		data[1].append(bottmp)

	#==
	
	nt = data[0][0].shape[0]; ndays = 30;
	time = ndays*86400*np.linspace(1, nt, nt)
	text_data = ptt.getTextData(time, 'month', 0.1, Y[-2], color='k')
	
	timeSeriesPlot = False
	if timeSeriesPlot:
		y0 = 90; y1 = 100
		for d in data[1]:
			plt.plot(time, np.mean(d[:,y0:y1], axis=-1))
			#plt.plot(d[0,90:100])
		plt.plot(time, np.ones(time.shape[0])*np.mean(botRef[y0:y1]))
		plt.show()
		quit()
		
	print('Look at SSH and theta changes, and in fully zonally uniform model')
	print('using timeSeriesPlot, clear that one has greater variance than other. Why?')
	print('Why do both lead to faster undercurrents?')
	
	ylabel = 'Y (km)'
	xlabel = 'Surface flow'
	constLineLabel = ['Default wind']
	constLineStyle = ['dashed']
	
	pt.animateLine(data[0], X=Y[1:-1], transposeAx=True, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[surfRef[1:-1]], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animSurfFlow_'+bathyName+'.mp4')
	
	#==
	
	xlabel = 'Bottom flow'
	
	pt.animateLine(data[1], X=Y[1:-1], transposeAx=True, vmin=-0.1, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[botRef[1:-1]], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animBotFlow_'+bathyName+'.mp4')
	

#==

sshAnim = 0
if sshAnim:

	path_root = '/home/michael/Documents/data/'
	
	#paths = ['MCS_154', ['MCS_301', 'MCS_302']]; bathyName = 'bathyS'; 	xx = [50,120]
	paths = ['PISOMIP_001', ['PISOMIP_001', 'PISOMIP_002']]; bathyName = 'bathyS'; xx = [50,120]
	#paths = ['MCS_303', ['MCS_304', 'MCS_305']]; bathyName = 'bathyUniform'; xx = [0,None]
	labels = ['Southward shift', 'Northward shift']
	
	# First get reference flow from IC run.
	path = path_root + paths[0] + '/run/'
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()
	
	VAR = 'oceTAUX'	
	vmin, vmax, cmap, title = getPlottingVars(VAR)	
	tmp = readVariable(VAR, path, file_format='nc', meta=False,  tt=[-12,None])
	
	# Time mean surface flow
	sshRef = np.mean(tmp[:,:,xx[0]:xx[1]], axis=(0,-1))
	
	#==
	
	# Now get time-dep. surface and bottom flows from perturbation experiments.
	
	paths = paths[1]

	data = [[],[]]
	for path_tmp in paths:
	
		sshtmp = readVariable(VAR, path_root+path_tmp+'/run/', file_format='nc', meta=False)
		
		# Surface flow
		sshtmp = np.mean(sshtmp[:,1:-1,xx[0]:xx[1]], axis=-1)
		sshtmp = tools.smooth3(sshtmp)
		data[0].append(sshtmp)

	#==
	
	nt = data[0][0].shape[0]; ndays = 30;
	time = ndays*86400*np.linspace(1, nt, nt)
	text_data = ptt.getTextData(time, 'month', Y[1], 0.15, color='k')
	
	xlabel = 'Y (km)'
	ylabel = 'SSH'
	constLineLabel = ['Default wind']
	constLineStyle = ['dashed']
	
	pt.animateLine(data[0], X=Y[1:-1], transposeAx=False, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[sshRef[1:-1]], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animSSH_'+bathyName+'.mp4')
	
	quit()
	
	#==
	
#==

isothermAnim = 0
if isothermAnim:

	paths = ['MCS_154', ['MCS_308', 'MCS_302']]; bathyName = 'bathyS'; xx = [50,120]
	#paths = ['MCS_154', ['MCS_306', 'MCS_307']]; bathyName = 'bathyS'; xx = [50,120]
	#paths = ['MCS_303', ['MCS_304', 'MCS_305']]; bathyName = 'bathyUniform'; xx = [0,None]

	labels = ['Southward shift', 'Northward shift']

	path_root = '/home/michael/Documents/data/'
	path_ref = paths[0]
	VAR = 'ThermZ_m05_' 
	
	# First get reference flow from IC run.
	path = path_root + path_ref + '/run/'
	grid = Grid(path)
	X = grid.XC / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()
	
	tmp = np.load(path+VAR+path_ref+'.npy')
	
	# Time mean surface flow
	isothermRef = np.mean(tmp[-120:,1:-1,xx[0]:xx[1]], axis=(0,-1))
	
	#==
	
	# Now get time-dep. surface and bottom flows from perturbation experiments.
	data = []
	for path_tmp in paths[1]:
	
		path = path_root+path_tmp+'/run/'
		tmp = np.load(path+VAR+path_tmp+'.npy')
		print(tmp.shape)
		tmp = np.mean(tmp[:,1:-1,xx[0]:xx[1]], axis=-1)
		#tmp = tools.smooth3(tmp)
		data.append(tmp[0:240])

	#==
	
	nt = data[1].shape[0]; ndays = 30
	time = ndays*86400*np.linspace(1, nt, nt)
	text_data = ptt.getTextData(time, 'month', Y[1], -220, color='k')
	
	xlabel = 'Y (km)'
	ylabel = '-0.5 deg. C isotherm depth'
	constLineLabel = ['Default wind']
	constLineStyle = ['dashed']
	
	vmin = -600; vmax = -150
	
	mask = grid.hFacC[:,:,120]
	mask = np.where(mask>0, 0, 1)
	mask = [Y, Z, mask]
	
	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[isothermRef], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animIsotherm_'+bathyName+'.mp4')
	
#==

isothermPlanAnim = 0
if isothermPlanAnim:

	run = 'MCS_308'

	path_root = '/home/michael/Documents/data/'
	VAR = 'ThermZ_m05_' 
	
	# First get reference flow from IC run.
	path = path_root + run + '/run/'
	grid = Grid(path)
	X = grid.XC[0,:] / 1000.
	Y = grid.YC[:,0] / 1000.
	Z = grid.RC.squeeze()
	
	data = np.load(path+VAR+run+'.npy')
	
	#==
	
	nt = data.shape[0]; ndays = 30
	time = ndays*86400*np.linspace(1, nt, nt)
	text_data = ptt.getTextData(time, 'month', X[1], Y[1], color='k')
	
	xlabel = 'X (km)'
	ylabel = 'Y (km)'
	title = '-0.5 deg. C isotherm depth'
	
	vmin = -600; vmax = -150
	cmap = 'YlOrRd'
	
	pt.animate1by1(data, X, Y, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=True, text_data=text_data, outname='isothermPlan_'+run+'.mp4', vmin=vmin, vmax=vmax, contour=grid.bathy)
	
	quit()

#==

heatContentTimeSeries = 0
if heatContentTimeSeries:

	paths = ['MCS_154', ['MCS_301', 'MCS_307']]; bathyName = 'bathyS'; xx = [50,120]; tt=[0,None]
	#paths = ['MCS_154', ['PISOMIP_001', 'PISOMIP_002']]; bathyName = 'bathyS'; xx = [50,120]; tt=[0,224]
	#paths = ['MCS_303', ['MCS_304', 'MCS_305']]; bathyName = 'bathyUniform'; xx = [0,None]; tt=[0,None]
	labels = ['Southward wind shift', 'Northward wind shift']
 	
	y0 = 90; z0 = 25
	yy = [0, y0]; zz = [0, z0]
 
	path_root = '/home/michael/Documents/data/'
	path = path_root + paths[0] + '/run/'
	grid = Grid(path)
	
	VAR = 'THETA'
	HCref = readVariable(VAR, path_root+paths[0]+'/run/', file_format='nc', meta=False, yy=yy, zz=zz)
	HCref = tools.heatContentShelfFast(HCref)
	Ntref = len(HCref)
	tref = np.linspace(1, Ntref, Ntref)
	
	HCs = []
	for path_tmp in paths[1]:
		print(path_tmp)
		tmp = readVariable(VAR, path_root+path_tmp+'/run/', file_format='nc', meta=False, yy=yy, zz=zz, tt=tt)
		print(tmp.shape)
		HCs.append(tools.heatContentShelfFast(tmp))
		
	Nt = len(HCs[1])
	t = np.linspace(Ntref+1, Ntref+1+Nt, Nt)
	
	plt.plot(tref, HCref, color='k', label='Default wind spin up')
	plt.plot(tref, HCs[0], label=labels[0])
	plt.plot(t, HCs[1], label=labels[1])
	plt.legend()
	plt.xlabel('Time (months)')
	plt.title('On shelf heat content (J), ' + bathyName)
	
	plt.grid()
	plt.show()
	
	quit()
	
#==

barotropicStreamfunction = 0
if barotropicStreamfunction:

	#path = '/home/michael/Documents/data/PAS_666/run/'
	#path = '/home/michael/Documents/data/PISOMIP_003/run/'
	path_root = '/home/michael/Documents/data/'
	
	run = 'MCS_301'
	
	tt = [0, None]
	
	path = path_root + run + '/run/'
	grid = Grid(path)
	bathy = grid.bathy
	
	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	xlabel = 'LON (km)'; ylabel = 'LAT (km)'
	ny, nx = bathy.shape
	
	plt.plot(bathy); plt.grid(); plt.show(); quit()
		
	#VAR = 'ISOTHERM'
	#VAR = 'BAROSTR'
	VAR = 'BOTVEL'

	if VAR == 'ISOTHERM':
	
		data = readVariable('ETAN', path, file_format='nc', tt=tt)

		
		data = np.load(path+'ThermZ_m05_'+run+'.npy')
		vmin = -600; vmax=-150;	cmap = 'YlOrRd'
		#vmin = -500; vmax = -350; cmap = 'jet'
		title = run + ' -0.5 deg. C isotherm depth'
		print(data.shape)

	#==
		
	elif VAR == 'BAROSTR':
			
		data = readVariable('UVEL', path, file_format='nc', tt=tt)
			
		vmin = -0.4e-5; vmax = 0.4e-5; cmap = 'jet'
		data = - 1.e-6 * tools.barotropicStreamfunction(data, grid, timeDep=True, norm=False)
		print(data.shape)
		
		title = 'Barotr. Streamfunction; ' + run
	
	#==
	
	elif VAR == 'BOTVEL':
		
		data = readVariable('UVEL', path, file_format='nc', tt=tt)
		data = tools.bottom(data, grid, 'u', shiftUp=1)
		
		title = 'u bottom; ' + run
				
		vmin = -0.1; vmax = 0.1; cmap = 'jet'
		print(data.shape)

	#==

	nt = data.shape[0]
	time_s = np.linspace(1,nt,nt)*86400.*30.	
	text_data = ptt.getTextData(time_s, 'month', X[1], Y[1], color='k')
	
	data = ptt.maskBathyXY(data, grid, 0, timeDep=True)
	
	pt.animate1by1(data, X, Y, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=True, text_data=text_data, outname='animate1by1Surf.mp4', vmin=vmin, vmax=vmax, contour=grid.bathy)

	quit()

#==

# Cross-shelf heat transport plots.
heatFluxes_heatRel = 1
if heatFluxes_heatRel:

	path_root = '/home/michael/Documents/data/'
	
	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)
	rho0 = 1030. # Units kg / m3
	
	dt_month = 30. * 86400.
	dt_rel = 1. / (100. * 86400.)
	
	# First experiment defined here.
	exp = 'MCS_308'
	
	path = path_root + exp + '/run/'#'MCS_158/run/'
	grid = Grid(path)

	X = grid.XC[1,:]/1000.
	Y = grid.YC[:,1]/1000.
	Z = grid.RC.squeeze()

	# Subregions for heat transport.

	lat = 88

	#v1 = np.mean(readVariable('VVEL', path, meta=False, tt=-1), axis=-1)
	#v2 = np.mean(readVariable('VVEL', path_root + exp3 + '/run/', meta=False, tt=-1), axis=-1)
	#v1 = ptt.maskBathyYZ(v1, grid, xi=120, timeDep=False); v2 = ptt.maskBathyYZ(v2, grid, xi=120, timeDep=False)
	#pt.plot1by2([v1,v2],X=[Y,Y],Y=[Z,Z], vmin=-0.001,vmax=0.001,mesh=True); quit()
	
	#==
	
	# 1. 
	# First get heat transport across shelf break

	T = readVariable('THETA', path, meta=False, yy=lat)
	Tf = -1.8 * np.ones(T.shape)
	v = np.mean(readVariable('VVEL', path, meta=False, yy=[lat,lat+2]), axis=-2)
	
	vT = v[-1]*T[-1]
	vT2 = readVariable('VVELTH', path, meta=False, yy=lat, tt=-1)
	
	pt.plot1by2([vT, vT2-vT]); quit()

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	T = rho0 * Cp * v * (T - Tf)

	Tshelf = np.ma.sum(ptt.maskBathyXZ(area*T, grid, yi=lat, timeDep=True), axis=(1,2)) * dt_month

	# Then get heat removal by relaxation in the south.
	nz = 25
	maskWidth = 4
	mask = tools.getRelMask(grid, maskWidth=maskWidth)
	mask = ptt.maskBathyAll(mask, grid, timeDep=False)
	
	T = readVariable('THETA', path, meta=False, yy=[1,maskWidth+1], zz=[0,nz])
	print(T.shape)
	Tf = -1.8 * np.ones(T.shape)
	vol = grid.volume(ylims=[1,maskWidth], zlims=[0,nz-1], limsIndex=True)
	print(vol.shape)
	Trel = rho0 * Cp * (T - Tf) * vol * mask[0:nz,1:maskWidth+1,:] 
	Trel = np.ma.sum(Trel, axis=(1,2,3)) * dt_rel * dt_month
	
	#==

	Tshelf = np.cumsum(Tshelf); Trel = np.cumsum(Trel)

	plt.plot(Tshelf, label='Tshelf')
	plt.plot(Trel, label='Trel')
	plt.plot(Tshelf-Trel, label='sum')
	plt.legend()
	plt.grid()
	plt.show()

#==

# Cross-shelf heat transport plots.
heatFluxes = 0
if heatFluxes:

	path_root = '/home/michael/Documents/data/'
	
	SMOOTH = True
	TF = True
	WARMFLUXES = True
	
	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)
	rho0 = 1030. # Units kg / m3
	normval = 1.e12
	
	# First experiment defined here.
	exp1 = 'MCS_154'; exp2 = 'MCS_301';	exp3 = 'MCS_302'
	#exp1 = 'MCS_303'; exp2 = 'MCS_304';	exp3 = 'MCS_305'
	
	path = path_root + exp1 + '/run/'#'MCS_158/run/'
	grid = Grid(path)

	# Subregions for heat transport.

	lat = 92

	# troughW
	lonsW = [100e3, 200e3]; depthW = [-10, -800]
	lonW = grid.getIndexFromLon(lonsW)
	depthW = grid.getIndexFromDepth(depthW)
	labelW = 'PITW'

	# troughE
	lonsE = [418e3, 520e3]; depthE = [-10, -800]
	lonE = grid.getIndexFromLon(lonsE)
	depthE = grid.getIndexFromDepth(depthE)
	labelE = 'PITE'

	# troughE AND sill
	lonsES = [418e3, 585e3]; depthES = [-10, -800]
	lonES = grid.getIndexFromLon(lonsES)
	depthES = grid.getIndexFromDepth(depthES)
	labelES = 'PITE + R'

	labelAll = 'All'
	labelU = 'Uniform lons'

	# Define colours for different lon ranges.
	cAll = 'k'
	cE = 'r'
	cU = 'b'
	cES = 'g'
	cW = 'c' 

	#==
	
	# 1. 
	# Get heat transport for first (bathyE) experiment.

	T = readVariable('THETA', path, meta=False, yy=lat)
	Tf = -1.8 * np.ones(T.shape)
	v = np.mean(readVariable('VVEL', path, meta=False, yy=[lat,lat+2]), axis=-2)

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	if TF:
		T = rho0 * Cp * v * (T - Tf)
	else:
		if WARMFLUXES:
			T = rho0 * Cp * v * np.where(T>0, T, 0)
		else:
			T = -rho0 * Cp * v * np.where(T<0, T, 0)
		
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)
	T *= area / normval
	
	# Heat transport for troughE, uniform lons and all lons.
	TES = np.ma.sum(T[:, depthES[0]:depthES[1], lonES[0]:lonES[1]], axis=(1,2))
	TAll = np.ma.sum(T, axis=(1,2))
	TU = TAll - TES
	
	if SMOOTH:
		TES = tools.smooth3(TES)
		TAll = tools.smooth3(TAll)
		TU = tools.smooth3(TU)

	av1 = np.mean(TAll[-120:])
	print(exp1 + ': ' + str(av1))
	
	Ts1 = [TES, TU, TAll]
	labels1 = [labelES, labelU, labelAll]
	colours1 = [cES, cU, cAll]
	title1 = '(a) Default wind'

	#==

	# 2. 
	# Repeat for bathyES experiment.
	path = path_root + exp2 + '/run/'
	grid = Grid(path)
	
	T = readVariable('THETA', path, meta=False, yy=lat)
	Tf = -1.8 * np.ones(T.shape)
	v = np.mean(readVariable('VVEL', path, meta=False, yy=[lat,lat+2]), axis=-2)

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	if TF:
		T = rho0 * Cp * v * (T - Tf)
	else:
		if WARMFLUXES:
			T = rho0 * Cp * v * np.where(T>0, T, 0)
		else:
			T = - rho0 * Cp * v * np.where(T<0, T, 0)
		
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)
	T *= area / normval

	# Heat transport for troughE, uniform lons and all lons.
	TES = np.ma.sum(T[:, depthES[0]:depthES[1], lonES[0]:lonES[1]], axis=(1,2))
	TAll = np.ma.sum(T, axis=(1,2))
	TU = TAll - TES
	
	plt.plot(30*86400.*normval*np.cumsum(TAll)); plt.show(); quit()
	
	if SMOOTH:
		TES = tools.smooth3(TES)
		TAll = tools.smooth3(TAll)
		TU = tools.smooth3(TU)

	av2 = np.mean(TAll[-120:])
	print(exp2 + ': ' + str(av2))
	
	Ts2 = [TES, TU, TAll]
	labels2 = [labelES, labelU, labelAll]
	colours2 = [cES, cU, cAll]
	title2 = '(b) Southward shift'
	
	#==
	
	# 3. 
	# Repeat for bathyWCES2 experiment.
	path = path_root + exp3 + '/run/'
	grid = Grid(path)
	
	T = readVariable('THETA', path, meta=False, yy=lat)
	Tf = -1.8 * np.ones(T.shape)
	v = np.mean(readVariable('VVEL', path, meta=False, yy=[lat,lat+2]), axis=-2)

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]
	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)

	if TF:
		T = rho0 * Cp * v * (T - Tf)
	else:
		if WARMFLUXES:
			T = rho0 * Cp * v * np.where(T>0, T, 0)
		else:
			T = - rho0 * Cp * v * np.where(T<0, T, 0)
		
	T = ptt.maskBathyXZ(T, grid, yi=lat, timeDep=True)
	T *= area / normval
	
	# Heat transport for troughE, uniform lons and all lons.
	TES = np.ma.sum(T[:, depthES[0]:depthES[1], lonES[0]:lonES[1]], axis=(1,2))
	TAll = np.ma.sum(T, axis=(1,2))
	TU = TAll - TES

	if SMOOTH:
		TES = tools.smooth3(TES)
		TAll = tools.smooth3(TAll)
		TU = tools.smooth3(TU)

	av3 = np.mean(TAll[-120:])
	print(exp3 + ': ' + str(av3))
	
	Ts3 = [TES, TU, TAll]
	labels3 = [labelES, labelU, labelAll]
	colours3 = [cES, cU, cAll]
	title3 = '(c) Northward shift'
	
	#==

	# Prepare plotting data.
	
	Ts = [Ts1, Ts2, Ts3]
	labels = [labels1, labels2, labels3]
	colours = [colours1, colours2, colours3]
	titles = [title1, title2, title3]
	
	ylabel = 'Heat transport (TW)'
	ylabels = [ylabel, ylabel, ylabel]
	
	xlabel = 'Time (months)'
	xlabels = [None, None, xlabel]
	
	xlims = [0, 360]
	if WARMFLUXES:
		ylims = [-.3, .1]
	else:
		#ylims = [-10,10]	
		ylims = [-3,3]
		
	xticks = np.linspace(0, 360, 10)
	xticks = [xticks]*3
	xticksvis = [False, False, True]
	
	yticks = np.linspace(ylims[0], ylims[1], 5)
	yticks = [yticks]*3
	yticksvis = [True, True, True]
	
	fdict = {'fontsize':12, 'color':'k'}
	text = ['{:.3f}'.format(av1), '{:.3f}'.format(av2), '{:.3f}'.format(av3)]
	text = ['Heat transport "all" average = ' + t + ' TW' for t in text]
	xloc = 100; yloc = 0.8*ylims[1]
	text_data = {'text':text, 'xloc':xloc, 'yloc':yloc, 'fontdict':fdict}
	
	#==

    # NOW PLOT
	pt.line1byN(Ts, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis, text_data=text_data, save=True, outname='heatFluxes.png')

#==

	
	
	
