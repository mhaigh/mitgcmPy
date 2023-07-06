import numpy as np



from grid import Grid

from grid_PAS import Grid as Grid_PAS



import matplotlib.pyplot as plt

import matplotlib.colors as cl

import matplotlib.patches as patches



import plotting as pt

import plotting_tools as ptt



import tools



from readData import *

from varDict import getPlottingVars, getTrefSref, getTitleData



import time



#==========================================================



LAT = 88

path_root = '/home/michael/Documents/data/'

	

animateQuivers = False

if animateQuivers:



	exp = 'MCS_314'

	

	path = path_root + exp + '/run/'

	titleData = getTitleData(exp)



	# Sample rate

	d = 8



	grid = Grid(path)



	level = grid.getIndexFromDepth(-450)



	X = grid.XC[1,:]/1000.

	Y = grid.YC[:,1]/1000.

	Z = grid.RC.squeeze()



	td = 20



	u = readVariable('UVEL', path, file_format='nc', meta=True)

	TIME = u['TIME'][:]

	utop = u['UVEL'][::td,0]

	ubot = tools.bottom(u['UVEL'][::td], grid, 'u')



	vbot = readVariable('VVEL', path, file_format='nc', meta=False)[::td]

	vtop = vbot[:,0]

	vbot = tools.bottom(vbot, grid, 'v')



	tscale = 86400. * 30

	tscale_name = 'month'

	t = TIME / tscale

	text = [tscale_name + ' ' + str(int(tt)) for tt in t]

	text_data = {'text':text, 'xloc':X[1], 'yloc':Y[1], 'fontdict':{'fontsize':14, 'color':'k'}}



	contour = grid.bathy

	contour = ptt.maskBathyXY(contour, grid, 0, timeDep=False)

	

	utop = ptt.maskBathyXY(utop, grid, level, timeDep=True)

	ubot = ptt.maskBathyXY(ubot, grid, level, timeDep=True)

	vtop = ptt.maskBathyXY(vtop, grid, level, timeDep=True)

	vbot = ptt.maskBathyXY(vbot, grid, level, timeDep=True)

	

	title = '(u, v) & bathymetry; ' + titleData



	udata = [utop, ubot]

	vdata = [vtop, vbot]



	pt.animate1by1quivers(udata, vdata, X, Y, d=d, contour=contour, contourf=False, text_data=text_data, title=title, figsize=(5,4), cmap='YlGn', scale=2., labels=['surf', 'bot'])

	

	quit()



#==

animYZ = 0

if animYZ:

	

	#==

	

	run = 'MCS_313'; xx = [50,190]; xi = 120#301 and 302

	#run = 'MCS_303'; xx = [0,None]; xi=120#304 and 305

	

	ylims = [88,118]

	

	path = '/home/michael/Documents/data/'+run+'/run/'

	grid = Grid(path)

	X = grid.XC / 1000.

	Y = grid.YC[:,0] / 1000.

	Z = grid.RC.squeeze()

	

	VAR = 'UVEL'	

	data = readVariable(VAR, path, file_format='nc', meta=False, xx=xx, tt=[0,None])

	vmin, vmax, cmap, title = getPlottingVars(VAR)	



	#==

	

	data = np.mean(data, axis=-1)

	#data = data[..., xi]

	data = ptt.maskBathyYZ(data, grid, 0, timeDep=True, xi=120)

			

	nt = data.shape[0]; ndays = 30;

	time = ndays*86400*np.linspace(0, nt-1, nt)

	text_data = ptt.getTextData(time, 'month', Y[1], Z[-1], color='k')

		

	#hlines = [grid.YC[92,0]/1000.]

	hlines = [Y[ylims[0]], Y[ylims[1]]] 

	

	data = tools.boundData(data, vmin, vmax)

	

	#if ylims is not None:

	#	data = data[...,ylims[0]:ylims[1]]

	#	Y = Y[ylims[0]:ylims[1]]

		

	pt.animate1by1(data, Y, Z, vmin=vmin, vmax=vmax, cmap='coolwarm', text_data=text_data, mesh=False, hlines=hlines, outname=VAR+'_'+run+'.mp4')

	

	quit()

	

#==



topMinusBotFlow = 0

if topMinusBotFlow:





	#paths = ['MCS_154', ['MCS_306', 'MCS_302']]; bathyName = 'bathyS'; xx = [50,190]

	paths = ['MCS_308', ['MCS_312', 'MCS_313']]; bathyName = 'bathyS'; xx = [0,None]

	

	#labels = ['Southward shift', 'Northward shift']

	labels = ['IPO pos', 'IPO neg']

		

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



crossCorr = 0

if crossCorr:



	run = 'MCS_320'; xi = 120

	

	zz = [20,30]

	yy = [84,105]

	xx = [1,190] 

	tt = [0, None]

	

	path = '/home/michael/Documents/data/'+run+'/run/'

	grid = Grid(path)

	bathy = grid.bathy

	X = grid.XC[0,:] / 1000.

	Y = grid.YC[:,0] / 1000.

	Z = grid.RC.squeeze()



	

	VAR = 'UVEL'	

	vmin, vmax, cmap, title = getPlottingVars(VAR)	

	# tmp = u(t, z, y, x)

	

	for ti in range(0,0):

		rect1 = patches.Rectangle((X[xx[0]], Y[yy[0]]), X[xx[1]]-X[xx[0]], Y[yy[1]]-Y[yy[0]], linewidth=1, edgecolor='k', facecolor='none', linestyle='dashed')

		rect2 = patches.Rectangle((X[xx[0]], Y[yy[0]]), X[xx[1]]-X[xx[0]], Y[yy[1]]-Y[yy[0]], linewidth=1, edgecolor='k', facecolor='none', linestyle='dashed')

		tmp = readVariable(VAR, path, file_format='nc', meta=False, tt=-ti)

		tmp1 = ptt.maskBathyXY(tmp[0], grid, zi=0)

		tmp2 = ptt.maskBathyXY(tmp[24], grid, zi=0)

		pt.plot1by2([tmp1, tmp2], X=[X,X], Y=[Y,Y], patches=[rect1,rect2], contour=[bathy,bathy])

		plt.clf()

	

	tmp = readVariable(VAR, path, file_format='nc', meta=False, tt=tt, yy=yy, xx=xx)

	

 	# Get surface flow in x-range

	# u0 = u0(t, x)

	u0 = tmp.copy()[:,0,:,:]

	u0 = np.min(u0, axis=1) # Max over y

	u0_av = np.mean(u0, axis=1) # Mean over x

	print(u0.shape, u0_av.shape)

	

	# Repeat for bottom flow

	#ub = tools.bottom(tmp.copy(), grid, 'u', shiftUp=1)[xx[0]:xx[1]]

	ub = tmp.copy()[:,zz[0]:zz[1],:,:]

	ub = np.max(ub, axis=(1,2)) # Max over z, y

	ub_av = np.mean(ub, axis=1) # Mean over x

	print(ub.shape, ub_av.shape)

	

	#==

	

	plt.plot(u0_av)

	plt.plot(ub_av)

	plt.show()	

	

	quit()

	

#==



BAROCL = 0

if BAROCL:



	#paths = ['MCS_154', ['MCS_306', 'MCS_302']]; bathyName = 'bathyS'; 	xx = [50,190]

	#paths = ['MCS_154', ['MCS_312', 'MCS_313']]; labels = ['IPO PH neg', 'IPO PH pos']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, None]; ylims = [95,110]

	#paths = ['MCS_154', ['MCS_315', 'MCS_314']]; labels = ['IPO unif neg', 'IPO unif pos']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, 154]; ylims = [95,110]

	#paths = ['MCS_154', ['MCS_319', 'MCS_318']]; labels = ['IPO sin neg, no rel', 'IPO sin pos, no rel']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, 205]; ylims = [95,110]

	#paths = ['MCS_154', ['MCS_320', 'MCS_316']]; labels = ['IPO sin exf', 'IPO sin pos']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, 208]; ylims = [95,110]

	#paths = ['MCS_154', ['MCS_321', 'MCS_316']]; labels = ['IPO sin pos sharp', 'IPO sin pos']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, 205]; ylims = [95,110]

	#paths = ['MCS_152', ['MCS_323', 'MCS_322']]; labels = ['IPO sin neg', 'IPO sin pos']; bathyName = 'bathy uniform'; xx = [0,-2];  tt = [0, 240]; ylims = [95,110]

	#paths = ['MCS_152', ['MCS_334', 'MCS_327']]; labels = ['IPO sin neg sharp', 'IPO sin pos sharp']; bathyName = 'bathy uniform'; xx = [0,-2];  tt = [0, 240]; ylims = [95,110]

	#paths = ['MCS_152', ['MCS_326', 'MCS_325']]; labels = ['IPO sin neg sharp (strong)', 'IPO sin pos sharp (strong)']; bathyName = 'bathy uniform'; xx = [0,-2];  tt = [0, 240]; ylims = [95,110]

	paths = ['MCS_340', ['PISOMIP_003', 'PISOMIP_004']]; labels = ['IPO sin neg sharp', 'IPO sin pos sharp']; bathyName = 'bathy uniform'; xx = [2,-2]; tt = [0, 240]; ylims = [90, 100]

		

	#labels = ['Southward shift', 'Northward shift']

	

	SMOOTH = False

	

	ttref = [0, None]



	# First get reference flow from IC run.

	path = path_root + paths[0] + '/run/'

	grid = Grid(path)

	X = grid.XC / 1000.

	Y = grid.YC[:,0] / 1000.

	Z = grid.RC.squeeze()

	

	#pt.plot1by1(grid.bathy); quit()

	

	VAR = 'UVEL'	

	vmin, vmax, cmap, title = getPlottingVars(VAR)	

	tmp = readVariable(VAR, path, file_format='nc', meta=False, tt=ttref)

	

	# Time mean surface flow

	surfRef = np.mean(tmp[:,0,:,xx[0]:xx[1]+1], axis=(0,-1))

	surfRefMean = np.mean(tmp[:,0,ylims[0]:ylims[1]+1,xx[0]:xx[1]+1], axis=(1,2))

	# Area mean surface flow

	

	# Time mean bottom flow

	ub = tools.bottom(tmp, grid, 'u', shiftUp=1)

	botRef = np.mean(ub[...,xx[0]:xx[1]+1], axis=(0,-1))

	

	# Area mean bottom flow

	botRefMean = np.mean(ub[:,ylims[0]:ylims[1]+1,xx[0]:xx[1]+1], axis=(-1,-2))

	

	baroclRef = np.mean(ub[:,ylims[0]:ylims[1]+1,:]-tmp[:,0,ylims[0]:ylims[1]+1,:], axis=(1,2))

	if SMOOTH:

		baroclRef = tools.smooth3(baroclRef)

	Ntref = len(baroclRef)

	tref = np.linspace(1, Ntref, Ntref)/12

	#plt.plot(baroclRef); plt.show(); print(baroclRef.shape); quit()



	#==

	

	# Now get time-dep. surface and bottom flows from perturbation experiments.

	

	wind = None

	data = [[],[]]

	barocl = []

	surf = []

	bot = []

	for path_tmp in paths[1]:

	

		print(path_tmp)

		

		if path_tmp == 'MCS_320' or path_tmp == 'MCS_312':# or path_tmp == 'MCS_326':

			wind = 2*readVariable('oceTAUX', path_root+path_tmp+'/run/', file_format='nc', meta=False, tt=tt)[:,100,0]

	

		tmp = readVariable(VAR, path_root+path_tmp+'/run/', file_format='nc', meta=False, tt=tt)

		#tmp = readAllnp('stateUvel', path_root+path_tmp+'/run/', (50, 200, 240))[:tt[1]]

			

		#tau = readVariable('oceTAUX', path_root+path_tmp+'/run/', file_format='nc', meta=False, tt=-1)

		#plt.plot(tau[:,0], Y); plt.grid(); plt.show()

		

		# Surface flow

		surftmp = np.mean(tmp[:,0,:,xx[0]:xx[1]], axis=-1)

		if SMOOTH:

			surftmp = tools.smooth3(surftmp)

		data[0].append(surftmp[:,1:-1])

		surf.append(np.mean(surftmp[:,ylims[0]:ylims[1]+1], axis=1))

		print(surftmp.shape)



		# Bottom flow

		bottmp = tools.bottom(tmp, grid, 'u', shiftUp=1)

		bottmp = np.mean(bottmp[...,xx[0]:xx[1]], axis=-1)

		if SMOOTH:

			bottmp = tools.smooth3(bottmp)

		data[1].append(bottmp[:,1:-1])

		bot.append(np.mean(bottmp[:,ylims[0]:ylims[1]+1], axis=1))

		

		barocl.append(np.mean(bottmp[:,ylims[0]:ylims[1]+1] - surftmp[:,ylims[0]:ylims[1]+1], axis=1))



	#==

	

	nt = data[0][0].shape[0]; ndays = 30;

	time = ndays*86400*np.linspace(1, nt, nt)

	text_data = ptt.getTextData(time, 'month', 0.1, Y[-2], color='k')

	

	#==

	

	Nt = len(barocl[1])

	t = np.linspace(Ntref+1, Ntref+1+Nt, Nt)/12

	

	plt.figure(figsize=(15,6))

	plt.subplot(121)

	plt.plot(tref, surfRefMean, color='k', label='Spin up surf. flow')

	plt.plot(t, surf[0], label=labels[0] + ', surf. flow', color='b')

	plt.plot(t, surf[1], label=labels[1] + ', surf. flow', color='r')

	plt.plot(tref, botRefMean, color='k', label='Spin up bot. flow', alpha=0.6)

	plt.plot(t, bot[0], label=labels[0] + ', bot. flow', color='b', alpha=0.6)

	plt.plot(t, bot[1], label=labels[1] + ', bot. flow', color='r', alpha=0.6)

	if wind is not None:

		plt.plot(t, wind, linestyle='dashed', color='k', alpha=0.8)

	plt.grid()

	plt.legend()

	plt.xlabel('Time (years)')

	plt.title('Surf. and bot. flow (m/s), ' + bathyName)

	

	plt.subplot(122)

	plt.plot(tref, baroclRef, color='k', label='Default wind spin up')

	plt.plot(t, barocl[0], label=labels[0], color='b')

	plt.plot(t, barocl[1], label=labels[1], color='r')

	if wind is not None:

		plt.plot(t, wind, linestyle='dashed', color='k', alpha=0.8)

	plt.legend()

	plt.xlabel('Time (years)')

	plt.title('Baroclinicity: bot. - surf. flow (m/s), ' + bathyName)

	

	plt.grid()

	plt.show()

	quit()

	

	#==



	ylabel = 'Y (km)'

	xlabel = 'Surface flow'

	constLineLabel = ['Default wind']

	constLineStyle = ['dashed']

	

	pt.animateLine(data[0], X=Y[1:-1], transposeAx=True, vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[surfRef[1:-1]], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animSurfFlow_'+bathyName+'.mp4')

	

	#==

	

	xlabel = 'Bottom flow'

	

	pt.animateLine(data[1], X=Y[1:-1], transposeAx=True, vmin=-0.1, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[botRef[1:-1]], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animBotFlow_'+bathyName+'.mp4')

	

	quit()



#==



surfBotFlowPlot = 1

if surfBotFlowPlot:



	path_root = '/home/michael/Documents/data/'



	#paths = ['MCS_313', 'MCS_154', 'MCS_312']; labels = ['IPO PH pos', 'Ref. wind', 'IPO PH neg']; tt = [-120, None]

	#paths = ['MCS_314', 'MCS_154', 'MCS_315'];  labels = ['IPO unif pos', 'Ref. wind', 'IPO unif neg'];	tt = [40, 160]

	#paths = ['MCS_316', 'MCS_154', 'MCS_317'];  labels = ['IPO sin pos', 'Ref. wind', 'IPO sin neg'];	tt = [-120, None]

	#paths = ['MCS_322', 'MCS_154', 'MCS_323'];  labels = ['IPO sin pos', 'Ref. wind', 'IPO sin neg'];	tt = [-120, None]

	paths = ['PISOMIP_001', 'MCS_340', 'PISOMIP_004'];  labels = ['IPO sin pos', 'Ref. wind', 'IPO sin neg'];	tt = [-120, None]

	ttref = [-120,None]

	

	xlims = [20, -20]

	ylims = [90,100]



	surfs = []

	bots = []

	stresses = []

	

	# For each simulation, get surface, bottom flow and difference.

	for exp in paths:

		

		path = path_root + exp + '/run/'

		grid = Grid(path)

		

		if exp == 'MCS_154':

			tt0 = ttref

		else:

			tt0 = tt

		

		# First get surf. and bottom flow		

		VAR = 'UVEL'	

		tmp = np.mean(readVariable(VAR, path, file_format='nc', meta=False, tt=tt0), axis=0)

		tmp_surf = ptt.maskBathyXY(tmp[0], grid, zi=0, timeDep=False)

		surfs.append(tmp_surf)		

	

		# Time mean bottom flow

		tmp_bot = tools.bottom(tmp, grid, 'u', shiftUp=1, timeDep=False)

		tmp_bot = ptt.maskBathyXY(tmp_bot, grid, zi=0, timeDep=False)

		bots.append(tmp_bot)

		

		# Second get surface stress.

		if exp == 'MCS_154':

			tmp = tools.get05cosWind(grid.Nx, grid.Ny)[:,0]

			tmp = 200 * tmp / np.max(np.abs(tmp)) + 300.

			stresses.append(tmp)

		else:

			VAR = 'oceTAUX'

			tmp = readVariable(VAR, path, file_format='nc', meta=False, xx=0, tt=0)

			tmp = 200 * tmp / np.max(np.abs(tmp)) + 300.

			stresses.append(tmp)

			

	#==

	

	# PLOT



	X = grid.XC[0,:] / 1000.

	Y = grid.YC[:,0] / 1000.

	Z = grid.RC.squeeze()

	

	rect = patches.Rectangle((X[xlims[0]], Y[ylims[0]]), X[xlims[1]]-X[xlims[0]], Y[ylims[1]]-Y[ylims[0]], linewidth=1, edgecolor='k', facecolor='none', linestyle='dashed')

	rect2 = patches.Rectangle((X[xlims[0]], Y[ylims[0]]), X[xlims[1]]-X[xlims[0]], Y[ylims[1]]-Y[ylims[0]], linewidth=1, edgecolor='k', facecolor='none', linestyle='dashed')

	rect3 = patches.Rectangle((X[xlims[0]], Y[ylims[0]]), X[xlims[1]]-X[xlims[0]], Y[ylims[1]]-Y[ylims[0]], linewidth=1, edgecolor='k', facecolor='none', linestyle='dashed')

		

	xticks = [0,100,200,300,400,500,600]

	yticks = [0,100,200,300,400,500]

	vmin = -0.2; vmax = 0.2

	bathy = ptt.maskBathyXY(grid.bathy, grid, zi=0, timeDep=False)	

	

	plt.figure(figsize=(21,15))

	for pi, path in enumerate(paths):

		

		plt.subplot(3,3,pi*3+1)

		plt.gca().patch.set_color('.5')

		plt.pcolormesh(X, Y, surfs[pi], cmap='bwr', vmin=vmin, vmax=vmax)

		plt.grid(); plt.xticks(xticks, labels=''); plt.yticks(yticks, labels='')

		plt.plot(stresses[pi], Y); plt.axvline(300)

		plt.ylabel(labels[pi])

		plt.colorbar()

		plt.contour(X, Y, bathy, levels=[-501, -401], colors='k', linestyles='solid')



		if pi == 0:

			plt.title('Surface flow')

			plt.gca().add_patch(rect)

			

		ax=plt.subplot(3,3,pi*3+2)

		plt.gca().patch.set_color('.5')

		plt.pcolormesh(X, Y, bots[pi], cmap='bwr', vmin=vmin, vmax=vmax)

		plt.grid(); plt.xticks(xticks, labels=''); plt.yticks(yticks, labels='')

		plt.colorbar()

		plt.contour(X, Y, bathy, levels=[-501, -401], colors='k', linestyles='solid')

		if pi == 0:

			plt.title('Bottom flow')

			ax.add_patch(rect2)

			

		plt.subplot(3,3,pi*3+3)

		plt.gca().patch.set_color('.5')

		plt.pcolormesh(X, Y, surfs[pi]-bots[pi], cmap='bwr', vmin=vmin, vmax=vmax)

		plt.grid(); plt.xticks(xticks, labels=''); plt.yticks(yticks, labels='')

		plt.colorbar()

		plt.contour(X, Y, bathy, levels=[-501, -401], colors='k', linestyles='solid')

		if pi == 0:

			plt.title('Surface flow - bottom flow')		

			plt.gca().add_patch(rect3)

				

	plt.tight_layout()

	plt.show()

	

	quit()

	

#==



sshAnim = 0

if sshAnim:



	path_root = '/home/michael/Documents/data/'

	

	#paths = ['MCS_154', ['MCS_301', 'MCS_302']]; bathyName = 'bathyS'; 	xx = [50,190]

	paths = ['MCS_308', ['MCS_309', 'MCS_310']]; bathyName = 'bathyS'; xx = [50,190]

	#paths = ['MCS_303', ['MCS_304', 'MCS_305']]; bathyName = 'bathyUniform'; xx = [0,None]

	labels = ['Southward shift', 'Northward shift']

	

	# First get reference flow from IC run.

	path = path_root + paths[0] + '/run/'

	grid = Grid(path)

	X = grid.XC / 1000.

	Y = grid.YC[:,0] / 1000.

	Z = grid.RC.squeeze()

	

	#VAR = 'oceTAUX'	

	VAR = 'ETAN'

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



	#paths = ['MCS_154', ['MCS_315', 'MCS_314']]; labels = ['IPO unif neg', 'IPO unif pos']; bathyName = 'bathyS'; xx = [50,190]; tt=160

	paths = ['MCS_154', ['MCS_318', 'MCS_319']]; labels = ['IPO sin neg', 'IPO sin pos']; bathyName = 'bathyS'; xx = [50,190]; tt=205

	#paths = ['MCS_154', ['MCS_312', 'MCS_313']]; labels = ['IPO PH neg', 'IPO PH pos']; bathyName = 'bathyS'; xx = [50,190]; tt = 240



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

		data.append(tmp[0:tt])



	#==

	

	nt = data[1].shape[0]; ndays = 30

	time = ndays*86400*np.linspace(1, nt, nt)

	text_data = ptt.getTextData(time, 'month', Y[1], -220, color='k')

	

	xlabel = 'Y (km)'

	ylabel = '-0.5 deg. C isotherm depth'

	constLineLabel = ['Default wind']

	constLineStyle = ['dashed']

	

	vmin = -500; vmax = -150

	

	mask = grid.hFacC[:,:,120]

	mask = np.where(mask>0, 0, 1)

	mask = [Y, Z, mask]

	

	pt.animateLine(data, X=Y[1:-1], vmin=vmin, vmax=vmax, xlabel=xlabel, ylabel=ylabel, labels=labels, text_data=text_data, constLine=[isothermRef], constLineLabel=constLineLabel, constLineStyle=constLineStyle, outname='animIsotherm_'+bathyName+'.mp4')

	

#==



isothermPlanAnim = 0

if isothermPlanAnim:



	run = 'MCS_313'



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

	

	vmin = -500; vmax = -200

	cmap = 'YlOrRd'

	

	data = tools.boundData(data, vmin, vmax)

	

	pt.animate1by1(data, X, Y, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=False, text_data=text_data, outname='isothermPlan_'+run+'.mp4', vmin=vmin, vmax=vmax, contour=grid.bathy)

	

	quit()



#==



heatContentTimeSeries = 0

if heatContentTimeSeries:



	paths = ['MCS_308', ['MCS_312', 'MCS_313']]; bathyName = 'bathyS'; tt=[0,None]

	#paths = ['MCS_154', ['PISOMIP_001', 'PISOMIP_002']]; bathyName = 'bathyS'; tt=[0,224]

	#paths = ['MCS_303', ['MCS_304', 'MCS_305']]; bathyName = 'bathyUniform'; tt=[0,None]



	#labels = ['Southward wind shift', 'Northward wind shift']

	labels = ['IPO pos', 'IPO neg']	

 	

	y0 = 85; z0 = 25

	yy = [0, y0]; zz = [0, z0]

 

	path_root = '/home/michael/Documents/data/'

	path = path_root + paths[0] + '/run/'

	grid = Grid(path)

	

	VAR = 'THETA'

	HCref = readVariable(VAR, path_root+paths[0]+'/run/', file_format='nc', meta=False, yy=yy, zz=zz)

	HCref = tools.heatContentShelfFast(HCref)

	Ntref = len(HCref)

	tref = np.linspace(1, Ntref, Ntref) / 12

	

	HCs = []

	for path_tmp in paths[1]:

		print(path_tmp)

		tmp = readVariable(VAR, path_root+path_tmp+'/run/', file_format='nc', meta=False, yy=yy, zz=zz, tt=tt)

		print(tmp.shape)

		HCs.append(tools.heatContentShelfFast(tmp))

		

	Nt = len(HCs[1])

	t = np.linspace(Ntref+1, Ntref+1+Nt, Nt) / 12

	

	plt.plot(tref, HCref, color='k', label='Default wind spin up')

	plt.plot(t, HCs[0], label=labels[0])

	plt.plot(t, HCs[1], label=labels[1])

	plt.legend()

	plt.xlabel('Time (years)')

	plt.title('On shelf heat content (J), ' + bathyName)

	

	plt.grid()

	plt.show()

	

	quit()

	

#==



barostr_timeSeries = 0

if barostr_timeSeries:



	#paths = ['MCS_154', ['MCS_315', 'MCS_314']]; labels = ['IPO unif neg', 'IPO unif pos']; bathyName = 'bathyS'; xx = [50,190]; tt=[0,160]

	#paths = ['MCS_154', ['MCS_317', 'MCS_316']]; labels = ['IPO sin neg', 'IPO sin pos']; bathyName = 'bathyS'; xx = [50,190]; tt=[0,205]

	paths = ['MCS_154', ['MCS_312', 'MCS_313']]; labels = ['IPO PH neg', 'IPO PH pos']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, None]

	

	ttRef = [0, None]

	d = 4

	

	path = path_root + paths[0] + '/run/'

	grid = Grid(path)

	bathy = grid.bathy

	

	X = grid.XC[1,:]/1000.

	Y = grid.YC[:,1]/1000.

	xlabel = 'LON (km)'; ylabel = 'LAT (km)'

	ny, nx = bathy.shape

	

	tmp = readVariable('UVEL', path, file_format='nc', tt=ttRef)[::d]

	tmp = - 1.e-6 * tools.barotropicStreamfunction(tmp, grid, timeDep=True, norm=False)

	

	print(tmp.shape)

	baroRef = np.abs(np.min(tmp[:,0:LAT,:], axis=(1,2)))

	Ntref = len(baroRef)

	tref = np.linspace(1, Ntref, Ntref) / 12

	print(Ntref)

	

	#==

	

	baros = []

	for path_tmp in paths[1]:

		

		print(path_tmp)

		tmp = readVariable('UVEL', path_root+path_tmp+'/run/', file_format='nc', meta=False, tt=tt)[::d]

		tmp = - 1.e-6 * tools.barotropicStreamfunction(tmp, grid, timeDep=True, norm=False)

		baros.append(np.abs(np.min(tmp[:,0:LAT,:], axis=(1,2))))

		

	#==

	

	Nt = len(baros[1])

	t = np.linspace(Ntref+1, Ntref+1+Nt, Nt) / 12

	

	plt.plot(tref, baroRef, color='k', label='Default wind spin up')

	plt.plot(t, baros[0], label=labels[0])

	plt.plot(t, baros[1], label=labels[1])

	plt.legend()

	plt.xlabel('Time (years)')

	plt.title('Barotropic streamfunction magnitude')

	

	plt.grid()

	plt.show()



	quit()

	

#==



barotropicStreamfunction = 0

if barotropicStreamfunction:



	#path = '/home/michael/Documents/data/PAS_666/run/'

	#path = '/home/michael/Documents/data/PISOMIP_003/run/'

	path_root = '/home/michael/Documents/data/'



	run = 'MCS_154'

	

	tt = [0, None]

	d = 1

	

	path = path_root + run + '/run/'

	grid = Grid(path)

	bathy = grid.bathy

	

	X = grid.XC[1,:]/1000.

	Y = grid.YC[:,1]/1000.

	xlabel = 'LON (km)'; ylabel = 'LAT (km)'

	ny, nx = bathy.shape

	

	#plt.plot(bathy); plt.grid(); plt.show(); quit()

		

	#VAR = 'ISOTHERM'

	VAR = 'BAROSTR'

	#VAR = 'BOTVEL'



	if VAR == 'ISOTHERM':

	

		data = readVariable('ETAN', path, file_format='nc', tt=tt)

		

		data = np.load(path+'ThermZ_m05_'+run+'.npy')

		vmin = -600; vmax=-150;	cmap = 'YlOrRd'

		#vmin = -500; vmax = -350; cmap = 'jet'

		title = run + ' -0.5 deg. C isotherm depth'

		print(data.shape)



	#==

		

	elif VAR == 'BAROSTR':

			

		d = 4

		data = readVariable('UVEL', path, file_format='nc', tt=tt)[::d]

			

		vmin = -0.6e-5; vmax = 0.6e-5; cmap = 'jet'

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



	nt = data.shape[0]*d

	time_s = np.linspace(d,nt+d,nt)*86400.*30.

	time_s = time_s[::d]

	text_data = ptt.getTextData(time_s, 'month', X[1], Y[1], color='k')

	

	data = ptt.maskBathyXY(data, grid, 0, timeDep=True)

	

	pt.animate1by1(data, X, Y, cmap=cmap, xlabel=xlabel, ylabel=ylabel, title=title, mesh=True, text_data=text_data, outname='animate1by1Surf.mp4', vmin=vmin, vmax=vmax, contour=grid.bathy)



	quit()



#==



heatTransport = 1

if heatTransport:



	#paths = ['MCS_154', ['MCS_315', 'MCS_314']]; labels = ['IPO unif neg', 'IPO unif pos']; bathyName = 'bathyS'; xx = [50,190]; tt=[0,160]

	#paths = ['MCS_154', ['MCS_317', 'MCS_316']]; labels = ['IPO sin neg', 'IPO sin pos']; bathyName = 'bathyS'; xx = [50,190]; tt=[0,205]

	paths = ['MCS_154', ['MCS_312', 'MCS_313']]; labels = ['IPO PH neg', 'IPO PH pos']; bathyName = 'bathyS'; xx = [50,190]; tt = [0, None]

	

	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)

	rho0 = 1030. # Units kg / m3

	

	lat = 50

	

	vRef = readVariable('VVEL', path_root+paths[0]+'/run/', file_format='nc', meta=False, yy=lat)

	thetaRef = readVariable('THETA', path_root+paths[0]+'/run/', file_format='nc', meta=False, yy=lat)

	vtRef = np.where(thetaRef>0, vRef*thetaRef, 0)

	vtRef = np.sum(vtRef[...,120:], axis=(1,2))

	Ntref = len(vtRef)

	tref = np.linspace(1, Ntref, Ntref) / 12

	

	vts = []

	for path_tmp in paths[1]:

		print(path_tmp)

		

		v = readVariable('VVEL', path_root+path_tmp+'/run/', file_format='nc', meta=False, yy=lat, tt=tt)

		theta = readVariable('THETA', path_root+path_tmp+'/run/', file_format='nc', meta=False, yy=lat, tt=tt)

		vt = np.where(theta>0, v*theta, 0)

		vts.append(np.sum(vt[...,120:], axis=(1,2)))

		

	Nt = len(vts[1])

	t = np.linspace(Ntref+1, Ntref+1+Nt, Nt) / 12

	

	plt.plot(tref, vtRef, color='k', label='Default wind spin up')

	plt.plot(t, vts[0], label=labels[0])

	plt.plot(t, vts[1], label=labels[1])

	plt.legend()

	plt.xlabel('Time (years)')

	plt.title('On shelf heat content (J), ' + bathyName)

	

	plt.grid()

	plt.show()

		

	

	



#==



# Cross-shelf heat transport plots.

heatFluxes_heatRel = 0

if heatFluxes_heatRel:



	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)

	rho0 = 1030. # Units kg / m3

	

	dt_month = 30. * 86400.

	dt_rel = 1. / (100. * 86400.)

	

	# First experiment defined here.

	exp = 'MCS_310'

	

	path = path_root + exp + '/run/'#'MCS_158/run/'

	grid = Grid(path)



	X = grid.XC[1,:]/1000.

	Y = grid.YC[:,1]/1000.

	Z = grid.RC.squeeze()



	lat = LAT



	#v1 = np.mean(readVariable('VVEL', path, meta=False, tt=-1), axis=-1)

	#v2 = np.mean(readVariable('VVEL', path_root + exp3 + '/run/', meta=False, tt=-1), axis=-1)

	#v1 = ptt.maskBathyYZ(v1, grid, xi=120, timeDep=False); v2 = ptt.maskBathyYZ(v2, grid, xi=120, timeDep=False)

	#pt.plot1by2([v1,v2],X=[Y,Y],Y=[Z,Z], vmin=-0.001,vmax=0.001,mesh=True); quit()

	

	#==

	

	# 1. 

	# First get heat transport across shelf break



	#T = readVariable('THETA', path, meta=False, yy=lat)

	#Tf = -1.8 * np.ones(T.shape)

	#v = np.mean(readVariable('VVEL', path, meta=False, yy=[lat,lat+2]), axis=-2)

	#T = rho0 * Cp * v * (T - Tf)

	

	T = rho0 * Cp * readVariable('VVELTH', path, meta=False, yy=lat)

	

	area = grid.DXG[lat] * grid.hFacS[:,lat] * grid.DRF[:,0]

	area = ptt.maskBathyXZ(area, grid, yi=lat, timeDep=False)



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



	Tshelf = np.ma.cumsum(Tshelf); Trel = np.ma.cumsum(Trel)



	plt.plot(Tshelf, label='Tshelf')

	#plt.plot(Trel, label='Trel')

	#plt.plot(Tshelf-Trel, label='sum')

	plt.legend()

	plt.grid()

	plt.show()



#==



# Cross-shelf heat transport plots.

heatFluxes = 0

if heatFluxes:



	path_root = '/home/michael/Documents/data/'

	SMOOTH = True

	TF = False

	WARMFLUXES = False

	

	Cp = 3974.0 # Units J / (kg C) = (kg m2 / s2) / (kg C) = m2 / (s2 C)

	rho0 = 1030. # Units kg / m3

	normval = 1.e12

	

	# First experiment defined here.

	#exp1 = 'MCS_154'; exp2 = 'MCS_301'; exp3 = 'MCS_302'

	#exp1 = 'MCS_303'; exp2 = 'MCS_304'; exp3 = 'MCS_305'

	exp1 = 'MCS_308'; exp2 = 'MCS_309';	exp3 = 'MCS_310'

	

	path = path_root + exp1 + '/run/'#'MCS_158/run/'

	grid = Grid(path)



	# Subregions for heat transport.



	lat = LAT



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

		ylims = [-.4, .2]

	else:

		#ylims = [-10,10]	

		ylims = [-12,12]

		

	xticks = np.linspace(0, 360, 10)

	xticks = [xticks]*3

	xticksvis = [False, False, True]

	

	yticks = np.linspace(ylims[0], ylims[1], 7)

	yticks = [yticks]*3

	yticksvis = [True, True, True]

	

	fdict = {'fontsize':12, 'color':'k'}

	text = ['{:.3f}'.format(av1), '{:.3f}'.format(av2), '{:.3f}'.format(av3)]

	text = ['Heat transport "all" average = ' + t + ' TW' for t in text]

	xloc = 100; yloc = 0.8*ylims[1]

	text_data = {'text':text, 'xloc':xloc, 'yloc':yloc, 'fontdict':fdict}

	

	#==



    # NOW PLOT

	pt.line1byN(Ts, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis, text_data=text_data, loc=1, save=True, outname='heatFluxes.png')



#==



	

	

	
