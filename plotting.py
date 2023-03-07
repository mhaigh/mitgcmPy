# plotting.py

import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

from plotting_tools import setText, getContourfLevels, maskBathyXY, maskDraftXY, makeList, doTicks, doTitle, doLabels

#==========================================================

def contourf(data):
	'''Quick contourf plot of data.'''
	
	plt.contourf(data)
	plt.colorbar()
	plt.show()	

#==

def plotBathymetry(grid, xlims=None, ylims=None):
	
	x = grid.XC; y = grid.YC
	if xlims is not None:
		x = x[xlims]
	if ylims is not None:
		y = y[ylims]
	
	plt.contourf(x, y, grid.bathy)
	plt.show()
	
#==

def timeSeries(data, TIME=None, Y=None, time_units='days', figsize=(6,3), labels=None, title=None, fontsize=14, ylabel=None, gridOn=True, save=False, outpath='', outname='timeSeries.png', show=True, dpi=200):
	'''1D time series of entries in listed data.'''
	
	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
	
	if time_units == 'days':
		T = TIME / 86400.
	elif time_units == 'months':
		T = TIME / (30*86400.)
	xlabel = 'Time (%s)' % time_units
				
	if TIME is not None:
		if labels is not None:
			for di in range(len(data)):
				plt.plot(T, data[di], label=labels[di])
		else:
			for d in data:
				plt.plot(T, d)			
	else:
		if labels is not None:
			for di in range(len(data)):
				plt.plot(data[di], label=labels[di])
		else:
			for d in data:
				plt.plot(d)	
			
	plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)
	
	if title is not None:
		plt.title(titles, fontsize=fontsize)
		
	if labels is not None:
		plt.legend()
		
	if gridOn:
		plt.grid()
		
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
	
#==

def plot1by1(data, X=None, Y=None, contour=None, contourlevels=None, figsize=(5,4), title=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabel=None, ylabel=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, hlines=None, xmin=[0], xmax=[1], vlines=None, xlim=None, ylim=None):
	
	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	plt.subplot(111)
	plt.gca().patch.set_color('.25')

	if mesh:
		if X is not None and Y is not None:
			cax = plt.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			cax = plt.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax)
	else:
		levels = getContourfLevels(vmin, vmax, contourfNlevels)
		if X is not None and Y is not None:
			cax = plt.contourf(X, Y, data, cmap=cmap, levels=levels)
		else: 
			cax = plt.contourf(data, cmap=cmap, levels=levels)

	if contour is not None:
		plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourlevels)
			
	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)
	
	if grid:
		plt.grid()

	if text_data is not None:
		setText(plt.gca(), text_data, set_invisible=False)

	plt.colorbar(cax, ax=ax)
	
	if title is not None:
		plt.title(title, fontsize=fontsize)

	#if yline is not None:
	#	plt.plot(yline, Y, color='k', linewidth=1.2)
	#	plt.axvline(x=X[X.shape[0]//2-1], color='k', linewidth=1.0, linestyle='--')

	if vlines is not None:
		for i in range(len(vlines)):
			plt.axvline(x=vlines[i], color='k', linewidth=1.0, linestyle='--')
	
	if hlines is not None:
		for i in range(len(hlines)):
			plt.axhline(y=hlines[i], xmin=xmin[i], xmax=xmax[i], color='k', linewidth=1.0, linestyle='--')
			
	if xlim is not None:
		plt.xlim(xlim)
	if ylim is not None:
		plt.ylim(ylim)
		
	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
	
#==

def line1byN(data, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis, fontsize=14, figsize=(9,3), save=False, outname='line1byN.png', show=True):

	N = len(data)
	if N == 3:
		figsize = (9,8)
	elif N == 4:
		figsize = (9,10.5)
		
	fig = plt.figure(figsize=figsize)
	
	#fig.suptitle('Meridional heat transport across shelf break')
	
	for pi in range(N):
		plt.subplot(N,1,pi+1)
		
		for i in range(len(data[pi])):
			plt.plot(data[pi][i], label=labels[pi][i], color=colours[pi][i])
			
		plt.xlim(xlims)
		plt.ylim(ylims)
		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])
		
		plt.title(titles[pi])
		plt.xlabel(xlabels[pi], fontsize=fontsize)
		plt.ylabel(ylabels[pi], fontsize=fontsize-2)
		plt.grid()
	
	plt.legend(loc=2) 
	
	plt.tight_layout()
		
	if save:
		plt.savefig(outname)
		
	if show:
		plt.show()

#==

def quiver1by1(u, v, Xd, Yd, scale=1, C=None, ccmap='bwr', contour=None, X=None, Y=None, contourf=True, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1.mp4', show=True, dpi=200, text_data=None):

	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	plt.subplot(111)
	plt.gca().patch.set_color('.25')

	if contour is not None:
		if contourf:
			plt.contourf(X, Y, contour, cmap=cmap)
		else:
			plt.pcolormesh(X, Y, contour, vmin=vmin, vmax=vmax, cmap=cmap)		
		plt.colorbar()

	if C is not None:
		cax = ax.quiver(Xd, Yd, u, v, C, cmap=ccmap, scale=scale)
		plt.colorbar(cax, ax=ax)
	else:
		cax = ax.quiver(Xd, Yd, u, v, scale=scale)

	ax.quiverkey(cax, 0.1, 0.1, .2, '0.2 m/s', labelpos='E', coordinates='axes')
			
	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)

	if text_data is not None:
		setText(plt.gca(), text_data, set_invisible=False)
	
	if title is not None:
		plt.title(title, fontsize=fontsize)

	#==
	
	plt.tight_layout()
		
	if show:
		plt.show()

	if save:
		plt.savefig(outpath + outname)
		plt.close()

#==

def quiver1byN(u, v, Xd, Yd, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, figsize=None, title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, text_data=None, width_ratios=None, labelData=None, cbar=True, grid=True, xticks=None, xticksvis=True, yticks=None, yticksvis=True, scale=2):

	N = len(u)

	Xd = makeList(Xd, N)
	Yd = makeList(Yd, N)

	if contourf is not None:
		X = makeList(X, N)
		Y = makeList(Y, N)

	vmin = makeList(vmin, N)
	vmax = makeList(vmax, N)

	cbar = makeList(cbar, N)
	grid = makeList(grid, N)
	title = makeList(title, N)

	xticks = makeList(xticks, N)
	xticksvis = makeList(xticksvis, N)

	yticks = makeList(yticks, N)
	yticksvis = makeList(yticksvis, N)

	xlabel = makeList(xlabel, N)
	ylabel = makeList(ylabel, N)

	if figsize is None:
		if N == 2:
			figsize = (8,3)
		elif N == 3:
			figsize = (11,3)
		else:
			figsize = (8,3)

	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)

	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=N, nrows=1, figure=fig, width_ratios=width_ratios)

	for pi in range(N):
		if width_ratios is not None:			
			ax = fig.add_subplot(gs[0,pi])
		else:
			plt.subplot(1,N,pi+1)
			ax = plt.gca()

		ax.patch.set_color('.5')

		if contourf is not None:
			if mesh:
				plt.pcolormesh(X[pi], Y[pi], contourf[pi], vmin=vmin[pi], vmax=vmax[pi], cmap=cmap)		
			else:
				levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)
				plt.contourf(X[pi], Y[pi], contourf[pi], cmap=cmap, levels=levels)

			if cbar[pi]:
				plt.colorbar()

		if contour is not None:
			if contourLevels is not None:
				plt.contour(X[pi], Y[pi], contour[pi], levels=contourLevels[pi], colors='k', linestyles='solid', linewidths=0.4)
			else:
				plt.contour(X[pi], Y[pi], contour[pi], levels=contourLevels[pi], colors='k', linestyles='solid', linewidths=0.4)

		if C is not None:
			cax = ax.quiver(Xd[pi], Yd[pi], u[pi], v[pi], C[pi], cmap=ccmap, scale=scale)
			if width_ratios is None:
				plt.colorbar(cax, ax=ax)
		else:
			cax = plt.quiver(Xd[pi], Yd[pi], u[pi], v[pi], scale=scale)
		ax.quiverkey(cax, 0.12, 0.03, .1, '0.1 m/s', labelpos='N', coordinates='axes')
				
		doLabels(xlabel[pi], ylabel[pi], fontsize=fontsize)
		doTitle(title[pi], fontsize=fontsize)

		if text_data is not None:
			setText(ax, text_data[pi], set_invisible=False)
		
		if labelData is not None:
			for li in labelData[pi]:
				plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
				plt.annotate(li['t'], li['x'])

		#ax.set_aspect('equal')

		if grid[pi]:
			plt.grid()

		if yticks[pi] is not None:
			if yticksvis[pi]:				
				plt.yticks(yticks[pi])
			else:			
				plt.yticks(yticks[pi], labels='')

	#==
	
	#fig.subplots_adjust(wspace=-1, hspace=0)

	plt.tight_layout()

	if save:
		plt.savefig(outpath + outname)

	if show:
		plt.show()

	plt.close()

#==

def quiver2by2(u, v, Xd, Yd, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, cbarShared=False, cbarSharedData=None, figsize=None, title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabels=None, ylabels=None, save=False, outpath='', outname='quiver2by2.png', show=True, dpi=200, text_data=None, width_ratios=[1,1], labelData=None, cbar=True, grid=True, xticks=None, xticksvis=True, yticks=None, yticksvis=True, scale=2):

	if len(u) != 2 or len(u[0]) != 2 or len(u[1]) != 2:
		print('Input data u and v must both be 2 by 2 lists.')
		return	
		
	M = 2
	N = 2

	Xd = makeList(Xd, M, N)
	Yd = makeList(Yd, M, N)

	if contourf is not None:
		X = makeList(X, M, N)
		Y = makeList(Y, M, N)

	vmin = makeList(vmin, M, N)
	vmax = makeList(vmax, M, N)
	C = makeList(C, M, N)
	
	cbar = makeList(cbar, M, N)
	grid = makeList(grid, M, N)
	title = makeList(title, M, N)

	xticks = makeList(xticks, M, N)
	xticksvis = makeList(xticksvis, M, N)
	yticks = makeList(yticks, M, N)
	yticksvis = makeList(yticksvis, M, N)

	xlabels = makeList(xlabels, M, N)
	ylabels = makeList(ylabels, M, N)


	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)

	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)

	for col in range(M):
		for row in range(N):

			ax = fig.add_subplot(gs[row,col])
			ax.patch.set_color('.5')

			if contourf[row][col] is not None:
				if mesh:
					im = plt.pcolormesh(X[row][col], Y[row][col], contourf[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap)		
				else:
					levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels)
					im = plt.contourf(X[row][col], Y[row][col], contourf[row][col], cmap=cmap, levels=levels)

			if cbar[row][col]:
				plt.colorbar()

			if C[row][col] is not None:
				cax = ax.quiver(Xd[row][col], Yd[row][col], u[row][col], v[row][col], C[row][col], cmap=ccmap, scale=scale)	
				plt.colorbar(cax, ax=ax)
			else:
				cax = plt.quiver(Xd[row][col], Yd[row][col], u[row][col], v[row][col], scale=scale)
			ax.quiverkey(cax, 0.12, 0.03, .1, '0.1 m/s', labelpos='N', coordinates='axes')
					
			doLabels(xlabels[row][col], ylabels[row][col], fontsize=fontsize)
			doTitle(title[row][col], fontsize=fontsize)
			doTicks(xticks[row][col], xticksvis[row][col], yticks[row][col], yticksvis[row][col])
			
			if text_data is not None:
				setText(ax, text_data[row][col], set_invisible=False)
			
			if labelData is not None:
				for li in labelData[row][col]:
					plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
					plt.annotate(li['t'], li['x'])

			#ax.set_aspect('equal')

			if grid[row][col]:
				plt.grid()
			
	if contourf:	
		if cbarShared:
			fig.subplots_adjust(right=0.8)
			cbar_ax = fig.add_axes(cbarSharedData[0])
			cbar = fig.colorbar(im, cax=cbar_ax)
			if cbarSharedData[1] is not None:
				cbar.set_label(cbarSharedData[1], rotation=270, labelpad=10)

	#==
	
	#fig.subplots_adjust(wspace=-1, hspace=0)

	#plt.tight_layout()

	if save:
		plt.savefig(outpath + outname)

	if show:
		plt.show()

	plt.close()

#==

def quiver1by2(u, v, Xd, Yd, qs=[0.1]*2, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, figsize=(8,3), title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabels=None, ylabels=None, save=False, outpath='', outname='quiver1by2.png', show=True, dpi=200, text_data=None,  width_ratios=None, labelData=None, scale=2, xticks=None, yticks=None, xticksvis=True, yticksvis=True):
	
	xticks = makeList(xticks, 2)
	yticks = makeList(yticks, 2)
	xticksvis = makeList(xticksvis, 2)
	yticksvis = makeList(yticksvis, 2)
	
	scale = makeList(scale, 2)
	
	fig = plt.figure(figsize=figsize, dpi=dpi)

	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)
		ax = fig.add_subplot(gs[0,0])
	else:
		plt.subplot(121)
		ax = plt.gca()

	ax.patch.set_color('.5')

	if contourf is not None:
		if mesh:
			plt.pcolormesh(X[0], Y[0], contourf[0], vmin=vmin[0], vmax=vmax[0], cmap=cmap)		
		else:
			levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)
			plt.contourf(X[0], Y[0], contourf[0], cmap=cmap, levels=levels)

		if width_ratios is None:
			plt.colorbar()

	if contour is not None:
		if contourLevels is not None:
			plt.contour(X[0], Y[0], contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.4)
		else:
			plt.contour(X[0], Y[0], contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.4)

	if C is not None:
		cax = ax.quiver(Xd[0], Yd[0], u[0], v[0], C[0], cmap=ccmap, scale=scale[0])
		if width_ratios is None:
			print('1')
			plt.colorbar(cax, ax=ax)
	else:
		cax = plt.quiver(Xd[0], Yd[0], u[0], v[0], scale=scale)
	ax.quiverkey(cax, 0.12, 0.03, qs[0], str(qs[0]) + ' m/s', labelpos='N', coordinates='axes')
			
		
	doLabels(xlabels[0], ylabels[0], fontsize=fontsize)
	doTicks(xticks[0], xticksvis[0], yticks[0], yticksvis[0])

	if text_data is not None:
		setText(ax, text_data[0], set_invisible=False)
	
	if title is not None:
		plt.title(title[0], fontsize=8)
		
	#plt.grid()

	#ax.set_aspect('equal')

	if labelData is not None:
		for li in labelData[0]:
			plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
			plt.annotate(li['t'], li['x'])

	#==

	plt.subplot(122)
	ax = plt.gca()
	ax.patch.set_color('.5')

	if contourf is not None:
		if mesh:
			plt.pcolormesh(X[1], Y[1], contourf[1], vmin=vmin[1], vmax=vmax[1], cmap=cmap)		
		else:
			levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)
			plt.contourf(X[1], Y[1], contourf[1], cmap=cmap, levels=levels)
		plt.colorbar()

	if C is not None:
		cax = plt.quiver(Xd[1], Yd[1], u[1], v[1], C[1], cmap=ccmap, scale=scale[1], linewidths=2)
		plt.colorbar(cax, ax=ax)
	else:
		cax = ax.quiver(Xd[1], Yd[1], u[1], v[1], scale=scale)
	ax.quiverkey(cax, 0.12, 0.03, qs[1], str(qs[1]) + ' m/s', labelpos='N', coordinates='axes')
			
	doLabels(xlabels[1], ylabels[1], fontsize=fontsize)
	doTicks(xticks[1], xticksvis[1], yticks[1], yticksvis[1])
	
	if text_data is not None:
		setText(ax, text_data[1], set_invisible=False)
	
	if title is not None:
		plt.title(title[1], fontsize=8)

	#plt.grid()
	#ax.set_aspect('equal')

	if labelData is not None:
		for li in labelData[1]:
			plt.scatter(li['x'][0], li['x'][1], s=1, color='k')
			plt.annotate(li['t'], li['x'])

	#==
	
	plt.tight_layout()
		
	if show:
		plt.show()

	if save:
		plt.savefig(outpath + outname)
		plt.close()
		
#==

def quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, figsize=(8,3), title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, save=False, outpath='', outname='quiver1by2.mp4', show=True, dpi=200, text_data=None, parallels=None, meridians=None, width_ratios=None, labelData=None, isf=[None, None]):
	
	from mpl_toolkits.basemap import Basemap
		
	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	fig = plt.figure(figsize=figsize, dpi=dpi)

	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)
		ax = fig.add_subplot(gs[0,0])
	else:
		plt.subplot(121)
		ax = plt.gca()

	ax.patch.set_color('.5')

	X0, Y0 = m(X[0],Y[0])
	Xd0, Yd0 = m(Xd[0],Yd[0])

	if contourf is not None:
		if mesh:
			m.pcolormesh(X0, Y0, contourf[0], vmin=vmin[0], vmax=vmax[0], cmap=cmap)		
		else:
			levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)
			m.contourf(X0, Y0, contourf[0], cmap=cmap, levels=levels)

		if width_ratios is None:
			plt.colorbar()

	if contour is not None:
		if contourLevels is not None:
			m.contour(X0, Y0, contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.4)
		else:
			m.contour(X0, Y0, contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.4)

	if C is not None:
		cax = ax.quiver(Xd0, Yd0, u[0], v[0], C[0], cmap=ccmap, scale=2)
		if width_ratios is None:
			plt.colorbar(cax, ax=ax)
	else:
		cax = m.quiver(Xd0, Yd0, u[0], v[0], scale=2.5)
	ax.quiverkey(cax, 0.12, 0.03, .2, '0.2 m/s', labelpos='E', coordinates='axes')
			
	if xlabel is not None:
		plt.xlabel(xlabel[0], fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel[0], fontsize=fontsize)
	
	if title is not None:
		plt.title(title[0], fontsize=fontsize)

	ax.set_aspect('equal')
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	if labelData is not None:
		for li in labelData[0]:
			plt.scatter(li['x'][0], li['x'][1], s=1.5, color='r')
			plt.annotate(li['t'], li['tx'], fontsize=12)

	if isf[0] is not None:
		extent = [X[0][0,0], X[0][0,-1], -Y[0][0,0], -Y[0][-1,0]]
		m.imshow(1-isf[0], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=1)
		
	if text_data is not None:
		setText(ax, text_data[0], set_invisible=False)
		
	#==

	plt.subplot(122)
	ax = plt.gca()
	ax.patch.set_color('.5')

	X0, Y0 = m(X[1],Y[1])
	Xd0, Yd0 = m(Xd[1],Yd[1])

	if contourf is not None:
		if mesh:
			m.pcolormesh(X0, Y0, contourf[1], vmin=vmin[1], vmax=vmax[1], cmap=cmap)		
		else:
			levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)
			m.contourf(X0, Y0, contourf[1], cmap=cmap, levels=levels)
		plt.colorbar()

	if C is not None:
		cax = m.quiver(Xd0, Yd0, u[1], v[1], C[1], cmap=ccmap, scale=1., linewidths=2)
		plt.colorbar(cax, ax=ax)
	else:
		cax = ax.quiver(Xd0, Yd0, u[1], v[1], scale=2)
	ax.quiverkey(cax, 0.22, 0.03, .2, '0.2 m/s', labelpos='E', coordinates='axes')
			
	if xlabel is not None:
		plt.xlabel(xlabel[1], fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel[1], fontsize=fontsize)
	
	if title is not None:
		plt.title(title[1], fontsize=fontsize)

	ax.set_aspect('equal')
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[False,False,False,False], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	if labelData is not None:
		for li in labelData[1]:
			plt.scatter(li['x'][0], li['x'][1], s=1.5, color='k')
			plt.annotate(li['t'], li['x'], fontsize=12)

	if isf[1] is not None:
		extent = [X[1][0,0], X[1][0,-1], -Y[1][0,0], -Y[1][-1,0]]
		m.imshow(1-isf[1], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=1)

	if text_data is not None:
		setText(ax, text_data[1], set_invisible=False)
				
	#==
	
	plt.tight_layout()
		
	if show:
		plt.show()

	if save:
		plt.savefig(outpath + outname)
		plt.close()

#==

def PAS_2by2_uvs(u, v, Xd, Yd, C=None, ccmap='bwr', contour=None, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1.mp4', show=True, dpi=200, text_data=None):

	fs = 8


	fig = plt.figure(figsize=figsize, dpi=dpi)
	fig.suptitle('PAS Seasonal flow and salinity')

	gs = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, width_ratios=[0.95,1.2])

	# Top left
	ax = fig.add_subplot(gs[0,0])
	plt.gca().patch.set_color('.25')

	plt.pcolormesh(X, Y, contour, vmin=vmin, vmax=vmax, cmap=cmap)		

	cax = ax.quiver(Xd, Yd, u[0][0], v[0][0], C[0][0], cmap=ccmap)

	ax.quiverkey(cax, 0.1, 0.1, .2, '0.2 m/s', labelpos='E', coordinates='axes')
			
	plt.xlabel(xlabel, fontsize=fontsize)

	ax.axes.xaxis.set_ticklabels([])
	plt.text(X[1], Y[1], 'SON', color='w', fontsize=fs)

	#==

	# Top right
	ax = fig.add_subplot(gs[0,1])
	plt.gca().patch.set_color('.25')

	plt.pcolormesh(X, Y, contour, vmin=vmin, vmax=vmax, cmap=cmap)		
	plt.colorbar()

	cax = plt.gca().quiver(Xd, Yd, u[0][1], v[0][1], C[0][1], cmap=ccmap)

	plt.gca().quiverkey(cax, 0.1, 0.1, .2, '0.2 m/s', labelpos='E', coordinates='axes')
	
	ax.axes.yaxis.set_ticklabels([])
	ax.axes.xaxis.set_ticklabels([])
	plt.text(X[1], Y[1], 'DJF', color='w', fontsize=fs)

	#==

	# Bottom left
	ax = fig.add_subplot(gs[1,0])
	plt.gca().patch.set_color('.25')

	plt.pcolormesh(X, Y, contour, vmin=vmin, vmax=vmax, cmap=cmap)

	cax = ax.quiver(Xd, Yd, u[1][0], v[1][0], C[1][0], cmap=ccmap)

	ax.quiverkey(cax, 0.1, 0.1, .2, '0.2 m/s', labelpos='E', coordinates='axes')
			
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.text(X[1], Y[1], 'MAM', color='w', fontsize=fs)

	#==

	# Bottom right
	ax = fig.add_subplot(gs[1,1])
	plt.gca().patch.set_color('.25')

	plt.pcolormesh(X, Y, contour, vmin=vmin, vmax=vmax, cmap=cmap)		
	
	cax = plt.gca().quiver(Xd, Yd, u[1][1], v[1][1], C[1][1], cmap=ccmap)
	plt.colorbar(cax, ax=ax)

	plt.gca().quiverkey(cax, 0.1, 0.1, .2, '0.2 m/s', labelpos='E', coordinates='axes')

			
	plt.xlabel(xlabel, fontsize=fontsize)

	plt.text(X[1], Y[1], 'JJA', color='w', fontsize=fs)

	ax.axes.yaxis.set_ticklabels([])

	#==
	
	plt.tight_layout()
		
	if show:
		plt.show()

	if save:
		plt.savefig(outpath + outname)
		plt.close()

#==
		
def plot1by2(data, X=[None]*2, Y=[None]*2, contour=[None]*2, contourlevels=[None]*2, figsize=(9,4), titles=[None]*2, fontsize=14, mesh=False, cmaps=['jet']*2, vmin=[None]*2, vmax=[None]*2, text_data=[None]*2, xlabels=[None]*2, ylabels=[None]*2, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*2, cbarLabel=[None]*2, xticks=[None]*2, yticks=[None]*2, xticksvis=[True]*2, yticksvis=[True]*2):
	
	vmin = makeList(vmin, 2)
	vmax = makeList(vmax, 2)
	cmaps = makeList(cmaps, 2)
	
	fig = plt.figure(figsize=figsize, dpi=dpi)
	
	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)

	for pi in range(2):

		if width_ratios is not None:
			ax = fig.add_subplot(gs[0,pi])
		else:
			plt.subplot(1,2,pi+1)
			ax = plt.gca()

		if pi == 1:
			ax1 = ax

		ax.patch.set_color('.5')

		if mesh:
			if X[pi] is not None and Y[pi] is not None:
				plt.pcolormesh(X[pi], Y[pi], data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])
			else: 
				plt.pcolormesh(data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])
		else:
			levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)
			if X[pi] is not None and Y[pi] is not None:
				plt.contourf(X[pi], Y[pi], data[pi], cmap=cmaps[pi], levels=levels)
			else: 
				plt.contourf(data[pi], cmap=cmaps[pi], levels=levels)

		if cbar[pi]:
			plt.colorbar()
			if cbarLabel[pi] is not None:
				cbar.set_label(cbarLabel[pi], rotation=270, labelpad=20, fontsize=10)

		if contour[pi] is not None:
			if X[pi] is not None and Y[pi] is not None:
				plt.contour(X[pi], Y[pi], contour[pi], levels=contourlevels[pi], colors='k', linestyles='solid', linewidths=0.8)
			else:
				plt.contour(X[pi], Y[pi], contour[pi], levels=contourlevels[pi], colors='k', linestyles='solid', linewidths=0.8)


		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)
		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])
		doTitle(titles[pi], fontsize=fontsize)
	
		if grid:
			plt.grid()

		if text_data[pi] is not None:
			setText(plt.gca(), text_data[pi], set_invisible=False)
		
	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==

def plot1by3(data, X=[None]*3, Y=[None]*3, contour=[None]*3, contourlevels=[None]*3, figsize=(11,3), titles=[None]*3, fontsize=14, mesh=False, cmaps=['jet']*3, vmin=[None]*3, vmax=[None]*3, text_data=[None]*3, xlabels=[None]*3, ylabels=[None]*3, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3):
	
	vmin = makeList(vmin, 3)
	vmax = makeList(vmax, 3)	
	
	fig = plt.figure(figsize=figsize, dpi=dpi)

	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig, width_ratios=width_ratios)

	for pi in range(3):

		if width_ratios is not None:
			ax = fig.add_subplot(gs[0,pi])
		else:
			plt.subplot(1,3,pi+1)
			ax = plt.gca()

		if pi == 1:
			ax1 = ax

		ax.patch.set_color('.5')

		if mesh:
			if X[pi] is not None and Y[pi] is not None:
				plt.pcolormesh(X[pi], Y[pi], data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])
			else: 
				plt.pcolormesh(data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])
		else:
			levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)
			if X[pi] is not None and Y[pi] is not None:
				plt.contourf(X[pi], Y[pi], data[pi], cmap=cmaps[pi], levels=levels)
			else: 
				plt.contourf(data[pi], cmap=cmaps[pi], levels=levels)

		if cbar[pi]:
			plt.colorbar()

		if contour[pi] is not None:
			if X[pi] is not None and Y[pi] is not None:
				plt.contour(X[pi], Y[pi], contour[pi], levels=contourlevels[pi], colors='k', linestyles='solid', linewidths=0.8)
			else:
				plt.contour(X[pi], Y[pi], contour[pi], levels=contourlevels[pi], colors='k', linestyles='solid', linewidths=0.8)


		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)
		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])
		doTitle(titles[pi], fontsize=fontsize)
	
		if grid:
			plt.grid()

		if text_data[pi] is not None:
			setText(plt.gca(), text_data[pi], set_invisible=False)


		
	box = ax1.get_position()
	box.x0 = box.x0 + 0.05
	ax1.set_position(box)


	#==

	#plt.tight_layout()
	
	if save:
		plt.savefig(outpath + outname, bbox_inches="tight")
	if show:
		plt.show()
		
#==

def plotMbyN(data, X=None, Y=None, figsize=(8,4), titles=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, cbar=True, cbarShared=False, cbarSharedData=None, text_data=None, xlabels=None, ylabels=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plotMbyN.png', show=True, dpi=200, width_ratios=None, xticks=None, yticks=None, xticksvis=True, yticksvis=True, hlines=None):

	if not isinstance(data, list):
		plot1by1(data, X=X, Y=Y, mesh=mesh, vmin=vmin, vmax=vmax)
		return
		
	N = len(data)
	M = len(data[0])
	
	X = makeList(X, M, N)
	Y = makeList(Y, M, N)

	hlines = makeList(hlines, M, N)
	
	vmin = makeList(vmin, M, N)
	vmax = makeList(vmax, M, N)	
	
	cbar = makeList(cbar, M, N)
	cmap = makeList(cmap, M, N)
	grid = makeList(grid, M, N)
	text_data = makeList(text_data, M, N)
	
	titles = makeList(titles, M, N)
	xlabels = makeList(xlabels, M, N)
	ylabels = makeList(ylabels, M, N)
	
	xticks = makeList(xticks, M, N)
	yticks = makeList(yticks, M, N)
	xticksvis = makeList(xticksvis, M, N)
	yticksvis = makeList(yticksvis, M, N)
			
	if width_ratios is None:
		width_ratios = [1]*M
	
	#==
	
	fig = plt.figure(figsize=figsize, dpi=dpi)
	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)
	
	for col in range(M):
		for row in range(N):

			ax = fig.add_subplot(gs[row,col])
			ax.patch.set_color('.5')

			if mesh:
				im = plt.pcolormesh(X[row][col], Y[row][col], data[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap[row][col])		
			else:
				levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels)
				im = plt.contourf(X[row][col], Y[row][col], data[row][col], cmap=cmap[row][col], levels=levels, extend='both')
			if cbar[row][col]:
				plt.colorbar()
				
			doTitle(titles[row][col], fontsize=fontsize)
			doTicks(xticks[row][col], xticksvis[row][col], yticks[row][col], yticksvis[row][col])
			doLabels(xlabels[row][col], ylabels[row][col], fontsize=fontsize)

			if hlines[row][col] is not None:
				plt.axhline(y=hlines[row][col], color='k', linewidth=0.7, linestyle='--')
		
			if text_data[row][col] is not None:
				setText(ax, text_data[row][col])
				
			if grid[row][col]:
				plt.grid()
				
	#==

	if cbarShared:
		fig.subplots_adjust(right=0.8)
		cbar_ax = fig.add_axes(cbarSharedData[0])
		cbar = fig.colorbar(im, cax=cbar_ax)
		if cbarSharedData[1] is not None:
			cbar.set_label(cbarSharedData[1], rotation=270, labelpad=20, fontsize=fontsize)
			
	#==
	
	#plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
			

#==

def quiver2plot1(u, v, Xd, Yd, data, X, Y, contourf=[None]*2, mesh=[True]*3, grid=[True]*3, figsize=(11,3), vmin=[None]*3, vmax=[None]*3, cmaps=['jet']*3, titles=[None]*3, fontsize=14, xlabels=[None]*3, ylabels=[None]*3, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3, contourlevels=[None]*3, contourfNlevels=9, qs=[1,1], scales=[1,1], text_data=[None]*3, qsfs=8):

#data, X=[None]*3, Y=[None]*3, contour=[None]*3, contourlevels=[None]*3, figsize=(11,3), titles=[None]*3, fontsize=14, mesh=False, cmaps=['jet']*3, vmin=[None]*3, vmax=[None]*3, text_data=[None]*3, xlabels=[None]*3, ylabels=[None]*3, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3):

	fig = plt.figure(figsize=figsize, dpi=dpi)

	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig, width_ratios=width_ratios)

	for pi in range(2):

		if width_ratios is not None:
			ax = fig.add_subplot(gs[0,pi])
		else:
			plt.subplot(1,3,pi+1)
			ax = plt.gca()

		if pi == 1:
			ax1 = ax

		ax.patch.set_color('.5')

		if contourf[pi] is not None:
			if mesh[pi]:
				im = plt.pcolormesh(X[pi], Y[pi], contourf[pi], vmin=vmin[pi], vmax=vmax[pi], cmap=cmaps[pi])		
			else:
				levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)
				im = plt.contourf(X[pi], Y[pi], contourf[pi], cmap=cmaps[pi], levels=levels)

		if cbar[pi]:
			plt.colorbar()

		cax = plt.quiver(Xd[pi], Yd[pi], u[pi], v[pi], scale=scales[pi])
		ax.quiverkey(cax, 0.1, 0.03, qs[pi], str(qs[pi]) + ' m/s', labelpos='N', coordinates='axes', fontproperties={'size':qsfs})
				
		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)
		doTitle(titles[pi], fontsize=fontsize)
		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])

		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)
		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])
		doTitle(titles[pi], fontsize=fontsize)
	
		if grid[pi]:
			plt.grid()

		if text_data[pi] is not None:
			setText(plt.gca(), text_data[pi], set_invisible=False)


	#==
	
	# 3rd panel
	
	pi = 2
	
	if width_ratios is not None:
		ax = fig.add_subplot(gs[0,pi])
	else:
		plt.subplot(1,3,pi+1)
		ax = plt.gca()
		
	ax.patch.set_color('.5')
				
	if mesh[pi]:
		plt.pcolormesh(X[pi], Y[pi], data, vmin=vmin[pi], vmax=vmax[pi], cmap=cmaps[pi], extend='both')		
	else:
		levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)
		plt.contourf(X[pi], Y[pi], data, cmap=cmaps[pi], levels=levels, extend='both')

	doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)
	doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])
	doTitle(titles[pi], fontsize=fontsize)
			
	if grid[pi]:
		plt.grid()
		
	if cbar[pi]:
		plt.colorbar()
		
	if text_data[pi] is not None:
		setText(plt.gca(), text_data[pi], set_invisible=False)
		
	#==
	
	#box = ax1.get_position()
	#box.x0 = box.x0 + 0.05
	#ax1.set_position(box)

	#==

	plt.tight_layout()
	
	if save:
		plt.savefig(outpath + outname, bbox_inches="tight")
	if show:
		plt.show()

#==

def animate1by1(data, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, contour=None, contourLevels=None):

	# Make animation
	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()
	plt.gca().patch.set_color('.25')

	if mesh:
		if X is not None and Y is not None:
			cax = ax.pcolormesh(X, Y, data[0], cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			cax = ax.pcolormesh(data[0], cmap=cmap, vmin=vmin, vmax=vmax)
		if text_data is not None:
			setText(ax, text_data, i=0)	

		if contour is not None:
			if contourLevels is None:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)
			else:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)

		plt.grid()

		def animate(i):
			cax.set_array(data[i].flatten())
			if text_data is not None:
				setText(ax, text_data, i=i, set_invisible=True)
			
		# End if mesh	

	else:
		nlevels = 9
		if X is not None and Y is not None:
			cax = ax.contourf(X, Y, data[0], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))
		else: 
			cax = ax.contourf(data[0], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))
		if text_data is not None:
			setText(ax, text_data, i=0)

		if contour is not None:
			if contourLevels is None:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)
			else:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)


		def animate(i):
			ax.clear()
			#plt.grid()
			#cax.set_data(data[i].flatten())			
			if text_data is not None:
				setText(ax, text_data, i=i, set_invisible=True)
			ax.contourf(X, Y, data[i], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))	
	#==
	
	fig.colorbar(cax, ax=ax)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.grid() 

	# Get number of timesteps/frames.
	if isinstance(data, list):
		Nt = len(data)
	else:
		Nt = data.shape[0]
		
	anim = animation.FuncAnimation(fig, animate, interval=50, frames=Nt)
	plt.tight_layout()
	plt.draw()
	
	if save:
		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)
		
	if show:
		plt.show()
		
#==

def animate1by1varCbar(data, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, contour=None, contourLevels=None):

	# Make animation
	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()
	plt.gca().patch.set_color('.25')

	if mesh:
		if X is not None and Y is not None:
			cax = ax.pcolormesh(X, Y, data[0], cmap=cmap)
		else: 
			cax = ax.pcolormesh(data[0], cmap=cmap)
			
		if text_data is not None:
			setText(ax, text_data, i=0)	

		if contour is not None:
			if contourLevels is None:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)
			else:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)

		plt.grid()

		def animate(i):
			print(i)
			ax.clear()
			cax.set_array(data[i].flatten())
			if text_data is not None:
				setText(ax, text_data, i=i, set_invisible=True)
			
		# End if mesh	

	else:
		nlevels = 9
		if X is not None and Y is not None:
			cax = ax.contourf(X, Y, data[0], cmap=cmap)
		else: 
			cax = ax.contourf(data[0], cmap=cmap)
		cb = fig.colorbar(cax, ax=ax)		
			
		if text_data is not None:
			setText(ax, text_data, i=0)

		if contour is not None:
			if contourLevels is None:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)
			else:
				plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)


		def animate(i):
			print(i)
			ax.clear()
			#plt.grid()
			#cax.set_data(data[i].flatten())			
			if text_data is not None:
				setText(ax, text_data, i=i, set_invisible=True)
			ax.contourf(X, Y, data[i], cmap=cmap)
			
	#==
	
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.grid() 

	# Get number of timesteps/frames.
	if isinstance(data, list):
		Nt = len(data)
	else:
		Nt = data.shape[0]
		
	anim = animation.FuncAnimation(fig, animate, interval=50, frames=Nt)
	plt.tight_layout()
	plt.draw()
	
	if save:
		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)
		
	if show:
		plt.show()


def animateLine(data, X=None, figsize=(5,4), title='', labels=None, fontsize=14, vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animateLine.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, constLine=None, constLineLabel=None, constLineStyle=['solid', 'dashed', 'dotted']):

	# Make animation
	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()

	if constLine is not None:
		for li in range(len(constLine)):
			if constLineLabel is not None:
				label = constLineLabel[li]
			else:
				label = None
			plt.plot(X, constLine[li], color='k', label=label, linestyle=constLineStyle[li])

	lines = []
	for di in range(len(data)):
		if labels is not None:
			line, = ax.plot(X, data[di][0], label=labels[di])
		else:
			line, = ax.plot(X, data[di][0])
		lines.append(line)

	
	if text_data is not None:
		setText(ax, text_data, i=0, set_invisible=True)

	if X is not None:
		plt.xlim(X[0], X[-1])
	plt.ylim(vmin, vmax)

	if labels is not None:
		plt.legend(loc=4)

	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.grid()
	

	def animate(i):
		for di in range(len(data)):
			lines[di].set_ydata(data[di][i])  # update the data.
		if text_data is not None:
			setText(ax, text_data, i=i, set_invisible=True)
		return lines

	anim = animation.FuncAnimation(fig, animate, interval=50, frames=data[0].shape[0])
	
	plt.tight_layout()
	plt.draw()
	
	if save:
		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)
		
	if show:
		plt.show()

#==

def animate1by1quiver(u, v, Xd, Yd, qlim=0.1, C=None, ccmap='coolwarm', contour=None, X=None, Y=None, cmap='viridis', vmin=None, vmax=None, contourf=True, grid=True, figsize=(5,4), title='', fontsize=14, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None):

	# Make animation
	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()
	plt.gca().patch.set_color('.25')

	if contour is not None:
		if len(contour.shape) == 2:
			ctr = None; cb = None
			if contourf:
				plt.contourf(X, Y, contour, cmap=cmap)
			else:
				plt.pcolormesh(X, Y, contour, vmin=vmin, vmax=vmax, cmap=cmap)		
			plt.colorbar()
		elif len(contour.shape) == 3:
			ctr = 1
			if contourf:
				C = ax.contourf(X, Y, contour[0], cmap=cmap)
			else:
				C = ax.pcolormesh(X, Y, contour[0], vmin=vmin, vmax=vmax, cmap=cmap)		
			plt.colorbar(C, ax=ax)

	if C is not None:
		Q = ax.quiver(Xd, Yd, u[0], v[0], C[0], cmap=ccmap)
		plt.colorbar(Q, ax=ax)
	else:
		Q = ax.quiver(Xd, Yd, u[0], v[0])
	ax.quiverkey(Q, 0.1, 0.1, qlim, str(qlim) + ' m/s', labelpos='E', coordinates='axes')

	if text_data is not None:
		setText(ax, text_data, i=0)

	if grid:
		plt.grid(linewidth=0.5)

	def animate(i):

		if text_data is not None:
			setText(ax, text_data, i=i, set_invisible=True)
		#if ctr is not None:
	#		if contourf:
	#			C = ax.contourf(X, Y, contour[i])
	#		else:
	#			C = ax.pcolormesh(X, Y, contour[i], vmin=vmin, vmax=vmax)
		if C is not None:
			Q.set_UVC(u[i], v[i], C[i])	
		else:		
			Q.set_UVC(u[i], v[i])
		# End if mesh	

	if title is not None:
		plt.title(title)

	anim = animation.FuncAnimation(fig, animate, interval=50, frames=u.shape[0])
	plt.tight_layout()
	plt.draw()

	#==
	
	if save:
		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)
		
	if show:
		plt.show()	

#==

def plot1by1Basemap(data, X, Y, lat_0, lon_0, contour=None, contourlevels=None, figsize=(5,4), title=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabel=None, ylabel=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, labelData=None):
	
	from mpl_toolkits.basemap import Basemap
		
	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	X,Y = m(X,Y)

	fig = plt.figure(figsize=figsize, dpi=dpi)

	ax = fig.add_subplot()
	plt.gca().patch.set_color('.5')

	if mesh:
		if X is not None and Y is not None:
			cax = m.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			cax = m.pcolormesh(data, cmap=cmap, vmin=vmin, vmax=vmax)
	else:
		levels = getContourfLevels(vmin, vmax, contourfNlevels)
		if X is not None and Y is not None:
			cax = m.contourf(X, Y, data, cmap=cmap, levels=levels)
		else: 
			cax = m.contourf(data, cmap=cmap, levels=levels)

	if contour is not None:
		m.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourlevels)
			
	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)
	
	if grid:
		plt.grid()

	if text_data is not None:
		setText(plt.gca(), text_data, set_invisible=False)

	plt.colorbar(cax, ax=ax)
	
	if title is not None:
		plt.title(title, fontsize=fontsize)

	if yline is not None:
		plt.plot(yline, Y, color='k', linewidth=1.2)
		plt.axvline(x=X[X.shape[0]//2-1], color='k', linewidth=1.0, linestyle='--')

	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	for li in labelData:
		plt.scatter(li['x'][0], li['x'][1], s=1, color='y')
		plt.annotate(li['t'], li['x'])

	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==
		
def plot1by2Basemap(data, X, Y, lat_0, lon_0, contour=None, contourlevels=None, figsize=(8,3), titles=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=[None,None], xlabels=None, ylabels=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, labelData=None):
	
	from mpl_toolkits.basemap import Basemap
		
	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	fig = plt.figure(figsize=figsize, dpi=dpi)

	plt.subplot(121)
	ax = plt.gca()
	ax.patch.set_color('.5')

	X0, Y0 = m(X[0],Y[0])

	if mesh:
		m.pcolormesh(X0, Y0, data[0], cmap=cmap, vmin=vmin[0], vmax=vmax[0])
	else:
		levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)
		m.contourf(X0, Y0, data[0], cmap=cmap, levels=levels)
			
	if xlabels is not None:
		plt.xlabel(xlabels[0], fontsize=fontsize)
	if ylabels is not None:
		plt.ylabel(ylabels[0], fontsize=fontsize)
	
	if grid:
		plt.grid()

	if text_data[0] is not None:
		setText(plt.gca(), text_data[0], set_invisible=False)

	plt.colorbar()
	
	if titles is not None:
		plt.title(titles[0], fontsize=fontsize)
	
	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])


	#==
	# Second Panel
	#==
	
	plt.subplot(122)
	ax = plt.gca()
	ax.patch.set_color('.5')


	X0, Y0 = m(X[1],Y[1])

	if mesh:
		m.pcolormesh(X0, Y0, data[1], cmap='jet', vmin=vmin[1], vmax=vmax[1])
	else:
		levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)
		m.contourf(X0, Y0, data[1], cmap=cmap, levels=levels)
	
	if xlabels is not None:
		plt.xlabel(xlabels[1], fontsize=fontsize)
	if ylabels is not None:
		plt.ylabel(ylabels[1], fontsize=fontsize)
		
	if grid:
		plt.grid()
	
	if text_data[1] is not None:
		setText(plt.gca(), text_data[1], set_invisible=False)

	plt.colorbar()
	
	if titles is not None:
		plt.title(titles[1], fontsize=fontsize)

	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==

# This is the same as above, but only the first panel uses Basemap 
def plot1by2Basemap1(data, X, Y, lat_0, lon_0, contour=[None,None], contourlevels=None, lws=[1,1], figsize=(7.5,3), titles=None, fontsize=14, tfontsize=9., mesh=False, cmaps=['jet', 'jet'], vmin=None, vmax=None, text_data=[None,None], xlabels=[None, None], ylabels=[None, None], grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, isf=[None, None], labelData=None):
	
	from mpl_toolkits.basemap import Basemap
		
	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	fig = plt.figure(figsize=figsize, dpi=dpi)

	plt.subplot(121)
	ax = plt.gca()
	ax.patch.set_color('.5')

	X0, Y0 = m(X[0],Y[0])

	if mesh:
		m.pcolormesh(X0, Y0, data[0], cmap=cmaps[0], vmin=vmin[0], vmax=vmax[0])
	else:
		levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)
		a = m.contourf(X0, Y0, data[0], cmap=cmaps[0], levels=levels, extend='both')
		#a.cmap.set_under('w')	
		#a.set_clim(vmin[0], vmax[0])
		#a.cmap.set_over('w')	
		#plt.colorbar(a)

	plt.colorbar()

	if contour[0] is not None:
		plt.contour(X0, Y0, contour[0], levels=contourlevels[0], colors='k', linewidths=lws[0], linestyles='solid')

	if xlabels[0] is not None:
		plt.xlabel(xlabels[0], fontsize=fontsize)
	if ylabels[0] is not None:
		plt.ylabel(ylabels[0], fontsize=fontsize)
	
	if grid:
		plt.grid()

	if text_data[0] is not None:
		setText(plt.gca(), text_data[0], set_invisible=False)

	if titles is not None:
		plt.title(titles[0], fontsize=tfontsize)

	if yline[0] is not None:
		x0, y0 = m(yline[0], 0)
		plt.axvline(x=x0, color='k', linewidth=1.0, linestyle='--')
	
	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	if isf[0] is not None:
		extent = [X[0][0,0], X[0][0,-1], -Y[0][0,0], -Y[0][-1,0]]
		m.imshow(1-isf[0], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)
		
	#==
	# Second Panel
	#==
	
	plt.subplot(122)
	ax = plt.gca()
	ax.patch.set_color('.5')

	if mesh:
		plt.pcolormesh(X[1], Y[1], data[1], cmaps=cmaps[1], vmin=vmin[1], vmax=vmax[1])
	else:
		levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)
		plt.contourf(X[1], Y[1], data[1], cmap=cmaps[1], levels=levels)
	plt.colorbar()
		
	if contour[1] is not None:
		plt.contour(X[1], Y[1], contour[1], levels=contourlevels[1], colors='k', linewidths=lws[1], linestyles='solid')


	if xlabels[1] is not None:
		plt.xlabel(xlabels[1], fontsize=fontsize)
	if ylabels[1] is not None:
		plt.ylabel(ylabels[1], fontsize=fontsize)
		
	#plt.xlim(-75,-70); plt.ylim(-1000,0)
	if grid:
		plt.grid()
	
	if text_data[1] is not None:
		setText(plt.gca(), text_data[1], set_invisible=False)


	
	if titles is not None:
		plt.title(titles[1], fontsize=tfontsize)
	
	if isf[1] is not None:
		extent = [X[1][0,0], X[1][0,-1], -Y[1][0,0], -Y[1][-1,0]]
		m.imshow(1-isf[1], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)
	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==

def quiver1by1Basemap(u, v, X, Y, d, lat_0, lon_0, width=2.8e6, height=1.7e6, contourf=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1Basemap.mp4', show=True, dpi=200, text_data=None, parallels=None, meridians=None, isf=None, outline=None):

	from mpl_toolkits.basemap import Basemap
	
	#m = Basemap(width=width,height=height, projection='gnom',lat_0=lat_0,lon_0=lon_0)
	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	Xd = X[::d,::d]; Yd = Y[::d,::d]
	X,Y = m(X,Y)
	Xd,Yd = m(Xd,Yd)

	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()
	plt.gca().patch.set_color('.5')

	if contourf is not None:
		if mesh:
			m.pcolormesh(X, Y, contourf, cmap=cmap)
			m.colorbar()
		else:
			m.contour(X, Y, contourf, levels=[-3000, -2000, -1000])
			m.colorbar()

	if isinstance(u, (list)):
		q = m.quiver(Xd,Yd,u[0][::d, ::d],v[0][::d, ::d], color='r', scale=2)
		#plt.quiverkey(q,0.1, -0.1, 0.1, r'0.1 N m$^{-2}$',labelpos='W')
		q = m.quiver(Xd,Yd,u[1][::d, ::d],v[1][::d, ::d], color='k', scale=2)
		#plt.quiverkey(q,0.1, -0.12, 0.1, r'0.1 N m$^{-2}$',labelpos='W')
		plt.quiverkey(q,0.94, 0.03, 0.1,r'$0.1$ N m$^{-2}$',labelpos='W')
	else:
		u = u[::d, ::d]; v = v[::d, ::d]
		q = m.quiver(Xd,Yd,u,v)
		plt.quiverkey(q,0.94, 0.03, 0.1, r'$0.1$ N m$^{-2}$',labelpos='W')

	#m.fillcontinents(color='#cc9955', zorder = 0)

	if xlabel is not None:
		plt.xlabel(xlabel, labelpad=20)
	if ylabel is not None:
		plt.ylabel(ylabel, labelpad=35)

	if title is not None:
		plt.title(title) 

	if text_data is not None:
		setText(ax, text_data, i=0)

	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])


	if isf is not None:
		#x, y = m(isf[0], isf[1])  # transform coordinates
		#m.scatter(x, y, color='k', s=1)
		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]
		m.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);
		#m.contour(X, Y, isf, levels=[1], colors='k')

	if outline is not None:
		lons = outline[0]
		lats = outline[1]
		lons, lats = m(lons, lats)

		m.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')
		m.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')
		m.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')
		m.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')

	plt.tight_layout()

	#==
	
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()	


#==

def animate1by1quiverBasemap(u, v, Xd, Yd, lat_0, lon_0, width=2.e6, height=1.7e6, contourf=None, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, parallels=None, meridians=None):

	from mpl_toolkits.basemap import Basemap

	m = Basemap(width=width,height=height, projection='gnom',lat_0=lat_0,lon_0=lon_0)

	X,Y = m(X,Y)
	Xd,Yd = m(Xd,Yd)

	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()
	#plt.gca().patch.set_color('.25')

	if contourf is not None:

		# Plot
		m.contourf(X, Y, contourf)
		plt.colorbar()

	q = m.quiver(Xd,Yd,u[0],v[0])

	#m.fillcontinents(color='#cc9955', zorder = 0)
	
	if parallels is not None:
		# Latitudes
		parallels = m.drawparallels(parallels)
		m.drawparallels(parallels,labels=[True,False,False,True])
	if meridians is not None:
		# Longitudes
		meridians = m.drawmeridians(meridians)
		m.drawmeridians(meridians,labels=[True,False,False,True])

	#plt.xlabel('Longitude')
	#plt.ylabel('Latitude')
	plt.quiverkey(q,0.9, 0.05, 6,'6 m/s',labelpos='W')

	if title is not None:
		plt.title(title) 

	if text_data is not None:
		setText(ax, text_data, i=0)

	def animate(i):
		q.set_UVC(u[i], v[i])
		if text_data is not None:
			setText(ax, text_data, i=i, set_invisible=True)
		
		# End if mesh	

	anim = animation.FuncAnimation(fig, animate, interval=50, frames=u.shape[0])
	plt.draw()

	#==
	
	if save:
		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)
		
	if show:
		plt.show()	

	
