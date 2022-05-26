# plotting.py

import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

from plotting_tools import setText, getContourfLevels, maskBathyXY, maskDraftXY

#==========================================================
	
#==

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

def plot1by1(data, X=None, Y=None, contour=None, contourlevels=None, figsize=(5,4), title=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabel=None, ylabel=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None):
	
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

	if yline is not None:
		plt.plot(yline, Y, color='k', linewidth=1.2)
		plt.axvline(x=X[X.shape[0]//2-1], color='k', linewidth=1.0, linestyle='--')

	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
	
#==

def quiver1by1(u, v, Xd, Yd, C=None, ccmap='bwr', contour=None, X=None, Y=None, contourf=True, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1.mp4', show=True, dpi=200, text_data=None):

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
		cax = ax.quiver(Xd, Yd, u, v, C, cmap=ccmap)
		plt.colorbar(cax, ax=ax)
	else:
		cax = ax.quiver(Xd, Yd, u, v)

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

def quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, figsize=(8,3), title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, save=False, outpath='', outname='quiver1by2.mp4', show=True, dpi=200, text_data=None, parallels=None, meridians=None, width_ratios=None, labelData=None):
	
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

	if text_data is not None:
		setText(ax, text_data[0], set_invisible=False)
	
	if title is not None:
		plt.title(title[0], fontsize=fontsize)

	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	if labelData is not None:
		for li in labelData[0]:
			plt.scatter(li['x'][0], li['x'][1], s=1, color='r')
			plt.annotate(li['t'], li['x'])

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

	if text_data is not None:
		setText(ax, text_data[1], set_invisible=False)
	
	if title is not None:
		plt.title(title[1], fontsize=fontsize)

	ax.set_aspect(1)
	if parallels is not None:
		# Latitudes
		m.drawparallels(parallels,labels=[False,False,False,False], color='k', linewidth=0.5,dashes=[2,1])
	if meridians is not None:
		# Longitudes
		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

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
		
def plot1by2(data, X=None, Y=None, figsize=(9,4), titles=None, fontsize=14, mesh=False, cmap='jet', vmin=[None,None], vmax=[None,None], text_data=[None,None], xlabels=None, ylabels=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by2.png', show=True, dpi=200):
	
	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	plt.subplot(121)
	plt.gca().patch.set_color('.5')

	if mesh:
		if X is not None and Y is not None:
			plt.pcolormesh(X[0], Y[0], data[0], cmap=cmap, vmin=vmin[0], vmax=vmax[0])
		else: 
			plt.pcolormesh(data[0], cmap=cmap, vmin=vmin[0], vmax=vmax[0])
	else:
		levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)
		if X is not None and Y is not None:
			plt.contourf(X[0], Y[0], data[0], cmap=cmap, levels=levels)
		else: 
			plt.contourf(data[0], cmap=cmap, levels=levels)
			
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
	

	#==
	# Second Panel
	#==
	
	plt.subplot(122)
	plt.gca().patch.set_color('.25')
	
	if mesh:
		if X is not None and Y is not None:
			plt.pcolormesh(X[1], Y[1], data[1], cmap='jet', vmin=vmin[1], vmax=vmax[1])
		else: 
			plt.pcolormesh(data[1], cmap=cmap,  vmin=vmin[1], vmax=vmax[1])
	else:
		levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)
		if X is not None and Y is not None:
			plt.contourf(X[1], Y[1], data[1], cmap=cmap, levels=levels)
		else: 
			plt.contourf(data[1], cmap=cmap, levels=levels)
	
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
		
	#==
	
	plt.tight_layout()
	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()

#==
	
		
def plot1by3(data, X=[None]*3, Y=[None]*3, contour=[None]*3, contourlevels=[None]*3, figsize=(11,3), titles=[None]*3, fontsize=14, mesh=False, cmaps=['jet']*3, vmin=[None]*3, vmax=[None]*3, text_data=[None]*3, xlabels=[None]*3, ylabels=[None]*3, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3):
	
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

		if xlabels[pi] is not None:
			plt.xlabel(xlabels[pi], fontsize=fontsize)
		if ylabels[pi] is not None:
			plt.ylabel(ylabels[pi], fontsize=fontsize)

		if xticks[pi] is not None:
			if xticksvis[pi]:				
				plt.xticks(xticks[pi])
			else:			
				plt.xticks(xticks[pi], labels='')
		if yticks[pi] is not None:
			if yticksvis[pi]:				
				plt.yticks(yticks[pi])
			else:			
				plt.yticks(yticks[pi], labels='')

		if grid:
			plt.grid()

		if text_data[pi] is not None:
			setText(plt.gca(), text_data[pi], set_invisible=False)

		if titles[pi] is not None:
			plt.title(titles[pi], fontsize=fontsize)
		
	box = ax1.get_position()
	box.x0 = box.x0 + 0.05
	ax1.set_position(box)


	#==

	#plt.tight_layout()
	
	if save:
		plt.savefig(outpath + outname, bbox_inches="tight")
	if show:
		plt.show()

	
		
def plot1by3_(data, X=[None]*3, Y=[None]*3, contour=[None]*3, contourlevels=[None]*3, figsize=(11,3), titles=[None]*3, fontsize=14, mesh=False, cmaps=['jet']*3, vmin=[None]*3, vmax=[None]*3, text_data=[None]*3, xlabels=[None]*3, ylabels=[None]*3, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3):
	
	fig = plt.figure(figsize=figsize, dpi=dpi)

	if width_ratios is not None:
		gs = gridspec.GridSpec(ncols=3, nrows=1, figure=fig, width_ratios=width_ratios)

	for pi in range(3):

		if width_ratios is not None:
			ax = fig.add_subplot(gs[0,pi])
		else:
			plt.subplot(1,3,pi+1)
			ax = plt.gca()

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

		if xlabels[pi] is not None:
			plt.xlabel(xlabels[pi], fontsize=fontsize)
		if ylabels[pi] is not None:
			plt.ylabel(ylabels[pi], fontsize=fontsize)

		if xticks[pi] is not None:
			if xticksvis[pi]:				
				plt.xticks(xticks[pi])
			else:			
				plt.xticks(xticks[pi], labels='')
		if yticks[pi] is not None:
			if yticksvis[pi]:				
				plt.yticks(yticks[pi])
			else:			
				plt.yticks(yticks[pi], labels='')

		if grid:
			plt.grid()

		if text_data[pi] is not None:
			setText(plt.gca(), text_data[pi], set_invisible=False)

		if titles[pi] is not None:
			plt.title(titles[pi], fontsize=fontsize)
	
		pos = [[0,0,0.33,0.33], [0.33, 0.33, 0.33, 0.33], [0.66, 0.66, 0.33, 0.33]]	

		ax.set_position(pos[pi])			

	#==

	plt.tight_layout()

	if save:
		plt.savefig(outpath + outname)
		
	if show:
		plt.show()
	

#==

def animate1by1(data, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, contour=None):

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
			plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)

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
			plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)

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
def plot1by2Basemap1(data, X, Y, lat_0, lon_0, contour=[None,None], contourlevels=None, figsize=(7.5,3), titles=None, fontsize=14, mesh=False, cmaps=['jet', 'jet'], vmin=None, vmax=None, text_data=[None,None], xlabels=[None, None], ylabels=[None, None], grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, labelData=None):
	
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
		a = m.contourf(X0, Y0, data[0], cmap=cmaps[0], levels=levels, extend='max')
		#a.cmap.set_under('w')	
		#a.set_clim(vmin[0], vmax[0])
		#a.cmap.set_over('w')	
		#plt.colorbar(a)

	plt.colorbar()

	if contour[0] is not None:
		plt.contour(X0, Y0, contour, levels=contourlevels, colors='k', linewidths=0.2, linestyles='solid')

	if xlabels[0] is not None:
		plt.xlabel(xlabels[0], fontsize=fontsize)
	if ylabels[0] is not None:
		plt.ylabel(ylabels[0], fontsize=fontsize)
	
	if grid:
		plt.grid()

	if text_data[0] is not None:
		setText(plt.gca(), text_data[0], set_invisible=False)

	if titles is not None:
		plt.title(titles[0], fontsize=fontsize)

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
	
	if xlabels[1] is not None:
		plt.xlabel(xlabels[1], fontsize=fontsize)
	if ylabels[1] is not None:
		plt.ylabel(ylabels[1], fontsize=fontsize)
		
	#plt.xlim(-75,-70); plt.ylim(-1000,0)
	if grid:
		plt.grid()
	
	if text_data[1] is not None:
		setText(plt.gca(), text_data[1], set_invisible=False)

	plt.colorbar()
	
	if titles is not None:
		plt.title(titles[1], fontsize=fontsize)
		
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
		plt.quiverkey(q,0.1, -0.12, 0.1, r'0.1 N m$^{-2}$',labelpos='W')

	else:
		u = u[::d, ::d]; v = v[::d, ::d]
		q = m.quiver(Xd,Yd,u,v)
		plt.quiverkey(q,0.0, -0.1, 0.1,'0.1 N m^-2',labelpos='W')
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

	
