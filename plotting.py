# plotting.py

import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec

from plotting_tools import setText, getContourfLevels
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
	plt.gca().patch.set_color('.25')

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


def animateLine(data, X=None, figsize=(5,4), title='', labels=None, fontsize=14, vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animateLine.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None):

	# Make animation
	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()

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

def quiver1by1Basemap(u, v, Xd, Yd, lat_0, lon_0, width=2.e6, height=1.7e6, contourf=None, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1Basemap.mp4', show=True, dpi=200, text_data=None, parallels=None, meridians=None):

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

	
