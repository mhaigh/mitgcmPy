# plotting.py

import sys

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
	
#==========================================================

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
		
def plot1by2(data, X=None, Y=None, figsize=(9,4), titles=None, fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabels=None, ylabels=None, save=False, outpath='', outname='plot1by2.png', show=True, dpi=200):
	
	fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

	plt.subplot(121)
	plt.gca().patch.set_color('.25')
	
	if mesh:
		if X is not None and Y is not None:
			plt.pcolormesh(X[0], Y[0], data[0], cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			plt.pcolormesh(data[0], cmap=cmap, vmin=vmin, vmax=vmax)
	else:
		if X is not None and Y is not None:
			plt.contourf(X[0], Y[0], data[0], cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			plt.contourf(data[0], cmap=cmap, vmin=vmin, vmax=vmax)
			
	if xlabels is not None:
		plt.xlabel(xlabels[0], fontsize=fontsize)
	if ylabels is not None:
		plt.ylabel(ylabels[0], fontsize=fontsize)
		
	plt.colorbar()
	
	if titles is not None:
		plt.title(titles[0], fontsize=fontsize)
	
	#==
	
	plt.subplot(122)
	plt.gca().patch.set_color('.25')
	
	if mesh:
		if X is not None and Y is not None:
			plt.pcolormesh(X[1], Y[1], data[1], cmap='jet', vmin=vmin, vmax=vmax)
		else: 
			plt.pcolormesh(data[1], cmap=cmap,  vmin=vmin, vmax=vmax)
	else:
		nlevels = 9
		if X is not None and Y is not None:
			plt.contourf(X[1], Y[1], data[1], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))
		else: 
			plt.contourf(data[1], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))
	
	if xlabels is not None:
		plt.xlabel(xlabels[1], fontsize=fontsize)
	if ylabels is not None:
		plt.ylabel(ylabels[1], fontsize=fontsize)
		
				
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

def animate1by1(data, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1):

	# Make animation
	fig = plt.figure(figsize=figsize, dpi=dpi)
	ax = fig.add_subplot()
	plt.gca().patch.set_color('.25')

	if mesh:
		if X is not None and Y is not None:
			cax = ax.pcolormesh(X, Y, data[0, ], cmap=cmap, vmin=vmin, vmax=vmax)
		else: 
			cax = ax.pcolormesh(data[0, ], cmap=cmap, vmin=vmin, vmax=vmax)
		def animate(i):
			cax.set_array(data[i,].flatten())
		# End if mesh	

	else:
		nlevels = 9
		if X is not None and Y is not None:
			cax = ax.contourf(X, Y, data[0, ], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))
		else: 
			cax = ax.contourf(data[0, ], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))
		def animate(i):
			ax.contourf(X, Y, data[i,], levels=np.linspace(vmin, vmax, nlevels))
			
	#==
	
	fig.colorbar(cax, ax=ax)
	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title, fontsize=fontsize)
	plt.grid() 

	anim = animation.FuncAnimation(fig, animate, interval=50, frames=data.shape[0])
	plt.tight_layout()
	plt.draw()
	
	if save:
		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)
		
	if show:
		plt.show()
		
	
	
