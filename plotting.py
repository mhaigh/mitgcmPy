

import sys



import numpy as np



import matplotlib.pyplot as plt

import matplotlib.animation as animation

import matplotlib.gridspec as gridspec

import matplotlib.patches as patches



from mpl_toolkits.axes_grid1 import make_axes_locatable

from mpl_toolkits.mplot3d import axes3d

from mpl_toolkits.basemap import Basemap

	

from plotting_tools import setText, getContourfLevels, maskBathyXY, maskDraftXY, makeList, doTicks, doTitle, doLabels, doStippling, doBox, makeBathyContour, doAntarcticInset



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



def plot1by1(data, X=None, Y=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', figsize=(5,4), title=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabel=None, ylabel=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, hlines=None, xmin=[0], xmax=[1], vlines=None, xlim=None, ylim=None, stippling=None, stipData=[0.05, 4, 4, .1], box=None):

	

	

	if contourLevels2 == None:

		contourLevels2 = contourLevels

	if contourLevels3 == None:

		contourLevels3 = contourLevels

		

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

		plt.contour(X, Y, contour, colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels)

	if contour2 is not None:

		plt.contour(X, Y, contour2, colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2)

	if contour3 is not None:

		plt.contour(X, Y, contour3, colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3)			

	

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

		

	if stippling is not None:

		doStippling(stippling, X=X, Y=Y, dataLim=stipData[0], dx=stipData[1], dy=stipData[2], s=stipData[3])

		

	doBox(box, X, Y)

		

	#==

	

	plt.tight_layout()

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()

	

#==



def line1byN(data, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis, fontsize=14, figsize=(9,3), text_data=None, loc=2, save=False, outname='line1byN.png', show=True):



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

		

		if text_data is not None:

			setText(plt.gca(), text_data, set_invisible=False, i=pi)

			

		plt.title(titles[pi])

		plt.xlabel(xlabels[pi], fontsize=fontsize)

		plt.ylabel(ylabels[pi], fontsize=fontsize-2)

		plt.grid()

	

	plt.legend(loc=loc) 

	

	plt.tight_layout()

		

	if save:

		plt.savefig(outname)

		

	if show:

		plt.show()



#==



def lineBar1byN(data, labels, colours, xlims, ylims, titles, xlabels, ylabels, xticks, xticksvis, yticks, yticksvis, t0, fontsize=14, figsize=(9,3), width_ratios=[1,0.5], save=False, outname='line1byN.png', show=True):



	N = len(data)

	if N == 3:

		figsize = (9,8)

	elif N == 4:

		figsize = (9,10.5)

		

	fig = plt.figure(figsize=figsize)

	gs = gridspec.GridSpec(ncols=2, nrows=N, figure=fig, width_ratios=width_ratios)

	

	for pi in range(N):

		#plt.subplot(N,2,2*pi+1)

		fig.add_subplot(gs[pi,0])

		

		ni = len(data[pi])

		for i in range(ni):

			plt.plot(data[pi][i], label=labels[pi][i], color=colours[pi][i])

			

		plt.xlim(xlims)

		plt.ylim(ylims)

		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])

		

		plt.title(titles[pi])

		plt.xlabel(xlabels[pi], fontsize=fontsize)

		plt.ylabel(ylabels[pi], fontsize=fontsize-2)

		plt.grid()

		

		if pi == 3:

			plt.legend(loc=2) 

		

		#==

		

		

		fig.add_subplot(gs[pi,1])

		

		for i in range(ni):

			plt.bar(i, np.mean(data[pi][ni-i-1][t0:]), label=labels[pi][ni-i-1], color=colours[pi][ni-i-1])

			

		#plt.xlim([-.5, 4.5])

		plt.grid()

		plt.ylim(ylims)

		doTicks(np.linspace(0.5,4.5,5), False, yticks[pi], False)

		

		#==

	

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



def quiver1byN(u, v, Xd, Yd, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, contour2=None, contourLevels2=None, figsize=None, title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, text_data=None, width_ratios=None, labelData=None, cbar=True, grid=True, xticks=None, xticksvis=True, yticks=None, yticksvis=True, scale=2, qs=0.1, qunits='m/s'):



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

	

	qs = makeList(qs, N)

	scale = makeList(scale, N)

	

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



		ax.patch.set_color('.6')



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

				

		if contour2[pi] is not None:

			plt.contour(X[pi], Y[pi], contour2[pi], colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2[pi])

			



		if C is not None:

			cax = ax.quiver(Xd[pi], Yd[pi], u[pi], v[pi], C[pi], cmap=ccmap, scale=scale[pi])

			if width_ratios is None:

				plt.colorbar(cax, ax=ax)

		else:

			cax = plt.quiver(Xd[pi], Yd[pi], u[pi], v[pi], scale=scale[pi])

		ax.quiverkey(cax, 0.15, 0.03, qs[pi], str(qs[pi]) + ' ' + qunits, labelpos='N', coordinates='axes')

				

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

		plt.savefig(outpath + outname, bbox_inches="tight")



	if show:

		plt.show()



	plt.close()



#==



def quiver1byN_Basemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, cbarTicks=None, contourfNlevels=9, X=None, Y=None, mesh=False, isf=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', DASH_NEG_CONTOURS=[False]*4, figsize=None, titles=None, fstitle=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, parallels=None, meridians=None, save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, text_data=None, width_ratios=1, labelData=None, cbar=True, cbarShrink=1., grid=True, scale=2, qs=0.1, qunits='m/s', labelpos='N', landscape=True, extend='neither', width=0.003, headwidth=5., headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03):



	N = len(u)



	Xd = makeList(Xd, N)

	Yd = makeList(Yd, N)



	if contourf is not None:

		X = makeList(X, N)

		Y = makeList(Y, N)



	vmin = makeList(vmin, N)

	vmax = makeList(vmax, N)



	contour = makeList(contour, N)

	contour2 = makeList(contour2, N)

	contour3 = makeList(contour3, N)

	contour4 = makeList(contour4, N)

	

	cmap = makeList(cmap, N)

	cbar = makeList(cbar, N)

	extend = makeList(extend, N)

	cbarTicks = makeList(cbarTicks, N)

	

	grid = makeList(grid, N)

	titles = makeList(titles, N)

	if fstitle is None:

		fstitle = [fontsize]*N

		

	width_ratios = makeList(width_ratios, N)



	xlabel = makeList(xlabel, N)

	ylabel = makeList(ylabel, N)



	qcolor = makeList(qcolor, N)

	qs = makeList(qs, N)

	scale = makeList(scale, N)

	labelpos = makeList(labelpos, N)

	

	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	

	if figsize is None:

		if N == 2:

			figsize = (8,3)

		elif N == 3:

			figsize = (11,3)

		else:

			figsize = (8,3)



	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)





	if landscape:

		gs = gridspec.GridSpec(ncols=N, nrows=1, figure=fig, width_ratios=width_ratios)

	else:

		gs = gridspec.GridSpec(ncols=1, nrows=N, figure=fig, height_ratios=width_ratios)



	for pi in range(N):



		if landscape:

			ax = fig.add_subplot(gs[0,pi])

		else:

			ax = fig.add_subplot(gs[pi,0])

			

		ax.patch.set_color('.6')



		X0, Y0 = m(X[pi],Y[pi])

		Xd0, Yd0 = m(Xd[pi],Yd[pi])

	

		if contourf is not None:

			if mesh:

				plt.pcolormesh(X0, Y0, contourf[pi], vmin=vmin[pi], vmax=vmax[pi], cmap=cmap[pi])		

			else:

				levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)

				plt.contourf(X0, Y0, contourf[pi], cmap=cmap[pi], levels=levels, extend=extend[pi])



			if cbar[pi]:

				cbar_ = plt.colorbar(ticks=cbarTicks[pi], shrink=cbarShrink)

				cbar_.ax.tick_params(labelsize=fontsize)

		

		if contour[pi] is not None:

			if DASH_NEG_CONTOURS[0]:

				ls = np.where(np.array(contourLevels[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour[pi], colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels[pi])

		if contour2[pi] is not None:

			if DASH_NEG_CONTOURS[1]:

				ls2 = np.where(np.array(contourLevels2[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour2[pi], colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2[pi],zorder=14)

			

		if contour3[pi] is not None:

			if DASH_NEG_CONTOURS[2]:

				ls3 = np.where(np.array(contourLevels3[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour3[pi], colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3[pi], zorder=14)

		if contour4[pi] is not None:

			if DASH_NEG_CONTOURS[3]:

				ls4 = np.where(np.array(contourLevels4[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour4[pi], colors=cc4, linestyles=ls4, linewidths=lw4, levels=contourLevels4[pi], zorder=14)



		if isf[pi] is not None:

			extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]

			m.imshow(1-isf[pi], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)

		

		if u[pi] is not None:

			if C is not None:

				cax = ax.quiver(Xd0, Yd0, u[pi], v[pi], C[pi], cmap=ccmap, scale=scale[pi])

				if width_ratios is None:

					plt.colorbar(cax, ax=ax)

			else:

				cax = plt.quiver(Xd0, Yd0, u[pi], v[pi], scale=scale[pi], width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, color=qcolor[pi])

			ax.quiverkey(cax, qlabelx, qlabely, qs[pi], str(qs[pi]) + ' ' + qunits, labelpos=labelpos[pi], coordinates='axes', labelcolor=qcolor[pi])

				

		doLabels(xlabel[pi], ylabel[pi], fontsize=fontsize)

		doTitle(titles[pi], fontsize=fstitle[pi])



		if text_data is not None:

			setText(ax, text_data[pi], set_invisible=False)

		

		if labelData is not None:

			for li in labelData[pi]:

				plt.scatter(li['x'][0], li['x'][1], s=1, color='r')

				plt.annotate(li['t'], li['x'])



		#ax.set_aspect('equal')



		if grid[pi]:

			plt.grid()

		

		ax.set_aspect('equal')

		if parallels[pi] is not None:

			# Latitudes

			m.drawparallels(parallels[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

		if meridians[pi] is not None:

			# Longitudes

			m.drawmeridians(meridians[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)



	#==

	

	#fig.subplots_adjust(wspace=-1, hspace=0)



	plt.tight_layout()



	if save:

		plt.savefig(outpath + outname, bbox_inches="tight")



	if show:

		plt.show()



	plt.close()



#==



def quiverMbyN_Basemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, cbarTicks=None, contourfNlevels=9, X=None, Y=None, mesh=False, isf=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', DASH_NEG_CONTOURS=[False]*4, figsize=None, titles=None, fstitle=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, parallels=None, meridians=None, save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, text_data=None, width_ratios=1, labelData=None, cbar=True, cbarShrink=1., grid=True, scale=2, qs=0.1, qunits='m/s', labelpos='E', landscape=True, extend='neither', width=0.003, headwidth=5., headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03, w_pad=0, h_pad=0):



	N = len(u)

	M = len(u[0])

	

	Xd = makeList(Xd, M, N)

	Yd = makeList(Yd, M, N)



	if contourf is not None:

		X = makeList(X, M, N)

		Y = makeList(Y, M, N)



	vmin = makeList(vmin, M, N)

	vmax = makeList(vmax, M, N)



	contour = makeList(contour, M, N)

	contour2 = makeList(contour2, M, N)

	contour3 = makeList(contour3, M, N)

	contour4 = makeList(contour4, M, N)

	

	cmap = makeList(cmap, M, N)

	cbar = makeList(cbar, M, N)

	extend = makeList(extend, M, N)

	cbarTicks = makeList(cbarTicks, M, N)

	contourfNlevels = makeList(contourfNlevels, M, N)

	

	grid = makeList(grid, M, N)

	titles = makeList(titles, M, N)

	if fstitle is None:

		fstitle = makeList(fontsize, M, N)

		

	width_ratios = makeList(width_ratios, M)



	xlabel = makeList(xlabel, M, N)

	ylabel = makeList(ylabel, M, N)



	qcolor = makeList(qcolor, M, N)

	qs = makeList(qs, M, N)

	scale = makeList(scale, M, N)

	labelpos = makeList(labelpos, M, N)

	

	m = Basemap(llcrnrlon=X[0][0][0,0],llcrnrlat=Y[0][0][0,0],urcrnrlon=X[0][0][-1,-1],urcrnrlat=Y[0][0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	

	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)

	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)

	

	for col in range(M):

		for row in range(N):



			ax = fig.add_subplot(gs[row,col])

			ax.patch.set_color('.6')



			X0, Y0 = m(X[row][col],Y[row][col])

			Xd0, Yd0 = m(Xd[row][col],Yd[row][col])

		

			if contourf is not None:

				if mesh:

					plt.pcolormesh(X0, Y0, contourf[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap[row][col])		

				else:

					levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels[row][col])

					plt.contourf(X0, Y0, contourf[row][col], cmap=cmap[row][col], levels=levels, extend=extend[row][col])



				if cbar[row][col]:

					cbar_ = plt.colorbar(ticks=cbarTicks[row][col], shrink=cbarShrink)

					cbar_.ax.tick_params(labelsize=fontsize)

			

			if contour[row][col] is not None:

				if DASH_NEG_CONTOURS[0]:

					ls = np.where(np.array(contourLevels[row][col]) > 0, "-", "--")

				plt.contour(X0, Y0, contour[row][col], colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels[row][col])

			if contour2[row][col] is not None:

				if DASH_NEG_CONTOURS[1]:

					ls2 = np.where(np.array(contourLevels2[row][col]) > 0, "-", "--")

				plt.contour(X0, Y0, contour2[row][col], colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2[row][col],zorder=14)

			if contour3[row][col] is not None:

				if DASH_NEG_CONTOURS[2]:

					ls3 = np.where(np.array(contourLevels3[row][col]) > 0, "-", "--")

				plt.contour(X0, Y0, contour3[row][col], colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3[row][col], zorder=14)

			if contour4[row][col] is not None:

				if DASH_NEG_CONTOURS[3]:

					ls4 = np.where(np.array(contourLevels4[row][col]) > 0, "-", "--")

				plt.contour(X0, Y0, contour4[row][col], colors=cc4, linestyles=ls4, linewidths=lw4, levels=contourLevels4[row][col], zorder=14)



			if isf[row][col] is not None:

				extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]

				m.imshow(1-isf[row][col], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)

			

			if u[row][col] is not None:

				if C is not None:

					cax = ax.quiver(Xd0, Yd0, u[row][col], v[row][col], C[row][col], cmap=ccmap, scale=scale[row][col])

					if width_ratios is None:

						plt.colorbar(cax, ax=ax)

				else:

					cax = plt.quiver(Xd0, Yd0, u[row][col], v[row][col], scale=scale[row][col], width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, color=qcolor[row][col])

				ax.quiverkey(cax, qlabelx, qlabely, qs[row][col], str(qs[row][col]) + ' ' + qunits, labelpos=labelpos[row][col], coordinates='axes', labelcolor=qcolor[row][col])

					

			doLabels(xlabel[row][col], ylabel[row][col], fontsize=fontsize)

			doTitle(titles[row][col], fontsize=fstitle[row][col])



			if text_data is not None:

				setText(ax, text_data[row][col], set_invisible=False)

			

			if labelData is not None:

				for li in labelData[row][col]:

					plt.scatter(li['x'][0], li['x'][1], s=1, color='r')

					plt.annotate(li['t'], li['x'])



			#ax.set_aspect('equal')



			if grid[row][col]:

				plt.grid()

			

			ax.set_aspect('equal')

			if parallels[row][col] is not None:

				# Latitudes

				m.drawparallels(parallels[row][col],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

			if meridians[row][col] is not None:

				# Longitudes

				m.drawmeridians(meridians[row][col],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)



		#==

	

	#fig.subplots_adjust(wspace=-1, hspace=0)



	plt.tight_layout(w_pad=w_pad, h_pad=h_pad)



	if save:

		plt.savefig(outpath + outname, bbox_inches="tight")



	if show:

		plt.show()



	plt.close()





#==



def quiver1byNBasemapAndSections(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, isf=None, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', DASH_NEG_CONTOURS=[False]*4, titles=None, fstitle=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabel=None, ylabel=None, parallels=None, meridians=None, text_data=None,  labelData=None, cbar=True,  cbarTicks=None, cbarLabels=None, cbarShrink=1., grid=True, scale=2, qs=0.1, qunits='m/s', labelpos='N', extend='neither', width=0.003, headwidth=5., headlength=5., headaxislength=3, qcolor='k', qlabelx=0.3, qlabely=0.03, \

sec_data=[None], sec_X=[None]*3, sec_Y=[None]*3, sec_contour=[None]*3, sec_contourLevels=[None]*3, sec_DASH_NEG_CONTOURS=False, sec_titles=[None]*3, sec_fstitle=None, sec_fontsize=14, sec_mesh=False, sec_cmaps=['jet']*3, sec_vmin=[None]*3, sec_vmax=[None]*3, sec_extend=['neither']*3, sec_text_data=[None]*3, sec_xlabels=[None]*3, sec_ylabels=[None]*3, sec_grid=True, sec_gridls='solid', sec_gridlw=0.8, sec_gridlc='k', sec_contourfNlevels=9, sec_clabel=[False]*3, sec_cbar=[True]*3, sec_cbarTicks=[None]*3, sec_xticks=[None]*3, sec_yticks=[None]*3, sec_xticksvis=[True]*3, sec_yticksvis=[True]*3, \

figsize=(8,8), save=False, outpath='', outname='quiver1byN.png', show=True, dpi=200, width_ratios=1):



	#=

	#=

	

	N = len(u)

	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)

	gs = gridspec.GridSpec(ncols=2, nrows=N, figure=fig, width_ratios=width_ratios)



	#=

	#=

	

	# Sections column

	

	sec_vmin = makeList(sec_vmin, 3)

	sec_vmax = makeList(sec_vmax, 3)

	sec_cmaps = makeList(sec_cmaps, 3)



	sec_contourfNlevels = makeList(sec_contourfNlevels, 3)

	

	sec_X = makeList(sec_X, 3)

	sec_Y = makeList(sec_Y, 3)	

	

	if sec_fstitle is None:

		sec_fstitle = [sec_fontsize]*3

	else:

		sec_fstitle = makeList(sec_fstitle, 3)



	for pi in range(N):



		ax = fig.add_subplot(gs[pi,0])

		ax.patch.set_color('.6')



		if sec_mesh:

			if sec_X[pi] is not None and sec_Y[pi] is not None:

				plt.pcolormesh(sec_X[pi], sec_Y[pi], sec_data[pi], cmap=sec_cmaps[pi], vmin=sec_vmin[pi], vmax=sec_vmax[pi])

			else: 

				plt.pcolormesh(sec_data[pi], cmap=sec_cmaps[pi], vmin=sec_vmin[pi], vmax=sec_vmax[pi])

		else:

			levels = getContourfLevels(sec_vmin[pi], sec_vmax[pi], sec_contourfNlevels[pi])

			if sec_X[pi] is not None and sec_Y[pi] is not None:

				plt.contourf(sec_X[pi], sec_Y[pi], sec_data[pi], cmap=sec_cmaps[pi], levels=levels, extend=sec_extend[pi])

			else: 

				plt.contourf(sec_data[pi], cmap=sec_cmaps[pi], levels=levels, extend=sec_extend[pi])



		if sec_cbar[pi]:

			cbar_ = plt.colorbar(ticks=sec_cbarTicks[pi])

			cbar_.formatter.set_useOffset(False)

			

		if sec_contour[pi] is not None:

			if sec_DASH_NEG_CONTOURS:

				linestyles = np.where(np.array(sec_contourLevels[pi]) > 0, "-", "--")

			else:

				linestyles = 'solid'

			if sec_X[pi] is not None and sec_Y[pi] is not None:

				CS = plt.contour(sec_X[pi], sec_Y[pi], sec_contour[pi], levels=sec_contourLevels[pi], colors='k', linestyles=linestyles, linewidths=0.8)

			else:

				CS = plt.contour(sec_X[pi], sec_Y[pi], sec_contour[pi], levels=sec_contourLevels[pi], colors='k', linestyles=linestyles, linewidths=0.8)

			if sec_clabel[pi]:

				plt.clabel(CS, fmt='%2.1f', colors='k', fontsize=7)



		doLabels(sec_xlabels[pi], sec_ylabels[pi], fontsize=sec_fontsize)

		doTicks(sec_xticks[pi], sec_xticksvis[pi], sec_yticks[pi], sec_yticksvis[pi])

		doTitle(sec_titles[pi], fontsize=sec_fstitle[pi])

	

		if grid:

			plt.grid(linestyle=sec_gridls, color=sec_gridlc, linewidth=sec_gridlw)



		if sec_text_data[pi] is not None:

			setText(ax, sec_text_data[pi], set_invisible=False)



	#=

	#=

	

	# Basemap column



	Xd = makeList(Xd, N)

	Yd = makeList(Yd, N)



	if contourf is not None:

		X = makeList(X, N)

		Y = makeList(Y, N)



	vmin = makeList(vmin, N)

	vmax = makeList(vmax, N)

	contourfNlevels = makeList(contourfNlevels, N)

	

	contour = makeList(contour, N)

	contour2 = makeList(contour2, N)

	contour3 = makeList(contour3, N)

	contour4 = makeList(contour4, N)

	

	cmap = makeList(cmap, N)

	cbar = makeList(cbar, N)

	extend = makeList(extend, N)

	cbarTicks = makeList(cbarTicks, N)

	cbarLabels = makeList(cbarLabels, N)

	

	grid = makeList(grid, N)

	titles = makeList(titles, N)

	if fstitle is None:

		fstitle = [fontsize]*N

		

	width_ratios = makeList(width_ratios, N)



	xlabel = makeList(xlabel, N)

	ylabel = makeList(ylabel, N)



	qcolor = makeList(qcolor, N)

	qs = makeList(qs, N)

	scale = makeList(scale, N)

	labelpos = makeList(labelpos, N)

	

	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)

	for pi in range(N):



		ax = fig.add_subplot(gs[pi,1])

			

		ax.patch.set_color('.6')



		X0, Y0 = m(X[pi],Y[pi])

		Xd0, Yd0 = m(Xd[pi],Yd[pi])

	

		if contourf is not None:

			if mesh:

				plt.pcolormesh(X0, Y0, contourf[pi], vmin=vmin[pi], vmax=vmax[pi], cmap=cmap[pi])		

			else:

				levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels[pi])

				plt.contourf(X0, Y0, contourf[pi], cmap=cmap[pi], levels=levels, extend=extend[pi])



			if cbar[pi]:

				cbar_ = plt.colorbar(ticks=cbarTicks[pi], shrink=cbarShrink)

				cbar_.ax.tick_params(labelsize=fontsize)

				cbar_.set_label(cbarLabels[pi], rotation=270)

				cbar_.ax.get_yaxis().labelpad = 13

				

		if contour[pi] is not None:

			if DASH_NEG_CONTOURS[0]:

				ls = np.where(np.array(contourLevels[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour[pi], colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels[pi])

		if contour2[pi] is not None:

			if DASH_NEG_CONTOURS[1]:

				ls2 = np.where(np.array(contourLevels2[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour2[pi], colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2[pi],zorder=14)

			

		if contour3[pi] is not None:

			if DASH_NEG_CONTOURS[2]:

				ls3 = np.where(np.array(contourLevels3[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour3[pi], colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3[pi], zorder=14)

		if contour4[pi] is not None:

			if DASH_NEG_CONTOURS[3]:

				ls4 = np.where(np.array(contourLevels4[pi]) > 0, "-", "--")

			plt.contour(X0, Y0, contour4[pi], colors=cc4, linestyles=ls4, linewidths=lw4, levels=contourLevels4[pi], zorder=14)



		if isf[pi] is not None:

			extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]

			m.imshow(1-isf[pi], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)

		

		if u[pi] is not None:

			if C is not None:

				cax = ax.quiver(Xd0, Yd0, u[pi], v[pi], C[pi], cmap=ccmap, scale=scale[pi])

				if width_ratios is None:

					plt.colorbar(cax, ax=ax)

			else:

				cax = plt.quiver(Xd0, Yd0, u[pi], v[pi], scale=scale[pi], width=width, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, color=qcolor[pi])

			ax.quiverkey(cax, qlabelx, qlabely, qs[pi], str(qs[pi]) + ' ' + qunits, labelpos=labelpos[pi], coordinates='axes', labelcolor=qcolor[pi])

				

		doLabels(xlabel[pi], ylabel[pi], fontsize=fontsize)

		doTitle(titles[pi], fontsize=fstitle[pi])



		if text_data is not None:

			setText(ax, text_data[pi], set_invisible=False)

		

		if labelData is not None:

			for li in labelData[pi]:

				plt.scatter(li['x'][0], li['x'][1], s=1, color='r')

				plt.annotate(li['t'], li['x'])



		#ax.set_aspect('equal')



		if grid[pi]:

			plt.grid()

		

		ax.set_aspect('equal')

		if parallels[pi] is not None:

			# Latitudes

			m.drawparallels(parallels[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

		if meridians[pi] is not None:

			# Longitudes

			m.drawmeridians(meridians[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)



	#==



	plt.tight_layout()



	if save:

		plt.savefig(outpath + outname, bbox_inches="tight")



	if show:

		plt.show()



	plt.close()



#==



def quiver2by2(u, v, Xd, Yd, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=None, contourLevels=None, cbarShared=False, cbarSharedData=None, figsize=None, title=None, fontsize=14, cmap='jet', vmin=None, vmax=None, xlabels=None, ylabels=None, save=False, outpath='', outname='quiver2by2.png', show=True, dpi=200, text_data=None, width_ratios=[1,1], labelData=None, cbar=True, grid=True, xticks=None, xticksvis=True, yticks=None, yticksvis=True, scale=2, linewidths=0.4):



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

	

	contour = makeList(contour, M, N)

	contourLevels = makeList(contourLevels, M, N)



	fig = plt.figure(figsize=figsize, dpi=dpi)#, constrained_layout=True)



	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)



	for col in range(M):

		for row in range(N):



			ax = fig.add_subplot(gs[row,col])

			ax.patch.set_color('.6')



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

					



			if contour is not None:

				if contourLevels is not None:

					plt.contour(X[row][col], Y[row][col], contour[row][col], levels=contourLevels[row][col], colors='k', linestyles='solid', linewidths=linewidths)

				else:

					plt.contour(X[row][col], Y[row][col], contour[row][col], levels=contourLevels[row][col], colors='k', linestyles='solid', linewidths=linewidths)

			

			

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



	ax.patch.set_color('.6')



	if contourf is not None:

		if mesh:

			plt.pcolormesh(X[0], Y[0], contourf[0], vmin=vmin[0], vmax=vmax[0], cmap=cmap)		

		else:

			levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)

			plt.contourf(X[0], Y[0], contourf[0], cmap=cmap, levels=levels)



		if width_ratios is None:

			plt.colorbar()



	if contour[0] is not None:

		plt.contour(X[0], Y[0], contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.3)



	if C is not None:

		cax = ax.quiver(Xd[0], Yd[0], u[0], v[0], C[0], cmap=ccmap, scale=scale[0])

		if width_ratios is None:

			print('1')

			plt.colorbar(cax, ax=ax)

	else:

		cax = plt.quiver(Xd[0], Yd[0], u[0], v[0], scale=scale[0])

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

	ax.patch.set_color('.6')



	if contourf is not None:

		if mesh:

			plt.pcolormesh(X[1], Y[1], contourf[1], vmin=vmin[1], vmax=vmax[1], cmap=cmap)		

		else:

			levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)

			plt.contourf(X[1], Y[1], contourf[1], cmap=cmap, levels=levels)

		plt.colorbar()

		

	if contour[1] is not None:

		plt.contour(X[1], Y[1], contour[1], levels=contourLevels[1], colors='k', linestyles='solid', linewidths=0.3)



	if C is not None:

		cax = plt.quiver(Xd[1], Yd[1], u[1], v[1], C[1], cmap=ccmap, scale=scale[1], linewidths=2)

		plt.colorbar(cax, ax=ax)

	else:

		cax = ax.quiver(Xd[1], Yd[1], u[1], v[1], scale=scale[1])

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



def quiver1by2Basemap(u, v, Xd, Yd, lat_0, lon_0, C=None, ccmap='bwr', contourf=None, contourfNlevels=9, X=None, Y=None, mesh=False, contour=[None]*2, contourLevels=None, contour2=[None]*2, contourLevels2=None, figsize=(8,3), title=None, fontsize=14, cmap='bwr', vmin=None, vmax=None, xlabel=None, ylabel=None, save=False, outpath='', outname='quiver1by2.mp4', show=True, dpi=200, text_data=[None]*2, parallels=None, meridians=None, width_ratios=None, labelData=None, isf=[None]*2, extend=['neither', 'neither'], qs=[0.2]*2, scale=1., qunits='m/s'):

	

		

	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	fig = plt.figure(figsize=figsize, dpi=dpi)



	if width_ratios is not None:

		gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)

		ax = fig.add_subplot(gs[0,0])

	else:

		plt.subplot(121)

		ax = plt.gca()



	ax.patch.set_color('.6')



	X0, Y0 = m(X[0],Y[0])

	Xd0, Yd0 = m(Xd[0],Yd[0])



	if contourf is not None:

		if mesh:

			m.pcolormesh(X0, Y0, contourf[0], vmin=vmin[0], vmax=vmax[0], cmap=cmap, extend=extend[0])		

		else:

			levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels)

			m.contourf(X0, Y0, contourf[0], cmap=cmap, levels=levels, extend=extend[0])



		if width_ratios is None:

			plt.colorbar()



	if contour[0] is not None:

		if contourLevels is not None:

			m.contour(X0, Y0, contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.4)

		else:

			m.contour(X0, Y0, contour[0], levels=contourLevels[0], colors='k', linestyles='solid', linewidths=0.4)



	if contour2[0] is not None:

		m.contour(X0, Y0, contour2[0], colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2[0])

			

	if C is not None:

		cax = ax.quiver(Xd0, Yd0, u[0], v[0], C[0], cmap=ccmap, scale=scale)

		if width_ratios is None:

			plt.colorbar(cax, ax=ax)

	else:

		cax = m.quiver(Xd0, Yd0, u[0], v[0], scale=scale)

	ax.quiverkey(cax, 0.12, 0.03, qs[0], str(qs[0]) + ' ' + qunits, labelpos='E', coordinates='axes')

			

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

		

	if text_data[0] is not None:

		setText(ax, text_data[0], set_invisible=False)



	#==



	plt.subplot(122)

	ax = plt.gca()

	ax.patch.set_color('.6')



	X0, Y0 = m(X[1],Y[1])

	Xd0, Yd0 = m(Xd[1],Yd[1])



	if contourf is not None:

		if mesh:

			m.pcolormesh(X0, Y0, contourf[1], vmin=vmin[1], vmax=vmax[1], cmap=cmap, extend=extend[1])		

		else:

			levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels)

			m.contourf(X0, Y0, contourf[1], cmap=cmap, levels=levels, extend=extend[1])

		plt.colorbar()

		

	if contour[1] is not None:

		if contourLevels is not None:

			m.contour(X0, Y0, contour[1], levels=contourLevels[1], colors='k', linestyles='solid', linewidths=0.4)

		else:

			m.contour(X0, Y0, contour[1], levels=contourLevels[1], colors='k', linestyles='solid', linewidths=0.4)



	if contour2[1] is not None:

		m.contour(X0, Y0, contour2[1], colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2[1])

		

	if C is not None:

		cax = m.quiver(Xd0, Yd0, u[1], v[1], C[1], cmap=ccmap, scale=1., linewidths=2)

		plt.colorbar(cax, ax=ax)

	else:

		cax = ax.quiver(Xd0, Yd0, u[1], v[1], scale=scale)

	ax.quiverkey(cax, 0.22, 0.03, qs[0], str(qs[0]) + ' ' + qunits, labelpos='E', coordinates='axes')

			

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

	

	if text_data[1] is not None:

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



def plot1by2(data, X=[None]*2, Y=[None]*2, contour=[None]*2, contourLevels=[None]*2, contour2=[None]*2, contourLevels2=[None]*2, DASH_NEG_CONTOURS=False, figsize=(9,4), titles=[None]*2, fontsize=14, titlefontsize=None, mesh=False, cmaps=['jet']*2, vmin=[None]*2, vmax=[None]*2, stippling=[None]*2, stipData=[0.05, 4, 4, .1], text_data=[None]*2, xlabels=[None]*2, ylabels=[None]*2, grid=True, contourfNlevels=[9]*2, save=False, outpath='', outname='plot1by2.png', show=True, dpi=200, width_ratios=None, cbar=[True]*2, cbarLabel=[None]*2, xticks=[None]*2, yticks=[None]*2, xticksvis=[True]*2, yticksvis=[True]*2, vlines=[None]*2, patches=[None]*2, extend=['neither']*2):

	

	X = makeList(X, 2)

	Y = makeList(Y, 2)

	contour = makeList(contour, 2)

	contour2 = makeList(contour2, 2)

	#contourLevels = makeList(contourLevels, 2, FORCE_MAKE_LIST=True)

	extend = makeList(extend, 2)

	vmin = makeList(vmin, 2)

	vmax = makeList(vmax, 2)

	cmaps = makeList(cmaps, 2)

	contourfNlevels = makeList(contourfNlevels, 2)

		

	if titlefontsize is None:

		titefontsize = fontsize

	

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



		ax.patch.set_color('.6')



		if mesh:

			if X[pi] is not None and Y[pi] is not None:

				plt.pcolormesh(X[pi], Y[pi], data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])

			else: 

				plt.pcolormesh(data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])

		else:

			levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels[pi])

			if X[pi] is not None and Y[pi] is not None:

				plt.contourf(X[pi], Y[pi], data[pi], cmap=cmaps[pi], levels=levels, extend=extend[pi])

			else: 

				plt.contourf(data[pi], cmap=cmaps[pi], levels=levels, extend=extend[pi])



		if cbar[pi]:

			cbar_ = plt.colorbar()

			if cbarLabel[pi] is not None:

				cbar_.set_label(cbarLabel[pi], rotation=270, labelpad=20, fontsize=10)



		if contour[pi] is not None:

			if DASH_NEG_CONTOURS:

				linestyles = np.where(np.array(contourLevels[pi]) > 0, "-", "--")

			else:

				linestyles = 'solid'

			if X[pi] is not None and Y[pi] is not None:

				plt.contour(X[pi], Y[pi], contour[pi], levels=contourLevels[pi], colors='k', linestyles=linestyles, linewidths=0.8)

			else:

				plt.contour(X[pi], Y[pi], contour[pi], levels=contourLevels[pi], colors='k', linestyles=linestyles, linewidths=0.8)

		

		if contour2[pi] is not None:

			plt.contour(X[pi], Y[pi], contour2[pi], colors='k', linestyles='solid', linewidths=1.2, levels=contourLevels2[pi])



		if vlines[pi] is not None:

			for vl in vlines[pi]:

				try:

					plt.axvline(x=vl[0], ymin=vl[1], ymax=vl[2], linewidth=1.0, linestyle='--', color='k')

				except:

					try:

						plt.axvline(x=vl, color='k', linewidth=1.0, linestyle='--')

					except:

						print('Error: plotting.plo1tby2. Entries of vlines must be list with 3 elements or float.')

						quit()

					

		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)

		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])

		doTitle(titles[pi], fontsize=titlefontsize)

		

		if stippling[pi] is not None:

			doStippling(stippling[pi], X=X[pi], Y=Y[pi], dataLim=stipData[0], dx=stipData[1], dy=stipData[2], s=stipData[3])

		

		if patches[pi] is not None:

			ax.add_patch(patches[pi])

	

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

		

	plt.close()

		



#==



def plot1by3(data, X=[None]*3, Y=[None]*3, contour=[None]*3, contourLevels=[None]*3, DASH_NEG_CONTOURS=False, figsize=(11,3), titles=[None]*3, fstitle=None, fontsize=14, mesh=False, cmaps=['jet']*3, vmin=[None]*3, vmax=[None]*3, extend=['neither']*3, text_data=[None]*3, xlabels=[None]*3, ylabels=[None]*3, grid=True, gridls='solid', gridlw=0.8, gridlc='k', contourfNlevels=9, clabel=[False]*3, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, cbarTicks=[None]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3, bbox=False):

	

	vmin = makeList(vmin, 3)

	vmax = makeList(vmax, 3)

	cmaps = makeList(cmaps, 3)



	contourfNlevels = makeList(contourfNlevels, 3)

	

	X = makeList(X, 3)

	Y = makeList(Y, 3)	

	

	if fstitle is None:

		fstitle = [fontsize]*3

	else:

		fstitle = makeList(fstitle, 3)

	

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



		ax.patch.set_color('.6')



		if mesh:

			if X[pi] is not None and Y[pi] is not None:

				plt.pcolormesh(X[pi], Y[pi], data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])

			else: 

				plt.pcolormesh(data[pi], cmap=cmaps[pi], vmin=vmin[pi], vmax=vmax[pi])

		else:

			levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels[pi])

			if X[pi] is not None and Y[pi] is not None:

				plt.contourf(X[pi], Y[pi], data[pi], cmap=cmaps[pi], levels=levels, extend=extend[pi])

			else: 

				plt.contourf(data[pi], cmap=cmaps[pi], levels=levels, extend=extend[pi])



		if cbar[pi]:

			cbar_ = plt.colorbar(ticks=cbarTicks[pi])

			cbar_.formatter.set_useOffset(False)

			

		if contour[pi] is not None:

			if DASH_NEG_CONTOURS:

				linestyles = np.where(np.array(contourLevels[pi]) > 0, "-", "--")

			else:

				linestyles = 'solid'

			if X[pi] is not None and Y[pi] is not None:

				CS = plt.contour(X[pi], Y[pi], contour[pi], levels=contourLevels[pi], colors='k', linestyles=linestyles, linewidths=0.8)

			else:

				CS = plt.contour(X[pi], Y[pi], contour[pi], levels=contourLevels[pi], colors='k', linestyles=linestyles, linewidths=0.8)

			if clabel[pi]:

				plt.clabel(CS, fmt='%2.1f', colors='k', fontsize=7)



		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)

		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])

		doTitle(titles[pi], fontsize=fstitle[pi])

	

		if grid:

			plt.grid(linestyle=gridls, color=gridlc, linewidth=gridlw)



		if text_data[pi] is not None:

			setText(plt.gca(), text_data[pi], set_invisible=False)



	if bbox:

		box = ax1.get_position()

		box.x0 = box.x0 + 0.05

		ax1.set_position(box)





	#==



	plt.tight_layout()

	

	if save:

		plt.savefig(outpath + outname, bbox_inches="tight")

	if show:

		plt.show()

		

#==



def plotMbyN(data, X=None, Y=None, figsize=(8,4), titles=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, cbar=True, cbarShared=False, cbarSharedData=None, text_data=None, xlabels=None, ylabels=None, grid=True, contourfNlevels=	13, save=False, outpath='', outname='plotMbyN.png', show=True, dpi=200, width_ratios=None, xticks=None, yticks=None, xticksvis=True, yticksvis=True, hlines=None, contour=None, contourLevels=None):



	if not isinstance(data, list):

		plot1by1(data, X=X, Y=Y, mesh=mesh, vmin=vmin, vmax=vmax)

		return

		

	N = len(data)

	M = len(data[0])

	

	if X is None:

		nX = data[0][0].shape[1]

		X = np.linspace(1, nX, nX)

	if Y is None:

		nY = data[0][0].shape[0]

		Y = np.linspace(1, nY, nY)	

	

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

			

	contour = makeList(contour, M, N)

	contourLevels = makeList(contourLevels, M, N, FORCE_MAKE_LIST=True)

				

	if width_ratios is None:

		width_ratios = [1]*M

	

	#==



	fig = plt.figure(figsize=figsize, dpi=dpi)

	gs = gridspec.GridSpec(ncols=M, nrows=N, figure=fig, width_ratios=width_ratios)

	

	for col in range(M):

		for row in range(N):



			ax = fig.add_subplot(gs[row,col])

			ax.patch.set_color('.6')



			if mesh:

				im = plt.pcolormesh(X[row][col], Y[row][col], data[row][col], vmin=vmin[row][col], vmax=vmax[row][col], cmap=cmap[row][col])		

			else:

				levels = getContourfLevels(vmin[row][col], vmax[row][col], contourfNlevels)

				im = plt.contourf(X[row][col], Y[row][col], data[row][col], cmap=cmap[row][col], levels=levels, extend='both')

			if cbar[row][col]:

				plt.colorbar()

				

			if contour[row][col] is not None:

				plt.contour(X[row][col], Y[row][col], contour[row][col], linewidths=.8, colors='k', levels=contourLevels[row][col], linestyles='solid')

				

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

	

	plt.tight_layout()

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()

			



#==



def quiver2plot1(u, v, Xd, Yd, data, X, Y, C=[None]*2, ccmap='bwr', contourf=[None]*2, mesh=[True]*3, grid=[True]*3, figsize=(11,3), vmin=[None]*3, vmax=[None]*3, cmaps=['jet']*3, titles=[None]*3, fontsize=14, xlabels=[None]*3, ylabels=[None]*3, save=False, outpath='', outname='plot1by3.png', show=True, dpi=200, width_ratios=None, cbar=[True]*3, xticks=[None]*3, yticks=[None]*3, xticksvis=[True]*3, yticksvis=[True]*3, contour=[None]*3, contourLevels=[None]*3, contourfNlevels=9, qs=[1,1], scales=[1,1], text_data=[None]*3, qsfs=8):





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



		ax.patch.set_color('.6')



		if contourf[pi] is not None:

			if mesh[pi]:

				im = plt.pcolormesh(X[pi], Y[pi], contourf[pi], vmin=vmin[pi], vmax=vmax[pi], cmap=cmaps[pi])		

			else:

				levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)

				im = plt.contourf(X[pi], Y[pi], contourf[pi], cmap=cmaps[pi], levels=levels)



		if cbar[pi]:

			plt.colorbar()



		if C[pi] is not None:

			cax = ax.quiver(Xd[pi], Yd[pi], u[pi], v[pi], C[pi], cmap=ccmap, scale=scales[pi])	

			plt.colorbar(cax, ax=ax)

		else:

			cax = plt.quiver(Xd[pi], Yd[pi], u[pi], v[pi], scale=scales[pi])

		ax.quiverkey(cax, 0.1, 0.03, qs[pi], str(qs[pi]) + ' m/s', labelpos='N', coordinates='axes', fontproperties={'size':qsfs})



		

		if contour[pi] is not None:

			plt.contour(X[pi], Y[pi], contour[pi], colors='k', linestyles='solid', linewidths=0.6, levels=contourLevels[pi])

				

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

		

	ax.patch.set_color('.6')

				

	if mesh[pi]:

		plt.pcolormesh(X[pi], Y[pi], data, vmin=vmin[pi], vmax=vmax[pi], cmap=cmaps[pi], extend='both')		

	else:

		levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels)

		plt.contourf(X[pi], Y[pi], data, cmap=cmaps[pi], levels=levels, extend='both')

		

	if cbar[pi]:

		plt.colorbar()

		

	if contour[pi] is not None:

		plt.contour(X[pi], Y[pi], contour[pi], colors='k', linestyles='solid', linewidths=0.6, levels=contourLevels[pi])

			

	doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)

	doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi])

	doTitle(titles[pi], fontsize=fontsize)

			

	if grid[pi]:

		plt.grid()

		

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



def animate1by1(data, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, text_data2=None, hlines=[], contourfNlevels=9, contourAnim=None, contourAnimLevels=None, contour=None, contourLevels=[-1000], contour2=None, contourLevels2=[-1000]):



	# Make animation

	fig = plt.figure(figsize=figsize, dpi=dpi)

	ax = fig.add_subplot()

	plt.gca().patch.set_color('.25')

	plt.title(title, fontsize=fontsize)

		

	if mesh:

		if X is not None and Y is not None:

			cax = ax.pcolormesh(X, Y, data[0], cmap=cmap, vmin=vmin, vmax=vmax)

		else: 

			cax = ax.pcolormesh(data[0], cmap=cmap, vmin=vmin, vmax=vmax)

			

		if text_data is not None:

			setText(ax, text_data, i=0)	

		if text_data2 is not None:

			setText(ax, text_data2, i=0)	

			

		if contour is not None:

			plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)



		if contour2 is not None:

			plt.contour(X, Y, contour2, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels2)

			

		for hline in hlines:

			plt.axvline(hline, color='k', linewidth=1.0, linestyle='--')



		plt.grid()

		

		fig.colorbar(cax, ax=ax)

		plt.xlabel(xlabel, fontsize=fontsize)

		plt.ylabel(ylabel, fontsize=fontsize)



		plt.grid() 

	

		def animate(i):

			cax.set_array(data[i].flatten())

			plt.gca().set_title(title)

			if text_data is not None:

				setText(ax, text_data, i=i, set_invisible=True)

			if text_data2 is not None:

				setText(ax, text_data2, i=i, set_invisible=False)

			for hline in hlines:

				plt.axvline(hline, color='k', linewidth=1.0, linestyle='--')



		# End if mesh	



	else:

		nlevels = contourfNlevels

		if X is not None and Y is not None:

			cax = ax.contourf(X, Y, data[0], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))

		else: 

			cax = ax.contourf(data[0], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))

			

		if contourAnim is not None:

			ax.contour(X, Y, contourAnim[0], linstyles=solid, linewidths=0.8, levels=contourAnimLevels)

			

		if text_data is not None:

			setText(ax, text_data, i=0)

		if text_data2 is not None:

			setText(ax, text_data2, i=0)

			

		if contour is not None:

			plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)



		if contour2 is not None:

			print('contour2')

			plt.contour(X, Y, contour2, colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2)

			

		for hline in hlines:

			plt.axvline(hline, color='k', linewidth=1.0, linestyle='--')



		fig.colorbar(cax, ax=ax)

		plt.xlabel(xlabel, fontsize=fontsize)

		plt.ylabel(ylabel, fontsize=fontsize)

		plt.title(title, fontsize=fontsize)

		plt.grid() 



		def animate(i):

			ax.clear()

			#plt.grid()

			#cax.set_data(data[i].flatten())

			plt.gca().set_title(title)	

			plt.gca().patch.set_color('.25')		

			if text_data is not None:

				setText(ax, text_data, i=i, set_invisible=True)

			if text_data2 is not None:

				setText(ax, text_data2, i=i, set_invisible=False)

			ax.contourf(X, Y, data[i], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))	

			if contourAnim is not None:

				ax.contour(X, Y, contourAnim[i], linstyles=solid, linewidths=0.8, levels=contourAnimLevels)

			if contour is not None:

				ax.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)

			if contour2 is not None:

				ax.contour(X, Y, contour2, colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2)

			for hline in hlines:

				plt.axvline(hline, color='k', linewidth=1.0, linestyle='--')

	#==

	



	

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



def animate1by1varCbar(data, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, contour=None, hline=None, vmin=None, vmax=None):



	# Make animation

	fig = plt.figure(figsize=figsize, dpi=dpi)

	ax = fig.add_subplot(111)

	ax.patch.set_color('.25')

	

	div = make_axes_locatable(ax)

	cax = div.append_axes('right', '5%', '5%')



	if mesh:

		#if X is not None and Y is not None:

		#	cax = ax.pcolormesh(X, Y, data[0], cmap=cmap)

		#else: 

		#	cax = ax.pcolormesh(data[0], cmap=cmap)

		#if text_data is not None:

		#	setText(ax, text_data, i=0)	





		ax.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4)



		#fig.colorbar(cax, ax=ax)

		

		def animate(i):

			#ax.cla()

			#cax.set_array(data[i].flatten())

			vmax = 0.5 * np.max(np.abs(data[i]))

			pc = ax.pcolormesh(X, Y, data[i], cmap=cmap, vmin=-vmax, vmax=vmax)

			if text_data is not None:

				setText(ax, text_data, i=i, set_invisible=True)

			if hline is not None:

				plt.axvline(hline, color='k', linewidth=1.0, linestyle='--')

			cax.cla()

			fig.colorbar(pc, cax=cax)

			ax.text(X[0,1], Y[-1,0]+15, title, fontsize=fontsize)



	#==

	





	

	# Get number of timesteps/frames.

	if isinstance(data, list):

		Nt = len(data)

	else:

		Nt = data.shape[0]

		

	anim = animation.FuncAnimation(fig, animate, interval=50, frames=Nt)

	plt.draw()

	

	if save:

		anim.save(outpath+outname,metadata={'artist':'Guido'},writer='ffmpeg',fps=fps,bitrate=bitrate)

		

	if show:

		plt.show()



#==





def animateLine(data, X=None, figsize=(5,4), title='', labels=None, fontsize=14, vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animateLine.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, constLine=None, constLineLabel=None, constLineStyle=['solid', 'dashed', 'dotted'], transposeAx=False):



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



def animate1by1quiver(u, v, Xd, Yd, qlim=0.1, qunits='m/s', C=None, ccmap='coolwarm', contourf=None, contourfNlevels=21, contour=None, contourLevels=[-1000, -500], contour2=None, contourLevels2=None, isf=None, X=None, Y=None, cmap='viridis', vmin=None, vmax=None, mesh=False, grid=True, figsize=(5,4), title='', fontsize=14, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, text_data2=None, scale=None):



	# Make animation

	fig = plt.figure(figsize=figsize, dpi=dpi)

	ax = fig.add_subplot()

	plt.gca().patch.set_color('.25')



	if contourf is not None:

		if not mesh:

			nlevels = contourfNlevels

			cax = ax.contourf(X, Y, contourf[0], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))

		else:

			cax = plt.pcolormesh(X, Y, contourf[0], vmin=vmin, vmax=vmax, cmap=cmap)		



		fig.colorbar(cax, ax=ax)



	if contour is not None:

		plt.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)

	if contour2 is not None:

		plt.contour(X, Y, contour2, colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2)

		

	plt.xlabel(xlabel, fontsize=fontsize)

	plt.ylabel(ylabel, fontsize=fontsize)



	if isf is not None:

		extent = [X[0], X[-1], Y[0], Y[-1]]

		ax.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)



	plt.grid() 



	if C is not None:

		Q = ax.quiver(Xd, Yd, u[0], v[0], C[0], cmap=ccmap, scale=scale)

		plt.colorbar(Q, ax=ax)

	else:

		Q = ax.quiver(Xd, Yd, u[0], v[0], scale=scale)

	ax.quiverkey(Q, 0.1, 0.1, qlim, str(qlim) + ' ' + qunits, labelpos='E', coordinates='axes', labelcolor='w')



	if text_data is not None:

		setText(ax, text_data, i=0)

	if text_data2 is not None:

		setText(ax, text_data2, i=0)

		

	if grid:

		plt.grid(linewidth=0.5)

	

	#==



	def animate(i):



		if text_data is not None:

				setText(ax, text_data, i=i, set_invisible=True)

		if text_data2 is not None:

			setText(ax, text_data2, i=i, set_invisible=False)

			

		if contourf is not None:

			if mesh:

				cax.set_array(contourf[i].flatten())

			else:

				ax.contourf(X, Y, contourf[i], cmap=cmap, levels=np.linspace(vmin, vmax, nlevels))			

		

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



def animate1by1quivers(udata, vdata, X, Y, d=8, qlim=0.1, C=None, ccmap='coolwarm', contour=None, cmap='viridis', vmin=None, vmax=None, contourf=True, grid=True, figsize=(5,4), title='', fontsize=14, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=300, fps=8, bitrate=-1, text_data=None, colors=['b', 'r'], scale=2, labels=['','']):

	'''The same as above, but can plot more than one vector field.'''



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

		elif len(contour.shape) == 3:

			ctr = 1

			if contourf:

				C = ax.contourf(X, Y, contour[0], cmap=cmap)

			else:

				C = ax.pcolormesh(X, Y, contour[0], vmin=vmin, vmax=vmax, cmap=cmap)		



	Qs = []

	for di, u in enumerate(udata):

		v = vdata[di]	

		Qs.append(ax.quiver(X[::d], Y[::d], u[0,::d,::d], v[0,::d,::d], color=colors[di], width=.002, scale=scale))

		ax.quiverkey(Qs[di], 0.05, 0.1+0.05*di, qlim, str(qlim) + ' m/s; ' + labels[di], labelpos='E', coordinates='axes')



	if text_data is not None:

		setText(ax, text_data, i=0)



	if grid:

		plt.grid(linewidth=0.5)



	def animate(i):



		for di, u in enumerate(udata):

			v = vdata[di]

			Qs[di].set_UVC(u[i,::d,::d], v[i,::d,::d])

		if text_data is not None:

			setText(ax, text_data, i=i, set_invisible=True)



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



def plot1by1Basemap(data, X, Y, lat_0, lon_0, contour=None, contourLevels=None, figsize=(5,4), title=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabel=None, ylabel=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, labelData=None, AntarcticInsetData=None):

	

	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	X,Y = m(X,Y)



	fig = plt.figure(figsize=figsize, dpi=dpi)



	ax = fig.add_subplot()

	plt.gca().patch.set_color('.6')



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

		m.contour(X, Y, contour, colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels)

			

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

		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1], fontsize=fontsize)

	if meridians is not None:

		# Longitudes

		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1], fontsize=fontsize)



	for li in labelData:

		plt.scatter(li['x'][0], li['x'][1], s=1, color='y')

		plt.annotate(li['t'], li['x'])



	if AntarcticInsetData is not None:

		doAntarcticInset(fig, AntarcticInsetData)

		

	#==

	

	plt.tight_layout()

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()



#==

		

def plot1by2Basemap(data, X, Y, lat_0, lon_0, contour=[None]*2, contourLevels=[None]*2, contour2=[None]*2, contourLevels2=[None]*2, isf=[None]*2, figsize=(8,3), titles=None, fstitle=None, fontsize=14, mesh=False, cmap=['jet']*2, vmin=[None]*2, vmax=[None]*2, text_data=[None,None], xlabels=None, ylabels=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1by2Basemap.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, labelData=None, cbar=[True]*2, cbarTicks=[None, None], cbarShrink=1., extend=['neither']*2, xticks=[None]*2, yticks=[None]*2, xticksvis=[True]*2, yticksvis=[True]*2, stippling=[None]*2, stipData=[0.05, 4, 4, .1], width_ratios=[1,1], landscape=True):

	

	X = makeList(X, 2)

	Y = makeList(Y, 2)

	vmin = makeList(vmin, 2)

	vmax = makeList(vmax, 2)

	

	titles = makeList(titles, 2)

	grid = makeList(grid, 2)

	xlabels = makeList(xlabels, 2)

	ylabels = makeList(ylabels, 2)	



	cmap = makeList(cmap, 2)

	contourfNlevels = makeList(contourfNlevels, 2)

	

	if fstitle is None:

		fstitle = [fontsize]*2

		

	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	fig = plt.figure(figsize=figsize, dpi=dpi)

	if landscape:

		gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)

	else:

		gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=width_ratios)



	for pi in range(2):



		if landscape:

			ax = fig.add_subplot(gs[0,pi])

		else:

			ax = fig.add_subplot(gs[pi,0])

					

		if pi == 1:

			ax1 = ax



		ax.patch.set_color('.6')

		

		X0, Y0 = m(X[pi],Y[pi])



		if mesh:

			cax = m.pcolormesh(X0, Y0, data[pi], cmap=cmap[pi], vmin=vmin[pi], vmax=vmax[pi])

		else:

			levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels[pi])

			cax = m.contourf(X0, Y0, data[pi], cmap=cmap[pi], levels=levels, extend=extend[pi])

	

		if cbar[pi]:

			plt.colorbar(ticks=cbarTicks[pi], shrink=cbarShrink)

		

		if contour[pi] is not None:

			m.contour(X0, Y0, contour[pi], colors='k', linestyles='solid', linewidths=0.4, levels=contourLevels[pi])

		

		if contour2[pi] is not None:

			m.contour(X0, Y0, contour2[pi], colors='k', linestyles='solid', linewidths=0.8, levels=contourLevels2[pi])

			

		if isf[pi] is not None:

			extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]

			m.imshow(1-isf[pi], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)

		

		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)

		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi], fontsize=fontsize)

		doTitle(titles[pi], fontsize=fstitle[pi])

		

		if grid:

			plt.grid()



		if text_data[pi] is not None:

			setText(plt.gca(), text_data[pi], set_invisible=False)

			

		if stippling[pi] is not None:

			doStippling(stippling[pi], X=X0, Y=Y0, dataLim=stipData[0], dx=stipData[1], dy=stipData[2], s=stipData[3])

		

		ax.set_aspect(1)

		if parallels is not None:

			# Latitudes

			m.drawparallels(parallels[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

		if meridians is not None:

			# Longitudes

			m.drawmeridians(meridians[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

			

	#==

	

	plt.tight_layout()

	

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()



#==



def plot1byNbasemap(data, X, Y, lat_0, lon_0, contour=None, contourLevels=None, ls='solid', lw=0.4, cc='k', contour2=None, contourLevels2=None, ls2='solid', lw2=0.4, cc2='grey', contour3=None, contourLevels3=None, ls3='solid', lw3=0.4, cc3='w', contour4=None, contourLevels4=None, ls4='solid', lw4=0.4, cc4='w', contour5=None, contourLevels5=None, ls5='solid', lw5=0.4, cc5='k', isf=None, figsize=(8,3), titles=None, fstitle=None, fontsize=14, mesh=False, cmap='jet', vmin=None, vmax=None, text_data=None, xlabels=None, ylabels=None, grid=True, contourfNlevels=9, save=False, outpath='', outname='plot1byNbasemap.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, labelData=None, cbar=True, cbarTicks=None, cbarShrink=1., extend='neither', xticks=None, yticks=None, xticksvis=True, yticksvis=True, stippling=None, stipData=[0.05, 4, 4, .1], width_ratios=[1,1], landscape=True):

	

	N = len(data)

	

	X = makeList(X, N)

	Y = makeList(Y, N)

	vmin = makeList(vmin, N)

	vmax = makeList(vmax, N)

	

	titles = makeList(titles, N)

	grid = makeList(grid, N)

	xlabels = makeList(xlabels, N)

	ylabels = makeList(ylabels, N)	

	text_data = makeList(text_data, N)

	

	contour = makeList(contour, N)

	contour2 = makeList(contour2, N)

	contour3 = makeList(contour3, N)

	isf = makeList(isf, N)

	

	cbar = makeList(cbar, N)

	cbarTicks = makeList(cbarTicks, N)

	extend = makeList(extend, N)

	

	xticks = makeList(xticks, N)

	xticksvis = makeList(xticksvis, N)

	yticks = makeList(yticks, N)

	yticksvis = makeList(yticksvis, N)

	

	stippling = makeList(stippling, N)



	cmap = makeList(cmap, N)

	contourfNlevels = makeList(contourfNlevels, N)

	

	if fstitle is None:

		fstitle = [fontsize]*N

		

	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	fig = plt.figure(figsize=figsize, dpi=dpi)

	if landscape:

		gs = gridspec.GridSpec(ncols=N, nrows=1, figure=fig, width_ratios=width_ratios)

	else:

		gs = gridspec.GridSpec(ncols=1, nrows=N, figure=fig, height_ratios=width_ratios)



	for pi in range(N):



		if landscape:

			ax = fig.add_subplot(gs[0,pi])

		else:

			ax = fig.add_subplot(gs[pi,0])

					

		if pi == 1:

			ax1 = ax



		ax.patch.set_color('.6')

		

		X0, Y0 = m(X[pi],Y[pi])



		if mesh:

			cax = m.pcolormesh(X0, Y0, data[pi], cmap=cmap[pi], vmin=vmin[pi], vmax=vmax[pi])

		else:

			levels = getContourfLevels(vmin[pi], vmax[pi], contourfNlevels[pi])

			cax = m.contourf(X0, Y0, data[pi], cmap=cmap[pi], levels=levels, extend=extend[pi])

	

		if cbar[pi]:

			cbar_ = plt.colorbar(ticks=cbarTicks[pi], shrink=cbarShrink)

			#cbar_.formatter.set_powerlimits((0, 0))

			#cbar_.formatter.set_useMathText(True)

			

		if contour[pi] is not None:

			plt.contour(X0, Y0, contour[pi], colors=cc, linestyles=ls, linewidths=lw, levels=contourLevels[pi])

		if contour2[pi] is not None:

			plt.contour(X0, Y0, contour2[pi], colors=cc2, linestyles=ls2, linewidths=lw2, levels=contourLevels2[pi],zorder=14)

		if contour3[pi] is not None:

			plt.contour(X0, Y0, contour3[pi], colors=cc3, linestyles=ls3, linewidths=lw3, levels=contourLevels3[pi], zorder=14)

		if contour4[pi] is not None:

			plt.contour(X0, Y0, contour4[pi], colors=cc4, linestyles=ls4, linewidths=lw4, levels=contourLevels4[pi], zorder=14)

		if contour5[pi] is not None:

			plt.contour(X0, Y0, contour5[pi], colors=cc5, linestyles=ls5, linewidths=lw5, levels=contourLevels5[pi], zorder=14)



		if isf[pi] is not None:

			extent = [X0[0,0], X0[0,-1], -Y0[0,0], -Y0[-1,0]]

			m.imshow(1-isf[pi], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=13)

		

		doLabels(xlabels[pi], ylabels[pi], fontsize=fontsize)

		doTicks(xticks[pi], xticksvis[pi], yticks[pi], yticksvis[pi], fontsize=fontsize)

		doTitle(titles[pi], fontsize=fstitle[pi])

		

		if grid:

			plt.grid()



		if text_data[pi] is not None:

			setText(plt.gca(), text_data[pi], set_invisible=False)

			

		if stippling[pi] is not None:

			doStippling(stippling[pi], X=X0, Y=Y0, dataLim=stipData[0], dx=stipData[1], dy=stipData[2], s=stipData[3])

		

		ax.set_aspect(1)

		if parallels is not None:

			# Latitudes

			m.drawparallels(parallels[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

		if meridians is not None:

			# Longitudes

			m.drawmeridians(meridians[pi],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], fontsize=fontsize)

			

	#==

	

	plt.tight_layout()

	

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()



#==



# This is the same as above, but only the first panel uses Basemap 

def plot1by2Basemap1(data, X, Y, lat_0, lon_0, contour=[None,None], contourLevels=None, lws=[1,1], figsize=(7.5,3), titles=None, fontsize=14, tfontsize=9., mesh=False, cmaps=['jet', 'jet'], vmin=None, vmax=None, text_data=[None,None], xlabels=[None, None], ylabels=[None, None], grid=True, contourfNlevels=[9,9], save=False, outpath='', outname='plot1by1.png', show=True, dpi=200, yline=None, parallels=None, meridians=None, isf=[None, None], labelData=None, clabel=[False,False]):

	

	m = Basemap(llcrnrlon=X[0][0,0],llcrnrlat=Y[0][0,0],urcrnrlon=X[0][-1,-1],urcrnrlat=Y[0][-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	fig = plt.figure(figsize=figsize, dpi=dpi)



	plt.subplot(121)

	ax = plt.gca()

	ax.patch.set_color('.6')



	X0, Y0 = m(X[0],Y[0])



	if mesh:

		m.pcolormesh(X0, Y0, data[0], cmap=cmaps[0], vmin=vmin[0], vmax=vmax[0])

	else:

		levels = getContourfLevels(vmin[0], vmax[0], contourfNlevels[0])

		a = m.contourf(X0, Y0, data[0], cmap=cmaps[0], levels=levels, extend='both')

		#a.cmap.set_under('w')	

		#a.set_clim(vmin[0], vmax[0])

		#a.cmap.set_over('w')	

		#plt.colorbar(a)



	plt.colorbar()



	if contour[0] is not None:

		CS = plt.contour(X0, Y0, contour[0], levels=contourLevels[0], colors='k', linewidths=lws[0], linestyles='solid')

		if clabel[0]:

			plt.clabel(CS, fmt='%2.1f', colors='k', fontsize=7)

				

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

	ax.patch.set_color('.6')



	if mesh:

		plt.pcolormesh(X[1], Y[1], data[1], cmaps=cmaps[1], vmin=vmin[1], vmax=vmax[1])

	else:

		levels = getContourfLevels(vmin[1], vmax[1], contourfNlevels[1])

		plt.contourf(X[1], Y[1], data[1], cmap=cmaps[1], levels=levels)

	plt.colorbar()

		

	if contour[1] is not None:

		CS = plt.contour(X[1], Y[1], contour[1], levels=contourLevels[1], colors='k', linewidths=lws[1], linestyles='solid')

		if clabel[1]:

			plt.clabel(CS, fmt='%2.1f', colors='k', fontsize=7)



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



def quiver1by1Basemap(u, v, X, Y, d, lat_0, lon_0, qx=0.3, qy=0.03, qs=0.05, qunits='m/s', scale=1, vmin=None, vmax=None, width=2.8e6, height=1.7e6, contourf=None, contourfNlevels=13, contour=None, contourLevels=None, contourColours=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1Basemap.png', show=True, dpi=200, text_data=None, parallels=None, meridians=None, isf=None, land=None, outline=None, cbarTicks=None, extend='neither', AntarcticInsetData=None, maskColour='.6', tight_layout=True):



	#m = Basemap(width=width,height=height, projection='gnom',lat_0=lat_0,lon_0=lon_0)

	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	Xd = X[::d,::d]; Yd = Y[::d,::d]

	X,Y = m(X,Y)

	Xd,Yd = m(Xd,Yd)



	fig = plt.figure(figsize=figsize, dpi=dpi)

	

	ax = fig.add_subplot()

	

	plt.gca().patch.set_color(maskColour)



	if contourf is not None:

		if mesh:

			m.pcolormesh(X, Y, contourf, cmap=cmap, vmin=vmin, vmax=vmax)

			m.colorbar(extend=extend, ticks=cbarTicks)

		else:

			levels = getContourfLevels(vmin, vmax, contourfNlevels)

			m.contourf(X, Y, contourf, levels=levels, cmap=cmap)

			m.colorbar(extend=extend, ticks=cbarTicks)

	

	if contour is not None:	

		for ci in range(len(contour)):

			cs = m.contour(X, Y, contour[ci], levels=contourLevels[ci], colors=contourColours[ci], linestyles='solid', linewidths=0.8)

			#plt.clabel(cs, fontsize=7)



	if isinstance(u, (list)):

		q = m.quiver(Xd,Yd,u[0][::d, ::d],v[0][::d, ::d], color='r', scale=scale)

		#plt.quiverkey(q,0.1, -0.1, 0.1, r'0.1 N m$^{-2}$',labelpos='W')

		q = m.quiver(Xd,Yd,u[1][::d, ::d],v[1][::d, ::d], color='k', scale=scale)

		plt.quiverkey(q, qx, qy, qs, str(qs) + ' ' + qunits ,labelpos='W', fontproperties={'size':11})

	else:

		u = u[::d, ::d]; v = v[::d, ::d]

		q = m.quiver(Xd, Yd, u, v, scale=scale, width=0.002, headwidth=5., headlength=5., headaxislength=3, zorder=14)

		plt.quiverkey(q, qx, qy, qs, str(qs) + ' ' + qunits ,labelpos='W', fontproperties={'size':11})

	#m.fillcontinents(color='#cc9955', zorder = 0)

	

	# width 0.005 

	# headwidth 3

	# headlength 5

	# headaxislength 4.5



	if xlabel is not None:

		plt.xlabel(xlabel, labelpad=20)

	if ylabel is not None:

		plt.ylabel(ylabel, labelpad=35)



	if title is not None:

		plt.title(title, fontsize=fontsize) 



	if text_data is not None:

		setText(ax, text_data)



	ax.set_aspect(1)

	if parallels is not None:

		# Latitudes

		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4)

	if meridians is not None:

		# Longitudes

		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4)



	if isf is not None:

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);

		#m.contourf(X, Y, 1-isf, cmap=plt.cm.gray, extent=extent, zorder=3, vmin=0, vmax=1);

		

	if land is not None:

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-land, cmap='gray', interpolation='nearest', extent=extent, zorder=2);

		

	if outline is not None:

		lons = outline[0]

		lats = outline[1]

		lons, lats = m(lons, lats)



		m.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')



	if AntarcticInsetData is not None:

		doAntarcticInset(fig, AntarcticInsetData)

	

	#==



	if tight_layout:

		plt.tight_layout()

	

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()	

		

#==



def quiver1by1Basemap_timeseries(u, v, X, Y, d, lat_0, lon_0, qx=0.3, qy=0.03, qs=0.05, qunits='m/s', scale=1, vmin=None, vmax=None, width=2.8e6, height=1.7e6, plotvlines=[], contourf=None, contourfNlevels=13, contour=None, lw=1.2, contourLevels=None, contourColours=None, figsize=(5,4), title='', fstitle=14, fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1Basemap.png', show=True, dpi=200, text_data=None, parallels=None, meridians=None, isf=None, land=None, outline=None, cbarTicks=None, extend='neither', AntarcticInsetData=None, arrows=[], maskColour='.6', height_ratios=None, ts_time=None, ts_data=None, ts_labels=None, ts_lines=None, ts_title=None, ts_xlabel=None, ts_ylabel=None, ts_xlim=None, ts_ylim=None, ts_xticks=None, ts_yticks=None, ts_legendLoc=None, ts_legendFontsize=8, ts_legendNcol=3):



	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	Xd = X[::d,::d]; Yd = Y[::d,::d]

	X,Y = m(X,Y)

	Xd,Yd = m(Xd,Yd)



	fig = plt.figure(figsize=figsize, dpi=dpi)

	gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=height_ratios)

	ax = fig.add_subplot(gs[0])

	plt.gca().patch.set_color(maskColour)



	if contourf is not None:

		if mesh:

			m.pcolormesh(X, Y, contourf, cmap=cmap, vmin=vmin, vmax=vmax)

			cbar = m.colorbar(extend=extend, ticks=cbarTicks)

			cbar.ax.tick_params(labelsize=fontsize+1)

		else:

			levels = getContourfLevels(vmin, vmax, contourfNlevels)

			m.contourf(X, Y, contourf, levels=levels, cmap=cmap)

			cbar = m.colorbar(extend=extend, ticks=cbarTicks)

			cbar.ax.tick_params(labelsize=fontsize+1)

	

	if contour is not None:	

		for ci in range(len(contour)):

			cs = m.contour(X, Y, contour[ci], levels=contourLevels[ci], colors=contourColours[ci], linestyles='solid', linewidths=lw)

			#plt.clabel(cs, fontsize=7)



	if isinstance(u, (list)):

		q = m.quiver(Xd,Yd,u[0][::d, ::d],v[0][::d, ::d], color='r', scale=scale)

		#plt.quiverkey(q,0.1, -0.1, 0.1, r'0.1 N m$^{-2}$',labelpos='W')

		q = m.quiver(Xd,Yd,u[1][::d, ::d],v[1][::d, ::d], color='k', scale=scale)

		plt.quiverkey(q, qx, qy, qs, str(qs) + ' ' + qunits ,labelpos='W', fontproperties={'size':20})

	else:

		u = u[::d, ::d]; v = v[::d, ::d]

		q = m.quiver(Xd, Yd, u, v, scale=scale, width=0.002, headwidth=5., headlength=5., headaxislength=3, zorder=14)

		plt.quiverkey(q, qx, qy, qs, str(qs) + ' ' + qunits ,labelpos='W', fontproperties={'size':20})

	#m.fillcontinents(color='#cc9955', zorder = 0)

	

	if xlabel is not None:

		plt.xlabel(xlabel, labelpad=20)

	if ylabel is not None:

		plt.ylabel(ylabel, labelpad=35)



	if title is not None:

		plt.title(title, fontsize=fstitle) 



	if text_data is not None:

		setText(ax, text_data)



	ax.set_aspect(1)

	if parallels is not None:

		# Latitudes

		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4, fontsize=fontsize+2)

	if meridians is not None:

		# Longitudes

		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4, fontsize=fontsize+2)



	if isf is not None:

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);

		#m.contourf(X, Y, 1-isf, cmap=plt.cm.gray, extent=extent, zorder=3, vmin=0, vmax=1);

		

	if land is not None:

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-land, cmap='gray', interpolation='nearest', extent=extent, zorder=2);

		

	if outline is not None:

		lons = outline[0]

		lats = outline[1]

		lons, lats = m(lons, lats)



		m.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')





	for vline in plotvlines:

		nvline = 40

		m.plot(vline[0]*np.ones(nvline), np.linspace(vline[1], vline[2], nvline), color=vline[3], linestyle='dashed', latlon=True, lw=vline[4])



	for arrow in arrows:

		plt.arrow(arrow[0], arrow[1], arrow[2], arrow[3], color=arrow[4], zorder=15, linewidth=3, head_width=1.5e4, head_length=1.5e4)	

			

	if AntarcticInsetData is not None:

		doAntarcticInset(fig, AntarcticInsetData)

		#doAntarcticInset(ax, AntarcticInsetData)



	#==

	

	# Now timeseries.

	ax = fig.add_subplot(gs[1])



	for ti, timeseries in enumerate(ts_data):

		plt.plot(ts_time[ti], timeseries, color=ts_lines[ti][0], linestyle=ts_lines[ti][1], label=ts_labels[ti]) 

		

	plt.xlim(ts_xlim)

	plt.ylim(ts_ylim)

	plt.xticks(ts_xticks)

	plt.yticks(ts_yticks)

	plt.title(ts_title, fontsize=fstitle)

	#ax.tick_params(axis='both', which='minor', labelsize=fontsize)

	plt.xticks(fontsize=fontsize)

	plt.yticks(fontsize=fontsize)

	plt.xlabel(ts_xlabel, fontsize=fontsize)

	plt.ylabel(ts_ylabel, fontsize=fontsize)

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)

	plt.legend(prop={'size': ts_legendFontsize}, ncol=ts_legendNcol, loc=ts_legendLoc)

	

	#==	



	plt.tight_layout()



	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()	

		

#==



def plot1by1Basemap_timeseries(data, X, Y, lat_0, lon_0, qx=0.3, qy=0.03, qs=0.05, qunits='m/s', scale=1, vmin=None, vmax=None, width=2.8e6, height=1.7e6, stippling=None, stipData=[0.05, 4, 4, .1], plotvlines=[], contourf=None, contourfNlevels=13, contour=None, lw=1.2, contourLevels=None, contourColours=None, contourLS=None, contourLW=None, figsize=(5,4), title='', fstitle=14, fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1Basemap.png', show=True, dpi=200, text_data=None, parallels=None, meridians=None, isf=None, land=None, outline=None, cbarTicks=None, extend='neither', AntarcticInsetData=None, maskColour='.6', height_ratios=None, ts_time=None, ts_data=None, ts_labels=None, ts_lines=None, ts_title=None, ts_xlabel=None, ts_ylabel=None, ts_xlim=None, ts_ylim=None, ts_fslegend=9):



	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	X,Y = m(X,Y)



	fig = plt.figure(figsize=figsize, dpi=dpi)

	gs = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, height_ratios=height_ratios)

	ax = fig.add_subplot(gs[0])

	plt.gca().patch.set_color(maskColour)



	if mesh:

		m.pcolormesh(X, Y, data, cmap=cmap, vmin=vmin, vmax=vmax)

		cbar = m.colorbar(extend=extend, ticks=cbarTicks)

		cbar.ax.tick_params(labelsize=fontsize+1)

	else:

		levels = getContourfLevels(vmin, vmax, contourfNlevels)

		m.contourf(X, Y, data, levels=levels, cmap=cmap)

		cbar = m.colorbar(extend=extend, ticks=cbarTicks)

		cbar.ax.tick_params(labelsize=fontsize+1)

	

	if contour is not None:	

		for ci in range(len(contour)):

			cs = m.contour(X, Y, contour[ci], levels=contourLevels[ci], colors=contourColours[ci], linestyles=contourLS[ci], linewidths=contourLW[ci])

			#plt.clabel(cs, fontsize=7)

	

	if stippling is not None:

		doStippling(stippling, X=X, Y=Y, dataLim=stipData[0], dx=stipData[1], dy=stipData[2], s=stipData[3])

		

	if xlabel is not None:

		plt.xlabel(xlabel, labelpad=20)

	if ylabel is not None:

		plt.ylabel(ylabel, labelpad=35)



	if title is not None:

		plt.title(title, fontsize=fstitle) 



	if text_data is not None:

		setText(ax, text_data)



	ax.set_aspect(1)

	if parallels is not None:

		# Latitudes

		m.drawparallels(parallels,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4, fontsize=fontsize+2)

	if meridians is not None:

		# Longitudes

		m.drawmeridians(meridians,labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[4,4], zorder=4, fontsize=fontsize+2)



	if isf is not None:

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-isf, cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);

		#m.contourf(X, Y, 1-isf, cmap=plt.cm.gray, extent=extent, zorder=3, vmin=0, vmax=1);

		

	if land is not None:

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-land, cmap='gray', interpolation='nearest', extent=extent, zorder=2);

		

	if outline is not None:

		lons = outline[0]

		lats = outline[1]

		lons, lats = m(lons, lats)



		m.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')





	for vline in plotvlines:

		nvline = 40

		m.plot(vline[0]*np.ones(nvline), np.linspace(vline[1], vline[2], nvline), color=vline[3], linestyle='dashed', latlon=True, lw=vline[4])



	if AntarcticInsetData is not None:

		doAntarcticInset(fig, AntarcticInsetData)

		#doAntarcticInset(ax, AntarcticInsetData)

		

	#==

	

	# Now timeseries.

	ax = fig.add_subplot(gs[1])



	for ti, timeseries in enumerate(ts_data):

		plt.plot(ts_time[ti], timeseries, color=ts_lines[ti][0], linestyle=ts_lines[ti][1], label=ts_labels[ti]) 

		

	plt.xlim(ts_xlim)

	plt.ylim(ts_ylim)

	plt.title(ts_title, fontsize=fstitle)

	#ax.tick_params(axis='both', which='minor', labelsize=fontsize)

	plt.xticks(fontsize=fontsize)

	plt.yticks(fontsize=fontsize)

	plt.xlabel(ts_xlabel, fontsize=fontsize)

	plt.ylabel(ts_ylabel, fontsize=fontsize)

	plt.grid(color='k', linestyle='dashed', linewidth=0.5)

	plt.legend(prop={'size': ts_fslegend}, ncol=1)

	

	#==	



	plt.tight_layout()



	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()	

		

#==



def quiverPlotSubr1by2Basemap(u, v, X, Y, d, lat_0, lon_0, Xsubr, Ysubr, bathySubr, latSubr_0, lonSubr_0, vmin=None, vmax=None, width=2.8e6, height=1.7e6, contourf=None, figsize=(9,4), title='', fontsize=14, mesh=True, cmap='jet', xlabel='', ylabel='', save=False, outpath='', outname='quiver1by1Basemap.png', show=True, dpi=200, text_data=None, parallels=None, meridians=None, isf=[None,None], outline=None, extend='neither', labelData=[None,None], fs=14, width_ratios=[1.6,1.], cbarTicks=None):



	fig = plt.figure(figsize=figsize, dpi=dpi)

	gs = gridspec.GridSpec(ncols=2, nrows=1, figure=fig, width_ratios=width_ratios)

		

	ax = fig.add_subplot(gs[0])

	ax.patch.set_color('.6')

	

	#m = Basemap(width=width,height=height, projection='gnom',lat_0=lat_0,lon_0=lon_0)

	m = Basemap(llcrnrlon=X[0,0],llcrnrlat=Y[0,0],urcrnrlon=X[-1,-1],urcrnrlat=Y[-1,-1], projection='merc',lat_0=lat_0+3,lon_0=lon_0-10)



	Xd = X[::d,::d]; Yd = Y[::d,::d]

	X1,Y1 = m(X,Y)

	Xd,Yd = m(Xd,Yd)



	if contourf is not None:

		if mesh:

			m.pcolormesh(X1, Y1, contourf, cmap=cmap, vmin=vmin[0], vmax=vmax[0])

			m.colorbar(extend=extend[0], ticks=cbarTicks)

		else:

			m.contour(X1, Y1, contourf, levels=[-3000, -2000, -1000])

			m.colorbar(extend=extend[0], ticks=cbarTicks)



	if isinstance(u, (list)):

		q = m.quiver(Xd,Yd,u[0][::d, ::d],v[0][::d, ::d], color='r', scale=2)

		#plt.quiverkey(q,0.1, -0.1, 0.1, r'0.1 N m$^{-2}$',labelpos='W')

		q = m.quiver(Xd,Yd,u[1][::d, ::d],v[1][::d, ::d], color='k', scale=2)

		plt.quiverkey(q,0.94, 0.03, 0.1, r'0.1 N m$^{-2}$',labelpos='W')



	else:

		u = u[::d, ::d]; v = v[::d, ::d]

		q = m.quiver(Xd,Yd,u,v)

		plt.quiverkey(q,0.94, 0.03, 0.1, r'0.1 N m$^{-2}$',labelpos='W')

	#m.fillcontinents(color='#cc9955', zorder = 0)



	if xlabel is not None:

		plt.xlabel(xlabel, labelpad=20)

	if ylabel is not None:

		plt.ylabel(ylabel, labelpad=35)



	if title is not None:

		plt.title(title[0], fontsize=fs) 



	if text_data is not None:

		setText(ax, text_data, i=0)



	ax.set_aspect(1)

	if parallels is not None:

		# Latitudes

		m.drawparallels(parallels[0],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	if meridians is not None:

		# Longitudes

		m.drawmeridians(meridians[0],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])



	if labelData[0] is not None:

		for li in labelData[0]:

			plt.scatter(li['x'][0], li['x'][1], s=1, color='r')

			plt.annotate(li['t'], li['x'])



	if isf[0] is not None:

		#x, y = m(isf[0], isf[1])  # transform coordinates

		#m.scatter(x, y, color='k', s=1)

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-isf[0], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);

		#m.contour(X, Y, isf, levels=[1], colors='k')



	if outline is not None:

		lons = outline[0]

		lats = outline[1]

		lons, lats = m(lons, lats)



		m.plot([lons[0], lons[0]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[1], lons[1]], [lats[0], lats[1]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[0], lats[0]], color='k', linestyle='--')

		m.plot([lons[0], lons[1]], [lats[1], lats[1]], color='k', linestyle='--')



	#==

	

	# Panel 2

	

	ax = fig.add_subplot(gs[1])

	ax.patch.set_color('.6')

	

	m = Basemap(llcrnrlon=Xsubr[0,0],llcrnrlat=Ysubr[0,0],urcrnrlon=Xsubr[-1,-1],urcrnrlat=Ysubr[-1,-1], projection='merc',lat_0=latSubr_0,lon_0=lonSubr_0)



	Xsubr, Ysubr = m(Xsubr, Ysubr)

	m.contourf(Xsubr, Ysubr, bathySubr, vmin=vmin[1], vmax=vmax[1], cmap='plasma', levels=np.linspace(vmin[1],vmax[1],11), extend=extend[1])

	m.colorbar()

	

	plt.title(title[1], fontsize=fs)

	

	if parallels is not None:

		# Latitudes

		m.drawparallels(parallels[1],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])

	if meridians is not None:

		# Longitudes

		m.drawmeridians(meridians[1],labels=[True,False,False,True], color='k', linewidth=0.5,dashes=[2,1])



	if labelData[1] is not None:

		for li in labelData[1]:

			plt.scatter(li['x'][0], li['x'][1], s=1, color='r')

			plt.annotate(li['t'], li['tx'], color=li['tc'], fontsize=li['ts'])



	if isf[1] is not None:

		#x, y = m(isf[0], isf[1])  # transform coordinates

		#m.scatter(x, y, color='k', s=1)

		extent = [X[0,0], X[0,-1], -Y[0,0], -Y[-1,0]]

		m.imshow(1-isf[1], cmap=plt.cm.gray, interpolation='nearest', extent=extent, zorder=3);

		#m.contour(X, Y, isf, levels=[1], colors='k')

				

	#==

	

	plt.tight_layout()

	

	if save:

		plt.savefig(outpath + outname)

		

	if show:

		plt.show()



#==



def animate1by1quiverBasemap(u, v, Xd, Yd, lat_0, lon_0, width=2.e6, height=1.7e6, contourf=None, X=None, Y=None, figsize=(5,4), title='', fontsize=14, mesh=True, cmap='jet', vmin=None, vmax=None, xlabel='', ylabel='', save=True, outpath='', outname='animate1by1.mp4', show=False, dpi=200, fps=8, bitrate=-1, text_data=None, parallels=None, meridians=None):



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



def plot3D(X, Y, data):

	'''Plot data in 3 dimensional graph.'''



	ax = plt.figure().add_subplot(projection='3d')

	ax.contour(X, Y, data, extend3d=True)

	plt.show()

	

#==



def plotRads(RADS, X, Y, grid, bathy, subr=False, lats=None, lons=None):

	'''Plot six surface radiation terms in one figure.'''



	# Check if RADs is list. If so, do first data in first dict minus data in second.

	if isinstance(RADS, list):

		diff = True

		LATENT = RADS[0]['LATENT'] - RADS[1]['LATENT']

		SENS = RADS[0]['SENS'] - RADS[1]['SENS']

		BB = RADS[0]['BB'] - RADS[1]['BB']

		LW = RADS[0]['LW'] - RADS[1]['LW']

		SW = RADS[0]['SW'] - RADS[1]['SW']

		FC = RADS[0]['FC'] - RADS[1]['FC']

		FMI = RADS[0]['FMI'] - RADS[1]['FMI']

	# Otherwise, assume RADS is single dictionary.

	else:

		diff = False

		LATENT = RADS['LATENT']

		SENS = RADS['SENS']

		BB = RADS['BB']

		LW = RADS['LW']

		SW = RADS['SW']

		FC = RADS['FC']

		FMI = RADS['FMI']

		

	#RADS  = {'LATENT':LATENT_out, 'SENS':SENS_out, 'BB':BB_out, 'LW':LW_out, 'SW':SW_out, 'FC':FC_out,...} 

	

	LATENT = maskBathyXY(LATENT, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	LATENT = maskDraftXY(LATENT, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	

	SENS = maskBathyXY(SENS, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	SENS = maskDraftXY(SENS, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	

	BB = maskBathyXY(BB, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	BB = maskDraftXY(BB, grid, zi=0, subregion=subr, lats=lats, lons=lons)



	LW = maskBathyXY(LW, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	LW = maskDraftXY(LW, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	

	SW = maskBathyXY(SW, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	SW = maskDraftXY(SW, grid, zi=0, subregion=subr, lats=lats, lons=lons)



	FC = maskBathyXY(FC, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	FC = maskDraftXY(FC, grid, zi=0, subregion=subr, lats=lats, lons=lons)



	FMI = maskBathyXY(FMI, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	FMI = maskDraftXY(FMI, grid, zi=0, subregion=subr, lats=lats, lons=lons)

	

	#==

	

	if diff:

		SWvmin = -0.02; SWvmax = 0.02

		LWvmin = -12; LWvmax = 12

		BBvmin = -12.; BBvmax = 12.

		LATvmin = -3.; LATvmax = 3.

		SENSvmin = -3.; SENSvmax = 3.

		FCvmin = -3.; FCvmax = 3.

		FMIvmin = -12.; FMIvmax = 12.

	else:

		SWvmin = 0.04; SWvmax = 0.24

		#LWvmin = 190; LWvmax = 240

		#BBvmin = 220; BBvmax = 270

		LWvmin = 190; LWvmax = 270

		BBvmin = 190; BBvmax = 270	

		LATvmin = 0; LATvmax = 6

		SENSvmin = -12; SENSvmax = 12

		FCvmin = 4; FCvmax = 20

		FMIvmin = -12.; FMIvmax = 12.

				

	#data = [[-SW, -LW, BB], [-LATENT, -SENS, FC]] #-(LAT+SENS+BB+LW+SW)

	#titles = [[r'SW $\downarrow$ (W/m^2)', r'LW $\downarrow$', r'BB $\uparrow$'], [r'LATENT $\downarrow$', r'SENS $\downarrow$', r'FC $\uparrow$']]

	#vmin = [[SWvmin, LWvmin, BBvmin], [LATvmin, SENSvmin, FCvmin]]

	#vmax = [[SWvmax, LWvmax, BBvmax], [LATvmax, SENSvmax, FCvmax]]

	

	data = [[-LW, BB, FMI], [-LATENT, -SENS, FC]] #-(LAT+SENS+BB+LW+SW)

	titles = [[r'LW $\downarrow$ (W/m^2)', r'BB $\uparrow$', r'FMI $\uparrow$'], [r'LATENT $\downarrow$', r'SENS $\downarrow$', r'FC $\uparrow$']]

	vmin = [[LWvmin, BBvmin, FMIvmin], [LATvmin, SENSvmin, FCvmin]]

	vmax = [[LWvmax, BBvmax, FMIvmax], [LATvmax, SENSvmax, FCvmax]]



	

	#vmin = None; vmax = None

	

	bathyC = makeBathyContour(bathy, X, Y)

	

	plotMbyN(data, X=X, Y=Y, xticksvis=False, yticksvis=False, xticks=[[[]]*3]*2, yticks=[[[]]*3]*2, titles=titles, vmin=vmin, vmax=vmax, contour=bathyC, contourLevels=[-1000], save=True)

	

#==

	

	





