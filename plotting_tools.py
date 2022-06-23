# plotting.py

import sys

import numpy as np

import matplotlib.pyplot as plt

#==========================================================

def maskBathyAll(data, grid, color='grey', timeDep=False):
	'''Mask bathymetry everywhere in data field.
	If timeDep=True, assumes data is time-dependent with
	time on the first axis.'''

	# Make bathy and z 3D fields.
	bathy = grid.bathy
	z = grid.RC
	bathy = np.tile(bathy, (z.shape[0], 1, 1))
	z = np.tile(z, (1, bathy.shape[1], bathy.shape[2]))

	if timeDep:
		print('Note: creating 4D arrays for masking can be slow.')
		Nt = data.shape[0]
		bathy = np.tile(bathy, (Nt,1,1,1))
		z = np.tile(bathy, (Nt,1,1,1))

	return np.ma.array(data, mask=grid.RC<bathy)

#==

def maskBathyAllT(data, grid):
	'''Calls above function in temporal for loop, to get around high mem. demands.'''

	Nt = data.shape[0]
	for ti in range(Nt):
		data[ti] = maskBathyAll(data[ti], grid)

	return data

#==

def maskBathyXY(data, grid, zi, color='grey', subregion=False, lats=[], lons=[], timeDep=False):
	'''Function to be called before plotting to mask the bathymetry in (X,Y)-slice.'''
		
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			XC = grid.XC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			bathy = grid.bathy[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
		except ValueError:
			print('Error: plotting_tools.maskBathyXY. If subregion set to True, both lons and lats need to be defined.')
			
	else:
		XC = grid.XC
		YC = grid.YC
		bathy = grid.bathy
	
	# The depth of the slice.
	z = grid.RC.squeeze()[zi]

	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		bathy = np.tile(bathy, (NT, 1, 1))
		
	return np.ma.array(data, mask=z<bathy)
	
#==

def maskBathyXZ(data, grid, color='grey', yi=0, subregion=False, lons=[], depths=[], timeDep=False):
	'''Function to be called before plotting to mask the bathymetry in (X,Z)-slice.'''
	
	if (len(data.shape) != 2 and timeDep == False):
		print('Error: plotting_toolts.maskBathyXZ. If data is more than 2-dimensional, timeDep must be set to True.')
		sys.exit()
		
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			XC = grid.XC[:, lons[0]:lons[1]+1]
			bathy = grid.bathy[:, lons[0]:lons[1]+1]
			
		except ValueError:
			print('Error: plotting_tools.maskBathyXZ. If subregion set to True, both lons and depths need to be defined.')
			
	else:
		RC = grid.RC.squeeze()
		XC = grid.XC
		bathy = grid.bathy
	
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (XC.shape[1], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	bathy = np.tile(bathy[yi, ], (RC.shape[0], 1))
	
	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		bathy = np.tile(bathy, (NT, 1, 1))
		
	return np.ma.array(data, mask=z<bathy)
	
	
#==

def maskBathyYZ(data, grid, color='grey', xi=0, subregion=False, lats=[], depths=[], timeDep=False):
	'''Function to be called before plotting to mask the bathymetry in (Y,Z)-slice.'''
	
	if (len(data.shape) != 2 and timeDep == False):
		print('Error: plotting_toolts.maskBathyYZ. If data is more than 2-dimensional, timeDep must be set to True.')
		sys.exit()
	
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1,]
			bathy = grid.bathy[lats[0]:lats[1]+1,]
			
		except ValueError:
			print('Error: plottingtools.maskBathyYZ. If subregion set to True, lats and depths need to be defined.')
			
	else:
		RC = grid.RC.squeeze()
		YC = grid.YC
		bathy = grid.bathy
	
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (YC.shape[0], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	bathy = np.tile(bathy[:,xi], (RC.shape[0], 1))

	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		bathy = np.tile(bathy, (NT, 1, 1))
		
	return np.ma.array(data, mask=z<bathy)

#==

def maskBathyYZ_allX(data, grid, timeDep=False):
	'''Calls maskBathyYZ for a range of x gridpoints.
	Makes the assumption that x is on the final axis.
	If timeDep this is axis=3, if not timeDep this is axis=2.'''

	data_masked = data.copy()

	if timeDep:
		Nx = data.shape[3]
	else:
		Nx = data.shape[2]

	for xi in range(Nx):
		data_masked[..., xi] = maskBathyYZ(data[..., xi], grid, xi=xi, timeDep=timeDep)

	return data_masked
	
#==


def maskDraftXY(data, grid, zi, color='grey', subregion=False, lats=[], lons=[], timeDep=False):
	'''Function to be called before plotting to mask the draft in (X,Y)-slice.'''
		
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			XC = grid.XC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
			draft = grid.draft[lats[0]:lats[1]+1, lons[0]:lons[1]+1]
		except ValueError:
			print('Error: plotting_tools.maskBathyXY. If subregion set to True, both lons and lats need to be defined.')
			
	else:
		XC = grid.XC
		YC = grid.YC
		draft = grid.draft
	
	# The depth of the slice.
	z = grid.RC.squeeze()[zi]

	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		draft = np.tile(draft, (NT, 1, 1))
		
	return np.ma.array(data, mask=z>draft)

#==

def maskDraftYZ(data, grid, color='grey', xi=10, subregion=False, lats=[], depths=[], timeDep=False):
	'''Function to be called before plotting to mask the ice shelf draft.
	If subregion is True, lats and depths lets grid know what the subregion is.'''
	
	if (len(data.shape) != 2 and timeDep == False):
		print('Error: plotting_toolts.maskDraftYZ. If data is more than 2-dimensional, timeDep must be set to True.')
		sys.exit()
	
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1,]
			draft = grid.draft[lats[0]:lats[1]+1,]
			
		except:
			print('Error: plottingtools.maskDraftYZ. If subregion set to True, lats and depths need to be defined.')
			sys.exit()
			
	else:  
		RC = grid.RC.squeeze()
		YC = grid.YC
		draft = grid.draft
			
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (YC.shape[0], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	draft = np.tile(draft[:,xi], (RC.shape[0], 1))
	
	# If data is time-dependent, assume time dimension is first and tile draft and z arrays.
	if timeDep:
		NT = data.shape[0]
		z = np.tile(z, (NT, 1, 1))
		draft = np.tile(draft, (NT, 1, 1))
	
	return np.ma.array(data, mask=z>draft)
	

#==

def getTextData(TIME, t_format, xloc, yloc, color='k', short_ctime=True):
	'''Get text data for animation plots.'''

	if t_format == 'month':
		tscale_val = 86400. * 30
		t = TIME / tscale_val
		text = [t_format + ' ' + str(int(tt)) for tt in t]

	elif t_format == 'ctime':
		from time import ctime
		text = [ctime(TIME[:][ti]) for ti in range(len(TIME[:]))]
		
		if short_ctime:
			text = [t[4:8] + t[-4:] for t in text]

	elif t_format == 'X':
		text = ['x = ' + str(time_/1000) + ' km' for time_ in TIME]

	text_data = {'text':text, 'xloc':xloc, 'yloc':yloc, 'fontdict':{'fontsize':14, 'color':color}}
	
	return text_data

#==

def setText(ax, text_data, i=None, set_invisible=False):
	'''Utility function for setting text on axis given text_data dictionary.'''

	if set_invisible:
		for t in ax.texts:
			t.set_visible(False)
	
	if i is None:
		ax.text(text_data['xloc'], text_data['yloc'], text_data['text'], fontdict=text_data['fontdict'])	
	else:
		ax.text(text_data['xloc'], text_data['yloc'], text_data['text'][i], fontdict=text_data['fontdict'])	

#==

def getContourfLevels(vmin, vmax, contourfNlevels):
	'''Get kwarg levels for contourf plot from given vmin and vmax.
	vmin and vmax can be None, in which case None is returned.'''

	if vmin == None or vmax == None:
		return None
	
	else:
		return np.linspace(vmin, vmax, contourfNlevels)
		
#==

def makeList(arg, m, n=None):
	'''Check if argument is list.
	If not list, make it a list of repeated argument.'''

	if not isinstance(arg, list):
		if n is None:
			arg = [arg]*m
		else:
			arg = [[arg]*m]*n
			
	return arg

#==

def doTicks(xticks, xticksvis, yticks, yticksvis):

	if xticks is not None:
		if xticksvis:				
			plt.xticks(xticks)
		else:			
			plt.xticks(xticks, labels='')
	if yticks is not None:
		if yticksvis:				
			plt.yticks(yticks)
		else:			
			plt.yticks(yticks, labels='')
			
#==

def doLabels(xlabel, ylabel, fontsize=12):

	if xlabel is not None:
		plt.xlabel(xlabel, fontsize=fontsize)
	if ylabel is not None:
		plt.ylabel(ylabel, fontsize=fontsize)
		
#==

def doTitle(title, fontsize=12):

	if title is not None:
		plt.title(title, fontsize=fontsize)
				
	
