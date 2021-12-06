# plotting.py

import sys

import numpy as np

import matplotlib.pyplot as plt

#==========================================================

def maskBathyXY(data, grid, zi, color='grey', subregion=False, lats=[], lons=[]):
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
		
	return np.ma.array(data, mask=z<bathy)
	
#==

def maskBathyXZ(data, grid, color='grey', yi=0, subregion=False, lons=[], depths=[]):
	'''Function to be called before plotting to mask the bathymetry in (X,Z)-slice.'''
	
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
		
	return np.ma.array(data, mask=z<bathy)
	
	
#==

def maskDraftYZ(data, grid, color='grey', xi=10, subregion=False, lats=[], depths=[]):
	'''Function to be called before plotting to mask the ice shelf draft.
	If subregion is True, lats and depths lets grid know what the subregion is.'''
	
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
	
	return np.ma.array(data, mask=z>draft)
	

#==

def maskBathyYZ(data, grid, color='grey', xi=0, subregion=False, lats=[], depths=[]):
	'''Function to be called before plotting to mask the bathymetry in (Y,Z)-slice.'''
	
	# If masking in subregion, grid needs to account for it.
	if subregion == True:
		try:
			RC = grid.RC.squeeze()[depths[0]:depths[1]+1]
			YC = grid.YC[lats[0]:lats[1]+1,]
			bathy = grid.bathy[lats[0]:lats[1]+1,]
			
		except ValueError:
			print('If subregion set to True, both lats and depths need to be defined.')
			
	else:
		RC = grid.RC.squeeze()
		YC = grid.YC
		bathy = grid.bathy
	
	# Make z a depth array with (y,z)-dimensions.
	z = np.tile(RC, (YC.shape[0], 1)).T
	
	# Make draft array with (y,z)-dimensions.
	bathy = np.tile(bathy[:,xi], (RC.shape[0], 1))
		
	return np.ma.array(data, mask=z<bathy)
	

	
