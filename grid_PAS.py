# grid.py

import sys

import numpy as np

from MITgcmutils import rdmds

#==========================================================

class Grid:

	'''
	Creates the grid from MITgcm binary files.
	----------
	ds : xarray.Dataset
		Dataset with gridinformation used to construct c-grid
	axis : str
		The appropriate xgcm axis. E.g. 'X' for longitudes.
	'''
	def __init__(self, path):
	
		# Read from MITgcm binary files.
		self.XC = rdmds(path+'XC')
		self.YC = rdmds(path+'YC')
		self.DXG = rdmds(path+'DXG')
		self.DYG = rdmds(path+'DYG')
		self.RAC = rdmds(path+'RAC')
		self.RC = rdmds(path+'RC')
		self.DRF = rdmds(path+'DRF')
		self.hFacC = rdmds(path+'hFacC')
		self.hFacS = rdmds(path+'hFacS')
		self.hFacW = rdmds(path+'hFacW')

		self.RF = self.RC + self.DRF/2
		self.RF = np.append(self.RF, [self.RF[-1]-self.DRF[-1]], 0)
		
		# Land boolean arrays.
		self.landC = np.sum(self.hFacC, axis=0) == 0
		self.landS = np.sum(self.hFacS, axis=0) == 0
		self.landW = np.sum(self.hFacW, axis=0) == 0
		
		# Ice shelf boolean arrays
		self.iceC = np.logical_and(np.sum(self.hFacC, axis=0) != 0, self.hFacC[0,] < 1)
		self.iceS = np.logical_and(np.sum(self.hFacS, axis=0) != 0, self.hFacS[0,] < 1)
		self.iceW = np.logical_and(np.sum(self.hFacW, axis=0) != 0, self.hFacW[0,] < 1)

		# Set draft everywhere to nan, to represent land or ice shelf.
		# For each level k, moving downwards, get grid points which are not entirely land AND draft is nan.
		# For these grid points at level k, set draft to z  - dz * (1-hFaC).  
		# Above will set draft to zero for all open ocean cells (hFacC = 1) for k=0.
		# At an ice shelf, first wet cell (hFacC > 0) will be given non-nan draft.
		# After all levels, aself.Nywhere with draft = nan has to be land.

		# Ice shelf draft
		self.Nz, self.Ny, self.Nx = self.hFacC.shape
		draft = np.zeros((self.Ny, self.Nx)) * np.nan
		for k in range(self.Nz):
			# Get indices of wet cells with draft not already computed.
			index = (self.hFacC[k,] != 0) * np.isnan(draft)
			# Set draft  = top of cell - dz * (fraction of cell not water).
			draft[index] = self.RF[k,0,0] - self.DRF[k,0,0] * (1 - self.hFacC[k][index])
    		# Set remaining NaNs to zero, representing land.
		self.draft = np.nan_to_num(draft)
		
		# Bathymetry (sum cell heights scaled by land fraction, subtract depth, subtract ice shelf contribution).
		self.bathy = np.sum(self.DRF * (1-self.hFacC), axis=0) + self.RF[-1,] + self.draft

		# Ice shelf draft
		self.Nz, self.Ny, self.Nx = self.hFacC.shape
		draft = np.zeros((self.Ny, self.Nx)) * np.nan
		for k in range(self.Nz):
			# Get indices of wet cells with draft not already computed.
			index = (self.hFacC[k,] != 0) * np.isnan(draft)
			# Set draft  = top of cell - dz * (fraction of cell not water).
			draft[index] = self.RF[k,0,0] - self.DRF[k,0,0] * (1 - self.hFacC[k][index])
    		# Set remaining NaNs to zero, representing land.
		self.draft = np.nan_to_num(draft)

		# Bathymetry (sum cell heights scaled by land fraction, subtract depth, subtract ice shelf contribution).
		self.bathy = np.sum(self.DRF * (1-self.hFacC), axis=0) + self.RF[-1,] + self.draft

	# Functions to get nearest index for given latitude, longitude, depth.
	# May need to be extended so is compatible with grid box corners.
	
	def getIndexFromLat(self, lats, xi=0):
		'''Given latitude in degrees, return closest grid point index.'''
		
		if isinstance(lats, list):
			return [np.argmin(np.abs(self.YC[:, xi] - lats[j])) for j in range(len(lats))]
			
		else:
			return np.argmin(np.abs(self.YC[:, xi] - lats))
			
	#==
	
	def getIndexFromLon(self, lons, yi=0):
		'''Given longitude in degrees, return closest grid point index.'''
		
		if isinstance(lons, list):
			return [np.argmin(np.abs(self.XC[yi, :] - lons[i])) for i in range(len(lons))]
			
		else:
			return np.argmin(np.abs(self.XC[yi, :] - lons))
		
	#==
	
	def getIndexFromDepth(self, depths):
		'''Given depth return closest grid point index.'''
		
		if isinstance(depths, list):
			return [np.argmin(np.abs(self.RC.squeeze() - depths[k])) for k in range(len(depths))]
			
		else:
			return np.argmin(np.abs(self.RC.squeeze() - depths))
	
			
	#==
	
	def XYsubr(self, lons, lats, data=None):
		'''Return 1D array of subregion of longitudes.
		Input lons should be list/tuple with two entries.'''
		
		if len(lons) != 2 or len(lats) != 2:
			print('Error: grid.Xsubr. lats and lons must have two entries.')
			sys.exit()
			
		i0 = 0

		ys = np.argmin(np.abs(self.YC[:, i0] - lats[0]))
		yn = np.argmin(np.abs(self.YC[:, i0] - lats[1]))

		iw = np.argmin(np.abs(self.XC[i0, :] - lons[0]))
		ie = np.argmin(np.abs(self.XC[i0, :] - lons[1]))
		
		if data is None:
			return self.XC[ys:yn+1, iw:ie+1].copy(), self.YC[ys:yn+1, iw:ie+1].copy()
		else:
			return self.XC[ys:yn+1, iw:ie+1].copy(), self.YC[ys:yn+1, iw:ie+1].copy(), data[ys:yn+1, iw:ie+1].copy()
	
	#==
	
	def Xsubr1D(self, lons, yi=0):
		'''Return 1D array of subregion of longitudes.
		Input lons should be list/tuple with two entries.'''
		
		if len(lons) != 2:
			print('Error: grid.Xsubr1D. lons must have two entries.')
			sys.exit()
			
		iw = np.argmin(np.abs(self.XC[yi, :] - lons[0]))
		ie = np.argmin(np.abs(self.XC[yi, :] - lons[1]))
		
		return self.XC[yi, iw:ie+1].copy()
		
	#==

	def Ysubr1D(self, lats, xi=0):
		'''Return 1D array of subregion of latitudes.
		Input lats should be list/tuple with two entries.'''
		
		if len(lats) != 2:
			print('Error: grid.Ysubr1D. lons must have two entries.')
			sys.exit()
			
		js = np.argmin(np.abs(self.YC[:, xi] - lats[0]))
		jn = np.argmin(np.abs(self.YC[:, xi] - lats[1]))
		
		return self.YC[js:jn+1, xi].copy()
		
		#==

	def Zsubr1D(self, depths):
		'''Return 1D array of subregion of depths.
		Input depths should be list/tuple with two entries.'''
		
		if len(depths) != 2:
			print('Error: grid.Zsubr1D. depths must have two entries.')
			sys.exit()
			
		kt = np.argmin(np.abs(self.RC.squeeze() - depths[0]))
		kb = np.argmin(np.abs(self.RC.squeeze() - depths[1]))
		
		return self.RC.squeeze()[kt:kb+1].copy()
	
	def getZFfromZC(self):
		
		zc = self.RC[:].squeeze().copy()
	

		
		

