# io.py

# For reading mitgcm outputs in either binary on netCDF format.
# Reading of grid data is handled by grid.py

#==========================================================

import sys

import numpy as np
#import varDict

#==========================================================

# Issues here:
# open_mdsdataset gets all variables in dir -- not sure if there's a way to specify one variable.
# In PISOMIP_001/run on Archer2 I had to hard copy mit2nc before I could execute it.

def readnp(inputfilename, dims, dtype='>f', rec=-1):
	''''''
	
	if len(dims) == 2:
		ny, nx = dims
		count = nx * ny
	elif len(dims) == 3:
		nz, ny, nx = dims
		count = nx * ny * nz
	else: 
		print('Try again: dims should be length 2 or 3')
		quit() 
		
	size = np.dtype(dtype).itemsize
	f = open(inputfilename, 'rb')
				
	# Read all records
	if rec == -1:
		data = np.fromfile(f, dtype=dtype, count=-1)
		data = data.reshape((int(data.size/(ny*nx)),ny,nx))
		
	# Read specific record
	else:
		f.seek(rec*size*count)
		data = np.fromfile(f, dtype=dtype, count=count)
		data = data.reshape(dims)	
	
	return data
	
#==

def readAllnp(fileHandle, path, dims, dtype='>f', reverse=False):
	''''''

	import os
	fnames = [filename for filename in os.listdir(path) if filename.startswith(fileHandle) and filename.endswith('.data')]
	fnames.sort(reverse=reverse)

	if len(fnames) == 0:
		print('Error: readData.realAllnp. No files found for VAR ' + fileHandle)
		
	if len(dims) == 2:
		data = np.zeros((len(fnames), dims[0], dims[1]))
	elif len(dims) == 3:
		data = np.zeros((len(fnames), dims[0], dims[1], dims[2]))
	else:
		print('Error: readData.readAllnp. dims must be length 2 or 3.')
		sys.exit()
		
	for fi, fname in enumerate(fnames):
		data[fi] = readnp(path+fname, dims, rec=0)
	
	return data
	
#==

def readVariable(VAR, path, file_format='nc', time_step=1, var2D=False, meta=False, tt=None, xx=None, yy=None, zz=None):
	'''Read mitgcm output for variable in given file format.
	Options are rdmds or nc, with netCDF default.
	If meta==True, returns additional metadata useful for plotting.'''
			
	from varDict import varDict
	

	# First check if sub-intervals are needed.
	if tt is not None:
		if tt == -1:
			tt = [-1, None]
		elif not isinstance(tt, list):
			tt = [tt, tt+1]
	else:
		tt = [None, None]
		
	if xx is not None:
		if not isinstance(xx, list):
			xx = [xx, xx+1]
	else:
		xx = [None, None]
		
	if yy is not None:
		if not isinstance(yy, list):
			yy = [yy, yy+1]
	else:
		yy = [None, None]
		
	if zz is not None:
		if not isinstance(zz, list):
			zz = [zz, zz+1]
	else:
		zz = [None, None]	
		
	if file_format == 'rdmds':
		from MITgcmutils import rdmds
		import numpy as np
		fname = varDict[VAR]['fname']
		return rdmds(path+fname, itrs=np.NaN)
	
	elif file_format == 'nc':
		from netCDF4 import Dataset
		from numpy import squeeze as sq
		fname = varDict[VAR]['fname'] + '.nc'
		
		if meta:
			return Dataset(path+fname, 'r')
		else:
			if fname == 'state2D.nc':
				return sq(Dataset(path+fname, 'r')[VAR][tt[0]:tt[1],yy[0]:yy[1],xx[0]:xx[1]])
			else:
				return sq(Dataset(path+fname, 'r')[VAR][tt[0]:tt[1],zz[0]:zz[1],yy[0]:yy[1],xx[0]:xx[1]])
		
		
	else:
		print('Error: readData.readVariable. file_format must be rdmds or nc')
		sys.exit()
	

# Could later implement xmitgcm compatibility...	
#return open_mdsdataset(path, delta_t=60, read_grid=False)
