# io.py

# For reading mitgcm outputs in either binary on netCDF format.
# Reading of grid data is handled by grid.py

#==========================================================

import sys
#import varDict

#==========================================================

# Issues here:
# open_mdsdataset gets all variables in dir -- not sure if there's a way to specify one variable.
# In PISOMIP_001/run on Archer2 I had to hard copy mit2nc before I could execute it.


def readVariable(VAR, path, file_format='nc', time_step=1, meta=False, interval=None):
	'''Read mitgcm output for variable in given file format.
	Options are rdmds or nc, with netCDF default.
	If meta==True, returns additional metadata useful for plotting.'''
			
	from varDict import varDict

	if file_format == 'rdmds':
		from MITgcmutils import rdmds
		import numpy as np
		fname = varDict[VAR]['fname']
		return rdmds(path+fname, itrs=np.NaN)
	
	elif file_format == 'nc':
		from netCDF4 import Dataset
		fname = varDict[VAR]['fname'] + '.nc'
		if meta:
			return Dataset(path+fname, 'r')
		else:
			if interval is not None:
				return Dataset(path+fname, 'r')[VAR][interval[0]:interval[1]+1,]
			else:
				return Dataset(path+fname, 'r')[VAR][:]
		
	else:
		print('Error: readData.readVariable. file_format must be rdmds or nc')
		sys.exit()
	

# Could later implement xmitgcm compatibility...	
#return open_mdsdataset(path, delta_t=60, read_grid=False)
