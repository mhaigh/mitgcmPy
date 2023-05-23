import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import plotting as pt
import plotting_tools as ptt
import tools 

from scipy.io import loadmat

#===========================================================================================


# Get PAS bathymetry for background	
#PASpath = '/home/michael/Documents/data/IdealisedAmundsenSea/data/PAS_851/run/'
#grid = Grid_PAS(PASpath)
#bathy = grid.bathy
#bathy = ptt.maskBathyXY(bathy, grid, zi=0)	
#land = np.where(bathy<0, 0, 1)

windPath = '/home/michael/Documents/data/WINDS/'
PASpath = '/home/michael/Documents/data/IdealisedAmundsenSea/data/PAS_851/run/'

#==

# Create lat and lon grids

nx = 1440; ny = 241; dims = (ny, nx)

uwind_lon0    = 0.0; uwind_lon_inc = 0.25
uwind_lat0    = -90.0; uwind_lat_inc = 0.25

lon = np.linspace(uwind_lon0, 360., nx)
lat = np.linspace(uwind_lat0, uwind_lat0+uwind_lat_inc*ny, ny)

# Zonal and meridional limits.
W = 360-130; E = 360-90;
S = -76; N = -64;

xw = np.argmin(np.abs(lon-W))
xe = np.argmin(np.abs(lon-E))
ys = np.argmin(np.abs(lat-S))
yn = np.argmin(np.abs(lat-N))
		
lat = lat[ys:yn]
lon = lon[xw:xe]

nySubr = lat.shape[0]
nxSubr = lon.shape[0]

#==

# Get PAS bathymetry for background	
grid = Grid_PAS(PASpath)
bathy = grid.bathy
bathy = ptt.maskBathyXY(bathy, grid, zi=0)	
land = np.where(bathy<0, 0, 1)

# Get PAS grid and land subrange.
lonPAS, latPAS, land = grid.XYsubr([W,E], [S,N], land)

print(bathy.shape)

# We want to show only xthe -1000m contour of bathy along the shelf break.
lonPASbathy, latPASbathy, bathy = grid.XYsubr([W,E], [-73,-70], bathy)

print(lonPASbathy.shape, latPASbathy.shape, bathy.shape)

Ny, Nx = latPAS.shape
latPAS_ = np.zeros((Ny+1,Nx))
latPAS_[1:,:] = latPAS; latPAS_[0,:] = -78
lonPAS_ = np.zeros((Ny+1,Nx))
lonPAS_[1:,:] = lonPAS; lonPAS_[0,:] = lonPAS[1,:]
land_ = np.zeros((Ny+1,Nx))
land_[1:,:] = land; land_[0,:] = 1

colors = ['white', 'cyan']; bounds = [0,0.5,1]
cmap = mpl.colors.ListedColormap(colors)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

#==

# Load winds
u = np.load(windPath + 'ERA5_uwind_av.npy')
v = np.load(windPath + 'ERA5_vwind_av.npy')

#==

# Readjust meridional limits.
S = -76; N = -66;

ys = np.argmin(np.abs(lat-S))
yn = np.argmin(np.abs(lat-N))
lat = lat[ys:yn]
u = u[ys:yn]
v = v[ys:yn]

#==

rhoa = 1.3 # kg / m3
ca = 0.0001
tauu = rhoa * ca * (u**2 + v**2)**0.5 * u 
tauv = rhoa * ca * (u**2 + v**2)**0.5 * v

u_av = np.mean(u, axis=1)
tauu_av = np.mean(tauu, axis=1)
tauRef = 0.025 * tools.get05cosWind(u.shape[1], u.shape[0])[:,0]

tauRef = 15 * tauRef / np.max(np.abs(tauRef)) + 250
u_av = 15 * u_av / np.max(np.abs(u_av)) + 250

#==

# PLOT

d = 4
latd = lat[::d]
lond = lon[::d]
ud = u[::d,::d]
vd = v[::d,::d]
tauud = tauu[::d,::d]
tauvd = tauv[::d,::d]
		
plt.contourf(lonPAS_, latPAS_, land_, cmap=cmap)
plt.contour(lonPASbathy, latPASbathy, bathy, levels=[-1000], colors='k', linestyles='solid')
plt.quiver(lond, latd, ud, vd)
plt.plot(tauRef, lat, label='Idealised wind (stress)')
plt.plot(u_av, lat, label='Zonal-mean wind')
plt.axvline(250, linestyle='dashed')
plt.legend()
	
plt.grid()
plt.xlim([W, E]); plt.ylim([S, N])
plt.title('ERA5 1979-2021 mean wind')
plt.show()
		
quit()

#####################################################################################################

# Code below here is for producing the mean wind file from ERA5 data on Archer2.

#####################################################################################################
#####################################################################################################

quit()

import os

import numpy as np

#===========================================================================================

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

outpath = ''

path = '/home/michael/Documents/data/WINDS/'
fileHandle = 'ERA5_uwind'
fnames = [filename for filename in os.listdir(path) if filename.startswith(fileHandle)]

# 6 hourly winds.
uwindstartdate1=19790101,
uwindstartdate2=000000,
uwindperiod=21600.0,

nx = 1440; ny = 241
dims = (ny, nx)

uwind_lon0    = 0.0
uwind_lon_inc = 0.25
uwind_lat0    = -90.0
uwind_lat_inc = 0.25

lon = np.linspace(uwind_lon0, 360., nx)
lat = np.linspace(uwind_lat0, uwind_lat0+uwind_lat_inc*ny, ny)

# Zonal and meridional limits.
W = 360-130; E = 360-90;
S = -76; N = -64;

xw = np.argmin(np.abs(lon-W))
xe = np.argmin(np.abs(lon-E))
ys = np.argmin(np.abs(lat-S))
yn = np.argmin(np.abs(lat-N))
		
lat = lat[ys:yn]
lon = lon[xw:xe]

nySubr = lat.shape[0]
nxSubr = lon.shape[0]

#==
	
data = np.zeros((len(fnames), nySubr, nxSubr))
for fi, fname in enumerate(fnames):
	data_tmp = readnp(path+fname, dims, rec=-1)
	data[fi] = np.mean(data_tmp[:,ys:yn, xw:xe], axis=0)
	
data = np.mean(data, axis=0)
np.save(fileHandle+'_av', data)



