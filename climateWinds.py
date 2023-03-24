import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import plotting as pt
import plotting_tools as ptt

from scipy.io import loadmat

#===========================================================================================

windPath = '/home/michael/Documents/data/WINDS/'
PASpath = '/home/michael/Documents/data/IdealisedAmundsenSea/data/PAS_851/run/'

# Zonal and meridional limits.
W = -130; E = -90;
S = -76; N = -64;

# Get PAS bathymetry for background	
grid = Grid_PAS(PASpath)
bathy = grid.bathy
bathy = ptt.maskBathyXY(bathy, grid, zi=0)	
land = np.where(bathy<0, 0, 1)

# Get PAS grid and land subrange.
lonPAS, latPAS, land = grid.XYsubr([W+360,E+360], [S,N], land); lonPAS -= 360

# We want to show only xthe -1000m contour of bathy along the shelf break.
lonPASbathy, latPASbathy, bathy = grid.XYsubr([W+360,E+360], [-73,-70], bathy); lonPASbathy -= 360
#pt.plot1by1(bathy, X=lonPASbathy, Y=latPASbathy); quit()

# Extend land array southwards, to match wind data
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

#extent = [lonPAS[0,0], lonPAS[0,-1], latPAS[0,0], latPAS[-1,0]]
extent = [W, E, S, N]

# Load Pauls' wind data.
f = loadmat(windPath + 'paneldata.mat')['paneldata']

d = [1,1,1,1,4,1]
titles = ['reconstruction', 'external', 'internal', 'GHG', 'IPO', 'ozone']

for pi in range(6):

	lat = f[pi,0]['lat'][0,0]
	lon = f[pi,0]['lon'][0,0]
	u = f[pi,0]['u'][0,0]
	v = f[pi,0]['v'][0,0]
	
	xw = np.argmin(np.abs(lon[0,:]-W))
	xe = np.argmin(np.abs(lon[0,:]-E))
	ys = np.argmin(np.abs(lat[:,0]-S))
	yn = np.argmin(np.abs(lat[:,0]-N))

	lat = lat[ys:yn, xw:xe]
	lon = lon[ys:yn, xw:xe]
	u = u[ys:yn, xw:xe]
	v = v[ys:yn, xw:xe]
	
	latd = lat[::d[pi],::d[pi]]
	lond = lon[::d[pi],::d[pi]]
	ud = u[::d[pi],::d[pi]]
	vd = v[::d[pi],::d[pi]]

	# Colour in land.
	
	uav = np.mean(u, axis=1)
	uav = uav * 0.3*(E-W) / np.max(np.abs(uav)) + (W+E)/2
	zeros = np.zeros(uav.shape[0]) + (W+E)/2
	
	#plt.imshow(land[::-1], interpolation='none', cmap=cmap, norm=norm, extent=extent, aspect=2)
	plt.plot(uav, lat[:,0], color='r')
	plt.plot(zeros, lat[:,0], linestyle='dashed', color='r')
	plt.contourf(lonPAS_, latPAS_, land_, cmap=cmap)
	#plt.contourf(lonPAS, latPAS, land, cmap=cmap)
	plt.contour(lonPASbathy, latPASbathy, bathy, levels=[-1000], colors='k', linestyles='solid')
	plt.quiver(lond, latd, ud, vd)
	
	plt.title(titles[pi])
	plt.xlim([W, E]); plt.ylim([S, N])
	plt.show()
	
	
	
