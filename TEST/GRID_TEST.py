
import numpy as np
import netCDF4 as nc

from ..grid_KN import Grid as Grid_KN
from ..grid import Grid

import matplotlib.pyplot as plt
import matplotlib.colors as cl

from MITgcmutils import rdmds


# Different variables defined at different parts of the grid/stencil.

# build_land_mask returns 2D boolean array of land.
##	 hfac is fraction of vertical cell occupied by water, so if zero, then its land. Full column=0 -> coastline.
# To get ice shelf, get ocean mask, multiply by (surface) regions where hfac < 1.

path = '/Users/mh115/Documents/BAS/data/PAS_666/'
grid = Grid(path)
grid_KN = Grid_KN(path) 



bathy = grid.bathy
draft = grid.draft

landC = grid.landC
iceC = grid.iceC




# Plot land and ice shelf.
if True:


	
	plt.subplot(121)
	plt.pcolormesh(draft)
	plt.colorbar()
	plt.subplot(122)
	plt.pcolormesh(grid_KN.draft-draft)
	plt.colorbar()
	plt.show()
	
	plt.subplot(121)
	plt.pcolormesh(bathy)
	plt.colorbar()
	plt.subplot(122)
	plt.pcolormesh(grid_KN.bathy-bathy)
	plt.colorbar()
	plt.show()
	
	quit()
	
	plt.subplot(121)
	plt.pcolormesh(draft)
	plt.colorbar()
	plt.subplot(122)
	plt.pcolormesh(grid_KN.draft)
	plt.colorbar()
	plt.show()
	
	quit()
	
	plt.subplot(121)
	plt.pcolormesh(grid.landC)
	plt.colorbar()
	plt.subplot(122)
	plt.pcolormesh(grid_KN.land_mask)
	plt.colorbar()
	plt.show()

	plt.subplot(121)
	plt.pcolormesh(grid.iceC)
	plt.colorbar()
	plt.subplot(122)       
	plt.pcolormesh(grid_KN.ice_mask)
	plt.colorbar()
	plt.show()       
