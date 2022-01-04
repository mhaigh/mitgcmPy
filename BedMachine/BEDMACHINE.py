# BEDMACHINE.py
# Read BedMachine netCDF file.

import numpy as np
import netCDF4 as nc

import matplotlib.pyplot as plt

#==========================================================

path = '/home/michai/Documents/data/BedMachine/'# '/Users/mh115/Documents/BAS/data/BedMachine/'
fname = 'BedMachineAntarctica_2020-07-15_v02.nc'
		
bed = 'bed'
mask = 'mask'
data = nc.Dataset(path+fname)

bed = data[bed][:][:,::-1]
mask = data[mask][:][:, ::-1]

xp = 400
xlims = [9250+xp,10750]
ylims = [7000+xp, 8500]
bed = bed[ylims[0]:ylims[1], xlims[0]:xlims[1]]
mask = mask[ylims[0]:ylims[1], xlims[0]:xlims[1]]

#variables(dimensions): |S1 mapping(), int32 x(x), int32 y(y), int8 mask(y, x), float32 firn(y, x), float32 surface(y, x), float32 thickness(y, x), float32 bed(y, x), int16 errbed(y, x), int8 source(y, x), int16 geoid(y, x)

bed = np.ma.array(bed, mask=mask>0)
lvls = [-2000, -1500, -1000, -500]

bed = bed.T[:, ::-1]

plt.figure(dpi=200)
plt.gca().patch.set_color('.25')
cax = plt.pcolormesh(bed, vmin=-1000, vmax=0, cmap='plasma')
plt.contour(bed, levels=lvls, colors='k', linestyles='solid', linewidths=0.8)
plt.title('Amundsen Sea bathymetry (units m)')
plt.colorbar(cax)
plt.savefig('fig1.png')
plt.show()
quit()

plt.contourf(mask)
plt.colorbar()
plt.show()
quit()





quit()




#print(bed)

		
