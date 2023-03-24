# BEDMACHINE.py
# Read BedMachine netCDF file.

import numpy as np
import netCDF4 as nc

import matplotlib.pyplot as plt

#==========================================================

fname = 'IBCSO_v2_bed_WGS84.nc'
		
bed = 'z'
data = nc.Dataset(fname)

print(data['lon'])

bed = data[bed][:][:,:]
x = data['lon'][:]# / 1000.
y = data['lat'][:]# / 1000.

print(x.shape)
#plt.plot(y[1:-1] - y[0:-2]); plt.show(); quit()

xp = 400
#xlims = [9250+xp,10750]
#ylims = [7000+xp, 8500]

xlims = [6000, 7700]
ylims = [1400, 1800]

#plt.plot(y); plt.grid(); plt.show();# quit()

bed = bed[ylims[0]:ylims[1], xlims[0]:xlims[1]]
x = x[xlims[0]:xlims[1]]
y = y[ylims[0]:ylims[1]]

bed = np.where(bed>0, 0, bed)
bed = np.where(bed<-1000, -1000, bed)

plt.contourf(x, y, bed, cmap='plasma')# vmin=-500, vmax=-400) 
plt.title('IBCSO V2 Amundsen Sea Bathymetry (m)')
plt.colorbar(); plt.show()
quit()

#print(bed)

		
