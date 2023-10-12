import numpy as np

from netCDF4 import Dataset



from grid import Grid

from grid_PAS import Grid as Grid_PAS



import matplotlib.pyplot as plt

import matplotlib.colors as cl



import plotting as pt

import plotting_tools as ptt



import tools



from readData import readVariable

from varDict import getPlottingVars



import time



#==========================================



path_root = '/data/oceans_output/shelf/pahol/mitgcm/'

run = 'PAS_851/run/'

# With PAS_851, Paul suggests looking at dates 1979-2012.

# Pre 1979 is spin up, post 2012 is cold anomaly.



path = path_root + run



# Load grid

grid = Grid_PAS(path)



lats = [-75.5, -70.5]; lons = [245, 262]#[230, 270]#

latsi = grid.getIndexFromLat(lats); lonsi = grid.getIndexFromLon(lons)



# Load UVEL.

u = readVariable('UVEL', path, file_format='nc', meta=True)



ts = 107; te = 502

ts = 24*12 + 5*12 - 8; te = -7



TIME = u['TIME'][:]

TIME = PAS_tools.getDecimalTime(TIME)[ts:te]



TIME = np.ma.filled(TIME, fill_value=0)

np.save('TIME_PAS851', TIME)



print(TIME.shape)

print(TIME[0])

print(TIME[-1])



WIND = False

if WIND:

	exf = Dataset(path+'stateExf.nc')

	uwind = np.mean(exf['EXFuwind'][ts:te+1], axis=0)

	vwind = np.mean(exf['EXFvwind'][ts:te+1], axis=0)

	uwind = np.ma.filled(uwind, fill_value=0)

	vwind = np.ma.filled(vwind, fill_value=0)

	np.save('uwindmean_PAS851', uwind)

	np.save('vwindmean_PAS851', vwind)



SSH = False

if SSH:

	ssh = Dataset(path+'state2D.nc')

	ssh = np.mean(ssh['ETAN'][ts:te+1], axis=0)

	ssh = np.ma.filled(ssh, fill_value=0)

	np.save('ssh_PAS851', ssh)

	

OCETAU = False

if OCETAU:



	state2D = Dataset(path+'state2D.nc')

	taux = np.mean(state2D['oceTAUX'][ts:te+1], axis=0)

	tauy = np.mean(state2D['oceTAUY'][ts:te+1], axis=0)

	taux = np.ma.filled(taux, fill_value=0)

	tauy = np.ma.filled(tauy, fill_value=0)

	np.save('tauxmean_PAS851', taux)

	np.save('tauymean_PAS851', tauy)



SALT = False

if SALT:

	# Repeat for SALT.

	S = readVariable('SALT', path, file_format='nc', meta=True)

	S1 = np.mean(S['SALT'][ts:te+1], axis=0)

	S1 = np.ma.filled(S1, fill_value=0)

	np.save('Smean_PAS851', S1)



	S = S['SALT'][:,:,latsi[0]:latsi[1]+1,lonsi[0]:lonsi[1]+1]

	print('S shape: ' + str(S.shape))

	Smean = np.mean(S[ts:te+1], axis=0)

	Smean = np.ma.filled(Smean, fill_value=0)

	np.save('Smean_PAS851', Smean)

	np.save('S0_PAS851', np.ma.filled(S[:,0], fill_value=0))

	quit()



UVEL = False

if UVEL:

	# Compute long-term mean.

	u = u['UVEL']#[:,:,latsi[0]:latsi[1]+1,lonsi[0]:lonsi[1]+1]

	print('u shape: ' + str(u.shape))

	umean = np.mean(u[ts:te+1], axis=0)

	umean = np.ma.filled(umean, fill_value=0)

	np.save('umean_PAS851', umean)

	#np.save('u0_PAS851', np.ma.filled(u[:,0], fill_value=0))



VVEL = False

if VVEL:

	# Repeat for VVEL.

	v = readVariable('VVEL', path, file_format='nc', meta=True)

	v = v['VVEL']#[:,:,latsi[0]:latsi[1]+1,lonsi[0]:lonsi[1]+1]

	print('v shape: ' + str(v.shape))

	vmean = np.mean(v[ts:te+1], axis=0)

	vmean = np.ma.filled(vmean, fill_value=0)

	np.save('vmean_PAS851', vmean)

	#np.save('v0_PAS851', np.ma.filled(v[:,0], fill_value=0))

	#w = readVariable('WVEL', path, file_format='nc', meta=False, interval=[ts,te])



THETA = False

if THETA:

	# Repeat for THETA.

	T = readVariable('THETA', path, file_format='nc', meta=True)

	T = T['THETA']#[:,:,latsi[0]:latsi[1]+1,lonsi[0]:lonsi[1]+1]

	print('T shape: ' + str(T.shape))



	THERM = None

	if THERM is not None:

		Z = grid.RC.squeeze()

		zi = grid.getIndexFromDepth(-1000)

		# Get z-indices of level with Theta closest to THERM.

		Tz = np.argmin(np.abs(T[:,:zi]-THERM), axis=1)

		ThermZ = Z[Tz]

		ThermZ = np.ma.filled(ThermZ, fill_value=0)

		np.save('Therm05Z_PAS851', ThermZ)	



	Tmean = np.mean(T[ts:te+1], axis=0)

	Tmean = np.ma.filled(Tmean, fill_value=0)

	np.save('Tmean_PAS851', Tmean)

	#np.save('T0_PAS851', np.ma.filled(T[:,0], fill_value=0))



#==



UVT_LEVEL = True

if UVT_LEVEL:



	level = 16

	

	u = u['UVEL'][ts:te,level,:,:]

	v = readVariable('VVEL', path, file_format='nc', meta=False, tt=[ts,te], zz=level)

	T = readVariable('THETA', path, file_format='nc', meta=False, tt=[ts,te], zz=level)

		

	print('u shape: ' + str(u.shape))

	print('v shape: ' + str(v.shape))

	

	u = np.ma.filled(np.mean(u, axis=0), fill_value=0)

	v = np.ma.filled(np.mean(v, axis=0), fill_value=0)

	T = np.ma.filled(np.mean(T, axis=0), fill_value=0)

	

	np.save('umean_PAS851_z'+str(level), u)

	np.save('vmean_PAS851_z'+str(level), v)

	np.save('Tmean_PAS851_z'+str(level), T)

	

#==



FLUXES = False

if FLUXES:

	# Heat fluxes.

	tauK = 0#273.15



	# First interpolate u,v to t points.

	u = tools.interp(u, 'u')

	v = tools.interp(v, 'v')



	uT = np.mean(u*(T+tauK), axis=0)

	vT = np.mean(v*(T+tauK), axis=0)

	#wT = np.mean(w*(T+tauK), axis=0)



	uT = np.ma.filled(uT, fill_value=0)

	vT = np.ma.filled(vT, fill_value=0)

	#wT = np.ma.filled(wT, fill_value=0)



	np.save('uT_PAS851', uT)

	np.save('vT_PAS851', vT)

	#np.save('wT_PAS851', wT)





