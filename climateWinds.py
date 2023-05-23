import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from grid import Grid
from grid_PAS import Grid as Grid_PAS

import plotting as pt
import plotting_tools as ptt
import tools 

from scipy.io import loadmat
from scipy.interpolate import interp1d

#===========================================================================================

windPath = '/home/michael/Documents/data/WINDS/'
PASpath = '/home/michael/Documents/data/IdealisedAmundsenSea/data/PAS_851/run/'

rhoa = 1.3 # kg / m3
ca = 0.0053
# ca values:
# 0053 Bett et al.
# 0.001 Holland et al. 2019
# 1.25e-3 Dotto et al. 2020
# 1.2e-3 u<11, (0.49+0.65u)e-3 11<u<25 m/s

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

#==

# Load Pauls' wind data.
f = loadmat(windPath + 'paneldata.mat')['paneldata']

d = [1,1,1,1,4,1]
titles = ['reconstruction', 'external', 'internal', 'GHG', 'IPO', 'ozone']

#==

INTERP = 1
if INTERP:

	nyears = 100
	
	def interp(data):
		cubic = interp1d(lat, data, kind='cubic')
		return cubic(latMCS)
			
	for pi in [1,4]:
	
		lat = f[pi,0]['lat'][0,0]
		lon = f[pi,0]['lon'][0,0]
		u = f[pi,0]['u'][0,0]
		v = f[pi,0]['v'][0,0]
		
		xw = np.argmin(np.abs(lon[0,:]-W))
		xe = np.argmin(np.abs(lon[0,:]-E))
		ys = np.argmin(np.abs(lat[:,0]-S))
		yn = np.argmin(np.abs(lat[:,0]-N))

		lat = lat[ys:yn, 0]
		lon = lon[0, xw:xe]
		u = u[ys:yn, xw:xe]
		v = v[ys:yn, xw:xe]
		
		nyMCS = 200
		latMCS = np.linspace(min(lat), max(lat), nyMCS)
		
		print(u.shape)
		
		u *= nyears
		v *= nyears

		#==
		
		# Things to do. Which order?
		# 1. Get Paul's wind trends 
		# 2. Get MCS reference stress, turn into wind.
		# 3. Interpolate onto MCS grid
		# 4. Take zonal average.
		# 5. Account for meridional trends in stress computation?
		
		# Compute wind stresses and wind speeds

		# Get reference wind stress on climate model grid and MCS grid.		
		tauRef = 0.025 * tools.get05cosWind(u.shape[1], u.shape[0])
		tauRefMCS = 0.025 * tools.get05cosWind(u.shape[1], 200)
		
		# Convert idealised ref. stress into wind speed.
		uRef = np.abs(tauRefMCS/(rhoa*ca))**0.5 * np.sign(tauRefMCS) 
		
		# Interpolate wind trend in y-direction onto MCS grid. Keep x-direction as climate model grid.
		u_interp = np.zeros(uRef.shape)
		v_interp = np.zeros(uRef.shape)
		for xi in range(u_interp.shape[1]):
			u_interp[:,xi] = interp(u[:,xi])
			v_interp[:,xi] = interp(v[:,xi])
		uPos = uRef + u_interp; uNeg = uRef - u_interp
		
		# Average winds first?
		AV_FIRST = False
		if AV_FIRST:

			uPos = np.mean(uPos, axis=1)
			uNeg = np.mean(uNeg, axis=1)
			
			lim = 1.; uPosAbs = np.abs(uPos); uNegAbs = np.abs(uNeg)
			#uPosAbs = np.where(np.abs(uPos)<lim, lim, np.abs(uPos)) 
			#uNegAbs = np.where(np.abs(uNeg)<lim, lim, np.abs(uNeg)) 
			
			tauPos = ca * rhoa * np.abs(uPosAbs) * uPos
			tauNeg = ca * rhoa * np.abs(uNegAbs) * uNeg
		
		else:
			
			lim = 1.; uPosAbs = np.abs(uPos); uNegAbs = np.abs(uNeg)
			#uPosAbs = np.where(np.abs(uPos)<lim, lim, np.abs(uPos)) 
			#uNegAbs = np.where(np.abs(uNeg)<lim, lim, np.abs(uNeg)) 
			#uPosAbs = (uPos**2 + v_interp**2)**0.5
			#uNegAbs = (uPos**2 + v_interp**2)**0.5
						
			tauPos = ca * rhoa * np.abs(uPosAbs) * uPos
			tauNeg = ca * rhoa * np.abs(uNegAbs) * uNeg

			uPos = np.mean(uPos, axis=1)
			uNeg = np.mean(uNeg, axis=1)
			tauPos = np.mean(tauPos, axis=1)
			tauNeg = np.mean(tauNeg, axis=1)
		
		#==
		
		uRef = uRef[:,0]
		u_interp = np.mean(u_interp, axis=1)
		tauRefMCS = tauRefMCS[:,0]
		tauPos = tools.smooth3(tauPos)
		tauNeg = tools.smooth3(tauNeg)
		
		outPos = 'taux_' + titles[pi][:3] + '_pos.bin'
		outNeg = 'taux_' + titles[pi][:3] + '_Neg.bin'
		
		PRINT = 0
		if PRINT:
			print('Printing tauPos')
			for val in tauPos:
				print(str(val)+', ', end='')
			print()	
			print('Printing tauNeg')
			for val in tauNeg:
				print(str(val)+', ', end='')	
		
		outPosFile = np.zeros((240,200), dtype='float32')
		for i in range(240):
			outPosFile[i] = np.float32(tauPos)
		with open(windPath+outPos, 'wb') as	file_:
			file_.write(outPosFile)
	
		#==

		
		# Plot
		lat = latMCS
		plt.figure(figsize=(21,6))
		plt.subplot(131)
		plt.title(titles[pi] + ': wind speed')
		plt.plot(uPos, lat, label='pos ' + titles[pi], color='b')
		plt.plot(uNeg, lat, label='neg ' + titles[pi], color='orange')
		plt.plot(u_interp, lat, linestyle='dashed', color='b')
		plt.plot(-u_interp, lat, linestyle='dashed', color='orange')
		plt.plot(uRef, lat, label='idealised reference', color='k')
		plt.grid(); plt.legend()
		plt.subplot(132)
		plt.title(titles[pi] + ': wind stress')
		plt.plot(tauPos, lat, label='pos ' + titles[pi], color='b')
		plt.plot(tauNeg, lat, label='neg ' + titles[pi], color='orange')
		plt.plot(tauRefMCS, lat, label='idealised reference', color='k')
		plt.grid(); plt.legend()
		plt.subplot(133)
		plt.title(titles[pi] + ': wind stress curl')
		plt.plot(tools.ddx(tauPos, 2.5e3), lat, label='pos ' + titles[pi], color='b')
		plt.plot(tools.ddx(tauNeg, 2.5e3), lat, label='neg ' + titles[pi], color='orange')
		plt.plot(tools.ddx(tauRefMCS, 2.5e3), lat, label='idealised reference', color='k')
		plt.grid(); plt.legend()
		plt.show()
#==

if False:

	nyears = 100
	
	for pi in [1,4]:
	
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
		
		print(u.shape)
		
		u *= nyears
		v *= nyears

		#==
		
		# Compute wind stresses and wind speeds
		
		tauRef = 0.025 * tools.get05cosWind(u.shape[1], u.shape[0])
		uRef = np.abs(tauRef/(rhoa*ca))**0.5 * np.sign(tauRef) 

		uPos = uRef + u; uNeg = uRef - u
		
		tauPos = ca * rhoa * (uPos**2 + v**2)**0.5 * uPos
		tauNeg = ca * rhoa * (uNeg**2 + v**2)**0.5 * uNeg
		
		uabs = (u**2+v**2)**0.5
		tau = ca * rhoa * uabs * u

		# Get zonal means of everything.
		tauRef = np.mean(tauRef, axis=1)
		tauPos = np.mean(tauPos, axis=1)
		tauNeg = np.mean(tauNeg, axis=1)
		uav = np.mean(u, axis=1)
		uRef = np.mean(uRef, axis=1)
		uPos = np.mean(uPos, axis=1)
		uNeg = np.mean(uNeg, axis=1)
	
		INTERP = 1
		if INTERP:
			cubic = interp1d(lat, uav, kind='cubic')
			uMCS = cubic(latMCS)
			
		#==
		
		# Plot

		plt.figure(figsize=(15,6))
		plt.subplot(121)
		plt.title(titles[pi] + ': wind speed')
		plt.plot(uPos, lat, label='pos ' + titles[pi], color='r')
		plt.plot(uNeg, lat, label='neg ' + titles[pi], color='b')
		plt.plot(uav, lat, linestyle='dashed', color='r')
		plt.plot(-uav, lat, linestyle='dashed', color='b')
		plt.plot(uRef, lat, label='idealised reference', color='k')
		plt.grid(); plt.legend()
		plt.subplot(122)
		plt.title(titles[pi] + ': wind stress')
		plt.plot(tauPos, lat, label='pos ' + titles[pi], color='r')
		plt.plot(tauNeg, lat, label='neg ' + titles[pi], color='b')
		plt.plot(tauRef, lat, label='idealised reference', color='k')
		plt.grid(); plt.legend()
		plt.show()

#==

if False:
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
		
	
	
