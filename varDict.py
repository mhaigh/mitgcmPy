# varDict.py

# Contains a dictionary for each outputted mitgcm variable.

#==========================================================

# Template:
# {'ncfile':'', 'var':'', 'vmin':, 'vmax':, cmap:''}

Eta = {'fname':'state2D', 'var':'ETAN', 'vmin':-10, 'vmax':10, 'cmap':'', 'title':'SSH'}
Rho = {'fname':'stateRho', 'var':'RHOAnoma', 'vmin':-2, 'vmax':-1, 'cmap':''}
Theta = {'fname':'stateTheta', 'var':'THETA', 'vmin':-1.8, 'vmax':1, 'cmap':'coolwarm', 'title':'THETA (deg. C)'}
Uvel = {'fname':'stateUvel', 'var':'UVEL', 'vmin':-0.02, 'vmax':0.02, 'cmap':'coolwarm', 'title':'U (m/s)'}
Vvel = {'fname':'stateUvel', 'var':'VVEL', 'vmin':-0.02, 'vmax':0.02, 'cmap':'coolwarm', 'title':'V (m/s)'}
Wvel = {'fname':'stateWvel', 'var':'WVEL', 'vmin':-1.e-6, 'vmax':1.e-6, 'cmap':'coolwarm', 'title':'W (m/s)'}

varDict={'2D':Eta, 'Rho':Rho, 'Theta':Theta, 'Uvel':Uvel, 'Vvel':Vvel, 'Wvel':Wvel}

def getPlottingVars(var):
	
	d = varDict[var]
	vmin = d['vmin']; vmax = d['vmax']
	cmap = d['cmap']
	title = d['title']

	return vmin, vmax, cmap, title

#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
#fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2

