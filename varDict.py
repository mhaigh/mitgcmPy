# varDict.py

# Contains a dictionary for each outputted mitgcm variable.

#==========================================================

# Template:
# {'ncfile':'', 'var':'', 'vmin':, 'vmax':, cmap:''}

TrefN=[-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.800e+00,-1.660e+00,-1.380e+00,-1.100e+00,-8.200e-01,-5.400e-01,-2.600e-01,2.000e-02,3.000e-01,5.800e-01,8.600e-01,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00,1.000e+00]
SrefN=[3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.350e+01,3.355e+01,3.365e+01,3.375e+01,3.385e+01,3.395e+01,3.405e+01,3.415e+01,3.425e+01,3.435e+01,3.445e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01,3.450e+01]



Eta = {'fname':'state2D', 'VAR':'ETAN', 'vmin':-0.2, 'vmax':0.2, 'cmap':'jet', 'title':'SSH anomaly (m)'}
Rho = {'fname':'stateRho', 'VAR':'RHOAnoma', 'vmin':-.2, 'vmax':.2, 'cmap':'jet', 'title':'Rho anom.'}
Theta = {'fname':'stateTheta', 'VAR':'THETA', 'vmin':-1.8, 'vmax':1, 'cmap':'coolwarm', 'title':'THETA (deg. C)','TrefN':TrefN}
Salt = {'fname':'stateSalt', 'VAR':'SALT', 'vmin':33.2, 'vmax':34.5, 'cmap':'jet', 'title':'Salinity (g/kg)', 'SrefN':SrefN}
Uvel = {'fname':'stateUvel', 'VAR':'UVEL', 'vmin':-0.2, 'vmax':0.2, 'cmap':'coolwarm', 'title':'U (m/s)'}
Vvel = {'fname':'stateVvel', 'VAR':'VVEL', 'vmin':-0.01, 'vmax':0.01, 'cmap':'coolwarm', 'title':'V (m/s)'}
Wvel = {'fname':'stateWvel', 'VAR':'WVEL', 'vmin':-1.e-5, 'vmax':1.e-5, 'cmap':'coolwarm', 'title':'W (m/s)'}

varDict={'ETAN':Eta, 'RHOAnoma':Rho, 'THETA':Theta, 'SALT':Salt, 'UVEL':Uvel, 'VVEL':Vvel, 'WVEL':Wvel}

varDict['DFrE_TH'] = {'fname':'stateTheta', 'vmin':-0.12, 'vmax':0.12, 'cmap':'coolwarm', 'title':'Pot. Temp. vertical flux'}
varDict['TOTTTEND'] = {'fname':'stateTheta', 'vmin':-1.e-3, 'vmax':1.e-3, 'cmap':'coolwarm', 'title':'Pot. Temp. Tendency'}
varDict['UVELTH'] = {'fname':'stateTheta', 'vmin':-1.e-1, 'vmax':1.e-1, 'cmap':'coolwarm', 'title':'Pot. Temp. U adv'}
varDict['VVELTH'] = {'fname':'stateTheta', 'vmin':-1.e-2, 'vmax':1.e-2, 'cmap':'coolwarm', 'title':'Pot. Temp. V adv'}
varDict['WVELTH'] = {'fname':'stateTheta', 'vmin':-1.e-5, 'vmax':1.e-5, 'cmap':'coolwarm', 'title':'Pot. Temp. W adv'}
varDict['PHIHYD'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Hydr. pressure pot. anom.'}
varDict['PHIBOT'] = {'fname':'state2D', 'vmin':-1, 'vmax':1, 'cmap':'coolwarm', 'title':'Bottom hydr. pressure pot. anom.'}
varDict['botTauX'] = {'fname':'state2D', 'vmin':-0.05, 'vmax':0.05, 'cmap':'coolwarm', 'title':'Zonal bottom stress'}
varDict['botTauY'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. bottom stress'}
varDict['oceTAUX'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Zonal surface stress'}
varDict['oceTAUY'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. surface stress'}
varDict['Um_dPhiX'] = {'fname':'stateUdpdx', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Zonal PGF'}
varDict['Vm_dPhiX'] = {'fname':'stateVdpdx', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. PGF'}

#==

def getPlottingVars(var):
	
	d = varDict[var]
	vmin = d['vmin']; vmax = d['vmax']
	cmap = d['cmap']
	title = d['title']

	return vmin, vmax, cmap, title

#==

def getTrefSref():

	return varDict['THETA']['TrefN'], varDict['SALT']['SrefN']

#==	

#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
#fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2

