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

Vvel = {'fname':'stateVvel', 'VAR':'VVEL', 'vmin':-0.02, 'vmax':0.02, 'cmap':'coolwarm', 'title':'V (m/s)'}

Wvel = {'fname':'stateWvel', 'VAR':'WVEL', 'vmin':-1.e-5, 'vmax':1.e-5, 'cmap':'coolwarm', 'title':'W (m/s)'}

varDict={'ETAN':Eta, 'RHOAnoma':Rho, 'THETA':Theta, 'SALT':Salt, 'UVEL':Uvel, 'VVEL':Vvel, 'WVEL':Wvel}



varDict['DFrE_TH'] = {'fname':'stateTheta', 'vmin':-0.12, 'vmax':0.12, 'cmap':'coolwarm', 'title':'Pot. Temp. vertical flux'}

varDict['TOTTTEND'] = {'fname':'stateTheta', 'vmin':-1.e-3, 'vmax':1.e-3, 'cmap':'coolwarm', 'title':'Pot. Temp. Tendency'}

varDict['UVELTH'] = {'fname':'stateTheta', 'vmin':-1.e-1, 'vmax':1.e-1, 'cmap':'coolwarm', 'title':'Pot. Temp. U adv'}

varDict['VVELTH'] = {'fname':'stateTheta', 'vmin':-1.e-2, 'vmax':1.e-2, 'cmap':'coolwarm', 'title':'Pot. Temp. V adv'}

varDict['WVELTH'] = {'fname':'stateTheta', 'vmin':-1.e-5, 'vmax':1.e-5, 'cmap':'coolwarm', 'title':'Pot. Temp. W adv'}

varDict['PHIHYD'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Hydr. pressure pot. anom.'}

varDict['PHIBOT'] = {'fname':'state2D', 'vmin':-1, 'vmax':1, 'cmap':'coolwarm', 'title':'Bottom hydr. pressure pot. anom.'}

varDict['oceSflux'] = {'fname':'state2D', 'vmin':-0.0005, 'vmax':0.0005, 'cmap':'coolwarm', 'title':'Surface salt flux'}

varDict['botTauX'] = {'fname':'state2D', 'vmin':-0.05, 'vmax':0.05, 'cmap':'coolwarm', 'title':'Zonal bottom stress'}

varDict['botTauY'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. bottom stress'}

varDict['oceTAUX'] = {'fname':'state2D', 'vmin':-0.004, 'vmax':0.004, 'cmap':'coolwarm', 'title':'Zonal surface stress'}

varDict['oceTAUY'] = {'fname':'state2D', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. surface stress'}

varDict['oceQnet'] = {'fname':'state2D', 'vmin':0, 'vmax':1, 'cmap':'jet', 'title':' Surf. specific humidity (kg/kg)'}

varDict['SIarea'] = {'fname':'state2D', 'vmin':0.0, 'vmax':1., 'cmap':'viridis', 'title':'Sea-ice coverage'}

varDict['SIuice'] = {'fname':'state2D', 'vmin':-.1, 'vmax':.1, 'cmap':'viridis', 'title':'Sea-ice zonal vel.'}

varDict['SIvice'] = {'fname':'state2D', 'vmin':-.1, 'vmax':.1, 'cmap':'viridis', 'title':'Sea-ice merid. vel.'}

varDict['SIheff'] = {'fname':'state2D', 'vmin':0, 'vmax':1, 'cmap':'viridis', 'title':'Sea-ice effective ice thickness'}

varDict['SIhsnow'] = {'fname':'state2D', 'vmin':0, 'vmax':1, 'cmap':'viridis', 'title':'Sea-ice effective snow thickness'}

varDict['oceFWflx'] = {'fname':'state2D', 'vmin':-1.e-4, 'vmax':1.e-4, 'cmap':'coolwarm', 'title':'Freshwater flux'}

varDict['SHIfwFlx'] = {'fname':'state2D', 'vmin':0., 'vmax':1., 'cmap':'coolwarm', 'title':'Ice shelf freshwater flux'}

varDict['Um_dPhiX'] = {'fname':'stateUdpdx', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Zonal PGF'}

varDict['Vm_dPhiY'] = {'fname':'stateVdpdy', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. PGF'}

varDict['ADTAUU'] = {'fname':'adxx_tauu', 'vmin':-0.3, 'vmax':0.3, 'cmap':'coolwarm', 'title':'Merid. PGF'}



varDict['EXFuwind'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':'Zonal 10 m wind'}

varDict['EXFvwind'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':'Merid. 10 m wind'}

varDict['EXFpress'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':''}

varDict['EXFatemp'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':''}

varDict['EXFpreci'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':''}

varDict['EXFevap'] = {'fname':'state2D', 'vmin':-1.e-5, 'vmax':1.e-5, 'cmap':'coolwarm', 'title':''}

varDict['EXFhs'] = {'fname':'state2D', 'vmin':-1, 'vmax':1, 'cmap':'coolwarm', 'title':''}

varDict['EXFhl'] = {'fname':'state2D', 'vmin':-1, 'vmax':1, 'cmap':'coolwarm', 'title':''}

varDict['EXFqnet'] = {'fname':'state2D', 'vmin':-1, 'vmax':1, 'cmap':'coolwarm', 'title':''}

varDict['EXFroff'] = {'fname':'stateExf', 'vmin':-1.e-4, 'vmax':0, 'cmap':'coolwarm', 'title':''}

varDict['EXFswdn'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':''}

varDict['EXFlwdn'] = {'fname':'state2D', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':''}

varDict['EXFaqh'] = {'fname':'stateExf', 'vmin':-4., 'vmax':4., 'cmap':'coolwarm', 'title':''}

	

#==





titleData = {}

titleData['MCS_313'] = 'IPO pos PH'

titleData['MCS_314'] = 'IPO neg PH'



#==



def getPlottingVars(var):

	

	d = varDict[var]

	vmin = d['vmin']; vmax = d['vmax']

	cmap = d['cmap']

	title = d['title']



	return vmin, vmax, cmap, title



#==



def getTitleData(var):



	return titleData[var]



#==



def getTrefSref():



	return varDict['THETA']['TrefN'], varDict['SALT']['SrefN']



#==	





#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1

#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'

#fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax

#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2

