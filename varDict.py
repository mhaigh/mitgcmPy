# varDict.py

# Contains a dictionary for each outputted mitgcm variable.
# Change with caution.  

#==========================================================

# Template:
# {'ncfile':'', 'var':'', 'vmin':, 'vmax':, cmap:''}


Rho = {'fname':'stateRho', 'var':'RHOAnoma', 'vmin':-2, 'vmax':-1, 'cmap':''}
Theta = {'fname':'stateTheta', 'var':'THETA', 'vmin':-2.5, 'vmax':2.5, 'cmap':'coolwarm'}

varDict={'Rho':Rho, 'Theta':Theta}

#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5; cmap='coolwarm'
#fname = 'stateUvel.nc'; var = 'UVEL'; cmap = 'coolwarm'; vmax = 0.1; vmin = -vmax
#fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2

