import numpy as np
import tools
from grid import Grid
import matplotlib.pyplot as plt

path = '/home/michael/Documents/data/MCS_155/run/'

grid = Grid(path)

Trel, Srel = tools.getTSrel(grid, salttoptop=33.5, salttop=33.5, saltbot=34.5, temptop=-1.8, tempbot=1.0, tclinewid=140., tclinetop=-200, hclinewid=140, hclinetop=-200, shift=0)
	
rho0 = 1030.
tAlpha=3.90e-5
sBeta =7.41e-4
g = 9.81
H = 200.
f0 = 1.4e-4
drho = 1.

rho =  rho0 * (1 - tAlpha * Trel + sBeta * Srel)

N = np.sqrt(g/rho0 * drho / H)
print(N)

Ld = N * H / f0

print(Ld)
 
 
#plt.plot(rho)
#plt.show()	


