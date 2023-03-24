# readSTDOUT.py

import os
import re
import matplotlib.pyplot as plt
import numpy as np

#==

f = open('STDOUT.0000', 'r')
data = []

for line in f:
	if re.search('ad_dynstat_aduvel_sd', line):
		data.append(float(line[-20:-1]))

nt = len(data)
dt = 7200./86400.
time = np.linspace(0,nt-1,nt)*dt
plt.plot(time, data)
plt.grid()
plt.show()


