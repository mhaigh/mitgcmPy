import numpy as np 
import matplotlib.pyplot as plt

#===========================================

N = 10000

x = np.zeros(N)
y = np.zeros(N)
pi = np.zeros(N)

inCircleX = []
inCircleY = []
inCircleCount = 0

centreX = 0.5
centreY = 0.5

#==

for i in range(N):

	x[i] = np.random.rand()
	y[i] = np.random.rand()
	
	if (x[i]-centreX)**2 + (y[i]-centreY)**2 < 0.25:
	
		inCircleX.append(x[i])
		inCircleY.append(y[i])
		inCircleCount = inCircleCount + 1				
		# If ends here
		
	pi[i] = 4 * inCircleCount / (i+1)

#==

plt.plot(pi)
plt.axhline(y=np.pi, linestyle='--')
plt.show()

#plt.plot(inCircleX, inCircleY, 'o')
#plt.show()
