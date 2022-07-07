import sys

import numpy as np

from grid import Grid
from readData import readVariable

#==========================================================

EXP = sys.argv[1]
VAR = sys.argv[2]

#path = '/data/oceans_output/shelf/michai/mitgcm/'+EXP+'/run/'
path = '/home/michael/Documents/data/'+EXP+'/run/'

depths = [-10,-450]
levels = []

ts = 0

#==

grid = Grid(path)
for depth in depths:
	levels.append(grid.getIndexFromDepth(depth))
		
data = readVariable(VAR, path, meta=False)[ts:,levels,]
data = np.ma.filled(data, fill_value=0)

outname = VAR + '_' + EXP + '_lvls'
for level in levels:
	outname += '_' + str(level)
	
np.save(outname, data)
	
#data1 = readVariable(VAR, path, meta=False)[ts:,]
#import plotting as pt
#pt.plot1by2([data[-1,0], data1[-1,levels[0]]])
#pt.plot1by2([data[-1,1], data1[-1,levels[1]]])


