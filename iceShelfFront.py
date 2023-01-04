import numpy as np

#==

# Build and return an open ocean mask for the given grid type.
def get_open_ocean_mask (grid):

    # Start with array of all ones
    open_ocean = np.ones([grid.Ny, grid.Nx])
    # Set to zero in land and ice shelf regions
    open_ocean[grid.iceC] = 0
    open_ocean[grid.landC] = 0

    return open_ocean


# Return the ice shelf mask for the given ice shelf and grid type.
def get_ice_mask (self, shelf='all', gtype='h'):

    # Select grid type
    if gtype == 'h':
        ice_mask_all = self.ice_mask
    elif gtype == 'u':
        ice_mask_all = sel
    elif gtype == 'v':
        ice_mask_all = self.ice_mask_v
    else:
        print(('Error (get_ice_mask): no mask exists for the ' + gtype + ' grid'))
        sys.exit()
    
    # Select ice shelf
    if shelf == 'all':
        return ice_mask_all
    else:
        return self.restrict_mask(ice_mask_all, shelf, gtype=gtype)

# Finds the value of the given array to the west, east, south, north of every point, as well as which neighbours are non-missing, and how many neighbours are non-missing.
# Can also do 1D arrays (so just neighbours to the left and right) if you pass use_1d=True.
def neighbours (data, missing_val=-9999, use_1d=False):

    # Find the value to the west, east, south, north of every point
    # Just copy the boundaries
    data_w = np.empty(data.shape)
    data_w[...,1:] = data[...,:-1]
    data_w[...,0] = data[...,0]
    data_e = np.empty(data.shape)
    data_e[...,:-1] = data[...,1:]
    data_e[...,-1] = data[...,-1]
    if not use_1d:
        data_s = np.empty(data.shape)
        data_s[...,1:,:] = data[...,:-1,:]
        data_s[...,0,:] = data[...,0,:]
        data_n = np.empty(data.shape)
        data_n[...,:-1,:] = data[...,1:,:]
        data_n[...,-1,:] = data[...,-1,:]     
    # Arrays of 1s and 0s indicating whether these neighbours are non-missing
    valid_w = (data_w != missing_val).astype(float)
    valid_e = (data_e != missing_val).astype(float)
    if use_1d:
        # Number of valid neighoburs of each point
        num_valid_neighbours = valid_w + valid_e
        # Finished
        return data_w, data_e, valid_w, valid_e, num_valid_neighbours
    valid_s = (data_s != missing_val).astype(float)
    valid_n = (data_n != missing_val).astype(float)
    num_valid_neighbours = valid_w + valid_e + valid_s + valid_n
    return data_w, data_e, data_s, data_n, valid_w, valid_e, valid_s, valid_n, num_valid_neighbours

#==

# Find all ice shelf front points and return them as a list.
# For a specific ice shelf, pass a special ice_mask 
def get_ice_shelf_front (grid, ice_mask, open_ocean, gtype='h', xmin=None, xmax=None, ymin=None, ymax=None, neighbourPts=1):

	iceMask = grid.iceC
	oce = get_open_ocean_mask(grid)

	# Set any remaining bounds
	lon = grid.XC
	lat = grid.YC
	if xmin is None:
		xmin = np.amin(lon)
	if xmax is None:
		xmax = np.amax(lon)
	if ymin is None:
		ymin = np.amin(lat)
	if ymax is None:
		ymax = np.amax(lat)

	fronta = np.nan * np.zeros((grid.Ny, grid.Nx))
	frontl = []
	for j in range(grid.Ny):
		for i in range(grid.Nx):
			if iceMask[j,i]:
				T = oce[j,i+1]+oce[j,i-1]+oce[j+1,i]+oce[j-1,i]
				if T > 0:
					frontl.append([lon[j,i],lat[j,i]])					
					fronta[j,i] = 1
					if neighbourPts == 1:
						fronta[j,i+1] = 1; fronta[j,i-1] = 1
						fronta[j+1,i] = 1; fronta[j-1,i] = 1
						fronta[j-1,i+1] = 1; fronta[j-1,i-1] = 1
						fronta[j+1,i+1] = 1; fronta[j+1,i-1] = 1
					elif neighbourPts == 2:
						fronta[j,i+1] = 1; fronta[j,i-1] = 1
						fronta[j+1,i] = 1; fronta[j-1,i] = 1

	# Prepare frontl so can it can be plotted using scatter.
	x, y = zip(*frontl)
	frontl = [x,y]

	return frontl, fronta

# Find all ice shelf front points and return them as a list.
# For a specific ice shelf, pass a special ice_mask 
def get_grounding_line(grid, ice_mask, open_ocean, gtype='h', xmin=None, xmax=None, ymin=None, ymax=None, neighbourPts=False):

	iceMask = grid.iceC
	draft = grid.draft
	hFac = grid.hFacC

	# Set any remaining bounds
	lon = grid.XC
	lat = grid.YC
	if xmin is None:
		xmin = np.amin(lon)
	if xmax is None:
		xmax = np.amax(lon)
	if ymin is None:
		ymin = np.amin(lat)
	if ymax is None:
		ymax = np.amax(lat)

	gla = np.nan * np.zeros((grid.Ny, grid.Nx))
	gll = []
	for j in range(grid.Ny):
		for i in range(grid.Nx):
			if draft[j,i] < 0:
				k = grid.getIndexFromDepth(draft[j,i])
				if hFac[k+2,j,i] <= 0.1:
					gll.append([lon[j,i],lat[j,i]])					
					gla[j,i] = 1
					if neighbourPts:
						gla[j,i+1] = 1; gla[j,i-1] = 1
						gla[j+1,i] = 1; gla[j-1,i] = 1
						gla[j-1,i+1] = 1; gla[j-1,i-1] = 1
						gla[j+1,i+1] = 1; gla[j+1,i-1] = 1

	# Prepare frontl so can it can be plotted using scatter.
	x, y = zip(*gll)
	gll = [x,y]

	return gll, gla


#==


