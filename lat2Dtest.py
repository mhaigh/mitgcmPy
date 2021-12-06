# test.py

# Test that ensures xarray and dependencies are functioning.

#==========================================================

import sys
import os 

import numpy as np
#import pandas as pd
#import xarray as xr
import netCDF4 as nc

from grid import Grid, choose_grid

import matplotlib.pyplot as plt
import matplotlib.colors as cl

#==========================================================

def polar_stereo (lon, lat, a=6378137., e=0.08181919, lat_c=-71, lon0=0):

    # Deep copies of arrays in case they are reused
    lon = np.copy(lon)
    lat = np.copy(lat)

    if lat_c < 0:
        # Southern hemisphere
        pm = -1
    else:
        # Northern hemisphere
        pm = 1

    # Prepare input
    lon = lon*pm*deg2rad
    lat = lat*pm*deg2rad
    lat_c = lat_c*pm*deg2rad
    lon0 = lon0*pm*deg2rad

    # Calculations
    t = np.tan(np.pi/4 - lat/2)/((1 - e*np.sin(lat))/(1 + e*np.sin(lat)))**(e/2)
    t_c = np.tan(np.pi/4 - lat_c/2)/((1 - e*np.sin(lat_c))/(1 + e*np.sin(lat_c)))**(e/2)
    m_c = np.cos(lat_c)/np.sqrt(1 - (e*np.sin(lat_c))**2)
    rho = a*m_c*t/t_c
    x = pm*rho*np.sin(lon - lon0)
    y = -pm*rho*np.cos(lon - lon0)

    return x, y
    
def get_x_y (lon, lat, pster=False):
    if pster:
        x, y = polar_stereo(lon, lat)
    else:
        x = lon
        y = lat
    return x, y


def cell_boundaries (data, grid, gtype='t', extrapolate=True, pster=False):

    # Inner function to pad the given array in the given direction(s), either extrapolating or copying.
    def extend_array (A, south=None, north=None, west=None, east=None):
        if south is not None:
            if south == 'extrapolate':
                s_bdry = 2*A[0,:]-A[1,:]
            elif south == 'copy':
                s_bdry = A[0,:]
            A = np.concatenate((s_bdry[None,:], A), axis=0)
        if north is not None:
            if north == 'extrapolate':
                n_bdry = 2*A[-1,:]-A[-2,:]
            elif north == 'copy':
                n_bdry = A[-1,:]
            A = np.concatenate((A, n_bdry[None,:]), axis=0)
        if west is not None:
            if west == 'extrapolate':
                w_bdry = 2*A[:,0]-A[:,1]
            elif west == 'copy':
                w_bdry = A[:,0]
            A = np.concatenate((w_bdry[:,None], A), axis=1)
        if east is not None:
            if east == 'extrapolate':
                e_bdry = 2*A[:,-1]-A[:,-2]
            elif east == 'copy':
                e_bdry = A[:,-1]
            A = np.concatenate((A, e_bdry[:,None]), axis=1)
        return A

    if gtype in ['t', 'w']:
        # Tracer grid: at centres of cells
        # Boundaries are corners of cells
        lon = grid.lon_corners_2d
        lat = grid.lat_corners_2d
        # Care about eastern and northern edges
        if extrapolate:
            lon = extend_array(lon, north='copy', east='extrapolate')
            lat = extend_array(lat, north='extrapolate', east='copy')
        else:
            data = data[...,:-1,:-1]
    elif gtype == 'u':
        # U-grid: on left edges of cells
        # Boundaries are centres of cells in X, corners of cells in Y
        lon = grid.lon_2d
        lat = grid.lat_corners_2d
        # Care about western and northern edges
        if extrapolate:
            lon = extend_array(lon, north='copy', west='extrapolate')
            lat = extend_array(lat, north='extrapolate', west='copy')
        else:
            data = data[...,:-1,1:]
    elif gtype == 'v':
        # V-grid: on bottom edges of cells
        # Boundaries are corners of cells in X, centres of cells in Y
        lon = grid.lon_corners_2d
        lat = grid.lat_2d
        # Care about eastern and southern edges
        if extrapolate:
            lon = extend_array(lon, south='copy', east='extrapolate')
            lat = extend_array(lat, south='extrapolate', east='copy')
        else:
            data = data[...,1:,:-1]
    elif gtype == 'psi':
        # Psi-grid: on southwest corners of cells
        # Boundaries are centres of cells
        lon = grid.lon_2d
        lat = grid.lat_2d
        # Care about western and southern edges
        if extrapolate:
            lon = extend_array(lon, south='copy', west='extrapolate')
            lat = extend_array(lat, south='extrapolate', west='copy')
        else:
            data = data[...,1:,1:]

    # Convert to polar stereographic if needed
    x, y = get_x_y(lon, lat, pster=pster)
    return x, y, data
    
def shade_mask (ax, mask, grid, gtype='t', pster=False, colour='grey', rasterized=False):

    # Properly mask all the False values, so that only True values are unmasked
    mask_plot = np.ma.masked_where(np.invert(mask), mask)
    # Prepare quadrilateral patches
    x, y, mask_plot = pcell_boundaries(mask_plot, grid, gtype=gtype, pster=pster)
    if colour == 'grey':
        rgb = (0.6, 0.6, 0.6)
    elif colour == 'white':
        rgb = (1, 1, 1)
    else:
        print(('Error (shade_mask): invalid colour ' + colour))
        sys.exit()
    # Add to plot        
    img = ax.pcolormesh(x, y, mask_plot, cmap=cl.ListedColormap([rgb]), linewidth=0, rasterized=rasterized)
    img.set_edgecolor('face')

    
def shade_land (ax, grid, gtype='t', pster=False, land_mask=None, rasterized=False):
    if land_mask is None:
        land_mask = grid.get_land_mask(gtype=gtype)
    shade_mask(ax, land_mask, grid, gtype=gtype, pster=pster, rasterized=rasterized)

    
def shade_land_ice (ax, grid, gtype='t', pster=False, land_mask=None, ice_mask=None, rasterized=False):
    if land_mask is None:
        land_mask = grid.get_land_mask(gtype=gtype)
    if ice_mask is None:
        ice_mask = grid.get_ice_mask(gtype=gtype)
    shade_mask(ax, land_mask+ice_mask, grid, gtype=gtype, pster=pster, rasterized=rasterized)
    
#==========================================================
#==========================================================

# CODE START


# Load data
path_grid = '/Users/mh115/Documents/BAS/data/PAS_666/'
path = '/Users/mh115/Documents/BAS/data/PAS_666/'
vmin = None; vmax = None

#fname = 'stateRho.nc'; var = 'RHOAnoma'; vmin = -2; vmax = - 1
#fname = 'stateTheta.nc'; var = 'THETA'; vmin = - 2.5; vmax = 2.5
#fname = 'stateUvel.nc'; var = 'UVEL'; vmin = -0.2; vmax = 0.2
fname = 'stateVvel.nc'; var = 'VVEL'; vmin = -0.2; vmax = 0.2

level = 0

ncfile = nc.Dataset(path+fname, 'r')
data = ncfile.variables[var][-1,:]

#plt.plot(data[:,100,100])
#plt.show()

# Build grid
grid = Grid(path_grid)
gtype='t'; pster=False; land_mask=None; ice_mask=None;rasterized=False
x, y, data_plot = cell_boundaries(data, grid, gtype=gtype, pster=pster)


# Create plot
fig, ax = plt.subplots()
#shade_land_ice(ax, grid, gtype=gtype, pster=pster, land_mask=land_mask, ice_mask=ice_mask, rasterized=rasterized)
shade_land_ice(ax, grid, gtype=gtype, pster=pster, land_mask=land_mask, rasterized=rasterized)
img = ax.pcolormesh(x, y, data[level], rasterized=rasterized, vmin=vmin, vmax=vmax)

plt.gca().set_aspect('equal', adjustable='box')

plt.title(var + ', level ' + str(level))
fig.colorbar(img, shrink=0.3)
plt.tight_layout()
plt.savefig(var+"_level"+str(level)+".png",bbox_inches='tight',dpi=200)
#plt.show()


ncfile.close()


