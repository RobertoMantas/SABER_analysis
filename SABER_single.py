"""
Example to be done with one single file.
"""

import os
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import functools
import operator
from pylab import *
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.cm import get_cmap

years = [2002] #2003,2004,2005,2006,2007,2008,2009,2010,2011,2012
file_path =  "/Volumes/GAS_HDD/SABER_instrument/20022/001/SABER_L2A_2005001_16603_02.0.nc"#SABER_L2A_2002026_00736_02.0.nc"


sites_lat_y = []
sites_lon_x = []
altitude_list = []
ktemp_list = []
O3_96_list = []
ncfile = Dataset(file_path,'r')
latitude = ((ncfile.variables['tplatitude'])[:].tolist())
altitude = ((ncfile.variables['tpaltitude'])[:].tolist())
ktemp = ((ncfile.variables['ktemp'])[:].tolist())
longitude = ((ncfile.variables['tplongitude'])[:].tolist())
tpAD_list = ((ncfile.variables['tpAD'])[:].tolist())
O3_96 = ((ncfile.variables['O3_96'])[:].tolist())
sites_lat_y.extend(functools.reduce(operator.iconcat, latitude, []))
sites_lon_x.extend(functools.reduce(operator.iconcat, longitude, []))
altitude_list.extend(functools.reduce(operator.iconcat, altitude, []))
ktemp_list.extend(functools.reduce(operator.iconcat, ktemp, []))
O3_96_list.extend(functools.reduce(operator.iconcat, O3_96, []))


########################plot function as a heatmap########################
def plot_heatmap(colormaps,data):
    """
    Helper function to plot data with associated colormap.
    """
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            constrained_layout=True, squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap="viridis", rasterized=True, vmin=160, vmax=280)
        # psm = ax.pcolormesh(data, cmap="viridis", rasterized=True, vmin=0, vmax=0.000012)
        fig.colorbar(psm, ax=ax)
    plt.show()

cmap = ListedColormap(["darkorchid","plum", "violet","cornflowerblue", "mediumblue", "lime", "yellow","gold", "orange", "red"])
########################plot function####################################


#The following parameters are extracted from the papers. It explains how is the data arranged.
no = 1 #14 #number of orbits per day
nlat = 10;step_lat=10 #ten groups of 10º lat bins
nhei = 18; step_hei = 5 #18 groups of 5km bins altitude
nlon = 24 ; step_lon = 15#24 groups of 15º bins longitude

lat = np.arange(-45.,-45.+nlat*step_lat,step_lat,float)
height = np.arange(22.5,22.5+nhei*step_hei,step_hei,float)
lon = np.arange(7.5,7.5+nlon*step_lon,step_lon,float)
zonal_mean_temp = np.zeros((nhei,nlat), float)

####################### None values are replaced with -999 (missing value) #############################
latitude_list = [-999.0 if v is None else v for v in sites_lat_y]
altitude_list = [-999.0 if v is None else v for v in altitude_list]
ktemp_list = [-999.0 if v is None else v for v in ktemp_list]
# longitude_list = [-999.0 if v is None else v for v in longitude_list]
# O3_96_list = [-999.0 if v is None else v for v in O3_96_list]
####################### Not values are replaced with -999 (missing value) #############################

for z in range(nhei):
    for y in range(nlat):
        for j in range(nlon):
            dum = np.ma.masked_where(latitude_list < (lat[y]-(step_lat/2)), ktemp_list) #take only the values insdide the corresponding bin
            dum = np.ma.masked_where(latitude_list > (lat[y]+(step_lat/2)), dum)#take only the values insdide the corresponding bin
            dum = np.ma.masked_where(altitude_list < (height[z]-(step_hei/2)), dum)#take only the values insdide the corresponding bin
            dum = np.ma.masked_where(altitude_list > (height[z]+(step_hei/2)), dum)#take only the values insdide the corresponding bin
            # dum = np.ma.masked_where(longitude_list < (lon[j]-(step_lon/2)), dum)
            # dum = np.ma.masked_where(longitude_list > (lon[j]+(step_lon/2)), dum)
            dum = dum.compressed()
            dum = ([x for x in dum if x != -999.0]) #Remove all the -999 values before the mean
            zonal_mean_temp[z,y] = np.mean(dum)
            # temp[z,y,j] = np.mean(dum)
# This will call the function to plot the heatmap. The value being plotted is the temperature. 
# You will see a Matrix plotted witht he temperature data inside.
plot_heatmap([cmap],zonal_mean_temp[:,:]) 