# from __future__ import print_function
import os
# from scipy.io import netcdf
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import functools
import operator
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.cm import get_cmap

def average (x):
    return sum(x)/len(x)

# from Sientific.IO.NetCDF import 
years = [2002] #2003,2004,2005,2006,2007,2008,2009,2010,2011,2012
in_dir =  "/Volumes/GAS_HDD/SABER_instrument/20022/001"#SABER_L2A_2002026_00736_02.0.nc"




# print(nc.dimensions)
sites_lat_y = []
sites_lon_x = []
altitude_list = []
ktemp_list = []
O3_96_list = []

# def WRF_extract_list(in_dir, out_dir):
for base, dirs, files in os.walk(in_dir):
    for file in files:
        if str(file).startswith('SABER') and str(file).endswith('0.nc'):
            file_path = base+'/'+file
            # with open(out_dir+ '/OUTPUT_'+ str(file)+".dat" ,'w') as output_file:
            # nc = netcdf.netcdf_file(file_path,'r')
            ncfile = Dataset(file_path,'r')
            latitude = ((ncfile.variables['tplatitude'])[:].tolist())
            altitude = ((ncfile.variables['tpaltitude'])[:].tolist())
            ktemp = ((ncfile.variables['ktemp'])[:].tolist())
            longitude = ((ncfile.variables['tplongitude'])[:].tolist())
            tpAD_list = ((ncfile.variables['tpAD'])[:].tolist())
            O3_96 = ((ncfile.variables['O3_96'])[:].tolist())
            #With the next lines you are concatenating the data from the different orbits into one list.
            sites_lat_y.extend(functools.reduce(operator.iconcat, latitude, []))
            altitude_list.extend(functools.reduce(operator.iconcat, altitude, []))
            sites_lon_x.extend(functools.reduce(operator.iconcat, longitude, []))
            ktemp_list.extend(functools.reduce(operator.iconcat, ktemp, []))
            O3_96_list.extend(functools.reduce(operator.iconcat, O3_96, []))

#### To see the orbits of one day around the glob and get a sense of what data are you working with do the following:###
fig = plt.figure(figsize=(12,9))
m = Basemap(projection='moll',resolution='c',area_thresh=10000.,lat_0=30.,lon_0=0.)
# m = Basemap(projection="mill", llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180, urcrnrlon=180, resolution='c')
m.drawcoastlines()
m.drawparallels(np.arange(-90,90,10),labels=[True,False,False,False])
m.drawmeridians(np.arange(-180,180,15), labels=[0,0,0,1])
sites_lat_y1 = [x for x in sites_lat_y if x is not None]
sites_lon_x1 = [x for x in  sites_lon_x if x is not None]
m.scatter(sites_lon_x1, sites_lat_y1, latlon = True)
plt.show()

# As told in the paper by Pancheva, a first 3 standard deaviation has to be done to the data to do avoid 
# data out of bounds. Obtain the third st. dev. from the lower bound and the higher one. Anything out this range
# has to be removed.
#Here remove the None values of the lists
ktemp_clean = list(filter(None.__ne__, ktemp_list))
    # ktemp_clean.append(Not_none_values)
min_max_deviation_ktemp = ([(np.mean(ktemp_clean) - 3 * (np.std(ktemp_clean)), (np.mean(ktemp_clean)) + 3 * (np.std(ktemp_clean)))])
#Now you have the range in which the data should be. You will have to remove the data outside this range.
print(min_max_deviation_ktemp)
"""

HERE THE WORK DONE IN THE SCRIPT Saber_single.py has to be added. It is better to first have a full understanding of a 
single orbit by its self and then introduce it into all the data.


"""

########################plot function as a heatmap########################
def plot_heatmap(colormaps,data):
    """
    Helper function to plot data with associated colormap.
    """
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            constrained_layout=True, squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap="viridis", rasterized=True, vmin=min_max_deviation_ktemp[0][0], vmax=min_max_deviation_ktemp[0][1])
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


#TODO
"""
Here is where the work is left. In the last step where the function is implemented with the data. 
The following equations have not been checked. I extracted it from Pancheva et al., 2009 paper section 3.
You have to extract the amplitude and phase from the planetary waves. 
The reasoning was the following: if you think about the slope winds, you have upslope flow during daytime and 
downslope at night. You can approximate this behaviour to a sinusoidal function and therefore we fit the model-predicted 
meridional wind to a sinusoidal function of the form A * COS (omg*t - phi) where A is the amplitude, omg = 2*PI/24h is the frequency, 
and PHI is the phase, defined as the time of maximum amplitude in hours.
The functions are defined as following COS (a + b) = COS (a) COS (b) - SIN (a) SIN (b). This is the theory, the actual formula is the one
in section 3 in panchevas 2009 paper. 

You can look on the internet how to do cross-sections and plot it. 
An example to extract the cross-sections:

frequency = 2*PI/360
model = p0 * COS( FREQ * X ) + p1 * SIN( FREQ * X )
return, (y-model)/err

FREQ is omg in the equation above, so p0 corresponds to A*COS (phi) and p1 to A*SIN(phi).

We ignore err here and I give as X the time (in our case 24 hours with a step of 1 hour) and Y the the meridional wind at each time (this is done at every model grid-point). The fit returns p[0] and p[1] that minimize the least-squares. 

You then extract your amplitude and phase as

amp = SQRT(p0^2+p1^2) 
pha = ACOS(p0/(SQRT(p0^2+p1^2)))*180.0/PI

For the amplitude, SQRT(p0^2 + p1^2) = SQRT(A^2 * COS^2(phi) +A^2 * SIN^2(phi)) = SQRT(A^2) = A

For the phase you just do arc cos (A * COS(phi)/A) = arc cos (COS (phi)) = phi. However, you need to add the conditions below because of the nature of arc cos. Thank you!

In different papers and thesis I have seen that it can be solved with Fouriers transformations. Here you can find it explained:
http://www.diva-portal.org/smash/get/diva2:142105/FULLTEXT01.pdf
"""
amplitude_A = "?"#TODO
amplitude_B = "?"#TODO
amplitude_C = "?"#TODO
phaseA = "?"#TODO
phaseB = "?"#TODO
phaseC = "?"#TODO

T = [24,16,11,5]
#traveling planetary wave function
for j in range(4): # Planetary wave period
    for s in range(-3,4): # Wavenumber: <0 westward propagation, >0 eastward propagation
        amplitude_A * np.cos(((2*np.pi)/(24*T[j]))*hour - ((2*np.pi)/360)* s * np.deg2rad(longitude) - ((2*np.pi/360)*phaseA))
#stationary planetary wave function
for s in range(3): # Wavenumber: <0 westward propagation, >0 eastward propagation
    amplitude_B * np.cos((2*np.pi/360)*s* np.deg2rad(longitude) - (2*np.pi/360)*phaseB)

#24- 12-h tides function
for k in range(1, 3):# number of diurnal harmonic 1-2
    for s in range(-4, 5): # Wavenumber: <0 westward propagation, >0 eastward propagation
        amplitude_C * np.cos(k*(2*np.pi/24)*hour - (2*np.pi/360)* s * np.deg2rad(longitude) - (2*np.pi/360)*phaseC) 
        