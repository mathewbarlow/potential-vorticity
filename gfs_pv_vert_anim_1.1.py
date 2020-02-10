#
# run on python 3.7
#
# Python code to calculate 3D tropopause from online GFS analysis data
# and animate over the elevation angle. As the data is accessed
# online and the calculation is poorly coded, the program can take a while to 
# run.  GFS output contains its own tropopause data but it often has problem
# points.
#
# The date and lat-lon range can be set below.  Not currently set up to
# consider cases where the longitude switches from negative to positive,
# although that would require only minor changes.
#
# (poorly) coded by Mathew Barlow
# initial release: 14 Nov 2017
# last updated: 9 Feb 2020
#
# This code has *not* been extensively tested and has been 
# awkwardly translated from other coding languages, so if you find
# any errors or have any suggestions or improvements, including for
# the plotting, please let me know at Mathew_Barlow@uml.edu . Thanks!
#
# Support from NSF AGS-1623912 is gratefully acknowledged
#

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from mpl_toolkits.mplot3d import axes3d
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import imageio
import os

from datetime import datetime


# VALUES TO SET *************************************************
# set date, lat-lon range, and PV-value definition of tropopause
mydate='20200206'
myhour='06'
(lat1,lat2)=(20,60)
(lon1,lon2)=(-130,-60)
tpdef=2   # definition of tropopause in PVU
#****************************************************************
if lon2 < 0:
    lon1=lon1+360
    lon2=lon2+360

#constants
re=6.37e6
g=9.81
cp=1004.5
r=2*cp/7
kap=r/cp
omega=7.292e-5
pi=3.14159265

# open dataset, retreive variables, close dataset

url='https://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs'+\
mydate+'/gfs_0p25_'+myhour+'z_anl'

file = netCDF4.Dataset(url)

lat_in  = file.variables['lat'][:]
lon_in  = file.variables['lon'][:]
lev = file.variables['lev'][:]

pres2pv_in = file.variables['pres2pv'][0,:,:]
pressfc_in = file.variables['pressfc'][0,:,:]

nlev = lev.size
nx = lon_in.size
ny = lat_in.size

u_in = np.full((nlev, ny, nx), None)
v_in = np.full((nlev, ny, nx), None)
t_in = np.full((nlev, ny, nx), None)
hgt_in = np.full((nlev, ny, nx), None)

#
# the following lines shouldn't really be necessary but reduce
# data access issues, at least for me
#

ilev = 0
while ilev < nlev:
    print(ilev)
    u_in[ilev, :, :] = file.variables['ugrdprs'][0, ilev, :, :]
    ilev = ilev + 1
    
ilev = 0
while ilev < nlev:
    v_in[ilev, :, :] = file.variables['vgrdprs'][0, ilev, :, :]
    ilev = ilev + 1
    
ilev = 0
while ilev < nlev:
    t_in[ilev, :, :] = file.variables['tmpprs'][0, ilev, :, :]
    ilev = ilev + 1

ilev = 0
while ilev < nlev:
    hgt_in[ilev, :, :] = file.variables['hgtprs'][0, ilev, :, :]
    ilev = ilev + 1

# if data access were working smoothly, these would be simpler:    
#t_in = file.variables['tmpprs'][0,:,:,:]
#u_in = file.variables['ugrdprs'][0,:,:,:]
#v_in = file.variables['vgrdprs'][0,:,:,:]
#hgt_in = file.variables['hgtprs'][0,:,:,:]

file.close()

# get array indices for lat-lon range
# specified above
iy1 = np.argmin( np.abs( lat_in - lat1 ) )
iy2 = np.argmin( np.abs( lat_in - lat2 ) ) 
ix1 = np.argmin( np.abs( lon_in - lon1 ) )
ix2 = np.argmin( np.abs( lon_in - lon2 ) )  

# select specified lat-lon range
t=t_in[:,iy1:iy2,ix1:ix2]
lon=lon_in[ix1:ix2]
lat=lat_in[iy1:iy2]
u=u_in[:,iy1:iy2,ix1:ix2]
v=v_in[:,iy1:iy2,ix1:ix2]
hgt=hgt_in[:,iy1:iy2,ix1:ix2]
pres2pv=pres2pv_in[iy1:iy2,ix1:ix2]
pressfc=pressfc_in[iy1:iy2,ix1:ix2]

# some prep work for derivatives
xlon,ylat=np.meshgrid(lon,lat)

# define potential temperature and Coriolis parameter
theta=t*(1.E5/(lev[:,np.newaxis,np.newaxis]*100))**kap
f=2*omega*np.sin(ylat*pi/180)

lon = np.array(lon, dtype='float')
lat = np.array(lat, dtype='float')
lev = np.array(lev, dtype='float')
u = np.array(u, dtype='float')
v = np.array(v, dtype='float')
hgt = np.array(hgt, dtype='float')
pres2pv = np.array(pres2pv, dtype='float')
pressfc = np.array(pressfc, dtype='float')
theta = np.array(theta, dtype='float')
f = np.array(f, dtype='float')

# calculate derivatives

def ddp(f):
# handle unevenly-spaced levels with 2nd order
# Lagrange interpolation
# except for top and bottom, where use forward diff
    lev3=lev.reshape(lev.size,1,1)*100
    dpp=lev3-np.roll(lev3,-1,axis=0)
    dpm=lev3-np.roll(lev3,1,axis=0)
    fp=np.roll(f,-1,axis=0)
    fm=np.roll(f,1,axis=0)
    ddp_f=(
        fm*dpp/( (dpp-dpm)*(-dpm) ) +
        f*(dpp+dpm)/( dpm*dpp ) +
        fp*dpm/( (dpm-dpp)*(-dpp) )
        )
    ddp_f[0,:,:]=(f[1,:,:]-f[0,:,:])/(lev3[1,:,:]-lev3[0,:,:])
    ddp_f[-1,:,:]=(f[-1,:,:]-f[-2,:,:])/(lev3[-2,:,:]-lev3[-1,:,:])
    return(ddp_f)
    
def ddx(f):
# use center-difference, assuming evenly spaced lon
# except for side-boundaries, where use forward diff
    x=(re*np.cos(ylat*np.pi/180)*np.pi/180)*lon
    x3=x.reshape(1,x.shape[0],x.shape[1])
    dx3=np.roll(x3,-1,axis=2)-np.roll(x3,1,axis=2)
    ddx_f=(np.roll(f,-1,axis=2)-np.roll(f,1,axis=2))/dx3
    ddx_f[:,:,0]=(f[:,:,1]-f[:,:,0])/(x3[:,:,1]-x3[:,:,0])
    ddx_f[:,:,-1]=(f[:,:,-2]-f[:,:,-1])/(x3[:,:,-2]-x3[:,:,-1])
    return(ddx_f)
    
def ddy(f):
# use center-difference, assuming evenly spaced lon
# except for N/S boundaries, where use forward diff
    y=(re*np.pi/180)*lat
    y3=y.reshape(1,y.shape[0],1)
    dy3=np.roll(y3,-1,axis=1)-np.roll(y3,1,axis=1)
    ddy_f=(np.roll(f,-1,axis=1)-np.roll(f,1,axis=1))/dy3
    ddy_f[:,0,:]=(f[:,1,:]-f[:,0,:])/(y3[:,1,:]-y3[:,0,:])
    ddy_f[:,-1,:]=(f[:,-2,:]-f[:,-1,:])/(y3[:,-2,:]-y3[:,-1,:])
    return(ddy_f)
    
# could also use numpy built-in gradient function but
# I abandoned this for some reason
#lev3=lev.reshape(lev.size,1,1)
#ddp_theta=np.gradient(theta,lev3*100,axis=0)
#ddx_theta=np.gradient(theta,axis=2)/dx
#ddy_theta=np.gradient(theta,axis=1)/dy

# smooth data a bit
gf=1

ddp_theta=ddp(theta)
ddp_u=ddp(gaussian_filter(u,sigma=gf))
ddp_v=ddp(gaussian_filter(v,sigma=gf))

ddx_theta=ddx(theta)
ddy_theta=ddy(theta)
ddx_v=ddx(gaussian_filter(v,sigma=gf))
ddy_ucos=ddy(gaussian_filter(u,sigma=gf)*np.cos(ylat*pi/180))

# calculate contributions to PV and PV
absvort=ddx_v-(1/np.cos(ylat*pi/180))*ddy_ucos+f
pv_one=g*absvort*(-ddp_theta)
pv_two=g*(ddp_v*ddx_theta-ddp_u*ddy_theta)
pv=pv_one+pv_two

# calculate pressure of tropopause, Fortran-style (alas!)
# as well as potential temperature (theta) and height
#
# starting from 10hPa and working down, to avoid
# more complicated vertical structure higher up
#
nx=ix2-ix1+1
ny=iy2-iy1+1
nz=lev.size
nzs=np.argwhere(lev==50.0)[0,0]
tp=np.empty((ny-1,nx-1))*np.nan   # initialize as undef
tp_theta=np.empty((ny-1,nx-1))*np.nan   # initialize as undef
tp_hgt=np.empty((ny-1,nx-1))*np.nan   # initialize as undef

for ix in range(0,nx-1):
    for iy in range(0,ny-1):
        for iz in range(nzs,0,-1):
            if pv[iz,iy,ix]/1e-6<=tpdef:
                if np.isnan(tp[iy,ix]):
                    tp[iy,ix]=(
                    (lev[iz]*(pv[iz+1,iy,ix]-tpdef*1e-6)
                    -lev[iz+1]*(pv[iz,iy,ix]-tpdef*1e-6))/
                    (pv[iz+1,iy,ix]-pv[iz,iy,ix])
                    )
    
                    tp_theta[iy,ix]=(
                    ((lev[iz]-tp[iy,ix])*theta[iz+1,iy,ix]+
                    (tp[iy,ix]-lev[iz+1])*theta[iz,iy,ix])/
                    (lev[iz]-lev[iz+1])
                    )
                    
                    tp_hgt[iy,ix]=(
                    ((lev[iz]-tp[iy,ix])*hgt[iz+1,iy,ix]+
                    (tp[iy,ix]-lev[iz+1])*hgt[iz,iy,ix])/
                    (lev[iz]-lev[iz+1])
                    )

# calculate PV on the 330K isentropic surface
# (also not in a pythonic way)
nx=ix2-ix1+1
ny=iy2-iy1+1
nz=lev.size
pv330=np.empty((ny-1,nx-1))*np.nan   # initialize as undef
for ix in range(0,nx-1):
    for iy in range(0,ny-1):
        for iz in range(nz-2,0,-1):
            if theta[iz,iy,ix]>=330:
                if theta[iz-1,iy,ix]<=330:
                    if np.isnan(pv330[iy,ix]):
                        pv330[iy,ix]=(
                        ((330-theta[iz-1,iy,ix])*pv[iz,iy,ix]+
                        (theta[iz,iy,ix]-330)*pv[iz-1,iy,ix])/
                        (theta[iz,iy,ix]-theta[iz-1,iy,ix])
                        )
                   

# slight smoothing of result
tp=gaussian_filter(tp,sigma=1)
tp_theta=gaussian_filter(tp_theta,sigma=1)
pv330=gaussian_filter(pv330,sigma=1)


# get date for plotting
fdate=datetime.strptime(mydate, '%Y%m%d').strftime('%d %b %Y')

# output location for files
dirout= '/Users/mathew_barlow/downloads/'

plt.close(fig='all')

print('got here')

images = []
filenames = []

# make individual still images at range of viewpoint elevations
nframe=30
iframe=0
while iframe<=nframe:
    plt.figure(iframe,figsize=plt.figaspect(0.5))

    pressfc_smooth=gaussian_filter(pressfc,sigma=1)
    ax=plt.gca(projection='3d')

    surf=ax.plot_surface(xlon,ylat,tp,cmap="coolwarm",alpha=1,
                       rstride=1,cstride=1,
                       linewidth=0, antialiased=False)

    ax.plot_surface(xlon,ylat,pressfc_smooth/100,color="lightgray",
                       rstride=1,cstride=1,
                       linewidth=0, antialiased=False)
                       
    ax.set_zlim(1000,100)
    ax.set_xlim(lon1,lon2)
    ax.set_ylim(lat1,lat2)
    ax.view_init(elev=90 - iframe*90/nframe,azim=-90)
    
    stime = myhour+'Z '+fdate

    plt.title('2PVU Dynamic Tropopause over topography\n'+stime)
    
    plt.colorbar(surf)
                    
    pdir = dirout
    filename=pdir +'temp'+ '{:04d}'.format(iframe)+'.png'
    plt.savefig(filename, bbox_inches='tight')
    
    images.append(imageio.imread(filename))
    filenames.append(filename)
    
    plt.clf()
                    
    iframe=iframe+1

# combine individual images into an animated gif    
imageio.mimsave(pdir+'tp_anom' + stime + '.gif', images)

# remove all individual images
for file in filenames:
    os.remove(file)

plt.close('all')
