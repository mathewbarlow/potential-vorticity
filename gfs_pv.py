#
# python code for some calculations related to the dynamic tropopause (DT)  -
# DT pressure, DT potential temperature, 330K PV, 
# and a cross-section of PV at the latitude where the tropopause is lowest -
# all based on the GFS analysis available online.  As the data is accessed
# online, the program can take a while to run.
#
# the date and lat-lon range can be set below
#
# (poorly) coded by Mathew Barlow
# last updated: 14 Nov 2017
#
# this code has *not* been extensively tested and has been 
# awkwardly translated from other coding languages, so if you find
# any errors or have any suggestions or improvements, including for
# the plotting, please let me know at Mathew_Barlow@uml.edu . Thanks!
#
# Support from NSF AGS-1623912 is gratefully acknowledged
#

import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
from scipy.ndimage import gaussian_filter
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime

# VALUES TO SET *************************************************
# set date, lat-lon range, and PV-value definition of tropopause
mydate='20171114'
myhour='06'
(lat1,lat2)=(25,60)
(lon1,lon2)=(-130,-60)
tpdef=2   # definition of tropopause in PVU
#****************************************************************

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

url='http://nomads.ncep.noaa.gov:80/dods/gfs_0p25/gfs'+\
mydate+'/gfs_0p25_'+myhour+'z_anl'

file = netCDF4.Dataset(url)

lat_in  = file.variables['lat'][:]
lon_in  = file.variables['lon'][:]
lev = file.variables['lev'][:]

pres2pv_in = file.variables['pres2pv'][0,:,:]

t_in = file.variables['tmpprs'][0,:,:,:]
u_in = file.variables['ugrdprs'][0,:,:,:]
v_in = file.variables['vgrdprs'][0,:,:,:]
hgt_in = file.variables['hgtprs'][0,:,:,:]

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

# some prep work for derivatives
xlon,ylat=np.meshgrid(lon,lat)
dlony,dlonx=np.gradient(xlon)
dlaty,dlatx=np.gradient(ylat)
dx=re*np.cos(ylat*pi/180)*dlonx*pi/180
dy=re*dlaty*pi/180

# define potential temperature and Coriolis parameter
theta=t*(1.E5/(lev[:,np.newaxis,np.newaxis]*100))**kap
f=2*omega*np.sin(ylat*pi/180)

# calculate derivatives
# (np.gradient can handle 1D uneven spacing,
# so build that in for p, but do dx and dy 
# external to the function since they are 2D)
ddp_theta=np.gradient(theta,lev*100,axis=0)
ddx_theta=np.gradient(theta,axis=2)/dx
ddy_theta=np.gradient(theta,axis=1)/dy
ddp_u=np.gradient(u,lev*100,axis=0)
ddp_v=np.gradient(v,lev*100,axis=0)
ddx_v=np.gradient(v,axis=2)/dx
ddy_ucos=np.gradient(u*np.cos(ylat*pi/180),axis=1)/dy

# calculate contributions to PV and PV
absvort=ddx_v-(1/np.cos(ylat*pi/180))*ddy_ucos+f
pv_one=g*absvort*(-ddp_theta)
pv_two=g*(ddp_v*ddx_theta-ddp_u*ddy_theta)
pv=pv_one+pv_two

# calculate pressure of tropopause, Fortran-style (alas!)
# as well as potential temperature (theta) and height
nx=ix2-ix1+1
ny=iy2-iy1+1
nz=lev.size
tp=np.empty((ny-1,nx-1))*np.nan   # initialize as undef
tp_theta=np.empty((ny-1,nx-1))*np.nan   # initialize as undef
tp_hgt=np.empty((ny-1,nx-1))*np.nan   # initialize as undef

for ix in range(0,nx-1):
    for iy in range(0,ny-1):
        for iz in range(nz-2,0,-1):
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
# (appears to work better than smoothing u,v,t first)
tp=gaussian_filter(tp,sigma=1)
tp_theta=gaussian_filter(tp_theta,sigma=1)
pv330=gaussian_filter(pv330,sigma=1)

# define spatial correlation function for testing results
def scorr(a,b):
    abar=np.mean(a)
    bbar=np.mean(b)
    covar=sum((a-abar)*(b-bbar))
    avar=sum((a-abar)**2)
    bvar=sum((b-bbar)**2)
    r=covar/np.sqrt(avar*bvar)
    return(r)

# identify latitude of lowest tropopause
maxloc=np.argwhere(tp==np.amax(tp))
latmax=lat[maxloc[0,0]]



# now make some plots - these badly need to be improved

states = NaturalEarthFeature(category='cultural', 
    scale='50m', facecolor='none', 
    name='admin_1_states_provinces_shp')

# get date for plotting
fdate=datetime.strptime(mydate, '%Y%m%d').strftime('%d %b %Y')

# plot of DT pressure
plt.figure(1)

ax = plt.axes(projection=ccrs.PlateCarree( ))
ax.set_extent([lon1,lon2,lat1,lat2],crs=ccrs.PlateCarree())
clevs=np.arange(50,800,50)
plt.contour(lon,lat,tp,clevs,transform=ccrs.PlateCarree(),colors='black',
linewidths=0.5)
cp=plt.contourf(lon,lat,tp,clevs,transform=ccrs.PlateCarree(),cmap='RdPu')
gl = ax.gridlines(draw_labels=True)
plt.contour(lon,lat,ylat,[latmax],transform=ccrs.PlateCarree(),colors='white',
linewidths=1,linestyles='dashed')
cbar = plt.colorbar(cp, ticks=clevs, orientation='horizontal')
cbar.set_label('hPa')
ax.add_feature(states, linewidth=0.8, color='gray')
ax.coastlines('50m', linewidth=0.8,color='gray')

gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator(np.arange(lon1-360,lon2-360+10,10))
gl.ylocator = mticker.FixedLocator(np.arange(lat1,lat2+5,5))
plt.title('Dynamic Tropopause (2PVU) Pressure\n'+myhour+'Z '+fdate)



plt.figure(2)

ax = plt.axes(projection=ccrs.PlateCarree( ))
ax.set_extent([lon1,lon2,lat1,lat2],crs=ccrs.PlateCarree())
clevs2=np.arange(260,400,10)
plt.contour(lon,lat,tp_theta,clevs2,transform=ccrs.PlateCarree(),
    colors='black',linewidths=0.5)
cp=plt.contourf(lon,lat,tp_theta,clevs2,transform=ccrs.PlateCarree(),
    cmap='RdBu_r')
cbar = plt.colorbar(cp, ticks=clevs2, orientation='horizontal')
cbar.set_label('K')

ax.add_feature(states, linewidth=0.8, color='gray')
ax.coastlines('50m', linewidth=0.8,color='gray')
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator(np.arange(lon1-360,lon2-360+10,10))
gl.ylocator = mticker.FixedLocator(np.arange(lat1,lat2+5,5))
plt.title('Dynamic Tropopause (2PVU) Potential Temperature\n'+myhour+'Z '+fdate)


plt.figure(3)

ax = plt.axes(projection=ccrs.PlateCarree( ))
ax.set_extent([lon1,lon2,lat1,lat2],crs=ccrs.PlateCarree())
clevs2=np.arange(-10,11,1)
plt.contour(lon,lat,pv330/1e-6,clevs2,transform=ccrs.PlateCarree(),
    colors='black',linewidths=0.5)
cp=plt.contourf(lon,lat,pv330/1e-6,clevs2,transform=ccrs.PlateCarree(),
    cmap='RdBu_r')
cbar = plt.colorbar(cp, ticks=clevs2, orientation='horizontal')
cbar.set_label('K')

ax.add_feature(states, linewidth=0.8, color='gray')
ax.coastlines('50m', linewidth=0.8,color='gray')
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlocator = mticker.FixedLocator(np.arange(lon1-360,lon2-360+10,10))
gl.ylocator = mticker.FixedLocator(np.arange(lat1,lat2+5,5))
plt.title('Potential Vorticity on the 330K Surface\n'+myhour+'Z '+fdate)

plt.figure(4)
# P-lon cross-section of PV at latitude
# of lowest tropopause
ax = plt.axes()

pv_smooth=gaussian_filter(pv,sigma=1)
theta_smooth=gaussian_filter(theta,sigma=1)

plt.ylim(lev[0],lev[20])
#plt.yscale('log')
clevs=np.arange(2,32,2)
plt.contour(lon-360,lev[0:21],pv_smooth[0:21,maxloc[0,0],:]/1e-6,clevs,
    colors='black')
cp=plt.contourf(lon-360,lev[0:21],pv_smooth[0:21,maxloc[0,0],:]/1e-6,clevs,
    cmap='RdPu')
clevs2=np.arange(260,490,10)
plt.contour(lon-360,lev[0:21],theta_smooth[0:21,maxloc[0,0],:],[330],
    colors='blue',linewidths=1.2)
cs=plt.contour(lon-360,lev[0:21],theta_smooth[0:21,maxloc[0,0],:],clevs2,
    colors='blue',linewidths=0.5)

plt.clabel(cs,inline=1,fontsize=8,fmt='%4.0f')
cbar = plt.colorbar(cp, ticks=clevs, orientation='horizontal')
cbar.set_label('PVU')

plt.title('LON-P Cross-section of PV (shading) and '+r'$\theta$'+
    ' (blue contours) at '+str(latmax)+'N\n'+myhour+'Z '+fdate)

plt.show()

