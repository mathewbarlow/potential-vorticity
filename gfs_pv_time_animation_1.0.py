#
# run on python 3.7
#
# python code to calculate the dynamic tropoapuse from online GFS data over
# a given date range (has to be within 7 days of the current time) and make
# a 3D animation.  As the data is accessed online and the calculations are 
# poorly coded, the program can take a while to run.
#
# NB:
# 1. Can take hours to run, due to both slow data access and poor coding
# 2. Does not fail gently if there is a data access problem (this could
#    be fixed)
# 3. Not currently able to handle a longitude range that crosses 0 (this 
#    could also be fixed, fairly easily)
# 4. Very poorly coded and not commented:  may cause nausea, vomiting, hair
#    loss, and uncontrollable crying.  You have been warned!
#
# The date and lat-lon range can be set below.
#
# (poorly) coded by Mathew Barlow
# initial release: 7 Feb 2020
# updated: 9 Feb 2020
#
# This code has *not* been extensively tested and has been, in part, 
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
from datetime import timedelta

real_current_time = datetime.now().strftime("%H:%M:%S")
print(real_current_time)


# VALUES TO SET *************************************************
# set date, lat-lon range, PV-value definition of tropopause
# and directory to put images
date_start='20200208'
hour_start='18'
date_end='20200209'
hour_end='06'
(lat1,lat2)=(20,60)
(lon1,lon2)=(-130,-60)
tpdef=2   # definition of tropopause in PVU
dirout= '/Users/mathew_barlow/downloads/'
#****************************************************************
if lon2<0:
    lon1=lon1+360
    lon2=lon2+360


start_date=datetime.strptime(date_start+hour_start, '%Y%m%d%H')
end_date=datetime.strptime(date_end+hour_end, '%Y%m%d%H')
diff = end_date-start_date
days, seconds = diff.days, diff.seconds
hours = days * 24 + seconds // 3600
nt = np.int((hours+6)/6)

dates = start_date
it=0
while(it<nt-1):
    lag_date=start_date+timedelta(hours=(it+1)*6)
    dates=np.append(dates,lag_date)
    it=it+1

images = []
filenames = []
    
it=nt-1
while(it>=0):
    
    mydate=dates[it].strftime('%Y%m%d')
    myhour=dates[it].strftime('%H')
    
    time_out=mydate+myhour

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
    
    # the following lines shouldn't be necessary but appear to help
    # with data access
    
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

    # without data access issues, this would be the way to go:
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
        
    
    # could also use the built-in gradient operator for derivatives,
    # can't remember why I abandoned this
    #lev3=lev.reshape(lev.size,1,1)
    #ddp_theta=np.gradient(theta,lev3*100,axis=0)
    #ddx_theta=np.gradient(theta,axis=2)/dx
    #ddy_theta=np.gradient(theta,axis=1)/dy
    
    # some spatial smoothing
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
    
    # get date for plotting
    fdate=datetime.strptime(mydate, '%Y%m%d').strftime('%d %b %Y')
    stime = myhour+'Z '+fdate
    
    plt.close(fig='all')
    
    print('done with calculations for '+stime+' '+np.str(it)+' out of '+np.str(nt))
    
    plt.figure(it,figsize=plt.figaspect(0.5))

    pressfc_smooth=gaussian_filter(pressfc,sigma=1)
    ax=plt.gca(projection='3d')

    surf=ax.plot_surface(xlon,ylat,tp,cmap="coolwarm",alpha=1,
                    rstride=1,cstride=1,
                    vmin=50,vmax=650,
                    linewidth=0, antialiased=False)

    ax.plot_surface(xlon,ylat,pressfc_smooth/100,color="lightgray",
                    rstride=1,cstride=1,
                    linewidth=0, antialiased=False)
                    
    ax.set_zlim(1000,100)
    ax.set_xlim(lon1,lon2)
    ax.set_ylim(lat1,lat2)
    ax.view_init(elev=80,azim=-90)


    plt.title('2PVU Dynamic Tropopause over topography\n'+stime)
    plt.colorbar(surf,shrink=0.5)
    
                    
    dirout = '/Users/mathew_barlow/downloads/'
    filename=pdir +'tp_3D'+ time_out +'.png'
    plt.savefig(filename, bbox_inches='tight')
    
    images.append(imageio.imread(filename))
    filenames.append(filename)
        
    plt.clf()

        
    it=it-1
    
imageio.mimsave(pdir+'tp_time_anim.gif', images)
#for file in filenames:
#    os.remove(file)

plt.close('all')

real_current_time = datetime.now().strftime("%H:%M:%S")
print(real_current_time)