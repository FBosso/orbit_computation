# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:41:11 2020

@author: Francesco Bosso
@title: Orbit computation from Almanac
@infomation: Almanac of satellite SVN 63, PRN 01 (Block IIR) for 2016/11/ 28
"""

# Space where import functions and libraries 
import numpy as NP
import math as MT
import matplotlib.pyplot as plt
import eta
import Fu as F
from mpl_toolkits.basemap import Basemap
import pandas as PA
from xyz2geo_2 import xyz2geo_2
#------------------------------------------------------------------------------------------------
# 1. vector time
#    it starts from 00:00 end at 23:59 (step 30 second) - TIME

t= NP.array(range(0,86400,30))

#------------------------------------------------------------------------------------------------
# 2. clock offset

dt0= -7.661711424589*10**(-5)
dt1= -3.183231456205*10**(-12)
dt2=  0.000000000000*10**(0)

# second order polinomial relation with the 3 parameters dtsat, asat, bsat
off = dt0 + dt1*t + dt2* t**2

# plot the clock offset
plt.plot(t,off)
plt.title('Clock offset [s]') 
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
plt.show() 

#------------------------------------------------------------------------------------------------
# 3. computation of cordinates

# ORS coordinates

# major semi-axis
a=  (5.153650835037*10**(3))**2 #(meters)

# eccentricy
e=   3.841053112410*10**(-3)

M0= 1.295004883409
GMe = 3.986005*10**(14)

n=MT.sqrt(GMe/a**3)

#vectorial computation
M_t= M0+n*t      # t0 = 0 (time starts at 00:00)

#creation of the ampty vector of eta
eta_t=NP.zeros(len(M_t))

# compuatation of eta for each t using the function eta.f
# filling the vector of eta at each one epoch
for i in range(len(M_t)):
    eta_t[i]=eta.f(M_t[i])
    
# computation of psi (true anomaly)
psi_t=NP.arctan2((NP.sqrt(1-e**2)*NP.sin(eta_t)),(NP.cos(eta_t)-e))

# computation coordintates in ORS
r_t= a*(1-e*NP.cos(eta_t))
#r_t= (a*(1-e**2))/(1+e*NP.cos(psi_t))
x_t = r_t * NP.cos(psi_t)
y_t = r_t * NP.sin(psi_t)

#print of the true anomaly
plt.plot(t,psi_t)
plt.title('true anomaly [rad]')
plt.show() 

# computation of the coordinates in ITRF 

# data in order to apply the linear relation
Omega0= -2.241692424630*10**(-1) #(radians)
Omegadot= -8.386063598924*10**(-9) #(radians/sec)
i0= 9.634782624741*10**(-1) #(radians)
idot= -7.286017777600*10**(-11)  #(radians/sec)
w0= 9.419793734505*10**(-1) #(radians)
wdot= 0.0  #(radians/sec)
OmegaEdot = 7.2921151467*10**(-5) #(radians)

# linear relation (vectorial, all the computation with the same relation)
Omega_t = Omega0+ (Omegadot-OmegaEdot)*t
w_t = w0 + wdot*t
i_t = i0 + idot*t


#------------------------------------------------------------------------------------------------
# 4. covertion from ITRF to Geodetic
# newaxis create an array composed by some bub-array
x_t= x_t[:, NP.newaxis]
coor_ors=NP.zeros((3,len(t)))
coor_itrf=NP.zeros((3,len(t)))
coor_itrf_geo=NP.zeros((3,len(t)))

for i in range(len(t)):
    #rotation matrix for each "epoch"
    R= (F.angle_to_R(-Omega_t[i], 3).dot(F.angle_to_R(-i_t[i], 1))).dot(F.angle_to_R(-w_t[i], 3))
    
    # coordinates ORS for each "epoch"
    coor_ors[0,i]=x_t[i]
    coor_ors[1,i]=y_t[i]
    
    #rotation of ORS in ITRF
    coor_itrf[:,i]=  R.dot(coor_ors[:,i])
    #convertion from geocentric to geodetic
    coor_itrf_geo[:,i] = xyz2geo_2(NP.matrix(coor_itrf[:,i]))
    
    
#------------------------------------------------------------------------------------------------
# 5. Satellite Ground Track
    
# height's avarage [km]    
avg_height=round(NP.mean(coor_itrf_geo[2,:])/1000,4)

# latitude and longitude vector
lat_sat=NP.rad2deg(coor_itrf_geo[0,:])
lon_sat=NP.rad2deg(coor_itrf_geo[1,:])

# plot the Earth surface
fig = plt.figure(figsize=(25, 21), edgecolor='w')
Earth = Basemap(projection='cyl', resolution='l',
                llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, suppress_ticks=False)

# check about reference system
#print(Earth.epsg)
# EPSG : 4326 -> WGS84
#cambiare epsg 3857


# first figure satellite position projected on Earth wrt each epoch
plt.figure(1)

F.draw_map(Earth,0.2)

plt.title("Satellite Ground Track") 

Earth.scatter(lon_sat[1:],lat_sat[1:],marker='.',s=0.5,color='red')

# initial point
Earth.scatter(lon_sat[0],lat_sat[0],marker='s',s=80,color='black') 

# second figure satellite height wrt each epoch
plt.figure(2,figsize=(20,6))


plt.title("ellipsoidic height variations [km] around mean height= %f km" %avg_height) 
h_sat=(coor_itrf_geo[2,:]/1000)-avg_height*NP.ones((1,len(t)))
plt.scatter(t,h_sat,marker='.',s=0.5,color='green')
plt.ticklabel_format(axis='x',style='sci',scilimits=(0,0))

#------------------------------------------------------------------------------------------------
# 6. Save data in a txt file

result=PA.DataFrame({'TIME [s]':t,'X_ORS[m]':NP.round(coor_ors[0,:]),'Y_ORS[m]':NP.round(coor_ors[1,:]),'Z_ORS[m]':NP.round(coor_ors[2,:]),
                     'X_ITRF[m]':NP.round(coor_itrf[0,:]),'Y_ITRF[m]':NP.round(coor_itrf[1,:]),'Z_ITRF[m]':NP.round(coor_itrf[2,:]),
                     'LAT[°]':NP.round(lat_sat,5),'LON[°]':NP.round(lon_sat,5),'H[m]':NP.round(coor_itrf_geo[2,:])})
result.to_csv('Result.txt',sep='|',index=False)

