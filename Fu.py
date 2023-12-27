"""
Created on Fri Mar 20 14:35:17 2020

@author: Marco

LISTA FUNZIONI DA IMPLEMENTARE
"""


import numpy as np
import math as m
from itertools import chain

#implementazione matrice di rotazione.
#axis has to contain the int 1,2 or3

def angle_to_R(alpha,axis):
    if axis==1:
        R=np.array([[1,0,0],
                    [0, m.cos(alpha), m.sin(alpha)],
                    [0,-m.sin(alpha), m.cos(alpha)]])
    if axis==2:
        R=np.array([[m.cos(alpha), 0, -m.sin(alpha)],
                    [     0,       1,       0],
                    [m.sin(alpha),-0, m.cos(alpha)]])
    if axis==3:
        R=np.array([[ m.cos(alpha),m.sin(alpha),0],
                    [-m.sin(alpha),m.cos(alpha),0],
                    [0,0,1]])
                   
   
    return R;



#conversione angoli

def dms2rad(dms):
    
    sex=dms[0]+dms[1]/60+dms[2]/3600
    rad=sex*m.pi/180
    
    return rad

def rad2dms(rad):
    
    sex=rad*180/m.pi
    d=m.floor(sex)
    M=m.floor((sex-d)*60)
    s=((((sex-d)*60)-M)*60)
    dms=np.array([d,M,s])
    
    return dms


#from GC to LC
def xyz2enu(xyz_0,xyz):
    
    #get phi and lambda coordinates of the origin of the Local RS
    geo=F.xyz2geo(xyz_0)
    phi=geo[0]
    lam=geo[1]
        
    #compute the Rotation Matrix
    R=np.array([[-m.sin(lam),               m.cos(lam),               0],
                [-m.sin(phi)*m.cos(lam),   -m.sin(phi)*m.sin(lam), m.cos(phi)],
                [ m.cos(phi)*m.cos(lam),    m.cos(phi)*m.sin(lam), m.sin(phi)]])
    
    
    #compute the LC thanks to the equation LC=R*(xyz-xyz_0)
    delta=xyz-xyz_0
    delta=delta[:, np.newaxis]
    
    enu= R.dot(delta)
    enu=np.squeeze(enu)
    return enu


# from LC to GC--> slide 66

def enu2xyz(xyz_0,enu):
    
    enu= enu[:,np.newaxis]
    xyz_0=xyz_0[:,np.newaxis]
    #get phi and lambda coordinates of the origin of the Local RS
    geo=F.xyz2geo(xyz_0)
    phi=geo[0]
    lam=geo[1]
    #compute the Rotation Matrix
    R=np.array([[-m.sin(lam),               m.cos(lam),               0],
                [-m.sin(phi)*m.cos(lam),   -m.sin(phi)*m.sin(lam), m.cos(phi)],
                [ m.cos(phi)*m.cos(lam),    m.cos(phi)*m.sin(lam), m.sin(phi)]])    
    
    #compute the GC thanks to the equation GC=R_t*LC + GC_0
    #(where GC_0 are the global cartesian coordinates of the origin of the local RS)
    xyz=(R.transpose()).dot(enu) + xyz_0
    
    xyz=np.squeeze(xyz)
    return xyz

# from GC to Geo--> slide 71
def xyz2geo(xyz):
    
    a=6378137.0000000
    f=1/(298.25722210882711243)

    b      = a*(1-f)
    e_2    = 1-(1-f)**2
    e_b_2  = 1/((1-f)**2) -1
    r      =m.sqrt( m.pow(xyz[0],2) +m.pow(xyz[1],2) )
    psi    =m.atan2( xyz[2] ,  r*m.sqrt(1-e_2))
    lam    =m.atan2( xyz[1] , xyz[0])
    phi    =m.atan2( xyz[2]+e_b_2*b*m.pow(m.sin(psi),3)  , r-e_2*a*m.pow(m.cos(psi),3))
    R_n    = a/(m.sqrt(1-e_2*m.pow(m.sin(phi),2)))
    h      = r/m.cos(phi) - R_n
    
    geo=np.array([phi,lam,h])
    return geo







def draw_map(m, scale):       
    
    # draw a shaded-relief image
        m.shadedrelief(scale=scale)
    
    # lats and longs are returned as a dictionary
        lats = m.drawparallels(np.linspace(-90, 90, 13))
        lons = m.drawmeridians(np.linspace(-180, 180, 13))
   
    # keys contain the plt.Line2D instances
        lat_lines = chain(*(tup[1][0] for tup in lats.items()))
        lon_lines = chain(*(tup[1][0] for tup in lons.items()))
        all_lines = chain(lat_lines, lon_lines)
    
    
    # cycle through these lines and set the desired style
        for line in all_lines:
        
            line.set(linestyle='-', alpha=0.3, color='w')           
   
        
        return
    
    
    
def xyz2geo_p(xyz):

    a=6378137.0000000

    e_2    =0.0818191908426215**2
    r      =m.sqrt( m.pow(xyz[0,0],2) +m.pow(xyz[0,1],2)+m.pow(xyz[0,2],2) )
    lam    =m.atan2( xyz[0,1] , xyz[0,0])
    phi_C  =m.atan(xyz[0,2]/(m.sqrt(m.pow(xyz[0,0],2) +m.pow(xyz[0,1],2))))
    psi    =m.atan( m.tan(phi_C)/m.sqrt(1-e_2))

    phi    =m.atan((r*m.sin(phi_C)+e_2*a/m.sqrt(1-e_2)*m.sin(psi)**3)/(r*m.cos(phi_C)-e_2*a*m.cos(psi)**3))
    
    N    = a/(m.sqrt(1-e_2*(m.sin(phi)**2)))
    h      = r*m.cos(phi_C)/m.cos(phi) - N
    
    geo=np.matrix([phi,lam,h])
    return geo


