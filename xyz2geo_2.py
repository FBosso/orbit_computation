import math as m
import numpy as np

def xyz2geo_2(xyz):
    lam = m.atan2(xyz[0,1],xyz[0,0])
    
    a=6378137.0000000
    e = m.sqrt(1-(1-1/(298.25722210882711243))**2)
    
    k = 2
    phi = 0
    phi_app = 1
    h = 0
    
    while round(phi,20) != round(phi_app,20):
        phi_app = phi
        k = k-1
        if k == 1:
            n = (xyz[0,2])
            d = m.sqrt(xyz[0,0]**2 + xyz[0,1]**2)
            phi = m.atan2(n,d)
            N =(a)/m.sqrt(1-(e**2 * (m.sin(phi))**2))
           
        else:
            n = (xyz[0,2] + e**2 * (N*m.sin(phi)))
            d = m.sqrt(xyz[0,0]**2 + xyz[0,1]**2)
            phi = m.atan2(n,d)
            N =(a)/m.sqrt(1-(e**2 * (m.sin(phi))**2))
            h = (m.sqrt(xyz[0,0]**2 + xyz[0,1]**2)/(m.cos(phi)))- N
    
    return np.matrix([phi,lam,h]) 