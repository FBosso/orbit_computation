# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:16:33 2020

@author: Marco
"""
import math as MT

# this function implements the iterative computation of the value eta,
# according to the formula : eta(t)= M(t) + e* sin(eta(t))
# e is considered constant.


def f(M_t):
    
    e = 3.841053112410*10**(-3)
    
    eta_f=M_t
    delta=1
    i=0
    
    while (delta > 10**(-12) and i<12) :
        eta_i=eta_f
        eta_f=M_t+e*MT.sin(eta_i)
        delta= (eta_f - eta_i) % 2*MT.pi
        i=i+1
       
    return eta_f
        