#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 13:08:55 2022

@author: brun1463
"""

########################## Calculate Pwave_bot and Uwave_rms #############################
# Calculate bottom orbital velocity using a parametric wave spectrum from Wiberg & 
# Sherwood 2008 Appendices D & E
#
# This is a python version of the matlab scripts in Appendices D & E of the Wiberg and 
# Sherwood 2008 paper - so everything in this script is taken exactly from that paper.
#
# The purpose of this script is to create the function from Wiberg and Sherwood that 
# calculates bottom orbital velocity and bottom wave period from significant wave height, 
# peak wave period, and water depth.
#
# Notes:
# - Matlab starts indexing at 1 (python starts at 0)
# - Should work for inputs that are arrays or single values
# - |* means check if it should be * or @
#  - All of these were checked an appear to be correct since only * works (due to 
#    incompatible sizes for matrix multiplication)
#############################################################################################


# Load in the packages 
import numpy as np
import math

# Define the function to find kh
def qkhfs(w, h):
    """
    This function is a quick iterative calculation of kh in dispersion 
    relationship
    
    kh = qkhfs(w, h)
    
    Input:
    - w: Angular wave frequency = 2*pi/T where T = wave period (1/s)
    - h: Water depth (m)
    
    Returns:
    - kh = wavenumber * depth 
    
    Either w or h can be a vector, but not both.
    Hard-wired for MKS units.
    Orbital velocities from kh are accurate to 3e-12 !
    
    RL Soulsby (2006) "Simplified calculation of wave orbital velocities"
    HR Wallingford Report TR 155, February 2006, Eqns. 12a - 14
    
    csherwood@usgs.gov
    Sept 10, 2006
    
    Converted to python here by Brianna Undzis, 10 August 2022
    """

    g = 9.80665
    x = (w**(2)*h)/g # |* only * works 
    
    y = np.sqrt(x) * (np.where(x < 1,1,0)) + x * (np.where(x >= 1,1,0)) 

    t = np.tanh(y)
    y = y - ((y * t  -x)/(t + y * (1. - t ** (2)))) 
    t = np.tanh(y)
    y = y - ((y * t - x)/(t + y * (1. - t ** (2)))) 
    t = np.tanh(y)
    y = y - ((y * t - x)/(t + y * (1. - t ** (2)))) 
    
    kh = y
    
    return kh


# Define the function to find bottom wave period and bottom orbital velocity
def ubspecpar(hs, tp, h, specform='D'):
    """
    Calculate ubr and Tbr from hs and tp using parametric spectrum
    
    The inputs can be given either as single values or numpy arrays.
    
    Input:
    - hs: Significant wave height (m)
    - tp: Peak period (s)
    - h: Water depth (m)
    - specform: spectral formulation to use
        specform = 'D' for Donelan spectrum (default)
        specform = 'J' for JONSWAP spectrum
    
    Returns:
    - ubr: representative bottom orbital velocity (m/s)
    - Tbr: representative bottom wave period (s)
    - iter: number of iterations if Donelan spectrum is chosen
        - If iter > 5 the calculation did not converge
        - iter = 0 for the JONSWAP spectrum
        
    Patricia Wiberg, UVa
    Last Modified 9 Mar 2007
    
    Converted to python here by Brianna Undzis, 10 August 2022
    """
    # Load in the packages
    import numpy as np
    import math
    
    # Check
    #print('got here 0')
    
    # Make it so this still works if only one value is entered
    if type(hs) == int:
        hs = np.full((1), hs)
        print('1 value hs')
    if type(tp) == int:
        tp = np.full((1), tp)
        print('1 value tp')
    if type(h) == int:
        h = np.full((1), h)
        print('1 value h')

    
    g = 9.81
    dffp = 0.01
    
    ffp = np.arange(0.2, 5+dffp, dffp) # an array
    
    # Check
    #print('got here 1')
    
    # Returns the length of the largest array dimension in X
    try:
        nt = max(tp.shape)
    except:
        nt = 1
    
    iterr = np.zeros((nt, 1)) # an array 
    
    ubr = np.zeros((nt,1)) # an array 
    
    Tbr = np.zeros((nt,1)) # an array
    
    # Check
    #print('got here 2')
    
    for i in range(nt):
        m0 = hs[i]**2 / 16
        fp = 1/tp[i]
        f = ffp * fp 
        df = dffp * fp 
        kh = qkhfs(2*np.pi*f, h[i]) # calling another function - defined in cell above
        
        # Check
        #print('got here 3')
        
        if specform == 'D':
            xi = 1
            tol = 1e-3
            m0sfD = 0
            
            while abs((m0-m0sfD)/m0 > tol):
                fpbar = ((m0 * fp**4)/((g**(2)*6.635e-6)))**(1/.7) # |* only * works  
                
                gam = 1.7
                
                if fpbar >= 0.05:
                    # alpha = 0.0165 * fpbar**(0.55)
                    if fpbar > 0.159:
                        gam = 6.5 + 2.606*math.log(fpbar) # |* only * works
                        
                    sig = 0.08 + 0.0013*((fpbar)**(-3))
                else:
                    # alpha = 0.0165 * 0.05**(0.55)
                    sig = 0.08 + 0.0013*((0.05)**(-3))
                
                eterm = - ((ffp - 1)**2)/(2*(sig**2))
                ee = np.exp(eterm)
                t2 = gam**(ee)
                t1 = - (ffp)**(-4)
                sffpn = (ffp**(xi)/(ffp**(5)))*np.exp(t1)*t2 # dimensionless spectrum
                XD = 1/(sum(sffpn * dffp)) # chi
                
                sf = (m0 * XD * (fp**(4))*1/((f**(4))*fp))*np.exp(t1)*t2 # eq. 20 
                
                m0sfD = sum(sf*df) 
                iterr[i] = iterr[i] + 1
                
                # Check
                #print('got here 4') 
                
                if iterr[i] > 5:
                    break # this is the same in python and matlab

                su = (2*np.pi*f/np.sinh(np.float128(kh)))**(2)*sf # eq. 5
                ubrD = np.sqrt(2*sum(su*df)) # eq. 6 |* only * works 
                frD = sum(su*f*df)/sum(su*df) # eq. 9
                fzD = np.sqrt(sum(f**(2)*su*df)/sum(su*df)) 
                TbrD = 1/frD
                TbzD = 1/fzD
                ubr[i] = ubrD
                Tbr[i] = TbrD
        
        if specform == 'J':
            gam = 3.3 # the value of gam can be changed
            xi = 0 
            #XJ = 3.283  
            
            jn = (ffp>1).nonzero() 
            
            sig = 0.07*np.ones_like((ffp), float) # |* only * works 
            print('sig shape: ', sig.shape)

            
            sig[jn] = 0.09
            eterm = - ((ffp-1)**2)/(2*sig**(2)) 
            ee = np.exp(eterm)
            t2 = gam**(ee)
            t1 = - 1.25*(ffp)**(-4) 
            sffpn = (ffp**(xi)/ffp**(5))*np.exp(t1)*t2 # dimensionless spectrum
            XJ = 1/sum(sffpn*dffp) # chi
            sf = ((m0*fp**(4)*XJ)/f**(5))*np.exp(t1)*t2 # eq. 20

            
            su = ((2*np.pi*f/np.sinh(np.float128(kh)))**(2))*sf # eq. 5
            ubrJ = np.sqrt(2*sum(su*df)) # eq. 6 |* only * works            
            frJ = sum(f*su*df)/sum(su*df) # eq. 9
            fzJ = np.sqrt(sum(f**(2)*su*df)/sum(su*df)) 
            TbrJ = 1/frJ
            TbzJ = 1/fzJ
            ubr[i] = ubrJ
            Tbr[i] = TbrJ
        
    # Check
    #print('Got here 30')
    
    # Return the output from this function
    return ubr[:,0], Tbr[:,0], iterr[:,0]




