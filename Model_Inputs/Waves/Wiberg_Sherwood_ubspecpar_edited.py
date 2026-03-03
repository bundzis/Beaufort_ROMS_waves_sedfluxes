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

    if type(w) == list:
        x = [item ** 2 * h / g for item in w]
    elif type(h) == list:
        x = [w ** 2 * item / g for item in h]
    else:
        x = w ** 2 * h / g

    
    lx = len(x)
    kh = [float(i) for i in x]
    kh = [np.sqrt(i) if i < 1 else i for i in x]
    for i in range(lx):
        y = kh[i]
        t = np.tanh(y)
        y = y - ((y * t - x[i]) / (t + y * (1. - t ** 2)))
        t = np.tanh(y)
        y = y - ((y * t - x[i]) / (t + y * (1. - t ** 2)))
        t = np.tanh(y)
        y = y - ((y * t - x[i]) / (t + y * (1. - t ** 2)))
        kh[i] = y
    
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
    #if type(hs) == float:
     #   hs = np.full((1), hs)
      #  print('1 value hs')
    #if type(tp) == float:
     #   tp = np.full((1), tp)
      #  print('1 value tp')
    #if type(h) == float:
     #   h = np.full((1), h)
      #  print('1 value h')

    
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

    iterr = np.zeros(nt) # an array OLD
    ubr = np.zeros(nt) # an array OLD
    Tbr = np.zeros(nt) # OLD
    #Tbz = Tbr # OLD

    fp = 1. / tp # NEW
    
    # Check
    #print('got here 2')
    
    for i in range(nt):
        m0 = hs[i] ** 2.0 / 16
        f = ffp * fp[i] 
        df = dffp * fp[i] 
        kh = qkhfs(2 * np.pi * f, h[i]) # calling another function - defined in cell above
        
        # Check
        #print('got here 3')

        itcnt = 0 # NEW
        
        if specform == 'D':
            xi = 1
            tol = 1e-3
            m0sfD = 0.0
            
            while abs((m0 - m0sfD) / m0 > tol):
                fpbar = ((m0 * fp[i] ** 4)/((g ** (2) * 6.635e-6))) ** (1 / 0.7) # |* only * works  
                
                gam = 1.7
                
                if fpbar >= 0.05:
                    if fpbar > 0.159:
                        gam = 6.5 + 2.606 * np.log(fpbar) 
                        
                    sig = 0.08 * 0.0013 * ((fpbar) ** (-3))
                else:
                    sig = 0.08 + 0.0013 * ((0.05) ** (-3))
                
                #eterm = - ((ffp - 1)**2)/(2*(sig**2)) # OLD
                a = - ((ffp - 1) ** 2.)
                b = 2.0 * sig ** 2
                fdiv = a/b
                #ee = np.exp(eterm) # OLD
                ee = np.exp(fdiv) # NEW
                t2 = gam ** (ee)
                t1 = - (ffp ** -4.)
                #sffpn = (ffp**(xi)/(ffp**(5)))*np.exp(t1)*t2 # dimensionless spectrum OLD                
                #XD = 1/(sum(sffpn * dffp)) # chi OLD
                a = ffp ** xi
                b = ffp ** 5.
                fdiv2 = a/b
                sffpn = fdiv2 * np.exp(t1) * t2
                dnm = np.sum(sffpn * dffp)
                XD = 1.0 / dnm
                a = m0 * fp[i] ** 4. * XD
                b = f ** 4. * fp[i]
                fdiv3 = a/b

                #sf = (m0 * XD * (fp**(4))*1/((f**(4))*fp))*np.exp(t1)*t2 # eq. 20 OLD
                sf = fdiv3 * np.exp(t1) * t2 # NEW
                
                m0sfD = np.sum(sf * df) 
                #iterr[i] = iterr[i] + 1  # OLD
                itcnt = itcnt + 1
                
                # Check
                #print('got here 4')
                
                if itcnt > 5: # NEW
                    print('Too many iteraions') # NEW
                    break # NEW
             
            ########################## BRI FIX ##############################
            if itcnt == 0:
                fpbar = (m0 * fp[i] ** 4. / (g ** 2 * 6.635e-06)) ** (1 / 0.7)
                gam = 1.7
                if fpbar >= 0.05:
                    if fpbar > 0.159: gam = 6.5 + 2.606 * np.log(fpbar)
                    sig = 0.08 + 0.0013 * fpbar ** -3
                else:
                    sig = 0.08 + 0.0013 * 0.05 ** -3
                a = - ((ffp - 1) ** 2.)
                b = 2.0 * sig ** 2.
                fdiv = a/b
                ee = np.exp(fdiv)
                t2 = gam ** ee
                t1 = -(ffp ** -4.)
                a = ffp ** xi
                b = ffp ** 5.
                fdiv2 = a/b
                sffpn = fdiv2 * np.exp(t1) * t2
                dnm = np.sum (sffpn * dffp)
                XD = 1.0 / dnm
                a = m0 * fp[i] ** 4. * XD
                b = f ** 4. * fp[i]
                fdiv3 = a/b
                sf = fdiv3 * np.exp(t1) * t2
            ######################## END BRI FIX ############################


            #su = (2*np.pi*f/np.sinh(kh))**(2)*sf # eq. 5 OLD
            a = 2 * np.pi * f
            b = np.sinh(kh)
            fdiv4 = (a/b) ** 2
            su = fdiv4 * sf # NEW
            ubrD = np.sqrt(2.0 * np.sum(su * df)) # eq. 6 |* only * works 
            frD = sum(su * f * df) / np.sum(su * df) # eq. 9
            fzD = np.sqrt(np.sum(f ** (2)* su * df) / np.sum(su * df)) 
            TbrD = 1.0 / frD
            TbzD = 1.0 / fzD
            ubr[i] = ubrD  #OLD 
            Tbr[i] = TbrD # OLD
            iterr[i] = itcnt # OLD

        
        if specform == 'J':
            gam = 3.3 # the value of gam can be changed
            xi = 0   
            
            jn = (ffp>1).nonzero() 
            
            sig = list(ffp)
            sig = np.array([0.07 if j < 1 else 0.09 for j in sig])
            #print('sig shape: ', sig.shape)

            a = - ((ffp - 1) ** 2)
            b = 2.0 * sig ** 2.
            fdiv = a/b
            ee = np.exp(fdiv)
            t2 = gam ** (ee)
            t1 = -1.25 * (ffp ** -4) 
            #sffpn = (ffp**(xi)/ffp**(5))*np.exp(t1)*t2 # dimensionless spectrum OLD
            a = ffp ** xi
            b = ffp ** 5.
            fdiv2 = a/b
            sffpn = fdiv2 * np.exp(t1) * t2
            dnm = np.sum(sffpn * dffp)
            #XJ = 1/sum(sffpn*dffp) # chi OLD
            XJ = 1.0 / dnm
            #sf = ((m0*fp**(4)*XJ)/f**(5))*np.exp(t1)*t2 # eq. 20 OLD
            a = m0 * fp[i] ** 4. * XJ
            b = f ** 5.
            fdiv3 = a/b
            sf = fdiv3 * np.exp(t1) * t2
            a = 2 * np.pi * f
            b = np.sinh(kh)
            fdiv4 = (a/b) ** 2
            su = fdiv4 * sf


            
            #su = ((2*np.pi*f/np.sinh(kh))**(2))*sf # eq. 5 OLD

            ubrJ = np.sqrt(2.0 * np.sum(su * df)) # eq. 6            
            frJ = np.sum(f * su * df) / np.sum(su * df) # eq. 9
            fzJ = np.sqrt(np.sum(f ** (2)* su * df) / np.sum(su * df)) 
            TbrJ = 1.0 / frJ
            TbzJ = 1.0 / fzJ
            ubr[i] = ubrJ # OLD
            Tbr[i] = TbrJ # OLD
            iterr[i] = itcnt # OLD
        
    # Check
    #print('Got here 30')
    
    # Return the output from this function
    return ubr, Tbr, iterr # NEW


# TEST
u, t, i = ubspecpar(2, 40, 15, 'D')
print('u: ', u)
print('t: ', t)
print('i: ', i)
