#fucntion to calculate the density based on the temperature
#and salinity of the water using the equation of state

#Adapted from FVCOM versin 4.3 (Chen et al. 2013) and Clark et al. 2020
import numpy as np

def dens2(S,T):
    #S is the salinity, T is the temperature, RHOF is the density

    SF = S

    TF = np.fmax(T,1.0) 

    RHOF = SF * SF * SF * \
           6.76786136E-6 - SF * SF * 4.8249614E-4 + \
           SF * 8.14876577E-1 - 0.22584586E0

    RHOF = RHOF * (TF*TF* \
           TF*1.667E-8-TF*TF*8.164E-7+ \
           TF*1.803E-5)

    RHOF = RHOF + 1. - TF * TF * \
           TF * 1.0843E-6 + TF * TF * \
           9.8185E-5 - TF * 4.786E-3

    RHOF = RHOF * (SF*SF* \
           SF*6.76786136E-6-SF*SF* \
           4.8249614E-4+SF*8.14876577E-1+3.895414E-2)

    RHOF = RHOF - (TF-3.98) ** 2 * ( \
           TF+283.0) / (503.57*(TF+67.26))

    RHOF=RHOF+1000.       
    return RHOF