#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:45:19 2021

@author: jbclark8
"""
#  main program for ICM-DOM-PD for python
#  a highly modularized code based on Cerco and Coel 1993 estuarine model
#  of Chesapeake Bay

#  Modified from Clark et al. 2020 to include more complex
#  light reactions and photochemistry
#  this main program contains a vertical solution to be calculated from two boundary conditions
#  for estuarine circulation
#  but Q can be set to zero to make it a stable system
#  See readme.txt for all the details on installing appropriate software and
#  some potential modifications to parameterization for different systems

#  Version 0.0.2

#  JB Clark, Aug 2021

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Density
import functions
import DefVars as dv
import Light
import Sediment as sed
import Nutrients as nt
import DOM
import POM
import Algae
import Oxygen
import time


# sys.path
# sys.path.append('/Users/jbclark8/Desktop/pyICM')
#############################################################################
# start main portion of code
# initialize all variables that are global, mostly parameters
errorcode = dv.InitializeVariables()
if(errorcode == 1):
    print('All variables were not defined, there is a problem and ICM is\
          stopping')
elif(errorcode == 0):
    print('All global variables have been defined')

# set the spatial dimensions
# Latitude, Longitude, and Year
# We are using a Climatology so year can be whatever
LAT = 38.
LON = -76.0
yy = 2020

# number of layers to split the depth
nLayers = 10
# depth of the box, m
dv.Z = 15.
# cross channel distance, m
dv.X = 1.e4
# along channel distance, m
# will be calculated from the difference in input coords
dv.Y = 2.e4
# now divide the box into layers that are evenly spaced
Z = np.ones(nLayers) * dv.Z / nLayers
# but need to split to have an odd number for integrating the light
ZZ = np.zeros(nLayers+1)
ZZ[0] = 0.
ZZ[1:nLayers + 1] = np.cumsum(Z)

# read in the time series data  of water quality inputs
WQ_fname = '/Users/jbclark8/Desktop/pyICM/Inputs/Chesapeake/CBP/CB1.1_Climatology.csv'
myTS_up = pd.read_csv(WQ_fname)

# read in spectral information for light attenuation
Light_fname = '/Users/jbclark8/Desktop/pyICM/Inputs/Chesapeake/CBay_Spectral.csv'
myLT = pd.read_csv(Light_fname)

# read in surfacec forcing information for light and weather
weather_fname = '/Users/jbclark8/Desktop/pyICM/Inputs/Chesapeake/Weather/CB1.1_Weather.csv'
mySFC = pd.read_csv(weather_fname)

# calculate the domain volume, specified in DefVars
Volume = functions.CalcVolume(dv.X, dv.Y, dv.Z)
print('Volume =', Volume*1e-9, 'km^3')

# now set up the time vector for integration, using the forcing file time

modtimein = myTS_up.Time
lasttime = modtimein[len(modtimein)-1]
print('Start and End Days are = ', modtimein[0]/86400, lasttime/86400)

# set the time step here
DT = 3600.
Nsteps = int(lasttime/DT)
# now set up a time array with the appropriate time step
modtime = np.linspace(modtimein[0], lasttime, Nsteps)
current_time = modtime[0]

print('My Time Step is = ', DT, 'Seconds')
# initialize and define variables

# phytoplankton 1 and 2 biomass (g C m^-3)
mylen = len(modtime)
B1 = np.zeros((mylen, nLayers))
B1[0] = 0.1
B2 = np.zeros((mylen, nLayers))
B2[0] = 0.1

# phytoplankton 1 and 2 chla
# convert carbon to chla in micrograms L6-1
chla1 = B1 / dv.cchla1 * 1e3
chla2 = B2 / dv.cchla2 * 1e3

chla = chla1 + chla2

# #ammonium and nitrate concentration (g N m^-3)
NH4 = np.zeros((mylen, nLayers))
NH4[0] = 0.05
rivNH4 = 0.0

NO3 = np.zeros((mylen, nLayers))
NO3[0] = 0.1
rivNO3 = 0.0

# dissolved oxygeen
DO2 = np.zeros((mylen, nLayers))
DO2[0] = 10.0

# temp var for DON remineralization
MNLDON = 0.00
##############################################################################
# calculate and set up things for light
# wavelengths, nm
WL = myLT.Lambda

# read in total downwelling shortwave flux and scale by 0.43 to remove IR
EdIn = mySFC.dswrf*0.43
EdUp = 0.
# set up an array for KD, PAR, and NP_Total to collect and pass back to main
# from the light attenuation function
KD = np.zeros((mylen, nLayers, len(WL)))

PAR = np.zeros((mylen, nLayers))
KD_PAR = np.zeros((mylen, nLayers))

NP_total = np.zeros((mylen, nLayers, len(WL)))

# absorption due to water, m^-1
aWater = myLT.aW
# mass specific absorption for each colored DOC, m^2 gC^-1
aC1 = myLT.aCDOC1
aC2 = myLT.aCDOC2
aC3 = myLT.aCDOC3
# take the averge for now
aCDOC = np.average(np.column_stack((aC1, aC2, aC3)), 1)
# mass specific absorption due to chla, m^2 mg chla^-1
aPhi = myLT.aPhi
# mass specific absorption due to particles, m^2 g^-1
aP = myLT.aP
# mass specific backscattering due to particles, m^2 g^-1
bbP = myLT.bbp
# spectral distribution of light, nm^-1
SpecDis = myLT.Spec_dist
# now find where WL == 400 and WL == 700 to mark off the PAR
dv.PAR_Id[0] = (np.abs(WL - 400.)).argmin()
dv.PAR_Id[1] = (np.abs(WL - 700.)).argmin()
#######################################################
# interpolate all forcing to the model time step
# Water Temperature
Tin = myTS_up.WTEMP
T = np.interp(modtime, modtimein, Tin)
# Water Salinity
Sin = myTS_up.SALINITY
S = np.interp(modtime, modtimein, Sin)

# wind velocity
Uin = mySFC.uwnd[0:365]
Vin = mySFC.vwnd[0:365]
Speed = np.sqrt(Uin ** 2 + Vin ** 2)
Uwind = np.interp(modtime, modtimein, Speed)

# Nitrogen
rivNH4in = myTS_up.NH4F
rivNO3in = myTS_up.NO23F

rivNH4 = np.interp(modtime, modtimein, rivNH4in)
rivNO3 = np.interp(modtime, modtimein, rivNO3in)


# Convert upstream chl a concentration into algae 1 and 2 carbon
RivAlgae1in = myTS_up.CHLA*0.5*dv.cchla1*1e-3
RivAlgae2in = myTS_up.CHLA*0.5*dv.cchla1*1e-3
RivAlgae1 = np.interp(modtime, modtimein, RivAlgae1in)
RivAlgae2 = np.interp(modtime, modtimein, RivAlgae2in)


# get the flushing rate with the flow time series and the box dimensions
Q = np.ones(mylen)*500.

# print('Flushing Rate = ',FlushingRate,'per day')

# calculate the change in concentration due to river inputs
RiverISSin = myTS_up.TSS
RiverISS = np.interp(modtime, modtimein, RiverISSin)

# #colored and non-colored DOC concentration (g C m^-3)
CDOC = np.ones((mylen, nLayers))
# NCDOC=np.zeros(mylen)
MNLDOC = 0.0
# #colored and non-colored DON concentration (g N m^-3)
# CDON=np.zeros(mylen)
# NCDON=np.zeros(mylen)

# #Labile and Refractory POC and PON in g C and g N m^-3
LPOC = np.zeros((mylen, nLayers))
# LPOC[0] = 0.5
RPOC = np.zeros((mylen, nLayers))
# RPOC[0] = 1.0
# LPON=np.zeros(mylen)
# RPON=np.zeros(mylen)

# inorganic suspended sediment (g m^-3)
ISS = np.zeros((mylen, nLayers))
ISS[0] = 10.

# TSS = np.zeros(mylen)
TSS = ISS + (LPOC + RPOC) * 2.5
# Algae1Out=pd.DataFrame(index=range(mylen),\
#                       columns=['NPP','GPP','DTB1','NH4A','NO3A','PN','NL','P1','FI'])

# Algae2Out=pd.DataFrame(index=range(mylen),\
#                       columns=['NPP','GPP','DTB2','NH4A','NO3A','PN','NL','P2','FI'])

i = 1

# DOY=np.zeros(mylen)
SZA = np.zeros(mylen-1)
# Sreflectance=np.zeros(mylen)
# while current_time < modtimein[1]:
# put diagnostic variables here
# many of these rates, mostly related to phytoplankton growth currently,
# are used in the other functions, especially DO2

FI1 = np.zeros((mylen, nLayers))
FI2 = np.zeros((mylen, nLayers))
NL1 = np.zeros((mylen, nLayers))
NL2 = np.zeros((mylen, nLayers))
NPP1 = np.zeros((mylen, nLayers))
NPP2 = np.zeros((mylen, nLayers))
P1 = np.zeros((mylen, nLayers))
P2 = np.zeros((mylen, nLayers))
PR1 = np.zeros((mylen, nLayers))
PR2 = np.zeros((mylen, nLayers))
PN1 = np.zeros((mylen, nLayers))
PN2 = np.zeros((mylen, nLayers))

NT = np.zeros((mylen, nLayers))

while current_time < modtime[len(modtime)-1]:
    start = time.time()
    JDAY = current_time/86400
    MyDay = np.floor(JDAY)

    # now loop over the layers
    k = 0
    
    # get the days downwelling irradiance
    EdTop = EdIn[MyDay]*SpecDis
    
    for k in range(nLayers):

        # calculate the total absorption and the total scattering
        aTotal = (aWater + chla[(i-1, k)] * aPhi +
                  CDOC[(i-1, k)] * aCDOC +
                  TSS[(i-1, k)] * aP)
        bTotal = TSS[(i-1, k)] * bbP
        # call the light attenuation functions
        # have to set it so it calls once for the surface layer and then uses
        # the previous calculated downwelling irradiance in the next layer
        Zin = (ZZ[k] + ZZ[k+1]) / 2
        if k == 0:
            EdTop = EdTop
        else:
            EdTop = EdUp        

        KD[(i-1, k, )], PAR[(i-1, k)],\
            NP_total[(i-1, k, )], KD_PAR[(i-1, k)], EdUp = \
            Light.Light_Attenuation(WL, EdTop, aTotal,
                                    bTotal, Zin,
                                    JDAY, LAT, LON, yy)

        FRate = functions.CalcFlushingRate(Q[i-1], Volume)
        # calculate the change in concentration due to biogeochemical processes

        # first algal growth and death
        Algae1 = Algae.DTB1(B1[(i-1, k)], NH4[(i-1, k)],
                            NO3[(i-1, k)], T[(i-1)],
                            PAR[(i-1, k)], FRate, RivAlgae1[(i-1)],
                            Q[(i-1)], Volume)
        Algae2 = Algae.DTB2(B2[(i-1, k)], NH4[(i-1, k)],
                            NO3[(i-1, k)], T[(i-1)],
                            PAR[(i-1, k)], FRate, RivAlgae1[(i-1)],
                            Q[(i-1)], Volume)
        # get some of the algae rates
        FI1[(i, k)] = Algae1.FI
        FI2[(i, k)] = Algae2.FI
        NL1[(i, k)] = Algae1.NL
        NL2[(i, k)] = Algae2.NL
        NPP1[(i, k)] = Algae1.NPP
        NPP2[(i, k)] = Algae2.NPP
        P1[(i, k)] = Algae1.P1
        P2[(i, k)] = Algae2.P2
        PR1[(i, k)] = Algae1.PR1
        PR2[(i, k)] = Algae2.PR2
        PN1[(i, k)] = Algae1.PN
        PN2[(i, k)] = Algae2.PN

        # next inorganic nitrogen
        NH4A = Algae1.NH4A+Algae2.NH4A
        NO3A = Algae1.NO3A+Algae2.NO3A

        DTNO3 = nt.DTNO3(NH4[(i-1, k)], NO3[(i-1, k)], T[(i-1)], DO2[(i-1, k)], NO3A,
                         Q[(i-1)], Volume, FRate, rivNO3[(i)])
        DTNH4 = nt.DTNH4(NH4[(i-1, k)], NO3[(i-1, k)], T[(i-1)], DO2[(i-1, k)], MNLDON,
                         NH4A, Q[(i-1)], Volume, FRate, rivNH4[(i)])
        # nitrification
        NT[i,k] = nt.Nitrification(NH4[(i-1, k)], T[(i-1)], DO2[(i-1, k)])

        # inorganic sediment
        DTISS = sed.deltaISS(Q[(i-1)], Volume, FRate,
                             RiverISS[(i-1)], ISS[(i-1, k)])

        # next calculate dissolved oxygen
        if k == 0:
            Uwindin = Uwind[i]
        else:
            Uwindin = 0.

        DTDO2 = Oxygen.DOXG(DO2[(i - 1, k)], T[(i)], S[(i)], Z[k],
                            PN1[(i, k)], PN2[(i, k)],
                            P1[(i, k)], P2[(i, k)],
                            PR1[(i, k)], PR2[(i, k)],
                            B1[(i, k)], B2[(i, k)],
                            NT[(i, k)], Uwindin, MNLDOC, 0)

        DO2[(i, k)] = DO2[(i-1, k)] + DTDO2 * DT

        #  update the concentrations
        ISS[(i, k)] = ISS[(i-1, k)]+DTISS*DT

        B1[(i, k)] = B1[(i-1, k)]+Algae1.DTB1*DT
        B2[(i, k)] = B2[(i-1, k)]+Algae2.DTB2*DT

        chla1[(i, k)] = B1[(i, k)] / dv.cchla1 * 1e3
        chla2[(i, k)] = B2[(i, k)] / dv.cchla2 * 1e3

        chla[(i, k)] = chla1[(i, k)] + chla2[(i, k)]

        NH4[(i, k)] = NH4[(i-1, k)]+DTNH4.DTNH4*DT
        NO3[(i, k)] = NO3[(i-1, k)]+DTNO3.DTNO3*DT

        TSS[(i, k)] = ISS[(i, k)] + B1[(i, k)] + B2[(i, k)]\
            + (LPOC[(i, k)] + RPOC[(i, k)]) * 2.5

        k = k + 1
# update the time
    current_time = modtime[i]
    i = i+1
    end = time.time()
    print(end-start)

##############################################################################
plot = plt.plot
contourf = plt.contourf
surf = plt.surface

# now make some plots
plt.plot(modtime/86400, TSS[:,2])
plt.xlabel('Time')
plt.ylabel('$ISS (g m^{-3})$')
plt.show()

plt.plot(modtime/86400, B1)
plt.xlabel('Time')
plt.ylabel('$Alg 1 (gC m^{-3})$')
plt.show()

plt.plot(modtime/86400, B2)
plt.xlabel('Time')
plt.ylabel('$Alg 2 (gC m^{-3})$')
plt.show()


plt.plot(modtime/86400, NH4)
plt.xlabel('Time')
plt.ylabel('$NH4 (gNm^{-3})$')
plt.show()

plt.plot(modtime/86400, NO3)
plt.xlabel('Time')
plt.ylabel('$NO3 (gNm^{-3})$')
plt.show()
