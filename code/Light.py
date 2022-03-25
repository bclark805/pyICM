#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:45:19 2021

@author: jbclark8
"""


# Module for calculating the underwater light field
# adapted from Clark et al. 2020, following Gallegos et al. 2006
# JB Clark, NASA GSFC, Aug 2021

import numpy as np
import DefVars as dv
from pysolar.solar import get_altitude
import datetime
import scipy.integrate as integrate
import OASIM


def Light_Attenuation(WL, Ed, aTotal, bTotal, Z,
                      JDAY, Latitude, Longitude, yy,T,S):
    # WL is an array of the wavelengths, nm
    # Ed is the spectral light intensity for this day, W m^-2 nm^-1
    # aTotal is the total spectral absorption of size (WL), m^-1
    # bTotal is the total backscatter of size (WL), m^-1
    # Z1 is the average depth of this layer
    # JDAY is the Julian Day (e.g. DOY)
    # Latitude and longitude are self explanatory (degrees)
    # T and S are temperature and salinity
    rlwn = np.zeros(len(WL))
    Rrs = np.zeros(len(WL))

    # KD_out is the spectral attenuation coefficient, m^-1
    # PAR is the number of photons/ day in the PAR range (400-700 nm)
    # NP_Total is the number of photons / second over the entire range per nm

    # get the solar zenith angle, the reflection (0-100), and DOY
    SZA, SREFLECT, DOY, CTHTAA = Declination(JDAY, Latitude, Longitude, yy)

    SZAa=180*CTHTAA/np.pi
    # now calculate the spectral KD using the Lee et al. 2013 Algo
    KD_out = KD_LEE(SZA,aTotal, bTotal)
    # scaled Ed  based on the reflection, which is 100
    # when the sun is below the horizon
    Ed_Top = Ed * (1. - SREFLECT / 100.)
    # optical depth for the surface layer
    OptDepth = KD_out * Z
    # now calculate the light at the bottom, based on exponential decay
    Ed_Bottom = Ed_Top * np.exp(-OptDepth)
    # now calculate the average light based on the optical depth
    Ed_Avg = (Ed_Top - Ed_Bottom) / OptDepth

    # convert Ed_Avg to total number of photons per second
    NP_Total = Ed_Avg * WL * dv.NP_sec / dv.Avo_N

    # integrate over the range of PAR
    PARin = integrate.simps(NP_Total[dv.PAR_Id[0]:dv.PAR_Id[1]],
                            WL[dv.PAR_Id[0]:dv.PAR_Id[1]])

    PARout = PARin * 86400

    KD_PAR = np.average(KD_out[dv.PAR_Id[0]:dv.PAR_Id[1]])
    # loop over each wavelength and
    # calculate the remote sensing reflectance (Rrs)
    # j = 0
    # while j <= len(WL)-1:
        # rlwn[j], Rrs[j] = OASIM.OASIM(WL[j], Ed_Top[j], T, S, aTotal[j],
        #                               bTotal[j])
        
    rlwn, Rrs = OASIM.OASIM(WL, Ed_Top, T, S, aTotal, bTotal)
        # j = j + 1

    return KD_out, PARout, NP_Total, KD_PAR, Ed_Avg, Rrs, SZA, SZAa


def Declination(JDAY, LAT, LON, yy):

    # function to calculate the declination of the sun
    # based on the time of day

    twopi = 2*np.pi
    DOY = np.remainder(JDAY, 365.)
    HOUR = 24.*(DOY-np.trunc(DOY))
    DOY = np.trunc(DOY)

    # SOLAR DATE
    SOLDATE = twopi*DOY/365.

    # DECLINATION
    DECL = twopi*(0.39637-(22.9133*np.cos(SOLDATE))+(4.02543*np.sin(SOLDATE)) -
                  (0.3872*np.cos(2.*SOLDATE))+(0.052*np.sin(2*SOLDATE)))/360.
    # reflectance
    # time of the day
    tau = twopi*HOUR/24.

   #COSIN OF THETA IN AIR
    CTHTAA = ((np.sin(LAT))*(np.sin(DECL))) - \
        ((np.cos(LAT))*(np.cos(DECL))*(np.cos(tau)))

    #INCIDENCE OF THETA IN AIR
    ITHTAA = np.arccos(CTHTAA)
    # THETA IN WATER
    #Wen Long, here 1.33 is the water refraction index of light
    THTAW = np.arcsin(np.sin(ITHTAA)/1.33)

    #REFLECTANCE, in fraction of 100
    if (ITHTAA > twopi/4.):
        SREFLECT = 100.
    else:
        SREFLECT = np.exp(0.6148+ITHTAA**3)

    #get the date time for input into the pysolar get_azimuth function
    mydate = datetime.datetime(yy, 1, 1, tzinfo=datetime.timezone.utc) +\
        datetime.timedelta(JDAY - 1)

    SZA = get_altitude(LAT, LON, mydate)

    return SZA, SREFLECT, DOY, CTHTAA


def KD_LEE(SZA, a_total, beta_tss):
    m0 = 0.0
    m1 = 4.259
    m2 = 0.52
    m3 = -10.80

    #Calculate KD
    m0 = 1.0+0.005*(90-SZA)
    m4 = m3*a_total
    m5 = 1.-m2*np.exp(m4)

    KD_out = m0 * a_total + m1 * m5 * beta_tss

    return KD_out
