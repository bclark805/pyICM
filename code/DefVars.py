#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:45:19 2021

@author: jbclark8
"""


# function to define global variables that will be used throughout the program
import numpy as np
import pandas as pd


def InitializeVariables():
    errorcode = 1
    # physical parameters
    global Z
    Z = 10.  # depth in meters
    global X
    X = 1.e4  # x width in meters
    global Y
    Y = 1.e4  # y width in meters
    global Q
    Q = 10.  # flow rate in m^3 s^-1 that will be forced into the model
    global mu_diff
    mu_diff = 1.e-6  # molecular diffusion rate
    global TurbDiff
    TurbDiff = 0.  # turbulence diffusion rate, will be calculated from surface mixing via wind velocity


##################################
# Below are coefficients that are fixed for each algae but will largely determine how and when
# They grow depending on the light, nutrient, and temperature conditions
    # half saturation of inorganic nitrogen uptake by algae

    # algal basal metabolic coefficient, same for both phytoplankton
    global ktb1
    ktb1 = 0.032
    global ktb2
    ktb2 = 0.032

    # reference temperature for algal metabolism
    global Tb1
    Tb1 = 20.
    global Tb2
    Tb2 = 20.

    # algal predation temperature coefficient
    global ktpr1
    ktpr1 = 0.032
    global ktpr2
    ktpr2 = 0.032
    # reference temperature for algal predation
    global Tpr1
    Tpr1 = 20.
    global Tpr2
    Tpr2 = 20.

    # maximum photosynthetic rate for algae 1 and algae 2 (g C g chla^-1 d^-1)
    global PM1
    PM1 = 300
    global PM2
    PM2 = 350

    # temperature scaled maximum photosynthetic rate, calculated
    global PP1
    PP1 = 0.
    global PP2
    PP2 = 0.

    # basal metabolic rate for algae 1 and 2
    global BMR1
    BMR1 = 0.01
    global BMR2
    BMR2 = 0.02

    # basal algal predation rate from Kimmet et al. 2006 for zooplankton in CBay
    global BPR1
    BPR1 = 1.5
    global BPR2
    BPR2 = 1.5

    # fraaction of predation that is lost as NH4
    global FNP1
    FNP1 = 0.1
    global FNP2
    FNP2 = 0.1

    # algaal nitrogeen to carbon ratio, fixed
    global ANC1
    ANC1 = 0.135
    global ANC2
    ANC2 = 0.175

    # slope of P vs I curve
    global alpha1
    alpha1 = 8.
    global alpha2
    alpha2 = 8.

    # nitrogen limitaation half saaturation
    global khn1
    khn1 = 0.025
    global khn2
    khn2 = 0.025

    # carbon to chla ratio g C g Chl^-1
    global cchla1
    cchla1 = 50.
    global cchla2
    cchla2 = 50.

    # algal respiration coefficient
    global PRSP1
    PRSP1 = 0.25
    global PRSP2
    PRSP2 = 0.25

    # algal sinking velocity (m/d)
    global WB1
    WB1 = 1.0
    global WB2
    WB2 = 1.0

    # algal temperature coefficients
    global kta1a
    kta1a = 0.0018  # sub-optimal temperature coefficient for algae 1 production
    global kta1b
    kta1b = 0.006  # super-optiamal temperature coefficient for algae 1 production
    global kta2a
    kta2a = 0.0035  # sub-optimal temperature coefficient for algae 2 production
    global kta2b
    kta2b = 0.000  # super-optiamal temperature coefficient for algae 2 production

    # reference temperature for algal production
    global TrefB1
    TrefB1 = 16.0
    global TrefB2
    TrefB2 = 35.0
##############################################################################
    # parameters for Carbon
    # algal DOC exudation fraction
    global FCD
    FCD = 0.1

    # oxygen to carbon ratio for production and respiration
    global AOCR
    AOCR = 2.666666667

    # Fraction of algal predation that is respiration
    global FDOP
    FDOP = 0.2

###############################################################################
    # parameters for light
    # the Ids of the Wavelength vector for 400 and 700 nm
    global PAR_Id
    PAR_Id = [0, 0]
    # Avogogadros Number
    global Avo_N
    Avo_N = 6.0221409e23
    # converting from Watts to Photons, coefficients
    # https://www.berthold.com/en/bio/how-do-i-convert-irradiance-photon-flux
    # Convert to number of photons per second
    global NP_sec
    NP_sec = 5.03e15
    # Conver to number of photons per day
    global NP_day
    NP_day = 4.3459e20

###############################################################################
    # parameters related to other things

    global WISS
    WISS = 0.

    errorcode = 0

###############################################################################

    # parameters related to nitrogen biogeochemical processes

    # nitrification
    global ktnt1
    ktnt1 = 0.0  # suboptimal nitrification coefficient
    global ktnt2
    ktnt2 = 0.0  # superoptimal nitrification coefficient
    global TMNT
    TMNT = 20.0  # temperature of maximum nitrification
    global KHONT
    KHONT = 0.00  # oxygen half saturation limitation for nitrification
    global KHNNT
    KHNNT = 0.0  # NH4 half saturation limitation on nitrifcation
    global NTM
    NTM = 0.5  # maximum nitrification rate

    # denitrification

    return errorcode
