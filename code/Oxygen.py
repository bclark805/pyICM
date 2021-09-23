#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Oxygen module to calculate Dissolved Oxygen dynamics 

from ICM-DOM-PD

for simplified estuarine and coastal ocean modeling

@author: jbclark8
"""
import numpy as np
import Nutrients as nt
import DefVars as dv
import functions as fct


def DOXG(DO2, T, S, Z,
         PN1, PN2, P1, P2, PR1, PR2, B1, B2,
         NT, Uwind, MNLDOC, SOD):
    # DO2 is dissolved oxygen concentration at previous time step
    # T is temperature, S is salinity, Z is layer depth
    # PN is algal preference for NH4
    # P1 and P2 are algal production rates
    # PR1 and PR2 are the algal predation
    # B1 and B2 are algal biomass
    # NT is the amount of nitrification
    # Uwind is the surface wind velocity, only called in surface layer
    # MNLDOC is the total mineralization of DOC
    # SOD is the sediment oxygen demand, only called in bottom layer

    # metabolic rate of algae
    BM1 = dv.BMR1 * fct.get_ft(dv.ktb1, dv.Tb1, T)
    BM2 = dv.BMR2 * fct.get_ft(dv.ktb2, dv.Tb2, T)

    # this is the fraction of algal predation and metabolism
    # that is lost to POC and DOC
    FRDO = 1 - dv.FCD

    # algal photorespiration and metabolism rates
    CP1 = P1 * dv.PRSP1 + BM1
    CP2 = P2 * dv.PRSP2 + BM2

    # production of oxygen from algal production
    DOR1 = ((1.3 - 0.3 * PN1) * P1 - FRDO * CP1) * dv.AOCR * B1
    DOR2 = ((1.3 - 0.3 * PN2) * P2 - FRDO * CP2) * dv.AOCR * B2

    # loss of oxygen from algal predation
    DOPR = dv.FDOP * (PR1 + PR2) * dv.AOCR

    # loss of oxygen from remineralization of DOC
    DDOC = MNLDOC * dv.AOCR

    # change in DO2 over time from biochemical processes
    DTDO2 = (DOR1 + DOR2 - DOPR - DDOC - NT) / 86400.

    # calculate the rearation, if Z is given
    if Z > 0.:
        RearDO2 = Rearation(DO2, T, S, Uwind, Z)
        DTDO2 = DTDO2 + RearDO2/86400.

    # calculate SOD if bottom layer and SOD is true
    if SOD:
        DTDO2 = DTDO2 - SOD

    return DTDO2


def Rearation(DO2, T, S, W, Z0):
    # functiom to calculate the reaeration of DO2 in
    # the surface layer from wind
    # velocity parameterized from Stumm and Morgan, 3rd ed., Ho et al. 2006,
    # Wanninkhof, 2014
    # W is wind velocity
    # Z0 is the surface layer thickness

    # reaeration coefficients
    a = 0.08
    b = 1.0
    c = 1.5

    F1 = a * (b * W) ** c

    # chlorine from salinity
    cl = S / 1.80655
    # ratio of kinematic viscosity of pure water at 20 C to
    # kinematic viscosity at modeled T and S
    Rnu = 0.54 + 0.7 * T / 30. - 0.07 * S / 35

    # reaeration velocity
    KRDO = F1 * Rnu

    # oxygen saturation at  given T and S
    DOS = 14.5532 + T * (0.0054258 * T - 0.38217) - cl\
        * (0.1665 + T * (9.796e-5 * T - 5.866e-3))

    RearDO = KRDO / Z0 * (DOS - DO2)

    return RearDO















    
    
    
    

